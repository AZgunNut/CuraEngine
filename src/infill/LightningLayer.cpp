// Copyright (c) 2022 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/LightningLayer.h" //The class we're implementing.

#include <algorithm>
#include <cmath>
#include <iterator> // advance
#include <limits>
#include <vector>

#include "geometry/OpenPolyline.h"
#include "infill/LightningDistanceField.h"
#include "infill/LightningTreeNode.h"
#include "sliceDataStorage.h"
#include "utils/SparsePointGridInclusive.h"
#include "utils/linearAlg2D.h"

using namespace cura;

namespace
{
constexpr size_t no_polyline = std::numeric_limits<size_t>::max();
constexpr double pi = 3.14159265358979323846;

struct ExperimentalSource
{
    Point2LL point;
    Point2LL outward_direction;
    size_t polyline_idx{ no_polyline };
    coord_t source_length{ 0 };
};

struct ExperimentalSegment
{
    Point2LL a;
    Point2LL b;
    size_t polyline_idx{ no_polyline };
};

struct ExperimentalVariantConfig
{
    bool allow_same_polyline{ false };
    bool actual_leaf_sources{ false };
    bool actual_tree_targets{ false };
    bool use_preclip_geometry{ false };
    bool angular_scoring{ false };
    bool use_arc{ false };
    bool dynamic_reach{ false };
    bool second_target_fallback{ false };
    coord_t max_reach{ 12000 };
    double minimum_distance_widths{ 2.0 };
};

ExperimentalVariantConfig getVariantConfig(const int variant)
{
    ExperimentalVariantConfig config;
    switch (variant)
    {
    case 1:
        config.allow_same_polyline = true;
        break;
    case 2:
        config.actual_leaf_sources = true;
        break;
    case 3:
        config.actual_tree_targets = true;
        break;
    case 4:
        config.use_preclip_geometry = true;
        break;
    case 5:
        config.angular_scoring = true;
        break;
    case 6:
        config.use_arc = true;
        break;
    case 7:
        config.max_reach = 20000;
        break;
    case 8:
        config.dynamic_reach = true;
        break;
    case 9:
        config.second_target_fallback = true;
        break;
    case 10:
        config.minimum_distance_widths = 0.75;
        break;
    case 11:
        config.actual_leaf_sources = true;
        config.actual_tree_targets = true;
        break;
    case 12:
        config.allow_same_polyline = true;
        config.angular_scoring = true;
        break;
    case 13:
        config.actual_tree_targets = true;
        config.angular_scoring = true;
        break;
    case 14:
        config.actual_leaf_sources = true;
        config.allow_same_polyline = true;
        break;
    case 15:
        config.actual_leaf_sources = true;
        config.use_preclip_geometry = true;
        break;
    case 16:
        config.allow_same_polyline = true;
        config.use_arc = true;
        break;
    case 17:
        config.actual_tree_targets = true;
        config.use_arc = true;
        break;
    case 18:
        config.angular_scoring = true;
        config.use_arc = true;
        break;
    case 19:
        config.dynamic_reach = true;
        config.angular_scoring = true;
        break;
    case 20:
        config.second_target_fallback = true;
        config.use_arc = true;
        break;
    default:
        break;
    }
    return config;
}

coord_t polylineLength(const OpenPolyline& line)
{
    coord_t total = 0;
    for (size_t idx = 1; idx < line.size(); ++idx)
    {
        total += vSize(line[idx] - line[idx - 1]);
    }
    return total;
}

std::vector<ExperimentalSegment> collectTreeSegments(const std::vector<LightningTreeNodeSPtr>& tree_roots)
{
    std::vector<ExperimentalSegment> segments;
    for (const LightningTreeNodeSPtr& tree : tree_roots)
    {
        tree->visitBranches(
            [&segments](const Point2LL& parent, const Point2LL& child)
            {
                segments.push_back({ parent, child, no_polyline });
            });
    }
    return segments;
}

std::vector<ExperimentalSource> collectActualLeaves(const std::vector<ExperimentalSegment>& tree_segments)
{
    std::vector<ExperimentalSource> leaves;
    for (const ExperimentalSegment& segment : tree_segments)
    {
        bool child_is_parent = false;
        for (const ExperimentalSegment& candidate : tree_segments)
        {
            if (candidate.a == segment.b)
            {
                child_is_parent = true;
                break;
            }
        }
        if (child_is_parent)
        {
            continue;
        }

        bool already_added = false;
        for (const ExperimentalSource& source : leaves)
        {
            if (source.point == segment.b)
            {
                already_added = true;
                break;
            }
        }
        if (! already_added)
        {
            leaves.push_back({ segment.b, segment.b - segment.a, no_polyline, vSize(segment.b - segment.a) });
        }
    }
    return leaves;
}

std::vector<ExperimentalSource> collectPolylineSources(const OpenLinesSet& lines)
{
    std::vector<ExperimentalSource> sources;
    for (size_t idx = 0; idx < lines.size(); ++idx)
    {
        const OpenPolyline& line = lines[idx];
        if (line.size() < 2)
        {
            continue;
        }
        sources.push_back({ line.front(), line.front() - line[1], idx, polylineLength(line) });
    }
    return sources;
}

std::vector<ExperimentalSegment> collectPolylineSegments(const OpenLinesSet& lines)
{
    std::vector<ExperimentalSegment> segments;
    for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx)
    {
        const OpenPolyline& line = lines[line_idx];
        for (size_t segment_idx = 1; segment_idx < line.size(); ++segment_idx)
        {
            segments.push_back({ line[segment_idx - 1], line[segment_idx], line_idx });
        }
    }
    return segments;
}

double alignmentCosine(const Point2LL& direction, const Point2LL& candidate_vector)
{
    const double ax = static_cast<double>(direction.X);
    const double ay = static_cast<double>(direction.Y);
    const double bx = static_cast<double>(candidate_vector.X);
    const double by = static_cast<double>(candidate_vector.Y);
    const double alen = std::hypot(ax, ay);
    const double blen = std::hypot(bx, by);
    if (alen < 1.0 || blen < 1.0)
    {
        return 0.0;
    }
    return (ax * bx + ay * by) / (alen * blen);
}

OpenLinesSet makeConnector(const Point2LL& source, const Point2LL& target, const Point2LL& outward_direction, const coord_t line_width, const bool use_arc)
{
    OpenLinesSet connector;
    if (! use_arc)
    {
        connector.addSegment(source, target);
        return connector;
    }

    const double dx = static_cast<double>(target.X - source.X);
    const double dy = static_cast<double>(target.Y - source.Y);
    const double chord = std::hypot(dx, dy);
    if (chord < 1.0)
    {
        return connector;
    }

    double nx = -dy / chord;
    double ny = dx / chord;
    const double cross = static_cast<double>(outward_direction.X) * dy - static_cast<double>(outward_direction.Y) * dx;
    if (cross < 0.0)
    {
        nx = -nx;
        ny = -ny;
    }
    const double bulge = std::min(chord * 0.25, static_cast<double>(line_width) * 4.0);
    const double control_x = (static_cast<double>(source.X) + static_cast<double>(target.X)) * 0.5 + nx * bulge;
    const double control_y = (static_cast<double>(source.Y) + static_cast<double>(target.Y)) * 0.5 + ny * bulge;

    OpenPolyline curve;
    constexpr int steps = 10;
    for (int step = 0; step <= steps; ++step)
    {
        const double t = static_cast<double>(step) / static_cast<double>(steps);
        const double omt = 1.0 - t;
        const double x = omt * omt * static_cast<double>(source.X) + 2.0 * omt * t * control_x + t * t * static_cast<double>(target.X);
        const double y = omt * omt * static_cast<double>(source.Y) + 2.0 * omt * t * control_y + t * t * static_cast<double>(target.Y);
        curve.push_back(Point2LL(static_cast<coord_t>(std::llround(x)), static_cast<coord_t>(std::llround(y))));
    }
    connector.push_back(std::move(curve));
    return connector;
}

OpenPolyline smoothPolyline(const OpenPolyline& source)
{
    OpenPolyline result;
    if (source.size() < 3)
    {
        result = source;
        return result;
    }

    for (size_t idx = 0; idx < source.size(); ++idx)
    {
        const size_t begin = idx > 2 ? idx - 2 : 0;
        const size_t end = std::min(source.size() - 1, idx + 2);
        long long sum_x = 0;
        long long sum_y = 0;
        size_t count = 0;
        for (size_t sample = begin; sample <= end; ++sample)
        {
            sum_x += source[sample].X;
            sum_y += source[sample].Y;
            ++count;
        }
        result.push_back(Point2LL(static_cast<coord_t>(sum_x / static_cast<long long>(count)), static_cast<coord_t>(sum_y / static_cast<long long>(count))));
    }
    result.front() = source.front();
    result.back() = source.back();
    return result;
}

bool deterministicReverseNeeded(const OpenPolyline& line)
{
    if (line.size() < 2)
    {
        return false;
    }
    if (line.front().X != line.back().X)
    {
        return line.front().X > line.back().X;
    }
    return line.front().Y > line.back().Y;
}

bool nearestOtherPoint(
    const Point2LL& point,
    const OpenLinesSet& original_lines,
    const size_t own_line,
    Point2LL& nearest,
    coord_t& nearest_distance)
{
    bool found = false;
    nearest_distance = std::numeric_limits<coord_t>::max();
    for (size_t line_idx = 0; line_idx < original_lines.size(); ++line_idx)
    {
        if (line_idx == own_line)
        {
            continue;
        }
        const OpenPolyline& line = original_lines[line_idx];
        for (size_t segment_idx = 1; segment_idx < line.size(); ++segment_idx)
        {
            const Point2LL candidate = LinearAlg2D::getClosestOnLineSegment(point, line[segment_idx - 1], line[segment_idx]);
            const coord_t distance = vSize(candidate - point);
            if (distance < nearest_distance)
            {
                nearest = candidate;
                nearest_distance = distance;
                found = true;
            }
        }
    }
    return found;
}

OpenLinesSet makeSmoothedContinuousPaths(const OpenLinesSet& original_lines, const Shape& limit_to_outline, const coord_t line_width, const bool sine_offset)
{
    OpenLinesSet generated;
    constexpr coord_t sine_wavelength = 10000; // 10 mm full wavelength.

    for (size_t line_idx = 0; line_idx < original_lines.size(); ++line_idx)
    {
        OpenPolyline guide = original_lines[line_idx];
        if (guide.size() < 2)
        {
            continue;
        }

        if (deterministicReverseNeeded(guide))
        {
            std::reverse(guide.begin(), guide.end());
        }
        guide = smoothPolyline(guide);

        std::vector<Point2LL> side_a;
        std::vector<Point2LL> side_b;
        side_a.reserve(guide.size());
        side_b.reserve(guide.size());

        double cumulative_distance = 0.0;
        double previous_nx = 0.0;
        double previous_ny = 0.0;
        bool have_previous_normal = false;

        for (size_t idx = 0; idx < guide.size(); ++idx)
        {
            if (idx > 0)
            {
                cumulative_distance += static_cast<double>(vSize(guide[idx] - guide[idx - 1]));
            }

            const size_t before_idx = idx > 2 ? idx - 2 : 0;
            const size_t after_idx = std::min(guide.size() - 1, idx + 2);
            const double dx = static_cast<double>(guide[after_idx].X - guide[before_idx].X);
            const double dy = static_cast<double>(guide[after_idx].Y - guide[before_idx].Y);
            const double tangent_length = std::hypot(dx, dy);

            double nx = 0.0;
            double ny = 0.0;
            if (tangent_length > 0.5)
            {
                nx = -dy / tangent_length;
                ny = dx / tangent_length;
            }
            if (! have_previous_normal)
            {
                if (nx < 0.0 || (std::abs(nx) < 0.000001 && ny < 0.0))
                {
                    nx = -nx;
                    ny = -ny;
                }
                previous_nx = nx;
                previous_ny = ny;
                have_previous_normal = true;
            }
            else
            {
                if (nx * previous_nx + ny * previous_ny < 0.0)
                {
                    nx = -nx;
                    ny = -ny;
                }
                previous_nx = nx;
                previous_ny = ny;
            }

            double desired_offset = static_cast<double>(line_width) * 4.0;
            if (sine_offset)
            {
                const double phase = 2.0 * pi * cumulative_distance / static_cast<double>(sine_wavelength);
                desired_offset = static_cast<double>(line_width) * 2.0 * (1.0 - std::cos(phase));
            }

            Point2LL a = guide[idx];
            Point2LL b(
                static_cast<coord_t>(std::llround(static_cast<double>(guide[idx].X) + nx * desired_offset)),
                static_cast<coord_t>(std::llround(static_cast<double>(guide[idx].Y) + ny * desired_offset)));

            // Cramped-area rule: when the two intended averaged tracks do not fit
            // between neighboring Lightning, put them side-by-side around the
            // median gap for this local area, then smoothly return to the normal
            // offset as space opens again.
            Point2LL other;
            coord_t other_distance = 0;
            if (nearestOtherPoint(guide[idx], original_lines, line_idx, other, other_distance)
                && other_distance > line_width && static_cast<double>(other_distance) < desired_offset + static_cast<double>(line_width) * 2.0)
            {
                const double gx = static_cast<double>(other.X - guide[idx].X);
                const double gy = static_cast<double>(other.Y - guide[idx].Y);
                const double gap = std::hypot(gx, gy);
                if (gap > 1.0)
                {
                    const double ux = gx / gap;
                    const double uy = gy / gap;
                    const double median_x = (static_cast<double>(guide[idx].X) + static_cast<double>(other.X)) * 0.5;
                    const double median_y = (static_cast<double>(guide[idx].Y) + static_cast<double>(other.Y)) * 0.5;
                    const double half_spacing = static_cast<double>(line_width) * 0.5;
                    a = Point2LL(
                        static_cast<coord_t>(std::llround(median_x - ux * half_spacing)),
                        static_cast<coord_t>(std::llround(median_y - uy * half_spacing)));
                    b = Point2LL(
                        static_cast<coord_t>(std::llround(median_x + ux * half_spacing)),
                        static_cast<coord_t>(std::llround(median_y + uy * half_spacing)));
                }
            }

            side_a.push_back(a);
            side_b.push_back(b);
        }

        OpenPolyline continuous;
        for (const Point2LL& point : side_a)
        {
            continuous.push_back(point);
        }
        for (auto it = side_b.rbegin(); it != side_b.rend(); ++it)
        {
            continuous.push_back(*it);
        }
        if (! continuous.empty() && continuous.back() != continuous.front())
        {
            continuous.push_back(continuous.front());
        }
        generated.push_back(std::move(continuous));
    }

    return limit_to_outline.intersection(generated);
}

} // namespace

coord_t LightningLayer::getWeightedDistance(const Point2LL& boundary_loc, const Point2LL& unsupported_location)
{
    return vSize(boundary_loc - unsupported_location);
}

Point2LL GroundingLocation::p() const
{
    if (tree_node != nullptr)
    {
        return tree_node->getLocation();
    }
    else
    {
        assert(boundary_location);
        return boundary_location->p();
    }
}

void LightningLayer::fillLocator(SparseLightningTreeNodeGrid& tree_node_locator)
{
    std::function<void(LightningTreeNodeSPtr)> add_node_to_locator_func = [&tree_node_locator](LightningTreeNodeSPtr node)
    {
        tree_node_locator.insert(node->getLocation(), node);
    };
    for (auto& tree : tree_roots)
    {
        tree->visitNodes(add_node_to_locator_func);
    }
}

void LightningLayer::generateNewTrees(
    const Shape& current_overhang,
    const Shape& current_outlines,
    const LocToLineGrid& outlines_locator,
    const coord_t supporting_radius,
    const coord_t wall_supporting_radius)
{
    LightningDistanceField distance_field(supporting_radius, current_outlines, current_overhang);

    SparseLightningTreeNodeGrid tree_node_locator(locator_cell_size);
    fillLocator(tree_node_locator);

    // Until no more points need to be added to support all:
    // Determine next point from tree/outline areas via distance-field
    Point2LL unsupported_location;
    while (distance_field.tryGetNextPoint(&unsupported_location))
    {
        GroundingLocation grounding_loc
            = getBestGroundingLocation(unsupported_location, current_outlines, outlines_locator, supporting_radius, wall_supporting_radius, tree_node_locator);

        LightningTreeNodeSPtr new_parent;
        LightningTreeNodeSPtr new_child;
        attach(unsupported_location, grounding_loc, new_child, new_parent);
        tree_node_locator.insert(new_child->getLocation(), new_child);
        if (new_parent)
        {
            tree_node_locator.insert(new_parent->getLocation(), new_parent);
        }

        // update distance field
        distance_field.update(grounding_loc.p(), unsupported_location);
    }
}

GroundingLocation LightningLayer::getBestGroundingLocation(
    const Point2LL& unsupported_location,
    const Shape& current_outlines,
    const LocToLineGrid& outline_locator,
    const coord_t supporting_radius,
    const coord_t wall_supporting_radius,
    const SparseLightningTreeNodeGrid& tree_node_locator,
    const LightningTreeNodeSPtr& exclude_tree)
{
    ClosestPointPolygon cpp = PolygonUtils::findClosest(unsupported_location, current_outlines);
    Point2LL node_location = cpp.p();
    const coord_t within_dist = vSize(node_location - unsupported_location);

    PolygonsPointIndex dummy;

    LightningTreeNodeSPtr sub_tree{ nullptr };
    coord_t current_dist = getWeightedDistance(node_location, unsupported_location);
    if (current_dist >= wall_supporting_radius) // Only reconnect tree roots to other trees if they are not already close to the outlines.
    {
        auto candidate_trees = tree_node_locator.getNearbyVals(unsupported_location, std::min(current_dist, within_dist));
        for (auto& candidate_wptr : candidate_trees)
        {
            auto candidate_sub_tree = candidate_wptr.lock();
            if ((candidate_sub_tree && candidate_sub_tree != exclude_tree) && ! (exclude_tree && exclude_tree->hasOffspring(candidate_sub_tree))
                && ! PolygonUtils::polygonCollidesWithLineSegment(unsupported_location, candidate_sub_tree->getLocation(), outline_locator, &dummy))
            {
                const coord_t candidate_dist = candidate_sub_tree->getWeightedDistance(unsupported_location, supporting_radius);
                if (candidate_dist < current_dist)
                {
                    current_dist = candidate_dist;
                    sub_tree = candidate_sub_tree;
                }
            }
        }
    }

    if (! sub_tree)
    {
        return GroundingLocation{ nullptr, cpp };
    }
    else
    {
        return GroundingLocation{ sub_tree, std::optional<ClosestPointPolygon>() };
    }
}

bool LightningLayer::attach(const Point2LL& unsupported_location, const GroundingLocation& grounding_loc, LightningTreeNodeSPtr& new_child, LightningTreeNodeSPtr& new_root)
{
    // Update trees & distance fields.
    if (grounding_loc.boundary_location)
    {
        new_root = LightningTreeNode::create(grounding_loc.p(), std::make_optional(grounding_loc.p()));
        new_child = new_root->addChild(unsupported_location);
        tree_roots.push_back(new_root);
        return true;
    }
    else
    {
        new_child = grounding_loc.tree_node->addChild(unsupported_location);
        return false;
    }
}

void LightningLayer::reconnectRoots(
    std::vector<LightningTreeNodeSPtr>& to_be_reconnected_tree_roots,
    const Shape& current_outlines,
    const LocToLineGrid& outline_locator,
    const coord_t supporting_radius,
    const coord_t wall_supporting_radius)
{
    constexpr coord_t tree_connecting_ignore_offset = 100;

    SparseLightningTreeNodeGrid tree_node_locator(locator_cell_size);
    fillLocator(tree_node_locator);

    const coord_t within_max_dist = outline_locator.getCellSize() * 2;
    for (auto root_ptr : to_be_reconnected_tree_roots)
    {
        auto old_root_it = std::find(tree_roots.begin(), tree_roots.end(), root_ptr);

        if (root_ptr->getLastGroundingLocation())
        {
            const Point2LL& ground_loc = root_ptr->getLastGroundingLocation().value();
            if (ground_loc != root_ptr->getLocation())
            {
                Point2LL new_root_pt;
                if (PolygonUtils::lineSegmentPolygonsIntersection(root_ptr->getLocation(), ground_loc, current_outlines, outline_locator, new_root_pt, within_max_dist))
                {
                    auto new_root = LightningTreeNode::create(new_root_pt, new_root_pt);
                    root_ptr->addChild(new_root);
                    new_root->reroot();

                    tree_node_locator.insert(new_root->getLocation(), new_root);
                    *old_root_it = std::move(new_root); // replace old root with new root
                    continue;
                }
            }
        }

        const coord_t tree_connecting_ignore_width
            = wall_supporting_radius - tree_connecting_ignore_offset; // Ideally, the boundary size in which the valence rule is ignored would be configurable.
        GroundingLocation ground
            = getBestGroundingLocation(root_ptr->getLocation(), current_outlines, outline_locator, supporting_radius, tree_connecting_ignore_width, tree_node_locator, root_ptr);
        if (ground.boundary_location)
        {
            if (ground.boundary_location.value().p() == root_ptr->getLocation())
            {
                continue; // Already on the boundary.
            }

            auto new_root = LightningTreeNode::create(ground.p(), ground.p());
            auto attach_ptr = root_ptr->closestNode(new_root->getLocation());
            attach_ptr->reroot();

            new_root->addChild(attach_ptr);
            tree_node_locator.insert(new_root->getLocation(), new_root);

            *old_root_it = std::move(new_root); // replace old root with new root
        }
        else
        {
            assert(ground.tree_node);
            assert(ground.tree_node != root_ptr);
            assert(! root_ptr->hasOffspring(ground.tree_node));
            assert(! ground.tree_node->hasOffspring(root_ptr));

            auto attach_ptr = root_ptr->closestNode(ground.tree_node->getLocation());
            attach_ptr->reroot();

            ground.tree_node->addChild(attach_ptr);

            // remove old root
            *old_root_it = std::move(tree_roots.back());
            tree_roots.pop_back();
        }
    }
}

// Returns 'added someting'.
OpenLinesSet LightningLayer::convertToLines(const Shape& limit_to_outline, const coord_t line_width) const
{
    OpenLinesSet raw_lines;
    if (tree_roots.empty())
    {
        return raw_lines;
    }

    for (const LightningTreeNodeSPtr& tree : tree_roots)
    {
        tree->convertToPolylines(raw_lines, line_width);
    }

    OpenLinesSet result_lines = limit_to_outline.intersection(raw_lines);

    // Variant 0 is the stock Cura control.
    if (experimental_variant <= 0 || experimental_variant > 22)
    {
        return result_lines;
    }

    // Variants 21 and 22 deliberately replace normal Lightning output with a
    // smoothed, paired, continuous out-and-back path. 22 modulates the offset
    // from 0..4 extrusion widths on a phase-locked 10 mm sine.
    if (experimental_variant == 21 || experimental_variant == 22)
    {
        return makeSmoothedContinuousPaths(result_lines, limit_to_outline, line_width, experimental_variant == 22);
    }

    const ExperimentalVariantConfig config = getVariantConfig(experimental_variant);
    const OpenLinesSet& working_lines = config.use_preclip_geometry ? raw_lines : result_lines;
    const std::vector<ExperimentalSegment> tree_segments = collectTreeSegments(tree_roots);
    const std::vector<ExperimentalSource> sources
        = config.actual_leaf_sources ? collectActualLeaves(tree_segments) : collectPolylineSources(working_lines);
    const std::vector<ExperimentalSegment> targets
        = config.actual_tree_targets ? tree_segments : collectPolylineSegments(working_lines);

    const coord_t minimum_distance = static_cast<coord_t>(std::llround(static_cast<double>(line_width) * config.minimum_distance_widths));

    for (const ExperimentalSource& source : sources)
    {
        coord_t source_max_reach = config.max_reach;
        if (config.dynamic_reach)
        {
            source_max_reach = std::max<coord_t>(6000, std::min<coord_t>(30000, source.source_length * 2));
        }

        Point2LL best_target;
        Point2LL second_target;
        double best_score = std::numeric_limits<double>::max();
        double second_score = std::numeric_limits<double>::max();
        coord_t best_distance = source_max_reach + 1;
        coord_t second_distance = source_max_reach + 1;
        double best_alignment = -1.0;
        bool found_best = false;
        bool found_second = false;

        for (const ExperimentalSegment& target_segment : targets)
        {
            if (! config.allow_same_polyline && source.polyline_idx != no_polyline && target_segment.polyline_idx == source.polyline_idx)
            {
                continue;
            }

            // When same-polyline targets are enabled, skip the immediate source
            // neighborhood so the connector cannot simply fold onto itself.
            if (config.allow_same_polyline && source.polyline_idx != no_polyline && target_segment.polyline_idx == source.polyline_idx
                && (target_segment.a == source.point || target_segment.b == source.point))
            {
                continue;
            }

            const Point2LL candidate = LinearAlg2D::getClosestOnLineSegment(source.point, target_segment.a, target_segment.b);
            const coord_t distance = vSize(candidate - source.point);
            if (distance < minimum_distance || distance > source_max_reach)
            {
                continue;
            }

            const double alignment = alignmentCosine(source.outward_direction, candidate - source.point);
            const double angular_multiplier = config.angular_scoring ? (1.0 + 0.75 * (1.0 - alignment) * 0.5) : 1.0;
            const double score = static_cast<double>(distance) * angular_multiplier;

            if (score < best_score)
            {
                second_score = best_score;
                second_target = best_target;
                second_distance = best_distance;
                found_second = found_best;

                best_score = score;
                best_target = candidate;
                best_distance = distance;
                best_alignment = alignment;
                found_best = true;
            }
            else if (score < second_score)
            {
                second_score = score;
                second_target = candidate;
                second_distance = distance;
                found_second = true;
            }
        }

        if (! found_best)
        {
            continue;
        }

        Point2LL selected_target = best_target;
        coord_t selected_distance = best_distance;
        if (config.second_target_fallback && best_alignment < 0.0 && found_second)
        {
            selected_target = second_target;
            selected_distance = second_distance;
        }
        if (selected_distance > source_max_reach)
        {
            continue;
        }

        OpenLinesSet connector = makeConnector(source.point, selected_target, source.outward_direction, line_width, config.use_arc);
        OpenLinesSet clipped_connector = limit_to_outline.intersection(connector);
        result_lines.push_back(std::move(clipped_connector));
    }

    return result_lines;
}
