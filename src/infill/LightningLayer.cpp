// Copyright (c) 2022 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/LightningLayer.h"

#include <algorithm>
#include <cmath>
#include <functional>
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
constexpr double pi = 3.14159265358979323846;
constexpr double sine_wavelength = 10000.0; // 10 mm in Cura coordinate units.

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

long double orientation(const Point2LL& a, const Point2LL& b, const Point2LL& c)
{
    const long double abx = static_cast<long double>(b.X) - static_cast<long double>(a.X);
    const long double aby = static_cast<long double>(b.Y) - static_cast<long double>(a.Y);
    const long double acx = static_cast<long double>(c.X) - static_cast<long double>(a.X);
    const long double acy = static_cast<long double>(c.Y) - static_cast<long double>(a.Y);
    return abx * acy - aby * acx;
}

bool properSegmentsCross(const Point2LL& a, const Point2LL& b, const Point2LL& c, const Point2LL& d)
{
    const long double o1 = orientation(a, b, c);
    const long double o2 = orientation(a, b, d);
    const long double o3 = orientation(c, d, a);
    const long double o4 = orientation(c, d, b);
    return ((o1 > 0.0L && o2 < 0.0L) || (o1 < 0.0L && o2 > 0.0L))
        && ((o3 > 0.0L && o4 < 0.0L) || (o3 < 0.0L && o4 > 0.0L));
}

bool segmentCrossesPolyline(const Point2LL& a, const Point2LL& b, const OpenPolyline& line)
{
    for (size_t idx = 1; idx < line.size(); ++idx)
    {
        if (properSegmentsCross(a, b, line[idx - 1], line[idx]))
        {
            return true;
        }
    }
    return false;
}

bool segmentCrossesLines(const Point2LL& a, const Point2LL& b, const OpenLinesSet& lines, const size_t skip_line = std::numeric_limits<size_t>::max())
{
    for (size_t line_idx = 0; line_idx < lines.size(); ++line_idx)
    {
        if (line_idx == skip_line)
        {
            continue;
        }
        if (segmentCrossesPolyline(a, b, lines[line_idx]))
        {
            return true;
        }
    }
    return false;
}

bool segmentCrossesPaths(const Point2LL& a, const Point2LL& b, const std::vector<OpenPolyline>& paths)
{
    for (const OpenPolyline& path : paths)
    {
        if (segmentCrossesPolyline(a, b, path))
        {
            return true;
        }
    }
    return false;
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

OpenPolyline smoothPolyline(const OpenPolyline& source, const double smoothing)
{
    if (source.size() < 3 || smoothing <= 0.0)
    {
        return source;
    }

    const double strength = std::max(0.0, std::min(1.0, smoothing / 100.0));
    const size_t max_radius = std::min<size_t>(50, source.size() - 1);
    const size_t radius = std::max<size_t>(1, static_cast<size_t>(std::llround(1.0 + strength * strength * static_cast<double>(max_radius - 1))));

    OpenPolyline result;
    result.reserve(source.size());
    for (size_t idx = 0; idx < source.size(); ++idx)
    {
        if (idx == 0 || idx + 1 == source.size())
        {
            result.push_back(source[idx]);
            continue;
        }

        const size_t begin = idx > radius ? idx - radius : 0;
        const size_t end = std::min(source.size() - 1, idx + radius);
        long double sum_x = 0.0L;
        long double sum_y = 0.0L;
        size_t count = 0;
        for (size_t sample = begin; sample <= end; ++sample)
        {
            sum_x += static_cast<long double>(source[sample].X);
            sum_y += static_cast<long double>(source[sample].Y);
            ++count;
        }

        const long double avg_x = sum_x / static_cast<long double>(count);
        const long double avg_y = sum_y / static_cast<long double>(count);
        const long double x = static_cast<long double>(source[idx].X) * (1.0L - strength) + avg_x * strength;
        const long double y = static_cast<long double>(source[idx].Y) * (1.0L - strength) + avg_y * strength;
        result.push_back(Point2LL(static_cast<coord_t>(std::llround(x)), static_cast<coord_t>(std::llround(y))));
    }
    return result;
}

OpenPolyline makeSafeSmoothedGuide(
    const OpenPolyline& original,
    const OpenLinesSet& original_lines,
    const size_t own_line,
    const double smoothing,
    const std::vector<OpenPolyline>& already_generated)
{
    OpenPolyline candidate = smoothPolyline(original, smoothing);
    if (candidate.size() < 2)
    {
        return candidate;
    }

    OpenPolyline safe;
    safe.reserve(candidate.size());
    safe.push_back(candidate.front());
    for (size_t idx = 1; idx < candidate.size(); ++idx)
    {
        Point2LL accepted = candidate[idx];
        if (segmentCrossesLines(safe.back(), accepted, original_lines, own_line)
            || segmentCrossesPaths(safe.back(), accepted, already_generated)
            || segmentCrossesPolyline(safe.back(), accepted, safe))
        {
            accepted = original[idx];
        }
        safe.push_back(accepted);
    }
    return safe;
}

Point2LL offsetPoint(const Point2LL& base, const double nx, const double ny, const double distance)
{
    return Point2LL(
        static_cast<coord_t>(std::llround(static_cast<double>(base.X) + nx * distance)),
        static_cast<coord_t>(std::llround(static_cast<double>(base.Y) + ny * distance)));
}

bool offsetCandidateIsSafe(
    const Point2LL& base,
    const Point2LL& previous,
    const Point2LL& candidate,
    const OpenLinesSet& original_lines,
    const size_t own_line,
    const OpenPolyline& side_a,
    const OpenPolyline& side_b,
    const std::vector<OpenPolyline>& already_generated)
{
    // Check the printable companion segment and the radial join. Proper endpoint
    // touches are allowed; true crossovers are not.
    if (segmentCrossesLines(previous, candidate, original_lines, own_line)
        || segmentCrossesLines(base, candidate, original_lines, own_line)
        || segmentCrossesPolyline(previous, candidate, side_a)
        || segmentCrossesPolyline(previous, candidate, side_b)
        || segmentCrossesPaths(previous, candidate, already_generated)
        || segmentCrossesPaths(base, candidate, already_generated))
    {
        return false;
    }
    return true;
}

std::vector<OpenPolyline> joinPathsSafely(std::vector<OpenPolyline> paths, const OpenLinesSet& original_lines)
{
    std::vector<OpenPolyline> result;
    while (! paths.empty())
    {
        OpenPolyline current = std::move(paths.front());
        paths.erase(paths.begin());
        if (current.empty())
        {
            continue;
        }

        bool joined = true;
        while (joined && ! paths.empty())
        {
            joined = false;
            size_t best_idx = 0;
            bool best_reverse = false;
            coord_t best_distance = std::numeric_limits<coord_t>::max();

            for (size_t idx = 0; idx < paths.size(); ++idx)
            {
                if (paths[idx].empty())
                {
                    continue;
                }
                for (int reverse = 0; reverse < 2; ++reverse)
                {
                    const Point2LL target = reverse ? paths[idx].back() : paths[idx].front();
                    const coord_t distance = vSize(target - current.back());
                    if (distance >= best_distance)
                    {
                        continue;
                    }
                    if (segmentCrossesLines(current.back(), target, original_lines)
                        || segmentCrossesPolyline(current.back(), target, current)
                        || segmentCrossesPaths(current.back(), target, paths))
                    {
                        continue;
                    }
                    best_distance = distance;
                    best_idx = idx;
                    best_reverse = reverse != 0;
                    joined = true;
                }
            }

            if (joined)
            {
                OpenPolyline next = std::move(paths[best_idx]);
                paths.erase(paths.begin() + static_cast<std::ptrdiff_t>(best_idx));
                if (best_reverse)
                {
                    std::reverse(next.begin(), next.end());
                }
                if (! next.empty() && current.back() != next.front())
                {
                    current.push_back(next.front());
                }
                for (size_t idx = 1; idx < next.size(); ++idx)
                {
                    current.push_back(next[idx]);
                }
            }
        }
        result.push_back(std::move(current));
    }
    return result;
}

OpenLinesSet makeControlledLightningPaths(
    const OpenLinesSet& original_lines,
    const Shape& limit_to_outline,
    const coord_t line_width,
    const double smoothing,
    const double offset_widths,
    const bool single_path)
{
    std::vector<OpenPolyline> generated_paths;
    const double max_offset = static_cast<double>(line_width) * std::max(0.0, std::min(100.0, offset_widths));

    for (size_t line_idx = 0; line_idx < original_lines.size(); ++line_idx)
    {
        OpenPolyline original = original_lines[line_idx];
        if (original.size() < 2)
        {
            continue;
        }
        if (deterministicReverseNeeded(original))
        {
            std::reverse(original.begin(), original.end());
        }

        OpenPolyline side_a = makeSafeSmoothedGuide(original, original_lines, line_idx, smoothing, generated_paths);
        if (side_a.size() < 2)
        {
            continue;
        }

        if (max_offset <= 0.0)
        {
            generated_paths.push_back(std::move(side_a));
            continue;
        }

        OpenPolyline side_b;
        side_b.reserve(side_a.size());
        double cumulative_distance = 0.0;
        double previous_nx = 0.0;
        double previous_ny = 0.0;
        bool have_previous_normal = false;

        for (size_t idx = 0; idx < side_a.size(); ++idx)
        {
            if (idx > 0)
            {
                cumulative_distance += static_cast<double>(vSize(side_a[idx] - side_a[idx - 1]));
            }

            const size_t before_idx = idx > 2 ? idx - 2 : 0;
            const size_t after_idx = std::min(side_a.size() - 1, idx + 2);
            const double dx = static_cast<double>(side_a[after_idx].X - side_a[before_idx].X);
            const double dy = static_cast<double>(side_a[after_idx].Y - side_a[before_idx].Y);
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

            const double phase = 2.0 * pi * cumulative_distance / sine_wavelength;
            double desired_offset = max_offset * 0.5 * (1.0 - std::cos(phase)); // exactly 0..requested max widths.

            Point2LL other;
            coord_t other_distance = 0;
            if (nearestOtherPoint(side_a[idx], original_lines, line_idx, other, other_distance))
            {
                const double safe_gap = std::max(0.0, static_cast<double>(other_distance) - static_cast<double>(line_width));
                desired_offset = std::min(desired_offset, safe_gap);
            }

            const Point2LL previous = side_b.empty() ? side_a[idx] : side_b.back();
            double low = 0.0;
            double high = desired_offset;
            Point2LL accepted = side_a[idx];
            for (int attempt = 0; attempt < 14; ++attempt)
            {
                const double test = attempt == 0 ? high : (low + high) * 0.5;
                const Point2LL candidate = offsetPoint(side_a[idx], nx, ny, test);
                if (offsetCandidateIsSafe(side_a[idx], previous, candidate, original_lines, line_idx, side_a, side_b, generated_paths))
                {
                    accepted = candidate;
                    low = test;
                    if (attempt == 0)
                    {
                        break;
                    }
                }
                else
                {
                    high = test;
                }
            }
            side_b.push_back(accepted);
        }

        if (single_path)
        {
            OpenPolyline continuous;
            continuous.reserve(side_a.size() + side_b.size() + 1);
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
            generated_paths.push_back(std::move(continuous));
        }
        else
        {
            generated_paths.push_back(std::move(side_a));
            generated_paths.push_back(std::move(side_b));
        }
    }

    if (single_path)
    {
        generated_paths = joinPathsSafely(std::move(generated_paths), original_lines);
    }

    OpenLinesSet generated;
    for (OpenPolyline& path : generated_paths)
    {
        if (path.size() >= 2)
        {
            generated.push_back(std::move(path));
        }
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
    assert(boundary_location);
    return boundary_location->p();
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
    if (current_dist >= wall_supporting_radius)
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
    return GroundingLocation{ sub_tree, std::optional<ClosestPointPolygon>() };
}

bool LightningLayer::attach(const Point2LL& unsupported_location, const GroundingLocation& grounding_loc, LightningTreeNodeSPtr& new_child, LightningTreeNodeSPtr& new_root)
{
    if (grounding_loc.boundary_location)
    {
        new_root = LightningTreeNode::create(grounding_loc.p(), std::make_optional(grounding_loc.p()));
        new_child = new_root->addChild(unsupported_location);
        tree_roots.push_back(new_root);
        return true;
    }
    new_child = grounding_loc.tree_node->addChild(unsupported_location);
    return false;
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
                    *old_root_it = std::move(new_root);
                    continue;
                }
            }
        }

        const coord_t tree_connecting_ignore_width = wall_supporting_radius - tree_connecting_ignore_offset;
        GroundingLocation ground
            = getBestGroundingLocation(root_ptr->getLocation(), current_outlines, outline_locator, supporting_radius, tree_connecting_ignore_width, tree_node_locator, root_ptr);
        if (ground.boundary_location)
        {
            if (ground.boundary_location.value().p() == root_ptr->getLocation())
            {
                continue;
            }

            auto new_root = LightningTreeNode::create(ground.p(), ground.p());
            auto attach_ptr = root_ptr->closestNode(new_root->getLocation());
            attach_ptr->reroot();
            new_root->addChild(attach_ptr);
            tree_node_locator.insert(new_root->getLocation(), new_root);
            *old_root_it = std::move(new_root);
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
            *old_root_it = std::move(tree_roots.back());
            tree_roots.pop_back();
        }
    }
}

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
    if (lightning_smoothing <= 0.0 && lightning_offset_widths <= 0.0 && ! lightning_single_path)
    {
        return result_lines;
    }

    return makeControlledLightningPaths(
        result_lines,
        limit_to_outline,
        line_width,
        lightning_smoothing,
        lightning_offset_widths,
        lightning_single_path);
}
