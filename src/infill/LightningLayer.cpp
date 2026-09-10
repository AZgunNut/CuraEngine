// Copyright (c) 2022 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#include "infill/LightningLayer.h" //The class we're implementing.

#include <iterator> // advance
#include <limits>

#include "geometry/OpenPolyline.h"
#include "infill/LightningDistanceField.h"
#include "infill/LightningTreeNode.h"
#include "sliceDataStorage.h"
#include "utils/SparsePointGridInclusive.h"
#include "utils/linearAlg2D.h"

using namespace cura;

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
    else
    {
        return GroundingLocation{ sub_tree, std::optional<ClosestPointPolygon>() };
    }
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

namespace
{
constexpr coord_t looped_lightning_max_distance = 50000; // 50 mm test reach.
constexpr size_t looped_lightning_curve_segments = 18;

struct LoopedLightningTarget
{
    Point2LL point;
    Point2LL tangent;
    coord_t distance;
    bool has_tangent;
};

/*!
 * Find the nearest usable point on another Lightning polyline, not merely the
 * nearest stored vertex.  This is important for the looped-lightning idea:
 * the middle of an existing branch is valid structure and should be available
 * as a landing point for a returning leaf.
 */
LoopedLightningTarget findLoopedLightningTarget(
    const OpenLinesSet& lines,
    const size_t source_line_idx,
    const Point2LL& source,
    const Shape& limit_to_outline)
{
    LoopedLightningTarget best{ PolygonUtils::findClosest(source, limit_to_outline).p(), Point2LL(0, 0), 0, false };
    best.distance = vSize(best.point - source);

    for (size_t candidate_line_idx = 0; candidate_line_idx < lines.size(); ++candidate_line_idx)
    {
        if (candidate_line_idx == source_line_idx)
        {
            continue;
        }

        const OpenPolyline& candidate_line = lines[candidate_line_idx];
        if (candidate_line.size() < 2)
        {
            continue;
        }

        for (size_t segment_idx = 1; segment_idx < candidate_line.size(); ++segment_idx)
        {
            const Point2LL a = candidate_line[segment_idx - 1];
            const Point2LL b = candidate_line[segment_idx];
            const Point2LL ab = b - a;
            const double length_squared = static_cast<double>(ab.X) * static_cast<double>(ab.X)
                                        + static_cast<double>(ab.Y) * static_cast<double>(ab.Y);
            if (length_squared <= 0.0)
            {
                continue;
            }

            const Point2LL ap = source - a;
            double t = (static_cast<double>(ap.X) * static_cast<double>(ab.X)
                      + static_cast<double>(ap.Y) * static_cast<double>(ab.Y)) / length_squared;
            t = std::max(0.0, std::min(1.0, t));
            const Point2LL candidate = a + ab * t;
            const coord_t candidate_distance = vSize(candidate - source);

            if (candidate_distance < best.distance)
            {
                best.point = candidate;
                best.tangent = ab;
                best.distance = candidate_distance;
                best.has_tangent = true;
            }
        }
    }

    return best;
}

/*!
 * Replace a dangling Lightning leaf with a flowing return to nearby existing
 * structure.  The return leaves tangent to the source branch and, when it lands
 * on another branch, arrives tangent to that branch as well.  The result is a
 * connected loop rather than a dead-ended twig plus a straight stitch.
 *
 * The generated return remains part of result_lines, so subsequent closure
 * passes can use an earlier return as a landing structure.  This lets the
 * experimental network grow branch-to-return-to-branch rather than limiting
 * connections to the original Lightning vertices.
 */
void addLoopedLightningClosures(OpenLinesSet& result_lines, const Shape& limit_to_outline, const coord_t line_width)
{
    if (result_lines.empty() || limit_to_outline.empty())
    {
        return;
    }

    const size_t original_line_count = result_lines.size();

    for (size_t line_idx = 0; line_idx < original_line_count; ++line_idx)
    {
        const OpenPolyline& source_line = result_lines[line_idx];
        if (source_line.size() < 2)
        {
            continue;
        }

        const Point2LL source = source_line.front();
        const LoopedLightningTarget target_info = findLoopedLightningTarget(result_lines, line_idx, source, limit_to_outline);
        const Point2LL target = target_info.point;
        const coord_t best_distance = target_info.distance;

        if (best_distance > looped_lightning_max_distance || best_distance < line_width * 2)
        {
            continue;
        }

        const Point2LL source_tangent = source - source_line[1];
        const coord_t source_tangent_length = vSize(source_tangent);
        if (source_tangent_length <= 0)
        {
            continue;
        }

        const Point2LL chord = target - source;
        const coord_t chord_length = vSize(chord);
        if (chord_length <= 0)
        {
            continue;
        }

        const double source_handle_scale = 0.45 * static_cast<double>(best_distance) / static_cast<double>(source_tangent_length);
        const Point2LL control1 = source + source_tangent * source_handle_scale;

        Point2LL control2;
        if (target_info.has_tangent)
        {
            Point2LL target_tangent = target_info.tangent;
            // Pick the tangent direction that makes the Bezier approach the
            // target rather than curl away from it.
            if (target_tangent.X * chord.X + target_tangent.Y * chord.Y < 0)
            {
                target_tangent = Point2LL(-target_tangent.X, -target_tangent.Y);
            }
            const coord_t target_tangent_length = vSize(target_tangent);
            const double target_handle_scale = 0.35 * static_cast<double>(best_distance) / static_cast<double>(target_tangent_length);
            control2 = target - target_tangent * target_handle_scale;
        }
        else
        {
            // Boundary landing: use a broad side sweep so a wall return is still
            // visibly curved instead of degenerating into a straight connector.
            const Point2LL perpendicular(-chord.Y, chord.X);
            const double side = (line_idx % 2 == 0) ? 1.0 : -1.0;
            control2 = target - chord * 0.20 + perpendicular * (side * 0.28);
        }

        OpenPolyline closure;
        for (size_t segment_idx = 0; segment_idx <= looped_lightning_curve_segments; ++segment_idx)
        {
            const double t = static_cast<double>(segment_idx) / static_cast<double>(looped_lightning_curve_segments);
            const double u = 1.0 - t;
            const Point2LL curve_point
                = source * (u * u * u)
                + control1 * (3.0 * u * u * t)
                + control2 * (3.0 * u * t * t)
                + target * (t * t * t);
            closure.push_back(curve_point);
        }

        OpenLinesSet clipped = limit_to_outline.intersection(OpenLinesSet{ closure });
        result_lines.push_back(std::move(clipped));
    }
}
} // namespace

// Returns 'added someting'.
OpenLinesSet LightningLayer::convertToLines(const Shape& limit_to_outline, const coord_t line_width) const
{
    OpenLinesSet result_lines;
    if (tree_roots.empty())
    {
        return result_lines;
    }

    for (const LightningTreeNodeSPtr& tree : tree_roots)
    {
        tree->convertToPolylines(result_lines, line_width);
    }
    result_lines = limit_to_outline.intersection(result_lines);
    addLoopedLightningClosures(result_lines, limit_to_outline, line_width);

    return result_lines;
}
