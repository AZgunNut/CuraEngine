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

namespace
{
constexpr coord_t looped_lightning_max_distance = 12000; // 12 mm in CuraEngine's micron coordinate system.
constexpr size_t looped_lightning_curve_segments = 6;

/*!
 * Close dangling Lightning polylines toward nearby infill or the island boundary.
 *
 * LightningTreeNode::convertToPolylines() emits polylines beginning at a leaf.
 * For each leaf we look for the nearest point on another Lightning polyline and
 * compare that with the nearest model boundary. If either is within 12 mm, add
 * a short quadratic Bezier connector. The initial tangent continues away from
 * the existing branch so the closure forms a flowing hook rather than a hard
 * V-shaped reversal.
 *
 * This is deliberately a small proof of concept. Later versions can replace
 * the vertex-only target search with nearest-point-on-segment indexing and make
 * the distance/curvature user settings.
 */
void addLoopedLightningClosures(OpenLinesSet& result_lines, const Shape& limit_to_outline, const coord_t line_width)
{
    if (result_lines.empty() || limit_to_outline.empty())
    {
        return;
    }

    const size_t original_line_count = result_lines.size();
    OpenLinesSet closures;

    for (size_t line_idx = 0; line_idx < original_line_count; ++line_idx)
    {
        const OpenPolyline& source_line = result_lines[line_idx];
        if (source_line.size() < 2)
        {
            continue;
        }

        const Point2LL source = source_line.front(); // Lightning polylines begin at leaves.
        Point2LL target = PolygonUtils::findClosest(source, limit_to_outline).p();
        coord_t best_distance = vSize(target - source);

        // Prefer another piece of Lightning if it is nearer than the wall.
        for (size_t candidate_line_idx = 0; candidate_line_idx < original_line_count; ++candidate_line_idx)
        {
            if (candidate_line_idx == line_idx)
            {
                continue;
            }

            const OpenPolyline& candidate_line = result_lines[candidate_line_idx];
            for (const Point2LL& candidate : candidate_line)
            {
                const coord_t candidate_distance = vSize(candidate - source);
                if (candidate_distance < best_distance)
                {
                    best_distance = candidate_distance;
                    target = candidate;
                }
            }
        }

        // Ignore long closures and tiny hooks that would only over-extrude a junction.
        if (best_distance > looped_lightning_max_distance || best_distance < line_width * 2)
        {
            continue;
        }

        const Point2LL branch_outward = source - source_line[1];
        const coord_t branch_length = vSize(branch_outward);

        Point2LL control = (source + target) / 2;
        if (branch_length > 0)
        {
            // Continue the leaf tangent for roughly one third of the closure length.
            const double tangent_scale = 0.35 * static_cast<double>(best_distance) / static_cast<double>(branch_length);
            control = source + branch_outward * tangent_scale;
        }

        OpenPolyline closure;
        for (size_t segment_idx = 0; segment_idx <= looped_lightning_curve_segments; ++segment_idx)
        {
            const double t = static_cast<double>(segment_idx) / static_cast<double>(looped_lightning_curve_segments);
            const double one_minus_t = 1.0 - t;
            const Point2LL curve_point
                = source * (one_minus_t * one_minus_t) + control * (2.0 * one_minus_t * t) + target * (t * t);
            closure.push_back(curve_point);
        }
        closures.push_back(std::move(closure), CheckNonEmptyParam::OnlyIfValid);
    }

    result_lines.push_back(std::move(closures));
    // Curves can bulge out of concave islands. Clip everything back to the valid infill region.
    result_lines = limit_to_outline.intersection(result_lines);
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
