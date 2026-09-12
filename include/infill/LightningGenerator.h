// Copyright (c) 2021 Ultimaker B.V.
// CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef LIGHTNING_GENERATOR_H
#define LIGHTNING_GENERATOR_H

#include <functional>
#include <memory>
#include <vector>

#include "../utils/polygonUtils.h"
#include "LightningLayer.h"

namespace cura
{
class SliceMeshStorage;
class SupportStorage;

/*!
 * Generates the Lightning Infill pattern.
 *
 * The lightning infill pattern is designed to use a minimal amount of material
 * to support the top skin of the print, while still printing with reasonably
 * consistently flowing lines. It sacrifices strength completely in favour of
 * top surface quality and reduced print time / material usage.
 *
 * Lightning Infill is so named because the patterns it creates resemble a
 * forked path with one main path and many small lines on the side. These paths
 * grow out from the sides of the model just below where the top surface needs
 * to be supported from the inside, so that minimal material is needed.
 *
 * This pattern is based on a paper called "Ribbed Support Vaults for 3D
 * Printing of Hollowed Objects" by Tricard, Claux and Lefebvre:
 * https://www.researchgate.net/publication/333808588_Ribbed_Support_Vaults_for_3D_Printing_of_Hollowed_Objects
 */
class LightningGenerator // "Just like Nicola used to make!"
{
public:
    LightningGenerator(const SliceMeshStorage& mesh);
    LightningGenerator(const SupportStorage& support);

    const LightningLayer& getTreesForLayer(const size_t& layer_id) const;

protected:
    void generate(
        const coord_t layer_thickness,
        const coord_t line_width,
        const coord_t wall_thickness,
        const coord_t line_distance,
        const AngleRadians& overhang_angle,
        const AngleRadians& prune_angle,
        const AngleRadians& straightening_angle,
        const std::vector<Shape>& areas_per_layer);

    void generateInitialInternalOverhangs(const coord_t infill_wall_thickness, const std::vector<Shape>& areas_per_layer);
    void generateTrees(const coord_t infill_wall_thickness, const std::vector<Shape>& areas_per_layer);

    // CuraLightning post-processing controls for ordinary mesh Lightning.
    double lightning_smoothing{ 0.0 };
    double lightning_offset_widths{ 0.0 };
    bool lightning_single_path{ false };

    coord_t supporting_radius;
    coord_t wall_supporting_radius;
    coord_t prune_length;
    coord_t straightening_max_distance;
    std::vector<Shape> overhang_per_layer;
    std::vector<LightningLayer> lightning_layers;
};

} // namespace cura

#endif // LIGHTNING_GENERATOR_H
