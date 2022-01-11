// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>

#include <deal.II/physics/transformations.h>
#include <deal.II/physics/vector_relations.h>


DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Auxiliary constructs for the pipe junction geometry.
   *
   * Please refer to the in-source documentation of the pipe_junction function
   * below for more information about the individual contents.
   */
  namespace PipeSegment
  {
    /**
     * Selection of pipe segment properties to calculate its height with the
     * function below.
     */
    struct AdditionalData
    {
      double skeleton_length;

      double cosecant_polar;
      double cotangent_polar;
      double cotangent_azimuth_half_right;
      double cotangent_azimuth_half_left;
    };



    /**
     * Calculate the height of a pipe segment, depending on the location in the
     * x-y plane.
     */
    inline double
    compute_z_expansion(const double          x,
                        const double          y,
                        const AdditionalData &data)
    {
      return
        // Scale the unit cylinder to the correct length.
        data.skeleton_length
        // Next, adjust for the polar angle. This part will be zero if all
        // openings and the bifurcation are located on a plane.
        + x * data.cotangent_polar
        // Last, adjust for the azimuth angle.
        - std::abs(y) * data.cosecant_polar *
            ((y > 0) ? data.cotangent_azimuth_half_right :
                       data.cotangent_azimuth_half_left);
    }
  } // namespace PipeSegment
} // namespace



namespace GridGenerator
{
  template <int dim, int spacedim>
  void
  pipe_junction(Triangulation<dim, spacedim> &,
                const std::vector<std::pair<Point<spacedim>, double>> &,
                const std::pair<Point<spacedim>, double> &,
                const double)
  {
    Assert(false, ExcNotImplemented());
  }



  // hide the template specialization from doxygen
#ifndef DOXYGEN

  template <>
  void
  pipe_junction(Triangulation<3, 3> &                           tria,
                const std::vector<std::pair<Point<3>, double>> &openings,
                const std::pair<Point<3>, double> &             bifurcation,
                const double                                    aspect_ratio)
  {
    constexpr unsigned int dim      = 3;
    constexpr unsigned int spacedim = 3;
    using vector                    = Tensor<1, spacedim, double>;

    constexpr unsigned int n_pipes   = 3;
    constexpr double       tolerance = 1.e-12;

#  ifdef DEBUG
    // Verify user input.
    Assert(bifurcation.second > 0,
           ExcMessage("Invalid input: negative radius."));
    Assert(openings.size() == n_pipes,
           ExcMessage("Invalid input: only 3 openings allowed."));
    for (const auto &opening : openings)
      Assert(opening.second > 0, ExcMessage("Invalid input: negative radius."));
#  endif

    // Each pipe segment will be identified by the index of its opening in the
    // parameter array. To determine the next and previous entry in the array
    // for a given index, we create auxiliary functions.
    const auto cyclic = [](const unsigned int i) -> unsigned int {
      constexpr unsigned int n_pipes = 3;
      return (i < (n_pipes - 1)) ? i + 1 : 0;
    };
    const auto anticyclic = [](const unsigned int i) -> unsigned int {
      constexpr unsigned int n_pipes = 3;
      return (i > 0) ? i - 1 : n_pipes - 1;
    };

    // Cartesian base represented by unit vectors.
    constexpr std::array<vector, spacedim> directions = {
      {vector({1., 0., 0.}), vector({0., 1., 0.}), vector({0., 0., 1.})}};

    // The skeleton corresponds to the axis of symmetry in the center of each
    // pipe segment. Each skeleton vector points from the associated opening to
    // the common bifurcation point. For convenience, we also compute length and
    // unit vector of every skeleton vector here.
    std::array<vector, n_pipes> skeleton;
    for (unsigned int p = 0; p < n_pipes; ++p)
      skeleton[p] = bifurcation.first - openings[p].first;

    std::array<double, n_pipes> skeleton_length;
    for (unsigned int p = 0; p < n_pipes; ++p)
      skeleton_length[p] = skeleton[p].norm();

    // In many assertions that come up below, we will verify the integrity of
    // the geometry. For this, we introduce a tolerance length which vectors
    // must exceed to avoid being considered "too short". We relate this length
    // to the longest pipe segment.
    const double tolerance_length =
      tolerance *
      *std::max_element(skeleton_length.begin(), skeleton_length.end());

    std::array<vector, n_pipes> skeleton_unit;
    for (unsigned int p = 0; p < n_pipes; ++p)
      {
        Assert(skeleton_length[p] > tolerance_length,
               ExcMessage("Invalid input: bifurcation matches opening."));
        skeleton_unit[p] = skeleton[p] / skeleton_length[p];
      }

    // To determine the orientation of the pipe segments to each other, we will
    // construct a plane: starting from the bifurcation point, we will move by
    // the magnitude one in each of the skeleton directions and span a plane
    // with the three points we reached.
    //
    // The normal vector of this particular plane then describes the edge at
    // which all pipe segments meet. If we would interpret the bifurcation as a
    // ball joint, the normal vector would correspond to the polar axis of the
    // ball.
    vector normal = cross_product_3d(skeleton_unit[1] - skeleton_unit[0],
                                     skeleton_unit[2] - skeleton_unit[0]);
    Assert(normal.norm() > tolerance_length,
           ExcMessage("Invalid input: all three openings "
                      "are located on one line."));
    normal /= normal.norm();

    // Projections of all skeleton vectors perpendicular to the normal vector,
    // or in other words, onto the plane described above.
    std::array<vector, n_pipes> skeleton_plane;
    for (unsigned int p = 0; p < n_pipes; ++p)
      {
        skeleton_plane[p] = skeleton[p] - (skeleton[p] * normal) * normal;
        Assert(std::abs(skeleton_plane[p] * normal) <
                 tolerance * skeleton_plane[p].norm(),
               ExcInternalError());
        Assert(skeleton_plane[p].norm() > tolerance_length,
               ExcMessage("Invalid input."));
      }

    // Create a hyperball domain in 2D that will act as the reference cross
    // section for each pipe segment.
    Triangulation<dim - 1, spacedim - 1> tria_base;
    GridGenerator::hyper_ball_balanced(tria_base,
                                       /*center=*/Point<spacedim - 1>(),
                                       /*radius=*/1.);

    // Now move on to actually build the pipe junction geometry!
    //
    // For each pipe segment, we create a separate triangulation object which
    // will be merged with the parameter triangulation in the end.
    Assert(tria.n_cells() == 0,
           ExcMessage("The output triangulation object needs to be empty."));
    for (unsigned int p = 0; p < n_pipes; ++p)
      {
        Triangulation<dim, spacedim> pipe;

        //
        // Step 1: create unit cylinder
        //
        // We create a unit cylinder by extrusion from the base cross section.
        // The number of layers depends on the ratio of the length of the
        // skeleton and half the minimal radius in the pipe segment. The latter
        // corresponds to the length in radial direction of the smallest cell in
        // the base cross section. Further, the aspect ratio of the extruded
        // cells can be set individually with a function parameter.
        const unsigned int n_slices =
          1 + static_cast<unsigned int>(std::ceil(
                aspect_ratio * skeleton_length[p] /
                (0.5 * std::min(openings[p].second, bifurcation.second))));
        GridGenerator::extrude_triangulation(tria_base,
                                             n_slices,
                                             /*height*/ 1.,
                                             pipe);

        // Set all material and manifold indicators on the unit cylinder, simply
        // because they are easier to handle in this geometry. We will set
        // boundary indicators at the end of the function. See general
        // documentation of this function.
        for (const auto &cell : pipe.active_cell_iterators())
          {
            cell->set_material_id(p);

            for (const auto &face : cell->face_iterators())
              if (face->at_boundary())
                {
                  const auto center_z = face->center()[2];

                  if (std::abs(center_z) < tolerance)
                    {
                      // opening cross section
                    }
                  else if (std::abs(center_z - 1.) < tolerance)
                    {
                      // bifurcation cross section
                    }
                  else
                    {
                      // cone mantle
                      face->set_all_manifold_ids(p);
                    }
                }
          }

        //
        // Step 2: transform unit cylinder to pipe segment
        //
        // For the given cylinder, we will interpret the base in the xy-plane as
        // the cross section of the opening, and the base at z=1 as the surface
        // where all pipe segments meet. On the latter surface, we assign the
        // section in positive y-direction to face the next (right/cyclic) pipe
        // segment, and allocate the domain in negative y-direction to border
        // the previous (left/anticyclic) pipe segment.
        //
        // In the end, the transformed pipe segment will look like this:
        //              z                   z
        //              ^                   ^
        //         left | right             |  /|
        //   anticyclic | cyclic            |/  |
        //             /|\                 /|   |
        //           /  |  \             /  |   |
        //          |   |   |           |   |   |
        //          |   |   |           |   |   |
        //        ------+----->y      ------+----->x

        // Before transforming the unit cylinder however, we compute angle
        // relations between the skeleton vectors viewed from the bifurcation
        // point. For this purpose, we interpret the bifurcation as a ball joint
        // as described above.
        //
        // In spherical coordinates, the polar angle describes the kink of the
        // skeleton vector with respect to the polar axis. If all openings and
        // the bifurcation are located on a plane, then this angle is pi/2 for
        // every pipe segment.
        const double polar_angle =
          Physics::VectorRelations::angle(skeleton[p], normal);
        Assert(std::abs(polar_angle) > tolerance &&
                 std::abs(polar_angle - numbers::PI) > tolerance,
               ExcMessage("Invalid input."));

        // Further, we compute the angles between this pipe segment to the other
        // two. The angle corresponds to the azimuthal direction if we stick to
        // the picture of the ball joint.
        const double azimuth_angle_right =
          Physics::VectorRelations::signed_angle(skeleton_plane[p],
                                                 skeleton_plane[cyclic(p)],
                                                 /*axis=*/normal);
        Assert(std::abs(azimuth_angle_right) > tolerance,
               ExcMessage("Invalid input: at least two openings located "
                          "in same direction from bifurcation"));

        const double azimuth_angle_left =
          Physics::VectorRelations::signed_angle(skeleton_plane[p],
                                                 skeleton_plane[anticyclic(p)],
                                                 /*axis=*/-normal);
        Assert(std::abs(azimuth_angle_left) > tolerance,
               ExcMessage("Invalid input: at least two openings located "
                          "in same direction from bifurcation"));

        // We compute some trigonometric relations with these angles, and store
        // them conveniently in a struct to be reused later.
        PipeSegment::AdditionalData data;
        data.skeleton_length = skeleton_length[p];
        data.cosecant_polar  = 1. / std::sin(polar_angle);
        data.cotangent_polar = std::cos(polar_angle) * data.cosecant_polar;
        data.cotangent_azimuth_half_right = std::cos(.5 * azimuth_angle_right) /
                                            std::sin(.5 * azimuth_angle_right);
        data.cotangent_azimuth_half_left =
          std::cos(.5 * azimuth_angle_left) / std::sin(.5 * azimuth_angle_left);

        // Now transform the cylinder as described above.
        const auto pipe_segment = [&](const Point<spacedim> &pt) {
          // We transform the cylinder in x- and y-direction to become a
          // truncated cone, similarly to GridGenerator::truncated_cone().
          const double r_factor =
            (bifurcation.second - openings[p].second) * pt[2] +
            openings[p].second;
          const double x_new = r_factor * pt[0];
          const double y_new = r_factor * pt[1];

          // Further, to be able to smoothly merge all pipe segments at the
          // bifurcation, we also need to transform in z-direction.
          const double z_factor =
            PipeSegment::compute_z_expansion(x_new, y_new, data);
          Assert(z_factor > 0,
                 ExcMessage("Invalid input: at least one pipe segment "
                            "is not long enough in this configuration"));
          const double z_new = z_factor * pt[2];

          // TODO: MSVC can't capture const or constexpr values in lambda
          // functions (due to either a missing implementation or a bug).
          // Instead, we duplicate the declaration here.
          //   See also: https://developercommunity.visualstudio.com/t/
          //             invalid-template-argument-expected-compile-time-co/187862
          constexpr unsigned int spacedim = 3;
          return Point<spacedim>(x_new, y_new, z_new);
        };
        GridTools::transform(pipe_segment, pipe);

        //
        // Step 3: rotate pipe segment to match skeleton direction
        //
        // The symmetry axis of the pipe segment in its current state points in
        // positive z-direction. We rotate the pipe segment that its symmetry
        // axis matches the direction of the skeleton vector. For this purpose,
        // we rotate the pipe segment around the axis that is described by the
        // cross product of both vectors.
        const double rotation_angle =
          Physics::VectorRelations::angle(directions[2], skeleton_unit[p]);
        const vector rotation_axis = [&]() {
          const vector rotation_axis =
            cross_product_3d(directions[2], skeleton_unit[p]);
          const double norm = rotation_axis.norm();
          if (norm < tolerance)
            return directions[1];
          else
            return rotation_axis / norm;
        }();
        const Tensor<2, spacedim, double> rotation_matrix =
          Physics::Transformations::Rotations::rotation_matrix_3d(
            rotation_axis, rotation_angle);
        GridTools::transform(
          [&](const Point<spacedim> &pt) { return rotation_matrix * pt; },
          pipe);

        //
        // Step 4: rotate laterally to align pipe segments
        //
        // On the unit cylinder, we find that the edge on which all pipe
        // segments meet is parallel to the x-axis. After the transformation to
        // the pipe segment, we notice that this statement still holds for the
        // projection of this edge onto the xy-plane, which corresponds to the
        // cross section of the opening.
        //
        // With the latest rotation however, this is no longer the case. We
        // rotate the unit vector in x-direction in the same fashion, which
        // gives us the current direction of the projected edge.
        const vector Rx = rotation_matrix * directions[0];

        // To determine how far we need to rotate, we also need to project the
        // polar axis of the bifurcation ball joint into the same plane.
        const vector normal_projected_on_opening =
          normal - (normal * skeleton_unit[p]) * skeleton_unit[p];

        // Both the projected normal and Rx must be in the opening plane.
        Assert(std::abs(skeleton_unit[p] * normal_projected_on_opening) <
                 tolerance,
               ExcInternalError());
        Assert(std::abs(skeleton_unit[p] * Rx) < tolerance, ExcInternalError());

        // Now we laterally rotate the pipe segment around its own symmetry axis
        // that the edge matches the polar axis.
        const double lateral_angle =
          Physics::VectorRelations::signed_angle(Rx,
                                                 normal_projected_on_opening,
                                                 /*axis=*/skeleton_unit[p]);
        GridTools::rotate(skeleton_unit[p], lateral_angle, pipe);

        //
        // Step 5: shift to final position
        //
        GridTools::shift(openings[p].first, pipe);

        GridGenerator::merge_triangulations(
          pipe, tria, tria, tolerance_length, /*copy_manifold_ids=*/true);
      }

    // Since GridGenerator::merge_triangulations() does not copy boundary IDs,
    // we need to set them after the final geometry is created. Luckily,
    // boundary IDs match with material IDs, so we simply translate them with
    // the help of manifold IDs to identify openings.
    for (const auto &cell : tria.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
        if (face->at_boundary())
          {
            if (face->manifold_id() == numbers::flat_manifold_id)
              // opening cross section
              face->set_boundary_id(cell->material_id());
            else
              // cone mantle
              face->set_boundary_id(n_pipes);
          }
  }

#endif

} // namespace GridGenerator


// explicit instantiations
#include "grid_generator_pipe_junction.inst"

DEAL_II_NAMESPACE_CLOSE
