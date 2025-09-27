// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

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



    /**
     * Pipe segment manifold description.
     *
     * The manifold class is too specific to being of any other use than for the
     * pipe junction geometry.
     */
    template <int dim, int spacedim = dim>
    class Manifold : public ChartManifold<dim, spacedim, 3>
    {
    public:
      /**
       * Constructor. The manifold described is a pipe segment whose central
       * axis points in @p direction and goes through the given
       * @p point_on_axis.
       *
       * The vector @p normal_direction needs to be perpendicular to
       * @p direction and acts as the reference for azimuth angles. In other
       * words, $\phi = 0$ corresponds to @p normal_direction. For the pipe
       * junction geometry, it needs to point in the same direction as the
       * projection of the bifurcation edge on the opening plane, as this is
       * the reference direction for azimuth angles by construction.
       *
       * Both @p normal_direction and @p direction have to be unit vectors.
       *
       * We need @p data to map the pipe segment to a regular cylinder and back.
       *
       * The @p tolerance value is used to validate both direction vectors and
       * to determine if a point in on the axis.
       */
      Manifold(const Tensor<1, spacedim> &normal_direction,
               const Tensor<1, spacedim> &direction,
               const Point<spacedim>     &point_on_axis,
               const AdditionalData      &data,
               const double               tolerance = 1e-10);

      /**
       * Make a clone of this Manifold object.
       */
      virtual std::unique_ptr<dealii::Manifold<dim, spacedim>>
      clone() const override;

      /**
       * Compute the cylindrical coordinates $(r, \phi, \lambda)$ for the given
       * space point and map them to a cylinder of height one, where $r$ denotes
       * the distance from the axis, $\phi$ the angle between the given point
       * and the normal direction, and $\lambda$ the axial position.
       */
      virtual Point<3>
      pull_back(const Point<spacedim> &space_point) const override;

      /**
       * Compute the Cartesian coordinates for a chart point given in
       * cylindrical coordinates $(r, \phi, \lambda)$ on a cylinder of height
       * one, where $r$ denotes the distance from the axis, $\phi$ the angle
       * between the given point and the normal direction, and $\lambda$ the
       * axial position.
       */
      virtual Point<spacedim>
      push_forward(const Point<3> &chart_point) const override;

    private:
      /**
       * A vector orthogonal to the normal direction.
       */
      const Tensor<1, spacedim> normal_direction;

      /**
       * The direction vector of the axis.
       */
      const Tensor<1, spacedim> direction;

      /**
       * An arbitrary point on the axis.
       */
      const Point<spacedim> point_on_axis;

      /**
       * Pipe segment properties to calculate its height.
       */
      const AdditionalData data;

      /**
       * Relative tolerance to measure zero distances.
       */
      const double tolerance;

      /**
       * The direction vector perpendicular to both direction and
       * normal_direction.
       */
      const Tensor<1, spacedim> dxn;
    };



    template <int dim, int spacedim>
    Manifold<dim, spacedim>::Manifold(
      const Tensor<1, spacedim> &normal_direction,
      const Tensor<1, spacedim> &direction,
      const Point<spacedim>     &point_on_axis,
      const AdditionalData      &data,
      const double               tolerance)
      : ChartManifold<dim, spacedim, 3>(Tensor<1, 3>({0, 2. * numbers::PI, 0}))
      , normal_direction(normal_direction)
      , direction(direction)
      , point_on_axis(point_on_axis)
      , data(data)
      , tolerance(tolerance)
      , dxn(cross_product_3d(direction, normal_direction))
    {
      Assert(spacedim == 3,
             ExcMessage(
               "PipeSegment::Manifold can only be used for spacedim==3!"));

      Assert(std::abs(normal_direction.norm() - 1) < tolerance,
             ExcMessage("Normal direction must be unit vector."));
      Assert(std::abs(direction.norm() - 1) < tolerance,
             ExcMessage("Direction must be unit vector."));
      Assert(normal_direction * direction < tolerance,
             ExcMessage(
               "Direction and normal direction must be perpendicular."));
    }



    template <int dim, int spacedim>
    std::unique_ptr<dealii::Manifold<dim, spacedim>>
    Manifold<dim, spacedim>::clone() const
    {
      return std::make_unique<Manifold<dim, spacedim>>(*this);
    }



    template <int dim, int spacedim>
    Point<3>
    Manifold<dim, spacedim>::pull_back(const Point<spacedim> &space_point) const
    {
      // First find the projection of the given point to the axis.
      const Tensor<1, spacedim> normalized_point = space_point - point_on_axis;
      double                    lambda           = normalized_point * direction;
      const Point<spacedim>     projection = point_on_axis + direction * lambda;
      const Tensor<1, spacedim> p_diff     = space_point - projection;
      const double              r          = p_diff.norm();

      Assert(r > tolerance * data.skeleton_length,
             ExcMessage(
               "This class won't handle points on the direction axis."));

      // Then compute the angle between the projection direction and
      // another vector orthogonal to the direction vector.
      const double phi =
        Physics::VectorRelations::signed_angle(normal_direction,
                                               p_diff,
                                               /*axis=*/direction);

      // Map the axial coordinate to a cylinder of height one.
      lambda /= compute_z_expansion(r * std::cos(phi), r * std::sin(phi), data);

      // Return distance from the axis, angle and signed distance on the axis.
      return {r, phi, lambda};
    }



    template <int dim, int spacedim>
    Point<spacedim>
    Manifold<dim, spacedim>::push_forward(const Point<3> &chart_point) const
    {
      // Rotate the orthogonal direction by the given angle.
      const double sine_r   = chart_point[0] * std::sin(chart_point[1]);
      const double cosine_r = chart_point[0] * std::cos(chart_point[1]);

      const Tensor<1, spacedim> intermediate =
        normal_direction * cosine_r + dxn * sine_r;

      // Map the axial coordinate back to the pipe segment.
      const double lambda =
        chart_point[2] * compute_z_expansion(cosine_r, sine_r, data);

      // Finally, put everything together.
      return point_on_axis + direction * lambda + intermediate;
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
    DEAL_II_NOT_IMPLEMENTED();
  }



  // hide the template specialization from doxygen
#ifndef DOXYGEN

  template <>
  void
  pipe_junction(Triangulation<3, 3>                            &tria,
                const std::vector<std::pair<Point<3>, double>> &openings,
                const std::pair<Point<3>, double>              &bifurcation,
                const double                                    aspect_ratio)
  {
    constexpr unsigned int dim      = 3;
    constexpr unsigned int spacedim = 3;
    using vector3d                  = Tensor<1, spacedim, double>;

    constexpr unsigned int n_pipes   = 3;
    constexpr double       tolerance = 1.e-12;

    if constexpr (running_in_debug_mode())
      {
        // Verify user input.
        Assert(bifurcation.second > 0,
               ExcMessage("Invalid input: negative radius."));
        Assert(openings.size() == n_pipes,
               ExcMessage("Invalid input: only 3 openings allowed."));
        for (const auto &opening : openings)
          Assert(opening.second > 0,
                 ExcMessage("Invalid input: negative radius."));
      }

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
    const std::array<vector3d, spacedim> directions = {
      {vector3d({1., 0., 0.}), vector3d({0., 1., 0.}), vector3d({0., 0., 1.})}};

    // The skeleton corresponds to the axis of symmetry in the center of each
    // pipe segment. Each skeleton vector points from the associated opening to
    // the common bifurcation point. For convenience, we also compute length and
    // unit vector of every skeleton vector here.
    std::array<vector3d, n_pipes> skeleton;
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

    std::array<vector3d, n_pipes> skeleton_unit;
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
    vector3d normal = cross_product_3d(skeleton_unit[1] - skeleton_unit[0],
                                       skeleton_unit[2] - skeleton_unit[0]);
    Assert(normal.norm() > tolerance_length,
           ExcMessage("Invalid input: all three openings "
                      "are located on one line."));
    normal /= normal.norm();

    // Projections of all skeleton vectors perpendicular to the normal vector,
    // or in other words, onto the plane described above.
    std::array<vector3d, n_pipes> skeleton_plane;
    for (unsigned int p = 0; p < n_pipes; ++p)
      {
        skeleton_plane[p] = skeleton[p] - (skeleton[p] * normal) * normal;
        Assert(std::abs(skeleton_plane[p] * normal) <
                 tolerance * skeleton_plane[p].norm(),
               ExcInternalError());
        Assert(skeleton_plane[p].norm() > tolerance_length,
               ExcMessage("Invalid input."));
      }

    // Create a hyperball domain in 2d that will act as the reference cross
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

    std::vector<PipeSegment::Manifold<dim, spacedim>> manifolds;
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
                      cell->set_all_manifold_ids(n_pipes);
                      face->set_all_manifold_ids(p);
                    }
                }
          }

        //
        // Step 2: transform unit cylinder to pipe segment
        //
        // For the given cylinder, we will interpret the base in the xy-plane as
        // the cross section of the opening, and the top at z=1 as the surface
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
        const auto pipe_segment_transform =
          [&](const Point<spacedim> &pt) -> Point<spacedim> {
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

          return {x_new, y_new, z_new};
        };
        GridTools::transform(pipe_segment_transform, pipe);

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
        const vector3d rotation_axis = [&]() {
          const vector3d rotation_axis =
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
        const vector3d Rx = rotation_matrix * directions[0];

        // To determine how far we need to rotate, we also need to project the
        // polar axis of the bifurcation ball joint into the same plane.
        const vector3d normal_projected_on_opening =
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
        // Step 5: shift to final position and merge this pipe into the entire
        // assembly
        //
        GridTools::shift(openings[p].first, pipe);

        // Create a manifold object for the mantle of this particular pipe
        // segment. Since GridGenerator::merge_triangulations() does not copy
        // manifold objects, but just IDs if requested, we will copy them to
        // the final triangulation later.
        manifolds.emplace_back(
          /*normal_direction=*/normal_projected_on_opening /
            normal_projected_on_opening.norm(),
          /*direction=*/skeleton_unit[p],
          /*point_on_axis=*/openings[p].first,
          data,
          tolerance);

        GridGenerator::merge_triangulations(
          pipe, tria, tria, tolerance_length, /*copy_manifold_ids=*/true);
      }

    for (unsigned int p = 0; p < n_pipes; ++p)
      tria.set_manifold(p, manifolds[p]);
    tria.set_manifold(n_pipes, FlatManifold<3>());

    // Since GridGenerator::merge_triangulations() does not copy boundary IDs
    // either, we need to set them after the final geometry is created. Luckily,
    // boundary IDs match with material IDs, so we simply translate them with
    // the help of manifold IDs to identify openings.
    for (const auto &cell : tria.active_cell_iterators())
      for (const auto &face : cell->face_iterators())
        if (face->at_boundary())
          {
            if (face->manifold_id() == numbers::flat_manifold_id ||
                face->manifold_id() == n_pipes)
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
#include "grid/grid_generator_pipe_junction.inst"

DEAL_II_NAMESPACE_CLOSE
