// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/particles/generators.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  namespace Generators
  {
    namespace
    {
      /**
       * Exception
       */
      template <int dim>
      DeclException1(
        ProbabilityFunctionNegative,
        Point<dim>,
        << "Your probability density function in the particle generator "
           "returned a negative value for the following position: "
        << arg1 << ". Please check your function expression.");


      // This function integrates the given probability density
      // function over all cells of the triangulation. For each
      // cell it stores the cumulative sum (of all previous
      // cells including the current cell) in a vector that is then
      // returned. Therefore the returned vector has as many entries
      // as active cells, the first entry being the integral over the
      // first cell, and the last entry the integral over the whole
      // locally owned domain. Cells that are not locally owned
      // simply store the same value as the cell before (equivalent
      // to assuming a probability density function value of 0).
      template <int dim, int spacedim = dim>
      std::vector<double>
      compute_local_cumulative_cell_weights(
        const Triangulation<dim, spacedim> &triangulation,
        const Mapping<dim, spacedim> &      mapping,
        const Function<spacedim> &          probability_density_function)
      {
        std::vector<double> cumulative_cell_weights(
          triangulation.n_active_cells());
        double cumulative_weight = 0.0;

        // Evaluate function at all cell midpoints
        const QMidpoint<dim> quadrature_formula;

        // In the simplest case we do not even need a FEValues object, because
        // using cell->center() and cell->measure() would be equivalent. This
        // fails however for higher-order mappings.
        FE_Nothing<dim, spacedim> alibi_finite_element;
        FEValues<dim, spacedim>   fe_values(mapping,
                                          alibi_finite_element,
                                          quadrature_formula,
                                          update_quadrature_points |
                                            update_JxW_values);

        // compute the integral weight
        for (const auto &cell : triangulation.active_cell_iterators())
          {
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                const double quadrature_point_weight =
                  probability_density_function.value(
                    fe_values.get_quadrature_points()[0]);

                AssertThrow(quadrature_point_weight >= 0.0,
                            ProbabilityFunctionNegative<dim>(
                              quadrature_formula.point(0)));

                // get_cell_weight makes sure to return positive values
                cumulative_weight += quadrature_point_weight * fe_values.JxW(0);
              }
            cumulative_cell_weights[cell->active_cell_index()] =
              cumulative_weight;
          }

        return cumulative_cell_weights;
      }
    } // namespace

    template <int dim, int spacedim>
    void
    regular_reference_locations(
      const Triangulation<dim, spacedim> &triangulation,
      const std::vector<Point<dim>> &     particle_reference_locations,
      ParticleHandler<dim, spacedim> &    particle_handler,
      const Mapping<dim, spacedim> &      mapping)
    {
      types::particle_index particle_index = 0;

#ifdef DEAL_II_WITH_MPI
      if (const auto tria =
            dynamic_cast<const parallel::Triangulation<dim, spacedim> *>(
              &triangulation))
        {
          const types::particle_index n_particles_to_generate =
            tria->n_locally_owned_active_cells() *
            particle_reference_locations.size();

          // The local particle start index is the number of all particles
          // generated on lower MPI ranks.
          MPI_Exscan(&n_particles_to_generate,
                     &particle_index,
                     1,
                     DEAL_II_PARTICLE_INDEX_MPI_TYPE,
                     MPI_SUM,
                     tria->get_communicator());
        }
#endif

      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto &reference_location :
                   particle_reference_locations)
                {
                  const Point<spacedim> position_real =
                    mapping.transform_unit_to_real_cell(cell,
                                                        reference_location);

                  const Particle<dim, spacedim> particle(position_real,
                                                         reference_location,
                                                         particle_index);
                  particle_handler.insert_particle(particle, cell);
                  ++particle_index;
                }
            }
        }

      particle_handler.update_cached_numbers();
    }



    template <int dim, int spacedim>
    Particle<dim, spacedim>
    random_particle_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
      const types::particle_index                                        id,
      std::mt19937 &                random_number_generator,
      const Mapping<dim, spacedim> &mapping)
    {
      // Uniform distribution on the interval [0,1]. This
      // will be used to generate random particle locations.
      std::uniform_real_distribution<double> uniform_distribution_01(0, 1);

      const BoundingBox<spacedim> cell_bounding_box(cell->bounding_box());
      const std::pair<Point<spacedim>, Point<spacedim>> cell_bounds(
        cell_bounding_box.get_boundary_points());

      // Generate random points in these bounds until one is within the cell
      unsigned int       iteration          = 0;
      const unsigned int maximum_iterations = 100;
      Point<spacedim>    particle_position;
      while (iteration < maximum_iterations)
        {
          for (unsigned int d = 0; d < spacedim; ++d)
            {
              particle_position[d] =
                uniform_distribution_01(random_number_generator) *
                  (cell_bounds.second[d] - cell_bounds.first[d]) +
                cell_bounds.first[d];
            }
          try
            {
              const Point<dim> p_unit =
                mapping.transform_real_to_unit_cell(cell, particle_position);
              if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                {
                  // Generate the particle
                  return Particle<dim, spacedim>(particle_position, p_unit, id);
                }
            }
          catch (typename Mapping<dim>::ExcTransformationFailed &)
            {
              // The point is not in this cell. Do nothing, just try again.
            }
          ++iteration;
        }
      AssertThrow(
        iteration < maximum_iterations,
        ExcMessage(
          "Couldn't generate a particle position within the maximum number of tries. "
          "The ratio between the bounding box volume in which the particle is "
          "generated and the actual cell volume is approximately: " +
          std::to_string(
            cell->measure() /
            (cell_bounds.second - cell_bounds.first).norm_square())));

      return Particle<dim, spacedim>();
    }



    template <int dim, int spacedim>
    void
    probabilistic_locations(
      const Triangulation<dim, spacedim> &triangulation,
      const Function<spacedim> &          probability_density_function,
      const bool                          random_cell_selection,
      const types::particle_index         n_particles_to_create,
      ParticleHandler<dim, spacedim> &    particle_handler,
      const Mapping<dim, spacedim> &      mapping,
      const unsigned int                  random_number_seed)
    {
      unsigned int combined_seed = random_number_seed;
      if (const auto tria =
            dynamic_cast<const parallel::Triangulation<dim, spacedim> *>(
              &triangulation))
        {
          const unsigned int my_rank =
            Utilities::MPI::this_mpi_process(tria->get_communicator());
          combined_seed += my_rank;
        }
      std::mt19937 random_number_generator(combined_seed);

      types::particle_index start_particle_id(0);
      types::particle_index n_local_particles(0);

      std::vector<types::particle_index> particles_per_cell(
        triangulation.n_active_cells(), 0);

      // First determine how many particles to generate in which cell
      {
        // Get the local accumulated probabilities for every cell
        const std::vector<double> cumulative_cell_weights =
          compute_local_cumulative_cell_weights(triangulation,
                                                mapping,
                                                probability_density_function);

        // Sum the local integrals over all nodes
        double local_weight_integral = cumulative_cell_weights.back();
        double global_weight_integral;

        if (const auto tria =
              dynamic_cast<const parallel::Triangulation<dim, spacedim> *>(
                &triangulation))
          {
            global_weight_integral =
              Utilities::MPI::sum(local_weight_integral,
                                  tria->get_communicator());
          }
        else
          {
            global_weight_integral = local_weight_integral;
          }

        AssertThrow(global_weight_integral > std::numeric_limits<double>::min(),
                    ExcMessage(
                      "The integral of the user prescribed probability "
                      "density function over the domain equals zero, "
                      "deal.II has no way to determine the cell of "
                      "generated particles. Please ensure that the "
                      "provided function is positive in at least a "
                      "part of the domain; also check the syntax of "
                      "the function."));

        // Determine the starting weight of this process, which is the sum of
        // the weights of all processes with a lower rank
        double local_start_weight = 0.0;

#ifdef DEAL_II_WITH_MPI
        if (const auto tria =
              dynamic_cast<const parallel::Triangulation<dim, spacedim> *>(
                &triangulation))
          {
            MPI_Exscan(&local_weight_integral,
                       &local_start_weight,
                       1,
                       MPI_DOUBLE,
                       MPI_SUM,
                       tria->get_communicator());
          }
#endif

        // Calculate start id
        start_particle_id =
          std::llround(static_cast<double>(n_particles_to_create) *
                       local_start_weight / global_weight_integral);

        // Calcualate number of local particles
        const types::particle_index end_particle_id =
          std::llround(static_cast<double>(n_particles_to_create) *
                       ((local_start_weight + local_weight_integral) /
                        global_weight_integral));
        n_local_particles = end_particle_id - start_particle_id;

        if (random_cell_selection)
          {
            // Uniform distribution on the interval [0,local_weight_integral).
            // This will be used to randomly select cells for all local
            // particles.
            std::uniform_real_distribution<double> uniform_distribution(
              0.0, local_weight_integral);

            // Loop over all particles to create locally and pick their cells
            for (types::particle_index current_particle_index = 0;
                 current_particle_index < n_local_particles;
                 ++current_particle_index)
              {
                // Draw the random number that determines the cell of the
                // particle
                const double random_weight =
                  uniform_distribution(random_number_generator);

                const std::vector<double>::const_iterator selected_cell =
                  std::lower_bound(cumulative_cell_weights.begin(),
                                   cumulative_cell_weights.end(),
                                   random_weight);
                const unsigned int cell_index =
                  std::distance(cumulative_cell_weights.begin(), selected_cell);

                ++particles_per_cell[cell_index];
              }
          }
        else
          {
            // Compute number of particles per cell according to the ratio
            // between their weight and the local weight integral
            types::particle_index particles_created = 0;

            for (const auto &cell : triangulation.active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  const types::particle_index cumulative_particles_to_create =
                    std::llround(
                      static_cast<double>(n_local_particles) *
                      cumulative_cell_weights[cell->active_cell_index()] /
                      local_weight_integral);

                  // Compute particles for this cell as difference between
                  // number of particles that should be created including this
                  // cell minus the number of particles already created.
                  particles_per_cell[cell->active_cell_index()] =
                    cumulative_particles_to_create - particles_created;
                  particles_created +=
                    particles_per_cell[cell->active_cell_index()];
                }
          }
      }

      // Now generate as many particles per cell as determined above
      {
        unsigned int current_particle_index = start_particle_id;

        std::multimap<
          typename Triangulation<dim, spacedim>::active_cell_iterator,
          Particle<dim, spacedim>>
          particles;

        for (const auto &cell : triangulation.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              for (unsigned int i = 0;
                   i < particles_per_cell[cell->active_cell_index()];
                   ++i)
                {
                  Particle<dim, spacedim> particle =
                    random_particle_in_cell(cell,
                                            current_particle_index,
                                            random_number_generator,
                                            mapping);
                  particles.emplace_hint(particles.end(),
                                         cell,
                                         std::move(particle));
                  ++current_particle_index;
                }
            }

        particle_handler.insert_particles(particles);
      }
    }

    template <int dim, int spacedim>
    void
    dof_support_points(const DoFHandler<dim, spacedim> &dof_handler,
                       const std::vector<std::vector<BoundingBox<spacedim>>>
                         &                             global_bounding_boxes,
                       ParticleHandler<dim, spacedim> &particle_handler,
                       const Mapping<dim, spacedim> &  mapping,
                       const ComponentMask &           components)
    {
      const auto &fe = dof_handler.get_fe();

      // Take care of components
      const ComponentMask mask =
        (components.size() == 0 ? ComponentMask(fe.n_components(), true) :
                                  components);

      std::map<types::global_dof_index, Point<spacedim>> support_points_map;

      DoFTools::map_dofs_to_support_points(mapping,
                                           dof_handler,
                                           support_points_map,
                                           mask);

      // Generate the vector of points from the map
      // Memory is reserved for efficiency reasons
      std::vector<Point<spacedim>> support_points_vec;
      support_points_vec.reserve(support_points_map.size());
      for (auto const &element : support_points_map)
        support_points_vec.push_back(element.second);

      particle_handler.insert_global_particles(support_points_vec,
                                               global_bounding_boxes);
    }


    template <int dim, int spacedim>
    void
    quadrature_points(
      const Triangulation<dim, spacedim> &triangulation,
      const Quadrature<dim> &             quadrature,
      // const std::vector<Point<dim>> &     particle_reference_locations,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &                             global_bounding_boxes,
      ParticleHandler<dim, spacedim> &particle_handler,
      const Mapping<dim, spacedim> &  mapping)
    {
      const std::vector<Point<dim>> &particle_reference_locations =
        quadrature.get_points();
      std::vector<Point<spacedim>> points_to_generate;

      //       Loop through cells and gather gauss points
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (const auto &reference_location :
                   particle_reference_locations)
                {
                  const Point<spacedim> position_real =
                    mapping.transform_unit_to_real_cell(cell,
                                                        reference_location);
                  points_to_generate.push_back(position_real);
                }
            }
        }
      particle_handler.insert_global_particles(points_to_generate,
                                               global_bounding_boxes);
    }
  } // namespace Generators
} // namespace Particles

#include "generators.inst"

DEAL_II_NAMESPACE_CLOSE
