// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_non_matching_dof_handler_coupling
#define dealii_non_matching_dof_handler_coupling

#include <deal.II/base/config.h>

#include <deal.II/base/iterator_range.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>


DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  /**
   * Handle coupling between two DoFHandler objects built on non-matching
   * Triangulation objects.
   *
   * @param cache1
   * @param cache2
   * @param dh1
   * @param dh2
   * @param comps1
   * @param comps2
   * @param mapping1
   * @param mapping2
   *
   * @author Luca Heltai, Bruno Blais, 2019
   */
  template <int dim, int spacedim>
  class DoFHandlerCoupling : public Subscriptor
  {
  public:
    DoFHandlerCoupling(
      const GridTools::Cache<spacedim> &     cache1,
      const GridTools::Cache<dim, spacedim> &cache2,
      const DoFHandler<spacedim> &           dh1,
      const DoFHandler<dim, spacedim> &      dh2,
      const ComponentMask &                  comps1_ = ComponentMask(),
      const ComponentMask &                  comps2_ = ComponentMask(),
      const Mapping<spacedim> &mapping1 = StaticMappingQ1<spacedim>::mapping,
      const Mapping<dim, spacedim> &mapping2 =
        StaticMappingQ1<dim, spacedim>::mapping);



    using CellsAndQuadratureIterator = typename std::vector<
      std::pair<typename DoFHandler<spacedim>::cell_iterator,
                ImmersedSurfaceQuadrature<spacedim>>>::iterator;

    /**
     * @brief DoFHandlerCoupling
     * @param dh1
     * @param dh2
     * @param comps1
     * @param comps2
     * @param mapping1
     * @param mapping2
     */
    DoFHandlerCoupling(
      const DoFHandler<spacedim> &     dh1,
      const DoFHandler<dim, spacedim> &dh2,
      const ComponentMask &            comps1 = ComponentMask(),
      const ComponentMask &            comps2 = ComponentMask(),
      const Mapping<spacedim> &mapping1 = StaticMappingQ1<spacedim>::mapping,
      const Mapping<dim, spacedim> &mapping2 =
        StaticMappingQ1<dim, spacedim>::mapping);



    /**
     * Assemble sparsity.
     *
     * @param[in] quad
     * @param[out] sparsity
     * @param[in] constraints
     * @param[in] reuse_internal_data_structures
     */
    template <class Sparsity, typename number = double>
    void
    create_interpolation_sparsity_pattern(
      Sparsity &                       sparsity,
      const AffineConstraints<number> &constraints =
        AffineConstraints<number>(),
      const bool reuse_internal_data_structures = false) const
    {
      possibly_generate_particle_handler(reuse_internal_data_structures);

      std::vector<types::global_dof_index> dof_indices1(fe1->dofs_per_cell);

      auto       particle = particle_handler->begin();
      const auto max_particles_per_cell =
        particle_handler->n_global_max_particles_per_cell();

      while (particle != particle_handler->end())
        {
          const auto &cell = particle->get_surrounding_cell(*tria1);
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, dh1);
          dh_cell->get_dof_indices(dof_indices1);

          const auto pic         = particle_handler->particles_in_cell(cell);
          const auto n_particles = particle_handler->n_particles_in_cell(cell);

          // local_matrix.reinit({n_particles, fe.dofs_per_cell});
          Assert(pic.begin() == particle, ExcInternalError());
          for (unsigned int i = 0; particle != pic.end(); ++particle, ++i)
            {
              const auto properties = particle->get_properties();
              for (unsigned int j = 0; j < fe1->dofs_per_cell; ++j)
                {
                  const auto comp_j =
                    gtl1[fe1->system_to_component_index(j).first];
                  if (comp_j != numbers::invalid_unsigned_int)
                    {
                      const auto cj = static_cast<types::global_dof_index>(
                        properties[1 + comp_j]);
                      constraints.add_entries_local_to_global({cj},
                                                              {dof_indices1[j]},
                                                              sparsity);
                    }
                }
              // [TODO]: when this works, use this:
              // constraints.add_entries_local_to_global(particle_indices,
              //                                         dof_indices,
              //                                         sparsity,
              //                                         dof_mask);
            }
        }
    }

    /**
     * Assemble sparsity for Nitsche method.
     *
     * @param[in] quadrature
     * @param[out] sparsity
     * @param[in] constraints
     * @param[in] reuse_internal_data_structures
     */
    template <class MatrixType, class VectorType, typename number = double>
    void
    create_nitsche_restriction(
      const Quadrature<dim> &          quadrature,
      const Function<spacedim> &       rhs_function,
      const double &                   penalty_term,
      MatrixType &                     matrix,
      VectorType &                     rhs,
      const AffineConstraints<number> &constraints =
        AffineConstraints<number>(),
      const bool reuse_internal_data_structures = false) const
    {
      if (penalty_term == 0.0)
        return;
      possibly_generate_particle_handler(reuse_internal_data_structures,
                                         quadrature);

      std::vector<types::global_dof_index> dof_indices1(fe1->dofs_per_cell);


      FullMatrix<double> local_matrix(fe1->dofs_per_cell, fe1->dofs_per_cell);
      dealii::Vector<double> local_rhs(fe1->dofs_per_cell);

      auto particle = quadrature_particle_handler->begin();
      while (particle != quadrature_particle_handler->end())
        {
          local_matrix     = 0;
          local_rhs        = 0;
          const auto &cell = particle->get_surrounding_cell(*tria1);
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, dh1);
          dh_cell->get_dof_indices(dof_indices1);

          const auto pic = quadrature_particle_handler->particles_in_cell(cell);

          for (const auto &p : pic)
            {
              const auto  ref_q      = p.get_reference_location();
              const auto  real_q     = p.get_location();
              const auto  properties = p.get_properties();
              const auto &JxW        = properties[0];
              for (unsigned int i = 0; i < fe1->dofs_per_cell; ++i)
                {
                  const auto comp_i = fe1->system_to_component_index(i).first;
                  if (comps1[comp_i])
                    {
                      for (unsigned int j = 0; j < fe1->dofs_per_cell; ++j)
                        {
                          const auto comp_j =
                            fe1->system_to_component_index(j).first;
                          if (comp_i == comp_j)
                            local_matrix(i, j) +=
                              penalty_term * fe1->shape_value(i, ref_q) *
                              fe1->shape_value(j, ref_q) * JxW;
                        }
                      local_rhs(i) += penalty_term *
                                      rhs_function.value(real_q, comp_i) *
                                      fe1->shape_value(i, ref_q) * JxW;
                    }
                }
            }
          constraints.distribute_local_to_global(
            local_matrix, local_rhs, dof_indices1, matrix, rhs);
          particle = pic.end();
        }
    }



    template <class Matrix, typename number = double>
    void
    create_interpolation_matrix(Matrix &                         matrix,
                                const AffineConstraints<number> &constraints =
                                  AffineConstraints<number>()) const
    {
      if (particle_handler->n_locally_owned_particles() == 0)
        {
          matrix.compress(VectorOperation::add);
          return; // nothing else to do here
        }

      AssertDimension(matrix.n(), dh1->n_dofs());
      AssertDimension(matrix.m(), dh2->n_dofs());

      auto       particle = particle_handler->begin();
      const auto max_particles_per_cell =
        particle_handler->n_global_max_particles_per_cell();

      std::vector<types::global_dof_index> dof_indices1(fe1->dofs_per_cell);
      std::vector<types::global_dof_index> dof_indices2(
        particle_handler->n_global_max_particles_per_cell());

      FullMatrix<double> local_matrix(max_particles_per_cell * n_comps,
                                      fe1->dofs_per_cell);


      while (particle != particle_handler->end())
        {
          const auto &cell = particle->get_surrounding_cell(*tria1);
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, dh1);
          dh_cell->get_dof_indices(dof_indices1);

          const auto pic         = particle_handler->particles_in_cell(cell);
          const auto n_particles = particle_handler->n_particles_in_cell(cell);
          dof_indices2.resize(n_particles * n_comps);
          local_matrix.reinit({n_particles * n_comps, fe1->dofs_per_cell});

          Assert(pic.begin() == particle, ExcInternalError());
          for (unsigned int i = 0; particle != pic.end(); ++particle, ++i)
            {
              const auto &reference_location =
                particle->get_reference_location();

              const auto properties = particle->get_properties();

              for (unsigned int d = 0; d < n_comps; ++d)
                dof_indices2[i * n_comps + d] =
                  static_cast<types::global_dof_index>(properties[1 + d]);

              for (unsigned int j = 0; j < fe1->dofs_per_cell; ++j)
                {
                  const auto comp_j =
                    gtl1[fe1->system_to_component_index(j).first];
                  if (comp_j != numbers::invalid_unsigned_int)
                    {
                      const auto li = i * n_comps + comp_j;
                      local_matrix(li, j) =
                        fe1->shape_value(j, reference_location);
                    }
                }
            }
          constraints.distribute_local_to_global(local_matrix,
                                                 dof_indices2,
                                                 dof_indices1,
                                                 matrix);
        }
      matrix.compress(VectorOperation::add);
    }

    const Particles::ParticleHandler<spacedim> &
    get_particle_handler() const
    {
      Assert(particle_handler,
             ExcMessage("You need to call possibly_generate_particle_handler "
                        "function first."));
      return *particle_handler;
    }

    const Particles::ParticleHandler<spacedim> &
    get_quadrature_particle_handler() const
    {
      Assert(quadrature_particle_handler,
             ExcMessage("You need to call possibly_generate_particle_handler "
                        "function first."));
      return *quadrature_particle_handler;
    }

    void
    possibly_generate_particle_handler(
      bool                   reuse_internal_data_structures,
      const Quadrature<dim> &quadrature = Quadrature<dim>()) const;

    void
    set_quadrature_particles_positions(
      const std::vector<Tensor<1, spacedim>> &positions_vector,
      const bool                              displace_particles = true) const
    {
      // There should be a position per per local particles
      AssertDimension(positions_vector.size(),
                      quadrature_particle_handler->n_local_particles());
      Particles::Utilities::set_particle_positions(positions_vector,
                                                   *quadrature_particle_handler,
                                                   displace_particles);
    }

    /**
     * Displace the particles of the quadrature_particle_handler using a
     * vector of Tensor<1,spacedim> of size
     * quadrature_particle_handler.n_local_particles()
     */
    void
    set_particles_positions(
      const std::vector<Tensor<1, spacedim>> &positions_vector,
      const bool                              displace_particles = true) const
    {
      // There should be a position per per local particles
      AssertDimension(positions_vector.size(),
                      particle_handler->n_local_particles());
      Particles::Utilities::set_particle_positions(positions_vector,
                                                   *particle_handler,
                                                   displace_particles);
    }

    /**
     * Displace the particles of the quadrature particle_handler using a
     * function
     */
    void
    set_quadrature_particles_positions(
      const Function<spacedim> &function,
      const bool                displace_particles = true) const
    {
      Particles::Utilities::set_particle_positions(function,
                                                   *quadrature_particle_handler,
                                                   displace_particles);
    }

    /**
     * Displace the particles of the particle_handler using a
     * function
     */
    void
    set_particles_positions(const Function<spacedim> &function,
                            const bool displace_particles = true) const
    {
      Particles::Utilities::set_particle_positions(function,
                                                   *quadrature_particle_handler,
                                                   displace_particles);
    }


  private:
    SmartPointer<const Mapping<spacedim>>      mapping1;
    SmartPointer<const Mapping<dim, spacedim>> mapping2;

    SmartPointer<const Triangulation<spacedim>>      tria1;
    SmartPointer<const Triangulation<dim, spacedim>> tria2;

    SmartPointer<const DoFHandler<spacedim>>      dh1;
    SmartPointer<const DoFHandler<dim, spacedim>> dh2;

    SmartPointer<const FiniteElement<spacedim>>      fe1;
    SmartPointer<const FiniteElement<dim, spacedim>> fe2;

    std::unique_ptr<const GridTools::Cache<spacedim>>      cache1_ptr;
    std::unique_ptr<const GridTools::Cache<dim, spacedim>> cache2_ptr;

    SmartPointer<const GridTools::Cache<spacedim>>      cache1;
    SmartPointer<const GridTools::Cache<dim, spacedim>> cache2;

    std::vector<std::pair<typename DoFHandler<spacedim>::cell_iterator,
                          ImmersedSurfaceQuadrature<spacedim>>>
      cells_and_quadratures;

    ComponentMask comps1;
    ComponentMask comps2;

    unsigned int n_comps;
    unsigned int mpi_proc;
    unsigned int n_mpi_procs;
    MPI_Comm     communicator;

    std::vector<unsigned int> gtl1;
    std::vector<unsigned int> gtl2;

    std::vector<std::vector<BoundingBox<spacedim>>> global_bounding_boxes;

    mutable std::unique_ptr<Particles::ParticleHandler<spacedim>>
      particle_handler;

    mutable std::unique_ptr<Particles::ParticleHandler<spacedim>>
      quadrature_particle_handler;
  };



  // Implementation
  template <int dim, int spacedim>
  DoFHandlerCoupling<dim, spacedim>::DoFHandlerCoupling(
    const GridTools::Cache<spacedim> &     cache1,
    const GridTools::Cache<dim, spacedim> &cache2,
    const DoFHandler<spacedim> &           dh1,
    const DoFHandler<dim, spacedim> &      dh2,
    const ComponentMask &                  comps1_,
    const ComponentMask &                  comps2_,
    const Mapping<spacedim> &              mapping1,
    const Mapping<dim, spacedim> &         mapping2)
    : mapping1(&mapping1)
    , mapping2(&mapping2)
    , tria1(&dh1.get_triangulation())
    , tria2(&dh2.get_triangulation())
    , dh1(&dh1)
    , dh2(&dh2)
    , fe1(&dh1.get_fe())
    , fe2(&dh2.get_fe())
    , cache1(&cache1)
    , cache2(&cache2)
    , comps1(comps1_ == ComponentMask() ?
               ComponentMask(fe1->n_components(), true) :
               comps1_)
    , comps2(comps2_ == ComponentMask() ?
               ComponentMask(fe2->n_components(), true) :
               comps2_)

  {
    // Check consistency of components
    AssertDimension(comps1.size(), fe1->n_components());
    AssertDimension(comps2.size(), fe2->n_components());


    gtl1.resize(fe1->n_components(), numbers::invalid_unsigned_int);
    gtl2.resize(fe2->n_components(), numbers::invalid_unsigned_int);

    for (unsigned int i = 0, j = 0; i < comps1.size(); ++i)
      if (comps1[i])
        gtl1[i] = j++;

    for (unsigned int i = 0, j = 0; i < comps2.size(); ++i)
      if (comps2[i])
        gtl2[i] = j++;

    n_comps =
      std::min(comps1.n_selected_components(), comps2.n_selected_components());

    if (auto tr =
          dynamic_cast<const parallel::distributed::Triangulation<spacedim> *>(
            &(*tria1)))
      {
        communicator = tr->get_communicator();
        mpi_proc     = Utilities::MPI::this_mpi_process(communicator);
        n_mpi_procs  = Utilities::MPI::n_mpi_processes(communicator);


        // Distribute the local points to the processor that owns
        // them on the triangulation
        auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
          *tria1, IteratorFilters::LocallyOwnedCell());

        global_bounding_boxes =
          Utilities::MPI::all_gather(communicator, my_bounding_box);
      }
    else
      {
        communicator = MPI_COMM_WORLD;
        mpi_proc     = 0;
        n_mpi_procs  = 1;

        global_bounding_boxes = {GridTools::compute_mesh_predicate_bounding_box(
          *tria1, IteratorFilters::LocallyOwnedCell())};
      }
  }



  template <int dim, int spacedim>
  DoFHandlerCoupling<dim, spacedim>::DoFHandlerCoupling(
    const DoFHandler<spacedim> &     dh1,
    const DoFHandler<dim, spacedim> &dh2,
    const ComponentMask &            comps1,
    const ComponentMask &            comps2,
    const Mapping<spacedim> &        mapping1,
    const Mapping<dim, spacedim> &   mapping2)
    : DoFHandlerCoupling(
        *(new GridTools::Cache<spacedim>(dh1.get_triangulation(), mapping1)),
        *(new GridTools::Cache<dim, spacedim>(dh2.get_triangulation(),
                                              mapping2)),
        dh1,
        dh2,
        comps1,
        comps2,
        mapping1,
        mapping2)
  {
    cache1_ptr.reset(cache1);
    cache2_ptr.reset(cache2);
  }



  template <int dim, int spacedim>
  void
  DoFHandlerCoupling<dim, spacedim>::possibly_generate_particle_handler(
    bool                   reuse_internal_data_structures,
    const Quadrature<dim> &quadrature) const
  {
    if (quadrature.size() == 0)
      {
        if (reuse_internal_data_structures == false || !particle_handler)
          {
            auto tr = dynamic_cast<
              const parallel::distributed::Triangulation<spacedim> *>(
              &(*tria1));

            particle_handler =
              std::make_unique<Particles::ParticleHandler<dim, spacedim>>(
                *tr, *mapping1, n_comps + 1);

            ComponentMask temp_mask(fe2->n_components(), false);
            unsigned int  first_comp = 0;
            for (; first_comp < fe2->n_components(); ++first_comp)
              if (comps2[first_comp])
                {
                  temp_mask.set(first_comp, true);
                  break;
                }

            std::map<types::global_dof_index, Point<spacedim>>
              support_points_map;
            DoFTools::map_dofs_to_support_points(*mapping2,
                                                 *dh2,
                                                 support_points_map,
                                                 temp_mask);

            const auto dh2_selected_indices =
              DoFTools::extract_dofs_per_component(
                *dh2, DoFTools::OwnershipType::owned, comps2);

            std::vector<Point<spacedim>> support_points_vec;
            support_points_vec.reserve(support_points_map.size());


            for (auto const &element : support_points_map)
              if (dh2_selected_indices[first_comp].is_element(element.first))
                support_points_vec.emplace_back(element.second);

            std::vector<double> properties(support_points_vec.size() *
                                           (n_comps + 1));

            std::vector<typename IndexSet::ElementIterator> elements;

            unsigned int j = 0;
            for (unsigned int c = first_comp; c < comps2.size() && j < n_comps;
                 ++c)
              if (comps2[c])
                elements.emplace_back(dh2_selected_indices[c].begin());

            for (unsigned int i = 0; i < support_points_vec.size(); ++i)
              {
                properties[i * (n_comps + 1)] = mpi_proc;
                for (unsigned int c = 0; c < n_comps; ++c)
                  properties[i * (n_comps + 1) + c + 1] = *(elements[c]++);
              }

            auto cpu_to_index =
              particle_handler->insert_global_particles(support_points_vec,
                                                        global_bounding_boxes,
                                                        properties);
          }
      }
    else
      {
        if (reuse_internal_data_structures == false ||
            !quadrature_particle_handler)
          {
            auto tr1 = dynamic_cast<
              const parallel::distributed::Triangulation<spacedim> *>(
              &(*tria1));
            Assert(tr1,
                   ExcMessage("Only available for parallel distributed "
                              "triangulations"));

            auto tr2 = dynamic_cast<
              const parallel::distributed::Triangulation<spacedim> *>(
              &(*tria2));

            Assert(tr2,
                   ExcMessage("Only available for parallel distributed "
                              "triangulations"));

            const unsigned int n_properties =
              (dim == spacedim) ? 1 : spacedim + 1;

            quadrature_particle_handler =
              std::make_unique<Particles::ParticleHandler<dim, spacedim>>(
                *tr1, *mapping1, n_properties);


            std::vector<Point<spacedim>> quadrature_points_vec(
              quadrature.size() * tr2->n_locally_owned_active_cells());

            std::vector<double> properties(quadrature.size() *
                                           tr2->n_locally_owned_active_cells() *
                                           n_properties);

            UpdateFlags flags = update_JxW_values | update_quadrature_points;
            if (spacedim > dim)
              flags |= update_normal_vectors;
            FEValues<dim, spacedim> fe_v(*mapping2, *fe2, quadrature, flags);

            unsigned int cell_index = 0;
            for (const auto &cell : dh2->active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  fe_v.reinit(cell);
                  const auto &points = fe_v.get_quadrature_points();
                  const auto &JxW    = fe_v.get_JxW_values();

                  for (unsigned int q = 0; q < points.size(); ++q)
                    {
                      quadrature_points_vec[cell_index * points.size() + q] =
                        points[q];
                      properties[cell_index * n_properties * points.size() +
                                 q * n_properties] = JxW[q];
                      if (dim < spacedim)
                        for (unsigned int d = 0; d < spacedim; ++d)
                          {
                            properties[cell_index * n_properties *
                                         points.size() +
                                       q * n_properties + 1 + d] =
                              fe_v.normal_vector(q)[d];
                          }
                    }
                  ++cell_index;
                }
            auto cpu_to_index =
              quadrature_particle_handler->insert_global_particles(
                quadrature_points_vec, global_bounding_boxes, properties);
          }
      }
  }

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif
