// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_particles_utilities
#define dealii_particles_utilities

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparsity_pattern_base.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/particles/particle_handler.h>

DEAL_II_NAMESPACE_OPEN



namespace Particles
{
  /**
   * A namespace for functions offering tools to handle ParticleHandler objects
   * and their coupling with DoFHandler objects.
   */
  namespace Utilities
  {
    /**
     * Create an interpolation sparsity pattern for particles.
     *
     * Given a triangulation representing the domain $\Omega$, a particle
     * handler of particles in $\Omega$, and a scalar finite element space
     * $V(\Omega) = \text{span}\{v_j\}_{j=0}^n$, compute the sparsity pattern
     * that would be necessary to assemble the matrix
     * \f[
     * M_{i,j} \dealcoloneq v_j(x_i) ,
     * \f]
     * where $V(\Omega)$ is the finite element space associated with the
     * `space_dh`, and the index `i` is given by the particle id whose position
     * is `x_i`.
     *
     * In the case of vector valued finite element spaces, the components on
     * which interpolation must be performed can be selected using a component
     * mask. Only primitive finite element spaces are supported.
     *
     * When selecting more than one component, the resulting sparsity will have
     * dimension equal to `particle_handler.n_global_particles() *
     * mask.n_selected_components()` times `space_dh.n_dofs()`, and the
     * corresponding matrix entries are given by
     * \f[
     *  M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j},
     * \f]
     * where `comp_j` is the only non zero component of the vector valued basis
     * function `v_j` (equal to `fe.system_to_component_index(j).first`),
     * `k` corresponds to its index within the selected components of the mask,
     * and $e_{comp_j}$ is the unit vector in the direction `comp_j`.
     *
     * The `sparsity` is filled by locating the position of the particle with
     * index `i` within the particle handler with respect to the embedding
     * triangulation $\Omega$, and coupling it with all the local degrees of
     * freedom specified in the component mask @p space_comps, following the
     * ordering in which they are selected in the mask @p space_comps.
     *
     * If a particle does not fall within $\Omega$, it is ignored, and the
     * corresponding rows of the sparsity will be empty.
     *
     * Constraints of the form supported by the AffineConstraints class may be
     * supplied with the @p constraints argument. The method
     * AffineConstraints::add_entries_local_to_global() is used to fill the
     * final sparsity pattern.
     */
    template <int dim, int spacedim, typename number = double>
    void
    create_interpolation_sparsity_pattern(
      const DoFHandler<dim, spacedim>                 &space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      SparsityPatternBase                             &sparsity,
      const AffineConstraints<number>                 &constraints =
        AffineConstraints<number>(),
      const ComponentMask &space_comps = {});

    /**
     * Create an interpolation matrix for particles.
     *
     * Given a triangulation representing the domains $\Omega$, a particle
     * handler of particles in $\Omega$, and a scalar finite element space
     * $V(\Omega) = \text{span}\{v_j\}_{j=0}^n$, compute the matrix
     * \f[
     * M_{ij} \dealcoloneq v_j(x_i) ,
     * \f]
     * where $V(\Omega)$ is the finite element space associated with the
     * `space_dh`, and the index `i` is given by the particle id whose position
     * is `x_i`.
     *
     * In the case of vector valued finite element spaces, the components on
     * which interpolation must be performed can be selected using a component
     * mask. Only primitive finite element spaces are supported.
     *
     * When selecting more than one component, the resulting sparsity will have
     * dimension equal to `particle_handler.n_global_particles() *
     * mask.n_selected_components()` times `space_dh.n_dofs()`, and the
     * corresponding matrix entries are given by
     * \f[
     *  M_{(i*n_comps+k),j} \dealcoloneq v_j(x_i) \cdot e_{comp_j},
     * \f]
     * where `comp_j` is the only non zero component of the vector valued basis
     * function `v_j` (equal to `fe.system_to_component_index(j).first`),
     * `k` corresponds to its index within the selected components of the mask,
     * and $e_{comp_j}$ is the unit vector in the direction `comp_j`.
     *
     * The matrix is filled by locating the position of the particle with
     * index `i` within the particle handler with respect to the embedding
     * triangulation $\Omega$, and coupling it with all the local degrees of
     * freedom specified in the component mask @p space_comps, following the
     * ordering in which they are selected in the mask @p space_comps.
     *
     * If a particle does not fall within $\Omega$, it is ignored, and the
     * corresponding rows of the matrix will be zero.
     *
     * Constraints of the form supported by the AffineConstraints class may be
     * supplied with the @p constraints argument. The method
     * AffineConstraints::distribute_local_to_global() is used to distribute
     * the entries of the matrix to respect the given constraints.
     */
    template <int dim, int spacedim, typename MatrixType>
    void
    create_interpolation_matrix(
      const DoFHandler<dim, spacedim>                 &space_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      MatrixType                                      &matrix,
      const AffineConstraints<typename MatrixType::value_type> &constraints =
        AffineConstraints<typename MatrixType::value_type>(),
      const ComponentMask &space_comps = {});

    /**
     * Given a DoFHandler and a ParticleHandler, interpolate a vector field
     * at the position of the particles. The result is stored in an output
     * vector whose size corresponds to the number of locally owned particles *
     * number of active components
     *
     * @param[in] field_dh The DOF Handler which was used to generate the
     * field vector that is to be interpolated.
     *
     * @param[in] particle_handler The particle handler whose particle serve as
     * the interpolation points.
     *
     * @param[in] field_vector The vector of the field to be interpolated. This
     * vector must be coherent with the dof_handler provided
     *
     * @param[in,out] interpolated_field The interpolated value of the field at
     * the position of the particles. The size of the vector must be
     * n_locally_owned_particles times the n_components
     *
     * @param[in] field_comps An optional component mask that decides which
     * subset of the vector fields are interpolated
     */
    template <int dim,
              int spacedim,
              typename InputVectorType,
              typename OutputVectorType>
    void
    interpolate_field_on_particles(
      const DoFHandler<dim, spacedim>                 &field_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      const InputVectorType                           &field_vector,
      OutputVectorType                                &interpolated_field,
      const ComponentMask                             &field_comps = {})
    {
      if (particle_handler.n_locally_owned_particles() == 0)
        {
          interpolated_field.compress(VectorOperation::add);
          return; // nothing else to do here
        }

      const auto &fe       = field_dh.get_fe();
      auto        particle = particle_handler.begin();

      // Take care of components
      const ComponentMask comps =
        (field_comps.size() == 0 ? ComponentMask(fe.n_components(), true) :
                                   field_comps);
      AssertDimension(comps.size(), fe.n_components());
      const auto n_comps = comps.n_selected_components();

      AssertDimension(field_vector.size(), field_dh.n_dofs());
      AssertDimension(interpolated_field.size(),
                      particle_handler.get_next_free_particle_index() *
                        n_comps);

      // Global to local indices
      std::vector<unsigned int> space_gtl(fe.n_components(),
                                          numbers::invalid_unsigned_int);
      for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
        if (comps[i])
          space_gtl[i] = j++;

      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

      while (particle != particle_handler.end())
        {
          const auto &cell = particle->get_surrounding_cell();
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, &field_dh);
          dh_cell->get_dof_indices(dof_indices);
          const auto pic = particle_handler.particles_in_cell(cell);

          Assert(pic.begin() == particle, ExcInternalError());
          for (unsigned int i = 0; particle != pic.end(); ++particle, ++i)
            {
              const Point<dim> reference_location =
                particle->get_reference_location();

              const auto id = particle->get_id();

              for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                {
                  const auto comp_j =
                    space_gtl[fe.system_to_component_index(j).first];
                  if (comp_j != numbers::invalid_unsigned_int)
                    interpolated_field[id * n_comps + comp_j] +=
                      fe.shape_value(j, reference_location) *
                      field_vector(dof_indices[j]);
                }
            }
        }
      interpolated_field.compress(VectorOperation::add);
    }


    /**
     * Given a DoFHandler and a ParticleHandler, interpolate a vector field
     * at the position of the particles. The result is stored in an output
     * vector whose size corresponds to the number of locally owned particles *
     * number of active components. This version relies on FePointEvaluator
     * to interpolate the field at the particle location efficiently. The
     * drawback is that this requires knowledge of the number of components at
     * compile time.
     *
     * @param[in] field_dh The DOF Handler which was used to generate the
     * field vector that is to be interpolated.
     *
     * @param[in] particle_handler The particle handler whose particle serve as
     * the interpolation points.
     *
     * @param[in] field_vector The vector of the field to be interpolated. This
     * vector must be coherent with the dof_handler provided
     *
     * @param[in,out] interpolated_field The interpolated value of the field at
     * the position of the particles. The size of the vector must be
     * n_locally_owned_particles times the n_components
     *
     * @param[in] field_comps An optional component mask that decides which
     * subset of the vector fields are interpolated
     *
     * @param[in] mapping The mapping used to transform the reference cell to
     *            the real cell. Defaults to the default mapping of the
     *            hypercube reference cell for the given dimension.
     *
     * @note This function is more efficient than interpolate_field_on_particles,
     *        but requires @p n_components to be specified as a
     *       compile-time template argument.
     */

    template <int n_components,
              int dim,
              int spacedim,
              typename InputVectorType,
              typename OutputVectorType>
    void
    interpolate_field_on_particles_fast(
      const DoFHandler<dim, spacedim>                 &field_dh,
      const Particles::ParticleHandler<dim, spacedim> &particle_handler,
      const InputVectorType                           &field_vector,
      OutputVectorType                                &interpolated_field,
      const ComponentMask                             &field_comps = {},
      const Mapping<dim, spacedim>                    &mapping =
        (ReferenceCells::get_hypercube<dim>()
#ifndef _MSC_VER
           .template get_default_linear_mapping<spacedim>()
#else
           .ReferenceCell<dim>::get_default_linear_mapping<spacedim>()
#endif
           ))
    {
      if (particle_handler.n_locally_owned_particles() == 0)
        {
          interpolated_field.compress(VectorOperation::add);
          return; // nothing else to do here
        }

      const auto &fe       = field_dh.get_fe();
      auto        particle = particle_handler.begin();

      // Take care of components
      const ComponentMask comps =
        (field_comps.size() == 0 ? ComponentMask(n_components, true) :
                                   field_comps);
      AssertDimension(comps.size(), n_components);
      const auto n_selected_comps = comps.n_selected_components();

      AssertDimension(field_vector.size(), field_dh.n_dofs());
      AssertDimension(interpolated_field.size(),
                      particle_handler.get_next_free_particle_index() *
                        n_selected_comps);

      // Global to local indices
      std::vector<unsigned int> space_gtl(fe.n_components(),
                                          numbers::invalid_unsigned_int);
      for (unsigned int i = 0, j = 0; i < space_gtl.size(); ++i)
        if (comps[i])
          space_gtl[i] = j++;

      std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

      // Generate an evaluator that will be used to interpolate the fields at
      // the particle location.
      FEPointEvaluation<n_components, dim, spacedim> evaluator(mapping,
                                                               fe,
                                                               update_values);
      std::vector<Point<dim>>            particle_reference_locations;
      std::vector<types::particle_index> particle_indices;
      Vector<double>                     local_dof_values(fe.dofs_per_cell);

      while (particle != particle_handler.end())
        {
          const auto &cell = particle->get_surrounding_cell();
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, &field_dh);
          dh_cell->get_dof_indices(dof_indices);
          const auto pic = particle_handler.particles_in_cell(cell);
          dh_cell->get_dof_values(field_vector, local_dof_values);

          Assert(pic.begin() == particle, ExcInternalError());

          // Gather the reference location  and ids of all particles
          particle_reference_locations.clear();
          particle_indices.clear();
          for (auto &p : pic)
            {
              particle_reference_locations.emplace_back(
                p.get_reference_location());
              particle_indices.emplace_back(p.get_id());
            }

          evaluator.reinit(cell, particle_reference_locations);
          evaluator.evaluate(make_array_view(local_dof_values),
                             EvaluationFlags::values);
          for (unsigned int particle_index = 0; particle != pic.end();
               ++particle, ++particle_index)
            {
              const unsigned id = particle_indices[particle_index];

              if constexpr (n_components == 1)
                {
                  interpolated_field[id] = evaluator.get_value(particle_index);
                }
              else
                {
                  unsigned int j_comp = 0;
                  for (unsigned int j = 0; j < n_components; ++j)
                    {
                      if (field_comps[j])
                        {
                          interpolated_field[id * n_selected_comps + j_comp] =
                            evaluator.get_value(particle_index)[j];
                          j_comp++;
                        }
                    }
                }
            }
        }
    }

  } // namespace Utilities
} // namespace Particles
DEAL_II_NAMESPACE_CLOSE

#endif
