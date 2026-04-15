// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_non_matching_dof_handler_coupling
#define dealii_non_matching_dof_handler_coupling

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
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
   * triangulations.
   *
   * Let $V_h = \mathrm{span}\{\varphi_i\}_{i=0}^{N_1-1}$ be the finite
   * element space associated with @p dh1 (the embedding/background space), and
   * let $W_h = \mathrm{span}\{\psi_a\}_{a=0}^{N_2-1}$ be the finite element
   * space associated with @p dh2 (the immersed/coupled space). The class builds
   * particle-based search and transfer data so that one can assemble coupling
   * operators on arbitrarily distributed and fully distributed non-matching
   * meshes.
   *
   * In particular, this class supports:
   * - interpolation/restriction operators of the form
   *   \f[
   *   u_h(x_a) = \sum_{i=0}^{N_1-1} U_i\,\varphi_i(x_a),
   *   \f]
   *   where $x_a$ are support points associated with selected DoFs of @p dh2;
   * - penalty-type Nitsche restriction terms on the immersed manifold/cells,
   *   resulting in a matrix
   *   \f[
   *   A^{\gamma}_{ij} = \gamma\int_{\Gamma_h} \varphi_i\,\varphi_j\,d\Gamma,
   *   \f]
   *   and right-hand side
   *   \f[
   *   b^{\gamma}_{i} = \gamma\int_{\Gamma_h} g\,\varphi_i\,d\Gamma,
   *   \f]
   *   where $\gamma$ is the penalty parameter and $g$ is the prescribed
   *   restriction data.
   *
   * These operators are central in immersed and mixed-dimensional couplings:
   * they encode, respectively, pointwise transfer from the background field to
   * immersed DoFs and weak enforcement of constraints through consistent
   * variational terms.
   */
  template <int dim, int spacedim>
  class DoFHandlerCoupling : public EnableObserverPointer
  {
  public:
    /**
     * Construct a coupling object using externally managed grid caches.
     *
     * @param cache1 Geometric search cache for the background triangulation.
     * @param cache2 Geometric search cache for the immersed triangulation.
     * @param dh1 DoFHandler defining the source/background space $V_h$.
     * @param dh2 DoFHandler defining the immersed/target space $W_h$.
     * @param comps1 Active components of @p dh1 used in coupling.
     * @param comps2 Active components of @p dh2 used in coupling.
     * @param mapping1 Mapping used to interpret geometry on @p dh1.
     * @param mapping2 Mapping used to interpret geometry on @p dh2.
     *
     * If component masks are omitted, all components are considered active.
     * Active components are matched in order: the first selected component in
     * @p comps1 couples with the first selected component in @p comps2, and so
     * on, up to the minimum number of selected components.
     */
    DoFHandlerCoupling(
      const GridTools::Cache<spacedim>      &cache1,
      const GridTools::Cache<dim, spacedim> &cache2,
      const DoFHandler<spacedim>            &dh1,
      const DoFHandler<dim, spacedim>       &dh2,
      const ComponentMask                   &comps1_ = ComponentMask(),
      const ComponentMask                   &comps2_ = ComponentMask(),
      const Mapping<spacedim> &mapping1 = StaticMappingQ1<spacedim>::mapping,
      const Mapping<dim, spacedim> &mapping2 =
        StaticMappingQ1<dim, spacedim>::mapping);

    /**
     * Construct a coupling object by creating internal grid caches.
     *
     * This constructor is equivalent to the cache-based constructor, but
     * allocates and owns the two GridTools::Cache objects internally.
     *
     * @param dh1 DoFHandler defining the source/background space $V_h$.
     * @param dh2 DoFHandler defining the immersed/target space $W_h$.
     * @param comps1 Active components of @p dh1 used in coupling.
     * @param comps2 Active components of @p dh2 used in coupling.
     * @param mapping1 Mapping used on the background mesh.
     * @param mapping2 Mapping used on the immersed mesh.
     */
    DoFHandlerCoupling(
      const DoFHandler<spacedim>      &dh1,
      const DoFHandler<dim, spacedim> &dh2,
      const ComponentMask             &comps1 = ComponentMask(),
      const ComponentMask             &comps2 = ComponentMask(),
      const Mapping<spacedim> &mapping1 = StaticMappingQ1<spacedim>::mapping,
      const Mapping<dim, spacedim> &mapping2 =
        StaticMappingQ1<dim, spacedim>::mapping);



    /**
     * Build the sparsity pattern for the interpolation/restriction matrix.
     *
     * The resulting matrix $I$ maps DoFs of @p dh1 to selected DoFs of @p dh2,
     * with entries
     * \f[
     * I_{a i} = \varphi_i(x_a),
     * \f]
     * where $x_a$ is the support point associated with the $a$-th selected DoF
     * of @p dh2 and $\varphi_i$ is the $i$-th shape function of @p dh1.
     *
     * Only nonzero couplings $(a,i)$ are inserted in @p sparsity, accounting
     * for component matching and for the provided affine constraints.
     *
     * @param[out] sparsity Sparsity pattern to be filled.
     * @param[in] constraints Constraints used while inserting matrix entries.
     * @param[in] reuse_internal_data_structures If true, reuse cached particle
     *   data when available.
     */
    template <class Sparsity, typename number = double>
    void
    create_interpolation_sparsity_pattern(
      Sparsity                        &sparsity,
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
          const auto &cell = particle->get_surrounding_cell();
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
     * Assemble penalty-type Nitsche restriction contributions.
     *
     * For each active component, this function adds to @p matrix and @p rhs the
     * terms
     * \f[
     * A^{\gamma}_{ij} = \gamma\int_{\Gamma_h} \varphi_i\,\varphi_j\,d\Gamma,
     * \qquad
     * b^{\gamma}_i = \gamma\int_{\Gamma_h} g\,\varphi_i\,d\Gamma,
     * \f]
     * where $\Gamma_h$ is the immersed integration set, $\gamma$ is
     * @p penalty_term, and $g$ is @p rhs_function. The basis functions are the
     * shape functions of @p dh1, evaluated at quadrature points on the immersed
     * mesh and integrated in the background mesh using the particle-based
     * quadrature representation.
     *
     * The matrix is the penalty block that enforces restriction/Dirichlet data
     * in a weak sense, while the right-hand side contains the prescribed target
     * values projected with the same basis.
     *
     * @param[in] quadrature Reference quadrature on the immersed cells.
     * @param[in] rhs_function Prescribed field $g$ in physical coordinates.
     * @param[in] penalty_term Penalty parameter $\gamma$.
     * @param[in,out] matrix Global matrix receiving Nitsche penalty entries.
     * @param[in,out] rhs Global right-hand side receiving penalty forcing.
     * @param[in] constraints Constraints used in local-to-global distribution.
     * @param[in] reuse_internal_data_structures If true, reuse cached particle
     *   data when available.
     */
    template <class MatrixType, class VectorType, typename number = double>
    void
    create_nitsche_restriction(
      const Quadrature<dim>           &quadrature,
      const Function<spacedim>        &rhs_function,
      const double                    &penalty_term,
      MatrixType                      &matrix,
      VectorType                      &rhs,
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
      Vector<double>     local_rhs(fe1->dofs_per_cell);

      auto particle = quadrature_particle_handler->begin();
      while (particle != quadrature_particle_handler->end())
        {
          local_matrix     = 0;
          local_rhs        = 0;
          const auto &cell = particle->get_surrounding_cell();
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



    /**
     * Assemble the interpolation/restriction matrix between @p dh1 and @p dh2.
     *
     * The assembled matrix $I$ satisfies
     * \f[
     * \mathbf{u}_{\mathrm{immersed}} = I\,\mathbf{u}_{\mathrm{background}},
     * \f]
     * and has entries
     * \f[
     * I_{a i} = \varphi_i(x_a),
     * \f]
     * with $x_a$ support points of selected DoFs in @p dh2.
     *
     * Significance: this matrix is the discrete trace/transfer operator from
     * the background field to immersed DoFs. It is typically used for
     * interpolation constraints, projection-like couplings, and
     * mixed-dimensional transfer operators.
     *
     * @param[in,out] matrix Global matrix with shape
     *   $(\text{n\_dofs}(dh2),\text{n\_dofs}(dh1))$.
     * @param[in] constraints Constraints applied during assembly.
     */
    template <class Matrix, typename number = double>
    void
    create_interpolation_matrix(Matrix                          &matrix,
                                const AffineConstraints<number> &constraints =
                                  AffineConstraints<number>()) const
    {
      if (particle_handler->n_locally_owned_particles() == 0)
        {
          matrix.compress(VectorOperation::add);
          return; // nothing else to do here
        }

      // Expect: matrix rows = dofs on dh2, matrix columns = dofs on dh1.
      AssertDimension(matrix.m(), dh2->n_dofs());
      AssertDimension(matrix.n(), dh1->n_dofs());

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
          const auto &cell = particle->get_surrounding_cell();
          const auto &dh_cell =
            typename DoFHandler<dim, spacedim>::cell_iterator(*cell, dh1);
          dh_cell->get_dof_indices(dof_indices1);

          const auto pic         = particle_handler->particles_in_cell(cell);
          const auto n_particles = particle_handler->n_particles_in_cell(cell);
          dof_indices2.resize(n_particles * n_comps);
          local_matrix.reinit({n_particles * n_comps, fe1->dofs_per_cell});
          local_matrix = 0;

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

    /**
     * Return the particle handler representing immersed support points mapped
     * into the background triangulation.
     *
     * Particles store, in their properties, ownership metadata and selected
     * DoF indices of @p dh2 used to build interpolation sparsity and matrix
     * entries.
     */
    const Particles::ParticleHandler<spacedim> &
    get_particle_handler() const
    {
      Assert(particle_handler,
             ExcMessage("You need to call possibly_generate_particle_handler "
                        "function first."));
      return *particle_handler;
    }

    /**
     * Return the quadrature particle handler used for immersed integration.
     *
     * Particles represent physical quadrature points on the immersed mesh,
     * located in the background mesh, and store at least $JxW$ in their
     * properties (and normals when `dim < spacedim`).
     */
    const Particles::ParticleHandler<spacedim> &
    get_quadrature_particle_handler() const
    {
      Assert(quadrature_particle_handler,
             ExcMessage("You need to call possibly_generate_particle_handler "
                        "function first."));
      return *quadrature_particle_handler;
    }

    /**
     * Lazily create particle-based internal data structures.
     *
     * If @p quadrature is empty, the function builds @ref particle_handler from
     * support points of selected @p dh2 DoFs. Otherwise, it builds
     * @ref quadrature_particle_handler from immersed quadrature points and stores
     * geometric weights (and normals in codimension one).
     *
     * @param[in] reuse_internal_data_structures If true, skip re-creation when
     *   handlers are already available.
     * @param[in] quadrature Optional quadrature rule used to generate
     *   integration particles.
     */
    void
    possibly_generate_particle_handler(
      bool                   reuse_internal_data_structures,
      const Quadrature<dim> &quadrature = Quadrature<dim>()) const;

    /**
     * Set positions of quadrature particles using an explicit vector.
     *
     * This updates quadrature-point particle positions used in immersed
     * integration. With @p displace_particles set to true, values are treated
     * as displacements; otherwise they are interpreted as absolute coordinates.
     *
     * @param[in] positions_vector New positions/displacements, one per local
     *   quadrature particle.
     * @param[in] displace_particles Interpret @p positions_vector as
     *   displacements (true) or absolute coordinates (false).
     */
    void
    set_quadrature_particles_positions(
      const std::vector<Tensor<1, spacedim>> &positions_vector,
      const bool                              displace_particles = true) const
    {
      // There should be a position per per local particles
      AssertDimension(positions_vector.size(),
                      quadrature_particle_handler->n_local_particles());
      quadrature_particle_handler->set_particle_positions(positions_vector,
                                                          displace_particles);
    }

    /**
     * Set positions of interpolation particles using an explicit vector.
     *
     * This updates support-point particles used by interpolation sparsity and
     * interpolation matrix assembly.
     *
     * @param[in] positions_vector New positions/displacements, one per local
     *   interpolation particle.
     * @param[in] displace_particles Interpret @p positions_vector as
     *   displacements (true) or absolute coordinates (false).
     */
    void
    set_particles_positions(
      const std::vector<Tensor<1, spacedim>> &positions_vector,
      const bool                              displace_particles = true) const
    {
      // There should be a position per per local particles
      AssertDimension(positions_vector.size(),
                      particle_handler->n_local_particles());
      particle_handler->set_particle_positions(positions_vector,
                                               displace_particles);
    }

    /**
     * Set quadrature-particle positions from a vector-valued function.
     *
     * @param[in] function Function prescribing displacements or absolute
     *   coordinates.
     * @param[in] displace_particles Interpret @p function as displacement field
     *   (true) or position field (false).
     */
    void
    set_quadrature_particles_positions(
      const Function<spacedim> &function,
      const bool                displace_particles = true) const
    {
      quadrature_particle_handler->set_particle_positions(function,
                                                          displace_particles);
    }

    /**
     * Set interpolation-particle positions from a vector-valued function.
     *
     * @param[in] function Function prescribing displacements or absolute
     *   coordinates.
     * @param[in] displace_particles Interpret @p function as displacement field
     *   (true) or position field (false).
     */
    void
    set_particles_positions(const Function<spacedim> &function,
                            const bool displace_particles = true) const
    {
      particle_handler->set_particle_positions(function, displace_particles);
    }


  private:
    /** Mapping used to evaluate shape functions on the background mesh. */
    ObserverPointer<const Mapping<spacedim>> mapping1;
    /** Mapping used to evaluate shape functions and quadrature on the immersed
     * mesh. */
    ObserverPointer<const Mapping<dim, spacedim>> mapping2;

    /** Background triangulation carrying space $V_h$. */
    ObserverPointer<const Triangulation<spacedim>> tria1;
    /** Immersed triangulation carrying space $W_h$. */
    ObserverPointer<const Triangulation<dim, spacedim>> tria2;

    /** DoFHandler for the background/source finite element space. */
    ObserverPointer<const DoFHandler<spacedim>> dh1;
    /** DoFHandler for the immersed/target finite element space. */
    ObserverPointer<const DoFHandler<dim, spacedim>> dh2;

    /** Finite element of @ref dh1 defining basis functions $\varphi_i$. */
    ObserverPointer<const FiniteElement<spacedim>> fe1;
    /** Finite element of @ref dh2 defining basis functions $\psi_a$. */
    ObserverPointer<const FiniteElement<dim, spacedim>> fe2;

    /** Optional ownership of internally created cache for @ref tria1. */
    std::unique_ptr<const GridTools::Cache<spacedim>> cache1_ptr;
    /** Optional ownership of internally created cache for @ref tria2. */
    std::unique_ptr<const GridTools::Cache<dim, spacedim>> cache2_ptr;

    /** Geometric search cache associated with @ref tria1. */
    ObserverPointer<const GridTools::Cache<spacedim>> cache1;
    /** Geometric search cache associated with @ref tria2. */
    ObserverPointer<const GridTools::Cache<dim, spacedim>> cache2;

    /**
     * Cached cell/quadrature pairs for immersed-surface integration kernels.
     */
    std::vector<std::pair<typename DoFHandler<spacedim>::cell_iterator,
                          ImmersedSurfaceQuadrature<spacedim>>>
      cells_and_quadratures;

    /** Active component mask for background DoFs. */
    ComponentMask comps1;
    /** Active component mask for immersed DoFs. */
    ComponentMask comps2;

    /** Number of coupled components:
     * $\min(\#\mathrm{comps1},\#\mathrm{comps2})$. */
    unsigned int n_comps;
    /** Rank of the current MPI process in @ref communicator. */
    unsigned int mpi_proc;
    /** Number of MPI processes participating in @ref communicator. */
    unsigned int n_mpi_procs;
    /** MPI communicator used by the background triangulation. */
    MPI_Comm communicator;

    /**
     * Background global-to-localized component map: for each FE component of
     * @ref fe1, stores the contiguous coupled-component index or
     * `numbers::invalid_unsigned_int` if inactive.
     */
    std::vector<unsigned int> gtl1;
    /**
     * Immersed global-to-localized component map analogous to @ref gtl1.
     */
    std::vector<unsigned int> gtl2;

    /**
     * Bounding boxes of locally owned background cells on all MPI ranks,
     * used to route particles to owning processes.
     */
    std::vector<std::vector<BoundingBox<spacedim>>> global_bounding_boxes;

    /**
     * Particle representation of immersed support points used for
     * interpolation matrix/sparsity construction.
     */
    mutable std::unique_ptr<Particles::ParticleHandler<spacedim>>
      particle_handler;

    /**
     * Particle representation of immersed quadrature points used for Nitsche
     * and other quadrature-based coupling terms.
     */
    mutable std::unique_ptr<Particles::ParticleHandler<spacedim>>
      quadrature_particle_handler;
  };



  // Implementation
  template <int dim, int spacedim>
  DoFHandlerCoupling<dim, spacedim>::DoFHandlerCoupling(
    const GridTools::Cache<spacedim>      &cache1,
    const GridTools::Cache<dim, spacedim> &cache2,
    const DoFHandler<spacedim>            &dh1,
    const DoFHandler<dim, spacedim>       &dh2,
    const ComponentMask                   &comps1_,
    const ComponentMask                   &comps2_,
    const Mapping<spacedim>               &mapping1,
    const Mapping<dim, spacedim>          &mapping2)
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
    const DoFHandler<spacedim>      &dh1,
    const DoFHandler<dim, spacedim> &dh2,
    const ComponentMask             &comps1,
    const ComponentMask             &comps2,
    const Mapping<spacedim>         &mapping1,
    const Mapping<dim, spacedim>    &mapping2)
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
              DoFTools::locally_owned_dofs_per_component(*dh2, comps2);

            std::vector<Point<spacedim>> support_points_vec;
            support_points_vec.reserve(support_points_map.size());


            for (const auto &element : support_points_map)
              if (dh2_selected_indices[first_comp].is_element(element.first))
                support_points_vec.emplace_back(element.second);

            std::vector<std::vector<double>> properties(
              support_points_vec.size(), std::vector<double>(n_comps + 1));

            std::vector<typename IndexSet::ElementIterator> elements;

            unsigned int j = 0;
            for (unsigned int c = first_comp; c < comps2.size() && j < n_comps;
                 ++c)
              if (comps2[c])
                elements.emplace_back(dh2_selected_indices[c].begin());

            for (unsigned int i = 0; i < support_points_vec.size(); ++i)
              {
                properties[i][0] = mpi_proc;
                for (unsigned int c = 0; c < n_comps; ++c)
                  properties[i][c + 1] = static_cast<double>(*(elements[c]++));
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
            std::vector<std::vector<double>> quad_properties(
              quadrature_points_vec.size(), std::vector<double>(n_properties));
            for (unsigned int i = 0; i < quadrature_points_vec.size(); ++i)
              for (unsigned int j = 0; j < n_properties; ++j)
                quad_properties[i][j] = properties[i * n_properties + j];

            auto cpu_to_index =
              quadrature_particle_handler->insert_global_particles(
                quadrature_points_vec, global_bounding_boxes, quad_properties);
          }
      }
  }

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif
