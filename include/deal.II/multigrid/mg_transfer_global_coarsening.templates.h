// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mg_transfer_global_coarsening_templates_h
#define dealii_mg_transfer_global_coarsening_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mpi_consensus_algorithms.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/repartitioning_policy_tools.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/cell_id_translator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/matrix_free/evaluation_kernels.h>
#include <deal.II/matrix_free/evaluation_template_factory.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tensor_product_kernels.h>
#include <deal.II/matrix_free/vector_access_internal.h>

#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>

#include <boost/algorithm/string/join.hpp>

#include <limits>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  /**
   * Helper class to select the right templated implementation.
   *
   * @note This class is similar to internal::FEEvaluationFactory
   */
  class CellTransferFactory
  {
  public:
    static const unsigned int max_degree = 9;

    CellTransferFactory(const unsigned int degree_fine,
                        const unsigned int degree_coarse)
      : degree_fine(degree_fine)
      , degree_coarse(degree_coarse)
    {}

    template <typename Fu, unsigned int deg = 1>
    bool
    run(Fu &fu)
    {
      if ((degree_fine == (2 * deg)) && (degree_coarse == deg))
        fu.template run<2 * deg, deg>(); // h-MG (FE_Q)
      else if ((degree_fine == (2 * deg + 1)) && (degree_coarse == deg))
        fu.template run<2 * deg + 1, deg>(); // h-MG
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg / 2, 1u)))
        {
          constexpr unsigned int degree_coarse_used = std::max(deg / 2u, 1u);
          fu.template run<deg, degree_coarse_used>(); // p-MG: bisection
        }
      else if ((degree_fine == deg) && (degree_coarse == deg))
        fu.template run<deg, deg>(); // identity (nothing to do)
      else if ((degree_fine == deg) && (degree_coarse == std::max(deg - 1, 1u)))
        {
          constexpr unsigned int degree_coarse_used = std::max(deg - 1u, 1u);
          fu.template run<deg, degree_coarse_used>(); // p-MG: --
        }
      else if ((degree_fine == deg) && (degree_coarse == 1))
        fu.template run<deg, 1>(); // p-MG: jump to 1
      else if (deg < max_degree)
        return run<Fu, std::min(deg + 1, max_degree)>(fu); // try next degree
      else
        {
          // no match -> slow path
          fu.template run<-1, -1>(degree_fine, degree_coarse);
          return false; // indicate that slow path has been taken
        }

      return true; // indicate that fast path has been taken
    }

  private:
    const unsigned int degree_fine;
    const unsigned int degree_coarse;
  };

  /**
   * Helper class containing the cell-wise prolongation operation.
   */
  template <int dim, typename Number, typename Number2>
  class CellProlongator
  {
  public:
    CellProlongator(const AlignedVector<Number> &prolongation_matrix,
                    const Number2               *evaluation_data_coarse,
                    Number2                     *evaluation_data_fine)
      : prolongation_matrix(prolongation_matrix)
      , evaluation_data_coarse(evaluation_data_coarse)
      , evaluation_data_fine(evaluation_data_fine)
    {}

    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int degree_fine_   = numbers::invalid_unsigned_int,
        const unsigned int degree_coarse_ = numbers::invalid_unsigned_int)
    {
      Assert(prolongation_matrix.size() > 0, ExcNotImplemented());

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        dim,
        degree_coarse + 1,
        degree_fine + 1>::do_forward(1,
                                     prolongation_matrix,
                                     evaluation_data_coarse,
                                     evaluation_data_fine,
                                     degree_coarse_ + 1,
                                     degree_fine_ + 1);
    }

    void
    run_full(const unsigned int n_dofs_fine, const unsigned int n_dofs_coarse)
    {
      AssertDimension(prolongation_matrix.size(), n_dofs_coarse * n_dofs_fine);

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        1,
        0,
        0>::do_forward(1,
                       prolongation_matrix,
                       evaluation_data_coarse,
                       evaluation_data_fine,
                       n_dofs_coarse,
                       n_dofs_fine);
    }

  private:
    const AlignedVector<Number> &prolongation_matrix;
    const Number2               *evaluation_data_coarse;
    Number2                     *evaluation_data_fine;
  };

  /**
   * Helper class containing the cell-wise restriction operation.
   */
  template <int dim, typename Number, typename Number2>
  class CellRestrictor
  {
  public:
    CellRestrictor(const AlignedVector<Number> &prolongation_matrix,
                   Number2                     *evaluation_data_fine,
                   Number2                     *evaluation_data_coarse)
      : prolongation_matrix(prolongation_matrix)
      , evaluation_data_fine(evaluation_data_fine)
      , evaluation_data_coarse(evaluation_data_coarse)
    {}

    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int degree_fine_   = numbers::invalid_unsigned_int,
        const unsigned int degree_coarse_ = numbers::invalid_unsigned_int)
    {
      Assert(prolongation_matrix.size() > 0, ExcNotImplemented());

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        dim,
        degree_coarse == 0 ? -1 : (degree_coarse + 1),
        degree_fine == 0 ? -1 : (degree_fine + 1)>::
        do_backward(1,
                    prolongation_matrix,
                    false,
                    evaluation_data_fine,
                    evaluation_data_coarse,
                    degree_coarse_ + 1,
                    degree_fine_ + 1);
    }

    void
    run_full(const unsigned int n_dofs_fine, const unsigned int n_dofs_coarse)
    {
      AssertDimension(prolongation_matrix.size(), n_dofs_coarse * n_dofs_fine);

      internal::FEEvaluationImplBasisChange<
        internal::evaluate_general,
        internal::EvaluatorQuantity::value,
        1,
        0,
        0>::do_backward(1,
                        prolongation_matrix,
                        false,
                        evaluation_data_fine,
                        evaluation_data_coarse,
                        n_dofs_coarse,
                        n_dofs_fine);
    }

  private:
    const AlignedVector<Number> &prolongation_matrix;
    Number2                     *evaluation_data_fine;
    Number2                     *evaluation_data_coarse;
  };

  class CellProlongatorTest
  {
  public:
    template <int degree_fine, int degree_coarse>
    void
    run(const unsigned int = numbers::invalid_unsigned_int,
        const unsigned int = numbers::invalid_unsigned_int)
    {}
  };

} // namespace internal



namespace internal
{
  template <typename MeshType, typename OPType>
  DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
  void loop_over_active_or_level_cells(const MeshType    &tria,
                                       const unsigned int level,
                                       const OPType      &op)
  {
    if (level == numbers::invalid_unsigned_int)
      {
        for (const auto &cell : tria.active_cell_iterators())
          if (cell->is_locally_owned())
            op(cell);
      }
    else
      {
        for (const auto &cell : tria.cell_iterators_on_level(level))
          if (cell->is_locally_owned_on_level())
            op(cell);
      }
  }

  template <int dim>
  unsigned int
  compute_shift_within_children(const unsigned int child,
                                const unsigned int fe_shift_1d,
                                const unsigned int fe_degree)
  {
    // we put the degrees of freedom of all child cells in lexicographic
    // ordering
    unsigned int c_tensor_index[dim];
    unsigned int tmp = child;
    for (unsigned int d = 0; d < dim; ++d)
      {
        c_tensor_index[d] = tmp % 2;
        tmp /= 2;
      }
    const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
    unsigned int       factor          = 1;
    unsigned int       shift           = fe_shift_1d * c_tensor_index[0];
    for (unsigned int d = 1; d < dim; ++d)
      {
        factor *= n_child_dofs_1d;
        shift = shift + factor * fe_shift_1d * c_tensor_index[d];
      }
    return shift;
  }

  template <int dim>
  void
  get_child_offset(const unsigned int         child,
                   const unsigned int         fe_shift_1d,
                   const unsigned int         fe_degree,
                   std::vector<unsigned int> &local_dof_indices)
  {
    const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
    const unsigned int shift =
      compute_shift_within_children<dim>(child, fe_shift_1d, fe_degree);
    const unsigned int n_components =
      local_dof_indices.size() / Utilities::fixed_power<dim>(fe_degree + 1);
    const unsigned int n_scalar_cell_dofs =
      Utilities::fixed_power<dim>(n_child_dofs_1d);
    for (unsigned int c = 0, m = 0; c < n_components; ++c)
      for (unsigned int k = 0; k < (dim > 2 ? (fe_degree + 1) : 1); ++k)
        for (unsigned int j = 0; j < (dim > 1 ? (fe_degree + 1) : 1); ++j)
          for (unsigned int i = 0; i < (fe_degree + 1); ++i, ++m)
            local_dof_indices[m] = c * n_scalar_cell_dofs +
                                   k * n_child_dofs_1d * n_child_dofs_1d +
                                   j * n_child_dofs_1d + i + shift;
  }

  template <int dim>
  std::vector<std::vector<unsigned int>>
  get_child_offsets(const unsigned int n_dofs_per_cell_coarse,
                    const unsigned int fe_shift_1d,
                    const unsigned int fe_degree)
  {
    std::vector<std::vector<unsigned int>> cell_local_children_indices(
      GeometryInfo<dim>::max_children_per_cell,
      std::vector<unsigned int>(n_dofs_per_cell_coarse));
    for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; c++)
      get_child_offset<dim>(c,
                            fe_shift_1d,
                            fe_degree,
                            cell_local_children_indices[c]);
    return cell_local_children_indices;
  }

  template <int dim>
  std::vector<std::vector<unsigned int>>
  get_child_offsets_general(const unsigned int n_dofs_per_cell_coarse)
  {
    std::vector<std::vector<unsigned int>> cell_local_children_indices(
      GeometryInfo<dim>::max_children_per_cell,
      std::vector<unsigned int>(n_dofs_per_cell_coarse));
    for (unsigned int c = 0, k = 0;
         c < GeometryInfo<dim>::max_children_per_cell;
         c++)
      for (unsigned int d = 0; d < n_dofs_per_cell_coarse; ++d, ++k)
        cell_local_children_indices[c][d] = k;
    return cell_local_children_indices;
  }

  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<1>>
  create_1D_fe(const FiniteElement<dim, spacedim> &fe)
  {
    std::string fe_name = fe.get_name();
    {
      const std::size_t template_starts = fe_name.find_first_of('<');
      Assert(fe_name[template_starts + 1] ==
               (dim == 1 ? '1' : (dim == 2 ? '2' : '3')),
             ExcInternalError());
      fe_name[template_starts + 1] = '1';
    }
    return FETools::get_fe_by_name<1, 1>(fe_name);
  }

  template <int dim, int spacedim>
  FullMatrix<double>
  get_restriction_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const unsigned int                  child,
    RefinementCase<dim> ref_case = RefinementCase<dim>::isotropic_refinement)
  {
    auto matrix = fe.get_restriction_matrix(child, ref_case);

    for (unsigned int c_other = 0; c_other < child; ++c_other)
      {
        const auto &matrix_other = fe.get_restriction_matrix(c_other, ref_case);
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          {
            if (fe.restriction_is_additive(i) == true)
              continue;

            bool do_zero = false;
            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
              if (std::fabs(matrix_other(i, j)) > 1e-12)
                do_zero = true;

            if (do_zero)
              for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
                matrix(i, j) = 0.0;
          }
      }
    return matrix;
  }

  template <int dim>
  bool
  use_fast_hanging_node_algorithm(const DoFHandler<dim> &dof_handler_coarse,
                                  const unsigned int     mg_level_coarse)
  {
    // algorithm is only needed on active levels
    bool use_fast_hanging_node_algorithm =
      mg_level_coarse == numbers::invalid_unsigned_int;

    // algorithm can be only used on meshes consisting of hypercube and
    // simplices
    if (use_fast_hanging_node_algorithm)
      {
        const auto &reference_cells =
          dof_handler_coarse.get_triangulation().get_reference_cells();
        use_fast_hanging_node_algorithm =
          std::all_of(reference_cells.begin(),
                      reference_cells.end(),
                      [](const auto &r) {
                        return r.is_hyper_cube() || r.is_simplex();
                      });
      }

    // local p-refinement is not supported
    if (use_fast_hanging_node_algorithm)
      {
        const auto &fes = dof_handler_coarse.get_fe_collection();

        use_fast_hanging_node_algorithm &=
          std::all_of(fes.begin(), fes.end(), [&fes](const auto &fe) {
            return fes[0].compare_for_domination(fe) ==
                   FiniteElementDomination::Domination::
                     either_element_can_dominate;
          });
      }

    // check that all components are either supported or not
    if (use_fast_hanging_node_algorithm)
      {
        const std::vector<std::vector<bool>> supported_components = internal::
          MatrixFreeFunctions::HangingNodes<dim>::compute_supported_components(
            dof_handler_coarse.get_fe_collection());

        use_fast_hanging_node_algorithm &=
          std::any_of(supported_components.begin(),
                      supported_components.end(),
                      [](const auto &supported_components_per_fe) {
                        return std::all_of(supported_components_per_fe.begin(),
                                           supported_components_per_fe.end(),
                                           [](const auto &a) {
                                             return a == true;
                                           });
                      });

        use_fast_hanging_node_algorithm &= std::all_of(
          supported_components.begin(),
          supported_components.end(),
          [](const auto &supported_components_per_fe) {
            return std::all_of(supported_components_per_fe.begin(),
                               supported_components_per_fe.end(),
                               [&supported_components_per_fe](const auto &a) {
                                 return a == supported_components_per_fe[0];
                               });
          });
      }

    return use_fast_hanging_node_algorithm;
  }



  /**
   * A class behaving like DoFCellAccessor. Intended to be used for locally
   * relevant cell as a wrapper around DoFCellAccessor and for other cells
   * behaving as if the cell would be available.
   */
  class FineDoFHandlerViewCell
  {
  public:
    /**
     * Constructor.
     */
    FineDoFHandlerViewCell(
      const std::function<bool()> &has_children_function,
      const std::function<void(std::vector<types::global_dof_index> &)>
                                          &get_dof_indices_function,
      const std::function<unsigned int()> &active_fe_index_function,
      const std::function<std::uint8_t()> &refinement_case_function)
      : has_children_function(has_children_function)
      , get_dof_indices_function(get_dof_indices_function)
      , active_fe_index_function(active_fe_index_function)
      , refinement_case_function(refinement_case_function)
    {}

    /**
     * Return if cell has child.
     */
    bool
    has_children() const
    {
      return has_children_function();
    }

    /**
     * Get DoF indices.
     */
    void
    get_dof_indices(std::vector<types::global_dof_index> &dof_indices) const
    {
      get_dof_indices_function(dof_indices);
    }
    /**
     * Get active FE index.
     */
    unsigned int
    active_fe_index() const
    {
      return active_fe_index_function();
    }
    /**
     * Get refinement choice.
     */
    std::uint8_t
    refinement_case() const
    {
      return refinement_case_function();
    }

  private:
    /**
     * Lambda function returning whether cell has children.
     */
    const std::function<bool()> has_children_function;

    /**
     * Lambda function returning DoF indices.
     */
    const std::function<void(std::vector<types::global_dof_index> &)>
      get_dof_indices_function;

    /**
     * Lambda function returning active FE index.
     */
    const std::function<unsigned int()> active_fe_index_function;

    /**
     * Lambda function returning the refinement choice.
     */
    const std::function<std::uint8_t()> refinement_case_function;
  };



  /**
   * Base class for a view on fine DoFHandler.
   *
   * Implementations include:
   *  - IdentityFineDoFHandlerView: all cells on the fine mesh are either
   *    locally owned or ghosted; useful for p-multigrid without repartitioning;
   *  - FirstChildPolicyFineDoFHandlerView: parent cells are owned by the first
   *    child cell; useful for local smoothing with fast setup;
   *  - PermutationFineDoFHandlerView: fine mesh has the same cells as the
   * coarse mesh but is partitioned differently; useful for p-multigrid with
   *    repartitioning;
   *  - GlobalCoarseningFineDoFHandlerView: cells on the coarse mesh are either
   *    refined or not; useful for global coarsening.
   */
  template <int dim>
  class FineDoFHandlerViewBase
  {
  public:
    virtual ~FineDoFHandlerViewBase() = default;

    /**
     * Return cell on fine DoFHandler.
     */
    virtual FineDoFHandlerViewCell
    get_cell_view(
      const typename DoFHandler<dim>::cell_iterator &cell) const = 0;

    /**
     * Return child of cell on fine DoFHandler.
     */
    virtual FineDoFHandlerViewCell
    get_cell_view(const typename DoFHandler<dim>::cell_iterator &cell,
                  const unsigned int                             c) const = 0;
  };



  template <int dim>
  class IdentityFineDoFHandlerView : public FineDoFHandlerViewBase<dim>
  {
  public:
    IdentityFineDoFHandlerView(const DoFHandler<dim> &dof_handler_fine,
                               const unsigned int     mg_level_fine)
      : dof_handler_fine(dof_handler_fine)
      , mg_level_fine(mg_level_fine)
    {}

    virtual ~IdentityFineDoFHandlerView() = default;

    FineDoFHandlerViewCell
    get_cell_view(
      const typename DoFHandler<dim>::cell_iterator &cell) const override
    {
      return FineDoFHandlerViewCell(
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        },
        [this, cell](auto &dof_indices) {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            cell->as_dof_handler_iterator(this->dof_handler_fine)
              ->get_dof_indices(dof_indices);
          else
            cell->as_dof_handler_level_iterator(this->dof_handler_fine)
              ->get_mg_dof_indices(dof_indices);
        },
        [this, cell]() {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            return cell->as_dof_handler_iterator(dof_handler_fine)
              ->active_fe_index();
          else
            return cell->as_dof_handler_level_iterator(this->dof_handler_fine)
              ->active_fe_index();
        },
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return 0;
        });
    }

    FineDoFHandlerViewCell
    get_cell_view(const typename DoFHandler<dim>::cell_iterator &,
                  const unsigned int) const override
    {
      DEAL_II_ASSERT_UNREACHABLE();

      return FineDoFHandlerViewCell(
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        },
        [](auto &) { DEAL_II_ASSERT_UNREACHABLE(); },
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return 0;
        },
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return 0;
        });
    }

  private:
    const DoFHandler<dim> &dof_handler_fine;
    const unsigned int     mg_level_fine;
  };



  template <int dim>
  class FirstChildPolicyFineDoFHandlerView : public FineDoFHandlerViewBase<dim>
  {
  public:
    FirstChildPolicyFineDoFHandlerView(const DoFHandler<dim> &dof_handler_fine,
                                       const unsigned int     mg_level_fine)
      : dof_handler_fine(dof_handler_fine)
      , mg_level_fine(mg_level_fine)
    {}

    virtual ~FirstChildPolicyFineDoFHandlerView() = default;

    FineDoFHandlerViewCell
    get_cell_view(
      const typename DoFHandler<dim>::cell_iterator &cell) const override
    {
      return FineDoFHandlerViewCell(
        [this, cell]() {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              // create fine cell in two steps, since the coarse cell and
              // the fine cell are associated to different Triangulation
              // objects
              const auto cell_id = cell->id();
              const auto cell_fine_raw =
                dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell_id);
              return cell_fine_raw->has_children();
            }
          else
            {
              return cell->has_children();
            }
        },
        [this, cell](auto &dof_indices) {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              const auto cell_id = cell->id();
              const auto cell_fine_raw =
                dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell_id);
              return cell_fine_raw->as_dof_handler_iterator(dof_handler_fine)
                ->get_dof_indices(dof_indices);
            }
          else
            {
              cell->get_mg_dof_indices(dof_indices);
            }
        },
        [this, cell]() {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              const auto cell_id = cell->id();
              const auto cell_fine_raw =
                dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell_id);
              return cell_fine_raw->as_dof_handler_iterator(dof_handler_fine)
                ->active_fe_index();
            }
          else
            {
              return cell->active_fe_index();
            }
        },
        []() {
          DEAL_II_ASSERT_UNREACHABLE(); // only children have the correct info
          return 0;
        });
    }

    FineDoFHandlerViewCell
    get_cell_view(const typename DoFHandler<dim>::cell_iterator &cell,
                  const unsigned int c) const override
    {
      return FineDoFHandlerViewCell(
        []() {
          DEAL_II_ASSERT_UNREACHABLE();
          return false;
        },
        [this, cell, c](auto &dof_indices) {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              const auto cell_id       = cell->id();
              const auto cell_fine_raw = dof_handler_fine.get_triangulation()
                                           .create_cell_iterator(cell_id)
                                           ->child(c);
              cell_fine_raw->as_dof_handler_iterator(dof_handler_fine)
                ->get_dof_indices(dof_indices);
            }
          else
            {
              cell->child(c)->get_mg_dof_indices(dof_indices);
            }
        },
        [this, cell, c]() {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              const auto cell_id       = cell->id();
              const auto cell_fine_raw = dof_handler_fine.get_triangulation()
                                           .create_cell_iterator(cell_id)
                                           ->child(c);
              return cell_fine_raw->as_dof_handler_iterator(dof_handler_fine)
                ->active_fe_index();
            }
          else
            {
              return cell->child(c)->active_fe_index();
            }
        },
        [this, cell]() {
          if (this->mg_level_fine == numbers::invalid_unsigned_int)
            {
              const auto cell_id = cell->id();
              const auto cell_raw =
                dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell_id);
              return cell_raw->refinement_case();
            }
          else
            {
              return cell->refinement_case();
            }
        });
    }

  private:
    const DoFHandler<dim> &dof_handler_fine;
    const unsigned int     mg_level_fine;
  };



  template <int dim>
  class BlackBoxFineDoFHandlerView : public FineDoFHandlerViewBase<dim>
  {
  public:
    BlackBoxFineDoFHandlerView(const DoFHandler<dim> &dof_handler_fine,
                               const DoFHandler<dim> &dof_handler_coarse,
                               const unsigned int     mg_level_fine)
      : dof_handler_fine(dof_handler_fine)
      , dof_handler_coarse(dof_handler_coarse)
      , mg_level_fine(mg_level_fine)
      , communicator(
          dof_handler_fine
            .get_mpi_communicator() /*TODO: fix for different comms*/)
      , cell_id_translator(
          dof_handler_fine.get_triangulation().n_global_coarse_cells(),
          dof_handler_fine.get_triangulation().n_global_levels())
    {
      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_coarse.get_triangulation().n_global_coarse_cells());
      AssertIndexRange(dof_handler_coarse.get_triangulation().n_global_levels(),
                       dof_handler_fine.get_triangulation().n_global_levels() +
                         1);
    }

    virtual ~BlackBoxFineDoFHandlerView() = default;

    void
    reinit(const IndexSet &is_dst_locally_owned,
           const IndexSet &is_dst_remote_input,
           const IndexSet &is_src_locally_owned,
           const bool      check_if_elements_in_is_dst_remote_exist = false)
    {
      IndexSet is_dst_remote = is_dst_remote_input;

      if (check_if_elements_in_is_dst_remote_exist)
        {
          IndexSet is_dst_remote_potentially_relevant = is_dst_remote;
          is_dst_remote.clear();

          is_dst_remote_potentially_relevant.subtract_set(is_dst_locally_owned);

          const std::vector<unsigned int> owning_ranks_of_ghosts =
            Utilities::MPI::compute_index_owner(
              is_dst_locally_owned,
              is_dst_remote_potentially_relevant,
              communicator);

          for (unsigned i = 0;
               i < is_dst_remote_potentially_relevant.n_elements();
               ++i)
            if (owning_ranks_of_ghosts[i] != numbers::invalid_unsigned_int)
              is_dst_remote.add_index(
                is_dst_remote_potentially_relevant.nth_index_in_set(i));
        }

      // determine owners of remote cells
      const auto [is_dst_remote_owners, targets_with_indexset] =
        Utilities::MPI::compute_index_owner_and_requesters(is_dst_locally_owned,
                                                           is_dst_remote,
                                                           communicator);


      this->is_dst_locally_owned = is_dst_locally_owned;
      this->is_dst_remote        = is_dst_remote;
      this->is_src_locally_owned = is_src_locally_owned;

#ifndef DEAL_II_WITH_MPI
      Assert(targets_with_indexset.empty(), ExcInternalError());
#else

      std::map<unsigned int, std::vector<types::global_dof_index>>
                               indices_to_be_sent;
      std::vector<MPI_Request> requests;
      requests.reserve(targets_with_indexset.size());
      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(communicator);

      {
        std::vector<types::global_dof_index> indices;

        for (const auto &i : targets_with_indexset)
          {
            // Skip communication in case we would send to ourselves or when
            // there are no indices to send (this can still happen in the run
            // of the consensus algorithms above if the index spaces are
            // sparse).
            if (i.first == my_rank || i.second.begin() == i.second.end())
              continue;

            indices_to_be_sent[i.first] = {};
            std::vector<types::global_dof_index> &buffer =
              indices_to_be_sent[i.first];

            for (auto cell_id : i.second)
              {
                typename DoFHandler<dim>::cell_iterator cell(
                  *dof_handler_fine.get_triangulation().create_cell_iterator(
                    cell_id_translator.to_cell_id(cell_id)),
                  &dof_handler_fine);

                indices.resize(cell->get_fe().n_dofs_per_cell());

                if (mg_level_fine == numbers::invalid_unsigned_int)
                  cell->get_dof_indices(indices);
                else
                  cell->get_mg_dof_indices(indices);

                buffer.push_back(cell->active_fe_index());
                if (cell->level() != 0)
                  buffer.push_back(cell->parent()->refinement_case());
                else
                  buffer.push_back(numbers::invalid_dof_index);
                buffer.insert(buffer.end(), indices.begin(), indices.end());
              }

            requests.emplace_back();

            const auto ierr_1 = MPI_Isend(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              i.first,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &requests.back());
            AssertThrowMPI(ierr_1);
          }
      }

      // process local cells
      {
        std::map<unsigned int, std::vector<types::global_dof_index>>
          rank_to_ids;
        for (unsigned int i = 0; i < is_dst_remote_owners.size(); ++i)
          rank_to_ids[is_dst_remote_owners[i]].push_back(
            is_dst_remote.nth_index_in_set(i));

        for (const auto &pair : rank_to_ids)
          {
            // above we skip messages sent to myself, so also skip MPI_Probe
            if (pair.first == my_rank)
              continue;

            MPI_Status status;
            const int  ierr_1 = MPI_Probe(
              MPI_ANY_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              &status);
            AssertThrowMPI(ierr_1);

            std::vector<types::global_dof_index> buffer;

            int       message_length;
            const int ierr_2 = MPI_Get_count(
              &status,
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              &message_length);
            AssertThrowMPI(ierr_2);

            buffer.resize(message_length);

            const int ierr_3 = MPI_Recv(
              buffer.data(),
              buffer.size(),
              Utilities::MPI::mpi_type_id_for_type<decltype(*buffer.data())>,
              status.MPI_SOURCE,
              Utilities::MPI::internal::Tags::fine_dof_handler_view_reinit,
              communicator,
              MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr_3);

            const unsigned int rank = status.MPI_SOURCE;

            const auto ids = rank_to_ids[rank];

            std::vector<types::global_dof_index> indices;

            for (unsigned int i = 0, k = 0; i < ids.size(); ++i)
              {
                const unsigned int active_fe_index = buffer[k++];
                const std::uint8_t refinement_case =
                  static_cast<std::uint8_t>(buffer[k++]);
                indices.resize(
                  dof_handler_fine.get_fe(active_fe_index).n_dofs_per_cell());

                for (unsigned int j = 0; j < indices.size(); ++j, ++k)
                  indices[j] = buffer[k];
                remote_cell_info_map[ids[i]] = {active_fe_index,
                                                indices,
                                                refinement_case};
              }
          }

        if (!requests.empty())
          {
            const int ierr_1 = MPI_Waitall(requests.size(),
                                           requests.data(),
                                           MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr_1);
          }
      }
#endif
    }

    FineDoFHandlerViewCell
    get_cell_view(
      const typename DoFHandler<dim>::cell_iterator &cell) const override
    {
      const auto id = this->cell_id_translator.translate(cell);

      const bool is_cell_locally_owned =
        this->is_dst_locally_owned.is_element(id);
      const bool is_cell_remotely_owned = this->is_dst_remote.is_element(id);

      const bool has_cell_any_children = [&]() {
        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
             ++i)
          {
            const auto j = this->cell_id_translator.translate(cell, i);

            if (this->is_dst_locally_owned.is_element(j))
              return true;

            if (this->is_dst_remote.is_element(j))
              return true;
          }

        AssertThrow(is_cell_locally_owned || is_cell_remotely_owned,
                    ExcInternalError());

        return false;
      }();

      return FineDoFHandlerViewCell(
        [has_cell_any_children]() { return has_cell_any_children; },
        [cell, is_cell_locally_owned, is_cell_remotely_owned, id, this](
          auto &dof_indices) {
          if (is_cell_locally_owned)
            {
              typename DoFHandler<dim>::cell_iterator cell_fine(
                *dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell->id()),
                &dof_handler_fine);
              if (mg_level_fine == numbers::invalid_unsigned_int)
                cell_fine->get_dof_indices(dof_indices);
              else
                cell_fine->get_mg_dof_indices(dof_indices);
            }
          else if (is_cell_remotely_owned)
            {
              dof_indices = std::get<1>(remote_cell_info_map.at(id));
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
            }
        },
        [cell, is_cell_locally_owned, is_cell_remotely_owned, id, this]()
          -> unsigned int {
          if (is_cell_locally_owned)
            {
              return (typename DoFHandler<dim>::cell_iterator(
                        *dof_handler_fine.get_triangulation()
                           .create_cell_iterator(cell->id()),
                        &dof_handler_fine))
                ->active_fe_index();
            }
          else if (is_cell_remotely_owned)
            {
              return std::get<0>(remote_cell_info_map.at(id));
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
              return 0;
            }
        },
        []() -> std::uint8_t {
          AssertThrow(false,
                      ExcNotImplemented()); // only children hold correct info
          return 0;
        });
    }

    FineDoFHandlerViewCell
    get_cell_view(const typename DoFHandler<dim>::cell_iterator &cell,
                  const unsigned int c) const override
    {
      const auto id = this->cell_id_translator.translate(cell, c);

      const bool is_cell_locally_owned =
        this->is_dst_locally_owned.is_element(id);
      const bool is_cell_remotely_owned = this->is_dst_remote.is_element(id);

      return FineDoFHandlerViewCell(
        [cell]() {
          AssertThrow(false, ExcNotImplemented()); // currently we do not need
                                                   // children of children

          return false;
        },
        [cell, is_cell_locally_owned, is_cell_remotely_owned, c, id, this](
          auto &dof_indices) {
          if (is_cell_locally_owned)
            {
              const auto cell_fine =
                (typename DoFHandler<dim>::cell_iterator(
                   *dof_handler_fine.get_triangulation().create_cell_iterator(
                     cell->id()),
                   &dof_handler_fine))
                  ->child(c);

              if (mg_level_fine == numbers::invalid_unsigned_int)
                cell_fine->get_dof_indices(dof_indices);
              else
                cell_fine->get_mg_dof_indices(dof_indices);
            }
          else if (is_cell_remotely_owned)
            {
              dof_indices = std::get<1>(remote_cell_info_map.at(id));
            }
          else
            {
              AssertThrow(false, ExcNotImplemented()); // should not happen!
            }
        },
        []() -> unsigned int {
          AssertThrow(false, ExcNotImplemented()); // currently we do not need
                                                   // active_fe_index() for
                                                   // children

          return 0;
        },
        [cell, is_cell_locally_owned, is_cell_remotely_owned, id, this]() {
          if (is_cell_locally_owned &&
              (mg_level_fine == numbers::invalid_unsigned_int))
            {
              const auto cell_fine = (typename DoFHandler<dim>::cell_iterator(
                *dof_handler_fine.get_triangulation().create_cell_iterator(
                  cell->id()),
                &dof_handler_fine));
              return cell_fine->refinement_case();
            }
          else if (is_cell_locally_owned)
            {
              return cell->refinement_case();
            }

          else if (is_cell_remotely_owned)
            {
              return RefinementCase<dim>(
                std::get<2>(remote_cell_info_map.at(id)));
            }
          else
            {
              AssertThrow(false, ExcNotImplemented());
            }
          // return no refinement
          return RefinementCase<dim>(0);
        });
    }

  private:
    const DoFHandler<dim> &dof_handler_fine;
    const DoFHandler<dim> &dof_handler_coarse;
    const unsigned int     mg_level_fine;

  protected:
    const MPI_Comm              communicator;
    const CellIDTranslator<dim> cell_id_translator;
    IndexSet                    is_dst_locally_owned;
    IndexSet                    is_dst_remote;

  private:
    IndexSet is_src_locally_owned;

    std::map<types::global_cell_index,
             std::tuple<unsigned int,
                        std::vector<types::global_dof_index>,
                        std::uint8_t>>
      remote_cell_info_map;
  };

  template <int dim>
  class GlobalCoarseningFineDoFHandlerView
    : public BlackBoxFineDoFHandlerView<dim>
  {
  public:
    GlobalCoarseningFineDoFHandlerView(const DoFHandler<dim> &dof_handler_dst,
                                       const DoFHandler<dim> &dof_handler_src,
                                       const unsigned int     mg_level_fine,
                                       const unsigned int     mg_level_coarse)
      : BlackBoxFineDoFHandlerView<dim>(dof_handler_dst,
                                        dof_handler_src,
                                        mg_level_fine)
    {
      Assert((mg_level_fine == numbers::invalid_unsigned_int &&
              mg_level_coarse == numbers::invalid_unsigned_int) ||
               (mg_level_coarse + 1 == mg_level_fine),
             ExcNotImplemented());

      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      const auto fine_operation = [&](const auto &cell) {
        is_dst_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
      };

      const auto coarse_operation = [&](const auto &cell) {
        is_src_locally_owned.add_index(
          this->cell_id_translator.translate(cell));

        // in the case of global coarsening identity transfer is possible
        if (mg_level_coarse == numbers::invalid_unsigned_int)
          is_dst_remote.add_index(this->cell_id_translator.translate(cell));

        if (cell->level() + 1u == tria_dst.n_global_levels())
          return;

        for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell;
             ++i)
          is_dst_remote.add_index(this->cell_id_translator.translate(cell, i));
      };

      loop_over_active_or_level_cells(tria_dst, mg_level_fine, fine_operation);
      loop_over_active_or_level_cells(tria_src,
                                      mg_level_coarse,
                                      coarse_operation);

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   true);

      // check if meshes are compatible
      if (mg_level_coarse == numbers::invalid_unsigned_int)
        {
          std::vector<std::string> not_found_cells_local;

          const auto coarse_operation_check = [&](const auto &cell) {
            bool flag = false;

            const auto index = this->cell_id_translator.translate(cell);

            flag |= this->is_dst_remote.is_element(index) ||
                    this->is_dst_locally_owned.is_element(index);

            if (cell->level() + 1u != tria_dst.n_global_levels())
              {
                for (unsigned int i = 0;
                     i < GeometryInfo<dim>::max_children_per_cell;
                     ++i)
                  {
                    const auto index =
                      this->cell_id_translator.translate(cell, i);

                    flag |= this->is_dst_remote.is_element(index) ||
                            this->is_dst_locally_owned.is_element(index);
                  }
              }

            if (!flag)
              not_found_cells_local.emplace_back(cell->id().to_string());
          };

          loop_over_active_or_level_cells(tria_src,
                                          mg_level_coarse,
                                          coarse_operation_check);

          auto not_found_cells =
            Utilities::MPI::reduce<std::vector<std::string>>(
              not_found_cells_local,
              this->communicator,
              [](const auto &a, const auto &b) {
                auto result = a;
                result.insert(result.end(), b.begin(), b.end());
                return result;
              },
              0);

          if (Utilities::MPI::this_mpi_process(this->communicator) == 0 &&
              !not_found_cells.empty())
            {
              std::sort(not_found_cells.begin(), not_found_cells.end());

              const std::string str =
                boost::algorithm::join(not_found_cells, ", ");

              AssertThrow(
                false,
                ExcMessage(
                  "Problem setting up two-level transfer operator, since coarse triangulation "
                  "seems to be obtainable by simple coarsening. Following coarse cells "
                  "or children cells could not be found in the fine mesh: " +
                  str + "."));
            }
        }
    }

    virtual ~GlobalCoarseningFineDoFHandlerView() = default;
  };



  template <int dim>
  class PermutationFineDoFHandlerView : public BlackBoxFineDoFHandlerView<dim>
  {
  public:
    PermutationFineDoFHandlerView(const DoFHandler<dim> &dof_handler_dst,
                                  const DoFHandler<dim> &dof_handler_src,
                                  const unsigned int     mg_level_fine,
                                  const unsigned int     mg_level_coarse)
      : BlackBoxFineDoFHandlerView<dim>(dof_handler_dst,
                                        dof_handler_src,
                                        mg_level_fine)
    {
      // get reference to triangulations
      const auto &tria_dst = dof_handler_dst.get_triangulation();
      const auto &tria_src = dof_handler_src.get_triangulation();

      // create index sets
      IndexSet is_dst_locally_owned(this->cell_id_translator.size());
      IndexSet is_dst_remote(this->cell_id_translator.size());
      IndexSet is_src_locally_owned(this->cell_id_translator.size());

      const auto fine_operation = [&](const auto &cell) {
        is_dst_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
      };

      const auto coarse_operation = [&](const auto &cell) {
        is_src_locally_owned.add_index(
          this->cell_id_translator.translate(cell));
        is_dst_remote.add_index(this->cell_id_translator.translate(cell));
      };

      loop_over_active_or_level_cells(tria_dst, mg_level_fine, fine_operation);
      loop_over_active_or_level_cells(tria_src,
                                      mg_level_coarse,
                                      coarse_operation);

      this->reinit(is_dst_locally_owned,
                   is_dst_remote,
                   is_src_locally_owned,
                   false);
    }

    virtual ~PermutationFineDoFHandlerView() = default;
  };



  template <int dim, int spacedim>
  bool
  p_transfer_involves_repartitioning(
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const unsigned int               mg_level_fine,
    const unsigned int               mg_level_coarse)
  {
    if (mg_level_fine != mg_level_coarse)
      return true;

    if (&dof_handler_fine.get_triangulation() !=
        &dof_handler_coarse.get_triangulation())
      return true;

    return false;
  }



  template <int dim, int spacedim>
  bool
  h_transfer_uses_first_child_policy(
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const unsigned int               mg_level_fine,
    const unsigned int               mg_level_coarse)
  {
    if (mg_level_fine == numbers::invalid_unsigned_int &&
        mg_level_coarse == numbers::invalid_unsigned_int)
      {
        // two DoFHandlers

        bool flag = true;

        loop_over_active_or_level_cells(
          dof_handler_coarse.get_triangulation(),
          mg_level_coarse,
          [&](const auto &cell) {
            const auto cell_id = cell->id();

            if (dof_handler_fine.get_triangulation().contains_cell(cell_id) ==
                false)
              {
                flag = false;
              }
            else
              {
                const auto cell_fine =
                  dof_handler_fine.get_triangulation().create_cell_iterator(
                    cell_id);

                if (cell_fine->is_active())
                  {
                    if (cell_fine->subdomain_id() != cell->subdomain_id())
                      flag = false;
                  }
                else
                  {
                    Assert(
                      cell_fine->child(0)->is_active(),
                      ExcMessage(
                        "In building a transfer operator, we are "
                        "expecting that a cell is not or at most once refined "
                        "on the other multigrid level. But it is refined more than once. "
                        "Are you trying to build a transfer operator across "
                        "more than one level of mesh refinement?"));
                    if (cell_fine->child(0)->subdomain_id() !=
                        cell->subdomain_id())
                      flag = false;
                  }
              }
          });

        return Utilities::MPI::min(static_cast<unsigned int>(flag),
                                   dof_handler_fine.get_mpi_communicator()) ==
               1;
      }
    else
      {
        // single DoFHandler
        if (mg_level_fine == numbers::invalid_unsigned_int ||
            mg_level_coarse == numbers::invalid_unsigned_int)
          return false;

        if (mg_level_coarse + 1 != mg_level_fine)
          return false;

        if (&dof_handler_fine != &dof_handler_coarse)
          return false;

        return true;
      }
  }



  template <int dim>
  bool
  fe_has_tp_structure(const FiniteElement<dim> &fe)
  {
    return fe.n_base_elements() == 1 &&
           ((dynamic_cast<const FE_Q<dim> *>(&fe.base_element(0)) != nullptr) ||
            (dynamic_cast<const FE_DGQ<dim> *>(&fe.base_element(0)) !=
             nullptr) ||
            (dynamic_cast<const FE_Q_iso_Q1<dim> *>(&fe.base_element(0)) !=
             nullptr));
  }



  class MGTwoLevelTransferImplementation
  {
    /**
     * Compute weights.
     */
    template <int dim, typename Number, typename MemorySpace>
    static void
    setup_weights(
      const dealii::AffineConstraints<Number> &constraints_fine,
      MGTwoLevelTransfer<
        dim,
        LinearAlgebra::distributed::Vector<Number, MemorySpace>> &transfer,
      const bool                                                  is_feq)
    {
      if (transfer.fine_element_is_continuous == false)
        return; // nothing to do

      // 1) compute weights globally
      LinearAlgebra::distributed::Vector<Number> weight_vector;

      weight_vector.reinit(transfer.partitioner_fine);

      // ... compute valence of DoFs
      for (const auto i : transfer.constraint_info_fine.dof_indices)
        weight_vector.local_element(i) += 1;
      weight_vector.compress(VectorOperation::add);

      // ... invert valence
      for (unsigned int i = 0; i < weight_vector.locally_owned_size(); ++i)
        if (weight_vector.local_element(i) > 0.)
          weight_vector.local_element(i) =
            Number(1.) / weight_vector.local_element(i);

      // ... clear constrained indices
      for (const auto &i : constraints_fine.get_lines())
        if (weight_vector.locally_owned_elements().is_element(i.index))
          weight_vector[i.index] = 0.0;

      weight_vector.update_ghost_values();

      // 2) store data cell-wise a DG format and try to compress
      const unsigned int n_lanes   = VectorizedArray<Number>::size();
      std::size_t        n_batches = 0;
      for (const auto &scheme : transfer.schemes)
        n_batches += (scheme.n_coarse_cells + n_lanes - 1) / n_lanes;

      transfer.weights.clear();
      transfer.weights.reserve(n_batches * Utilities::pow(3, dim));
      transfer.weights_start.clear();
      transfer.weights_start.reserve(n_batches + 1);
      transfer.weights_are_compressed.clear();
      transfer.weights_are_compressed.resize(n_batches);

      AlignedVector<VectorizedArray<Number>> local_weights;
      std::array<VectorizedArray<Number>, Utilities::pow(3, dim)>
                   local_weights_compressed;
      unsigned int offset = 0;

      // ... loop over cells
      for (const auto &scheme : transfer.schemes)
        {
          local_weights.resize(scheme.n_dofs_per_cell_fine);
          for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
               cell += n_lanes)
            {
              transfer.weights_start.push_back(transfer.weights.size());

              const unsigned int n_lanes_filled =
                (cell + n_lanes > scheme.n_coarse_cells) ?
                  (scheme.n_coarse_cells - cell) :
                  n_lanes;

              for (unsigned int v = 0; v < n_lanes_filled;
                   ++v, offset += scheme.n_dofs_per_cell_fine)
                for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine; ++i)
                  local_weights[i][v] = weight_vector.local_element(
                    transfer.constraint_info_fine.dof_indices[offset + i]);

              const bool can_compress_weights =
                is_feq && compute_weights_fe_q_dofs_by_entity<dim, -1>(
                            local_weights.data(),
                            transfer.n_components,
                            scheme.degree_fine + 1,
                            local_weights_compressed.data());

              if (can_compress_weights && scheme.degree_fine > 1)
                {
                  for (const VectorizedArray<Number> weight :
                       local_weights_compressed)
                    transfer.weights.push_back(weight);
                  transfer
                    .weights_are_compressed[transfer.weights_start.size() - 1] =
                    1;
                }
              else
                for (const VectorizedArray<Number> weight : local_weights)
                  transfer.weights.push_back(weight);
            }
        }
      AssertDimension(transfer.weights_start.size(), n_batches);
    }



  public:
    template <int dim, typename Number>
    static std::shared_ptr<const Utilities::MPI::Partitioner>
    create_coarse_partitioner(
      const DoFHandler<dim>                   &dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_coarse)
    {
      IndexSet locally_relevant_dofs =
        (mg_level_coarse == numbers::invalid_unsigned_int) ?
          DoFTools::extract_locally_active_dofs(dof_handler_coarse) :
          DoFTools::extract_locally_active_level_dofs(dof_handler_coarse,
                                                      mg_level_coarse);

      std::vector<types::global_dof_index> locally_relevant_dofs_temp;

      for (const auto i : locally_relevant_dofs)
        {
          if (locally_relevant_dofs.is_element(i) == false)
            locally_relevant_dofs_temp.emplace_back(i);

          const auto constraints = constraints_coarse.get_constraint_entries(i);

          if (constraints)
            for (const auto &p : *constraints)
              if (locally_relevant_dofs.is_element(p.first) == false)
                locally_relevant_dofs_temp.emplace_back(p.first);
        }

      std::sort(locally_relevant_dofs_temp.begin(),
                locally_relevant_dofs_temp.end());
      locally_relevant_dofs.add_indices(locally_relevant_dofs_temp.begin(),
                                        locally_relevant_dofs_temp.end());

      return std::make_shared<Utilities::MPI::Partitioner>(
        mg_level_coarse == numbers::invalid_unsigned_int ?
          dof_handler_coarse.locally_owned_dofs() :
          dof_handler_coarse.locally_owned_mg_dofs(mg_level_coarse),
        locally_relevant_dofs,
        dof_handler_coarse.get_mpi_communicator());
    }



    template <int dim, typename Number, typename MemorySpace>
    static void
    reinit_geometric_transfer(
      const DoFHandler<dim>                   &dof_handler_fine,
      const DoFHandler<dim>                   &dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_fine,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_fine,
      const unsigned int                       mg_level_coarse,
      MGTwoLevelTransfer<
        dim,
        LinearAlgebra::distributed::Vector<Number, MemorySpace>> &transfer)
    {
      Assert((mg_level_fine == numbers::invalid_unsigned_int &&
              mg_level_coarse == numbers::invalid_unsigned_int) ||
               (mg_level_coarse + 1 == mg_level_fine),
             ExcNotImplemented());

      AssertDimension(constraints_fine.n_inhomogeneities(), 0);
      AssertDimension(constraints_coarse.n_inhomogeneities(), 0);

      transfer.dof_handler_fine = &dof_handler_fine;
      transfer.mg_level_fine    = mg_level_fine;

      std::unique_ptr<FineDoFHandlerViewBase<dim>> dof_handler_fine_view;

      if (internal::h_transfer_uses_first_child_policy(dof_handler_fine,
                                                       dof_handler_coarse,
                                                       mg_level_fine,
                                                       mg_level_coarse))
        dof_handler_fine_view =
          std::make_unique<FirstChildPolicyFineDoFHandlerView<dim>>(
            dof_handler_fine, mg_level_fine);
      else
        dof_handler_fine_view =
          std::make_unique<GlobalCoarseningFineDoFHandlerView<dim>>(
            dof_handler_fine,
            dof_handler_coarse,
            mg_level_fine,
            mg_level_coarse);

      // gather ranges for active FE indices on both fine and coarse dofhandlers
      std::array<unsigned int, 2> min_active_fe_indices = {
        {std::numeric_limits<unsigned int>::max(),
         std::numeric_limits<unsigned int>::max()}};
      std::array<unsigned int, 2> max_active_fe_indices = {{0, 0}};

      loop_over_active_or_level_cells(
        dof_handler_fine, mg_level_fine, [&](const auto &cell) {
          min_active_fe_indices[0] =
            std::min<unsigned int>(min_active_fe_indices[0],
                                   cell->active_fe_index());
          max_active_fe_indices[0] =
            std::max<unsigned int>(max_active_fe_indices[0],
                                   cell->active_fe_index());
        });

      loop_over_active_or_level_cells(
        dof_handler_coarse, mg_level_coarse, [&](const auto &cell) {
          min_active_fe_indices[1] =
            std::min<unsigned int>(min_active_fe_indices[1],
                                   cell->active_fe_index());
          max_active_fe_indices[1] =
            std::max<unsigned int>(max_active_fe_indices[1],
                                   cell->active_fe_index());
        });

      const auto comm = dof_handler_fine.get_mpi_communicator();

      Assert(comm == dof_handler_coarse.get_mpi_communicator(),
             ExcNotImplemented());

      ArrayView<unsigned int> temp_min(min_active_fe_indices);
      ArrayView<unsigned int> temp_max(max_active_fe_indices);
      Utilities::MPI::min(temp_min, comm, temp_min);
      Utilities::MPI::max(temp_max, comm, temp_max);

      // make sure that hp is used neither on the coarse nor on the fine
      // dofhandler
      AssertDimension(min_active_fe_indices[0], max_active_fe_indices[0]);
      AssertDimension(min_active_fe_indices[1], max_active_fe_indices[1]);

      const auto reference_cell = dof_handler_fine.get_fe(0).reference_cell();
      // set up mg-schemes
      //   (0) no refinement -> identity
      //   (1) h-refinement
      //   (2) h-refinement choice II <- e.g. for Tets
      //    .
      //    .
      //    .
      transfer.schemes.resize(1 +
                              reference_cell.n_isotropic_refinement_choices());

      const unsigned int fe_index_fine   = min_active_fe_indices[0];
      const unsigned int fe_index_coarse = min_active_fe_indices[1];

      const auto &fe_fine   = dof_handler_fine.get_fe(fe_index_fine);
      const auto &fe_coarse = dof_handler_coarse.get_fe(fe_index_coarse);

      // extract number of components
      AssertDimension(fe_fine.n_components(), fe_coarse.n_components());

      transfer.n_components = fe_fine.n_components();

      // helper function: to process the fine level cells; function @p fu_non_refined is
      // performed on cells that are not refined and @fu_refined is performed on
      // children of cells that are refined
      const auto process_cells = [&](const auto &fu_non_refined,
                                     const auto &fu_refined) {
        loop_over_active_or_level_cells(
          dof_handler_coarse, mg_level_coarse, [&](const auto &cell_coarse) {
            if (mg_level_coarse == numbers::invalid_unsigned_int)
              {
                // get a reference to the equivalent cell on the fine
                // triangulation
                const auto cell_coarse_on_fine_mesh =
                  dof_handler_fine_view->get_cell_view(cell_coarse);

                // check if cell has children
                if (cell_coarse_on_fine_mesh.has_children())
                  // ... cell has children -> process children
                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       c++)
                    fu_refined(cell_coarse,
                               dof_handler_fine_view->get_cell_view(cell_coarse,
                                                                    c),
                               c);
                else // ... cell has no children -> process cell
                  fu_non_refined(cell_coarse, cell_coarse_on_fine_mesh);
              }
            else
              {
                // check if cell has children
                if (cell_coarse->has_children())
                  // ... cell has children -> process children
                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       c++)
                    fu_refined(cell_coarse,
                               dof_handler_fine_view->get_cell_view(cell_coarse,
                                                                    c),
                               c);
              }
          });
      };

      // check if FE is the same
      AssertDimension(fe_coarse.n_dofs_per_cell(), fe_fine.n_dofs_per_cell());


      const bool is_feq = fe_fine.n_base_elements() == 1 &&
                          ((dynamic_cast<const FE_Q<dim> *>(
                              &fe_fine.base_element(0)) != nullptr));

      const bool has_tp_structure = fe_has_tp_structure(fe_fine);

      for (auto &scheme : transfer.schemes)
        {
          // number of dofs on coarse and fine cells
          scheme.n_dofs_per_cell_coarse = fe_coarse.n_dofs_per_cell();
          scheme.n_dofs_per_cell_fine =
            is_feq ? (fe_fine.n_components() *
                      Utilities::pow(2 * fe_fine.degree + 1, dim)) :
                     (fe_coarse.n_dofs_per_cell() *
                      GeometryInfo<dim>::max_children_per_cell);

          // degree of FE on coarse and fine cell
          scheme.degree_coarse = fe_coarse.degree;
          scheme.degree_fine =
            is_feq ? (fe_coarse.degree * 2) : (fe_coarse.degree * 2 + 1);

          // reset number of coarse cells
          scheme.n_coarse_cells = 0;
        }
      // correct for first scheme
      transfer.schemes[0].n_dofs_per_cell_fine = fe_coarse.n_dofs_per_cell();
      transfer.schemes[0].degree_fine          = fe_coarse.degree;

      // continuous or discontinuous
      transfer.fine_element_is_continuous  = fe_fine.n_dofs_per_vertex() > 0;
      std::uint8_t current_refinement_case = static_cast<std::uint8_t>(-1);

      // count coarse cells for each scheme (0, 1, ...)
      {
        // count by looping over all coarse cells
        process_cells(
          [&](const auto &, const auto &) {
            transfer.schemes[0].n_coarse_cells++;
          },
          [&](const auto &, const auto &cell_fine, const auto c) {
            std::uint8_t refinement_case = cell_fine.refinement_case();
            // Assert triggers if cell has no children
            if (reference_cell == ReferenceCells::Tetrahedron)
              Assert(RefinementCase<dim>(refinement_case) ==
                         RefinementCase<dim>(static_cast<std::uint8_t>(
                           IsotropicRefinementChoice::cut_tet_68)) ||
                       RefinementCase<dim>(refinement_case) ==
                         RefinementCase<dim>(static_cast<std::uint8_t>(
                           IsotropicRefinementChoice::cut_tet_57)) ||
                       RefinementCase<dim>(refinement_case) ==
                         RefinementCase<dim>(static_cast<std::uint8_t>(
                           IsotropicRefinementChoice::cut_tet_49)),
                     ExcNotImplemented());
            else
              {
                Assert(RefinementCase<dim>(refinement_case) ==
                         RefinementCase<dim>::isotropic_refinement,
                       ExcNotImplemented());
                refinement_case = 1;
              }

            if (c == 0)
              {
                transfer.schemes[refinement_case].n_coarse_cells++;
                current_refinement_case = refinement_case;
              }
            else
              // Check that all children have the same refinement case
              AssertThrow(current_refinement_case == refinement_case,
                          ExcNotImplemented());
          });
      }


      const auto cell_local_children_indices =
        (has_tp_structure) ?
          get_child_offsets<dim>(transfer.schemes[0].n_dofs_per_cell_coarse,
                                 is_feq ? fe_fine.degree : (fe_fine.degree + 1),
                                 fe_fine.degree) :
          get_child_offsets_general<dim>(
            transfer.schemes[0].n_dofs_per_cell_coarse);

      std::vector<unsigned int> n_dof_indices_fine(transfer.schemes.size() + 1);
      std::vector<unsigned int> n_dof_indices_coarse(transfer.schemes.size() +
                                                     1);

      for (unsigned int i = 0; i < transfer.schemes.size(); ++i)
        {
          n_dof_indices_fine[i + 1] = transfer.schemes[i].n_dofs_per_cell_fine *
                                      transfer.schemes[i].n_coarse_cells;
          n_dof_indices_coarse[i + 1] =
            transfer.schemes[i].n_dofs_per_cell_coarse *
            transfer.schemes[i].n_coarse_cells;
        }

      for (unsigned int i = 0; i < transfer.schemes.size(); ++i)
        {
          n_dof_indices_fine[i + 1] += n_dof_indices_fine[i];
          n_dof_indices_coarse[i + 1] += n_dof_indices_coarse[i];
        }

      // indices
      {
        std::vector<types::global_dof_index> local_dof_indices(
          transfer.schemes[0].n_dofs_per_cell_coarse);

        // ---------------------- lexicographic_numbering ----------------------
        std::vector<unsigned int> lexicographic_numbering_fine;
        std::vector<unsigned int> lexicographic_numbering_coarse;
        if (has_tp_structure)
          {
            const Quadrature<1> dummy_quadrature(
              std::vector<Point<1>>(1, Point<1>()));
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering_fine = shape_info.lexicographic_numbering;
            shape_info.reinit(dummy_quadrature, fe_coarse, 0);
            lexicographic_numbering_coarse = shape_info.lexicographic_numbering;
          }
        else
          {
            const auto dummy_quadrature =
              reference_cell.template get_gauss_type_quadrature<dim>(1);
            internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
            shape_info.reinit(dummy_quadrature, fe_fine, 0);
            lexicographic_numbering_fine = shape_info.lexicographic_numbering;
            shape_info.reinit(dummy_quadrature, fe_coarse, 0);
            lexicographic_numbering_coarse = shape_info.lexicographic_numbering;
          }

        // ------------------------------ indices ------------------------------
        std::vector<types::global_dof_index> level_dof_indices_coarse(
          transfer.schemes[0].n_dofs_per_cell_fine);
        std::vector<types::global_dof_index> level_dof_indices_fine(
          transfer.schemes[1].n_dofs_per_cell_fine);


        unsigned int n_coarse_cells_total = 0;
        for (const auto &scheme : transfer.schemes)
          n_coarse_cells_total += scheme.n_coarse_cells;

        transfer.constraint_info_coarse.reinit(
          dof_handler_coarse,
          n_coarse_cells_total,
          constraints_coarse.n_constraints() > 0 &&
            use_fast_hanging_node_algorithm(dof_handler_coarse,
                                            mg_level_coarse));
        transfer.constraint_info_coarse.set_locally_owned_indices(
          (mg_level_coarse == numbers::invalid_unsigned_int) ?
            dof_handler_coarse.locally_owned_dofs() :
            dof_handler_coarse.locally_owned_mg_dofs(mg_level_coarse));

        transfer.constraint_info_fine.reinit(n_coarse_cells_total);
        transfer.constraint_info_fine.set_locally_owned_indices(
          (mg_level_fine == numbers::invalid_unsigned_int) ?
            dof_handler_fine.locally_owned_dofs() :
            dof_handler_fine.locally_owned_mg_dofs(mg_level_fine));


        std::vector<unsigned int> cell_no(transfer.schemes.size(), 0);
        for (unsigned int i = 1; i < transfer.schemes.size(); ++i)
          cell_no[i] = cell_no[i - 1] + transfer.schemes[i - 1].n_coarse_cells;

        process_cells(
          [&](const auto &cell_coarse, const auto &cell_fine) {
            // first process cells with scheme 0
            // parent
            {
              transfer.constraint_info_coarse.read_dof_indices(
                cell_no[0],
                mg_level_coarse,
                cell_coarse,
                constraints_coarse,
                {});
            }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[0].n_dofs_per_cell_coarse;
                   i++)
                level_dof_indices_coarse[i] =
                  local_dof_indices[lexicographic_numbering_fine[i]];

              transfer.constraint_info_fine.read_dof_indices(
                cell_no[0], level_dof_indices_coarse, {});
            }

            // move pointers
            {
              ++cell_no[0];
            }
          },
          [&](const auto &cell_coarse, const auto &cell_fine, const auto c) {
            // process rest of cells
            const std::uint8_t refinement_case =
              reference_cell == ReferenceCells::Tetrahedron ?
                cell_fine.refinement_case() :
                1;
            // parent (only once at the beginning)
            if (c == 0)
              {
                transfer.constraint_info_coarse.read_dof_indices(
                  cell_no[refinement_case],
                  mg_level_coarse,
                  cell_coarse,
                  constraints_coarse,
                  {});

                level_dof_indices_fine.assign(level_dof_indices_fine.size(),
                                              numbers::invalid_dof_index);
              }

            // child
            {
              cell_fine.get_dof_indices(local_dof_indices);
              for (unsigned int i = 0;
                   i < transfer.schemes[refinement_case].n_dofs_per_cell_coarse;
                   ++i)
                {
                  const auto index =
                    local_dof_indices[lexicographic_numbering_fine[i]];
                  Assert(
                    level_dof_indices_fine[cell_local_children_indices[c][i]] ==
                        numbers::invalid_dof_index ||
                      level_dof_indices_fine[cell_local_children_indices[c]
                                                                        [i]] ==
                        index,
                    ExcInternalError());

                  level_dof_indices_fine[cell_local_children_indices[c][i]] =
                    index;
                }
            }

            // move pointers (only once at the end)
            if (c + 1 == GeometryInfo<dim>::max_children_per_cell)
              {
                transfer.constraint_info_fine.read_dof_indices(
                  cell_no[refinement_case], level_dof_indices_fine, {});

                ++cell_no[refinement_case];
              }
          });
      }

      {
        transfer.partitioner_coarse = transfer.constraint_info_coarse.finalize(
          dof_handler_coarse.get_mpi_communicator());
        transfer.vec_coarse.reinit(transfer.partitioner_coarse);

        transfer.partitioner_fine = transfer.constraint_info_fine.finalize(
          dof_handler_fine.get_mpi_communicator());
        transfer.vec_fine.reinit(transfer.partitioner_fine);
      }


      // ------------- prolongation matrix (0) -> identity matrix --------------

      // nothing to do since for identity prolongation matrices a short-cut
      // code path is used during prolongation/restriction

      // -------------------prolongation matrix (i = 1 ... n)-------------------
      {
        AssertDimension(fe_fine.n_base_elements(), 1);
        for (unsigned int transfer_scheme_index = 1;
             transfer_scheme_index < transfer.schemes.size();
             ++transfer_scheme_index)
          {
            if (has_tp_structure)
              {
                const auto fe = create_1D_fe(fe_fine.base_element(0));

                std::vector<unsigned int> renumbering(fe->n_dofs_per_cell());
                {
                  AssertIndexRange(fe->n_dofs_per_vertex(), 2);
                  renumbering[0] = 0;
                  for (unsigned int i = 0; i < fe->dofs_per_line; ++i)
                    renumbering[i + fe->n_dofs_per_vertex()] =
                      GeometryInfo<1>::vertices_per_cell *
                        fe->n_dofs_per_vertex() +
                      i;
                  if (fe->n_dofs_per_vertex() > 0)
                    renumbering[fe->n_dofs_per_cell() -
                                fe->n_dofs_per_vertex()] =
                      fe->n_dofs_per_vertex();
                }

                // TODO: data structures are saved in form of DG data structures
                // here
                const unsigned int shift =
                  is_feq ? (fe->n_dofs_per_cell() - fe->n_dofs_per_vertex()) :
                           (fe->n_dofs_per_cell());
                const unsigned int n_child_dofs_1d =
                  is_feq ?
                    (fe->n_dofs_per_cell() * 2 - fe->n_dofs_per_vertex()) :
                    (fe->n_dofs_per_cell() * 2);

                {
                  transfer.schemes[transfer_scheme_index]
                    .prolongation_matrix.resize(fe->n_dofs_per_cell() *
                                                n_child_dofs_1d);

                  for (unsigned int c = 0;
                       c < GeometryInfo<1>::max_children_per_cell;
                       ++c)
                    for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
                      for (unsigned int j = 0; j < fe->n_dofs_per_cell(); ++j)
                        transfer.schemes[transfer_scheme_index]
                          .prolongation_matrix[i * n_child_dofs_1d + j +
                                               c * shift] =
                          fe->get_prolongation_matrix(c)(renumbering[j],
                                                         renumbering[i]);
                }
                {
                  transfer.schemes[transfer_scheme_index]
                    .restriction_matrix.resize(fe->n_dofs_per_cell() *
                                               n_child_dofs_1d);

                  for (unsigned int c = 0;
                       c < GeometryInfo<1>::max_children_per_cell;
                       ++c)
                    {
                      const auto matrix = get_restriction_matrix(*fe, c);
                      for (unsigned int i = 0; i < fe->n_dofs_per_cell(); ++i)
                        for (unsigned int j = 0; j < fe->n_dofs_per_cell(); ++j)
                          transfer.schemes[transfer_scheme_index]
                            .restriction_matrix[i * n_child_dofs_1d + j +
                                                c * shift] +=
                            matrix(renumbering[i], renumbering[j]);
                    }
                }
              }
            else
              {
                const auto        &fe              = fe_fine.base_element(0);
                const unsigned int n_dofs_per_cell = fe.n_dofs_per_cell();

                {
                  transfer.schemes[transfer_scheme_index]
                    .prolongation_matrix.resize(
                      n_dofs_per_cell * n_dofs_per_cell *
                      GeometryInfo<dim>::max_children_per_cell);

                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       ++c)
                    {
                      const auto matrix =
                        reference_cell == ReferenceCells::Tetrahedron ?
                          fe.get_prolongation_matrix(
                            c, RefinementCase<dim>(transfer_scheme_index)) :
                          fe.get_prolongation_matrix(c);


                      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                          transfer.schemes[transfer_scheme_index]
                            .prolongation_matrix
                              [i * n_dofs_per_cell *
                                 GeometryInfo<dim>::max_children_per_cell +
                               j + c * n_dofs_per_cell] = matrix(j, i);
                    }
                }
                {
                  transfer.schemes[transfer_scheme_index]
                    .restriction_matrix.resize(
                      n_dofs_per_cell * n_dofs_per_cell *
                      GeometryInfo<dim>::max_children_per_cell);

                  for (unsigned int c = 0;
                       c < GeometryInfo<dim>::max_children_per_cell;
                       ++c)
                    {
                      const auto matrix =
                        reference_cell == ReferenceCells::Tetrahedron ?
                          get_restriction_matrix(
                            fe, c, RefinementCase<dim>(transfer_scheme_index)) :
                          get_restriction_matrix(fe, c);
                      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                          transfer.schemes[transfer_scheme_index]
                            .restriction_matrix
                              [i * n_dofs_per_cell *
                                 GeometryInfo<dim>::max_children_per_cell +
                               j + c * n_dofs_per_cell] += matrix(i, j);
                    }
                }
              }
          }
      }

      // ------------------------------- weights -------------------------------
      if (transfer.fine_element_is_continuous)
        setup_weights(constraints_fine, transfer, is_feq);
    }



    template <int dim, typename Number>
    static void
    compute_prolongation_and_restriction_matrices(
      const FiniteElement<dim> &fe_fine,
      const FiniteElement<dim> &fe_coarse,
      typename MGTwoLevelTransfer<
        dim,
        LinearAlgebra::distributed::Vector<Number>>::MGTransferScheme &scheme)
    {
      AssertDimension(fe_fine.n_base_elements(), 1);
      AssertDimension(fe_coarse.n_base_elements(), 1);

      const bool has_tp_structure = fe_has_tp_structure(fe_fine);

      const FiniteElement<dim> &fe_fine_scalar   = fe_fine.base_element(0);
      const FiniteElement<dim> &fe_coarse_scalar = fe_coarse.base_element(0);

      Assert(fe_fine.reference_cell() == fe_coarse.reference_cell(),
             ExcNotImplemented());

      if (has_tp_structure && (fe_coarse != fe_fine) &&
          (fe_coarse.n_dofs_per_cell() != 0 && fe_fine.n_dofs_per_cell() != 0))
        {
          const auto fe_fine_1d = create_1D_fe(fe_fine_scalar);

          std::vector<unsigned int> renumbering_fine(
            fe_fine_1d->n_dofs_per_cell());
          AssertIndexRange(fe_fine_1d->n_dofs_per_vertex(), 2);
          renumbering_fine[0] = 0;
          for (unsigned int i = 0; i < fe_fine_1d->dofs_per_line; ++i)
            renumbering_fine[i + fe_fine_1d->n_dofs_per_vertex()] =
              GeometryInfo<1>::vertices_per_cell *
                fe_fine_1d->n_dofs_per_vertex() +
              i;
          if (fe_fine_1d->n_dofs_per_vertex() > 0)
            renumbering_fine[fe_fine_1d->n_dofs_per_cell() -
                             fe_fine_1d->n_dofs_per_vertex()] =
              fe_fine_1d->n_dofs_per_vertex();

          const auto fe_coarse_1d = create_1D_fe(fe_coarse_scalar);

          std::vector<unsigned int> renumbering_coarse(
            fe_coarse_1d->n_dofs_per_cell());
          AssertIndexRange(fe_coarse_1d->n_dofs_per_vertex(), 2);
          renumbering_coarse[0] = 0;
          for (unsigned int i = 0; i < fe_coarse_1d->dofs_per_line; ++i)
            renumbering_coarse[i + fe_coarse_1d->n_dofs_per_vertex()] =
              GeometryInfo<1>::vertices_per_cell *
                fe_coarse_1d->n_dofs_per_vertex() +
              i;
          if (fe_coarse_1d->n_dofs_per_vertex() > 0)
            renumbering_coarse[fe_coarse_1d->n_dofs_per_cell() -
                               fe_coarse_1d->n_dofs_per_vertex()] =
              fe_coarse_1d->n_dofs_per_vertex();

          FullMatrix<double> matrix(fe_fine_1d->n_dofs_per_cell(),
                                    fe_coarse_1d->n_dofs_per_cell());
          FETools::get_projection_matrix(*fe_coarse_1d, *fe_fine_1d, matrix);
          scheme.prolongation_matrix.resize(fe_fine_1d->n_dofs_per_cell() *
                                            fe_coarse_1d->n_dofs_per_cell());

          for (unsigned int i = 0, k = 0; i < fe_coarse_1d->n_dofs_per_cell();
               ++i)
            for (unsigned int j = 0; j < fe_fine_1d->n_dofs_per_cell();
                 ++j, ++k)
              scheme.prolongation_matrix[k] =
                matrix(renumbering_fine[j], renumbering_coarse[i]);

          matrix.reinit(fe_coarse_1d->n_dofs_per_cell(),
                        fe_fine_1d->n_dofs_per_cell());
          FETools::get_projection_matrix(*fe_fine_1d, *fe_coarse_1d, matrix);
          scheme.restriction_matrix.resize(fe_fine_1d->n_dofs_per_cell() *
                                           fe_coarse_1d->n_dofs_per_cell());

          for (unsigned int i = 0, k = 0; i < fe_coarse_1d->n_dofs_per_cell();
               ++i)
            for (unsigned int j = 0; j < fe_fine_1d->n_dofs_per_cell();
                 ++j, ++k)
              scheme.restriction_matrix[k] =
                matrix(renumbering_coarse[i], renumbering_fine[j]);
        }
      else if ((!has_tp_structure) && (fe_fine != fe_coarse) &&
               (fe_fine.n_dofs_per_cell() != 0 &&
                fe_coarse.n_dofs_per_cell() != 0))
        {
          FullMatrix<double> matrix(fe_fine_scalar.n_dofs_per_cell(),
                                    fe_coarse_scalar.n_dofs_per_cell());
          FETools::get_projection_matrix(fe_coarse_scalar,
                                         fe_fine_scalar,
                                         matrix);
          scheme.prolongation_matrix.resize(matrix.m() * matrix.n());

          for (unsigned int i = 0, k = 0; i < matrix.n(); ++i)
            for (unsigned int j = 0; j < matrix.m(); ++j, ++k)
              scheme.prolongation_matrix[k] = matrix(j, i);

          matrix.reinit(fe_coarse_scalar.n_dofs_per_cell(),
                        fe_fine_scalar.n_dofs_per_cell());
          FETools::get_projection_matrix(fe_fine_scalar,
                                         fe_coarse_scalar,
                                         matrix);
          scheme.restriction_matrix.resize(matrix.m() * matrix.n());

          for (unsigned int i = 0, k = 0; i < matrix.m(); ++i)
            for (unsigned int j = 0; j < matrix.n(); ++j, ++k)
              scheme.restriction_matrix[k] = matrix(i, j);
        }
    }



    template <int dim, typename Number>
    static void
    reinit_polynomial_transfer(
      const DoFHandler<dim>                   &dof_handler_fine,
      const DoFHandler<dim>                   &dof_handler_coarse,
      const dealii::AffineConstraints<Number> &constraints_fine,
      const dealii::AffineConstraints<Number> &constraints_coarse,
      const unsigned int                       mg_level_fine,
      const unsigned int                       mg_level_coarse,
      MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
        &transfer)
    {
      Assert(mg_level_fine == numbers::invalid_unsigned_int ||
               mg_level_fine <= MGTools::max_level_for_coarse_mesh(
                                  dof_handler_fine.get_triangulation()),
             ExcMessage("Polynomial transfer is only allowed on the active "
                        "level (numbers::invalid_unsigned_int) or on "
                        "refinement levels without hanging nodes."));
      Assert(mg_level_coarse == numbers::invalid_unsigned_int ||
               mg_level_coarse <= MGTools::max_level_for_coarse_mesh(
                                    dof_handler_coarse.get_triangulation()),
             ExcMessage("Polynomial transfer is only allowed on the active "
                        "level (numbers::invalid_unsigned_int) or on "
                        "refinement levels without hanging nodes."));

      AssertDimension(constraints_fine.n_inhomogeneities(), 0);
      AssertDimension(constraints_coarse.n_inhomogeneities(), 0);

      transfer.dof_handler_fine = &dof_handler_fine;
      transfer.mg_level_fine    = mg_level_fine;

      std::unique_ptr<FineDoFHandlerViewBase<dim>> dof_handler_fine_view;

      if (internal::p_transfer_involves_repartitioning(dof_handler_fine,
                                                       dof_handler_coarse,
                                                       mg_level_fine,
                                                       mg_level_coarse))
        dof_handler_fine_view =
          std::make_unique<PermutationFineDoFHandlerView<dim>>(
            dof_handler_fine,
            dof_handler_coarse,
            mg_level_fine,
            mg_level_coarse);
      else
        dof_handler_fine_view =
          std::make_unique<IdentityFineDoFHandlerView<dim>>(dof_handler_fine,
                                                            mg_level_fine);

      // TODO: adjust assert
      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_active_cells(),
        dof_handler_coarse.get_triangulation().n_global_active_cells());

      // extract number of components
      AssertDimension(dof_handler_fine.get_fe_collection().n_components(),
                      dof_handler_coarse.get_fe_collection().n_components());

      transfer.n_components =
        dof_handler_fine.get_fe_collection().n_components();

      transfer.fine_element_is_continuous =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) {
                      return fe.n_dofs_per_cell() == 0 ||
                             fe.n_dofs_per_vertex() > 0;
                    });

      if constexpr (running_in_debug_mode())
        {
          const bool fine_element_is_discontinuous =
            std::all_of(dof_handler_fine.get_fe_collection().begin(),
                        dof_handler_fine.get_fe_collection().end(),
                        [](const auto &fe) {
                          return fe.n_dofs_per_cell() == 0 ||
                                 fe.n_dofs_per_vertex() == 0;
                        });

          Assert(transfer.fine_element_is_continuous !=
                   fine_element_is_discontinuous,
                 ExcNotImplemented());
        }
      const bool is_feq =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) {
                      return fe.n_base_elements() == 1 &&
                             (dynamic_cast<const FE_Q<dim> *>(
                                &fe.base_element(0)) != nullptr);
                    });

      const bool has_tp_structure =
        std::all_of(dof_handler_fine.get_fe_collection().begin(),
                    dof_handler_fine.get_fe_collection().end(),
                    [](const auto &fe) { return fe_has_tp_structure(fe); });

      const auto process_cells = [&](const auto &fu) {
        loop_over_active_or_level_cells(
          dof_handler_coarse, mg_level_coarse, [&](const auto &cell_coarse) {
            const auto cell_coarse_on_fine_mesh =
              dof_handler_fine_view->get_cell_view(cell_coarse);
            fu(cell_coarse, &cell_coarse_on_fine_mesh);
          });
      };

      std::map<std::pair<unsigned int, unsigned int>, unsigned int>
        fe_index_pairs;

      process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
        fe_index_pairs.emplace(
          std::pair<unsigned int, unsigned int>(cell_coarse->active_fe_index(),
                                                cell_fine->active_fe_index()),
          0);
      });

      unsigned int counter = 0;
      for (auto &f : fe_index_pairs)
        f.second = counter++;

      transfer.schemes.resize(counter);

      // extract number of coarse cells
      {
        for (auto &scheme : transfer.schemes)
          scheme.n_coarse_cells = 0;
        process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
          transfer
            .schemes[fe_index_pairs[std::pair<unsigned int, unsigned int>(
              cell_coarse->active_fe_index(), cell_fine->active_fe_index())]]
            .n_coarse_cells++;
        });
      }

      for (const auto &fe_index_pair : fe_index_pairs)
        {
          transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first)
              .n_dofs_per_cell();
          transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second)
              .n_dofs_per_cell();

          transfer.schemes[fe_index_pair.second].degree_coarse =
            dof_handler_coarse.get_fe(fe_index_pair.first.first).degree;
          transfer.schemes[fe_index_pair.second].degree_fine =
            dof_handler_fine.get_fe(fe_index_pair.first.second).degree;
        }

      std::vector<unsigned int> n_dof_indices_fine(fe_index_pairs.size() + 1);
      std::vector<unsigned int> n_dof_indices_coarse(fe_index_pairs.size() + 1);
      std::vector<unsigned int> cell_no(fe_index_pairs.size() + 1, 0);

      {
        std::vector<std::vector<unsigned int>> lexicographic_numbering_fine(
          fe_index_pairs.size());
        std::vector<std::vector<unsigned int>> lexicographic_numbering_coarse(
          fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_coarse(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_coarse_lex(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_fine(fe_index_pairs.size());
        std::vector<std::vector<types::global_dof_index>>
          local_dof_indices_fine_lex(fe_index_pairs.size());

        for (const auto &fe_index_pair : fe_index_pairs)
          {
            local_dof_indices_coarse[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse);
            local_dof_indices_coarse_lex[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse);
            local_dof_indices_fine[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine);
            local_dof_indices_fine_lex[fe_index_pair.second].resize(
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine);

            n_dof_indices_fine[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_fine *
              transfer.schemes[fe_index_pair.second].n_coarse_cells;
            n_dof_indices_coarse[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_dofs_per_cell_coarse *
              transfer.schemes[fe_index_pair.second].n_coarse_cells;
            cell_no[fe_index_pair.second + 1] =
              transfer.schemes[fe_index_pair.second].n_coarse_cells;

            const auto reference_cell =
              dof_handler_fine.get_fe(fe_index_pair.first.second)
                .reference_cell();

            Assert(reference_cell ==
                     dof_handler_coarse.get_fe(fe_index_pair.first.first)
                       .reference_cell(),
                   ExcNotImplemented());

            // ------------------- lexicographic_numbering  --------------------
            if (has_tp_structure)
              {
                const Quadrature<1> dummy_quadrature(
                  std::vector<Point<1>>(1, Point<1>()));
                internal::MatrixFreeFunctions::ShapeInfo<
                  VectorizedArray<Number>>
                  shape_info;
                shape_info.reinit(dummy_quadrature,
                                  dof_handler_fine.get_fe(
                                    fe_index_pair.first.second),
                                  0);
                lexicographic_numbering_fine[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;

                shape_info.reinit(dummy_quadrature,
                                  dof_handler_coarse.get_fe(
                                    fe_index_pair.first.first),
                                  0);
                lexicographic_numbering_coarse[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;
              }
            else
              {
                const auto dummy_quadrature =
                  reference_cell.template get_gauss_type_quadrature<dim>(1);

                internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
                shape_info.reinit(dummy_quadrature,
                                  dof_handler_fine.get_fe(
                                    fe_index_pair.first.second),
                                  0);
                lexicographic_numbering_fine[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;

                shape_info.reinit(dummy_quadrature,
                                  dof_handler_coarse.get_fe(
                                    fe_index_pair.first.first),
                                  0);
                lexicographic_numbering_coarse[fe_index_pair.second] =
                  shape_info.lexicographic_numbering;
              }
          }

        for (unsigned int i = 0; i < fe_index_pairs.size(); ++i)
          {
            n_dof_indices_fine[i + 1] += n_dof_indices_fine[i];
            n_dof_indices_coarse[i + 1] += n_dof_indices_coarse[i];
            cell_no[i + 1] += cell_no[i];
          }

        // ------------------------------ indices  -----------------------------

        transfer.constraint_info_coarse.reinit(
          dof_handler_coarse,
          cell_no.back(),
          constraints_coarse.n_constraints() > 0 &&
            use_fast_hanging_node_algorithm(dof_handler_coarse,
                                            mg_level_coarse));
        transfer.constraint_info_coarse.set_locally_owned_indices(
          (mg_level_coarse == numbers::invalid_unsigned_int) ?
            dof_handler_coarse.locally_owned_dofs() :
            dof_handler_coarse.locally_owned_mg_dofs(mg_level_coarse));

        transfer.constraint_info_fine.reinit(cell_no.back());
        transfer.constraint_info_fine.set_locally_owned_indices(
          (mg_level_fine == numbers::invalid_unsigned_int) ?
            dof_handler_fine.locally_owned_dofs() :
            dof_handler_fine.locally_owned_mg_dofs(mg_level_fine));

        process_cells([&](const auto &cell_coarse, const auto &cell_fine) {
          const auto fe_pair_no =
            fe_index_pairs[std::pair<unsigned int, unsigned int>(
              cell_coarse->active_fe_index(), cell_fine->active_fe_index())];

          // parent
          {
            transfer.constraint_info_coarse.read_dof_indices(
              cell_no[fe_pair_no],
              mg_level_coarse,
              cell_coarse,
              constraints_coarse,
              {});
          }

          // child
          {
            cell_fine->get_dof_indices(local_dof_indices_fine[fe_pair_no]);

            for (unsigned int i = 0;
                 i < transfer.schemes[fe_pair_no].n_dofs_per_cell_fine;
                 ++i)
              local_dof_indices_fine_lex[fe_pair_no][i] = local_dof_indices_fine
                [fe_pair_no][lexicographic_numbering_fine[fe_pair_no][i]];

            transfer.constraint_info_fine.read_dof_indices(
              cell_no[fe_pair_no], local_dof_indices_fine_lex[fe_pair_no], {});
          }

          // move pointers
          {
            cell_no[fe_pair_no]++;
          }
        });
      }

      {
        transfer.partitioner_coarse = transfer.constraint_info_coarse.finalize(
          dof_handler_coarse.get_mpi_communicator());
        transfer.vec_coarse.reinit(transfer.partitioner_coarse);

        transfer.partitioner_fine = transfer.constraint_info_fine.finalize(
          dof_handler_fine.get_mpi_communicator());
        transfer.vec_fine.reinit(transfer.partitioner_fine);
      }

      // ------------------------- prolongation matrix -------------------------
      for (const auto &fe_index_pair_ : fe_index_pairs)
        compute_prolongation_and_restriction_matrices<dim, Number>(
          dof_handler_fine.get_fe(fe_index_pair_.first.second),
          dof_handler_coarse.get_fe(fe_index_pair_.first.first),
          transfer.schemes[fe_index_pair_.second]);

      // ------------------------------- weights -------------------------------
      if (transfer.fine_element_is_continuous)
        setup_weights(constraints_fine, transfer, is_feq);
    }
  };



  template <typename Number>
  struct SimpleVectorDataExchange
  {
    SimpleVectorDataExchange(
      const std::shared_ptr<const Utilities::MPI::Partitioner>
                            &embedded_partitioner,
      AlignedVector<Number> &buffer)
      : embedded_partitioner(embedded_partitioner)
      , buffer(buffer)
    {}

    template <typename VectorType>
    void
    update_ghost_values(const VectorType &vec) const
    {
      update_ghost_values_start(vec);
      update_ghost_values_finish(vec);
    }

    template <typename VectorType>
    void
    update_ghost_values_start(const VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      buffer.resize_fast(embedded_partitioner->n_import_indices());

      embedded_partitioner
        ->template export_to_ghosted_array_start<Number, MemorySpace::Host>(
          0,
          dealii::ArrayView<const Number>(
            vec.begin(), embedded_partitioner->locally_owned_size()),
          dealii::ArrayView<Number>(buffer.begin(), buffer.size()),
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    update_ghost_values_finish(const VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      embedded_partitioner
        ->template export_to_ghosted_array_finish<Number, MemorySpace::Host>(
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);

      vec.set_ghost_state(true);
#endif
    }

    template <typename VectorType>
    void
    compress(VectorType &vec) const
    {
      compress_start(vec);
      compress_finish(vec);
    }

    template <typename VectorType>
    void
    compress_start(VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      buffer.resize_fast(embedded_partitioner->n_import_indices());

      embedded_partitioner
        ->template import_from_ghosted_array_start<Number, MemorySpace::Host>(
          VectorOperation::add,
          0,
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          dealii::ArrayView<Number>(buffer.begin(), buffer.size()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    compress_finish(VectorType &vec) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)vec;
#else
      const auto &vector_partitioner = vec.get_partitioner();

      embedded_partitioner
        ->template import_from_ghosted_array_finish<Number, MemorySpace::Host>(
          VectorOperation::add,
          dealii::ArrayView<const Number>(buffer.begin(), buffer.size()),
          dealii::ArrayView<Number>(vec.begin(),
                                    embedded_partitioner->locally_owned_size()),
          dealii::ArrayView<Number>(
            const_cast<Number *>(vec.begin()) +
              embedded_partitioner->locally_owned_size(),
            vector_partitioner->n_ghost_indices()),
          requests);
#endif
    }

    template <typename VectorType>
    void
    zero_out_ghost_values(const VectorType &vec) const
    {
      const auto &vector_partitioner = vec.get_partitioner();

      ArrayView<Number> ghost_array(const_cast<VectorType &>(vec).begin() +
                                      vector_partitioner->locally_owned_size(),
                                    vector_partitioner->n_ghost_indices());

      for (const auto &my_ghosts :
           embedded_partitioner->ghost_indices_within_larger_ghost_set())
        for (unsigned int j = my_ghosts.first; j < my_ghosts.second; ++j)
          ghost_array[j] = 0.;

      vec.set_ghost_state(false);
    }

  private:
    const std::shared_ptr<const Utilities::MPI::Partitioner>
                                     embedded_partitioner;
    dealii::AlignedVector<Number>   &buffer;
    mutable std::vector<MPI_Request> requests;
  };

} // namespace internal



namespace MGTransferGlobalCoarseningTools
{
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &fine_triangulation_in)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

    coarse_grid_triangulations.back().reset(&fine_triangulation_in, [](auto *) {
      // empty deleter, since fine_triangulation_in is an external field
      // and its destructor is called somewhere else
    });

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    Assert(
      (dynamic_cast<
         const parallel::fullydistributed::Triangulation<dim, spacedim> *>(
         &fine_triangulation_in) == nullptr),
      ExcMessage(
        "Triangulations of type parallel::fullydistributed::Triangulation are "
        "not supported by this function!"));

    const auto create_new_empty_triangulation =
      [&]() -> std::shared_ptr<Triangulation<dim, spacedim>> {
#ifdef DEAL_II_WITH_P4EST
      if (const auto fine_triangulation = dynamic_cast<
            const parallel::distributed::Triangulation<dim, spacedim> *>(
            &fine_triangulation_in))
        return std::make_shared<
          parallel::distributed::Triangulation<dim, spacedim>>(
          fine_triangulation->get_mpi_communicator());
      else
#endif
#ifdef DEAL_II_WITH_MPI
        if (const auto fine_triangulation = dynamic_cast<
              const parallel::shared::Triangulation<dim, spacedim> *>(
              &fine_triangulation_in))
        return std::make_shared<parallel::shared::Triangulation<dim, spacedim>>(
          fine_triangulation->get_mpi_communicator(),
          Triangulation<dim, spacedim>::none,
          fine_triangulation->with_artificial_cells());
      else
#endif
        return std::make_shared<Triangulation<dim, spacedim>>();
    };

    const unsigned int max_level = fine_triangulation_in.n_global_levels() - 1;

    // clear 'eliminate_unrefined_islands' from MeshSmoothing flags
    // to prevent unintentional refinement during coarsen_global()
    const auto mesh_smoothing =
      static_cast<typename Triangulation<dim, spacedim>::MeshSmoothing>(
        fine_triangulation_in.get_mesh_smoothing() &
        ~(Triangulation<dim, spacedim>::eliminate_unrefined_islands));

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // copy triangulation
        auto new_tria = create_new_empty_triangulation();
        new_tria->copy_triangulation(*coarse_grid_triangulations[l]);
        new_tria->set_mesh_smoothing(mesh_smoothing);

        // coarsen mesh
        new_tria->coarsen_global();

        // save mesh
        coarse_grid_triangulations[l - 1] = new_tria;
      }

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim>                         &fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool keep_fine_triangulation,
    const bool repartition_fine_triangulation)
  {
    std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
      coarse_grid_triangulations(fine_triangulation_in.n_global_levels());

#ifndef DEAL_II_WITH_P4EST
    DEAL_II_NOT_IMPLEMENTED();
    (void)policy;
    (void)keep_fine_triangulation;
    (void)repartition_fine_triangulation;
#else
    const auto fine_triangulation =
      dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
        &fine_triangulation_in);

    Assert(fine_triangulation, ExcNotImplemented());

    const auto comm = fine_triangulation->get_mpi_communicator();

    if (keep_fine_triangulation == true &&
        repartition_fine_triangulation == false)
      {
        coarse_grid_triangulations.back().reset(&fine_triangulation_in,
                                                [](auto *) {
                                                  // empty deleter, since
                                                  // fine_triangulation_in is an
                                                  // external field and its
                                                  // destructor is called
                                                  // somewhere else
                                                });
      }
    else
      {
        // create triangulation description
        const auto construction_data =
          repartition_fine_triangulation ?
            TriangulationDescription::Utilities::
              create_description_from_triangulation(
                *fine_triangulation, policy.partition(*fine_triangulation)) :
            TriangulationDescription::Utilities::
              create_description_from_triangulation(*fine_triangulation, comm);

        // create new triangulation
        const auto new_fine_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            new_fine_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        new_fine_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations.back() = new_fine_triangulation;
      }

    // for a single level nothing has to be done
    if (fine_triangulation_in.n_global_levels() == 1)
      return coarse_grid_triangulations;

    parallel::distributed::Triangulation<dim, spacedim> temp_triangulation(
      comm);

    if (keep_fine_triangulation == true)
      temp_triangulation.copy_triangulation(*fine_triangulation);

    auto *temp_triangulation_ptr =
      keep_fine_triangulation ? &temp_triangulation : fine_triangulation;

    // clear 'eliminate_unrefined_islands' from MeshSmoothing flags
    // to prevent unintentional refinement during coarsen_global()
    const auto mesh_smoothing =
      static_cast<typename Triangulation<dim, spacedim>::MeshSmoothing>(
        temp_triangulation_ptr->get_mesh_smoothing() &
        ~(Triangulation<dim, spacedim>::eliminate_unrefined_islands));
    temp_triangulation_ptr->set_mesh_smoothing(mesh_smoothing);

    const unsigned int max_level = fine_triangulation->n_global_levels() - 1;

    // create coarse meshes
    for (unsigned int l = max_level; l > 0; --l)
      {
        // coarsen mesh
        temp_triangulation_ptr->coarsen_global();

        // create triangulation description
        const auto construction_data = TriangulationDescription::Utilities::
          create_description_from_triangulation(
            *temp_triangulation_ptr, policy.partition(*temp_triangulation_ptr));

        // create new triangulation
        const auto level_triangulation = std::make_shared<
          parallel::fullydistributed::Triangulation<dim, spacedim>>(comm);

        for (const auto i : fine_triangulation->get_manifold_ids())
          if (i != numbers::flat_manifold_id)
            level_triangulation->set_manifold(
              i, fine_triangulation->get_manifold(i));

        level_triangulation->create_triangulation(construction_data);

        // save mesh
        coarse_grid_triangulations[l - 1] = level_triangulation;
      }

    // recover MeshSmoothing flags in case we used the fine_triangulation
    // to build the sequence
    if (keep_fine_triangulation == false)
      fine_triangulation->set_mesh_smoothing(
        coarse_grid_triangulations.back()->get_mesh_smoothing());
#endif

    AssertDimension(coarse_grid_triangulations[0]->n_global_levels(), 1);

    return coarse_grid_triangulations;
  }



  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim>                   &fine_triangulation_in,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation)
  {
    // remove const and convert it to flag
    return create_geometric_coarsening_sequence(
      const_cast<Triangulation<dim, spacedim> &>(fine_triangulation_in),
      policy,
      true,
      repartition_fine_triangulation);
  }

} // namespace MGTransferGlobalCoarseningTools



template <typename VectorType>
MGTwoLevelTransferBase<VectorType>::MGTwoLevelTransferBase()
  : vec_fine_needs_ghost_update(true)
{}



template <typename VectorType>
void
MGTwoLevelTransferBase<VectorType>::prolongate_and_add(
  VectorType       &dst,
  const VectorType &src) const
{
  const bool  use_dst_inplace = this->vec_fine.size() == 0;
  auto *const vec_fine_ptr    = use_dst_inplace ? &dst : &this->vec_fine;
  Assert(vec_fine_ptr->get_partitioner().get() == this->partitioner_fine.get(),
         ExcInternalError());

  const bool        use_src_inplace = this->vec_coarse.size() == 0;
  const auto *const vec_coarse_ptr = use_src_inplace ? &src : &this->vec_coarse;
  Assert(vec_coarse_ptr->get_partitioner().get() ==
           this->partitioner_coarse.get(),
         ExcInternalError());

  const bool src_ghosts_have_been_set = src.has_ghost_elements();

  if (use_src_inplace == false)
    this->vec_coarse.copy_locally_owned_data_from(src);

  if ((use_src_inplace == false) || (src_ghosts_have_been_set == false))
    this->update_ghost_values(*vec_coarse_ptr);

  if (use_dst_inplace == false)
    *vec_fine_ptr = Number(0.);

  this->prolongate_and_add_internal(*vec_fine_ptr, *vec_coarse_ptr);

  if (this->vec_fine_needs_ghost_update || use_dst_inplace == false)
    this->compress(*vec_fine_ptr, VectorOperation::add);

  if (use_dst_inplace == false)
    dst += this->vec_fine;

  if (use_src_inplace && (src_ghosts_have_been_set == false))
    this->zero_out_ghost_values(*vec_coarse_ptr);
}



namespace internal
{
  // Helper class to compute correct weights, which works by simply using
  // the degrees of freedom stored in MatrixFree, bypassing the hanging node
  // interpolation matrices.
  template <int dim, typename Number>
  class FEEvaluationNoConstraints : public FEEvaluation<dim, -1, 0, 1, Number>
  {
  public:
    FEEvaluationNoConstraints(const MatrixFree<dim, Number> &data,
                              const unsigned int             dof_index,
                              const unsigned int             component = 0)
      : FEEvaluation<dim, -1, 0, 1, Number>(data, dof_index, 0, component)
    {}

    template <typename VectorType>
    void
    read_dof_values_unconstrained(VectorType &src)
    {
      std::array<VectorType *, 1> src_vector{{&src}};
      internal::VectorReader<Number, VectorizedArray<Number>> reader;
      this->template read_write_operation<VectorType>(
        reader,
        src_vector,
        {},
        std::bitset<VectorizedArray<Number>::size()>().flip(),
        false);
    }

    template <typename VectorType>
    void
    distribute_local_to_global_unconstrained(VectorType &dst)
    {
      std::array<VectorType *, 1> dst_vector{{&dst}};
      internal::VectorDistributorLocalToGlobal<Number, VectorizedArray<Number>>
        writer;
      this->template read_write_operation<VectorType>(
        writer,
        dst_vector,
        {},
        std::bitset<VectorizedArray<Number>::size()>().flip(),
        false);
    }
  };
} // namespace internal

template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::prolongate_and_add_internal(
  VectorType       &dst,
  const VectorType &src) const
{
  if (matrix_free_data.get() != nullptr)
    {
      matrix_free_data->matrix_free_fine->template cell_loop<VectorType, int>(
        [&](const MatrixFree<dim, Number> &data,
            VectorType                    &dst,
            const int &,
            const std::pair<unsigned int, unsigned int> &range) {
          for (unsigned int comp = 0; comp < n_components; ++comp)
            {
              internal::FEEvaluationNoConstraints<dim, Number> eval_fine(
                data, matrix_free_data->dof_handler_index_fine, comp);
              FEEvaluation<dim, -1, 0, 1, Number> eval_coarse(
                *matrix_free_data->matrix_free_coarse,
                matrix_free_data->dof_handler_index_coarse,
                0,
                comp);

              internal::CellTransferFactory cell_transfer(
                eval_fine.get_shape_info().data[0].fe_degree,
                eval_coarse.get_shape_info().data[0].fe_degree);
              for (unsigned int cell = range.first; cell < range.second; ++cell)
                {
                  eval_fine.reinit(cell);
                  eval_coarse.reinit(
                    matrix_free_data->cell_list_fine_to_coarse[cell]);

                  eval_coarse.read_dof_values(src);
                  if (schemes[0].prolongation_matrix.empty() == false)
                    {
                      internal::
                        CellProlongator<dim, double, VectorizedArrayType>
                          cell_prolongator(schemes[0].prolongation_matrix,
                                           eval_coarse.begin_dof_values(),
                                           eval_fine.begin_dof_values());

                      if (schemes[0].prolongation_matrix.size() <
                          eval_fine.dofs_per_cell * eval_coarse.dofs_per_cell)
                        cell_transfer.run(cell_prolongator);
                      else
                        cell_prolongator.run_full(eval_fine.dofs_per_cell,
                                                  eval_coarse.dofs_per_cell);
                    }
                  else
                    for (unsigned int i = 0; i < eval_coarse.dofs_per_cell; ++i)
                      eval_fine.begin_dof_values()[i] =
                        eval_coarse.begin_dof_values()[i];

                  if (weights.size() > 0)
                    {
                      const VectorizedArrayType *cell_weights =
                        weights.data() + weights_start[cell];
                      if (weights_are_compressed[cell])
                        internal::weight_fe_q_dofs_by_entity<
                          dim,
                          -1,
                          VectorizedArrayType>(
                          cell_weights,
                          1,
                          eval_fine.get_shape_info().data[0].fe_degree + 1,
                          eval_fine.begin_dof_values());
                      else
                        for (unsigned int i = 0; i < eval_fine.dofs_per_cell;
                             ++i)
                          eval_fine.begin_dof_values()[i] *= cell_weights[i];
                    }
                  eval_fine.distribute_local_to_global_unconstrained(dst);
                }
            }
        },
        dst,
        0);
    }
  else
    {
      const unsigned int n_lanes = VectorizedArrayType::size();

      AlignedVector<VectorizedArrayType> evaluation_data_fine;
      AlignedVector<VectorizedArrayType> evaluation_data_coarse;

      unsigned int cell_counter  = 0;
      unsigned int batch_counter = 0;

      for (const auto &scheme : schemes)
        {
          if (scheme.n_coarse_cells == 0)
            continue;

          const bool needs_interpolation =
            scheme.prolongation_matrix.empty() == false;

          evaluation_data_fine.clear();
          evaluation_data_coarse.clear();

          const unsigned int max_n_dofs_per_cell =
            std::max(scheme.n_dofs_per_cell_fine,
                     scheme.n_dofs_per_cell_coarse);
          evaluation_data_fine.resize(max_n_dofs_per_cell);
          evaluation_data_coarse.resize(max_n_dofs_per_cell);

          internal::CellTransferFactory cell_transfer(scheme.degree_fine,
                                                      scheme.degree_coarse);

          const unsigned int n_scalar_dofs_fine =
            scheme.n_dofs_per_cell_fine / n_components;
          const unsigned int n_scalar_dofs_coarse =
            scheme.n_dofs_per_cell_coarse / n_components;

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
               cell += n_lanes, ++batch_counter)
            {
              const unsigned int n_lanes_filled =
                (cell + n_lanes > scheme.n_coarse_cells) ?
                  (scheme.n_coarse_cells - cell) :
                  n_lanes;

              // read from src vector (similar to
              // FEEvaluation::read_dof_values())
              internal::VectorReader<Number, VectorizedArrayType> reader;
              constraint_info_coarse.read_write_operation(
                reader,
                src,
                evaluation_data_coarse.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_coarse,
                true);
              constraint_info_coarse.apply_hanging_node_constraints(
                cell_counter, n_lanes_filled, false, evaluation_data_coarse);

              // ---------------------------- coarse ---------------------------
              if (needs_interpolation)
                for (int c = n_components - 1; c >= 0; --c)
                  {
                    internal::CellProlongator<dim, double, VectorizedArrayType>
                      cell_prolongator(scheme.prolongation_matrix,
                                       evaluation_data_coarse.begin() +
                                         c * n_scalar_dofs_coarse,
                                       evaluation_data_fine.begin() +
                                         c * n_scalar_dofs_fine);

                    if (scheme.prolongation_matrix.size() <
                        n_scalar_dofs_fine * n_scalar_dofs_coarse)
                      cell_transfer.run(cell_prolongator);
                    else
                      cell_prolongator.run_full(n_scalar_dofs_fine,
                                                n_scalar_dofs_coarse);
                  }
              else
                evaluation_data_fine = evaluation_data_coarse; // TODO
              // ------------------------------ fine ---------------------------

              // weight
              if (weights.size() > 0)
                {
                  const VectorizedArrayType *cell_weights =
                    weights.data() + weights_start[batch_counter];
                  if (weights_are_compressed[batch_counter])
                    internal::
                      weight_fe_q_dofs_by_entity<dim, -1, VectorizedArrayType>(
                        cell_weights,
                        n_components,
                        scheme.degree_fine + 1,
                        evaluation_data_fine.begin());
                  else
                    for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine;
                         ++i)
                      evaluation_data_fine[i] *= cell_weights[i];
                }

              // add into dst vector
              internal::VectorDistributorLocalToGlobal<Number,
                                                       VectorizedArrayType>
                writer;
              constraint_info_fine.read_write_operation(
                writer,
                dst,
                evaluation_data_fine.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_fine,
                false);

              cell_counter += n_lanes_filled;
            }
        }
    }
}



template <typename VectorType>
void
MGTwoLevelTransferBase<VectorType>::restrict_and_add(
  VectorType       &dst,
  const VectorType &src) const
{
  const bool        use_src_inplace = this->vec_fine.size() == 0;
  const auto *const vec_fine_ptr    = use_src_inplace ? &src : &this->vec_fine;
  Assert(vec_fine_ptr->get_partitioner().get() == this->partitioner_fine.get(),
         ExcInternalError());

  const bool  use_dst_inplace = this->vec_coarse.size() == 0;
  auto *const vec_coarse_ptr  = use_dst_inplace ? &dst : &this->vec_coarse;
  Assert(vec_coarse_ptr->get_partitioner().get() ==
           this->partitioner_coarse.get(),
         ExcInternalError());

  const bool src_ghosts_have_been_set = src.has_ghost_elements();

  if (use_src_inplace == false)
    this->vec_fine.copy_locally_owned_data_from(src);

  if ((use_src_inplace == false) ||
      (vec_fine_needs_ghost_update && (src_ghosts_have_been_set == false)))
    this->update_ghost_values(*vec_fine_ptr);

  if (use_dst_inplace == false)
    *vec_coarse_ptr = Number(0.0);

  // since we might add into the ghost values and call compress
  this->zero_out_ghost_values(*vec_coarse_ptr);

  this->restrict_and_add_internal(*vec_coarse_ptr, *vec_fine_ptr);

  // clean up related to update_ghost_values()
  if (vec_fine_needs_ghost_update == false && use_src_inplace == false)
    this->zero_out_ghost_values(*vec_fine_ptr); // internal vector (DG)
  else if (vec_fine_needs_ghost_update && use_src_inplace == false)
    vec_fine_ptr->set_ghost_state(false); // internal vector (CG)
  else if (vec_fine_needs_ghost_update && (src_ghosts_have_been_set == false))
    this->zero_out_ghost_values(*vec_fine_ptr); // external vector

  this->compress(*vec_coarse_ptr, VectorOperation::add);

  if (use_dst_inplace == false)
    dst += this->vec_coarse;
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::restrict_and_add_internal(
  VectorType       &dst,
  const VectorType &src) const
{
  if (matrix_free_data.get() != nullptr)
    {
      int dummy = 0;
      matrix_free_data->matrix_free_fine->template cell_loop<int, VectorType>(
        [&](const MatrixFree<dim, Number> &data,
            int &,
            const VectorType                            &src,
            const std::pair<unsigned int, unsigned int> &range) {
          for (unsigned int comp = 0; comp < n_components; ++comp)
            {
              internal::FEEvaluationNoConstraints<dim, Number> eval_fine(
                data, matrix_free_data->dof_handler_index_fine, comp);
              FEEvaluation<dim, -1, 0, 1, Number> eval_coarse(
                *matrix_free_data->matrix_free_coarse,
                matrix_free_data->dof_handler_index_coarse,
                0,
                comp);

              internal::CellTransferFactory cell_transfer(
                eval_fine.get_shape_info().data[0].fe_degree,
                eval_coarse.get_shape_info().data[0].fe_degree);
              for (unsigned int cell = range.first; cell < range.second; ++cell)
                {
                  eval_fine.reinit(cell);
                  eval_coarse.reinit(
                    matrix_free_data->cell_list_fine_to_coarse[cell]);

                  eval_fine.read_dof_values_unconstrained(src);
                  if (weights.size() > 0)
                    {
                      const VectorizedArrayType *cell_weights =
                        weights.data() + weights_start[cell];
                      if (weights_are_compressed[cell])
                        internal::weight_fe_q_dofs_by_entity<
                          dim,
                          -1,
                          VectorizedArrayType>(
                          cell_weights,
                          1,
                          eval_fine.get_shape_info().data[0].fe_degree + 1,
                          eval_fine.begin_dof_values());
                      else
                        for (unsigned int i = 0; i < eval_fine.dofs_per_cell;
                             ++i)
                          eval_fine.begin_dof_values()[i] *= cell_weights[i];
                    }

                  if (schemes[0].prolongation_matrix.empty() == false)
                    {
                      internal::CellRestrictor<dim, double, VectorizedArrayType>
                        cell_restrictor(schemes[0].prolongation_matrix,
                                        eval_fine.begin_dof_values(),
                                        eval_coarse.begin_dof_values());

                      if (schemes[0].prolongation_matrix.size() <
                          eval_fine.dofs_per_cell * eval_coarse.dofs_per_cell)
                        cell_transfer.run(cell_restrictor);
                      else
                        cell_restrictor.run_full(eval_fine.dofs_per_cell,
                                                 eval_coarse.dofs_per_cell);
                    }
                  else
                    for (unsigned int i = 0; i < eval_fine.dofs_per_cell; ++i)
                      eval_coarse.begin_dof_values()[i] =
                        eval_fine.begin_dof_values()[i];

                  eval_coarse.distribute_local_to_global(dst);
                }
            }
        },
        dummy,
        src);
    }
  else
    {
      const unsigned int n_lanes = VectorizedArrayType::size();

      AlignedVector<VectorizedArrayType> evaluation_data_fine;
      AlignedVector<VectorizedArrayType> evaluation_data_coarse;

      unsigned int cell_counter  = 0;
      unsigned int batch_counter = 0;

      for (const auto &scheme : schemes)
        {
          if (scheme.n_coarse_cells == 0)
            continue;

          const bool needs_interpolation =
            scheme.prolongation_matrix.empty() == false;

          evaluation_data_fine.clear();
          evaluation_data_coarse.clear();

          const unsigned int max_n_dofs_per_cell =
            std::max(scheme.n_dofs_per_cell_fine,
                     scheme.n_dofs_per_cell_coarse);
          evaluation_data_fine.resize(max_n_dofs_per_cell);
          evaluation_data_coarse.resize(max_n_dofs_per_cell);

          internal::CellTransferFactory cell_transfer(scheme.degree_fine,
                                                      scheme.degree_coarse);

          const unsigned int n_scalar_dofs_fine =
            scheme.n_dofs_per_cell_fine / n_components;
          const unsigned int n_scalar_dofs_coarse =
            scheme.n_dofs_per_cell_coarse / n_components;

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
               cell += n_lanes, ++batch_counter)
            {
              const unsigned int n_lanes_filled =
                (cell + n_lanes > scheme.n_coarse_cells) ?
                  (scheme.n_coarse_cells - cell) :
                  n_lanes;

              // read from source vector
              internal::VectorReader<Number, VectorizedArrayType> reader;
              constraint_info_fine.read_write_operation(
                reader,
                src,
                evaluation_data_fine.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_fine,
                false);

              // weight
              if (weights.size() > 0)
                {
                  const VectorizedArrayType *cell_weights =
                    weights.data() + weights_start[batch_counter];
                  if (weights_are_compressed[batch_counter])
                    internal::
                      weight_fe_q_dofs_by_entity<dim, -1, VectorizedArrayType>(
                        cell_weights,
                        n_components,
                        scheme.degree_fine + 1,
                        evaluation_data_fine.data());
                  else
                    for (unsigned int i = 0; i < scheme.n_dofs_per_cell_fine;
                         ++i)
                      evaluation_data_fine[i] *= cell_weights[i];
                }

              // ------------------------------ fine ---------------------------
              if (needs_interpolation)
                for (int c = n_components - 1; c >= 0; --c)
                  {
                    internal::CellRestrictor<dim, double, VectorizedArrayType>
                      cell_restrictor(scheme.prolongation_matrix,
                                      evaluation_data_fine.begin() +
                                        c * n_scalar_dofs_fine,
                                      evaluation_data_coarse.begin() +
                                        c * n_scalar_dofs_coarse);

                    if (scheme.prolongation_matrix.size() <
                        n_scalar_dofs_fine * n_scalar_dofs_coarse)
                      cell_transfer.run(cell_restrictor);
                    else
                      cell_restrictor.run_full(n_scalar_dofs_fine,
                                               n_scalar_dofs_coarse);
                  }
              else
                evaluation_data_coarse = evaluation_data_fine; // TODO
              // ----------------------------- coarse --------------------------

              // write into dst vector (similar to
              // FEEvaluation::distribute_global_to_local())
              internal::VectorDistributorLocalToGlobal<Number,
                                                       VectorizedArrayType>
                writer;
              constraint_info_coarse.apply_hanging_node_constraints(
                cell_counter, n_lanes_filled, true, evaluation_data_coarse);
              constraint_info_coarse.read_write_operation(
                writer,
                dst,
                evaluation_data_coarse.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_coarse,
                true);

              cell_counter += n_lanes_filled;
            }
        }
    }
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::interpolate(VectorType       &dst,
                                                 const VectorType &src) const
{
  if (matrix_free_data.get() != nullptr)
    {
      src.update_ghost_values();
      for (unsigned int comp = 0; comp < n_components; ++comp)
        {
          internal::FEEvaluationNoConstraints<dim, Number> eval_fine(
            *matrix_free_data->matrix_free_fine,
            matrix_free_data->dof_handler_index_fine,
            comp);
          FEEvaluation<dim, -1, 0, 1, Number> eval_coarse(
            *matrix_free_data->matrix_free_coarse,
            matrix_free_data->dof_handler_index_coarse,
            0,
            comp);

          internal::CellTransferFactory cell_transfer(
            eval_fine.get_shape_info().data[0].fe_degree,
            eval_coarse.get_shape_info().data[0].fe_degree);
          for (unsigned int cell = 0;
               cell < matrix_free_data->matrix_free_fine->n_cell_batches();
               ++cell)
            {
              eval_fine.reinit(cell);
              eval_coarse.reinit(
                matrix_free_data->cell_list_fine_to_coarse[cell]);

              eval_fine.read_dof_values_unconstrained(src);

              if (schemes[0].restriction_matrix.empty() == false)
                {
                  internal::CellRestrictor<dim, double, VectorizedArrayType>
                    cell_restrictor(schemes[0].restriction_matrix,
                                    eval_fine.begin_dof_values(),
                                    eval_coarse.begin_dof_values());

                  if (schemes[0].prolongation_matrix.size() <
                      eval_fine.dofs_per_cell * eval_coarse.dofs_per_cell)
                    cell_transfer.run(cell_restrictor);
                  else
                    cell_restrictor.run_full(eval_fine.dofs_per_cell,
                                             eval_coarse.dofs_per_cell);
                }
              else
                for (unsigned int i = 0; i < eval_fine.dofs_per_cell; ++i)
                  eval_coarse.begin_dof_values()[i] =
                    eval_fine.begin_dof_values()[i];

              eval_coarse.set_dof_values(dst);
            }
        }
      src.zero_out_ghost_values();
      dst.zero_out_ghost_values();
    }
  else
    {
      const unsigned int n_lanes = VectorizedArrayType::size();

      const bool        use_src_inplace = this->vec_fine.size() == 0;
      const auto *const vec_fine_ptr = use_src_inplace ? &src : &this->vec_fine;
      Assert(vec_fine_ptr->get_partitioner().get() ==
               this->partitioner_fine.get(),
             ExcInternalError());

      const bool  use_dst_inplace = this->vec_coarse.size() == 0;
      auto *const vec_coarse_ptr  = use_dst_inplace ? &dst : &this->vec_coarse;
      Assert(vec_coarse_ptr->get_partitioner().get() ==
               this->partitioner_coarse.get(),
             ExcInternalError());

      const bool src_ghosts_have_been_set = src.has_ghost_elements();

      if (use_src_inplace == false)
        this->vec_fine.copy_locally_owned_data_from(src);

      if ((use_src_inplace == false) || (this->vec_fine_needs_ghost_update &&
                                         (src_ghosts_have_been_set == false)))
        this->update_ghost_values(*vec_fine_ptr);

      *vec_coarse_ptr = Number(0.0);

      AlignedVector<VectorizedArrayType> evaluation_data_fine;
      AlignedVector<VectorizedArrayType> evaluation_data_coarse;

      unsigned int cell_counter = 0;

      for (const auto &scheme : schemes)
        {
          if (scheme.n_coarse_cells == 0)
            continue;

          if (scheme.n_dofs_per_cell_fine == 0 ||
              scheme.n_dofs_per_cell_coarse == 0)
            {
              cell_counter += scheme.n_coarse_cells;

              continue;
            }

          const bool needs_interpolation =
            scheme.restriction_matrix.empty() == false;

          // general case -> local restriction is needed
          evaluation_data_fine.resize(scheme.n_dofs_per_cell_fine);
          evaluation_data_coarse.resize(scheme.n_dofs_per_cell_fine);

          internal::CellTransferFactory cell_transfer(scheme.degree_fine,
                                                      scheme.degree_coarse);

          const unsigned int n_scalar_dofs_fine =
            scheme.n_dofs_per_cell_fine / n_components;
          const unsigned int n_scalar_dofs_coarse =
            scheme.n_dofs_per_cell_coarse / n_components;

          for (unsigned int cell = 0; cell < scheme.n_coarse_cells;
               cell += n_lanes)
            {
              const unsigned int n_lanes_filled =
                (cell + n_lanes > scheme.n_coarse_cells) ?
                  (scheme.n_coarse_cells - cell) :
                  n_lanes;

              // read from source vector
              internal::VectorReader<Number, VectorizedArrayType> reader;
              constraint_info_fine.read_write_operation(
                reader,
                *vec_fine_ptr,
                evaluation_data_fine.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_fine,
                false);

              // ------------------------------ fine ---------------------------
              if (needs_interpolation)
                for (int c = n_components - 1; c >= 0; --c)
                  {
                    internal::CellRestrictor<dim, double, VectorizedArrayType>
                      cell_restrictor(scheme.restriction_matrix,
                                      evaluation_data_fine.begin() +
                                        c * n_scalar_dofs_fine,
                                      evaluation_data_coarse.begin() +
                                        c * n_scalar_dofs_coarse);

                    if (scheme.restriction_matrix.size() <
                        n_scalar_dofs_fine * n_scalar_dofs_coarse)
                      cell_transfer.run(cell_restrictor);
                    else
                      cell_restrictor.run_full(n_scalar_dofs_fine,
                                               n_scalar_dofs_coarse);
                  }
              else
                evaluation_data_coarse = evaluation_data_fine; // TODO
              // ----------------------------- coarse --------------------------

              // write into dst vector (similar to
              // FEEvaluation::set_dof_values_plain())
              internal::VectorSetter<Number, VectorizedArrayType> writer;
              constraint_info_coarse.read_write_operation(
                writer,
                *vec_coarse_ptr,
                evaluation_data_coarse.data(),
                cell_counter,
                n_lanes_filled,
                scheme.n_dofs_per_cell_coarse,
                false);

              cell_counter += n_lanes_filled;
            }
        }

      // clean up related to update_ghost_values()
      if (use_src_inplace == false)
        vec_fine_ptr->set_ghost_state(false); // internal vector
      else if (this->fine_element_is_continuous &&
               (src_ghosts_have_been_set == false))
        this->zero_out_ghost_values(*vec_fine_ptr); // external vector

      if (use_dst_inplace == false)
        dst.copy_locally_owned_data_from(this->vec_coarse);
    }
}


namespace internal
{
  inline bool
  is_partitioner_contained(
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner)
  {
    // no external partitioner has been given
    if (external_partitioner.get() == nullptr)
      return false;

    // check if locally owned ranges are the same
    if (external_partitioner->size() != partitioner->size())
      return false;

    if (external_partitioner->locally_owned_range() !=
        partitioner->locally_owned_range())
      return false;

    const int ghosts_locally_contained =
      ((external_partitioner->ghost_indices() & partitioner->ghost_indices()) ==
       partitioner->ghost_indices()) ?
        1 :
        0;

    // check if ghost values are contained in external partititioner
    return Utilities::MPI::min(ghosts_locally_contained,
                               partitioner->get_mpi_communicator()) == 1;
  }

  inline std::shared_ptr<Utilities::MPI::Partitioner>
  create_embedded_partitioner(
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &larger_partitioner)
  {
    auto embedded_partitioner = std::make_shared<Utilities::MPI::Partitioner>(
      larger_partitioner->locally_owned_range(),
      larger_partitioner->get_mpi_communicator());

    embedded_partitioner->set_ghost_indices(
      partitioner->ghost_indices(), larger_partitioner->ghost_indices());

    return embedded_partitioner;
  }
} // namespace internal



template <typename VectorType>
template <int dim, std::size_t width, typename IndexType>
std::pair<bool, bool>
MGTwoLevelTransferBase<VectorType>::
  internal_enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
         &external_partitioner_fine,
    bool &vec_fine_needs_ghost_update,
    internal::MatrixFreeFunctions::
      ConstraintInfo<dim, VectorizedArray<Number, width>, IndexType>
                              &constraint_info_coarse,
    std::vector<unsigned int> &dof_indices_fine)
{
  std::pair<bool, bool> success_flags = {false, false};

  if (this->partitioner_coarse->is_globally_compatible(
        *external_partitioner_coarse))
    {
      this->vec_coarse.reinit(0);
      this->partitioner_coarse = external_partitioner_coarse;
      success_flags.first      = true;
    }
  else if (internal::is_partitioner_contained(this->partitioner_coarse,
                                              external_partitioner_coarse))
    {
      this->vec_coarse.reinit(0);

      for (auto &i : constraint_info_coarse.dof_indices)
        i = external_partitioner_coarse->global_to_local(
          this->partitioner_coarse->local_to_global(i));

      for (auto &i : constraint_info_coarse.plain_dof_indices)
        i = external_partitioner_coarse->global_to_local(
          this->partitioner_coarse->local_to_global(i));

      this->partitioner_coarse_embedded =
        internal::create_embedded_partitioner(this->partitioner_coarse,
                                              external_partitioner_coarse);

      this->partitioner_coarse = external_partitioner_coarse;
      success_flags.first      = true;
    }

  vec_fine_needs_ghost_update =
    Utilities::MPI::max(this->partitioner_fine->ghost_indices().n_elements(),
                        this->partitioner_fine->get_mpi_communicator()) != 0;

  if (this->partitioner_fine->is_globally_compatible(
        *external_partitioner_fine))
    {
      this->vec_fine.reinit(0);
      this->partitioner_fine = external_partitioner_fine;
      success_flags.second   = true;
    }
  else if (internal::is_partitioner_contained(this->partitioner_fine,
                                              external_partitioner_fine))
    {
      this->vec_fine.reinit(0);

      for (auto &i : dof_indices_fine)
        i = external_partitioner_fine->global_to_local(
          this->partitioner_fine->local_to_global(i));

      this->partitioner_fine_embedded =
        internal::create_embedded_partitioner(this->partitioner_fine,
                                              external_partitioner_fine);

      this->partitioner_fine = external_partitioner_fine;
      success_flags.second   = true;
    }

  return success_flags;
}



template <int dim, typename VectorType>
std::pair<bool, bool>
MGTwoLevelTransfer<dim, VectorType>::enable_inplace_operations_if_possible(
  const std::shared_ptr<const Utilities::MPI::Partitioner>
    &external_partitioner_coarse,
  const std::shared_ptr<const Utilities::MPI::Partitioner>
    &external_partitioner_fine)
{
  if (matrix_free_data.get() != nullptr)
    return std::make_pair(true, true);
  else
    return this->internal_enable_inplace_operations_if_possible(
      external_partitioner_coarse,
      external_partitioner_fine,
      this->vec_fine_needs_ghost_update,
      constraint_info_coarse,
      constraint_info_fine.dof_indices);
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::reinit_geometric_transfer(
  const DoFHandler<dim>           &dof_handler_fine,
  const DoFHandler<dim>           &dof_handler_coarse,
  const AffineConstraints<Number> &constraints_fine,
  const AffineConstraints<Number> &constraints_coarse,
  const unsigned int               mg_level_fine,
  const unsigned int               mg_level_coarse)
{
  matrix_free_data.reset();
  internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraints_fine,
    constraints_coarse,
    mg_level_fine,
    mg_level_coarse,
    *this);
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::reinit_polynomial_transfer(
  const DoFHandler<dim>           &dof_handler_fine,
  const DoFHandler<dim>           &dof_handler_coarse,
  const AffineConstraints<Number> &constraints_fine,
  const AffineConstraints<Number> &constraints_coarse,
  const unsigned int               mg_level_fine,
  const unsigned int               mg_level_coarse)
{
  matrix_free_data.reset();
  internal::MGTwoLevelTransferImplementation::reinit_polynomial_transfer(
    dof_handler_fine,
    dof_handler_coarse,
    constraints_fine,
    constraints_coarse,
    mg_level_fine,
    mg_level_coarse,
    *this);
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::reinit(
  const DoFHandler<dim>           &dof_handler_fine,
  const DoFHandler<dim>           &dof_handler_coarse,
  const AffineConstraints<Number> &constraints_fine,
  const AffineConstraints<Number> &constraints_coarse,
  const unsigned int               mg_level_fine,
  const unsigned int               mg_level_coarse)
{
  // determine if polynomial transfer can be performed via the following two
  // criteria:
  // 1) multigrid levels can be only used with polynomial transfer
  bool do_polynomial_transfer =
    (mg_level_fine != numbers::invalid_unsigned_int) ||
    (mg_level_coarse != numbers::invalid_unsigned_int);

  // 2) the meshes are identical
  if (do_polynomial_transfer == false)
    {
      const internal::CellIDTranslator<dim> cell_id_translator(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_fine.get_triangulation().n_global_levels());

      AssertDimension(
        dof_handler_fine.get_triangulation().n_global_coarse_cells(),
        dof_handler_coarse.get_triangulation().n_global_coarse_cells());
      AssertIndexRange(dof_handler_coarse.get_triangulation().n_global_levels(),
                       dof_handler_fine.get_triangulation().n_global_levels() +
                         1);

      IndexSet is_locally_owned_fine(cell_id_translator.size());
      IndexSet is_locally_owned_coarse(cell_id_translator.size());

      for (const auto &cell : dof_handler_fine.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        is_locally_owned_fine.add_index(cell_id_translator.translate(cell));

      for (const auto &cell : dof_handler_coarse.active_cell_iterators() |
                                IteratorFilters::LocallyOwnedCell())
        is_locally_owned_coarse.add_index(cell_id_translator.translate(cell));

      const MPI_Comm communicator = dof_handler_fine.get_mpi_communicator();

      const std::vector<unsigned int> owning_ranks =
        Utilities::MPI::compute_index_owner(is_locally_owned_fine,
                                            is_locally_owned_coarse,
                                            communicator);

      bool all_cells_found = true;

      for (unsigned i = 0; i < is_locally_owned_coarse.n_elements(); ++i)
        all_cells_found &= (owning_ranks[i] != numbers::invalid_unsigned_int);

      do_polynomial_transfer =
        Utilities::MPI::min(static_cast<unsigned int>(all_cells_found),
                            communicator) == 1;
    }

  if (do_polynomial_transfer)
    internal::MGTwoLevelTransferImplementation::reinit_polynomial_transfer(
      dof_handler_fine,
      dof_handler_coarse,
      constraints_fine,
      constraints_coarse,
      mg_level_fine,
      mg_level_coarse,
      *this);
  else
    internal::MGTwoLevelTransferImplementation::reinit_geometric_transfer(
      dof_handler_fine,
      dof_handler_coarse,
      constraints_fine,
      constraints_coarse,
      mg_level_fine,
      mg_level_coarse,
      *this);
}



template <int dim, typename VectorType>
void
MGTwoLevelTransfer<dim, VectorType>::reinit(
  const MatrixFree<dim, Number> &matrix_free_fine,
  const unsigned int             dof_no_fine,
  const MatrixFree<dim, Number> &matrix_free_coarse,
  const unsigned int             dof_no_coarse)
{
  matrix_free_data = std::make_unique<MatrixFreeRelatedData>();

  MatrixFreeRelatedData &data = *matrix_free_data;
  data.matrix_free_fine       = &matrix_free_fine;
  data.matrix_free_coarse     = &matrix_free_coarse;
  AssertIndexRange(dof_no_fine, matrix_free_fine.n_components());
  data.dof_handler_index_fine = dof_no_fine;
  AssertIndexRange(dof_no_coarse, matrix_free_coarse.n_components());
  data.dof_handler_index_coarse = dof_no_coarse;

  const DoFHandler<dim> &dof_fine =
    matrix_free_fine.get_dof_handler(dof_no_fine);
  const DoFHandler<dim> &dof_coarse =
    matrix_free_coarse.get_dof_handler(dof_no_coarse);

  Assert(&dof_fine.get_triangulation() == &dof_coarse.get_triangulation(),
         ExcMessage("You can only use this class if both MatrixFree objects "
                    "use the same underlying triangulation!"));
  AssertDimension(matrix_free_fine.get_mg_level(),
                  matrix_free_coarse.get_mg_level());

  AssertDimension(matrix_free_fine.n_physical_cells(),
                  matrix_free_coarse.n_physical_cells());

  Assert(dof_fine.get_fe_collection().size() == 1,
         ExcMessage("hp finite element objects are currently not supported "
                    "by this class."));
  Assert(dof_coarse.get_fe_collection().size() == 1,
         ExcMessage("hp finite element objects are currently not supported "
                    "by this class."));

  this->dof_handler_fine = &dof_fine;
  this->mg_level_fine    = matrix_free_fine.get_mg_level();

  const FiniteElement<dim> &fe_fine = dof_fine.get_fe();

  AssertDimension(fe_fine.n_components(), dof_coarse.get_fe().n_components());
  this->n_components = fe_fine.n_components();

  // Compute the link between cells on matrix_free_coarse and
  // matrix_free_fine, with the latter controlling the loop
  const Triangulation<dim> &tria = dof_fine.get_triangulation();
  std::vector<unsigned int> coarse_cell_indices(tria.n_active_cells());
  for (unsigned int cell = 0; cell < matrix_free_coarse.n_cell_batches();
       ++cell)
    for (unsigned int v = 0;
         v < matrix_free_coarse.n_active_entries_per_cell_batch(cell);
         ++v)
      {
        coarse_cell_indices[matrix_free_coarse.get_cell_iterator(cell, v)
                              ->active_cell_index()] =
          cell * VectorizedArrayType::size() + v;
      }

  std::array<unsigned int, VectorizedArrayType::size()> default_entry;
  default_entry.fill(numbers::invalid_unsigned_int);
  data.cell_list_fine_to_coarse.resize(matrix_free_fine.n_cell_batches(),
                                       default_entry);
  for (unsigned int cell = 0; cell < matrix_free_fine.n_cell_batches(); ++cell)
    for (unsigned int v = 0;
         v < matrix_free_fine.n_active_entries_per_cell_batch(cell);
         ++v)
      {
        data.cell_list_fine_to_coarse[cell][v] =
          coarse_cell_indices[matrix_free_fine.get_cell_iterator(cell, v)
                                ->active_cell_index()];
      }

  // Compute weights
  if (fe_fine.dofs_per_vertex > 0)
    {
      VectorType weight_vector;
      matrix_free_fine.initialize_dof_vector(weight_vector, dof_no_fine);

      internal::FEEvaluationNoConstraints<dim, Number> evaluator(
        matrix_free_fine, dof_no_fine);
      for (unsigned int cell = 0; cell < matrix_free_fine.n_cell_batches();
           ++cell)
        {
          evaluator.reinit(cell);
          for (unsigned int i = 0; i < evaluator.dofs_per_cell; ++i)
            evaluator.begin_dof_values()[i] = 1.;
          evaluator.distribute_local_to_global_unconstrained(weight_vector);
        }
      weight_vector.compress(VectorOperation::add);

      for (unsigned int i = 0; i < weight_vector.locally_owned_size(); ++i)
        if (weight_vector.local_element(i) > 0.)
          weight_vector.local_element(i) =
            Number(1.) / weight_vector.local_element(i);
      for (const unsigned int index :
           matrix_free_fine.get_constrained_dofs(dof_no_fine))
        weight_vector.local_element(index) = 0;
      weight_vector.update_ghost_values();

      weights.clear();
      weights.reserve(matrix_free_fine.n_cell_batches() *
                      Utilities::pow(3, dim));
      weights_start.clear();
      weights_start.resize(matrix_free_fine.n_cell_batches());
      weights_are_compressed.clear();
      weights_are_compressed.resize(weights_start.size());
      std::array<VectorizedArrayType, Utilities::pow(3, dim)>
        weights_compressed;
      for (unsigned int cell = 0; cell < matrix_free_fine.n_cell_batches();
           ++cell)
        {
          weights_start[cell] = weights.size();

          weights_compressed.fill(VectorizedArray<Number>());
          evaluator.reinit(cell);
          evaluator.read_dof_values_unconstrained(weight_vector);

          const bool can_compress_weights =
            internal::compute_weights_fe_q_dofs_by_entity<dim, -1>(
              evaluator.begin_dof_values(),
              1,
              fe_fine.degree + 1,
              weights_compressed.data());
          if (can_compress_weights && fe_fine.degree > 1)
            {
              for (const VectorizedArrayType weight : weights_compressed)
                weights.push_back(weight);
              weights_are_compressed[cell] = 1;
            }
          else
            for (unsigned int i = 0; i < evaluator.dofs_per_cell; ++i)
              weights.push_back(evaluator.begin_dof_values()[i]);
        }
    }

  // Compute interpolation matrices, possibly in the 1D version
  schemes.resize(1);
  internal::MGTwoLevelTransferImplementation::
    compute_prolongation_and_restriction_matrices<dim, Number>(
      fe_fine, dof_coarse.get_fe(), schemes[0]);

  // We internally call the ghost updates through the matrix-free framework
  this->vec_fine_needs_ghost_update = false;

  this->partitioner_fine = matrix_free_fine.get_vector_partitioner(dof_no_fine);
  this->partitioner_coarse =
    matrix_free_coarse.get_vector_partitioner(dof_no_coarse);
}



template <int dim, typename VectorType>
bool
MGTwoLevelTransfer<dim, VectorType>::fast_polynomial_transfer_supported(
  const unsigned int fe_degree_fine,
  const unsigned int fe_degree_coarse)
{
  internal::CellTransferFactory cell_transfer(fe_degree_fine, fe_degree_coarse);
  internal::CellProlongatorTest cell_transfer_test;

  return cell_transfer.run(cell_transfer_test);
}



template <int dim, typename VectorType>
std::size_t
MGTwoLevelTransfer<dim, VectorType>::memory_consumption() const
{
  std::size_t size = 0;

  for (const auto &scheme : schemes)
    {
      size += scheme.prolongation_matrix.memory_consumption();
      size += scheme.restriction_matrix.memory_consumption();
    }

  if (matrix_free_data.get() != nullptr)
    size += MemoryConsumption::memory_consumption(
      matrix_free_data->cell_list_fine_to_coarse);

  size += this->partitioner_fine->memory_consumption();
  size += this->partitioner_coarse->memory_consumption();
  size += this->vec_fine.memory_consumption();
  size += this->vec_coarse.memory_consumption();
  size += weights.memory_consumption();
  size += MemoryConsumption::memory_consumption(weights_start);
  size += MemoryConsumption::memory_consumption(weights_are_compressed);
  size += constraint_info_coarse.memory_consumption();
  size += constraint_info_fine.memory_consumption();

  return size;
}



template <typename VectorType>
void
MGTwoLevelTransferBase<VectorType>::update_ghost_values(
  const VectorType &vec) const
{
  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .update_ghost_values(vec);
  else if ((vec.get_partitioner().get() == this->partitioner_fine.get()) &&
           (this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .update_ghost_values(vec);
  else
    vec.update_ghost_values();
}



template <typename VectorType>
void
MGTwoLevelTransferBase<VectorType>::compress(
  VectorType                   &vec,
  const VectorOperation::values op) const
{
  Assert(op == VectorOperation::add, ExcNotImplemented());

  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .compress(vec);
  else if ((vec.get_partitioner().get() == this->partitioner_fine.get()) &&
           (this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .compress(vec);
  else
    vec.compress(op);
}



template <typename VectorType>
void
MGTwoLevelTransferBase<VectorType>::zero_out_ghost_values(
  const VectorType &vec) const
{
  if ((vec.get_partitioner().get() == this->partitioner_coarse.get()) &&
      (this->partitioner_coarse_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(
      this->partitioner_coarse_embedded, this->buffer_coarse_embedded)
      .zero_out_ghost_values(vec);
  else if ((vec.get_partitioner().get() == (this->partitioner_fine.get()) &&
            this->partitioner_fine_embedded != nullptr))
    internal::SimpleVectorDataExchange<Number>(this->partitioner_fine_embedded,
                                               this->buffer_fine_embedded)
      .zero_out_ghost_values(vec);
  else
    vec.zero_out_ghost_values();
}



template <int dim, typename Number, typename MemorySpace>
MGTransferMF<dim, Number, MemorySpace>::MGTransferMF()
{
  this->transfer.clear();
  this->internal_transfer.clear();
}



template <int dim, typename Number, typename MemorySpace>
MGTransferMF<dim, Number, MemorySpace>::MGTransferMF(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->transfer.clear();
  this->internal_transfer.clear();
  this->initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::initialize_constraints(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->mg_constrained_dofs = &mg_constrained_dofs;
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::initialize_internal_transfer(
  const DoFHandler<dim>                          &dof_handler,
  const ObserverPointer<const MGConstrainedDoFs> &mg_constrained_dofs)
{
  const unsigned int min_level = 0;
  const unsigned int max_level =
    dof_handler.get_triangulation().n_global_levels() - 1;

  MGLevelObject<AffineConstraints<typename VectorType::value_type>> constraints(
    min_level, max_level);

  if (mg_constrained_dofs)
    for (unsigned int l = min_level; l <= max_level; ++l)
      mg_constrained_dofs->merge_constraints(
        constraints[l],
        l,
        /*add_boundary_indices*/ true,
        /*add_refinement_edge_indices*/ false,
        /*add_level_constraints*/ true,
        /*add_user_constraints*/ true);

  this->internal_transfer.resize(min_level, max_level);

  for (unsigned int l = min_level; l < max_level; ++l)
    internal_transfer[l + 1].reinit_geometric_transfer(
      dof_handler, dof_handler, constraints[l + 1], constraints[l], l + 1, l);
}



template <int dim, typename Number, typename MemorySpace>
std::pair<const DoFHandler<dim> *, unsigned int>
MGTransferMF<dim, Number, MemorySpace>::get_dof_handler_fine() const
{
  if (this->transfer.n_levels() <= 1)
    // single level: the information cannot be retrieved
    return {nullptr, numbers::invalid_unsigned_int};

  if (const auto t = dynamic_cast<const MGTwoLevelTransfer<
        dim,
        LinearAlgebra::distributed::Vector<Number, MemorySpace>> *>(
        this->transfer[this->transfer.max_level()].get()))
    {
      return {t->dof_handler_fine, t->mg_level_fine};
    }

  if constexpr (std::is_same_v<MemorySpace, ::dealii::MemorySpace::Host>)
    {
      // MGTwoLevelTransferNonNested transfer is only instantiated for Host
      // memory:
      if (const auto t = dynamic_cast<const MGTwoLevelTransferNonNested<
            dim,
            LinearAlgebra::distributed::Vector<Number, MemorySpace>> *>(
            this->transfer[this->transfer.max_level()].get()))
        {
          return {t->dof_handler_fine, t->mg_level_fine};
        }
    }

  {
    DEAL_II_NOT_IMPLEMENTED();
    return {nullptr, numbers::invalid_unsigned_int};
  }
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::
  fill_and_communicate_copy_indices_global_coarsening(
    const DoFHandler<dim> &dof_handler_out)
{
  const auto dof_handler_and_level_in = get_dof_handler_fine();
  const auto dof_handler_in           = dof_handler_and_level_in.first;
  const auto level_in                 = dof_handler_and_level_in.second;

  if ((dof_handler_in == nullptr) || (dof_handler_in == &dof_handler_out))
    return; // nothing to do

  this->copy_indices.resize(1);
  this->copy_indices[0].reinit(2, dof_handler_out.n_locally_owned_dofs());

  std::vector<types::global_dof_index> dof_indices_in;
  std::vector<types::global_dof_index> dof_indices_out;

  this->perform_plain_copy = true;

  const auto &is_out = (level_in == numbers::invalid_unsigned_int) ?
                         dof_handler_out.locally_owned_dofs() :
                         dof_handler_out.locally_owned_mg_dofs(level_in);

  const auto &is_in = (level_in == numbers::invalid_unsigned_int) ?
                        dof_handler_in->locally_owned_dofs() :
                        dof_handler_in->locally_owned_mg_dofs(level_in);

  internal::loop_over_active_or_level_cells(
    dof_handler_in->get_triangulation(), level_in, [&](const auto &cell) {
      const auto cell_id = cell->id();

      Assert(
        dof_handler_out.get_triangulation().contains_cell(cell_id),
        ExcMessage(
          "DoFHandler instances used for set up of MGTransferMF and copy_to_mg(), "
          "copy_from_mg(), or interpolate_to_mg() are not compatible."));

      if (level_in == numbers::invalid_unsigned_int)
        {
          const auto cell_in  = cell->as_dof_handler_iterator(*dof_handler_in);
          const auto cell_out = dof_handler_out.get_triangulation()
                                  .create_cell_iterator(cell_id)
                                  ->as_dof_handler_iterator(dof_handler_out);

          AssertDimension(cell_in->get_fe().n_dofs_per_cell(),
                          cell_out->get_fe().n_dofs_per_cell());

          dof_indices_in.resize(cell_in->get_fe().n_dofs_per_cell());
          dof_indices_out.resize(cell_out->get_fe().n_dofs_per_cell());

          cell_in->get_dof_indices(dof_indices_in);
          cell_out->get_dof_indices(dof_indices_out);
        }
      else
        {
          const auto cell_in =
            cell->as_dof_handler_level_iterator(*dof_handler_in);
          const auto cell_out =
            dof_handler_out.get_triangulation()
              .create_cell_iterator(cell_id)
              ->as_dof_handler_level_iterator(dof_handler_out);

          AssertDimension(cell_in->get_fe().n_dofs_per_cell(),
                          cell_out->get_fe().n_dofs_per_cell());

          dof_indices_in.resize(cell_in->get_fe().n_dofs_per_cell());
          dof_indices_out.resize(cell_out->get_fe().n_dofs_per_cell());

          cell_in->get_mg_dof_indices(dof_indices_in);
          cell_out->get_mg_dof_indices(dof_indices_out);
        }

      this->perform_plain_copy &= (dof_indices_in == dof_indices_out);

      for (unsigned int i = 0; i < dof_indices_in.size(); ++i)
        if (is_out.is_element(dof_indices_out[i]))
          this->copy_indices[0](1,
                                is_out.index_within_set(dof_indices_out[i])) =
            is_in.index_within_set(dof_indices_in[i]);
    });


  this->perform_plain_copy =
    Utilities::MPI::max(this->perform_plain_copy ? 1 : 0,
                        dof_handler_out.get_mpi_communicator()) != 0;

  if (this->perform_plain_copy)
    {
      this->copy_indices.clear();
    }
  else
    {
      this->perform_renumbered_plain_copy = true;
      this->solution_copy_indices         = this->copy_indices;
    }
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::build(
  const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    &external_partitioners)
{
  this->external_partitioners = external_partitioners;

  if (this->external_partitioners.size() > 0)
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();

      AssertDimension(this->external_partitioners.size(), transfer.n_levels());

      for (unsigned int l = min_level + 1; l <= max_level; ++l)
        transfer[l]->enable_inplace_operations_if_possible(
          this->external_partitioners[l - 1 - min_level],
          this->external_partitioners[l - min_level]);
    }
  else
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();

      for (unsigned int l = min_level + 1; l <= max_level; ++l)
        {
          if (l == min_level + 1)
            this->external_partitioners.push_back(
              transfer[l]->partitioner_coarse);

          this->external_partitioners.push_back(transfer[l]->partitioner_fine);
        }
    }

  this->perform_plain_copy            = true;
  this->perform_renumbered_plain_copy = false;
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::build(
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  if (initialize_dof_vector)
    {
      const unsigned int min_level = transfer.min_level();
      const unsigned int max_level = transfer.max_level();
      const unsigned int n_levels  = transfer.n_levels();

      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        external_partitioners(n_levels);

      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          LinearAlgebra::distributed::Vector<typename VectorType::value_type,
                                             MemorySpace>
            vector;
          initialize_dof_vector(l, vector);
          external_partitioners[l - min_level] = vector.get_partitioner();
        }

      this->build(external_partitioners);
    }
  else
    {
      this->build();
    }
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::build(
  const DoFHandler<dim> &dof_handler,
  const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
    &external_partitioners)
{
  const bool use_local_smoothing =
    this->transfer.n_levels() == 0 || this->internal_transfer.n_levels() > 0;

  if (use_local_smoothing)
    {
      this->initialize_internal_transfer(dof_handler,
                                         this->mg_constrained_dofs);
      this->initialize_transfer_references(internal_transfer);
    }

  this->build(external_partitioners);

  if (use_local_smoothing)
    this->fill_and_communicate_copy_indices(dof_handler);
  else
    this->fill_and_communicate_copy_indices_global_coarsening(dof_handler);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::build(
  const DoFHandler<dim> &dof_handler,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
{
  const bool use_local_smoothing =
    this->transfer.n_levels() == 0 || this->internal_transfer.n_levels() > 0;

  if (use_local_smoothing)
    {
      this->initialize_internal_transfer(dof_handler,
                                         this->mg_constrained_dofs);
      this->initialize_transfer_references(internal_transfer);
    }

  this->build(initialize_dof_vector);

  if (use_local_smoothing)
    this->fill_and_communicate_copy_indices(dof_handler);
  else
    this->fill_and_communicate_copy_indices_global_coarsening(dof_handler);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::prolongate(const unsigned int to_level,
                                                   VectorType        &dst,
                                                   const VectorType  &src) const
{
  dst = Number(0.0);
  prolongate_and_add(to_level, dst, src);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::prolongate_and_add(
  const unsigned int to_level,
  VectorType        &dst,
  const VectorType  &src) const
{
  this->transfer[to_level]->prolongate_and_add(dst, src);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::restrict_and_add(
  const unsigned int from_level,
  VectorType        &dst,
  const VectorType  &src) const
{
  this->transfer[from_level]->restrict_and_add(dst, src);
}



template <int dim, typename Number, typename MemorySpace>
void
MGTransferMF<dim, Number, MemorySpace>::assert_dof_handler(
  const DoFHandler<dim> &dof_handler_out) const
{
  (void)dof_handler_out;
  if constexpr (running_in_debug_mode())
    {
      const auto dof_handler_and_level_in = get_dof_handler_fine();
      const auto dof_handler_in           = dof_handler_and_level_in.first;
      const auto level_in                 = dof_handler_and_level_in.second;

      if ((dof_handler_out.n_dofs() == 0) ||  // dummy DoFHandler
          (dof_handler_in == nullptr) ||      // single level
          (dof_handler_in == &dof_handler_out // same DoFHandler
           ))
        return; // nothing to do

      if (this->perform_plain_copy)
        {
          // global-coarsening path: compare indices of cells

          std::vector<types::global_dof_index> dof_indices_in;
          std::vector<types::global_dof_index> dof_indices_out;

          internal::loop_over_active_or_level_cells(
            dof_handler_in->get_triangulation(),
            level_in,
            [&](const auto &cell) {
              const auto cell_id = cell->id();

              Assert(
                dof_handler_out.get_triangulation().contains_cell(cell_id),
                ExcMessage(
                  "DoFHandler instances used for set up of MGTransferMF and copy_to_mg(), "
                  "copy_from_mg(), or interpolate_to_mg() are not compatible."));

              if (level_in == numbers::invalid_unsigned_int)
                {
                  const auto cell_in =
                    cell->as_dof_handler_iterator(*dof_handler_in);
                  const auto cell_out =
                    dof_handler_out.get_triangulation()
                      .create_cell_iterator(cell_id)
                      ->as_dof_handler_iterator(dof_handler_out);

                  AssertDimension(cell_in->get_fe().n_dofs_per_cell(),
                                  cell_out->get_fe().n_dofs_per_cell());

                  dof_indices_in.resize(cell_in->get_fe().n_dofs_per_cell());
                  dof_indices_out.resize(cell_out->get_fe().n_dofs_per_cell());

                  cell_in->get_dof_indices(dof_indices_in);
                  cell_out->get_dof_indices(dof_indices_out);
                }
              else
                {
                  const auto cell_in =
                    cell->as_dof_handler_level_iterator(*dof_handler_in);
                  const auto cell_out =
                    dof_handler_out.get_triangulation()
                      .create_cell_iterator(cell_id)
                      ->as_dof_handler_level_iterator(dof_handler_out);

                  AssertDimension(cell_in->get_fe().n_dofs_per_cell(),
                                  cell_out->get_fe().n_dofs_per_cell());

                  dof_indices_in.resize(cell_in->get_fe().n_dofs_per_cell());
                  dof_indices_out.resize(cell_out->get_fe().n_dofs_per_cell());

                  cell_in->get_mg_dof_indices(dof_indices_in);
                  cell_out->get_mg_dof_indices(dof_indices_out);
                }

              Assert(
                dof_indices_in == dof_indices_out,
                ExcMessage(
                  "DoFHandler instances used for set up of MGTransferMF and copy_to_mg(), "
                  "copy_from_mg(), or interpolate_to_mg() are not compatible."));
            });
        }
      else if (this->perform_renumbered_plain_copy)
        {
          // nothing to do
        }
    }
}



template <int dim, typename Number, typename MemorySpace>
std::size_t
MGTransferMF<dim, Number, MemorySpace>::memory_consumption() const
{
  std::size_t size = 0;

  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  for (unsigned int l = min_level + 1; l <= max_level; ++l)
    size += this->transfer[l]->memory_consumption();

  return size;
}



template <int dim, typename Number, typename MemorySpace>
inline unsigned int
MGTransferMF<dim, Number, MemorySpace>::min_level() const
{
  return transfer.min_level();
}



template <int dim, typename Number, typename MemorySpace>
inline unsigned int
MGTransferMF<dim, Number, MemorySpace>::max_level() const
{
  return transfer.max_level();
}


template <int dim, typename Number, typename MemorySpace>
inline void
MGTransferMF<dim, Number, MemorySpace>::clear()
{
  MGLevelGlobalTransfer<VectorType>::clear();

  internal_transfer.clear();
  transfer.clear();
  external_partitioners.clear();
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>
    &transfer_operator)
  : MGTransferBlockMatrixFreeBase<
      dim,
      Number,
      MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>(true)
{
  this->transfer_operators = {&transfer_operator};
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const MGConstrainedDoFs &mg_constrained_dofs)
  : MGTransferBlockMatrixFreeBase<
      dim,
      Number,
      MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>(true)
{
  initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number>
MGTransferBlockMF<dim, Number>::MGTransferBlockMF(
  const std::vector<MGConstrainedDoFs> &mg_constrained_dofs)
  : MGTransferBlockMatrixFreeBase<
      dim,
      Number,
      MGTransferMF<dim, Number, ::dealii::MemorySpace::Host>>(false)
{
  initialize_constraints(mg_constrained_dofs);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::initialize_constraints(
  const MGConstrainedDoFs &mg_constrained_dofs)
{
  this->transfer_operators_internal.clear();
  this->transfer_operators.clear();

  Assert(this->same_for_all,
         ExcMessage("This object was initialized with support for usage with "
                    "one DoFHandler for each block, but this method assumes "
                    "that the same DoFHandler is used for all the blocks!"));

  this->transfer_operators_internal.emplace_back(mg_constrained_dofs);
  this->transfer_operators = {&transfer_operators_internal.back()};
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::initialize_constraints(
  const std::vector<MGConstrainedDoFs> &mg_constrained_dofs)
{
  this->transfer_operators_internal.clear();
  this->transfer_operators.clear();

  Assert(!this->same_for_all,
         ExcMessage("This object was initialized with support for using "
                    "the same DoFHandler for all the blocks, but this "
                    "method assumes that there is a separate DoFHandler "
                    "for each block!"));

  for (const auto &dofs : mg_constrained_dofs)
    this->transfer_operators_internal.emplace_back(dofs);

  for (const auto &transfer : this->transfer_operators_internal)
    this->transfer_operators.emplace_back(&transfer);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::build(const DoFHandler<dim> &dof_handler)
{
  AssertDimension(transfer_operators.size(), 1);
  this->transfer_operators_internal[0].build(dof_handler);
}



template <int dim, typename Number>
void
MGTransferBlockMF<dim, Number>::build(
  const std::vector<const DoFHandler<dim> *> &dof_handler)
{
  AssertDimension(transfer_operators.size(), dof_handler.size());
  AssertDimension(transfer_operators_internal.size(), dof_handler.size());

  for (unsigned int i = 0; i < dof_handler.size(); ++i)
    this->transfer_operators_internal[i].build(*dof_handler[i]);
}



template <int dim, typename Number>
const MGTransferMF<dim, Number, ::dealii::MemorySpace::Host> &
MGTransferBlockMF<dim, Number>::get_matrix_free_transfer(
  const unsigned int b) const
{
  AssertIndexRange(b, transfer_operators.size());
  return *transfer_operators[b];
}

namespace internal
{
  template <int dim, typename Number>
  std::shared_ptr<NonMatching::MappingInfo<dim, dim, Number>>
  fill_mapping_info(const Utilities::MPI::RemotePointEvaluation<dim> &rpe)
  {
    const auto &cell_data = rpe.get_cell_data();

    std::vector<typename Triangulation<dim>::active_cell_iterator>
                                         cell_iterators;
    std::vector<std::vector<Point<dim>>> unit_points_vector;

    for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
      {
        typename Triangulation<dim>::active_cell_iterator cell(
          &rpe.get_triangulation(),
          cell_data.cells[i].first,
          cell_data.cells[i].second);

        const ArrayView<const Point<dim>> unit_points(
          cell_data.reference_point_values.data() +
            cell_data.reference_point_ptrs[i],
          cell_data.reference_point_ptrs[i + 1] -
            cell_data.reference_point_ptrs[i]);

        cell_iterators.emplace_back(cell);
        unit_points_vector.emplace_back(unit_points.begin(), unit_points.end());
      }

    typename NonMatching::MappingInfo<dim, dim, Number>::AdditionalData ad;
    ad.store_cells = true;

    auto mapping_info =
      std::make_shared<NonMatching::MappingInfo<dim, dim, Number>>(
        rpe.get_mapping(), update_values, ad);
    mapping_info->reinit_cells(cell_iterators, unit_points_vector);

    return mapping_info;
  }

  /**
   * This function provides information which DoF index is associated with
   * a support point.
   *
   * @param[in] dof_handler DoFHandler with @c FE_DGQ or @c FE_Q elements
   * providing DoF indices which are collected at support points.
   * @param[in] dof_handler_support_points DoFHandler with one component used
   * to determine support point indices (the underlying finite element is @c
   * FE_Q or @c FE_DGQ in case of polynomial degree 0).
   * @param[in] constraint AffineConstrains associated with @p dof_handler.
   *   Only unconstrained DoFs are considered
   * @return a tuple containing 0) local support point indices,
   *   1) pointers to global DoF indices, and 2) global DoF indices.
   */
  template <int dim, int spacedim, typename Number>
  std::tuple<std::vector<unsigned int>,
             std::vector<unsigned int>,
             std::vector<types::global_dof_index>>
  support_point_indices_to_dof_indices(
    const DoFHandler<dim, spacedim>         &dof_handler,
    const DoFHandler<dim, spacedim>         &dof_handler_support_points,
    const dealii::AffineConstraints<Number> &constraint)
  {
    // in case a FE_DGQ space of order 0 is provided, DoFs indices are always
    // uniquely assigned to support points (they are always defined in the
    // center of the element) and are never shared at vertices or faces.
    Assert((dynamic_cast<const FE_DGQ<dim, spacedim> *>(
              &dof_handler.get_fe().base_element(0)) != nullptr) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(
                &dof_handler.get_fe().base_element(0)) != nullptr) ||
             (dynamic_cast<const FE_SimplexP<dim, spacedim> *>(
                &dof_handler.get_fe().base_element(0)) != nullptr) ||
             (dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(
                &dof_handler.get_fe().base_element(0)) != nullptr),
           ExcMessage("Function expects FE_DGQ, FE_Q, FE_SimplexP, or "
                      "FE_SimplexDGP in dof_handler."));

    Assert(
      (dynamic_cast<const FE_Q<dim, spacedim> *>(
         &dof_handler_support_points.get_fe().base_element(0)) != nullptr ||
       dynamic_cast<const FE_SimplexP<dim, spacedim> *>(
         &dof_handler_support_points.get_fe().base_element(0)) != nullptr) ||
        ((dynamic_cast<const FE_DGQ<dim, spacedim> *>(
            &dof_handler_support_points.get_fe().base_element(0)) != nullptr ||
          dynamic_cast<const FE_SimplexDGP<dim, spacedim> *>(
            &dof_handler_support_points.get_fe().base_element(0)) != nullptr) &&
         dof_handler_support_points.get_fe().degree == 0),
      ExcMessage("Function expects (FE_DGQ||FE_SimplexDGP)&&degree==0 or "
                 "(FE_Q||FE_SimplexP) in dof_handler_support_points."));

    Assert(
      dof_handler_support_points.get_fe().n_components() == 1,
      ExcMessage(
        "dof_handler_support_points needs element with exactly one component."));

    Assert(&dof_handler.get_triangulation() ==
             &dof_handler_support_points.get_triangulation(),
           ExcMessage("DoFHandlers need the same underlying triangulation."));

    Assert(dof_handler.get_fe().degree ==
             dof_handler_support_points.get_fe().degree,
           ExcMessage("DoFHandlers need the same degree."));

    Assert(dof_handler.get_fe().is_primitive(),
           ExcMessage("Only primitive elements are allowed."));

    const auto degree        = dof_handler.get_fe().degree;
    const auto dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
    const auto support_points_per_cell =
      dof_handler_support_points.get_fe().n_dofs_per_cell();

    std::vector<std::pair<unsigned int, types::global_dof_index>>
      support_point_dofs;
    support_point_dofs.reserve(dof_handler.n_locally_owned_dofs());

    const unsigned int n_components = dof_handler.get_fe().n_components();

    // fill support_point_dofs
    {
      // Support points have a hierarchic numbering, L2 DoFs have
      // lexicographic numbering. Therefore, we need to convert the DoF
      // indices if DoFHandler is L2 conforming and has degree > 0.
      const bool needs_conversion =
        dof_handler.get_fe().conforming_space ==
          FiniteElementData<dim>::Conformity::L2 &&
        (dof_handler.get_fe().degree > 0) &&
        dof_handler.get_fe().reference_cell().is_hyper_cube();
      std::vector<unsigned int> lexicographic_to_hierarchic;
      if (needs_conversion)
        lexicographic_to_hierarchic =
          FETools::lexicographic_to_hierarchic_numbering<dim>(degree);

      const Utilities::MPI::Partitioner partitioner_support_points(
        dof_handler_support_points.locally_owned_dofs(),
        dof_handler_support_points.get_mpi_communicator());

      const Utilities::MPI::Partitioner partitioner_dof(
        dof_handler.locally_owned_dofs(),
        DoFTools::extract_locally_relevant_dofs(dof_handler),
        dof_handler.get_mpi_communicator());

      std::vector<bool> dof_processed(partitioner_dof.locally_owned_size() +
                                        partitioner_dof.n_ghost_indices(),
                                      false);


      std::vector<types::global_dof_index> support_point_indices(
        support_points_per_cell);
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      std::vector<std::pair<unsigned int, types::global_dof_index>>
        support_point_dofs_comp;
      support_point_dofs_comp.reserve(n_components);

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned() || cell->is_ghost())
            {
              const auto cell_support_point =
                cell->as_dof_handler_iterator(dof_handler_support_points);

              cell_support_point->get_dof_indices(support_point_indices);
              cell->get_dof_indices(dof_indices);

              // collect unconstrained DoFs for support point. In case of DG
              // elements with polynomial degree > 0 or continuous elements
              // with multiple components, more DoFs are associated to the
              // same support point.
              for (unsigned int i = 0; i < support_point_indices.size(); ++i)
                if (partitioner_support_points.in_local_range(
                      support_point_indices[i]))
                  {
                    for (unsigned int c = 0; c < n_components; ++c)
                      {
                        const auto global_dof_idx =
                          needs_conversion ?
                            dof_indices
                              [dof_handler.get_fe().component_to_system_index(
                                c, lexicographic_to_hierarchic[i])] :
                            dof_indices[dof_handler.get_fe()
                                          .component_to_system_index(c, i)];

                        const auto local_dof_idx =
                          partitioner_dof.global_to_local(global_dof_idx);

                        AssertIndexRange(local_dof_idx, dof_processed.size());

                        if (dof_processed[local_dof_idx] == false)
                          {
                            if (!constraint.is_constrained(global_dof_idx))
                              support_point_dofs_comp.emplace_back(
                                partitioner_support_points.global_to_local(
                                  support_point_indices[i]),
                                global_dof_idx);

                            dof_processed[local_dof_idx] = true;
                          }
                      }

                    Assert(support_point_dofs_comp.size() == 0 ||
                             support_point_dofs_comp.size() == n_components,
                           ExcNotImplemented());

                    if (support_point_dofs_comp.empty() == false)
                      support_point_dofs.insert(support_point_dofs.end(),
                                                support_point_dofs_comp.begin(),
                                                support_point_dofs_comp.end());

                    support_point_dofs_comp.clear();
                  }
            }
        }
    }

    // sort for support points (stable sort needed for multiple components)
    std::stable_sort(support_point_dofs.begin(),
                     support_point_dofs.end(),
                     [](const auto &a, const auto &b) {
                       return a.first < b.first;
                     });

    // convert to CRS format
    std::vector<types::global_dof_index> dof_indices;
    dof_indices.reserve(support_point_dofs.size());
    std::vector<unsigned int> dof_ptrs;
    dof_ptrs.reserve(dof_handler_support_points.n_locally_owned_dofs() + 1);
    dof_ptrs.push_back(0);
    std::vector<unsigned int> support_point_indices;
    support_point_indices.reserve(
      dof_handler_support_points.n_locally_owned_dofs());

    auto it = support_point_dofs.begin();
    while (it != support_point_dofs.end())
      {
        const unsigned int index = std::get<0>(*it);
        while (it != support_point_dofs.end() && it->first == index)
          {
            dof_indices.push_back(it->second);
            ++it;
          }
        support_point_indices.push_back(index);
        dof_ptrs.push_back(dof_indices.size() / n_components);
      }

    return std::make_tuple(std::move(support_point_indices),
                           std::move(dof_ptrs),
                           std::move(dof_indices));
  }


  /**
   * Create DoFHandler with unique support points.
   */
  template <int dim, int spacedim>
  std::shared_ptr<const DoFHandler<dim, spacedim>>
  create_support_point_dof_handler(const DoFHandler<dim, spacedim> &dof_handler)
  {
    const auto &fe           = dof_handler.get_fe();
    const auto &tria         = dof_handler.get_triangulation();
    const auto  degree       = fe.degree;
    const auto  n_components = fe.n_components();

    if (n_components == 1 &&
        ((fe.reference_cell().is_hyper_cube() ||
          fe.reference_cell().is_simplex()) &&
         (fe.conforming_space == FiniteElementData<dim>::Conformity::H1 ||
          degree == 0)))
      {
        // in case a DG space of order 0 is provided, DoFs indices are always
        // uniquely assigned to support points (they are always defined in the
        // center of the element) and are never shared at vertices or faces.
        return std::shared_ptr<const DoFHandler<dim, spacedim>>(&dof_handler,
                                                                [](auto *) {});
      }
    else
      {
        // Create dummy dof handler for support point numbering.
        // Unique support points are generally numbered according to FE_Q with
        // n_components==1. If degree==0 we use FE_DGQ which ensures a unique
        // support point numbering since the support point is located in the
        // center of the cell.
        auto dof_handler_support_points =
          std::make_shared<DoFHandler<dim, spacedim>>(tria);

        if (fe.reference_cell().is_simplex() && (degree == 0))
          dof_handler_support_points->distribute_dofs(
            FE_SimplexDGP<dim, spacedim>(degree));
        else if (fe.reference_cell().is_simplex())
          dof_handler_support_points->distribute_dofs(
            FE_SimplexP<dim, spacedim>(degree));
        else if (degree == 0)
          dof_handler_support_points->distribute_dofs(
            FE_DGQ<dim, spacedim>(degree));
        else
          dof_handler_support_points->distribute_dofs(
            FE_Q<dim, spacedim>(degree));

        return dof_handler_support_points;
      }
  }

  // Loop over cells and collect unique set of points
  template <int dim, typename Number>
  std::tuple<std::vector<Point<dim>>,
             std::vector<unsigned int>,
             std::vector<types::global_dof_index>>
  collect_unconstrained_unique_support_points(
    const DoFHandler<dim>                   &dof_handler,
    const Mapping<dim>                      &mapping,
    const dealii::AffineConstraints<Number> &constraint)
  {
    AssertThrow(dof_handler.get_fe().has_support_points(), ExcNotImplemented());

    // create DoFHandler for support points
    const auto dof_handler_support_points =
      create_support_point_dof_handler(dof_handler);

    // compute mapping: index of locally owned support points to (global) DoF
    // indices
    const auto support_point_dofs_crs =
      support_point_indices_to_dof_indices(dof_handler,
                                           *dof_handler_support_points,
                                           constraint);

    const std::vector<unsigned int> &local_support_point_indices =
      std::get<0>(support_point_dofs_crs);

    // compute locally owned support points
    std::vector<Point<dim>> points;
    points.resize(local_support_point_indices.size());

    const auto locally_owned_support_point =
      dof_handler_support_points->locally_owned_dofs();
    std::vector<unsigned int> indices_state(
      locally_owned_support_point.n_elements(), numbers::invalid_unsigned_int);

    AssertIndexRange(local_support_point_indices.size(),
                     indices_state.size() + 1);

    for (unsigned int i = 0; i < local_support_point_indices.size(); ++i)
      indices_state[local_support_point_indices[i]] = i;

    const auto   &fe_support_point = dof_handler_support_points->get_fe();
    FEValues<dim> fe_values(mapping,
                            fe_support_point,
                            Quadrature<dim>(
                              fe_support_point.get_unit_support_points()),
                            update_quadrature_points);

    std::vector<types::global_dof_index> dof_indices(
      fe_support_point.n_dofs_per_cell());

    for (const auto &cell :
         dof_handler_support_points->active_cell_iterators() |
           IteratorFilters::LocallyOwnedCell())
      {
        fe_values.reinit(cell);
        cell->get_dof_indices(dof_indices);

        for (const unsigned int q : fe_values.quadrature_point_indices())
          if (locally_owned_support_point.is_element(dof_indices[q]))
            {
              const auto index =
                locally_owned_support_point.index_within_set(dof_indices[q]);

              if (indices_state[index] != numbers::invalid_unsigned_int)
                {
                  points[indices_state[index]] = fe_values.quadrature_point(q);
                  indices_state[index]         = numbers::invalid_unsigned_int;
                }
            }
      }

    return std::make_tuple(
      std::move(points),
      std::move(std::get<1>(support_point_dofs_crs)),  // global_dofs_ptrs
      std::move(std::get<2>(support_point_dofs_crs))); // global_dofs_indices
  }

} // namespace internal



template <int dim, typename VectorType>
MGTwoLevelTransferNonNested<dim, VectorType>::MGTwoLevelTransferNonNested(
  const AdditionalData &data)
  : additional_data(data)
  , rpe(typename Utilities::MPI::RemotePointEvaluation<dim>::AdditionalData(
      data.tolerance,
      false,
      data.rtree_level,
      {}))
{}

template <int dim, typename VectorType>
void
MGTwoLevelTransferNonNested<dim, VectorType>::reinit(
  const DoFHandler<dim>           &dof_handler_fine,
  const DoFHandler<dim>           &dof_handler_coarse,
  const Mapping<dim>              &mapping_fine,
  const Mapping<dim>              &mapping_coarse,
  const AffineConstraints<Number> &constraint_fine,
  const AffineConstraints<Number> &constraint_coarse)
{
  AssertThrow(dof_handler_coarse.get_fe().has_support_points(),
              ExcNotImplemented());
  Assert(dof_handler_coarse.get_fe().n_components() > 0 &&
           dof_handler_fine.get_fe().n_components() > 0,
         ExcNotImplemented());

  this->dof_handler_fine = &dof_handler_fine;
  this->mg_level_fine    = numbers::invalid_unsigned_int;

  this->fine_element_is_continuous =
    dof_handler_fine.get_fe().n_dofs_per_vertex() > 0;

  // collect points, ptrs, and global indices
  const auto points_ptrs_indices =
    internal::collect_unconstrained_unique_support_points(dof_handler_fine,
                                                          mapping_fine,
                                                          constraint_fine);
  const auto &global_dof_indices = std::get<2>(points_ptrs_indices);

  // create partitioners and internal vectors
  {
    this->partitioner_coarse =
      internal::MGTwoLevelTransferImplementation::create_coarse_partitioner(
        dof_handler_coarse, constraint_coarse, numbers::invalid_unsigned_int);

    this->vec_coarse.reinit(this->partitioner_coarse);
  }
  {
    // in case a DG space of order 0 is provided, DoFs indices are never defined
    // on element faces or vertices and therefore, the partitioner is fine
    IndexSet locally_relevant_dofs(dof_handler_fine.n_dofs());
    if (!this->fine_element_is_continuous &&
        dof_handler_fine.get_fe().degree != 0)
      locally_relevant_dofs.add_indices(global_dof_indices.begin(),
                                        global_dof_indices.end());

    this->partitioner_fine.reset(
      new Utilities::MPI::Partitioner(dof_handler_fine.locally_owned_dofs(),
                                      locally_relevant_dofs,
                                      dof_handler_fine.get_mpi_communicator()));

    this->vec_fine.reinit(this->partitioner_fine);
  }

  const auto &points = std::get<0>(points_ptrs_indices);

  // using level_dof_indices_fine_ptrs always works but in case of CG or DG
  // with degree==0 and n_components==1 support points to dof mapping is unique
  // and we don't need it.
  if (dof_handler_fine.get_fe().n_components() == 1 &&
      (this->fine_element_is_continuous ||
       dof_handler_fine.get_fe().degree == 0))
    this->level_dof_indices_fine_ptrs.clear();
  else
    this->level_dof_indices_fine_ptrs = std::get<1>(points_ptrs_indices);

  // fill level_dof_indices_fine with local indices
  this->level_dof_indices_fine.resize(global_dof_indices.size());
  for (unsigned int i = 0; i < global_dof_indices.size(); ++i)
    this->level_dof_indices_fine[i] =
      this->partitioner_fine->global_to_local(global_dof_indices[i]);

  // hand points over to RPE
  rpe.reinit(points, dof_handler_coarse.get_triangulation(), mapping_coarse);

  AssertThrow(
    !additional_data.enforce_all_points_found || rpe.all_points_found(),
    ExcMessage(
      "You requested that all points should be found, but this didn'thappen."
      " You can change this option through the AdditionalData struct in the constructor."));

  // set up MappingInfo for easier data access
  mapping_info = internal::fill_mapping_info<dim, Number>(rpe);

  // set up constraints
  const auto &cell_data = rpe.get_cell_data();

  constraint_info.reinit(dof_handler_coarse,
                         cell_data.cells.size(),
                         false /*TODO*/);

  for (unsigned int i = 0; i < cell_data.cells.size(); ++i)
    {
      typename DoFHandler<dim>::active_cell_iterator cell(
        &rpe.get_triangulation(),
        cell_data.cells[i].first,
        cell_data.cells[i].second,
        &dof_handler_coarse);

      constraint_info.read_dof_indices(i,
                                       numbers::invalid_unsigned_int,
                                       cell,
                                       constraint_coarse,
                                       this->partitioner_coarse);
    }

  constraint_info.finalize();

  const auto        &fe_base      = dof_handler_coarse.get_fe().base_element(0);
  const unsigned int n_components = dof_handler_coarse.get_fe().n_components();

  if (const auto fe = dynamic_cast<const FE_Q<dim> *>(&fe_base))
    fe_coarse = std::make_unique<FESystem<dim>>(FE_DGQ<dim>(fe->get_degree()),
                                                n_components);
  else if (const auto fe = dynamic_cast<const FE_SimplexP<dim> *>(&fe_base))
    fe_coarse =
      std::make_unique<FESystem<dim>>(FE_SimplexDGP<dim>(fe->get_degree()),
                                      n_components);
  else if (dynamic_cast<const FE_DGQ<dim> *>(&fe_base) ||
           dynamic_cast<const FE_SimplexP<dim> *>(&fe_base))
    fe_coarse = dof_handler_coarse.get_fe().clone();
  else
    AssertThrow(false, ExcMessage(dof_handler_coarse.get_fe().get_name()));
}


// Access utilities for scalar or vector valued fields computed with
// FEPointEvaluation. Mimicking the members of
// internal::FEPointEvaluation::EvaluatorTypeTraits
namespace internal
{
  template <typename T>
  const T &
  access(const T &value, const unsigned int)
  {
    return value;
  }

  template <typename T>
  T &
  access(T &value, const unsigned int)
  {
    return value;
  }

  template <int dim, typename T>
  const T &
  access(const Tensor<1, dim, T> &value, const unsigned int c)
  {
    return value[c];
  }

  template <int dim, typename T>
  T &
  access(Tensor<1, dim, T> &value, const unsigned int c)
  {
    return value[c];
  }
} // namespace internal



template <int dim, typename VectorType>
template <int n_components>
void
MGTwoLevelTransferNonNested<dim, VectorType>::prolongate_and_add_internal_comp(
  VectorType       &dst,
  const VectorType &src) const
{
  const auto evaluation_function = [&](auto &values, const auto &cell_data) {
    this->signals_non_nested.prolongation_cell_loop(true);
    std::vector<Number> solution_values;

    FEPointEvaluation<n_components, dim, dim, Number> evaluator(*mapping_info,
                                                                *fe_coarse);

    const auto &send_permutation = rpe.get_send_permutation();

    for (unsigned int cell = 0; cell < cell_data.cells.size(); ++cell)
      {
        solution_values.resize(fe_coarse->n_dofs_per_cell());

        // gather and resolve constraints
        internal::VectorReader<Number, VectorizedArrayType> reader;
        constraint_info.read_write_operation(
          reader,
          src,
          reinterpret_cast<VectorizedArrayType *>(solution_values.data()),
          cell,
          1,
          solution_values.size(),
          true);

        // evaluate and scatter
        evaluator.reinit(cell);

        evaluator.evaluate(solution_values, dealii::EvaluationFlags::values);

        for (const auto q : evaluator.quadrature_point_indices())
          {
            const unsigned int index =
              send_permutation[q + cell_data.reference_point_ptrs[cell]];

            for (unsigned int c = 0; c < n_components; ++c)
              values[index * n_components + c] =
                internal::access(evaluator.get_value(q), c);
          }
      }
    this->signals_non_nested.prolongation_cell_loop(false);
  };

  this->signals_non_nested.prolongation(true);

  std::vector<Number> &evaluation_point_results = this->rpe_input_output;
  std::vector<Number> &buffer                   = this->rpe_buffer;

  rpe.template evaluate_and_process<Number, n_components>(
    evaluation_point_results, buffer, evaluation_function, false);
  this->signals_non_nested.prolongation(false);

  // Keep a vector of typical inverse touch counts that avoid divisions, all
  // other cases are handled by a division in the code below.
  std::array<Number, 8> typical_weights;
  for (unsigned int i = 0; i < typical_weights.size(); ++i)
    typical_weights[i] = Number(1) / Number(i + 1);

  const bool  must_interpolate     = (rpe.is_map_unique() == false);
  const auto &ptr                  = rpe.get_point_ptrs();
  const auto &recv_permutation_inv = rpe.get_inverse_recv_permutation();
  for (unsigned int j = 0; j < ptr.size() - 1; ++j)
    {
      std::array<Number, n_components> result;
      // Weight operator in case some points are owned by multiple cells.
      if (must_interpolate)
        {
          const unsigned int n_entries = ptr[j + 1] - ptr[j];

          if (n_entries > 1)
            {
              result = {};
              for (unsigned int k = 0; k < n_entries; ++k)
                {
                  const unsigned int index = recv_permutation_inv[ptr[j] + k];

                  for (unsigned int c = 0; c < n_components; ++c)
                    result[c] +=
                      evaluation_point_results[index * n_components + c];
                }
              if (n_entries <= typical_weights.size())
                for (unsigned int c = 0; c < n_components; ++c)
                  result[c] *= typical_weights[n_entries - 1];
              else
                for (unsigned int c = 0; c < n_components; ++c)
                  result[c] *= Number(1) / Number(n_entries);
            }
          else if (n_entries == 1)
            {
              const unsigned int index = recv_permutation_inv[ptr[j]];

              for (unsigned int c = 0; c < n_components; ++c)
                result[c] = evaluation_point_results[index * n_components + c];
            }
          else
            result = {};
        }
      else
        for (unsigned int c = 0; c < n_components; ++c)
          result[c] = evaluation_point_results[j * n_components + c];

      if (level_dof_indices_fine_ptrs.empty())
        {
          for (unsigned int c = 0; c < n_components; ++c)
            {
              AssertIndexRange(n_components * j + c,
                               this->level_dof_indices_fine.size());
              dst.local_element(
                this->level_dof_indices_fine[n_components * j + c]) +=
                result[c];
            }
        }
      else
        {
          for (unsigned int i = this->level_dof_indices_fine_ptrs[j];
               i < this->level_dof_indices_fine_ptrs[j + 1];
               ++i)
            for (unsigned int c = 0; c < n_components; ++c)
              {
                AssertIndexRange(n_components * i + c,
                                 this->level_dof_indices_fine.size());
                dst.local_element(
                  this->level_dof_indices_fine[n_components * i + c]) +=
                  result[c];
              }
        }
    }
}



template <int dim, typename VectorType>
void
MGTwoLevelTransferNonNested<dim, VectorType>::prolongate_and_add_internal(
  VectorType       &dst,
  const VectorType &src) const
{
  if (this->fe_coarse->n_components() == 1)
    prolongate_and_add_internal_comp<1>(dst, src);
  else if (this->fe_coarse->n_components() == dim)
    prolongate_and_add_internal_comp<dim>(dst, src);
  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename VectorType>
template <int n_components>
void
MGTwoLevelTransferNonNested<dim, VectorType>::restrict_and_add_internal_comp(
  VectorType       &dst,
  const VectorType &src) const
{
  std::vector<Number> &evaluation_point_results = this->rpe_input_output;
  std::vector<Number> &buffer                   = this->rpe_buffer;

  std::array<Number, 8> typical_weights;
  for (unsigned int i = 0; i < typical_weights.size(); ++i)
    typical_weights[i] = Number(1) / Number(i + 1);

  const bool  must_interpolate     = (rpe.is_map_unique() == false);
  const auto &ptr                  = rpe.get_point_ptrs();
  const auto &recv_permutation_inv = rpe.get_inverse_recv_permutation();
  evaluation_point_results.resize(ptr.back() * n_components);

  for (unsigned int j = 0; j < ptr.size() - 1; ++j)
    {
      std::array<Number, n_components> result;
      if (level_dof_indices_fine_ptrs.empty())
        {
          for (unsigned int c = 0; c < n_components; ++c)
            {
              AssertIndexRange(n_components * j + c,
                               this->level_dof_indices_fine.size());

              result[c] = src.local_element(
                this->level_dof_indices_fine[n_components * j + c]);
            }
        }
      else
        {
          result = {};

          for (unsigned int i = this->level_dof_indices_fine_ptrs[j];
               i < this->level_dof_indices_fine_ptrs[j + 1];
               ++i)
            {
              for (unsigned int c = 0; c < n_components; ++c)
                {
                  AssertIndexRange(n_components * i + c,
                                   this->level_dof_indices_fine.size());
                  result[c] += src.local_element(
                    this->level_dof_indices_fine[n_components * i + c]);
                }
            }
        }
      if (must_interpolate)
        {
          const unsigned int n_entries = ptr[j + 1] - ptr[j];
          if (n_entries > 1)
            {
              if (n_entries <= typical_weights.size())
                for (unsigned int c = 0; c < n_components; ++c)
                  result[c] *= typical_weights[n_entries - 1];
              else
                for (unsigned int c = 0; c < n_components; ++c)
                  result[c] /= Number(n_entries);
            }
        }

      for (unsigned int i = ptr[j]; i < ptr[j + 1]; ++i)
        {
          const unsigned int index = recv_permutation_inv[i];

          for (unsigned int c = 0; c < n_components; ++c)
            evaluation_point_results[index * n_components + c] = result[c];
        }
    }

  const auto evaluation_function = [&](const auto &values,
                                       const auto &cell_data) {
    this->signals_non_nested.restriction_cell_loop(true);
    std::vector<Number>                               solution_values;
    FEPointEvaluation<n_components, dim, dim, Number> evaluator(*mapping_info,
                                                                *fe_coarse);

    const auto &send_permutation = rpe.get_send_permutation();

    for (unsigned int cell = 0; cell < cell_data.cells.size(); ++cell)
      {
        solution_values.resize(fe_coarse->n_dofs_per_cell());

        // gather and integrate
        evaluator.reinit(cell);

        for (const auto q : evaluator.quadrature_point_indices())
          {
            typename internal::FEPointEvaluation::
              EvaluatorTypeTraits<dim, dim, n_components, Number>::value_type
                value;

            const unsigned int index =
              send_permutation[q + cell_data.reference_point_ptrs[cell]];

            for (unsigned int c = 0; c < n_components; ++c)
              internal::access(value, c) = values[index * n_components + c];

            evaluator.submit_value(value, q);
          }

        evaluator.test_and_sum(solution_values, EvaluationFlags::values);

        // resolve constraints and scatter
        internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
          writer;
        constraint_info.read_write_operation(
          writer,
          dst,
          reinterpret_cast<VectorizedArrayType *>(solution_values.data()),
          cell,
          1,
          solution_values.size(),
          true);
      }
    this->signals_non_nested.restriction_cell_loop(false);
  };

  this->signals_non_nested.restriction(true);
  rpe.template process_and_evaluate<Number, n_components>(
    evaluation_point_results, buffer, evaluation_function, false);
  this->signals_non_nested.restriction(false);
}



template <int dim, typename VectorType>
void
MGTwoLevelTransferNonNested<dim, VectorType>::restrict_and_add_internal(
  VectorType       &dst,
  const VectorType &src) const
{
  if (this->fe_coarse->n_components() == 1)
    restrict_and_add_internal_comp<1>(dst, src);
  else if (this->fe_coarse->n_components() == dim)
    restrict_and_add_internal_comp<dim>(dst, src);
  else
    AssertThrow(false, ExcNotImplemented());
}



template <int dim, typename VectorType>
void
MGTwoLevelTransferNonNested<dim, VectorType>::interpolate(
  VectorType       &dst,
  const VectorType &src) const
{
  AssertThrow(false, ExcNotImplemented());
  (void)dst;
  (void)src;
}



template <int dim, typename VectorType>
std::pair<bool, bool>
MGTwoLevelTransferNonNested<dim, VectorType>::
  enable_inplace_operations_if_possible(
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_coarse,
    const std::shared_ptr<const Utilities::MPI::Partitioner>
      &external_partitioner_fine)
{
  return this->internal_enable_inplace_operations_if_possible(
    external_partitioner_coarse,
    external_partitioner_fine,
    this->vec_fine_needs_ghost_update,
    constraint_info,
    this->level_dof_indices_fine);
}



template <int dim, typename VectorType>
std::size_t
MGTwoLevelTransferNonNested<dim, VectorType>::memory_consumption() const
{
  std::size_t size = 0;

  size += this->partitioner_coarse->memory_consumption();
  size += this->vec_coarse.memory_consumption();
  size += MemoryConsumption::memory_consumption(this->level_dof_indices_fine);
  // TODO: add consumption for rpe, mapping_info and constraint_info.

  return size;
}



template <int dim, typename VectorType>
boost::signals2::connection
MGTwoLevelTransferNonNested<dim, VectorType>::connect_prolongation_cell_loop(
  const std::function<void(const bool)> &slot)
{
  return signals_non_nested.prolongation_cell_loop.connect(slot);
}



template <int dim, typename VectorType>
boost::signals2::connection
MGTwoLevelTransferNonNested<dim, VectorType>::connect_restriction_cell_loop(
  const std::function<void(const bool)> &slot)
{
  return signals_non_nested.restriction_cell_loop.connect(slot);
}



template <int dim, typename VectorType>
boost::signals2::connection
MGTwoLevelTransferNonNested<dim, VectorType>::connect_prolongation(
  const std::function<void(const bool)> &slot)
{
  return signals_non_nested.prolongation.connect(slot);
}



template <int dim, typename VectorType>
boost::signals2::connection
MGTwoLevelTransferNonNested<dim, VectorType>::connect_restriction(
  const std::function<void(const bool)> &slot)
{
  return signals_non_nested.restriction.connect(slot);
}


DEAL_II_NAMESPACE_CLOSE

#endif
