/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2010 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Guido Kanschat, Texas A&M University, 2009
 */


// The include files for the linear algebra: A regular SparseMatrix, which in
// turn will include the necessary files for SparsityPattern and Vector
// classes.
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/block_vector.h>

// Include files for setting up the mesh
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>

// Include files for FiniteElement classes and DoFHandler.
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>

// The include file for using MeshWorker::mesh_loop()
#include <deal.II/meshworker/mesh_loop.h>

// The include file for local integrators associated with the Laplacian
#include <deal.II/integrators/laplace.h>

// Support for multigrid methods
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>

// Finally, we take our exact solution from the library as well as quadrature
// and additional tools.
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <array>
#include <fstream>
#include <iostream>

// All classes of the deal.II library are in the namespace dealii. In order to
// save typing, we tell the compiler to search names in there as well.
namespace Step39
{
  using namespace dealii;

  // This is the function we use to set the boundary values and also the exact
  // solution we compare to.
  Functions::SlitSingularityFunction<2> exact_solution;

  // To begin with, let us recall how MeshWorker::mesh_loop() works. Much like
  // the WorkStream mechanism discussed in step-9, it separates the traversal
  // of the mesh from the computation performed on each cell, boundary face, or
  // interior face. The loop visits these objects one at a time, runs worker
  // functions that compute purely local contributions, and then hands the
  // results to a copier function that accumulates them into global matrices,
  // vectors, or error indicators.
  //
  // This organization is what makes parallel assembly practical: the worker
  // stage must not touch shared global data directly, and it should also avoid
  // repeatedly allocating expensive temporary objects. Consequently, each
  // worker receives a scratch object that owns reusable FEValues-like data and
  // other temporary arrays, together with a copy object that stores the local
  // contributions produced on the current cell or face. If you would like to
  // see the general WorkStream idea in a simpler setting, step-9 gives a more
  // introductory discussion before we apply the same pattern here to cells,
  // faces, and subfaces.

  // We start with a class that stores scratch data for assembling the global
  // and multigrid matrices. This is the most elaborate scratch object because
  // matrix assembly needs FEValues on cells, on boundary faces, on regular
  // interior faces, and on subfaces at refinement edges.
  template <int dim>
  struct MatrixScratchData
  {
    MatrixScratchData(const Mapping<dim>        &mapping,
                      const FiniteElement<dim>  &fe,
                      const Quadrature<dim>     &cell_quadrature,
                      const Quadrature<dim - 1> &face_quadrature,
                      const UpdateFlags          cell_update_flags,
                      const UpdateFlags          face_update_flags)
      : fe_values(mapping, fe, cell_quadrature, cell_update_flags)
      , boundary_fe_values(mapping, fe, face_quadrature, face_update_flags)
      , face_fe_values(mapping, fe, face_quadrature, face_update_flags)
      , subface_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_face_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_subface_values(mapping, fe, face_quadrature, face_update_flags)
    {}


    MatrixScratchData(const MatrixScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , boundary_fe_values(scratch_data.boundary_fe_values.get_mapping(),
                           scratch_data.boundary_fe_values.get_fe(),
                           scratch_data.boundary_fe_values.get_quadrature(),
                           scratch_data.boundary_fe_values.get_update_flags())
      , face_fe_values(scratch_data.face_fe_values.get_mapping(),
                       scratch_data.face_fe_values.get_fe(),
                       scratch_data.face_fe_values.get_quadrature(),
                       scratch_data.face_fe_values.get_update_flags())
      , subface_values(scratch_data.subface_values.get_mapping(),
                       scratch_data.subface_values.get_fe(),
                       scratch_data.subface_values.get_quadrature(),
                       scratch_data.subface_values.get_update_flags())
      , neighbor_face_values(
          scratch_data.neighbor_face_values.get_mapping(),
          scratch_data.neighbor_face_values.get_fe(),
          scratch_data.neighbor_face_values.get_quadrature(),
          scratch_data.neighbor_face_values.get_update_flags())
      , neighbor_subface_values(
          scratch_data.neighbor_subface_values.get_mapping(),
          scratch_data.neighbor_subface_values.get_fe(),
          scratch_data.neighbor_subface_values.get_quadrature(),
          scratch_data.neighbor_subface_values.get_update_flags())
    {}


    FEValues<dim>        fe_values;
    FEFaceValues<dim>    boundary_fe_values;
    FEFaceValues<dim>    face_fe_values;
    FESubfaceValues<dim> subface_values;
    FEFaceValues<dim>    neighbor_face_values;
    FESubfaceValues<dim> neighbor_subface_values;
  };


  // Next comes the scratch object used for assembling the right-hand side.
  // Here we only need access to boundary faces because the inhomogeneous terms
  // come from the Nitsche boundary contributions.
  template <int dim>
  struct RightHandSideScratchData
  {
    RightHandSideScratchData(const Mapping<dim>        &mapping,
                             const FiniteElement<dim>  &fe,
                             const Quadrature<dim - 1> &face_quadrature,
                             const UpdateFlags          face_update_flags)
      : boundary_fe_values(mapping, fe, face_quadrature, face_update_flags)
      , boundary_values(face_quadrature.size())
    {}


    RightHandSideScratchData(const RightHandSideScratchData<dim> &scratch_data)
      : boundary_fe_values(scratch_data.boundary_fe_values.get_mapping(),
                           scratch_data.boundary_fe_values.get_fe(),
                           scratch_data.boundary_fe_values.get_quadrature(),
                           scratch_data.boundary_fe_values.get_update_flags())
      , boundary_values(scratch_data.boundary_values.size())
    {}


    FEFaceValues<dim>   boundary_fe_values;
    std::vector<double> boundary_values;
  };


  // The error estimator, in turn, needs the following scratch class. Besides
  // FEValues-like objects on the relevant geometric entities, this object
  // stores reusable buffers for solution values, gradients, Hessians, and
  // exact boundary values so that the worker lambdas do not allocate these
  // arrays repeatedly.
  template <int dim>
  struct EstimatorScratchData
  {
    EstimatorScratchData(const Mapping<dim>        &mapping,
                         const FiniteElement<dim>  &fe,
                         const Quadrature<dim>     &cell_quadrature,
                         const Quadrature<dim - 1> &boundary_quadrature,
                         const Quadrature<dim - 1> &face_quadrature,
                         const UpdateFlags          cell_update_flags,
                         const UpdateFlags          boundary_update_flags,
                         const UpdateFlags          face_update_flags)
      : fe_values(mapping, fe, cell_quadrature, cell_update_flags)
      , boundary_fe_values(mapping,
                           fe,
                           boundary_quadrature,
                           boundary_update_flags)
      , face_fe_values(mapping, fe, face_quadrature, face_update_flags)
      , subface_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_face_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_subface_values(mapping, fe, face_quadrature, face_update_flags)
      , cell_hessians(cell_quadrature.size())
      , boundary_solution_values(boundary_quadrature.size())
      , boundary_exact_values(boundary_quadrature.size())
      , face_solution_values(face_quadrature.size())
      , neighbor_face_solution_values(face_quadrature.size())
      , face_solution_gradients(face_quadrature.size())
      , neighbor_face_solution_gradients(face_quadrature.size())
    {}


    EstimatorScratchData(const EstimatorScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , boundary_fe_values(scratch_data.boundary_fe_values.get_mapping(),
                           scratch_data.boundary_fe_values.get_fe(),
                           scratch_data.boundary_fe_values.get_quadrature(),
                           scratch_data.boundary_fe_values.get_update_flags())
      , face_fe_values(scratch_data.face_fe_values.get_mapping(),
                       scratch_data.face_fe_values.get_fe(),
                       scratch_data.face_fe_values.get_quadrature(),
                       scratch_data.face_fe_values.get_update_flags())
      , subface_values(scratch_data.subface_values.get_mapping(),
                       scratch_data.subface_values.get_fe(),
                       scratch_data.subface_values.get_quadrature(),
                       scratch_data.subface_values.get_update_flags())
      , neighbor_face_values(
          scratch_data.neighbor_face_values.get_mapping(),
          scratch_data.neighbor_face_values.get_fe(),
          scratch_data.neighbor_face_values.get_quadrature(),
          scratch_data.neighbor_face_values.get_update_flags())
      , neighbor_subface_values(
          scratch_data.neighbor_subface_values.get_mapping(),
          scratch_data.neighbor_subface_values.get_fe(),
          scratch_data.neighbor_subface_values.get_quadrature(),
          scratch_data.neighbor_subface_values.get_update_flags())
      , cell_hessians(scratch_data.cell_hessians.size())
      , boundary_solution_values(scratch_data.boundary_solution_values.size())
      , boundary_exact_values(scratch_data.boundary_exact_values.size())
      , face_solution_values(scratch_data.face_solution_values.size())
      , neighbor_face_solution_values(
          scratch_data.neighbor_face_solution_values.size())
      , face_solution_gradients(scratch_data.face_solution_gradients.size())
      , neighbor_face_solution_gradients(
          scratch_data.neighbor_face_solution_gradients.size())
    {}


    FEValues<dim>               fe_values;
    FEFaceValues<dim>           boundary_fe_values;
    FEFaceValues<dim>           face_fe_values;
    FESubfaceValues<dim>        subface_values;
    FEFaceValues<dim>           neighbor_face_values;
    FESubfaceValues<dim>        neighbor_subface_values;
    std::vector<Tensor<2, dim>> cell_hessians;
    std::vector<double>         boundary_solution_values;
    std::vector<double>         boundary_exact_values;
    std::vector<double>         face_solution_values;
    std::vector<double>         neighbor_face_solution_values;
    std::vector<Tensor<1, dim>> face_solution_gradients;
    std::vector<Tensor<1, dim>> neighbor_face_solution_gradients;
  };


  // To compute the actual error norms, we use another scratch class. As above,
  // we keep both FEValues-like objects and temporary storage for function
  // values and gradients in one per-thread object that can be reused for many
  // cells and faces.
  template <int dim>
  struct ErrorScratchData
  {
    ErrorScratchData(const Mapping<dim>        &mapping,
                     const FiniteElement<dim>  &fe,
                     const Quadrature<dim>     &cell_quadrature,
                     const Quadrature<dim - 1> &boundary_quadrature,
                     const Quadrature<dim - 1> &face_quadrature,
                     const UpdateFlags          cell_update_flags,
                     const UpdateFlags          boundary_update_flags,
                     const UpdateFlags          face_update_flags)
      : fe_values(mapping, fe, cell_quadrature, cell_update_flags)
      , boundary_fe_values(mapping,
                           fe,
                           boundary_quadrature,
                           boundary_update_flags)
      , face_fe_values(mapping, fe, face_quadrature, face_update_flags)
      , subface_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_face_values(mapping, fe, face_quadrature, face_update_flags)
      , neighbor_subface_values(mapping, fe, face_quadrature, face_update_flags)
      , cell_solution_values(cell_quadrature.size())
      , cell_solution_gradients(cell_quadrature.size())
      , cell_exact_values(cell_quadrature.size())
      , cell_exact_gradients(cell_quadrature.size())
      , boundary_solution_values(boundary_quadrature.size())
      , boundary_exact_values(boundary_quadrature.size())
      , face_solution_values(face_quadrature.size())
      , neighbor_face_solution_values(face_quadrature.size())
    {}


    ErrorScratchData(const ErrorScratchData<dim> &scratch_data)
      : fe_values(scratch_data.fe_values.get_mapping(),
                  scratch_data.fe_values.get_fe(),
                  scratch_data.fe_values.get_quadrature(),
                  scratch_data.fe_values.get_update_flags())
      , boundary_fe_values(scratch_data.boundary_fe_values.get_mapping(),
                           scratch_data.boundary_fe_values.get_fe(),
                           scratch_data.boundary_fe_values.get_quadrature(),
                           scratch_data.boundary_fe_values.get_update_flags())
      , face_fe_values(scratch_data.face_fe_values.get_mapping(),
                       scratch_data.face_fe_values.get_fe(),
                       scratch_data.face_fe_values.get_quadrature(),
                       scratch_data.face_fe_values.get_update_flags())
      , subface_values(scratch_data.subface_values.get_mapping(),
                       scratch_data.subface_values.get_fe(),
                       scratch_data.subface_values.get_quadrature(),
                       scratch_data.subface_values.get_update_flags())
      , neighbor_face_values(
          scratch_data.neighbor_face_values.get_mapping(),
          scratch_data.neighbor_face_values.get_fe(),
          scratch_data.neighbor_face_values.get_quadrature(),
          scratch_data.neighbor_face_values.get_update_flags())
      , neighbor_subface_values(
          scratch_data.neighbor_subface_values.get_mapping(),
          scratch_data.neighbor_subface_values.get_fe(),
          scratch_data.neighbor_subface_values.get_quadrature(),
          scratch_data.neighbor_subface_values.get_update_flags())
      , cell_solution_values(scratch_data.cell_solution_values.size())
      , cell_solution_gradients(scratch_data.cell_solution_gradients.size())
      , cell_exact_values(scratch_data.cell_exact_values.size())
      , cell_exact_gradients(scratch_data.cell_exact_gradients.size())
      , boundary_solution_values(scratch_data.boundary_solution_values.size())
      , boundary_exact_values(scratch_data.boundary_exact_values.size())
      , face_solution_values(scratch_data.face_solution_values.size())
      , neighbor_face_solution_values(
          scratch_data.neighbor_face_solution_values.size())
    {}


    FEValues<dim>               fe_values;
    FEFaceValues<dim>           boundary_fe_values;
    FEFaceValues<dim>           face_fe_values;
    FESubfaceValues<dim>        subface_values;
    FEFaceValues<dim>           neighbor_face_values;
    FESubfaceValues<dim>        neighbor_subface_values;
    std::vector<double>         cell_solution_values;
    std::vector<Tensor<1, dim>> cell_solution_gradients;
    std::vector<double>         cell_exact_values;
    std::vector<Tensor<1, dim>> cell_exact_gradients;
    std::vector<double>         boundary_solution_values;
    std::vector<double>         boundary_exact_values;
    std::vector<double>         face_solution_values;
    std::vector<double>         neighbor_face_solution_values;
  };


  // The first copy-data structure represents the local data associated with
  // one face contribution. These data are kept separately because each
  // interior face contributes four blocks coupling the two adjacent cells.
  struct FaceCopyData
  {
    FaceCopyData()
      : level_1(numbers::invalid_unsigned_int)
      , level_2(numbers::invalid_unsigned_int)
    {}

    unsigned int level_1;
    unsigned int level_2;

    FullMatrix<double> matrix_11;
    FullMatrix<double> matrix_12;
    FullMatrix<double> matrix_21;
    FullMatrix<double> matrix_22;

    std::vector<types::global_dof_index> dof_indices_1;
    std::vector<types::global_dof_index> dof_indices_2;
  };


  // Matrix assembly then uses the following copy-data class. Each worker first
  // fills the cell matrix and, if necessary, appends additional face
  // contributions; the copier then transfers all of this local data into the
  // global sparse matrices.
  template <int dim>
  struct MatrixCopyData
  {
    MatrixCopyData(const unsigned int dofs_per_cell = 0)
      : level(numbers::invalid_unsigned_int)
      , cell_matrix(dofs_per_cell, dofs_per_cell)
      , local_dof_indices(dofs_per_cell)
    {}


    template <typename CellIterator>
    void reinit(const CellIterator &cell, const bool use_level_dofs = false)
    {
      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

      level = cell->level();
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      if (use_level_dofs)
        cell->get_mg_dof_indices(local_dof_indices);
      else
        cell->get_dof_indices(local_dof_indices);

      face_data.clear();
    }


    template <typename CellIterator>
    FaceCopyData &emplace_face_data(const CellIterator &cell_1,
                                    const CellIterator &cell_2,
                                    const bool          use_level_dofs = false)
    {
      const unsigned int dofs_per_cell_1 = cell_1->get_fe().n_dofs_per_cell();
      const unsigned int dofs_per_cell_2 = cell_2->get_fe().n_dofs_per_cell();

      face_data.emplace_back();
      FaceCopyData &face_copy = face_data.back();

      face_copy.level_1 = cell_1->level();
      face_copy.level_2 = cell_2->level();
      face_copy.matrix_11.reinit(dofs_per_cell_1, dofs_per_cell_1);
      face_copy.matrix_12.reinit(dofs_per_cell_1, dofs_per_cell_2);
      face_copy.matrix_21.reinit(dofs_per_cell_2, dofs_per_cell_1);
      face_copy.matrix_22.reinit(dofs_per_cell_2, dofs_per_cell_2);
      face_copy.dof_indices_1.resize(dofs_per_cell_1);
      face_copy.dof_indices_2.resize(dofs_per_cell_2);

      if (use_level_dofs)
        {
          cell_1->get_mg_dof_indices(face_copy.dof_indices_1);
          cell_2->get_mg_dof_indices(face_copy.dof_indices_2);
        }
      else
        {
          cell_1->get_dof_indices(face_copy.dof_indices_1);
          cell_2->get_dof_indices(face_copy.dof_indices_2);
        }

      return face_copy;
    }


    unsigned int level;

    FullMatrix<double>                   cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<FaceCopyData>            face_data;
  };


  // For the right-hand side, a simpler copy-data class is sufficient. In this
  // case the local result is just one cell vector together with the
  // corresponding global DoF indices.
  struct RightHandSideCopyData
  {
    RightHandSideCopyData(const unsigned int dofs_per_cell = 0)
      : cell_rhs(dofs_per_cell)
      , local_dof_indices(dofs_per_cell)
    {}


    template <typename CellIterator>
    void reinit(const CellIterator &cell)
    {
      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();

      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
    }


    Vector<double>                       cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;
  };


  // Finally, the estimator and error computations share the following
  // copy-data class. Instead of matrices or vectors, the local results are a
  // small fixed number of scalar indicators attached to the current cell and
  // to the faces touching it. The template argument is the number of such
  // indicators, which is known at compile time. Specifically, we will use
  // this class below with `n_values=1` for the error estimator, and with
  // `n_values=2` for the error computation (where we compute both the $L_2$
  // and $H^1$ errors).
  template <unsigned int n_values>
  struct ErrorCopyData
  {
    struct FaceContribution
    {
      FaceContribution()
        : cell_index_1(numbers::invalid_unsigned_int)
        , cell_index_2(numbers::invalid_unsigned_int)
      {
        values.fill(0.);
      }

      unsigned int                 cell_index_1;
      unsigned int                 cell_index_2;
      std::array<double, n_values> values;
    };


    ErrorCopyData()
      : cell_index(numbers::invalid_unsigned_int)
    {
      cell_values.fill(0.);
    }


    template <typename CellIterator>
    void reinit(const CellIterator &cell)
    {
      cell_index = cell->active_cell_index();
      cell_values.fill(0.);
      face_data.clear();
    }


    template <typename CellIterator>
    FaceContribution &emplace_face_data(const CellIterator &cell_1,
                                        const CellIterator &cell_2)
    {
      face_data.emplace_back();
      FaceContribution &face_contribution = face_data.back();
      face_contribution.cell_index_1      = cell_1->active_cell_index();
      face_contribution.cell_index_2      = cell_2->active_cell_index();
      face_contribution.values.fill(0.);
      return face_contribution;
    }


    unsigned int                  cell_index;
    std::array<double, n_values>  cell_values;
    std::vector<FaceContribution> face_data;
  };

  // @sect3{The local integrators}

  // The MeshWorker::mesh_loop() function separates the local integration on
  // cells and faces from the traversal of the mesh. In the functions below, we
  // therefore use a scratch object holding the FEValues-like data needed on
  // the current cell or face, together with copy objects that store local
  // matrices, vectors, or per-cell indicators before they are transferred to
  // the global objects by a copier.

  // The first namespace defining local integrators is responsible for
  // assembling the global matrix as well as the level matrices.
  // On each cell, we integrate the Dirichlet form as well as the
  // Nitsche boundary conditions and the interior penalty fluxes between
  // cells.
  //
  // The boundary and flux terms need a penalty parameter, which should be
  // adjusted to the cell size and the polynomial degree. We compute it
  // in two steps: First, we compute on each cell
  // $K_i$ the value $P_i = p_i(p_i+1)/h_i$, where
  // $p_i$ is the polynomial degree on cell $K_i$ and $h_i$ is the length of
  // $K_i$ orthogonal to the current face. Second, if exactly one of the two
  // cells adjacent to the face has children, its penalty is multiplied
  // by two (to account for the fact that the mesh size $h_i$ there is
  // only half that previously computed); it is possible that both adjacent
  // cells are refined, in which case we are integrating over a non-active
  // face and no adjustment is necessary. Finally, we return the average
  // of the two penalty values.
  namespace MatrixIntegrator
  {
    template <int dim, typename CellIterator>
    double ip_penalty_factor(const CellIterator &cell1,
                             const unsigned int  face1,
                             const unsigned int  deg1,
                             const CellIterator &cell2,
                             const unsigned int  face2,
                             const unsigned int  deg2)
    {
      const unsigned int normal1 =
        GeometryInfo<dim>::unit_normal_direction[face1];
      const unsigned int normal2 =
        GeometryInfo<dim>::unit_normal_direction[face2];
      const unsigned int deg1sq = (deg1 == 0) ? 1 : deg1 * (deg1 + 1);
      const unsigned int deg2sq = (deg2 == 0) ? 1 : deg2 * (deg2 + 1);

      double penalty1 = deg1sq / cell1->extent_in_direction(normal1);
      double penalty2 = deg2sq / cell2->extent_in_direction(normal2);
      if (cell1->has_children() && !cell2->has_children())
        penalty1 *= 2;
      else if (!cell1->has_children() && cell2->has_children())
        penalty2 *= 2;

      const double penalty = 0.5 * (penalty1 + penalty2);
      return penalty;
    }


    template <int dim>
    void cell(const FEValues<dim> &fe_values, FullMatrix<double> &M)
    {
      for (unsigned int k = 0; k < fe_values.n_quadrature_points; ++k)
        {
          const double dx = fe_values.JxW(k);

          for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i)
            {
              const double Mii =
                fe_values.shape_grad(i, k) * fe_values.shape_grad(i, k) * dx;

              M(i, i) += Mii;

              for (unsigned int j = i + 1; j < fe_values.dofs_per_cell; ++j)
                {
                  const double Mij = fe_values.shape_grad(j, k) *
                                     fe_values.shape_grad(i, k) * dx;

                  M(i, j) += Mij;
                  M(j, i) += Mij;
                }
            }
        }
    }


    // Boundary faces use the Nitsche method to impose boundary values:
    template <int dim, typename CellIterator>
    void boundary(const CellIterator      &cell,
                  const unsigned int       face_no,
                  const FEFaceValues<dim> &fe_face_values,
                  FullMatrix<double>      &M)
    {
      AssertDimension(M.n(), fe_face_values.dofs_per_cell);
      AssertDimension(M.m(), fe_face_values.dofs_per_cell);

      const unsigned int polynomial_degree =
        fe_face_values.get_fe().tensor_degree();

      const double ip_penalty = ip_penalty_factor<dim>(
        cell, face_no, polynomial_degree, cell, face_no, polynomial_degree);

      for (unsigned int k = 0; k < fe_face_values.n_quadrature_points; ++k)
        {
          const double          dx = fe_face_values.JxW(k);
          const Tensor<1, dim> &n  = fe_face_values.normal_vector(k);

          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i)
            for (unsigned int j = 0; j < fe_face_values.dofs_per_cell; ++j)
              M(i, j) += (2. * fe_face_values.shape_value(i, k) * ip_penalty *
                            fe_face_values.shape_value(j, k) -
                          (n * fe_face_values.shape_grad(i, k)) *
                            fe_face_values.shape_value(j, k) -
                          (n * fe_face_values.shape_grad(j, k)) *
                            fe_face_values.shape_value(i, k)) *
                         dx;
        }
    }

    // Interior faces use the interior penalty method:
    template <int dim, typename CellIterator>
    void face(const CellIterator &cell_1,
              const unsigned int  face_no_1,
              const FEFaceValuesBase<dim>
                &fe_face_values_1, // These can be FEFaceValues or
                                   // FESubfaceValues objects.
              const CellIterator          &cell_2,
              const unsigned int           face_no_2,
              const FEFaceValuesBase<dim> &fe_face_values_2, // Same here
              FullMatrix<double>          &M11,
              FullMatrix<double>          &M12,
              FullMatrix<double>          &M21,
              FullMatrix<double>          &M22)
    {
      AssertDimension(M11.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M11.m(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M12.n(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M12.m(), fe_face_values_2.dofs_per_cell);
      AssertDimension(M21.n(), fe_face_values_2.dofs_per_cell);
      AssertDimension(M21.m(), fe_face_values_1.dofs_per_cell);
      AssertDimension(M22.n(), fe_face_values_2.dofs_per_cell);
      AssertDimension(M22.m(), fe_face_values_2.dofs_per_cell);

      const unsigned int polynomial_degree =
        fe_face_values_1.get_fe().tensor_degree();
      const double ip_penalty = ip_penalty_factor<dim>(cell_1,
                                                       face_no_1,
                                                       polynomial_degree,
                                                       cell_2,
                                                       face_no_2,
                                                       polynomial_degree);

      const double nui = 1.;
      const double nue = 1.;
      const double nu  = .5 * (nui + nue);

      for (unsigned int k = 0; k < fe_face_values_1.n_quadrature_points; ++k)
        {
          const double          dx = fe_face_values_1.JxW(k);
          const Tensor<1, dim> &n  = fe_face_values_1.normal_vector(k);

          for (unsigned int i = 0; i < fe_face_values_1.dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < fe_face_values_1.dofs_per_cell; ++j)
                {
                  const double vi   = fe_face_values_1.shape_value(i, k);
                  const double dnvi = n * fe_face_values_1.shape_grad(i, k);
                  const double ve   = fe_face_values_2.shape_value(i, k);
                  const double dnve = n * fe_face_values_2.shape_grad(i, k);
                  const double ui   = fe_face_values_1.shape_value(j, k);
                  const double dnui = n * fe_face_values_1.shape_grad(j, k);
                  const double ue   = fe_face_values_2.shape_value(j, k);
                  const double dnue = n * fe_face_values_2.shape_grad(j, k);

                  M11(i, j) += (-.5 * nui * dnvi * ui - .5 * nui * dnui * vi +
                                nu * ip_penalty * ui * vi) *
                               dx;
                  M12(i, j) += (.5 * nui * dnvi * ue - .5 * nue * dnue * vi -
                                nu * ip_penalty * vi * ue) *
                               dx;
                  M21(i, j) += (-.5 * nue * dnve * ui + .5 * nui * dnui * ve -
                                nu * ip_penalty * ui * ve) *
                               dx;
                  M22(i, j) += (.5 * nue * dnve * ue + .5 * nue * dnue * ve +
                                nu * ip_penalty * ue * ve) *
                               dx;
                }
            }
        }
    }
  } // namespace MatrixIntegrator

  // The second set of local integrators builds the right hand side. In our
  // example, the right hand side function is zero, such that only the boundary
  // condition is set here in weak form.
  namespace RHSIntegrator
  {
    template <int dim, typename CellIterator>
    void boundary(const CellIterator        &cell,
                  const unsigned int         face_no,
                  const FEFaceValues<dim>   &fe,
                  const std::vector<double> &boundary_values,
                  Vector<double>            &local_vector)
    {
      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             cell->face(face_no)->measure() / cell->measure();

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
          local_vector(i) +=
            (-penalty * fe.shape_value(i, k) // (-sigma * v_i(x_k)
             +
             fe.normal_vector(k) * fe.shape_grad(i, k)) // + n * grad v_i(x_k))
            * boundary_values[k] * fe.JxW(k);           // u^D(x_k) * dx
    }


  } // namespace RHSIntegrator

  // The third local integrator is responsible for the contributions to the
  // error estimate. This is the standard energy estimator due to Karakashian
  // and Pascal (2003).
  // The cell contribution is the Laplacian of the discrete solution, since
  // the right hand side is zero.
  namespace Estimator
  {
    template <int dim, typename CellIterator>
    double cell(const CellIterator                &cell,
                const FEValues<dim>               &fe,
                const std::vector<Tensor<2, dim>> &DDuh)
    {
      double value = 0.;
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double t = cell->diameter() * trace(DDuh[k]);
          value += t * t * fe.JxW(k);
        }
      return std::sqrt(value);
    }

    // At the boundary, we use simply a weighted form of the boundary residual,
    // namely the norm of the difference between the finite element solution and
    // the correct boundary condition.
    template <int dim, typename CellIterator>
    double boundary(const CellIterator        &cell,
                    const unsigned int         face_no,
                    const FEFaceValues<dim>   &fe,
                    const std::vector<double> &uh,
                    const std::vector<double> &boundary_values)
    {
      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             cell->face(face_no)->measure() / cell->measure();

      double value = 0.;
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = boundary_values[k] - uh[k];
          value += penalty * diff * diff * fe.JxW(k);
        }
      return std::sqrt(value);
    }


    // Finally, on interior faces, the estimator consists of the jumps of the
    // solution and its normal derivative, weighted appropriately.
    template <int dim, typename CellIterator>
    double face(const CellIterator &cell_1,
                const unsigned int  face_no_1,
                const FEFaceValuesBase<dim>
                  &fe, // This can be an FEFaceValues or FESubfaceValues object.
                const std::vector<double>         &uh1,
                const std::vector<Tensor<1, dim>> &Duh1,
                const CellIterator                &cell_2,
                const unsigned int                 face_no_2,
                const std::vector<double>         &uh2,
                const std::vector<Tensor<1, dim>> &Duh2)
    {
      const unsigned int degree   = fe.get_fe().tensor_degree();
      const double       penalty1 = degree * (degree + 1) *
                              cell_1->face(face_no_1)->measure() /
                              cell_1->measure();
      const double penalty2 = degree * (degree + 1) *
                              cell_2->face(face_no_2)->measure() /
                              cell_2->measure();
      const double penalty = penalty1 + penalty2;
      const double h       = cell_1->face(face_no_1)->measure();

      double value = 0.;
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff1 = uh1[k] - uh2[k];
          const double diff2 =
            fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
          value += (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k);
        }
      return std::sqrt(value);
    }
  } // namespace Estimator

  // Finally we have an integrator for the error. Since the energy norm for
  // discontinuous Galerkin problems not only involves the difference of the
  // gradient inside the cells, but also the jump terms across faces and at
  // the boundary, we cannot just use VectorTools::integrate_difference().
  // Instead, we use MeshWorker::mesh_loop() to compute the error ourselves.

  // There are several different ways to define this energy norm, but all of
  // them are equivalent to each other uniformly with mesh size (some not
  // uniformly with polynomial degree). Here, we choose @f[ \|u\|_{1,h} =
  // \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
  // 4\sigma_F\|\average{ u \mathbf n}\|^2_F + \sum_{F \in F_h^b}
  // 2\sigma_F\|u\|^2_F @f]
  //
  // Below, the first function is, as always, the integration on cells. The
  // exact solution is evaluated directly in the quadrature points and then
  // compared against the discrete solution.

  // Additionally to computing the error in the energy norm, we also compute
  // the <i>L<sup>2</sup></i>-error in the same loop. Obviously, this one does
  // not have any jump terms and only appears in the integration on cells.

  namespace ErrorIntegrator
  {
    template <int dim>
    std::array<double, 2>
    cell(const FEValues<dim>               &fe,
         const std::vector<double>         &uh,
         const std::vector<Tensor<1, dim>> &Duh,
         const std::vector<double>         &exact_values,
         const std::vector<Tensor<1, dim>> &exact_gradients)
    {
      std::array<double, 2> values = {};

      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          double sum = 0;
          for (unsigned int d = 0; d < dim; ++d)
            {
              const double diff = exact_gradients[k][d] - Duh[k][d];
              sum += diff * diff;
            }
          const double diff = exact_values[k] - uh[k];
          values[0] += sum * fe.JxW(k);
          values[1] += diff * diff * fe.JxW(k);
        }
      values[0] = std::sqrt(values[0]);
      values[1] = std::sqrt(values[1]);
      return values;
    }


    template <int dim, typename CellIterator>
    double boundary(const CellIterator        &cell,
                    const unsigned int         face_no,
                    const FEFaceValues<dim>   &fe,
                    const std::vector<double> &uh,
                    const std::vector<double> &exact_values)
    {
      const unsigned int degree  = fe.get_fe().tensor_degree();
      const double       penalty = 2. * degree * (degree + 1) *
                             cell->face(face_no)->measure() / cell->measure();

      double value = 0.;
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = exact_values[k] - uh[k];
          value += penalty * diff * diff * fe.JxW(k);
        }
      return std::sqrt(value);
    }


    template <int dim, typename CellIterator>
    double face(const CellIterator &cell_1,
                const unsigned int  face_no_1,
                const FEFaceValuesBase<dim>
                  &fe, // This can be an FEFaceValues or FESubfaceValues object.
                const std::vector<double> &uh1,
                const CellIterator        &cell_2,
                const unsigned int         face_no_2,
                const std::vector<double> &uh2)
    {
      const unsigned int degree   = fe.get_fe().tensor_degree();
      const double       penalty1 = degree * (degree + 1) *
                              cell_1->face(face_no_1)->measure() /
                              cell_1->measure();
      const double penalty2 = degree * (degree + 1) *
                              cell_2->face(face_no_2)->measure() /
                              cell_2->measure();
      const double penalty = penalty1 + penalty2;

      double value = 0.;
      for (unsigned k = 0; k < fe.n_quadrature_points; ++k)
        {
          const double diff = uh1[k] - uh2[k];
          value += penalty * diff * diff * fe.JxW(k);
        }
      return std::sqrt(value);
    }
  } // namespace ErrorIntegrator


  // @sect3{The main class}

  // This class does the main job, like in previous examples. For a
  // description of the functions declared here, please refer to the
  // implementation below.
  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    InteriorPenaltyProblem();

    void run(unsigned int n_steps);

  private:
    void   setup_system();
    void   assemble_matrix();
    void   assemble_mg_matrix();
    void   assemble_right_hand_side();
    void   error();
    double estimate();
    void   solve();
    void   output_results(const unsigned int cycle) const;

    // The member objects related to the discretization are here.
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;
    const FE_DGQ<2>      fe;
    DoFHandler<dim>      dof_handler;

    // Then, we have the matrices and vectors related to the global discrete
    // system.
    SparsityPattern      sparsity;
    SparseMatrix<double> matrix;
    Vector<double>       solution;
    Vector<double>       right_hand_side;
    BlockVector<double>  estimates;

    // Finally, we have a group of sparsity patterns and sparse matrices
    // related to the multilevel preconditioner.  First, we have a level
    // matrix and its sparsity pattern.
    MGLevelObject<SparsityPattern>      mg_sparsity;
    MGLevelObject<SparseMatrix<double>> mg_matrix;

    // When we perform multigrid with local smoothing on locally refined
    // meshes, additional matrices are required; see Kanschat (2004). Here is
    // the sparsity pattern for these edge matrices. We only need one, because
    // the pattern of the up matrix is the transpose of that of the down
    // matrix. These matrices are filled in the multigrid mesh loop below.
    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
    // The flux matrix at the refinement edge, coupling fine level degrees of
    // freedom to coarse level.
    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_down;
    // The transpose of the flux matrix at the refinement edge, coupling
    // coarse level degrees of freedom to fine level.
    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_up;
  };


  // The constructor simply sets up the coarse grid and the DoFHandler.
  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem()
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
    , mapping()
    , fe(3)
    , dof_handler(triangulation)
    , estimates(1)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  // In this function, we set up the dimension of the linear system and the
  // sparsity patterns for the global matrix as well as the level matrices.
  template <int dim>
  void InteriorPenaltyProblem<dim>::setup_system()
  {
    // First, we use the finite element to distribute degrees of freedom over
    // the mesh and number them.
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs();
    unsigned int n_dofs = dof_handler.n_dofs();
    // Then, we already know the size of the vectors representing finite
    // element functions.
    solution.reinit(n_dofs);
    right_hand_side.reinit(n_dofs);

    // Next, we set up the sparsity pattern for the global matrix. Since we do
    // not know the row sizes in advance, we first fill a temporary
    // DynamicSparsityPattern object and copy it to the regular
    // SparsityPattern once it is complete.
    DynamicSparsityPattern dsp(n_dofs);
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
    sparsity.copy_from(dsp);
    matrix.reinit(sparsity);

    const unsigned int n_levels = triangulation.n_levels();
    // The global system is set up, now we attend to the level matrices. We
    // resize all matrix objects to hold one matrix per level.
    mg_matrix.resize(0, n_levels - 1);
    mg_matrix.clear_elements();
    mg_matrix_dg_up.resize(0, n_levels - 1);
    mg_matrix_dg_up.clear_elements();
    mg_matrix_dg_down.resize(0, n_levels - 1);
    mg_matrix_dg_down.clear_elements();
    // It is important to update the sparsity patterns after <tt>clear()</tt>
    // was called for the level matrices, since the matrices lock the sparsity
    // pattern through the ObserverPointer and
    // EnableObserverPointer mechanism.
    mg_sparsity.resize(0, n_levels - 1);
    mg_sparsity_dg_interface.resize(0, n_levels - 1);

    // Now all objects are prepared to hold one sparsity pattern or matrix per
    // level. What's left is setting up the sparsity patterns on each level.
    for (unsigned int level = mg_sparsity.min_level();
         level <= mg_sparsity.max_level();
         ++level)
      {
        // These are roughly the same lines as above for the global matrix,
        // now for each level.
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level);
        mg_sparsity[level].copy_from(dsp);
        mg_matrix[level].reinit(mg_sparsity[level]);

        // Additionally, we need to initialize the transfer matrices at the
        // refinement edge between levels. They are stored at the index
        // referring to the finer of the two indices, thus there is no such
        // object on level 0.
        if (level > 0)
          {
            DynamicSparsityPattern dsp;
            dsp.reinit(dof_handler.n_dofs(level - 1),
                       dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, dsp, level);
            mg_sparsity_dg_interface[level].copy_from(dsp);
            mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
            mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
          }
      }
  }


  // In this function, we assemble the global system matrix, where by global we
  // indicate that this is the matrix of the discrete system we solve and it
  // is covering the whole mesh.
  //
  // The call to MeshWorker::mesh_loop() below follows the usual WorkStream
  // pattern. It traverses all active cells and, depending on the flags we pass
  // in, also visits boundary faces and interior faces. For each such geometric
  // object, a worker function computes the corresponding local contribution:
  // the cell worker fills the cell matrix, the boundary worker adds the
  // Nitsche terms on Dirichlet boundaries, and the interior face worker
  // computes the four coupling blocks associated with one face.
  //
  // None of these workers writes into the global sparse matrix directly.
  // Instead, they only populate the copy object associated with the current
  // cell. Once the local work is complete, the copier takes the data stored in
  // that copy object and inserts it into the global matrix. This separation is
  // what allows MeshWorker::mesh_loop() to parallelize the local integration
  // safely.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_matrix()
  {
    using CellIterator = typename DoFHandler<dim>::active_cell_iterator;

    const MatrixScratchData<dim> scratch(mapping,
                                         fe,
                                         QGauss<dim>(fe.degree + 1),
                                         QGauss<dim - 1>(fe.degree + 1),
                                         update_gradients | update_JxW_values,
                                         update_values | update_gradients |
                                           update_JxW_values |
                                           update_normal_vectors);
    const MatrixCopyData<dim>    copy_data(fe.n_dofs_per_cell());

    MeshWorker::mesh_loop(
      dof_handler.begin_active(),
      dof_handler.end(),
      /* cell worker: */
      [&](const CellIterator     &cell,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        copy.reinit(cell);
        scratch_data.fe_values.reinit(cell);
        MatrixIntegrator::cell<dim>(scratch_data.fe_values, copy.cell_matrix);
      },
      /* copier: */
      [&](const MatrixCopyData<dim> &copy) {
        matrix.add(copy.local_dof_indices, copy.cell_matrix);
        for (const auto &face : copy.face_data)
          {
            matrix.add(face.dof_indices_1, face.dof_indices_1, face.matrix_11);
            matrix.add(face.dof_indices_1, face.dof_indices_2, face.matrix_12);
            matrix.add(face.dof_indices_2, face.dof_indices_1, face.matrix_21);
            matrix.add(face.dof_indices_2, face.dof_indices_2, face.matrix_22);
          }
      },
      /* scratch and copy objects: */
      scratch,
      copy_data,
      /* where and how we want to integrate: */
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces |
        MeshWorker::assemble_own_interior_faces_once,
      /* boundary face worker: */
      [&](const CellIterator     &cell,
          const unsigned int      face_no,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        scratch_data.boundary_fe_values.reinit(cell, face_no);
        MatrixIntegrator::boundary<dim>(cell,
                                        face_no,
                                        scratch_data.boundary_fe_values,
                                        copy.cell_matrix);
      },
      /* interior face worker: */
      [&](const CellIterator     &cell,
          const unsigned int      face_no,
          const unsigned int      subface_no,
          const CellIterator     &neighbor,
          const unsigned int      neighbor_face_no,
          const unsigned int      neighbor_subface_no,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        FaceCopyData &face_copy = copy.emplace_face_data(cell, neighbor);
        if (subface_no == numbers::invalid_unsigned_int)
          {
            scratch_data.face_fe_values.reinit(cell, face_no);
            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                MatrixIntegrator::face<dim>(cell,
                                            face_no,
                                            scratch_data.face_fe_values,
                                            neighbor,
                                            neighbor_face_no,
                                            scratch_data.neighbor_face_values,
                                            face_copy.matrix_11,
                                            face_copy.matrix_12,
                                            face_copy.matrix_21,
                                            face_copy.matrix_22);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                MatrixIntegrator::face<dim>(
                  cell,
                  face_no,
                  scratch_data.face_fe_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_subface_values,
                  face_copy.matrix_11,
                  face_copy.matrix_12,
                  face_copy.matrix_21,
                  face_copy.matrix_22);
              }
          }
        else
          {
            scratch_data.subface_values.reinit(cell, face_no, subface_no);
            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                MatrixIntegrator::face<dim>(cell,
                                            face_no,
                                            scratch_data.subface_values,
                                            neighbor,
                                            neighbor_face_no,
                                            scratch_data.neighbor_face_values,
                                            face_copy.matrix_11,
                                            face_copy.matrix_12,
                                            face_copy.matrix_21,
                                            face_copy.matrix_22);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                MatrixIntegrator::face<dim>(
                  cell,
                  face_no,
                  scratch_data.subface_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_subface_values,
                  face_copy.matrix_11,
                  face_copy.matrix_12,
                  face_copy.matrix_21,
                  face_copy.matrix_22);
              }
          }
      });
  }


  // Now, we do the same for the level matrices. Not too surprisingly, this
  // function looks like a twin of the previous one. The mesh loop again
  // traverses cells, boundary faces, and interior faces, and the worker
  // functions compute local matrix contributions without touching global data.
  // The scratch object provides the FEValues-like data needed for those local
  // computations, and the copy object collects the local cell matrix together
  // with the face blocks.
  //
  // The main difference lies in the copier. Rather than assembling into a
  // single global matrix, it dispatches the local data into the appropriate
  // level matrix and, on refinement edges, into the up- and down-transfer
  // matrices that are required by the multigrid algorithm.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_mg_matrix()
  {
    using CellIterator = typename DoFHandler<dim>::level_cell_iterator;

    const MatrixScratchData<dim> scratch(mapping,
                                         fe,
                                         QGauss<dim>(fe.degree + 1),
                                         QGauss<dim - 1>(fe.degree + 1),
                                         update_gradients | update_JxW_values,
                                         update_values | update_gradients |
                                           update_JxW_values |
                                           update_normal_vectors);
    const MatrixCopyData<dim>    copy_data(fe.n_dofs_per_cell());

    MeshWorker::mesh_loop(
      dof_handler.begin_mg(),
      dof_handler.end_mg(),
      /* cell worker: */
      [&](const CellIterator     &cell,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        copy.reinit(cell, true);
        scratch_data.fe_values.reinit(cell);
        MatrixIntegrator::cell<dim>(scratch_data.fe_values, copy.cell_matrix);
      },
      /* copier: */
      [&](const MatrixCopyData<dim> &copy) {
        mg_matrix[copy.level].add(copy.local_dof_indices, copy.cell_matrix);

        for (const auto &face : copy.face_data)
          if (face.level_1 == face.level_2)
            {
              mg_matrix[face.level_1].add(face.dof_indices_1,
                                          face.dof_indices_1,
                                          face.matrix_11);
              mg_matrix[face.level_1].add(face.dof_indices_1,
                                          face.dof_indices_2,
                                          face.matrix_12);
              mg_matrix[face.level_1].add(face.dof_indices_2,
                                          face.dof_indices_1,
                                          face.matrix_21);
              mg_matrix[face.level_1].add(face.dof_indices_2,
                                          face.dof_indices_2,
                                          face.matrix_22);
            }
          else
            {
              Assert(face.level_1 > face.level_2, ExcInternalError());

              mg_matrix[face.level_1].add(face.dof_indices_1,
                                          face.dof_indices_1,
                                          face.matrix_11);

              for (unsigned int j = 0; j < face.dof_indices_2.size(); ++j)
                for (unsigned int k = 0; k < face.dof_indices_1.size(); ++k)
                  {
                    mg_matrix_dg_up[face.level_1].add(face.dof_indices_2[j],
                                                      face.dof_indices_1[k],
                                                      face.matrix_12(k, j));
                    mg_matrix_dg_down[face.level_1].add(face.dof_indices_2[j],
                                                        face.dof_indices_1[k],
                                                        face.matrix_21(j, k));
                  }
            }
      },
      /* scratch and copy objects: */
      scratch,
      copy_data,
      /* where and how we want to integrate: */
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces |
        MeshWorker::assemble_own_interior_faces_once,
      /* boundary face worker: */
      [&](const CellIterator     &cell,
          const unsigned int      face_no,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        scratch_data.boundary_fe_values.reinit(cell, face_no);
        MatrixIntegrator::boundary<dim>(cell,
                                        face_no,
                                        scratch_data.boundary_fe_values,
                                        copy.cell_matrix);
      },
      /* interior face worker: */
      [&](const CellIterator     &cell,
          const unsigned int      face_no,
          const unsigned int      subface_no,
          const CellIterator     &neighbor,
          const unsigned int      neighbor_face_no,
          const unsigned int      neighbor_subface_no,
          MatrixScratchData<dim> &scratch_data,
          MatrixCopyData<dim>    &copy) {
        FaceCopyData &face_copy = copy.emplace_face_data(cell, neighbor, true);
        if (subface_no == numbers::invalid_unsigned_int)
          {
            scratch_data.face_fe_values.reinit(cell, face_no);
            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                MatrixIntegrator::face<dim>(cell,
                                            face_no,
                                            scratch_data.face_fe_values,
                                            neighbor,
                                            neighbor_face_no,
                                            scratch_data.neighbor_face_values,
                                            face_copy.matrix_11,
                                            face_copy.matrix_12,
                                            face_copy.matrix_21,
                                            face_copy.matrix_22);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                MatrixIntegrator::face<dim>(
                  cell,
                  face_no,
                  scratch_data.face_fe_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_subface_values,
                  face_copy.matrix_11,
                  face_copy.matrix_12,
                  face_copy.matrix_21,
                  face_copy.matrix_22);
              }
          }
        else
          {
            scratch_data.subface_values.reinit(cell, face_no, subface_no);
            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                MatrixIntegrator::face<dim>(cell,
                                            face_no,
                                            scratch_data.subface_values,
                                            neighbor,
                                            neighbor_face_no,
                                            scratch_data.neighbor_face_values,
                                            face_copy.matrix_11,
                                            face_copy.matrix_12,
                                            face_copy.matrix_21,
                                            face_copy.matrix_22);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                MatrixIntegrator::face<dim>(
                  cell,
                  face_no,
                  scratch_data.subface_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_subface_values,
                  face_copy.matrix_11,
                  face_copy.matrix_12,
                  face_copy.matrix_21,
                  face_copy.matrix_22);
              }
          }
      });
  }


  // Here we have another clone of the assembly function. The difference to
  // assembling the system matrix consists in that we assemble a vector here.
  //
  // The mesh loop still uses the same worker/copier split. The cell worker is
  // only responsible for initializing the copy object for the current cell,
  // whereas the actual local work happens on boundary faces: there, the worker
  // evaluates the inhomogeneous boundary values and accumulates the associated
  // Nitsche terms into the local right-hand-side vector. The copier then adds
  // that local vector to the global right-hand side.
  template <int dim>
  void InteriorPenaltyProblem<dim>::assemble_right_hand_side()
  {
    using CellIterator = typename DoFHandler<dim>::active_cell_iterator;

    const RightHandSideScratchData<dim> scratch(mapping,
                                                fe,
                                                QGauss<dim - 1>(fe.degree + 1),
                                                update_quadrature_points |
                                                  update_values |
                                                  update_gradients |
                                                  update_JxW_values |
                                                  update_normal_vectors);
    const RightHandSideCopyData         copy_data(fe.n_dofs_per_cell());

    MeshWorker::mesh_loop(
      dof_handler.begin_active(),
      dof_handler.end(),
      /* cell worker: */
      [&](const CellIterator &cell,
          RightHandSideScratchData<dim> &,
          RightHandSideCopyData &copy) { copy.reinit(cell); },
      /* copier: */
      [&](const RightHandSideCopyData &copy) {
        right_hand_side.add(copy.local_dof_indices, copy.cell_rhs);
      },
      /* scratch and copy objects: */
      scratch,
      copy_data,
      /* where and how we want to integrate: */
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces,
      /* boundary face worker: */
      [&](const CellIterator            &cell,
          const unsigned int             face_no,
          RightHandSideScratchData<dim> &scratch_data,
          RightHandSideCopyData         &copy) {
        scratch_data.boundary_fe_values.reinit(cell, face_no);
        exact_solution.value_list(
          scratch_data.boundary_fe_values.get_quadrature_points(),
          scratch_data.boundary_values);
        RHSIntegrator::boundary<dim>(cell,
                                     face_no,
                                     scratch_data.boundary_fe_values,
                                     scratch_data.boundary_values,
                                     copy.cell_rhs);
      });

    right_hand_side *= -1.;
  }


  // Now that we have coded all functions building the discrete linear system,
  // it is about time that we actually solve it.
  template <int dim>
  void InteriorPenaltyProblem<dim>::solve()
  {
    // The solver of choice is conjugate gradient.
    SolverControl            control(1000, 1.e-12);
    SolverCG<Vector<double>> solver(control);

    // Now we are setting up the components of the multilevel
    // preconditioner. First, we need transfer between grid levels. The object
    // we are using here generates sparse matrices for these transfers.
    MGTransferPrebuilt<Vector<double>> mg_transfer;
    mg_transfer.build(dof_handler);

    // Then, we need an exact solver for the matrix on the coarsest level.
    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from(mg_matrix[0]);
    MGCoarseGridHouseholder<double, Vector<double>> mg_coarse;
    mg_coarse.initialize(coarse_matrix);

    // While transfer and coarse grid solver are pretty much generic, more
    // flexibility is offered for the smoother. First, we choose Gauss-Seidel
    // as our smoothing method.
    GrowingVectorMemory<Vector<double>> mem;
    using RELAXATION = PreconditionSOR<SparseMatrix<double>>;
    mg::SmootherRelaxation<RELAXATION, Vector<double>> mg_smoother;
    RELAXATION::AdditionalData                         smoother_data(1.);
    mg_smoother.initialize(mg_matrix, smoother_data);

    // Do two smoothing steps on each level.
    mg_smoother.set_steps(2);
    // Since the SOR method is not symmetric, but we use conjugate gradient
    // iteration below, here is a trick to make the multilevel preconditioner
    // a symmetric operator even for nonsymmetric smoothers.
    mg_smoother.set_symmetric(true);
    // The smoother class optionally implements the variable V-cycle, which we
    // do not want here.
    mg_smoother.set_variable(false);

    // Finally, we must wrap our matrices in an object having the required
    // multiplication functions.
    mg::Matrix<Vector<double>> mgmatrix(mg_matrix);
    mg::Matrix<Vector<double>> mgdown(mg_matrix_dg_down);
    mg::Matrix<Vector<double>> mgup(mg_matrix_dg_up);

    // Now, we are ready to set up the V-cycle operator and the multilevel
    // preconditioner.
    Multigrid<Vector<double>> mg(
      mgmatrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
    // Let us not forget the edge matrices needed because of the adaptive
    // refinement.
    mg.set_edge_flux_matrices(mgdown, mgup);

    // After all preparations, wrap the Multigrid object into another object,
    // which can be used as a regular preconditioner,
    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>>
      preconditioner(dof_handler, mg, mg_transfer);
    // and use it to solve the system.
    solver.solve(matrix, solution, right_hand_side, preconditioner);

    std::cout << "Converged in " << control.last_step() << " iterations"
              << std::endl;
  }


  // The next function estimates the error. The big difference to the previous
  // mesh loop functions is that we now also read from the discrete solution
  // vector. The results of the estimator are stored in a vector with one entry
  // per cell.
  //
  // As before, MeshWorker::mesh_loop() separates the local work from the
  // accumulation into the global output vector. The workers evaluate the
  // current solution on cells and faces, compute the cell, boundary, and jump
  // contributions to the estimator, and store them in the copy object. The
  // copier then distributes these local indicators to the per-cell entries of
  // the global estimator vector, splitting face terms evenly between the two
  // cells that share the face.
  template <int dim>
  double InteriorPenaltyProblem<dim>::estimate()
  {
    estimates.block(0).reinit(triangulation.n_active_cells());
    using CellIterator = typename DoFHandler<dim>::active_cell_iterator;

    const unsigned int n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;
    const EstimatorScratchData<dim> scratch(mapping,
                                            fe,
                                            QGauss<dim>(n_gauss_points),
                                            QGauss<dim - 1>(n_gauss_points + 1),
                                            QGauss<dim - 1>(n_gauss_points),
                                            update_hessians | update_JxW_values,
                                            update_quadrature_points |
                                              update_values | update_gradients |
                                              update_JxW_values |
                                              update_normal_vectors,
                                            update_values | update_gradients |
                                              update_JxW_values |
                                              update_normal_vectors);
    const ErrorCopyData<1>          copy_data;

    MeshWorker::mesh_loop(
      dof_handler.begin_active(),
      dof_handler.end(),
      /* cell worker: */
      [&](const CellIterator        &cell,
          EstimatorScratchData<dim> &scratch_data,
          ErrorCopyData<1>          &copy) {
        copy.reinit(cell);
        scratch_data.fe_values.reinit(cell);
        const FEValues<dim> &fe_values = scratch_data.fe_values;
        fe_values.get_function_hessians(solution, scratch_data.cell_hessians);
        copy.cell_values[0] =
          Estimator::cell<dim>(cell, fe_values, scratch_data.cell_hessians);
      },
      /* copier: */
      [&](const ErrorCopyData<1> &copy) {
        estimates.block(0)(copy.cell_index) += copy.cell_values[0];
        for (const auto &face : copy.face_data)
          {
            estimates.block(0)(face.cell_index_1) += 0.5 * face.values[0];
            estimates.block(0)(face.cell_index_2) += 0.5 * face.values[0];
          }
      },
      /* scratch and copy objects: */
      scratch,
      copy_data,
      /* where and how we want to integrate: */
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces |
        MeshWorker::assemble_own_interior_faces_once,
      /* boundary face worker: */
      [&](const CellIterator        &cell,
          const unsigned int         face_no,
          EstimatorScratchData<dim> &scratch_data,
          ErrorCopyData<1>          &copy) {
        scratch_data.boundary_fe_values.reinit(cell, face_no);
        const FEFaceValues<dim> &fe_face_values =
          scratch_data.boundary_fe_values;
        fe_face_values.get_function_values(
          solution, scratch_data.boundary_solution_values);
        exact_solution.value_list(fe_face_values.get_quadrature_points(),
                                  scratch_data.boundary_exact_values);
        copy.cell_values[0] +=
          Estimator::boundary<dim>(cell,
                                   face_no,
                                   fe_face_values,
                                   scratch_data.boundary_solution_values,
                                   scratch_data.boundary_exact_values);
      },
      /* interior face worker: */
      [&](const CellIterator        &cell,
          const unsigned int         face_no,
          const unsigned int         subface_no,
          const CellIterator        &neighbor,
          const unsigned int         neighbor_face_no,
          const unsigned int         neighbor_subface_no,
          EstimatorScratchData<dim> &scratch_data,
          ErrorCopyData<1>          &copy) {
        auto &face_data = copy.emplace_face_data(cell, neighbor);

        if (subface_no == numbers::invalid_unsigned_int)
          {
            scratch_data.face_fe_values.reinit(cell, face_no);
            const FEFaceValuesBase<dim> &fe_face_values =
              scratch_data.face_fe_values;

            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_face_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);
                fe_face_values.get_function_gradients(
                  solution, scratch_data.face_solution_gradients);
                neighbor_fe_face_values.get_function_gradients(
                  solution, scratch_data.neighbor_face_solution_gradients);

                face_data.values[0] = Estimator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  scratch_data.face_solution_gradients,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values,
                  scratch_data.neighbor_face_solution_gradients);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_subface_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);
                fe_face_values.get_function_gradients(
                  solution, scratch_data.face_solution_gradients);
                neighbor_fe_face_values.get_function_gradients(
                  solution, scratch_data.neighbor_face_solution_gradients);

                face_data.values[0] = Estimator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  scratch_data.face_solution_gradients,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values,
                  scratch_data.neighbor_face_solution_gradients);
              }
          }
        else
          {
            scratch_data.subface_values.reinit(cell, face_no, subface_no);
            const FEFaceValuesBase<dim> &fe_face_values =
              scratch_data.subface_values;

            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_face_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);
                fe_face_values.get_function_gradients(
                  solution, scratch_data.face_solution_gradients);
                neighbor_fe_face_values.get_function_gradients(
                  solution, scratch_data.neighbor_face_solution_gradients);

                face_data.values[0] = Estimator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  scratch_data.face_solution_gradients,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values,
                  scratch_data.neighbor_face_solution_gradients);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_subface_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);
                fe_face_values.get_function_gradients(
                  solution, scratch_data.face_solution_gradients);
                neighbor_fe_face_values.get_function_gradients(
                  solution, scratch_data.neighbor_face_solution_gradients);

                face_data.values[0] = Estimator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  scratch_data.face_solution_gradients,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values,
                  scratch_data.neighbor_face_solution_gradients);
              }
          }
      });

    return estimates.block(0).l2_norm();
  }

  // Here we compare our finite element solution with the known exact solution
  // and compute the mean quadratic error of the gradient and the function
  // itself. This function is a close relative of the estimation function right
  // above: the mesh loop again visits cells, boundary faces, and interior
  // faces; the workers evaluate local quantities with the help of the scratch
  // object; and the copier writes the resulting indicators into global data
  // structures only after the local computation is finished.
  //
  // Since we compute the error in the energy and the
  // <i>L<sup>2</sup></i>-norm, respectively, our block vector needs two
  // blocks here. Consequently, the copy object stores two local values per
  // cell and per face contribution, and the copier accumulates each of them
  // into the corresponding block.
  template <int dim>
  void InteriorPenaltyProblem<dim>::error()
  {
    BlockVector<double> errors(2);
    errors.block(0).reinit(triangulation.n_active_cells());
    errors.block(1).reinit(triangulation.n_active_cells());
    using CellIterator = typename DoFHandler<dim>::active_cell_iterator;

    const unsigned int n_gauss_points =
      dof_handler.get_fe().tensor_degree() + 1;
    const ErrorScratchData<dim> scratch(mapping,
                                        fe,
                                        QGauss<dim>(n_gauss_points),
                                        QGauss<dim - 1>(n_gauss_points + 1),
                                        QGauss<dim - 1>(n_gauss_points),
                                        update_quadrature_points |
                                          update_values | update_gradients |
                                          update_JxW_values,
                                        update_quadrature_points |
                                          update_values | update_JxW_values,
                                        update_values | update_JxW_values);
    const ErrorCopyData<2>      copy_data;

    MeshWorker::mesh_loop(
      dof_handler.begin_active(),
      dof_handler.end(),
      /* cell worker: */
      [&](const CellIterator    &cell,
          ErrorScratchData<dim> &scratch_data,
          ErrorCopyData<2>      &copy) {
        copy.reinit(cell);
        scratch_data.fe_values.reinit(cell);
        const FEValues<dim> &fe_values = scratch_data.fe_values;
        fe_values.get_function_values(solution,
                                      scratch_data.cell_solution_values);
        fe_values.get_function_gradients(solution,
                                         scratch_data.cell_solution_gradients);
        exact_solution.value_list(fe_values.get_quadrature_points(),
                                  scratch_data.cell_exact_values);
        exact_solution.gradient_list(fe_values.get_quadrature_points(),
                                     scratch_data.cell_exact_gradients);
        copy.cell_values =
          ErrorIntegrator::cell<dim>(fe_values,
                                     scratch_data.cell_solution_values,
                                     scratch_data.cell_solution_gradients,
                                     scratch_data.cell_exact_values,
                                     scratch_data.cell_exact_gradients);
      },
      /* copier: */
      [&](const ErrorCopyData<2> &copy) {
        errors.block(0)(copy.cell_index) += copy.cell_values[0];
        errors.block(1)(copy.cell_index) += copy.cell_values[1];
        for (const auto &face : copy.face_data)
          {
            errors.block(0)(face.cell_index_1) += 0.5 * face.values[0];
            errors.block(0)(face.cell_index_2) += 0.5 * face.values[0];
          }
      },
      /* scratch and copy objects: */
      scratch,
      copy_data,
      /* where and how we want to integrate: */
      MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces |
        MeshWorker::assemble_own_interior_faces_once,
      /* boundary face worker: */
      [&](const CellIterator    &cell,
          const unsigned int     face_no,
          ErrorScratchData<dim> &scratch_data,
          ErrorCopyData<2>      &copy) {
        scratch_data.boundary_fe_values.reinit(cell, face_no);
        const FEFaceValues<dim> &fe_face_values =
          scratch_data.boundary_fe_values;
        fe_face_values.get_function_values(
          solution, scratch_data.boundary_solution_values);
        exact_solution.value_list(fe_face_values.get_quadrature_points(),
                                  scratch_data.boundary_exact_values);
        copy.cell_values[0] +=
          ErrorIntegrator::boundary<dim>(cell,
                                         face_no,
                                         fe_face_values,
                                         scratch_data.boundary_solution_values,
                                         scratch_data.boundary_exact_values);
      },
      /* interior face worker: */
      [&](const CellIterator    &cell,
          const unsigned int     face_no,
          const unsigned int     subface_no,
          const CellIterator    &neighbor,
          const unsigned int     neighbor_face_no,
          const unsigned int     neighbor_subface_no,
          ErrorScratchData<dim> &scratch_data,
          ErrorCopyData<2>      &copy) {
        auto &face_data = copy.emplace_face_data(cell, neighbor);

        if (subface_no == numbers::invalid_unsigned_int)
          {
            scratch_data.face_fe_values.reinit(cell, face_no);
            const FEFaceValuesBase<dim> &fe_face_values =
              scratch_data.face_fe_values;

            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_face_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);

                face_data.values[0] = ErrorIntegrator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_subface_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);

                face_data.values[0] = ErrorIntegrator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values);
              }
          }
        else
          {
            scratch_data.subface_values.reinit(cell, face_no, subface_no);
            const FEFaceValuesBase<dim> &fe_face_values =
              scratch_data.subface_values;

            if (neighbor_subface_no == numbers::invalid_unsigned_int)
              {
                scratch_data.neighbor_face_values.reinit(neighbor,
                                                         neighbor_face_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_face_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);

                face_data.values[0] = ErrorIntegrator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values);
              }
            else
              {
                scratch_data.neighbor_subface_values.reinit(
                  neighbor, neighbor_face_no, neighbor_subface_no);
                const FEFaceValuesBase<dim> &neighbor_fe_face_values =
                  scratch_data.neighbor_subface_values;

                fe_face_values.get_function_values(
                  solution, scratch_data.face_solution_values);
                neighbor_fe_face_values.get_function_values(
                  solution, scratch_data.neighbor_face_solution_values);

                face_data.values[0] = ErrorIntegrator::face<dim>(
                  cell,
                  face_no,
                  fe_face_values,
                  scratch_data.face_solution_values,
                  neighbor,
                  neighbor_face_no,
                  scratch_data.neighbor_face_solution_values);
              }
          }
      });

    std::cout << "energy-error: " << errors.block(0).l2_norm() << std::endl;
    std::cout << "L2-error:     " << errors.block(1).l2_norm() << std::endl;
  }


  // Create graphical output. We produce the filename by collating the
  // name from its various components, including the refinement cycle
  // that we output with two digits.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::output_results(const unsigned int cycle) const
  {
    const std::string filename =
      "sol-" + Utilities::int_to_string(cycle, 2) + ".gnuplot";

    std::cout << "Writing solution to <" << filename << ">..." << std::endl
              << std::endl;
    std::ofstream gnuplot_output(filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "u");
    data_out.add_data_vector(estimates.block(0), "est");

    data_out.build_patches();

    data_out.write_gnuplot(gnuplot_output);
  }

  // And finally the adaptive loop, more or less like in previous examples.
  template <int dim>
  void InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    std::cout << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s = 0; s < n_steps; ++s)
      {
        std::cout << "Step " << s << std::endl;
        if (estimates.block(0).empty())
          triangulation.refine_global(1);
        else
          {
            GridRefinement::refine_and_coarsen_fixed_fraction(
              triangulation, estimates.block(0), 0.5, 0.0);
            triangulation.execute_coarsening_and_refinement();
          }

        std::cout << "Triangulation " << triangulation.n_active_cells()
                  << " cells, " << triangulation.n_levels() << " levels"
                  << std::endl;

        setup_system();
        std::cout << "DoFHandler " << dof_handler.n_dofs()
                  << " dofs, level dofs";
        for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
          std::cout << ' ' << dof_handler.n_dofs(l);
        std::cout << std::endl;

        std::cout << "Assemble matrix" << std::endl;
        assemble_matrix();
        std::cout << "Assemble multilevel matrix" << std::endl;
        assemble_mg_matrix();
        std::cout << "Assemble right hand side" << std::endl;
        assemble_right_hand_side();
        std::cout << "Solve" << std::endl;
        solve();
        error();
        std::cout << "Estimate " << estimate() << std::endl;
        output_results(s);
      }
  }
} // namespace Step39



int main()
{
  try
    {
      using namespace Step39;

      InteriorPenaltyProblem<2> test1;
      test1.run(12);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
