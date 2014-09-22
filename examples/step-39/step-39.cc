/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Guido Kanschat, Texas A&M University, 2009
 */


// The include files for the linear algebra: A regular SparseMatrix, which in
// turn will include the necessary files for SparsityPattern and Vector
// classes.
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
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

// The include files for using the MeshWorker framework
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

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

#include <iostream>
#include <fstream>

// All classes of the deal.II library are in the namespace dealii. In order to
// save typing, we tell the compiler to search names in there as well.
namespace Step39
{
  using namespace dealii;

  // This is the function we use to set the boundary values and also the exact
  // solution we compare to.
  Functions::SlitSingularityFunction<2> exact_solution;

  // @sect3{The local integrators}

  // MeshWorker separates local integration from the loops over cells and
  // faces. Thus, we have to write local integration classes for generating
  // matrices, the right hand side and the error estimator.

  // All these classes have the same three functions for integrating over
  // cells, boundary faces and interior faces, respectively. All the
  // information needed for the local integration is provided by
  // MeshWorker::IntegrationInfo<dim>. Note that the signature of the
  // functions cannot be changed, because it is expected by
  // MeshWorker::integration_loop().

  // The first class defining local integrators is responsible for computing
  // cell and face matrices. It is used to assemble the global matrix as well
  // as the level matrices.
  template <int dim>
  class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo,
              typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                  typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  // On each cell, we integrate the Dirichlet form. We use the library of
  // ready made integrals in LocalIntegrators to avoid writing these loops
  // ourselves. Similarly, we implement Nitsche boundary conditions and the
  // interior penalty fluxes between cells.
  //
  // The boundary and flux terms need a penalty parameter, which should be
  // adjusted to the cell size and the polynomial degree. A safe choice of
  // this parameter for constant coefficients can be found in
  // LocalIntegrators::Laplace::compute_penalty() and we use this below.
  template <int dim>
  void MatrixIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values());
  }


  template <int dim>
  void MatrixIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::nitsche_matrix(
      dinfo.matrix(0,false).matrix, info.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, deg, deg));
  }

  // Interior faces use the interior penalty method
  template <int dim>
  void MatrixIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &dinfo1,
    MeshWorker::DoFInfo<dim> &dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
    LocalIntegrators::Laplace::ip_matrix(
      dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
      dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
      info1.fe_values(0), info2.fe_values(0),
      LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
  }

  // The second local integrator builds the right hand side. In our example,
  // the right hand side function is zero, such that only the boundary
  // condition is set here in weak form.
  template <int dim>
  class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  template <int dim>
  void RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &, typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  template <int dim>
  void RHSIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();
    Vector<double> &local_vector = dinfo.vector(0).block(0);

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        local_vector(i) += (- fe.shape_value(i,k) * penalty * boundary_values[k]
                            + (fe.normal_vector(k) * fe.shape_grad(i,k)) * boundary_values[k])
                           * fe.JxW(k);
  }


  template <int dim>
  void RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
                                MeshWorker::DoFInfo<dim> &,
                                typename MeshWorker::IntegrationInfo<dim> &,
                                typename MeshWorker::IntegrationInfo<dim> &) const
  {}


  // The third local integrator is responsible for the contributions to the
  // error estimate. This is the standard energy estimator due to Karakashian
  // and Pascal (2003).
  template <int dim>
  class Estimator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };


  // The cell contribution is the Laplacian of the discrete solution, since
  // the right hand side is zero.
  template <int dim>
  void Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    const std::vector<Tensor<2,dim> > &DDuh = info.hessians[0][0];
    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        const double t = dinfo.cell->diameter() * trace(DDuh[k]);
        dinfo.value(0) +=  t*t * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }

  // At the boundary, we use simply a weighted form of the boundary residual,
  // namely the norm of the difference between the finite element solution and
  // the correct boundary condition.
  template <int dim>
  void Estimator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> boundary_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      dinfo.value(0) += penalty * (boundary_values[k] - uh[k]) * (boundary_values[k] - uh[k])
                        * fe.JxW(k);
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  // Finally, on interior faces, the estimator consists of the jumps of the
  // solution and its normal derivative, weighted appropriately.
  template <int dim>
  void Estimator<dim>::face(MeshWorker::DoFInfo<dim> &dinfo1,
                            MeshWorker::DoFInfo<dim> &dinfo2,
                            typename MeshWorker::IntegrationInfo<dim> &info1,
                            typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &fe = info1.fe_values();
    const std::vector<double> &uh1 = info1.values[0][0];
    const std::vector<double> &uh2 = info2.values[0][0];
    const std::vector<Tensor<1,dim> > &Duh1 = info1.gradients[0][0];
    const std::vector<Tensor<1,dim> > &Duh2 = info2.gradients[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty1 = deg * (deg+1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 = deg * (deg+1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;
    const double h = dinfo1.face->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double diff1 = uh1[k] - uh2[k];
        double diff2 = fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k];
        dinfo1.value(0) += (penalty * diff1*diff1 + h * diff2*diff2)
                           * fe.JxW(k);
      }
    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
    dinfo2.value(0) = dinfo1.value(0);
  }

  // Finally we have an integrator for the error. Since the energy norm for
  // discontinuous Galerkin problems not only involves the difference of the
  // gradient inside the cells, but also the jump terms across faces and at
  // the boundary, we cannot just use VectorTools::integrate_difference().
  // Instead, we use the MeshWorker interface to compute the error ourselves.

  // There are several different ways to define this energy norm, but all of
  // them are equivalent to each other uniformly with mesh size (some not
  // uniformly with polynomial degree). Here, we choose @f[ \|u\|_{1,h} =
  // \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
  // 4\sigma_F\|\{\!\{ u \mathbf n\}\!\}\|^2_F + \sum_{F \in F_h^b}
  // 2\sigma_F\|u\|^2_F @f]

  template <int dim>
  class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim>
  {
  public:
    void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
    void face(MeshWorker::DoFInfo<dim> &dinfo1,
              MeshWorker::DoFInfo<dim> &dinfo2,
              typename MeshWorker::IntegrationInfo<dim> &info1,
              typename MeshWorker::IntegrationInfo<dim> &info2) const;
  };

  // Here we have the integration on cells. There is currently no good
  // interface in MeshWorker that would allow us to access values of regular
  // functions in the quadrature points. Thus, we have to create the vectors
  // for the exact function's values and gradients inside the cell
  // integrator. After that, everything is as before and we just add up the
  // squares of the differences.

  // Additionally to computing the error in the energy norm, we use the
  // capability of the mesh worker to compute two functionals at the same time
  // and compute the <i>L<sup>2</sup></i>-error in the same loop. Obviously,
  // this one does not have any jump terms and only appears in the integration
  // on cells.
  template <int dim>
  void ErrorIntegrator<dim>::cell(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();
    std::vector<Tensor<1,dim> > exact_gradients(fe.n_quadrature_points);
    std::vector<double> exact_values(fe.n_quadrature_points);

    exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<Tensor<1,dim> > &Duh = info.gradients[0][0];
    const std::vector<double> &uh = info.values[0][0];

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double sum = 0;
        for (unsigned int d=0; d<dim; ++d)
          {
            const double diff = exact_gradients[k][d] - Duh[k][d];
            sum += diff*diff;
          }
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) +=  sum * fe.JxW(k);
        dinfo.value(1) +=  diff*diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
    dinfo.value(1) = std::sqrt(dinfo.value(1));
  }


  template <int dim>
  void ErrorIntegrator<dim>::boundary(
    MeshWorker::DoFInfo<dim> &dinfo,
    typename MeshWorker::IntegrationInfo<dim> &info) const
  {
    const FEValuesBase<dim> &fe = info.fe_values();

    std::vector<double> exact_values(fe.n_quadrature_points);
    exact_solution.value_list(fe.get_quadrature_points(), exact_values);

    const std::vector<double> &uh = info.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        const double diff = exact_values[k] - uh[k];
        dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
      }
    dinfo.value(0) = std::sqrt(dinfo.value(0));
  }


  template <int dim>
  void ErrorIntegrator<dim>::face(
    MeshWorker::DoFInfo<dim> &dinfo1,
    MeshWorker::DoFInfo<dim> &dinfo2,
    typename MeshWorker::IntegrationInfo<dim> &info1,
    typename MeshWorker::IntegrationInfo<dim> &info2) const
  {
    const FEValuesBase<dim> &fe = info1.fe_values();
    const std::vector<double> &uh1 = info1.values[0][0];
    const std::vector<double> &uh2 = info2.values[0][0];

    const unsigned int deg = fe.get_fe().tensor_degree();
    const double penalty1 = deg * (deg+1) * dinfo1.face->measure() / dinfo1.cell->measure();
    const double penalty2 = deg * (deg+1) * dinfo2.face->measure() / dinfo2.cell->measure();
    const double penalty = penalty1 + penalty2;

    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
      {
        double diff = uh1[k] - uh2[k];
        dinfo1.value(0) += (penalty * diff*diff)
                           * fe.JxW(k);
      }
    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
    dinfo2.value(0) = dinfo1.value(0);
  }



  // @sect3{The main class}

  // This class does the main job, like in previous examples. For a
  // description of the functions declared here, please refer to the
  // implementation below.
  template <int dim>
  class InteriorPenaltyProblem
  {
  public:
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    InteriorPenaltyProblem(const FiniteElement<dim> &fe);

    void run(unsigned int n_steps);

  private:
    void setup_system ();
    void assemble_matrix ();
    void assemble_mg_matrix ();
    void assemble_right_hand_side ();
    void error ();
    double estimate ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    // The member objects related to the discretization are here.
    Triangulation<dim>        triangulation;
    const MappingQ1<dim>      mapping;
    const FiniteElement<dim> &fe;
    DoFHandler<dim>           dof_handler;

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
    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparseMatrix<double> > mg_matrix;

    // When we perform multigrid with local smoothing on locally refined
    // meshes, additional matrices are required; see Kanschat (2004). Here is
    // the sparsity pattern for these edge matrices. We only need one, because
    // the pattern of the up matrix is the transpose of that of the down
    // matrix. Actually, we do not care too much about these details, since
    // the MeshWorker is filling these matrices.
    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
    // The flux matrix at the refinement edge, coupling fine level degrees of
    // freedom to coarse level.
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_down;
    // The transpose of the flux matrix at the refinement edge, coupling
    // coarse level degrees of freedom to fine level.
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_up;
  };


  // The constructor simply sets up the coarse grid and the DoFHandler. The
  // FiniteElement is provided as a parameter to allow flexibility.
  template <int dim>
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem(const FiniteElement<dim> &fe)
    :
    mapping(),
    fe(fe),
    dof_handler(triangulation),
    estimates(1)
  {
    GridGenerator::hyper_cube_slit(triangulation, -1, 1);
  }


  // In this function, we set up the dimension of the linear system and the
  // sparsity patterns for the global matrix as well as the level matrices.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::setup_system()
  {
    // First, we use the finite element to distribute degrees of freedom over
    // the mesh and number them.
    dof_handler.distribute_dofs(fe);
    dof_handler.distribute_mg_dofs(fe);
    unsigned int n_dofs = dof_handler.n_dofs();
    // Then, we already know the size of the vectors representing finite
    // element functions.
    solution.reinit(n_dofs);
    right_hand_side.reinit(n_dofs);

    // Next, we set up the sparsity pattern for the global matrix. Since we do
    // not know the row sizes in advance, we first fill a temporary
    // CompressedSparsityPattern object and copy it to the regular
    // SparsityPattern once it is complete.
    CompressedSparsityPattern c_sparsity(n_dofs);
    DoFTools::make_flux_sparsity_pattern(dof_handler, c_sparsity);
    sparsity.copy_from(c_sparsity);
    matrix.reinit(sparsity);

    const unsigned int n_levels = triangulation.n_levels();
    // The global system is set up, now we attend to the level matrices. We
    // resize all matrix objects to hold one matrix per level.
    mg_matrix.resize(0, n_levels-1);
    mg_matrix.clear();
    mg_matrix_dg_up.resize(0, n_levels-1);
    mg_matrix_dg_up.clear();
    mg_matrix_dg_down.resize(0, n_levels-1);
    mg_matrix_dg_down.clear();
    // It is important to update the sparsity patterns after <tt>clear()</tt>
    // was called for the level matrices, since the matrices lock the sparsity
    // pattern through the Smartpointer and Subscriptor mechanism.
    mg_sparsity.resize(0, n_levels-1);
    mg_sparsity_dg_interface.resize(0, n_levels-1);

    // Now all objects are prepared to hold one sparsity pattern or matrix per
    // level. What's left is setting up the sparsity patterns on each level.
    for (unsigned int level=mg_sparsity.min_level();
         level<=mg_sparsity.max_level(); ++level)
      {
        // These are roughly the same lines as above for the global matrix,
        // now for each level.
        CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
        mg_sparsity[level].copy_from(c_sparsity);
        mg_matrix[level].reinit(mg_sparsity[level]);

        // Additionally, we need to initialize the transfer matrices at the
        // refinement edge between levels. They are stored at the index
        // referring to the finer of the two indices, thus there is no such
        // object on level 0.
        if (level>0)
          {
            CompressedSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level-1), dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, ci_sparsity, level);
            mg_sparsity_dg_interface[level].copy_from(ci_sparsity);
            mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
            mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
          }
      }
  }


  // In this function, we assemble the global system matrix, where by global
  // we indicate that this is the matrix of the discrete system we solve and
  // it is covering the whole mesh.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_matrix()
  {
    // First, we need t set up the object providing the values we
    // integrate. This object contains all FEValues and FEFaceValues objects
    // needed and also maintains them automatically such that they always
    // point to the current cell. To this end, we need to tell it first, where
    // and what to compute. Since we are not doing anything fancy, we can rely
    // on their standard choice for quadrature rules.
    //
    // Since their default update flags are minimal, we add what we need
    // additionally, namely the values and gradients of shape functions on all
    // objects (cells, boundary and interior faces). Afterwards, we are ready
    // to initialize the container, which will create all necessary
    // FEValuesBase objects for integration.
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    // This is the object into which we integrate local data. It is filled by
    // the local integration routines in MatrixIntegrator and then used by the
    // assembler to distribute the information into the global matrix.
    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Furthermore, we need an object that assembles the local matrix into the
    // global matrix. These assembler objects have all the knowledge
    // of the structures of the target object, in this case a
    // SparseMatrix, possible constraints and the mesh structure.
    MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(matrix);

    // Now comes the part we coded ourselves, the local
    // integrator. This is the only part which is problem dependent.
    MatrixIntegrator<dim> integrator;
    // Now, we throw everything into a MeshWorker::loop(), which here
    // traverses all active cells of the mesh, computes cell and face matrices
    // and assembles them into the global matrix. We use the variable
    // <tt>dof_handler</tt> here in order to use the global numbering of
    // degrees of freedom.
    MeshWorker::integration_loop<dim, dim>(
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);
  }


  // Now, we do the same for the level matrices. Not too surprisingly, this
  // function looks like a twin of the previous one. Indeed, there are only
  // two minor differences.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_mg_matrix()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Obviously, the assembler needs to be replaced by one filling level
    // matrices. Note that it automatically fills the edge matrices as well.
    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(mg_matrix);
    assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down);

    MatrixIntegrator<dim> integrator;
    // Here is the other difference to the previous function: we run
    // over all cells, not only the active ones. And we use functions
    // ending on <code>_mg</code> since we need the degrees of freedom
    // on each level, not the global numbering.
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_mg(), dof_handler.end_mg(),
      dof_info, info_box,
      integrator, assembler);
  }


  // Here we have another clone of the assemble function. The difference to
  // assembling the system matrix consists in that we assemble a vector here.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::assemble_right_hand_side()
  {
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // Since this assembler allows us to fill several vectors, the interface is
    // a little more complicated as above. The pointers to the vectors have to
    // be stored in a NamedData object. While this seems to cause two extra
    // lines of code here, it actually comes handy in more complex
    // applications.
    MeshWorker::Assembler::ResidualSimple<Vector<double> > assembler;
    NamedData<Vector<double>* > data;
    Vector<double> *rhs = &right_hand_side;
    data.add(rhs, "RHS");
    assembler.initialize(data);

    RHSIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim>(
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    right_hand_side *= -1.;
  }


  // Now that we have coded all functions building the discrete linear system,
  // it is about time that we actually solve it.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::solve()
  {
    // The solver of choice is conjugate gradient.
    SolverControl control(1000, 1.e-12);
    SolverCG<Vector<double> > solver(control);

    // Now we are setting up the components of the multilevel
    // preconditioner. First, we need transfer between grid levels. The object
    // we are using here generates sparse matrices for these transfers.
    MGTransferPrebuilt<Vector<double> > mg_transfer;
    mg_transfer.build_matrices(dof_handler);

    // Then, we need an exact solver for the matrix on the coarsest level.
    FullMatrix<double> coarse_matrix;
    coarse_matrix.copy_from (mg_matrix[0]);
    MGCoarseGridHouseholder<double, Vector<double> > mg_coarse;
    mg_coarse.initialize(coarse_matrix);

    // While transfer and coarse grid solver are pretty much generic, more
    // flexibility is offered for the smoother. First, we choose Gauss-Seidel
    // as our smoothing method.
    GrowingVectorMemory<Vector<double> > mem;
    typedef PreconditionSOR<SparseMatrix<double> > RELAXATION;
    mg::SmootherRelaxation<RELAXATION, Vector<double> >
    mg_smoother;
    RELAXATION::AdditionalData smoother_data(1.);
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
    MGMatrix<SparseMatrix<double>, Vector<double> > mgmatrix(&mg_matrix);
    MGMatrix<SparseMatrix<double>, Vector<double> > mgdown(&mg_matrix_dg_down);
    MGMatrix<SparseMatrix<double>, Vector<double> > mgup(&mg_matrix_dg_up);

    // Now, we are ready to set up the V-cycle operator and the multilevel
    // preconditioner.
    Multigrid<Vector<double> > mg(dof_handler, mgmatrix,
                                  mg_coarse, mg_transfer,
                                  mg_smoother, mg_smoother);
    // Let us not forget the edge matrices needed because of the adaptive
    // refinement.
    mg.set_edge_flux_matrices(mgdown, mgup);

    // After all preparations, wrap the Multigrid object into another object,
    // which can be used as a regular preconditioner,
    PreconditionMG<dim, Vector<double>,
                   MGTransferPrebuilt<Vector<double> > >
                   preconditioner(dof_handler, mg, mg_transfer);
    // and use it to solve the system.
    solver.solve(matrix, solution, right_hand_side, preconditioner);
  }


  // Another clone of the assemble function. The big difference to the
  // previous ones is here that we also have an input vector.
  template <int dim>
  double
  InteriorPenaltyProblem<dim>::estimate()
  {
    // The results of the estimator are stored in a vector with one entry per
    // cell. Since cells in deal.II are not numbered, we have to create our
    // own numbering in order to use this vector.
    //
    // On the other hand, somebody might have used the user indices
    // already. So, let's be good citizens and save them before tampering with
    // them.
    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);

    estimates.block(0).reinit(triangulation.n_active_cells());
    unsigned int i=0;
    for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell,++i)
      cell->set_user_index(i);

    // This starts like before,
    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
    info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);

    // but now we need to notify the info box of the finite element function we
    // want to evaluate in the quadrature points. First, we create a NamedData
    // object with this vector, which is the solution we just computed.
    NamedData<Vector<double>* > solution_data;
    solution_data.add(&solution, "solution");

    // Then, we tell the Meshworker::VectorSelector for cells, that we need
    // the second derivatives of this solution (to compute the
    // Laplacian). Therefore, the Boolean arguments selecting function values
    // and first derivatives a false, only the last one selecting second
    // derivatives is true.
    info_box.cell_selector.add("solution", false, false, true);
    // On interior and boundary faces, we need the function values and the
    // first derivatives, but not second derivatives.
    info_box.boundary_selector.add("solution", true, true, false);
    info_box.face_selector.add("solution", true, true, false);

    // And we continue as before, with the exception that the default update
    // flags are already adjusted to the values and derivatives we requested
    // above.
    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    // The assembler stores one number per cell, but else this is the same as
    // in the computation of the right hand side.
    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    NamedData<BlockVector<double>* > out_data;
    BlockVector<double> *est = &estimates;
    out_data.add(est, "cells");
    assembler.initialize(out_data, false);

    Estimator<dim> integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    // Right before we return the result of the error estimate, we restore the
    // old user indices.
    triangulation.load_user_indices(old_user_indices);
    return estimates.block(0).l2_norm();
  }

  // Here we compare our finite element solution with the (known) exact
  // solution and compute the mean quadratic error of the gradient and the
  // function itself. This function is a clone of the estimation function
  // right above.

  // Since we compute the error in the energy and the
  // <i>L<sup>2</sup></i>-norm, respectively, our block vector needs two
  // blocks here.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::error()
  {
    BlockVector<double> errors(2);
    errors.block(0).reinit(triangulation.n_active_cells());
    errors.block(1).reinit(triangulation.n_active_cells());
    unsigned int i=0;
    for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell,++i)
      cell->set_user_index(i);

    MeshWorker::IntegrationInfoBox<dim> info_box;
    const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
    info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);

    NamedData<Vector<double>* > solution_data;
    solution_data.add(&solution, "solution");

    info_box.cell_selector.add("solution", true, true, false);
    info_box.boundary_selector.add("solution", true, false, false);
    info_box.face_selector.add("solution", true, false, false);

    info_box.add_update_flags_cell(update_quadrature_points);
    info_box.add_update_flags_boundary(update_quadrature_points);
    info_box.initialize(fe, mapping, solution_data);

    MeshWorker::DoFInfo<dim> dof_info(dof_handler);

    MeshWorker::Assembler::CellsAndFaces<double> assembler;
    NamedData<BlockVector<double>* > out_data;
    BlockVector<double> *est = &errors;
    out_data.add(est, "cells");
    assembler.initialize(out_data, false);

    ErrorIntegrator<dim> integrator;
    MeshWorker::integration_loop<dim, dim> (
      dof_handler.begin_active(), dof_handler.end(),
      dof_info, info_box,
      integrator, assembler);

    deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl;
    deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl;
  }


  // Some graphical output
  template <int dim>
  void InteriorPenaltyProblem<dim>::output_results (const unsigned int cycle) const
  {
    // Output of the solution in gnuplot format.
    char *fn = new char[100];
    sprintf(fn, "sol-%02d", cycle);

    std::string filename(fn);
    filename += ".gnuplot";
    deallog << "Writing solution to <" << filename << ">..."
            << std::endl << std::endl;
    std::ofstream gnuplot_output (filename.c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "u");
    data_out.add_data_vector (estimates.block(0), "est");

    data_out.build_patches ();

    data_out.write_gnuplot(gnuplot_output);
  }

  // And finally the adaptive loop, more or less like in previous examples.
  template <int dim>
  void
  InteriorPenaltyProblem<dim>::run(unsigned int n_steps)
  {
    deallog << "Element: " << fe.get_name() << std::endl;
    for (unsigned int s=0; s<n_steps; ++s)
      {
        deallog << "Step " << s << std::endl;
        if (estimates.block(0).size() == 0)
          triangulation.refine_global(1);
        else
          {
            GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
                                                               estimates.block(0),
                                                               0.5, 0.0);
            triangulation.execute_coarsening_and_refinement ();
          }

        deallog << "Triangulation "
                << triangulation.n_active_cells() << " cells, "
                << triangulation.n_levels() << " levels" << std::endl;

        setup_system();
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs";
        for (unsigned int l=0; l<triangulation.n_levels(); ++l)
          deallog << ' ' << dof_handler.n_dofs(l);
        deallog << std::endl;

        deallog << "Assemble matrix" << std::endl;
        assemble_matrix();
        deallog << "Assemble multilevel matrix" << std::endl;
        assemble_mg_matrix();
        deallog << "Assemble right hand side" << std::endl;
        assemble_right_hand_side();
        deallog << "Solve" << std::endl;
        solve();
        error();
        deallog << "Estimate " << estimate() << std::endl;
        output_results(s);
      }
  }
}



int main()
{
  try
    {
      using namespace dealii;
      using namespace Step39;

      std::ofstream logfile("deallog");
      deallog.attach(logfile);
      FE_DGQ<2> fe1(3);
      InteriorPenaltyProblem<2> test1(fe1);
      test1.run(12);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
