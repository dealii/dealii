/* $Id$ */
/* Author: Guido Kanschat, Texas A&M University, 2009 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // The include files for the linear
				 // algebra: A regular SparseMatrix,
				 // which in turn will include the
				 // necessary files for
				 // SparsityPattern and Vector classes.
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/precondition_block.h>

				 // Include files for setting up the
				 // mesh
#include <grid/grid_generator.h>

				 // Include files for FiniteElement
				 // classes and DoFHandler.
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_handler.h>

				 // The include files for using the
				 // MeshWorker framework
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_loop.h>

				 // Support for multigrid methods
#include <multigrid/mg_tools.h>
#include <multigrid/multigrid.h>
#include <multigrid/mg_matrix.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>

				 // Finally, we take our exact
				 // solution from the library as well
				 // as quadrature and additional tools.
#include <base/function_lib.h>
#include <base/quadrature_lib.h>
#include <numerics/vectors.h>
#include <numerics/data_out.h>

#include <iostream>
#include <fstream>

				 // All classes of the deal.II library
				 // are in the namespace dealii. In
				 // order to save typing, we tell the
				 // compiler to search names in there
				 // as well.
using namespace dealii;

				 // This is the function we use to set
				 // the boundary values and also the
				 // exact solution we compare to.
//Functions::LSingularityFunction exact_solution;
Functions::Q1WedgeFunction<2> exact_solution;

				 // @sect3{The local integrators}

				 // MeshWorker separates local
				 // integration from the loops over
				 // cells and faces. Thus, we have to
				 // write local integration classes
				 // for generating matrices, the right
				 // hand side and the error
				 // estimator.

				 // All these classes have the same
				 // three functions for integrating
				 // over cells, boundary faces and
				 // interior faces, respectively. All
				 // the information needed for the
				 // local integration is provided by
				 // MeshWorker::IntegrationWorker<dim>::CellInfo
				 // and
				 // MeshWorker::IntegrationWorker<dim>::FaceInfo. Note
				 // that this public interface cannot
				 // be changed, because it is expected
				 // by MeshWorker::integration_loop().
template <int dim>
class MatrixIntegrator : public Subscriptor
{
  public:
    void cell(typename MeshWorker::IntegrationWorker<dim>::CellInfo& info) const;
    void bdry(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info) const;
    void face(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info1,
	      typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info2) const;
};


template <int dim>
void MatrixIntegrator<dim>::cell(typename MeshWorker::IntegrationWorker<dim>::CellInfo& info) const
{
  const FEValuesBase<dim>& fe = info.fe();
  FullMatrix<double>& local_matrix = info.M1[0].matrix;
  
  for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	local_matrix(i,j) += (fe.shape_grad(i,k) * fe.shape_grad(j,k))
			     * fe.JxW(k);
}


template <int dim>
void MatrixIntegrator<dim>::bdry(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info) const
{
  const FEFaceValuesBase<dim>& fe = info.fe();
  FullMatrix<double>& local_matrix = info.M1[0].matrix;
  
  const unsigned int deg = fe.get_fe().tensor_degree();
  const double penalty = 2. * deg * (deg+1) * info.face->measure() / info.cell->measure();
  
  for (unsigned k=0;k<fe.n_quadrature_points;++k)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	local_matrix(i,j) += (fe.shape_value(i,k) * penalty * fe.shape_value(j,k)
			      - (fe.normal_vector(k) * fe.shape_grad(i,k)) * fe.shape_value(j,k)
			      - (fe.normal_vector(k) * fe.shape_grad(j,k)) * fe.shape_value(i,k))
			     * fe.JxW(k);
}


template <int dim>
void MatrixIntegrator<dim>::face(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info1,
				 typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info2) const
{
  const FEFaceValuesBase<dim>& fe1 = info1.fe();
  const FEFaceValuesBase<dim>& fe2 = info2.fe();
  FullMatrix<double>& matrix_v1u1 = info1.M1[0].matrix;
  FullMatrix<double>& matrix_v1u2 = info1.M2[0].matrix;
  FullMatrix<double>& matrix_v2u1 = info2.M2[0].matrix;
  FullMatrix<double>& matrix_v2u2 = info2.M1[0].matrix;
  
  const unsigned int deg = fe1.get_fe().tensor_degree();
  const double penalty1 = deg * (deg+1) * info1.face->measure() / info1.cell->measure();
  const double penalty2 = deg * (deg+1) * info2.face->measure() / info2.cell->measure();
  const double penalty = penalty1 + penalty2;
  
  for (unsigned k=0;k<fe1.n_quadrature_points;++k)
    for (unsigned int i=0; i<fe1.dofs_per_cell; ++i)
      for (unsigned int j=0; j<fe1.dofs_per_cell; ++j)
	{
	  matrix_v1u1(i,j) += (fe1.shape_value(i,k) * penalty * fe1.shape_value(j,k)
			       - (fe1.normal_vector(k) * fe1.shape_grad(i,k)) * fe1.shape_value(j,k)
			       - (fe1.normal_vector(k) * fe1.shape_grad(j,k)) * fe1.shape_value(i,k)
	  ) * .5 * fe1.JxW(k);
	  matrix_v1u2(i,j) += (-fe1.shape_value(i,k) * penalty * fe2.shape_value(j,k)
			       + (fe1.normal_vector(k) * fe1.shape_grad(i,k)) * fe2.shape_value(j,k)
			       - (fe1.normal_vector(k) * fe2.shape_grad(j,k)) * fe1.shape_value(i,k)
	  ) * .5 * fe1.JxW(k);
	  matrix_v2u1(i,j) += (-fe2.shape_value(i,k) * penalty * fe1.shape_value(j,k)
			       - (fe1.normal_vector(k) * fe2.shape_grad(i,k)) * fe1.shape_value(j,k)
			       + (fe1.normal_vector(k) * fe1.shape_grad(j,k)) * fe2.shape_value(i,k)
	  ) * .5 * fe1.JxW(k);
	  matrix_v2u2(i,j) += (fe2.shape_value(i,k) * penalty * fe2.shape_value(j,k)
			       + (fe1.normal_vector(k) * fe2.shape_grad(i,k)) * fe2.shape_value(j,k)
			       + (fe1.normal_vector(k) * fe2.shape_grad(j,k)) * fe2.shape_value(i,k)
	  ) * .5 * fe1.JxW(k);
	}
}


template <int dim>
class RHSIntegrator : public Subscriptor
{
  public:
    void cell(typename MeshWorker::IntegrationWorker<dim>::CellInfo& info) const;
    void bdry(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info) const;
    void face(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info1,
	      typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info2) const;
};


template <int dim>
void RHSIntegrator<dim>::cell(typename MeshWorker::IntegrationWorker<dim>::CellInfo&) const
{}


template <int dim>
void RHSIntegrator<dim>::bdry(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info) const
{
  const FEFaceValuesBase<dim>& fe = info.fe();
  Vector<double>& local_vector = info.R[0].block(0);
  
  std::vector<double> boundary_values(fe.n_quadrature_points);
  exact_solution.value_list(fe.get_quadrature_points(), boundary_values);
  
  const unsigned int deg = fe.get_fe().tensor_degree();
  const double penalty = 2. * deg * (deg+1) * info.face->measure() / info.cell->measure();
  
  for (unsigned k=0;k<fe.n_quadrature_points;++k)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      local_vector(i) += (- fe.shape_value(i,k) * penalty * boundary_values[k]
			  + (fe.normal_vector(k) * fe.shape_grad(i,k)) * boundary_values[k])
			 * fe.JxW(k);
}


template <int dim>
void RHSIntegrator<dim>::face(typename MeshWorker::IntegrationWorker<dim>::FaceInfo&,
			      typename MeshWorker::IntegrationWorker<dim>::FaceInfo&) const
{}


				 // @sect3{The main class}
template <int dim>
class Step39
{
  public:
    Step39(const FiniteElement<dim>& fe);

    void run(unsigned int n_steps);
    
  private:
    void setup_system ();
    void assemble_matrix ();
    void assemble_mg_matrix ();
    void assemble_right_hand_side ();
    void error ();
    void estimate ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;
    
    Triangulation<dim> triangulation;
    const MappingQ1<dim> mapping;
    const FiniteElement<dim>& fe;
    MGDoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity;
    SparseMatrix<double> matrix;
    Vector<double>       solution;
    Vector<double>       right_hand_side;

    MGLevelObject<SparsityPattern> mg_sparsity;
    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface;
    MGLevelObject<SparseMatrix<double> > mg_matrix;
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_up;
    MGLevelObject<SparseMatrix<double> > mg_matrix_dg_down;
};


template <int dim>
Step39<dim>::Step39(const FiniteElement<dim>& fe)
		:
		fe(fe),
		dof_handler(triangulation)
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
}


template <int dim>
void
Step39<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);
  unsigned int n_dofs = dof_handler.n_dofs();
  
  CompressedSparsityPattern c_sparsity(n_dofs);
  const DoFHandler<dim>& dof = dof_handler;
  DoFTools::make_flux_sparsity_pattern(dof, c_sparsity);
  sparsity.copy_from(c_sparsity);
  
  matrix.reinit(sparsity);
  
  solution.reinit(n_dofs);
  right_hand_side.reinit(n_dofs);

  const unsigned int n_levels = triangulation.n_levels();
  
  mg_matrix.resize(0, n_levels-1);
  mg_matrix.clear();
  mg_matrix_dg_up.resize(0, n_levels-1);
  mg_matrix_dg_up.clear();
  mg_matrix_dg_down.resize(0, n_levels-1);
  mg_matrix_dg_down.clear();

  mg_sparsity.resize(0, n_levels-1);
  mg_sparsity_dg_interface.resize(0, n_levels-1);
  
  for (unsigned int level=mg_sparsity.get_minlevel();
       level<=mg_sparsity.get_maxlevel();++level)
    {
      CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(level));
      CompressedSparsityPattern ci_sparsity;
      if (level>0)
	ci_sparsity.reinit(dof_handler.n_dofs(level-1), dof_handler.n_dofs(level));
      
      MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
      if (level>0)
	MGTools::make_flux_sparsity_pattern_edge(dof_handler, ci_sparsity, level);
      
      mg_sparsity[level].copy_from(c_sparsity);
      mg_matrix[level].reinit(mg_sparsity[level]);
      if (level>0)
	{
	  mg_sparsity_dg_interface[level].copy_from(ci_sparsity);
	  mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]);
	  mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]);
	}
    }
}


template <int dim>
void
Step39<dim>::assemble_matrix()
{
  const MatrixIntegrator<dim> local;
  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> >, MatrixIntegrator<dim> >
    integrator(local);
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  integrator.initialize_gauss_quadrature(n_gauss_points, n_gauss_points, n_gauss_points);
  UpdateFlags update_flags = update_values | update_gradients;
  integrator.add_update_flags(update_flags, true, true, true, true);

  integrator.initialize(matrix);
  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(integrator, fe, mapping);
  MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(), info_box, integrator);
}


template <int dim>
void
Step39<dim>::assemble_mg_matrix()
{
  const MatrixIntegrator<dim> local;
  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> >, MatrixIntegrator<dim> >
    integrator(local);
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  integrator.initialize_gauss_quadrature(n_gauss_points, n_gauss_points, n_gauss_points);
  UpdateFlags update_flags = update_values | update_gradients;
  integrator.add_update_flags(update_flags, true, true, true, true);

  integrator.initialize(mg_matrix);
  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(integrator, fe, mapping);
  MeshWorker::integration_loop(dof_handler.begin(), dof_handler.end(), info_box, integrator);
}


template <int dim>
void
Step39<dim>::assemble_right_hand_side()
{
  const RHSIntegrator<dim> local;
  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::ResidualSimple<Vector<double> >, RHSIntegrator<dim> >
  integrator(local);
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  integrator.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  integrator.add_update_flags(update_flags, true, true, true, true);

  NamedData<Vector<double>* > data;
  Vector<double>* rhs = &right_hand_side;
  data.add(rhs, "RHS");
  integrator.initialize(data);
  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(integrator, fe, mapping);
  MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(), info_box, integrator);
  
  right_hand_side *= -1.;
}


template <int dim>
void
Step39<dim>::solve()
{
  SolverControl control(1000, 1.e-12);
  SolverCG<Vector<double> > cg(control);

  GrowingVectorMemory<Vector<double> > mem;
  MGTransferPrebuilt<Vector<double> > mg_transfer;
  mg_transfer.build_matrices(dof_handler);
  FullMatrix<double> coarse_matrix;
  coarse_matrix.copy_from (mg_matrix[0]);
  MGCoarseGridHouseholder<double, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);
  typedef PreconditionSOR<SparseMatrix<double> > RELAXATION;
  MGSmootherRelaxation<SparseMatrix<double>, RELAXATION, Vector<double> >
    mg_smoother(mem);
  RELAXATION::AdditionalData smoother_data(1.);
//  RELAXATION::AdditionalData smoother_data(fe.dofs_per_cell, 1.);
  mg_smoother.initialize(mg_matrix, smoother_data);
  
				   // Do two smoothing steps per level
  mg_smoother.set_steps(2);
				   // Since the SOR method is not
				   // symmetric, but we use conjugate
				   // gradient iteration below, here
				   // is a trick to make the
				   // multilevel preconditioner a
				   // symmetric operator even for
				   // nonsymmetric smoothers.
  mg_smoother.set_symmetric(true);
  mg_smoother.set_variable(false);

				   // We must wrap our matrices in an
				   // object having the required
				   // multiplication functions.
  MGMatrix<SparseMatrix<double>, Vector<double> > mgmatrix(&mg_matrix);
  MGMatrix<SparseMatrix<double>, Vector<double> > mgdown(&mg_matrix_dg_down);
  MGMatrix<SparseMatrix<double>, Vector<double> > mgup(&mg_matrix_dg_up);

  
				   // Now, we are ready to set up the
				   // V-cycle operator and the
				   // multilevel preconditioner.
  Multigrid<Vector<double> > mg(dof_handler, mgmatrix,
				mg_coarse, mg_transfer,
				mg_smoother, mg_smoother);
  mg.set_edge_matrices(mgdown, mgup);
  mg.set_debug(0);
  mg_smoother.set_debug(0);
  
  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
    preconditioner(dof_handler, mg, mg_transfer);
  
  cg.solve(matrix, solution, right_hand_side, preconditioner);
}


template <int dim>
void
Step39<dim>::error()
{
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+2;
  Vector<double> cell_errors(triangulation.n_active_cells());
  
  QGauss<dim> quadrature(n_gauss_points);
  VectorTools::integrate_difference(mapping, dof_handler, solution, exact_solution,
				    cell_errors, quadrature, VectorTools::L2_norm);
  deallog << "Error " << cell_errors.l2_norm() << std::endl;
}


template <int dim>
void Step39<dim>::output_results (const unsigned int cycle) const
{
				   // Output of the solution in
				   // gnuplot format.
  std::string filename = "sol-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gnuplot";
  std::cout << "Writing solution to <" << filename << ">..."
	    << std::endl << std::endl;
  std::ofstream gnuplot_output (filename.c_str());
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");

  data_out.build_patches ();
  
  data_out.write_gnuplot(gnuplot_output);
}


template <int dim>
void
Step39<dim>::run(unsigned int n_steps)
{
  for (unsigned int s=0;s<n_steps;++s)
    {
      deallog << "Step " << s << std::endl;
      if (s != 0)
	triangulation.refine_global(1);
      
      deallog << "Triangulation "
	      << triangulation.n_active_cells() << " cells, "
	      << triangulation.n_levels() << " levels" << std::endl;
      
      setup_system();
      deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, levels";
      for (unsigned int l=0;l<triangulation.n_levels();++l)
	deallog << ' ' << dof_handler.n_dofs(l);
      deallog << std::endl;

      deallog << "Assemble matrix" << std::endl;
      assemble_matrix();
      deallog << "Assemble multilevel matrix" << std::endl;
      assemble_mg_matrix();
      deallog << "Assemble right hand side" << std::endl;
      assemble_right_hand_side();
      
      solve();
      error();
      output_results(s);
    }
}


int main()
{
  FE_DGQ<2> dgq1(1);
  Step39<2> test1(dgq1);
  test1.run(7);
}
