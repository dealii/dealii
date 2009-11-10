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

				 // Include files for setting up the
				 // mesh
#include <grid/grid_generator.h>

				 // Include files for FinitElement
				 // classes and DoFHandler.
#include <fe/fe_dgq.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_handler.h>

				 // The include files for using the
				 // MeshWorker framework
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_loop.h>

				 // Finally, we take our exact
				 // solution from the library.
#include <base/function_lib.h>

#include <iostream>
#include <fstream>

				 // All classes of the deal.II library
				 // are in the namespace dealii. In
				 // order to save typing, we tell the
				 // compiler to search names in there
				 // as well.
using namespace dealii;


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
	local_matrix(i,j) += fe.shape_grad(i,k) * fe.shape_grad(j,k)
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
  const FEFaceValuesBase<dim>& fe = info1.fe();
  FullMatrix<double>& matrix_v1u1 = info1.M1[0].matrix;
  FullMatrix<double>& matrix_v1u2 = info1.M2[0].matrix;
  FullMatrix<double>& matrix_v2u1 = info2.M2[0].matrix;
  FullMatrix<double>& matrix_v2u2 = info2.M1[0].matrix;
  
  const unsigned int deg = fe.get_fe().tensor_degree();
  const double penalty1 = deg * (deg+1) * info1.face->measure() / info1.cell->measure();
  const double penalty2 = deg * (deg+1) * info2.face->measure() / info2.cell->measure();
  const double penalty = penalty1 + penalty2;
  
  for (unsigned k=0;k<fe.n_quadrature_points;++k)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	{
	  matrix_v1u1(i,j) += (fe.shape_value(i,k) * penalty * fe.shape_value(j,k)
			       - (fe.normal_vector(k) * fe.shape_grad(i,k)) * fe.shape_value(j,k)
			       - (fe.normal_vector(k) * fe.shape_grad(j,k)) * fe.shape_value(i,k))
			      * .5 * fe.JxW(k);
	  matrix_v1u2(i,j) += (-fe.shape_value(i,k) * penalty * fe.shape_value(j,k)
			       + (fe.normal_vector(k) * fe.shape_grad(i,k)) * fe.shape_value(j,k)
			       - (fe.normal_vector(k) * fe.shape_grad(j,k)) * fe.shape_value(i,k))
			      * .5 * fe.JxW(k);
	  matrix_v2u1(i,j) += (-fe.shape_value(i,k) * penalty * fe.shape_value(j,k)
			       - (fe.normal_vector(k) * fe.shape_grad(i,k)) * fe.shape_value(j,k)
			       + (fe.normal_vector(k) * fe.shape_grad(j,k)) * fe.shape_value(i,k))
			      * .5 * fe.JxW(k);
	  matrix_v2u2(i,j) += (fe.shape_value(i,k) * penalty * fe.shape_value(j,k)
			       + (fe.normal_vector(k) * fe.shape_grad(i,k)) * fe.shape_value(j,k)
			       + (fe.normal_vector(k) * fe.shape_grad(j,k)) * fe.shape_value(i,k))
			      * .5 * fe.JxW(k);
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
  
  static Functions::LSingularityFunction lshaped;
  std::vector<double> boundary_values(fe.n_quadrature_points);
  lshaped.value_list(fe.get_quadrature_points(), boundary_values);
  
  const unsigned int deg = fe.get_fe().tensor_degree();
  const double penalty = 2. * deg * (deg+1) * info.face->measure() / info.cell->measure();
  
  for (unsigned k=0;k<fe.n_quadrature_points;++k)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      local_vector(i) += (-fe.shape_value(i,k) * penalty * boundary_values[k]
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
    void assemble_right_hand_side ();
    double estimate ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;
    
    Triangulation<dim> triangulation;
    const MappingQ1<dim> mapping;
    const FiniteElement<dim>& fe;
    MGDoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> matrix;
    Vector<double>       solution;
    Vector<double>       right_hand_side;
};


template <int dim>
Step39<dim>::Step39(const FiniteElement<dim>& fe)
		:
		fe(fe),
		dof_handler(triangulation)
{
  GridGenerator::hyper_L(triangulation);
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
  sparsity_pattern.copy_from(c_sparsity);
  
  matrix.reinit(sparsity_pattern);
  
  solution.reinit(n_dofs);
  right_hand_side.reinit(n_dofs);
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
Step39<dim>::assemble_right_hand_side()
{
  const RHSIntegrator<dim> local;
  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::ResidualSimple<Vector<double> >, RHSIntegrator<dim> >
  integrator(local);
  const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
  integrator.initialize_gauss_quadrature(1, n_gauss_points, n_gauss_points);
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  integrator.add_update_flags(update_flags, true, true, true, true);

  NamedData<Vector<double>* > data;
  Vector<double>* rhs = &right_hand_side;
  data.add(rhs, "RHS");
  integrator.initialize(data);
  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(integrator, fe, mapping);
  MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(), info_box, integrator);
}


template <int dim>
void
Step39<dim>::solve()
{
  SolverControl control(1000, 1.e-12);
  SolverCG<Vector<double> > cg(control);
  cg.solve(matrix, solution, right_hand_side, PreconditionIdentity());
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

      assemble_matrix();
      assemble_right_hand_side();
      solve();
    }
}


int main()
{
  FE_DGQ<2> dgq1(1);
  Step39<2> test1(dgq1);
  test1.run(7);
}
