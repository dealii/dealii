// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Multigrid for continuous finite elements using MeshWorker

#include "../tests.h"
#include <deal.II/base/logstream.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/output.h>
#include <deal.II/meshworker/loop.h>

#include <fstream>
#include <sstream>

using namespace dealii;
using namespace LocalIntegrators;

template <int dim>
class LaplaceMatrix : public MeshWorker::LocalIntegrator<dim>
{
public:
  LaplaceMatrix();
  virtual void cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void boundary(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
		    MeshWorker::IntegrationInfo<dim>& info1, MeshWorker::IntegrationInfo<dim>& info2) const;
};


template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix()
		:
		MeshWorker::LocalIntegrator<dim>(true, false, false)
{}


template <int dim>
void LaplaceMatrix<dim>::cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const
{
  AssertDimension (dinfo.n_matrices(), 1);  
  Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0));
}


template <int dim>
void LaplaceMatrix<dim>::boundary(MeshWorker::DoFInfo<dim>& /*dinfo*/,
				  typename MeshWorker::IntegrationInfo<dim>& /*info*/) const
{
//  const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
//  Laplace::nitsche_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0),
//  			  Laplace::compute_penalty(dinfo, dinfo, deg, deg));
}


template <int dim>
void LaplaceMatrix<dim>::face(
  MeshWorker::DoFInfo<dim>& /*dinfo1*/, MeshWorker::DoFInfo<dim>& /*dinfo2*/,
  MeshWorker::IntegrationInfo<dim>& /*info1*/, MeshWorker::IntegrationInfo<dim>& /*info2*/) const
{
//  const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
//  Laplace::ip_matrix(dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix, 
//		     dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
//		     info1.fe_values(0), info2.fe_values(0),
//		     Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
}

template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem (const unsigned int deg);
  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void assemble_multigrid (const bool &use_mw);
  void solve ();
  void refine_grid (const std::string& reftype);
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;
  FE_Q<dim>            fe;
  MGDoFHandler<dim>    mg_dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  ConstraintMatrix     hanging_node_constraints;
  ConstraintMatrix     constraints;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  const unsigned int degree;
  LaplaceMatrix<dim> matrix_integrator;

  MGLevelObject<SparsityPattern>       mg_sparsity_patterns;
  MGLevelObject<SparseMatrix<double> > mg_matrices;
  MGLevelObject<SparseMatrix<double> > mg_interface_in;
  MGLevelObject<SparseMatrix<double> > mg_interface_out;
  MGConstrainedDoFs                    mg_constrained_dofs;
};


template <int dim>
class Coefficient : public Function<dim>
{
public:
  Coefficient () : Function<dim>() {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  virtual void value_list (const std::vector<Point<dim> > &points,
                           std::vector<double>            &values,
                           const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
                                const unsigned int) const
{
//  if (p.square() < 0.5*0.5)
//    return 20;
//  else
    return 1;
}



template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
                                   std::vector<double>            &values,
                                   const unsigned int              component) const
{
  const unsigned int n_points = points.size();

  Assert (values.size() == n_points,
          ExcDimensionMismatch (values.size(), n_points));

  Assert (component == 0,
          ExcIndexRange (component, 0, 1));

  for (unsigned int i=0; i<n_points; ++i)
    values[i] = Coefficient<dim>::value (points[i]);
}


template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree)
  :
  triangulation (Triangulation<dim>::
                 limit_level_difference_at_vertices),
  fe (degree),
  mg_dof_handler (triangulation),
  degree(degree),
  matrix_integrator()
{}


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);
  deallog << "Number of degrees of freedom: "
          << mg_dof_handler.n_dofs();

  for (unsigned int l=0; l<triangulation.n_levels(); ++l)
    deallog << "   " << 'L' << l << ": "
            << mg_dof_handler.n_dofs(l);
  deallog  << std::endl;

  sparsity_pattern.reinit (mg_dof_handler.n_dofs(),
                           mg_dof_handler.n_dofs(),
                           mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (
    static_cast<const DoFHandler<dim>&>(mg_dof_handler),
    sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

  constraints.clear ();
  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (mg_dof_handler, constraints);
  DoFTools::make_hanging_node_constraints (mg_dof_handler, hanging_node_constraints);
  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  MappingQ1<dim> mapping;
  VectorTools::interpolate_boundary_values (mapping,
                                            mg_dof_handler,
                                            dirichlet_boundary,
                                            constraints);
  constraints.close ();
  hanging_node_constraints.close ();
  constraints.condense (sparsity_pattern);
  sparsity_pattern.compress();
  system_matrix.reinit (sparsity_pattern);

  mg_constrained_dofs.clear();
  mg_constrained_dofs.initialize(mg_dof_handler, dirichlet_boundary);
  const unsigned int n_levels = triangulation.n_levels();

  mg_interface_in.resize(0, n_levels-1);
  mg_interface_in.clear ();
  mg_interface_out.resize(0, n_levels-1);
  mg_interface_out.clear ();
  mg_matrices.resize(0, n_levels-1);
  mg_matrices.clear ();
  mg_sparsity_patterns.resize(0, n_levels-1);

  for (unsigned int level=0; level<n_levels; ++level)
    {
      CompressedSparsityPattern csp;
      csp.reinit(mg_dof_handler.n_dofs(level),
                 mg_dof_handler.n_dofs(level));
      MGTools::make_sparsity_pattern(mg_dof_handler, csp, level);

      mg_sparsity_patterns[level].copy_from (csp);

      mg_matrices[level].reinit(mg_sparsity_patterns[level]);
      mg_interface_in[level].reinit(mg_sparsity_patterns[level]);
      mg_interface_out[level].reinit(mg_sparsity_patterns[level]);
    }
}


template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  const QGauss<dim>  quadrature_formula(degree+1);

  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  typename MGDoFHandler<dim>::active_cell_iterator
  cell = mg_dof_handler.begin_active(),
  endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      coefficient.value_list (fe_values.get_quadrature_points(),
                              coefficient_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (coefficient_values[q_point] *
                                   fe_values.shape_grad(i,q_point) *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            1.0 *
                            fe_values.JxW(q_point));
          }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                              local_dof_indices,
                                              system_matrix, system_rhs);
    }
}


template <int dim>
void LaplaceProblem<dim>::assemble_multigrid (const bool& use_mw)
{
  if(use_mw == true)
  {
    mg_matrices = 0.;

    MappingQ1<dim> mapping;
    MeshWorker::IntegrationInfoBox<dim> info_box;
    UpdateFlags update_flags = update_values | update_gradients | update_hessians;
    info_box.add_update_flags_all(update_flags);
    info_box.initialize(fe, mapping);

    MeshWorker::DoFInfo<dim> dof_info(mg_dof_handler);

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
    assembler.initialize(mg_constrained_dofs);
    assembler.initialize(mg_matrices);
    assembler.initialize_interfaces(mg_interface_in, mg_interface_out);

    MeshWorker::integration_loop<dim, dim> (
        mg_dof_handler.begin(), mg_dof_handler.end(),
        dof_info, info_box, matrix_integrator, assembler);

    const unsigned int nlevels = triangulation.n_levels();
    for (unsigned int level=0;level<nlevels;++level)
    {
      for(unsigned int i=0; i<mg_dof_handler.n_dofs(level); ++i)
        if(mg_matrices[level].diag_element(i)==0)
          mg_matrices[level].set(i,i,1.);
    }
  }
  else
  {
    QGauss<dim>  quadrature_formula(1+degree);

    FEValues<dim> fe_values (fe, quadrature_formula,
        update_values   | update_gradients |
        update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const Coefficient<dim> coefficient;
    std::vector<double>    coefficient_values (n_q_points);

    std::vector<std::vector<bool> > interface_dofs
      = mg_constrained_dofs.get_refinement_edge_indices ();
    std::vector<std::vector<bool> > boundary_interface_dofs
      = mg_constrained_dofs.get_refinement_edge_boundary_indices ();

    std::vector<ConstraintMatrix> boundary_constraints (triangulation.n_levels());
    std::vector<ConstraintMatrix> boundary_interface_constraints (triangulation.n_levels());
    for (unsigned int level=0; level<triangulation.n_levels(); ++level)
    {
      boundary_constraints[level].add_lines (interface_dofs[level]);
      boundary_constraints[level].add_lines (mg_constrained_dofs.get_boundary_indices()[level]);
      boundary_constraints[level].close ();

      boundary_interface_constraints[level]
        .add_lines (boundary_interface_dofs[level]);
      boundary_interface_constraints[level].close ();
    }

    typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
             endc = mg_dof_handler.end();

    for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      fe_values.reinit (cell);

      coefficient.value_list (fe_values.get_quadrature_points(),
          coefficient_values);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (coefficient_values[q_point] *
                fe_values.shape_grad(i,q_point) *
                fe_values.shape_grad(j,q_point) *
                fe_values.JxW(q_point));

      cell->get_mg_dof_indices (local_dof_indices);

      boundary_constraints[cell->level()]
        .distribute_local_to_global (cell_matrix,
            local_dof_indices,
            mg_matrices[cell->level()]);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          if ( !(interface_dofs[cell->level()][local_dof_indices[i]]==true &&
                interface_dofs[cell->level()][local_dof_indices[j]]==false))
            cell_matrix(i,j) = 0;

      boundary_interface_constraints[cell->level()]
        .distribute_local_to_global (cell_matrix,
            local_dof_indices,
            mg_interface_in[cell->level()]);
    }
  }
}


template <int dim>
void LaplaceProblem<dim>::solve ()
{
  MGTransferPrebuilt<Vector<double> > mg_transfer(hanging_node_constraints, mg_constrained_dofs);
  mg_transfer.build_matrices(mg_dof_handler);

  FullMatrix<double> coarse_matrix;
  coarse_matrix.copy_from (mg_matrices[0]);
  MGCoarseGridHouseholder<> coarse_grid_solver;
  coarse_grid_solver.initialize (coarse_matrix);

  typedef PreconditionSOR<SparseMatrix<double> > Smoother;
  MGSmootherRelaxation<SparseMatrix<double>, Smoother, Vector<double> >
  mg_smoother;
  mg_smoother.initialize(mg_matrices);
  mg_smoother.set_steps(2);
  mg_smoother.set_symmetric(true);

  MGMatrix<> mg_matrix(&mg_matrices);
  MGMatrix<> mg_interface_up(&mg_interface_in);
  MGMatrix<> mg_interface_down(&mg_interface_in);

  Multigrid<Vector<double> > mg(mg_dof_handler,
                                mg_matrix,
                                coarse_grid_solver,
                                mg_transfer,
                                mg_smoother,
                                mg_smoother);
  mg.set_edge_matrices(mg_interface_down, mg_interface_up);

  PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double> > >
  preconditioner(mg_dof_handler, mg, mg_transfer);

  SolverControl solver_control (1000, 1e-12);
  SolverCG<>    cg (solver_control);

  solution = 0;

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);
  constraints.distribute (solution);

  deallog << "   " << solver_control.last_step()
          << " CG iterations needed to obtain convergence."
          << std::endl;
}


template <int dim>
void LaplaceProblem<dim>::refine_grid (const std::string& reftype)
{
  bool cell_refined = false;
  if (reftype == "center" || !cell_refined)
    {
      for (typename Triangulation<dim>::active_cell_iterator
	     cell = triangulation.begin_active();
	   cell != triangulation.end(); ++cell)
	for (unsigned int vertex=0;
	     vertex < GeometryInfo<dim>::vertices_per_cell;
	     ++vertex)
	  {
	    {
	      const Point<dim> p = cell->vertex(vertex);
	      const Point<dim> origin = (dim == 2 ?
					 Point<dim>(0,0) :
					 Point<dim>(0,0,0));
	      const double dist = p.distance(origin);
	      if(dist<0.25/numbers::PI)
		{
		  cell->set_refine_flag ();
		  cell_refined = true;
		  break;
		}
	    }
	  }
    }
  if (reftype=="global" || !cell_refined)
    triangulation.refine_global(1);
  else
    triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
           << cycle
           << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}


template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube (triangulation, -1, 1);

//          static const HyperBallBoundary<dim> boundary;
//          triangulation.set_boundary (0, boundary);

          triangulation.refine_global (1);
        }
      else
        refine_grid ("center");


      deallog << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

      setup_system ();

      deallog << "   Number of degrees of freedom: "
              << mg_dof_handler.n_dofs()
              << " (by level: ";
      for (unsigned int level=0; level<triangulation.n_levels(); ++level)
        deallog << mg_dof_handler.n_dofs(level)
                << (level == triangulation.n_levels()-1
                    ? ")" : ", ");
      deallog << std::endl;

      assemble_system ();
      assemble_multigrid (true);

      solve ();
//      output_results (cycle);
    }
}


int main ()
{
  initlog();

  try
    {
      deallog.depth_console (0);

      LaplaceProblem<2> laplace_problem(1);
      laplace_problem.run ();
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
