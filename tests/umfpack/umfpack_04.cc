// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


// test the umfpack sparse direct solver on a block mass matrix

#include "../tests.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

template <int dim, typename MATRIX, typename VECTOR>
void assemble_laplace (MATRIX &B, VECTOR &bb, DoFHandler<dim> &dof_handler, FiniteElement<dim> &fe)
{
  QGauss<dim>  quadrature_formula(2);
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values | update_gradients | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

      unsigned int comp_i, comp_j;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          comp_i = fe.system_to_component_index(i).first;
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              comp_j = fe.system_to_component_index(j).first;
              if (comp_i == comp_j)
                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                  cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                       fe_values.shape_grad (j, q_point) *
                                       fe_values.JxW (q_point));
            }
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                              1 *
                              fe_values.JxW (q_point));
        }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          B.add (local_dof_indices[i],
                 local_dof_indices[j],
                 cell_matrix(i,j));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        bb(local_dof_indices[i]) += cell_rhs(i);
    }
}



template <int dim>
void test ()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria,0,1);
  tria.refine_global (1);

  // destroy the uniformity of the matrix by
  // refining one cell
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global(7-2*dim);

  FESystem<dim> fe (FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);
  DoFRenumbering::component_wise (dof_handler);

  deallog << "Number of dofs = " << dof_handler.n_dofs() << std::endl;

  std::vector<types::global_dof_index> size (dim);
  DoFTools::count_dofs_per_component(dof_handler, size);

  BlockSparsityPattern b_sparsity_pattern;
  SparsityPattern sparsity_pattern;

  b_sparsity_pattern.reinit(size.size(),size.size());
  for (unsigned int k=0; k<size.size(); ++k)
    for (unsigned int l=0; l<size.size(); ++l)
      b_sparsity_pattern.block(k,l).reinit(size[k],
                                           size[l],
                                           dof_handler.max_couplings_between_dofs());
  b_sparsity_pattern.collect_sizes();

  sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.max_couplings_between_dofs());

  DoFTools::make_sparsity_pattern (dof_handler, b_sparsity_pattern);
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

  b_sparsity_pattern.compress();
  sparsity_pattern.compress ();

  BlockSparseMatrix<double> Bb;
  Bb.reinit (b_sparsity_pattern);

  SparseMatrix<double> B;
  B.reinit(sparsity_pattern);

  BlockVector<double> bb(size), bx(size);
  Vector<double> b(dof_handler.n_dofs()), x(dof_handler.n_dofs()),
         ub(dof_handler.n_dofs()), ubb(dof_handler.n_dofs());

  assemble_laplace(Bb, bb, dof_handler, fe);
  assemble_laplace(B, b, dof_handler, fe);

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(size.size()),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values, Bb, bx, bb);
  MatrixTools::apply_boundary_values (boundary_values, B, x, b);

  ub = b;
  ubb = bb;

  SolverControl control (1000, 1e-13);

  SolverCG<BlockVector<double> > bcg(control, SolverCG<BlockVector<double> >::AdditionalData());
  bcg.solve(Bb, bx, bb, PreconditionIdentity());

  SolverCG<Vector<double> > cg(control, SolverCG<Vector<double> >::AdditionalData());
  cg.solve(B, x, b, PreconditionIdentity());

  deallog << "Sparse Factorization" << std::endl;
  SparseDirectUMFPACK umfpack, umfpackb;
  umfpack.factorize(B);
  umfpack.solve(ub);
  deallog << "absolute norms = "
          << x.l2_norm() << ' ' << ub.l2_norm()
          << std::endl;
  x -= ub;

  deallog << "relative norm distance = "
          << x.l2_norm() / ub.l2_norm()
          << std::endl;
  Assert (x.l2_norm() / ub.l2_norm() < 1e-12,
          ExcInternalError());

  deallog << "Block Sparse Factorization" << std::endl;
  umfpackb.factorize(Bb);
  // umfpackb.factorize(B);
  umfpackb.solve(ubb);
  deallog << "absolute norms = "
          << bx.l2_norm() << ' ' << ubb.l2_norm()
          << std::endl;
  b = bx;
  b -= ubb;

  deallog << "relative norm distance = "
          << b.l2_norm() / ubb.l2_norm()
          << std::endl;
  Assert (b.l2_norm() / ubb.l2_norm() < 1e-12,
          ExcInternalError());
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();
}
