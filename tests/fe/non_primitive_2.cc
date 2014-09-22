// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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


// assemble a matrix for the stokes equation in different ways for
// non-primitive elements. this is the counterpart to the non_primitive_1
// program, which did the same thing for actually primitive finite elements
//
// of course, it makes absolutely no sense to work the Stokes equation
// with a Nedelec element, but this is just to test the library, no?

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>


// note: create_stokes_matrix_1 from non_primitive_1 does not work here,
// since the elements in use are actually non-primitive

// create the matrix in the simple
// way that is necessary when you
// want to use non-primitive shape
// functions
template <int dim>
void
create_stokes_matrix_2 (const DoFHandler<dim> &dof_handler,
                        SparseMatrix<double>  &A)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);

  QGauss<dim> quadrature (3);
  const unsigned int n_q_points = quadrature.size();

  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients |
                           update_JxW_values);

  const double nu = 3.14159265358e-2;

  for (; cell!=endc; ++cell)
    {
      local_matrix = 0;
      fe_values.reinit (cell);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int comp_i=0; comp_i<fe.n_components(); ++comp_i)
          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
            for (unsigned int comp_j=0; comp_j<fe.n_components(); ++comp_j)
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  // velocity-velocity coupling?
                  if ((comp_i<dim) && (comp_j<dim))
                    if (comp_i == comp_j)
                      local_matrix(i,j)
                      += (nu *
                          (fe_values.shape_grad_component(i,q,comp_i) *
                           fe_values.shape_grad_component(j,q,comp_j)  ) *
                          fe_values.JxW(q));

                  // velocity-pressure coupling
                  if ((comp_i<dim) && (comp_j==dim))
                    local_matrix(i,j)
                    += (-fe_values.shape_grad_component(i,q,comp_i)[comp_i] *
                        fe_values.shape_value_component(j,q,comp_j)        *
                        fe_values.JxW(q));

                  // pressure-velocity coupling
                  if ((comp_i==dim) && (comp_j<dim))
                    local_matrix(i,j)
                    += (fe_values.shape_value_component(i,q,comp_i) *
                        fe_values.shape_grad_component(j,q,comp_j)[comp_j] *
                        fe_values.JxW(q));
                };


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          A.add (local_dof_indices[i],
                 local_dof_indices[j],
                 local_matrix(i,j));
    };
}



// create the matrix in a way that is
// necessary when you want to use
// non-primitive shape
// functions. compared to the second
// possibility used above, use some
// optimizations
template <int dim>
void
create_stokes_matrix_3 (const DoFHandler<dim> &dof_handler,
                        SparseMatrix<double>  &A)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);

  QGauss<dim> quadrature (3);
  const unsigned int n_q_points = quadrature.size();

  FEValues<dim> fe_values (fe, quadrature,
                           update_values | update_gradients |
                           update_JxW_values);

  const double nu = 3.14159265358e-2;

  for (; cell!=endc; ++cell)
    {
      local_matrix = 0;
      fe_values.reinit (cell);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int comp_i=0; comp_i<fe.n_components(); ++comp_i)
          if (fe.get_nonzero_components(i)[comp_i] == true)
            for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
              for (unsigned int comp_j=0; comp_j<fe.n_components(); ++comp_j)
                if (fe.get_nonzero_components(j)[comp_j] == true)
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      // velocity-velocity coupling?
                      if ((comp_i<dim) && (comp_j<dim))
                        if (comp_i == comp_j)
                          local_matrix(i,j)
                          += (nu *
                              (fe_values.shape_grad_component(i,q,comp_i) *
                               fe_values.shape_grad_component(j,q,comp_j)  ) *
                              fe_values.JxW(q));

                      // velocity-pressure coupling
                      if ((comp_i<dim) && (comp_j==dim))
                        local_matrix(i,j)
                        += (-fe_values.shape_grad_component(i,q,comp_i)[comp_i] *
                            fe_values.shape_value_component(j,q,comp_j)        *
                            fe_values.JxW(q));

                      // pressure-velocity coupling
                      if ((comp_i==dim) && (comp_j<dim))
                        local_matrix(i,j)
                        += (fe_values.shape_value_component(i,q,comp_i) *
                            fe_values.shape_grad_component(j,q,comp_j)[comp_j] *
                            fe_values.JxW(q));
                    };


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          A.add (local_dof_indices[i],
                 local_dof_indices[j],
                 local_matrix(i,j));
    };
}




template <int dim>
void
test (const unsigned int p)
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube (triangulation);

  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();

  deallog << "dim=" << dim
          << ", n_cells=" << triangulation.n_active_cells()
          << std::endl;

  FESystem<dim> fe (FE_Nedelec<dim>(p), 1,
                    FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);

  SparsityPattern sparsity(dof_handler.n_dofs(),
                           dof_handler.n_dofs());

  // have a mask that indicates that
  // for the stokes equation the
  // pressure does not couple to
  // itself
//   std::vector<std::vector<bool> > mask (dim+1, std::vector<bool> (dim+1, true));
//   mask[dim][dim] = false;

  DoFTools::make_sparsity_pattern (dof_handler, /*mask,*/ sparsity);
  sparsity.compress ();

  SparseMatrix<double> A2 (sparsity);
  SparseMatrix<double> A3 (sparsity);

  create_stokes_matrix_2 (dof_handler, A2);
  create_stokes_matrix_3 (dof_handler, A3);

  // write out the contents of the
  // matrix and compare for equality
  // with the other matrices. to
  // reduce the amount of data
  // written out a little bit, only
  // write every so-many-th element
  for (unsigned int i=0; i<A2.n_nonzero_elements(); ++i)
    {
      if (i % (dim*dim*dim) == 0)
        deallog << i << ' ' << A2.global_entry(i) << std::endl;
      Assert (A3.global_entry(i) == A2.global_entry(i),
              ExcInternalError());
    };
}



int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  deallog << "Degree 0: " << std::endl;
  test<2> (0);
  test<3> (0);
  deallog << "Degree 1: " << std::endl;
  test<2> (1);
  test<3> (1);
}




