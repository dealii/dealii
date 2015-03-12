// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2015 by the deal.II authors
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



/* Author: Wolfgang Bangerth, University of Heidelberg, 2001 */



#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>


template<int dim>
class MySquareFunction : public Function<dim>
{
public:
  MySquareFunction () : Function<dim>(2) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return (component+1)*p.square();
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
    values(1) = value(p,1);
  }
};




template <int dim>
void
check_boundary (const DoFHandler<dim> &dof,
                const Mapping<dim>    &mapping)
{
  MySquareFunction<dim> coefficient;
  typename FunctionMap<dim>::type function_map;
  function_map[0] = &coefficient;

  QGauss<dim-1> face_quadrature(6);

  std::vector<types::global_dof_index> dof_to_boundary_mapping;
  DoFTools::map_dof_to_boundary_indices (dof,
                                         dof_to_boundary_mapping);

  SparsityPattern sparsity(dof.n_boundary_dofs(function_map),
                           dof.max_couplings_between_boundary_dofs());
  DoFTools::make_boundary_sparsity_pattern (dof,
                                            function_map,
                                            dof_to_boundary_mapping,
                                            sparsity);
  sparsity.compress ();

  SparseMatrix<double> matrix;
  matrix.reinit (sparsity);

  Vector<double> rhs (dof.n_boundary_dofs(function_map));
  MatrixTools::
  create_boundary_mass_matrix (mapping, dof,
                               face_quadrature, matrix,
                               function_map, rhs,
                               dof_to_boundary_mapping,
                               &coefficient);

  // since we only generate
  // output with two digits after
  // the dot, and since matrix
  // entries are usually in the
  // range of 1 or below,
  // multiply matrix by 100 to
  // make test more sensitive
  matrix *= 100;

  // finally write out matrix
  matrix.print (deallog.get_file_stream());
}



void
check_boundary (const DoFHandler<1> &,
                const Mapping<1> &)
{}




template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  // create a system element composed
  // of one Q1 and one Q2 element
  FESystem<dim> element(FE_Q<dim>(1), 1,
                        FE_Q<dim>(2), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  MappingQ<dim> mapping (3);
  QGauss<dim> quadrature(6);

  // create sparsity pattern. note
  // that different components should
  // not couple, so use pattern
  SparsityPattern sparsity (dof.n_dofs(), dof.n_dofs());
  Table<2,DoFTools::Coupling> mask (2, 2);
  mask(0,0) = mask(1,1) = DoFTools::always;
  mask(0,1) = mask(1,0) = DoFTools::none;
  DoFTools::make_sparsity_pattern (dof, mask, sparsity);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof, constraints);
  constraints.close ();
  constraints.condense (sparsity);
  sparsity.compress ();

  SparseMatrix<double> matrix;

  Functions::ExpFunction<dim> coefficient;

  typename FunctionMap<dim>::type function_map;
  function_map[0] = &coefficient;

  for (unsigned int test=0; test<2; ++test)
    {
      matrix.reinit(sparsity);
      switch (test)
        {
        case 0:
          MatrixTools::
          create_mass_matrix (mapping, dof,
                              quadrature, matrix, &coefficient);
          break;
        case 1:
          MatrixTools::
          create_laplace_matrix (mapping, dof,
                                 quadrature, matrix, &coefficient);
          break;
        default:
          Assert (false, ExcInternalError());
        };

      // since we only generate
      // output with two digits after
      // the dot, and since matrix
      // entries are usually in the
      // range of 1 or below,
      // multiply matrix by 100 to
      // make test more sensitive
      for (SparseMatrix<double>::const_iterator p=matrix.begin();
	   p!=matrix.end(); ++p)
	deallog.get_file_stream() << p->value() * 100
				  << std::endl;
    };

  if (dim > 1)
    check_boundary (dof, mapping);
}



int main ()
{
  std::ofstream logfile ("output");
  logfile << std::setprecision (2);
  logfile << std::fixed;
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
