// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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



// like hp/matrices, but with different fe objects


#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>



template <int dim>
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
check_boundary (const hp::DoFHandler<dim> &dof,
                const hp::MappingCollection<dim>    &mapping)
{
  MySquareFunction<dim> coefficient;
  typename FunctionMap<dim>::type function_map;
  function_map[0] = &coefficient;

  hp::QCollection<dim-1> face_quadrature;
  for (unsigned int i=1; i<7-dim; ++i)
    face_quadrature.push_back (QGauss<dim-1>(3+i));

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
check_boundary (const hp::DoFHandler<1> &,
                const hp::MappingCollection<1> &)
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
  tr.reset_manifold(0);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  // create a system element composed
  // of one Q1 and one Q2 element
  hp::FECollection<dim> element;
  for (unsigned int i=1; i<7-dim; ++i)
    element.push_back (FESystem<dim> (FE_Q<dim>(QIterated<1>(QTrapez<1>(),i)), 1,
                                      FE_Q<dim>(QIterated<1>(QTrapez<1>(),i+1)), 1));
  hp::DoFHandler<dim> dof(tr);
  for (typename hp::DoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active(); cell!=dof.end(); ++cell)
    cell->set_active_fe_index (Testing::rand() % element.size());

  dof.distribute_dofs(element);

  // use a more complicated mapping
  // of the domain and a quadrature
  // formula suited to the elements
  // we have here
  hp::MappingCollection<dim> mapping;
  for (unsigned int i=1; i<7-dim; ++i)
    mapping.push_back (MappingQ<dim>(i+1));

  hp::QCollection<dim> quadrature;
  for (unsigned int i=1; i<7-dim; ++i)
    quadrature.push_back (QGauss<dim>(3+i));

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
  initlog();
  deallog.get_file_stream() << std::setprecision(12);

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
