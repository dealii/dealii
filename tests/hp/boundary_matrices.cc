// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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



// Like hp/matrices, but only the boundary mass matrices for vector-valued
// (non-primitive) hp objects are being tested.


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
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
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
  MySquareFunction () : Function<dim>(dim+1) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return (component+1)*p.square();
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    for (unsigned int i=0; i<dim+1; ++i)
      values(i) = value(p,i);
  }
};



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

  // Create a system element composed
  // of one RT1 and one DGQ1
  hp::FECollection<dim> element;
  element.push_back (FESystem<dim> (FE_RaviartThomasNodal<dim>(1), 1,
                                    FE_DGQ<dim>(1), 1));
  hp::DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  MySquareFunction<dim> coefficient;
  typename FunctionMap<dim>::type function_map;
  function_map[0] = &coefficient;

  hp::QCollection<dim-1> face_quadrature;
  face_quadrature.push_back (QGauss<dim-1>(6));

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
  create_boundary_mass_matrix (dof,
                               face_quadrature, matrix,
                               function_map, rhs,
                               dof_to_boundary_mapping,
                               &coefficient);

  // Multiply matrix by 100 to
  // make test more sensitive
  matrix *= 100;

  // Write out matrix
  matrix.print (deallog.get_file_stream());
}



int main ()
{
  std::ofstream logfile ("output");
  logfile.precision (2);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
