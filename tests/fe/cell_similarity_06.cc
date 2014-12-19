// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// since early 2009, the FEValues objects try to be more efficient by only
// recomputing things like gradients of shape functions if the cell on which
// we are is not a translation of the previous one. in this series of tests we
// make sure that this actually works the way it's supposed to be; in
// particular, if we create a mesh of two identical cells but one has a curved
// boundary, then they are the same if we use a Q1 mapping, but not a Q2
// mapping. so we test that the mass matrix we get from each of these cells is
// actually different in the latter case, but the same in the former
//
// the various tests cell_similarity_?? differ in the mappings and finite
// elements in use


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_c1.h>

#include <fstream>


bool equal (const FullMatrix<double> &m1,
            const FullMatrix<double> &m2)
{
  double d = 0, s = 0;
  for (unsigned int i=0; i<m1.m(); ++i)
    for (unsigned int j=0; j<m1.n(); ++j)
      {
        d += (m1(i,j)-m2(i,j)) * (m1(i,j)-m2(i,j));
        s += m1(i,j) * m1(i,j);
      }
  return (d<1e-8*s);
}


template<int dim>
void test (const Triangulation<dim> &tr)
{
  FE_Q<dim> fe(2);
  deallog << "FE=" << fe.get_name() << std::endl;

  MappingQ<dim> mapping(2);
  deallog << "Mapping=Q2" << std::endl;


  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (mapping, fe, quadrature,
                           update_values | update_gradients |
                           update_JxW_values);

  FullMatrix<double> mass_matrix[2], laplace_matrix[2];
  mass_matrix[0].reinit(fe.dofs_per_cell, fe.dofs_per_cell);
  mass_matrix[1].reinit(fe.dofs_per_cell, fe.dofs_per_cell);
  laplace_matrix[0].reinit(fe.dofs_per_cell, fe.dofs_per_cell);
  laplace_matrix[1].reinit(fe.dofs_per_cell, fe.dofs_per_cell);

  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      fe_values.reinit (cell);

      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
          for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
            {
              mass_matrix[cell->index()](i,j) += fe_values.shape_value (i,q) *
                                                 fe_values.shape_value (j,q) *
                                                 fe_values.JxW(q);
              laplace_matrix[cell->index()](i,j) += fe_values.shape_grad (i,q) *
                                                    fe_values.shape_grad (j,q) *
                                                    fe_values.JxW(q);
            }
    }

  // check what we expect for this mapping
  // about equality or inequality of the
  // matrices
  deallog << "Mass matrices "
          << (equal (mass_matrix[0], mass_matrix[1]) ? "are" : "are not")
          << " equal."
          << std::endl;
  deallog << "Laplace matrices "
          << (equal (laplace_matrix[0], laplace_matrix[1]) ? "are" : "are not")
          << " equal."
          << std::endl;

  for (unsigned int cell=0; cell<2; ++cell)
    {
      deallog << "cell=" << cell << std::endl;
      deallog << "mass_matrix:" << std::endl;
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
            deallog << mass_matrix[cell](i,j) << ' ';
          deallog << std::endl;
        }

      deallog << "laplace_matrix:" << std::endl;
      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
            deallog << laplace_matrix[cell](i,j) << ' ';
          deallog << std::endl;
        }
    }
}



template <int dim>
void test()
{
  Triangulation<dim> tr;
  Point<dim> p1 = Point<dim>();
  Point<dim> p2 = (dim == 2 ?
                   Point<dim>(2,1) :
                   Point<dim>(2,1,1));
  std::vector<unsigned int> subdivisions (dim, 1);
  subdivisions[0] = 2;
  GridGenerator::subdivided_hyper_rectangle(tr, subdivisions, p1, p2);

  static const HyperBallBoundary<dim> boundary(tr.begin_active()->center());
  tr.set_boundary (1, boundary);

  // set boundary id on cell 1
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    if (tr.begin_active()->at_boundary(f))
      tr.begin_active()->face(f)->set_boundary_indicator (1);

  test(tr);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test<2>();
  test<3>();
}
