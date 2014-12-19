// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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



// Controls that the covariant matrix is calculated properly. It uses
// a Q1 finite element to calculate the scalar product of the gradient
// of a projected function (a monomial) with the tangential to the
// cell surface taken in the cell midpoint.  The result obtained is
// compared with the exact one in the <1,2> case.

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <string>

// all include files needed for the program

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <cmath>



std::ofstream logfile("output");

template <int dim, int spacedim>
void test(std::string filename)

{
  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  const QMidpoint<dim> q_midpoint;


  // finite elements used for the
  // projection
  const FE_Q<dim,spacedim> fe (1);
  DoFHandler<dim,spacedim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  FEValues<dim,spacedim> fe_values (fe, q_midpoint,
                                    update_values |
                                    update_JxW_values |
                                    update_cell_normal_vectors |
                                    update_gradients);

  // finite elements used for the
  // graphical representation with
  // data_out
  const FE_DGQ<dim,spacedim> fe_help(0);
  DoFHandler<dim,spacedim> dof_handler_help (triangulation);
  dof_handler_help.distribute_dofs (fe_help);
  FEValues<dim,spacedim> fe_values_help (fe_help, q_midpoint,
                                         update_cell_normal_vectors);

  deallog
      << "no. of cells "<< triangulation.n_cells() <<std::endl;
  deallog
      << "no. of dofs "<< dof_handler.n_dofs()<< std::endl;
  deallog
      << "no. of dofs per cell "<< fe.dofs_per_cell<< std::endl;
  deallog
      << "no. of help dofs "<< dof_handler_help.n_dofs()<< std::endl;
  deallog
      << "no. of help dofs per cell "<< fe_help.dofs_per_cell<< std::endl;



  //  definition of the exact function
  //  and calculation of the projected
  //  one
  Vector<double> projected_one(dof_handler.n_dofs());

//  Functions::CosineFunction<spacedim> cosine;

  Tensor<1,spacedim> exp;
  exp[0]=1;
  exp[1]=0;
  if (spacedim==3)
    exp[2]=0;
  Functions::Monomial<spacedim> monomial(exp);

  const QGauss<dim> quad(5);
  ConstraintMatrix constraints;
  constraints.close();
  VectorTools::project(dof_handler, constraints, quad, monomial, projected_one);


  // calculate its gradient

  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector< Point<spacedim> > cell_normals(q_midpoint.size());
  std::vector< Point<spacedim> > cell_tangentials(q_midpoint.size());
  std::vector<double> shape_directional_derivative(dofs_per_cell);
  Vector<double> projected_directional_derivative(triangulation.n_cells());

  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);

  typename DoFHandler<dim, spacedim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {

      fe_values.reinit(cell);
      cell-> get_dof_indices (local_dof_indices);
      cell_normals = fe_values.get_cell_normal_vectors();

      // The cell tangential is calculated
      // in the midpoint of the cell. For
      // the <2,3> case the tangential with
      // z component = 0 is chosen out of
      // the plane tangential to the
      // element surface.
      cell_tangentials[0][0] = cell_normals[0][1]
                               / sqrt( pow(cell_normals[0][0],2) + pow(cell_normals[0][1],2) );
      cell_tangentials[0][1] = -cell_normals[0][0]
                               / sqrt( pow(cell_normals[0][0],2) + pow(cell_normals[0][1],2) );
      if (spacedim == 3)
        cell_tangentials[0][2]=0.;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          shape_directional_derivative[i]=
            contract(
              fe_values.shape_grad(i,0),
              cell_tangentials[0]);

          // notice that the dof_index for
          // fe_dgq(0) is the same as that of
          // the cell
          projected_directional_derivative(cell->index())
          +=
            projected_one(local_dof_indices[i])
            *
            shape_directional_derivative[i];
        }

      deallog
          << "cell no. "
          << cell->index()<< "; "
          << "dir.deriv. "
          << projected_directional_derivative(cell->index())<< "; "<<std::endl;
      if (spacedim == 2)
        deallog
            << "exact solution "
            << cos( 2*numbers::PI*
                    (cell->index()+.5) / triangulation.n_cells() )
            << std::endl;


    }

  //  write graphical output
  DataOut<dim, DoFHandler<dim,spacedim> > dataout;
  dataout.attach_triangulation(triangulation);
  dataout.add_data_vector(projected_directional_derivative, "derivative");
  dataout.build_patches();
  dataout.write_vtk(logfile);
}



int main ()
{
  logfile.precision (4);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-12);

  deallog<<"Test <1,2>"<<std::endl;
  test<1,2>(SOURCE_DIR "/grids/circle_4.inp");

  deallog<<std::endl;

  deallog<<"Test <2,3>"<<std::endl;
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp");

//     test<2,3>(SOURCE_DIR "/grids/sphere_2.inp");
//     test<2,3>(SOURCE_DIR "/grids/sphere_3.inp");
//     test<2,3>(SOURCE_DIR "/grids/sphere_4.inp");

  return 0;
}

