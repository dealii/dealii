// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


// Refinement base: perform REF/dim refinements
// REF 21 needs 9GBytes on a 64 Bit architecture

#ifdef DEBUG
#define REF 6
#else
#define REF 21
#endif

#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;

template <int dim>
void indices (const DoFHandler<dim> &dof, unsigned int repeat)
{
  typedef typename DoFHandler<dim>::active_cell_iterator I;

  std::vector<unsigned int> dofs(dof.get_fe().dofs_per_cell);
  const I end = dof.end();

  for (unsigned int k=0; k<repeat; ++k)
    for (I i=dof.begin_active(); i!=end; ++i)
      i->get_dof_indices(dofs);
}


template <int dim>
void fevalues (const DoFHandler<dim> &dof,
               UpdateFlags updates)
{
  typedef typename DoFHandler<dim>::active_cell_iterator I;
  const I end = dof.end();

  QGauss<dim> quadrature(5);
  MappingQ1<dim> mapping;
  FEValues<dim> fe(mapping, dof.get_fe(), quadrature, updates);

  for (I i=dof.begin_active(); i!=end; ++i)
    fe.reinit(i);
}


template <int dim>
void fefacevalues (const DoFHandler<dim> &dof,
                   UpdateFlags updates)
{
  typedef typename DoFHandler<dim>::active_cell_iterator I;
  const I end = dof.end();

  QGauss<dim-1> quadrature(5);
  MappingQ1<dim> mapping;
  FEFaceValues<dim> fe(mapping, dof.get_fe(), quadrature, updates);

  for (I i=dof.begin_active(); i!=end; ++i)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      fe.reinit(i, f);
}


template <int dim>
void check_mapping (const DoFHandler<dim> &dof)
{
  deallog.push("cell");
  fevalues(dof,  update_q_points);
  deallog << "qpoints" << std::endl;
  fevalues(dof, update_JxW_values);
  deallog << "JxW" << std::endl;
  fevalues(dof,  update_q_points | update_JxW_values);
  deallog << "qpoints|JxW" << std::endl;
  deallog.pop();

  deallog.push("face");
  fefacevalues(dof,  update_quadrature_points);
  deallog << "qpoints" << std::endl;
  fefacevalues(dof, update_JxW_values);
  deallog << "JxW" << std::endl;
  fefacevalues(dof, update_normal_vectors);
  deallog << "normals" << std::endl;
  fefacevalues(dof,  update_q_points | update_JxW_values | update_normal_vectors);
  deallog << "qpoints|JxW|normals" << std::endl;
  deallog.pop();

}


template <int dim>
void check_values (const DoFHandler<dim> &dof)
{
  indices(dof, 100);
  deallog << "Index*100" << std::endl;
  deallog.push("cell");
  fevalues(dof,  update_values);
  deallog << "values" << std::endl;
  fevalues(dof,  update_gradients);
  deallog << "gradients" << std::endl;
  fevalues(dof,  update_values | update_JxW_values);
  deallog << "values|JxW" << std::endl;
  fevalues(dof,  update_values | update_gradients
           | update_q_points | update_JxW_values);
  deallog << "values|gradients|qpoints|JxW" << std::endl;
//  fevalues(dof,  update_values | update_gradients | update_second_derivatives
//     | update_q_points | update_JxW_values);
//  deallog << "values|gradients|2nds|qpoints|JxW" << std::endl;
  deallog.pop();

  deallog.push("face");
  fefacevalues(dof,  update_values);
  deallog << "values" << std::endl;
  fefacevalues(dof,  update_gradients);
  deallog << "gradients" << std::endl;
  fefacevalues(dof,  update_values | update_JxW_values);
  deallog << "values|JxW" << std::endl;
  fefacevalues(dof,  update_values | update_gradients
               | update_q_points | update_JxW_values | update_normal_vectors);
  deallog << "values|gradients|qpoints|JxW|normals" << std::endl;
//  fefacevalues(dof,  update_values | update_gradients | update_second_derivatives
//     | update_q_points | update_JxW_values);
//  deallog << "values|gradients|2nds|qpoints|JxW" << std::endl;
  deallog.pop();
}


template <int dim>
void check_q ()
{
  deallog << "Mesh" << std::endl;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(REF/dim);

  for (unsigned int i=1; i<5; ++i)
    {
      FE_Q<dim> q(i);
      deallog.push(q.get_name());
      deallog << "Dofs per cell " << q.dofs_per_cell << std::endl;
      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(q);
      deallog << "Dofs " << dof.n_dofs() << std::endl;

      if (i==1)
        check_mapping(dof);
      check_values(dof);
      deallog.pop();
    }
}


template <int dim>
void check_sys ()
{
  deallog << "Mesh" << std::endl;
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(REF/dim);

  FE_Q<dim> q(1);

  for (unsigned int i=1; i<5; ++i)
    {
      FESystem<dim> fe(q,1<<i);
      deallog.push(fe.get_name());
      deallog << "Dofs per cell " << fe.dofs_per_cell << std::endl;
      DoFHandler<dim> dof(tr);
      dof.distribute_dofs(fe);
      deallog << "Dofs " << dof.n_dofs() << std::endl;

      if (i==1)
        check_mapping(dof);
      check_values(dof);
      deallog.pop();
    }
}


int main()
{
  std::ofstream out("dof_handler_timing/output");
  deallog.attach(out);
  deallog.log_execution_time(true);
  deallog.log_time_differences(true);
  check_q<2>();
  check_sys<2>();
  check_q<3>();
  check_sys<3>();
}

