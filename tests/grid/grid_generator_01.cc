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



// Test grid generation functions  in GridGenerator.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


template<int dim>
void test(std::ostream &out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  if (dim>2) p1[2] = 0.;
  Point<dim> p2;
  p2[0] = 3.;
  if (dim>1) p2[1] = 2.;
  if (dim>2) p2[2] = 4.;
  Point<dim> p3;
  p3[0] = 2.;
  if (dim>1) p3[1] = 1.;
  if (dim>2) p3[2] = 4.;

  GridOut go;
  GridOutFlags::XFig xfig_flags;
  xfig_flags.level_color = false;
  xfig_flags.fill_style = 25;

  go.set_flags(xfig_flags);

  GridOut::OutputFormat format = GridOut::gnuplot;
  if (dim==2) format = GridOut::xfig;
  if (dim==3) format = GridOut::dx;

  if (true)
    {
      deallog << "hyper_cube" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_cube(tr, 3., 7.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "subdivided_hyper_cube" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::subdivided_hyper_cube(tr, 3, 1., 7.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "hyper_rectangle" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_rectangle(tr, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim> tr;
      std::vector<unsigned int> sub(dim);
      sub[0] = 2;
      if (dim>1) sub[1] = 3;
      if (dim>2) sub[2] = 4;
      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, p2, true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim==2)
    {
      deallog << "parallelogram" << std::endl;
      Triangulation<dim> tr;
      Tensor<2,dim> corners;
      corners[0] = p1;
      if (dim>1) corners[1] = p2;
      if (dim>2) corners[2] = p3;
      GridGenerator::parallelogram(tr, corners, true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>1)
    {
      deallog << "enclosed_hyper_cube" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::enclosed_hyper_cube(tr, 3., 7., 1., true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>1)
    {
      deallog << "hyper_ball" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_ball(tr, p1, 3.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>1)
    {
      deallog << "cylinder" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::cylinder(tr, 1., 3.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>1)
    {
      deallog << "hyper_L" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_L(tr, -1., 1.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim==2)
    {
      deallog << "hyper_cube_slit" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_cube_slit(tr, -2., 2., true);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim==2)
    {
      deallog << "hyper_shell" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::hyper_shell(tr, p1, 4., 6.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>2)
    {
      deallog << "cylinder_shell" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::cylinder_shell(tr, 2., 5., 6.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim>1)
    {
      deallog << "half_hyper_ball" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::half_hyper_ball(tr, p1, 3.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
  if (dim==2)
    {
      deallog << "half_hyper_shell" << std::endl;
      Triangulation<dim> tr;
      GridGenerator::half_hyper_shell(tr, p1, 4., 6.);
      if (tr.n_cells() > 0)
        go.write(tr, out, format);
    }
}


int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("1d");
  test<1>(logfile);
  deallog.pop();
  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
