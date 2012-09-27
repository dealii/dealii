//--------------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------------------

// Output the results of all begin and end functions.

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>

// One function for testing all begin, last and end functions for
// cells, faces, lines, quads and hexes, respectively.

template <int dim>
void test_cell(const Triangulation<dim>& tria, const DoFHandler<dim>& dof)
{
  LogStream::Prefix lp("cell");
  deallog.push("Regular");
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin(); cell != tria.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(); cell != dof.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;

  deallog << "begin\t\t"; dof.begin().print(deallog); deallog << '\t'; tria.begin().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last() .print(deallog); deallog << '\t'; tria.last() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; dof.end()  .print(deallog); deallog << '\t'; tria.end()  .print(deallog); deallog << std::endl;

  for (unsigned int l=0;l<tria.n_levels();++l)
    {
      deallog << "begin\t" << l << '\t'; dof.begin(l).print(deallog); deallog << '\t'; tria.begin(l).print(deallog); deallog << std::endl;
      deallog << "last \t" << l << '\t'; dof.last(l) .print(deallog); deallog << '\t'; tria.last(l) .print(deallog); deallog << std::endl;
      deallog << "end  \t" << l << '\t'; dof.end(l)  .print(deallog); deallog << '\t'; tria.end(l)  .print(deallog); deallog << std::endl;
    }

  deallog.pop();
  deallog << std::endl;
  deallog.push("Active ");
  for (typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(); cell != tria.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active(); cell != dof.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;

  
  deallog << "begin\t\t"; dof.begin_active().print(deallog); deallog << '\t'; tria.begin_active().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_active() .print(deallog); deallog << '\t'; tria.last_active() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  for (unsigned int l=0;l<tria.n_levels();++l)
    {
      deallog << "begin\t" << l << '\t'; dof.begin_active(l).print(deallog); deallog << '\t'; tria.begin_active(l).print(deallog); deallog << std::endl;
      deallog << "last \t" << l << '\t'; dof.last_active(l) .print(deallog); deallog << '\t'; tria.last_active(l) .print(deallog); deallog << std::endl;
      deallog << "end  \t" << l << '\t'; dof.end_active(l)  .print(deallog); deallog << '\t'; tria.end_active(l)  .print(deallog); deallog << std::endl;
    }
  deallog.pop();
  deallog << std::endl;
  deallog.push("Raw ");
  for (typename Triangulation<dim>::raw_cell_iterator cell = tria.begin_raw(); cell != tria.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::raw_cell_iterator cell = dof.begin_raw(); cell != dof.end();++cell)
    {
      deallog << '\t'; cell.print(deallog);
    }
  deallog << std::endl;

  
  deallog << "begin\t\t"; dof.begin_raw().print(deallog); deallog << '\t'; tria.begin_raw().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_raw() .print(deallog); deallog << '\t'; tria.last_raw() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  for (unsigned int l=0;l<tria.n_levels();++l)
    {
      deallog << "begin\t" << l << '\t'; dof.begin_raw(l).print(deallog); deallog << '\t'; tria.begin_raw(l).print(deallog); deallog << std::endl;
      deallog << "last \t" << l << '\t'; dof.last_raw(l) .print(deallog); deallog << '\t'; tria.last_raw(l) .print(deallog); deallog << std::endl;
      deallog << "end  \t" << l << '\t'; dof.end_raw(l)  .print(deallog); deallog << '\t'; tria.end_raw(l)  .print(deallog); deallog << std::endl;
    }
  deallog.pop();
}


template <int dim>
void test_face(const Triangulation<dim>& tria, const DoFHandler<dim>& dof)
{
  LogStream::Prefix lp("face");
  deallog.push("Regular");
  for (typename Triangulation<dim>::face_iterator face = tria.begin_face(); face != tria.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::face_iterator face = dof.begin_face(); face != dof.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;

  deallog << "begin\t\t"; dof.begin_face().print(deallog); deallog << '\t'; tria.begin_face().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_face() .print(deallog); deallog << '\t'; tria.last_face() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; dof.end_face()  .print(deallog); deallog << '\t'; tria.end_face()  .print(deallog); deallog << std::endl;

  deallog.pop();
  deallog << std::endl;
  deallog.push("Active ");
    for (typename Triangulation<dim>::active_face_iterator face = tria.begin_active_face(); face != tria.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::active_face_iterator face = dof.begin_active_face(); face != dof.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;

  deallog << "begin\t\t"; dof.begin_active_face().print(deallog); deallog << '\t'; tria.begin_active_face().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_active_face() .print(deallog); deallog << '\t'; tria.last_active_face() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  deallog.pop();
  deallog << std::endl;
  deallog.push("Raw ");
  for (typename Triangulation<dim>::raw_face_iterator face = tria.begin_raw_face(); face != tria.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::raw_face_iterator face = dof.begin_raw_face(); face != dof.end_face();++face)
    {
      deallog << '\t'; face.print(deallog);
    }
  deallog << std::endl;


  deallog << "begin\t\t"; dof.begin_raw_face().print(deallog); deallog << '\t'; tria.begin_raw_face().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_raw_face() .print(deallog); deallog << '\t'; tria.last_raw_face() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;
  
  deallog.pop();
}


template <int dim>
void test_line(const Triangulation<dim>& tria, const DoFHandler<dim>& dof)
{
  LogStream::Prefix lp("line");
  for (typename Triangulation<dim>::line_iterator line = tria.begin_line(); line != tria.end_line();++line)
    {
      deallog << '\t'; line.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::line_iterator line = dof.begin_line(); line != dof.end_line();++line)
    {
      deallog << '\t'; line.print(deallog);
    }
  deallog << std::endl;

  deallog.push("Regular");

  deallog << "begin\t\t"; dof.begin_line().print(deallog); deallog << '\t'; tria.begin_line().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_line() .print(deallog); deallog << '\t'; tria.last_line() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; dof.end_line()  .print(deallog); deallog << '\t'; tria.end_line()  .print(deallog); deallog << std::endl;

  if (dim == 1)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_line(l).print(deallog); deallog << '\t'; tria.begin_line(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_line(l) .print(deallog); deallog << '\t'; tria.last_line(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_line(l)  .print(deallog); deallog << '\t'; tria.end_line(l)  .print(deallog); deallog << std::endl;
      }

  deallog.pop();
  deallog << std::endl;
  deallog.push("Active ");
  
  deallog << "begin\t\t"; dof.begin_active_line().print(deallog); deallog << '\t'; tria.begin_active_line().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_active_line() .print(deallog); deallog << '\t'; tria.last_active_line() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  if (dim == 1)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_active_line(l).print(deallog); deallog << '\t'; tria.begin_active_line(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_active_line(l) .print(deallog); deallog << '\t'; tria.last_active_line(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_active_line(l)  .print(deallog); deallog << '\t'; tria.end_active_line(l)  .print(deallog); deallog << std::endl;
    }
  deallog.pop();
  deallog << std::endl;
  deallog.push("Raw ");
  
  deallog << "begin\t\t"; dof.begin_raw_line().print(deallog); deallog << '\t'; tria.begin_raw_line().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_raw_line() .print(deallog); deallog << '\t'; tria.last_raw_line() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;
  
  if (dim == 1)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_raw_line(l).print(deallog); deallog << '\t'; tria.begin_raw_line(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_raw_line(l) .print(deallog); deallog << '\t'; tria.last_raw_line(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_raw_line(l)  .print(deallog); deallog << '\t'; tria.end_raw_line(l)  .print(deallog); deallog << std::endl;
      }
  deallog.pop();
}


template <int dim>
void test_quad(const Triangulation<dim>& tria, const DoFHandler<dim>& dof)
{
  LogStream::Prefix lp("quad");
  for (typename Triangulation<dim>::quad_iterator quad = tria.begin_quad(); quad != tria.end_quad();++quad)
    {
      deallog << '\t';
      quad.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::quad_iterator quad = dof.begin_quad(); quad != dof.end_quad();++quad)
    {
      deallog << '\t';
      quad.print(deallog);
    }
  deallog << std::endl;

  deallog.push("Regular");

  deallog << "begin\t\t"; dof.begin_quad().print(deallog); deallog << '\t'; tria.begin_quad().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_quad() .print(deallog); deallog << '\t'; tria.last_quad() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; dof.end_quad()  .print(deallog); deallog << '\t'; tria.end_quad()  .print(deallog); deallog << std::endl;

  if (dim == 2)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_quad(l).print(deallog); deallog << '\t'; tria.begin_quad(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_quad(l) .print(deallog); deallog << '\t'; tria.last_quad(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_quad(l)  .print(deallog); deallog << '\t'; tria.end_quad(l)  .print(deallog); deallog << std::endl;
      }

  deallog.pop();
  deallog << std::endl;
  deallog.push("Active ");
  
  deallog << "begin\t\t"; dof.begin_active_quad().print(deallog); deallog << '\t'; tria.begin_active_quad().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_active_quad() .print(deallog); deallog << '\t'; tria.last_active_quad() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  if (dim == 2)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_active_quad(l).print(deallog); deallog << '\t'; tria.begin_active_quad(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_active_quad(l) .print(deallog); deallog << '\t'; tria.last_active_quad(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_active_quad(l)  .print(deallog); deallog << '\t'; tria.end_active_quad(l)  .print(deallog); deallog << std::endl;
    }
  deallog.pop();
  deallog << std::endl;
  deallog.push("Raw ");
  
  deallog << "begin\t\t"; dof.begin_raw_quad().print(deallog); deallog << '\t'; tria.begin_raw_quad().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_raw_quad() .print(deallog); deallog << '\t'; tria.last_raw_quad() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;
  
  if (dim == 2)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_raw_quad(l).print(deallog); deallog << '\t'; tria.begin_raw_quad(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_raw_quad(l) .print(deallog); deallog << '\t'; tria.last_raw_quad(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_raw_quad(l)  .print(deallog); deallog << '\t'; tria.end_raw_quad(l)  .print(deallog); deallog << std::endl;
      }
  deallog.pop();
}

template <int dim>
void test_hex(const Triangulation<dim>& tria, const DoFHandler<dim>& dof)
{
  LogStream::Prefix lp("hex");
  for (typename Triangulation<dim>::hex_iterator hex = tria.begin_hex(); hex != tria.end_hex();++hex)
    {
      deallog << '\t';
      hex.print(deallog);
    }
  deallog << std::endl;
  
  for (typename DoFHandler<dim>::hex_iterator hex = dof.begin_hex(); hex != dof.end_hex();++hex)
    {
      deallog << '\t';
      hex.print(deallog);
    }
  deallog << std::endl;

  deallog.push("Regular");

  deallog << "begin\t\t"; dof.begin_hex().print(deallog); deallog << '\t'; tria.begin_hex().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_hex() .print(deallog); deallog << '\t'; tria.last_hex() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; dof.end_hex()  .print(deallog); deallog << '\t'; tria.end_hex()  .print(deallog); deallog << std::endl;

  if (dim == 3)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_hex(l).print(deallog); deallog << '\t'; tria.begin_hex(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_hex(l) .print(deallog); deallog << '\t'; tria.last_hex(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_hex(l)  .print(deallog); deallog << '\t'; tria.end_hex(l)  .print(deallog); deallog << std::endl;
      }

  deallog.pop();
  deallog << std::endl;
  deallog.push("Active ");
  
  deallog << "begin\t\t"; dof.begin_active_hex().print(deallog); deallog << '\t'; tria.begin_active_hex().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_active_hex() .print(deallog); deallog << '\t'; tria.last_active_hex() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;

  if (dim == 3)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_active_hex(l).print(deallog); deallog << '\t'; tria.begin_active_hex(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_active_hex(l) .print(deallog); deallog << '\t'; tria.last_active_hex(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_active_hex(l)  .print(deallog); deallog << '\t'; tria.end_active_hex(l)  .print(deallog); deallog << std::endl;
    }
  deallog.pop();
  deallog << std::endl;
  deallog.push("Raw ");
  
  deallog << "begin\t\t"; dof.begin_raw_hex().print(deallog); deallog << '\t'; tria.begin_raw_hex().print(deallog); deallog << std::endl;
  deallog << "last \t\t"; dof.last_raw_hex() .print(deallog); deallog << '\t'; tria.last_raw_hex() .print(deallog); deallog << std::endl;
  deallog << "end  \t\t"; deallog << std::endl;
  
  if (dim == 3)
    for (unsigned int l=0;l<tria.n_levels();++l)
      {
	deallog << "begin\t" << l << '\t'; dof.begin_raw_hex(l).print(deallog); deallog << '\t'; tria.begin_raw_hex(l).print(deallog); deallog << std::endl;
	deallog << "last \t" << l << '\t'; dof.last_raw_hex(l) .print(deallog); deallog << '\t'; tria.last_raw_hex(l) .print(deallog); deallog << std::endl;
	deallog << "end  \t" << l << '\t'; dof.end_raw_hex(l)  .print(deallog); deallog << '\t'; tria.end_raw_hex(l)  .print(deallog); deallog << std::endl;
      }
  deallog.pop();
}

template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof (tria);
  dof.distribute_dofs (fe);
  
  test_cell(tria, dof);
  if (dim>1) test_face(tria, dof);
  if (dim >= 1) test_line(tria, dof);
  if (dim >= 2) test_quad(tria, dof);
  if (dim >= 3) test_hex(tria, dof);
}



int main ()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
//deallog.depth_console(0);

  deallog.push("1D");
  test<1> ();
  deallog.pop();
  deallog.push("2D");
  test<2> ();
  deallog.pop();
  deallog.push("3D");
  test<3> ();
  deallog.pop();

  return 0;
}
