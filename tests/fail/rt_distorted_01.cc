// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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


/*
 * Snippet to demonstrate some properties of the deal.II implementation of
 * the RT spaces.
 */



#include "../tests.h"
#include <deal.II/base/logstream.h>

#define PRECISION 2

#include <fstream>

std::ofstream logfile ("output");

#include <fstream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <fstream>


template <int dim>
class TestMap1 : public Function<dim>
{
public:
  TestMap1 (const unsigned int n_components) :
    Function<dim> (n_components)
  {}

  virtual ~TestMap1 () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  void vector_value (const Point<dim> &p,
                     Vector<double>   &return_value) const;
};



template <int dim>
double
TestMap1<dim>::value (const Point<dim>   &/*p*/,
                      const unsigned int  /*component*/) const
{
  return (1);
}


template <int dim>
void TestMap1<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));

  // Parabolic inflow profile
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}

///-----------------------------------------------------------------------

template <int dim>
class TestDef1 : public Function<dim>
{
private:
  const double phi;

public:
  TestDef1 (const unsigned int n_components, const double ph) :
    Function<dim> (n_components),
    phi (ph)
  {}

  virtual ~TestDef1 () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  void vector_value (const Point<dim> &p,
                     Vector<double>   &return_value) const;
};



template <int dim>
double
TestDef1<dim>::value (const Point<dim>   &p,
                      const unsigned int  component) const
{
  Point<2> center;
  center(0) = 0.5;
  center(1) = 0.5;
  double rad = p.distance (center),
         phi_p = atan2 (p(0) - center(0), p(1) - center(1));

  if (component == 0)
    return rad * (sin (phi + phi_p) - sin (phi_p));
  else
    return rad * (cos (phi + phi_p) - cos (phi_p));
}


template <int dim>
void TestDef1<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}


///-----------------------------------------------------------------------


template <int dim>
class TestDef2 : public Function<dim>
{
private:
  const double scale;

public:
  TestDef2 (const unsigned int n_components, const double sc) :
    Function<dim> (n_components),
    scale (sc)
  {}

  virtual ~TestDef2 () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  void vector_value (const Point<dim> &p,
                     Vector<double>   &return_value) const;
};


template <int dim>
double
TestDef2<dim>::value (const Point<dim>   &p,
                      const unsigned int  component) const
{
  double x = p(0),
         y = p(1);

  if (component == 0)
    return scale * x;
  else
    return scale * y;
}


template <int dim>
void TestDef2<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}


///-----------------------------------------------------------------------
// testDef3 implements parallelograms ...


template <int dim>
class TestDef3 : public Function<dim>
{
private:
  const double scale;

public:
  TestDef3 (const unsigned int n_components, const double sc) :
    Function<dim> (n_components),
    scale (sc)
  {}

  virtual ~TestDef3 () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const;

  void vector_value (const Point<dim> &p,
                     Vector<double>   &return_value) const;
};


template <int dim>
double
TestDef3<dim>::value (const Point<dim>   &p,
                      const unsigned int  component) const
{
  double y = p(1);

  if (component == 0)
    return scale * y;
  else
    return 0;
}


template <int dim>
void TestDef3<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
          ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}



/*
 * Integrate the function value over the element.
 */

double EvaluateArea (Mapping<2> &mapping,
                     DoFHandler<2> *dof_handler,
                     Vector<double> &solution)
{
  // Use a high order quadrature.
  QGauss<2> quad (6);
  FEValues<2> fe_values (mapping, dof_handler->get_fe (), quad,
                         UpdateFlags(update_values    |
                                     update_q_points  |
                                     update_JxW_values));

  const unsigned int   n_q_points    = quad.size();
  const unsigned int   n_components   = dof_handler->get_fe().n_components();

  // Cell iterators
  DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active(),
                                      endc = dof_handler->end();
  double result_u = 0,
         result_v = 0;

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      // Get values from solution vector (For Trap.Rule)
      std::vector<Vector<double> > this_value
      (n_q_points, Vector<double>(n_components));
      fe_values.get_function_values (solution, this_value);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
          double JxW = fe_values.JxW (q_point);
          result_u += this_value[q_point](0) * JxW;
          result_v += this_value[q_point](1) * JxW;
        }
    }

  return (result_v);
}


int main (int /*argc*/, char **/*argv*/)
{
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria_test;
  DoFHandler<2> *dof_handler,
             *dof_handler_def;
  Point<2> p1 (0,0),
        p2 (1, 1);

  GridGenerator::hyper_rectangle (tria_test, p1, p2);
  tria_test.refine_global (1);

  // Uncommenting the following line, demonstrates the problem
  // of RT elements on distorted Quads!
  tria_test.distort_random (0.4);

  // Create a DoFHandler for the RT space
  FE_RaviartThomas<2> fe (2);
  dof_handler = new DoFHandler<2> (tria_test);
  dof_handler->distribute_dofs (fe);

  // Create an deformation object for the Eulerian mapping
  FESystem<2> fe_def (FE_Q<2>(1), 2);
  dof_handler_def = new DoFHandler<2> (tria_test);
  dof_handler_def->distribute_dofs (fe_def);

  // Alloc some DoFs
  Vector<double> solution,
         solution_q,
         deformation;
  solution.reinit (dof_handler->n_dofs ());
  solution_q.reinit (dof_handler_def->n_dofs ());
  deformation.reinit (dof_handler_def->n_dofs ());

  // Project solution onto RT FE field
  ConstraintMatrix     hn_constraints;
  hn_constraints.clear ();
  DoFTools::make_hanging_node_constraints (*dof_handler,
                                           hn_constraints);
  hn_constraints.close ();
  VectorTools::project (*dof_handler, hn_constraints,
                        QGauss<2> (6), TestMap1<2>(2),
                        solution);

  // Project reference solution onto RT FE field
  ConstraintMatrix     hn_constraints_def;
  hn_constraints_def.clear ();
  DoFTools::make_hanging_node_constraints (*dof_handler_def,
                                           hn_constraints_def);
  hn_constraints_def.close ();

  VectorTools::project (*dof_handler_def, hn_constraints_def,
                        QGauss<2> (6), TestMap1<2>(2),
                        solution_q);

  MappingQ1Eulerian<2> *mapping_euler
    = new MappingQ1Eulerian<2> (deformation, *dof_handler_def);

  char buf[1000];
  sprintf (buf,
           "FE_RT Area %e  FE_Q Area %e\n",
           EvaluateArea (*mapping_euler, dof_handler, solution),
           EvaluateArea (*mapping_euler, dof_handler_def, solution_q));
  deallog << buf;

  unsigned int test_out = 0;
  // Try rotating the elements
  for (double rotat = 0; rotat < 2 * M_PI; rotat += 0.25 * M_PI)
    {
      // Rotate element
      VectorTools::project (*dof_handler_def, hn_constraints_def,
                            QGauss<2> (6), TestDef1<2>(2, rotat),
                            deformation);

      // Project 1 function to element
      VectorTools::project (*mapping_euler, *dof_handler, hn_constraints,
                            QGauss<2> (6), TestMap1<2>(2),
                            solution);

      // Write output files
      DataOut<2> *data_out = new DataOut<2>;
      data_out->attach_dof_handler (*dof_handler);
      data_out->add_data_vector (solution, "solution");
      data_out->build_patches (*mapping_euler, 8);

      data_out->write_gnuplot (deallog.get_file_stream());
      test_out++;

      delete data_out;

      double area_rt = EvaluateArea (*mapping_euler, dof_handler, solution);
      double area_q = EvaluateArea (*mapping_euler, dof_handler_def, solution_q);

      char buf[100];
      sprintf (buf, "phi = %e FE_RT Area %e  FE_Q Area %e\n",
               rotat,
               area_rt, area_q);
      deallog << buf;
    }

  // Try resizing the elements
  for (double scale = -0.75; scale < 4.0; scale += 0.25)
    {
      VectorTools::project (*dof_handler_def, hn_constraints_def,
                            QGauss<2> (6), TestDef2<2>(2, scale),
                            deformation);

      // Project 1 function to element
      VectorTools::project (*mapping_euler, *dof_handler, hn_constraints,
                            QGauss<2> (6), TestMap1<2>(2),
                            solution);

      char buf[1000];
      sprintf (buf, "Scale = %e FE_RT Area %e  FE_Q Area %e\n",
               scale,
               EvaluateArea (*mapping_euler, dof_handler, solution),
               EvaluateArea (*mapping_euler, dof_handler_def, solution_q));
      deallog << buf;
    }


  // Try parallelograms
  for (double scale = -1.0; scale < 1.0; scale += 0.25)
    {
      VectorTools::project (*dof_handler_def, hn_constraints_def,
                            QGauss<2> (6), TestDef3<2>(2, scale),
                            deformation);

      // Project 1 function to element
      VectorTools::project (*mapping_euler, *dof_handler, hn_constraints,
                            QGauss<2> (6), TestMap1<2>(2),
                            solution);

      char buf[1000];
      sprintf (buf, "Parallelogram = %e FE_RT Area %e  FE_Q Area %e\n",
               scale,
               EvaluateArea (*mapping_euler, dof_handler, solution),
               EvaluateArea (*mapping_euler, dof_handler_def, solution_q));
      deallog << buf;
    }

  delete (mapping_euler);

  delete (dof_handler);
  delete (dof_handler_def);

  return (0);
}
