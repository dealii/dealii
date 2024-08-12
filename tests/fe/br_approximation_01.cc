// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// projects the shape functions from the Bernardi-Raugel elements
// into different polynomial spaces and computes the error of the
// approximation. The error should be close to machine precision
// for polynomials of degree 2 or higher

// this is almost a verbatim copy of rt_approximation_01.cc and
// may needs additional modifications to make it more robust

#include "../tests.h"

#define PRECISION 8


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bernardi_raugel.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>



template <int dim>
class TestMap1 : public Function<dim>
{
public:
  TestMap1(const unsigned int n_components)
    : Function<dim>(n_components)
  {}

  virtual ~TestMap1()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  void
  vector_value(const Point<dim> &p, Vector<double> &return_value) const;
};



template <int dim>
double
TestMap1<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  if (component == 0)
    return (p(0) * p(0));
  if (component == 1)
    return (p(1) * p(1));

  return (0);
}


template <int dim>
void
TestMap1<dim>::vector_value(const Point<dim> &p,
                            Vector<double>   &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));

  for (unsigned int iCount = 0; iCount < this->n_components; ++iCount)
    return_value(iCount) = value(p, iCount);
}



///-----------------------------------------------------------------------

template <int dim>
class TestDef1 : public Function<dim>
{
private:
  const double phi;

public:
  TestDef1(const unsigned int n_components, const double ph)
    : Function<dim>(n_components)
    , phi(ph)
  {}

  virtual ~TestDef1()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  void
  vector_value(const Point<dim> &p, Vector<double> &return_value) const;
};



template <int dim>
double
TestDef1<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  Point<2> center;
  center[0]    = 0.5;
  center[1]    = 0.5;
  double rad   = p.distance(center),
         phi_p = atan2(p[0] - center[0], p[1] - center[1]);

  if (component == 0)
    return rad * (sin(phi + phi_p) - sin(phi_p));
  else
    return rad * (cos(phi + phi_p) - cos(phi_p));
}


template <int dim>
void
TestDef1<dim>::vector_value(const Point<dim> &p,
                            Vector<double>   &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; ++iCount)
    return_value(iCount) = value(p, iCount);
}


///-----------------------------------------------------------------------


template <int dim>
class TestDef2 : public Function<dim>
{
private:
  const double scale;

public:
  TestDef2(const unsigned int n_components, const double sc)
    : Function<dim>(n_components)
    , scale(sc)
  {}

  virtual ~TestDef2()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  void
  vector_value(const Point<dim> &p, Vector<double> &return_value) const;
};


template <int dim>
double
TestDef2<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  double x = p[0], y = p[1];

  if (component == 0)
    return scale * x;
  else
    return scale * y;
}


template <int dim>
void
TestDef2<dim>::vector_value(const Point<dim> &p,
                            Vector<double>   &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; ++iCount)
    return_value(iCount) = value(p, iCount);
}


///-----------------------------------------------------------------------
// testDef3 implements parallelograms ...


template <int dim>
class TestDef3 : public Function<dim>
{
private:
  const double scale;

public:
  TestDef3(const unsigned int n_components, const double sc)
    : Function<dim>(n_components)
    , scale(sc)
  {}

  virtual ~TestDef3()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  void
  vector_value(const Point<dim> &p, Vector<double> &return_value) const;
};


template <int dim>
double
TestDef3<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  double y = p[1];

  if (component == 0)
    return scale * y;
  else
    return 0;
}


template <int dim>
void
TestDef3<dim>::vector_value(const Point<dim> &p,
                            Vector<double>   &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; ++iCount)
    return_value(iCount) = value(p, iCount);
}

template <int dim>
class TestPoly : public Function<dim>
{
private:
  std::vector<Polynomials::Polynomial<double>> polys;

public:
  TestPoly(unsigned int deg)
    : Function<dim>(2)
  {
    std::vector<double> coeff(deg, 0.0);
    for (unsigned int p = 0; p < 4; ++p)
      {
        for (unsigned int i = 0; i < deg; ++i)
          {
            double c = ((double)Testing::rand()) / ((double)RAND_MAX + 1);
            coeff[i] = c;
          }
        polys.push_back(Polynomials::Polynomial<double>(coeff));
      }
  }

  virtual ~TestPoly()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  void
  vector_value(const Point<dim> &p, Vector<double> &return_value) const;
};


template <int dim>
double
TestPoly<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  double x = p[0], y = p[1];

  if (component == 0)
    return polys[0].value(x) + polys[1].value(y);
  else
    return polys[2].value(x) + polys[3].value(y);
}


template <int dim>
void
TestPoly<dim>::vector_value(const Point<dim> &p,
                            Vector<double>   &return_value) const
{
  Assert(return_value.size() == this->n_components,
         ExcDimensionMismatch(return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; ++iCount)
    return_value(iCount) = value(p, iCount);
}



/*
 * Check the value of the derivative field.
 */

double
TestProjection(Mapping<2> &mapping, DoFHandler<2> *dof_handler)
{
  Vector<double> solution;
  solution.reinit(dof_handler->n_dofs());

  // Create and Project test function to H_div space and evaluate the error
  // If the error is in the range of machine precision, this polynomial
  // degree can obviously be represented by projected field.

  for (unsigned int deg = 1; deg < 4; ++deg)
    {
      TestPoly<2> pol(deg);

      // Project solution onto BR FE field
      AffineConstraints<double> hn_constraints;
      hn_constraints.clear();
      DoFTools::make_hanging_node_constraints(*dof_handler, hn_constraints);
      hn_constraints.close();
      VectorTools::project(
        mapping, *dof_handler, hn_constraints, QGauss<2>(6), pol, solution);

      // Now evaluate error ...
      // Use a high order quadrature.
      QGauss<2>   quad(6);
      FEValues<2> fe_values(mapping,
                            dof_handler->get_fe(),
                            quad,
                            UpdateFlags(update_values |
                                        update_quadrature_points |
                                        update_gradients | update_JxW_values |
                                        update_contravariant_transformation));

      const unsigned int n_q_points   = quad.size();
      const unsigned int n_components = dof_handler->get_fe().n_components();

      // Cell iterators
      DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active(),
                                          endc = dof_handler->end();

      double err_u = 0, err_v = 0;

      for (; cell != endc; ++cell)
        {
          fe_values.reinit(cell);

          // Get values from solution vector (For Trap.Rule)
          std::vector<Vector<double>> this_value(n_q_points,
                                                 Vector<double>(n_components));
          fe_values.get_function_values(solution, this_value);

          for (const auto q_point : fe_values.quadrature_point_indices())
            {
              double   u = this_value[q_point](0), v = this_value[q_point](1);
              Point<2> p = fe_values.quadrature_point(q_point);

              double u_ex = pol.value(p, 0), v_ex = pol.value(p, 1);

              double JxW = fe_values.JxW(q_point);

              err_u += (u - u_ex) * (u - u_ex) * JxW;
              err_v += (v - v_ex) * (v - v_ex) * JxW;
            }
        }

      deallog << dof_handler->get_fe().get_name() << ", testfun(" << deg
              << "), error_u=" << err_u << ", error_v=" << err_v << std::endl;
    }


  // Test the core functionality
  DataOut<2> *data_out = new DataOut<2>;
  data_out->attach_dof_handler(*dof_handler);
  data_out->add_data_vector(solution, "solution");
  data_out->build_patches(mapping, 4);

  data_out->write_gnuplot(deallog.get_file_stream());

  delete data_out;

  return (0.0);
}


int
main()
{
  initlog();
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.get_file_stream() << std::setprecision(PRECISION);
  deallog.get_file_stream() << std::fixed;

  Triangulation<2> tria_test;
  DoFHandler<2>   *dof_handler, *dof_handler_def;
  Point<2>         p1(0, 0), p2(1, 1);

  GridGenerator::hyper_rectangle(tria_test, p1, p2);
  //  tria_test.refine_global (1);
  //   GridTools::distort_random (0.4, tria_test);

  // Create a DoFHandler for the BR space
  FE_BernardiRaugel<2> fe(1);
  dof_handler = new DoFHandler<2>(tria_test);
  dof_handler->distribute_dofs(fe);

  QGauss<2> quad_temp(6);
  deallog << "DoFs per quad: " << fe.dofs_per_quad
          << ", dofs per line: " << fe.dofs_per_line
          << ", dofs per vertex: " << fe.dofs_per_vertex << std::endl;
  deallog << "n_q_points=" << quad_temp.size() << std::endl;

  // Create a deformation object for
  // the mapping
  FESystem<2> fe_def(FE_Q<2>(1), 2);
  dof_handler_def = new DoFHandler<2>(tria_test);
  dof_handler_def->distribute_dofs(fe_def);

  // Alloc some DoFs
  Vector<double> solution, solution_q, deformation;
  solution.reinit(dof_handler->n_dofs());
  solution_q.reinit(dof_handler->n_dofs());
  deformation.reinit(dof_handler_def->n_dofs());

  AffineConstraints<double> hn_constraints_def;
  hn_constraints_def.clear();
  DoFTools::make_hanging_node_constraints(*dof_handler_def, hn_constraints_def);
  hn_constraints_def.close();

  {
    MappingQ1Eulerian<2> mapping_euler(*dof_handler_def, deformation);

    // Try rotating the elements
    for (double rotate = 0; rotate < 2 * numbers::PI;
         rotate += 0.25 * numbers::PI)
      {
        // Rotate element
        VectorTools::project(*dof_handler_def,
                             hn_constraints_def,
                             QGauss<2>(6),
                             TestDef1<2>(2, rotate),
                             deformation);
        deallog << "phi = " << rotate << std::endl;
        TestProjection(mapping_euler, dof_handler);
      }

    // Try resizing the elements
    for (double scale = -0.75; scale < 4.0; scale += 0.25)
      {
        VectorTools::project(*dof_handler_def,
                             hn_constraints_def,
                             QGauss<2>(6),
                             TestDef2<2>(2, scale),
                             deformation);
        deallog << "scale = " << scale << std::endl;
        TestProjection(mapping_euler, dof_handler);
      }

    // Try paralellogramming the elements
    for (double scale = -1.0; scale < 1.0; scale += 0.25)
      {
        VectorTools::project(*dof_handler_def,
                             hn_constraints_def,
                             QGauss<2>(6),
                             TestDef3<2>(2, scale),
                             deformation);
        deallog << "scale = " << scale << std::endl;
        TestProjection(mapping_euler, dof_handler);
      }

    // Try arbitrary deformation ...
    for (unsigned int i = 0; i < deformation.size(); ++i)
      {
        double c       = ((double)Testing::rand()) / ((double)RAND_MAX + 1);
        deformation(i) = 0.35 * (c - 0.5);
      }

    deallog << "Arbitrary\n" << std::endl;
    TestProjection(mapping_euler, dof_handler);
  }

  delete (dof_handler);
  delete (dof_handler_def);

  return (0);
}
