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



/* Author: Denis Davydov, University of Erlangen-Nuremberg, 2015 */
// make sure that the modified Kelly error estimator with
// strategy = face_diameter_over_twice_max_degree
// returns correct results.
// To that end consider 3 cases at domain [0,2]^dim and one at domain [0,1]^dim:
//
// 1)
// Neuman BC g(x) = c; // constant
// single element with p =3; u=0;
// \eta = h * A * c^2 / p;
//
//
// 2)
//
// ---------------
// |      |      |
// |      |      |
// |  p1  |  p2  |
// |      |      |
// |      |      |
// ---------------
//
// u (left)  = 0;
// u (right) = kx;
//
// \eta = h * A * k^2 / 2 max(p1,p2)
//
//
// 3)
// ----------------------
// |     |      |       |
// |  p3 | p3  f_1      |
// |------------|  p1   |
// |  p2 | p2  f_2      |
// |     |      |       |
// ----------------------
// solution is the same as above;
//
// \eta_1 = h * A * k^2 / 2 max(p3,p1)
// \eta_2 = h * A * k^2 / 2 max(p2,p1)
// \eta_3 = \eta_1 + \eta_2
//
// 4) just arbitrary mesh
// -------------------
// |        |        |
// |    1   |   1    |
// |        |        |
// |------------------
// |  3 | 2 |        |
// |--------|    1   |
// |  3 | 2 |        |
// ------------------
//
// and evaluate error of the interpolated to it function.



#include "../tests.h"

// dealii
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/error_estimator.h>

// c++
#include <fstream>
#include <iostream>

using namespace dealii;

template<int dim>
class MyFunction : public dealii::Function<dim>
{
public:
  MyFunction(const double k);

  virtual double value(const dealii::Point<dim> &point,
                       const unsigned int component = 0 ) const;

  double get_k() const;

private:
  const double k;
};

template<int dim>
MyFunction<dim>::MyFunction(const double k)
  :
  dealii::Function<dim>(1),
  k(k)
{

}

template<int dim>
double MyFunction<dim>::value(const dealii::Point<dim> &point,
                              const unsigned int ) const
{
  const double x = point[0]-1.0;

  if (x < 0)
    return 0.0;
  else
    return k * x;
}

template<int dim>
double MyFunction<dim>::get_k() const
{
  return k;
}

// neuman bc
template<int dim>
class NeumanBC : public dealii::Function<dim>
{
public:
  NeumanBC(const double c);

  virtual double value(const dealii::Point<dim> &point,
                       const unsigned int component = 0 ) const;

  double get_c() const;

private:
  const double c;
};

template<int dim>
NeumanBC<dim>::NeumanBC(const double c)
  :
  dealii::Function<dim>(1),
  c(c)
{
}

template<int dim>
double NeumanBC<dim>::value(const dealii::Point<dim> &point,
                            const unsigned int ) const
{
  return c;
}

template<int dim>
double NeumanBC<dim>::get_c() const
{
  return c;
}

// helper function to get diagonal and
// area of the squared element with lenght h
template<int dim>
void get_h_area(double &h, double &a, const double L);

template<>
void get_h_area<2>(double &h, double &a, const double L)
{
  h = L;
  a = L;
}

template<>
void get_h_area<3>(double &h, double &a, const double L)
{
  h = std::sqrt(2.0)*L;
  a = L*L;
}

// helper function to get diagonal and area of the
// h-refined face.
template<int dim>
void get_h_area_sub(double &h, double &a, const double L);

template<>
void get_h_area_sub<2>(double &h, double &a, const double L)
{
  h = L/2;
  a = L/2;
}

template<>
void get_h_area_sub<3>(double &h, double &a, const double L)
{
  h = std::sqrt(2.0)*L/2;
  a = L*L/4.0;
}

// output for inspection
template<int dim>
void output(const std::string          name,
            const Triangulation<dim>  &triangulation,
            const hp::DoFHandler<dim> &dof_handler,
            const Vector<double>      &values,
            const Vector<float>       &error)
{
  dealii::Vector<double> fe_degrees(triangulation.n_active_cells());
  {
    typename dealii::hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      fe_degrees(index) = dof_handler.get_fe()[cell->active_fe_index()].degree;
  }

  dealii::DataOut<dim,dealii::hp::DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector(values,
                           std::string("function_interpolation"));
  data_out.add_data_vector(fe_degrees,
                           std::string("fe_degree"));
  data_out.add_data_vector(error,
                           std::string("error"));
  data_out.build_patches ();

  std::ofstream output (name.c_str());
  data_out.write_vtu (output);
}

// case 1)
template<int dim>
void test_neumann(const NeumanBC<dim> &func)
{
  deallog << "NeumanBC case:"<<std::endl;
  deallog << "--------------"<<std::endl;
  Triangulation<dim> triangulation;
  hp::DoFHandler<dim> dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_formula;
  hp::QCollection<dim-1> face_quadrature_formula;
  ConstraintMatrix constraints;

  const unsigned int p = 3;

  fe_collection.push_back(dealii::FE_Q<dim>(QIterated<1>(QTrapez<1>(),p)));
  quadrature_formula.push_back(dealii::QGauss<dim>(p+5));
  face_quadrature_formula.push_back(dealii::QGauss<dim-1>(p+5));

  const double L = 2.0;

  // set-up domain
  {
    GridGenerator::hyper_cube (triangulation,.0,L,/*colorize*/true);
  }

  dof_handler.distribute_dofs (fe_collection);

  // constraints
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  // interpolate some function
  dealii::Vector<double> values(dof_handler.n_dofs());
  values = 0.0;

  dealii::deallog << "dof values:"<<std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    dealii::deallog << " " << values[i];
  dealii::deallog << std::endl;

  // call Kelly
  typename dealii::FunctionMap<dim>::type function_map;
  function_map[0] = &func;

  dealii::Vector<float> error(dof_handler.n_dofs());
  dealii::KellyErrorEstimator<dim>::estimate(dof_handler,
                                             face_quadrature_formula,
                                             function_map,
                                             values,
                                             error,
                                             dealii::ComponentMask(),
                                             0,
                                             dealii::numbers::invalid_unsigned_int,
                                             dealii::numbers::invalid_subdomain_id,
                                             dealii::numbers::invalid_material_id,
                                             dealii::KellyErrorEstimator<dim>::face_diameter_over_twice_max_degree);

  dealii::deallog <<"error:"<<std::endl;
  for (unsigned int i = 0; i < error.size(); i++)
    dealii::deallog << " " << error[i];
  dealii::deallog << std::endl;

  //output("neuman.vtu",
  //       triangulation,
  //       dof_handler,
  //       values,
  //       error);

  double h,A;
  get_h_area<dim>(h,A,L);
  const double expected_value_squared = h*A*std::pow(func.get_c(),2)/p;
  dealii::deallog << "expected:"<< std::endl <<" "<< std::sqrt(expected_value_squared) << std::endl;

  AssertThrow (std::fabs(std::sqrt(expected_value_squared) - error[0] ) < 1e-5, dealii::ExcInternalError());

  dof_handler.clear();
}

// case 2)
template<int dim>
void test_regular(const MyFunction<dim> &func)
{
  deallog << std::endl;
  deallog << "Regular face:"<<std::endl;
  deallog << "-------------"<<std::endl;
  Triangulation<dim> triangulation;
  hp::DoFHandler<dim> dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_formula;
  hp::QCollection<dim-1> face_quadrature_formula;
  ConstraintMatrix constraints;

  const unsigned int p1 = 1;
  const unsigned int p2 = 2;
  std::vector<unsigned int> p_degree;
  p_degree.push_back(p1);
  p_degree.push_back(p2);

  for (unsigned int i=0; i<p_degree.size(); i++)
    {
      const unsigned int &p = p_degree[i];
      fe_collection.push_back(dealii::FE_Q<dim>(QIterated<1>(QTrapez<1>(),p)));
      quadrature_formula.push_back(dealii::QGauss<dim>(p+5));
      face_quadrature_formula.push_back(dealii::QGauss<dim-1>(p+5));
    }

  const double L = 2.0;

  // set-up domain
  {
    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = 2;
    dealii::Point<dim> p1;
    dealii::Point<dim> p2;
    for (unsigned int d = 0; d < dim; d++)
      {
        p1[d] = 0.0;
        p2[d] = L;
      }
    GridGenerator::subdivided_hyper_rectangle   (triangulation,
                                                 repetitions,
                                                 p1,
                                                 p2,
                                                 /*colorize*/ false);

    typename dealii::hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell != endc; cell++)
      if (cell->center()[0] > 1.0)
        {
          cell->set_active_fe_index(1);
          break;
        }

  }

  dof_handler.distribute_dofs (fe_collection);

  // constraints
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  // interpolate some function
  dealii::Vector<double> values(dof_handler.n_dofs());
  dealii::VectorTools::interpolate (dof_handler,
                                    func,
                                    values);

  dealii::deallog << "dof values:"<<std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    dealii::deallog << " " << values[i];
  dealii::deallog << std::endl;

  // call Kelly
  dealii::Vector<float> error(dof_handler.n_dofs());
  dealii::KellyErrorEstimator<dim>::estimate(dof_handler,
                                             face_quadrature_formula,
                                             typename dealii::FunctionMap<dim>::type(),
                                             values,
                                             error,
                                             dealii::ComponentMask(),
                                             0,
                                             dealii::numbers::invalid_unsigned_int,
                                             dealii::numbers::invalid_subdomain_id,
                                             dealii::numbers::invalid_material_id,
                                             dealii::KellyErrorEstimator<dim>::face_diameter_over_twice_max_degree);

  dealii::deallog <<"error:"<<std::endl;
  for (unsigned int i = 0; i < error.size(); i++)
    dealii::deallog << " " << error[i];
  dealii::deallog << std::endl;

  //output("regular.vtu",
  //       triangulation,
  //       dof_handler,
  //       values,
  //       error);

  double h,A;
  get_h_area<dim>(h,A,L);
  const double expected_value_squared = h*A*std::pow(func.get_k(),2)/2.0/std::max(p1,p2);
  dealii::deallog << "expected:"<< std::endl <<" "<< std::sqrt(expected_value_squared) << std::endl;
  for (unsigned int i = 0; i < error.size(); i++)
    AssertThrow (std::fabs(std::sqrt(expected_value_squared) - error[i] ) < 1e-6, dealii::ExcInternalError());

  dof_handler.clear();
}

// case 3)
template<int dim>
void test_irregular(const MyFunction<dim> &func)
{
  deallog << std::endl;
  deallog << "Irregular face:"<<std::endl;
  deallog << "---------------"<<std::endl;
  Triangulation<dim> triangulation;
  hp::DoFHandler<dim> dof_handler(triangulation);
  hp::FECollection<dim> fe_collection;
  hp::QCollection<dim> quadrature_formula;
  hp::QCollection<dim-1> face_quadrature_formula;
  ConstraintMatrix constraints;

  const unsigned int p1 = 1;
  const unsigned int p2 = 2;
  const unsigned int p3 = 3;
  std::vector<unsigned int> p_degree;
  p_degree.push_back(p1);
  p_degree.push_back(p2);
  p_degree.push_back(p3);

  for (unsigned int i=0; i<p_degree.size(); i++)
    {
      const unsigned int &p = p_degree[i];
      fe_collection.push_back(dealii::FE_Q<dim>(QIterated<1>(QTrapez<1>(),p)));
      quadrature_formula.push_back(dealii::QGauss<dim>(p+5));
      face_quadrature_formula.push_back(dealii::QGauss<dim-1>(p+5));
    }

  const double L = 2.0;

  // set-up domain
  {
    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = 2;
    dealii::Point<dim> p1;
    dealii::Point<dim> p2;
    for (unsigned int d = 0; d < dim; d++)
      {
        p1[d] = 0.0;
        p2[d] = L;
      }
    GridGenerator::subdivided_hyper_rectangle   (triangulation,
                                                 repetitions,
                                                 p1,
                                                 p2,
                                                 /*colorize*/ false);
    // refine left side
    {
      typename dealii::hp::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active();
      cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }

    typename dealii::hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell != endc; cell++)
      if (cell->center()[0] > 1.0) // right
        {
          cell->set_active_fe_index(0);
        }
      else if (cell->center()[1] > 1.0) // top
        {
          cell->set_active_fe_index(2);
        }
      else
        {
          cell->set_active_fe_index(1);
        }
  }

  dof_handler.distribute_dofs (fe_collection);

  // constraints
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  // interpolate some function
  dealii::Vector<double> values(dof_handler.n_dofs());
  dealii::VectorTools::interpolate (dof_handler,
                                    func,
                                    values);

  dealii::deallog << "dof values:"<<std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    dealii::deallog << " " << values[i];
  dealii::deallog << std::endl;

  // call Kelly
  dealii::Vector<float> error(dof_handler.n_dofs());
  dealii::KellyErrorEstimator<dim>::estimate(dof_handler,
                                             face_quadrature_formula,
                                             typename dealii::FunctionMap<dim>::type(),
                                             values,
                                             error,
                                             dealii::ComponentMask(),
                                             0,
                                             dealii::numbers::invalid_unsigned_int,
                                             dealii::numbers::invalid_subdomain_id,
                                             dealii::numbers::invalid_material_id,
                                             dealii::KellyErrorEstimator<dim>::face_diameter_over_twice_max_degree);

  dealii::deallog <<"error:"<<std::endl;
  for (unsigned int i = 0; i < error.size(); i++)
    dealii::deallog << " " << error[i];
  dealii::deallog << std::endl;

  //output("irregular.vtu",
  //       triangulation,
  //       dof_handler,
  //       values,
  //       error);

  double h,A;
  get_h_area_sub<dim>(h,A,L);
  //
  const double expected_squared_1 = h*A*std::pow(func.get_k(),2)/2.0/std::max(p3,p1);
  const double expected_squared_2 = h*A*std::pow(func.get_k(),2)/2.0/std::max(p2,p1);
  const double expected_squared_3 = (dim==2) ?
                                    expected_squared_1 +   expected_squared_2:
                                    2*expected_squared_1 + 2*expected_squared_2;

  std::vector<double> expected_error(error.size(),0.0);

  expected_error[0] = std::sqrt(expected_squared_3);
  expected_error[2] = std::sqrt(expected_squared_2);
  expected_error[4] = std::sqrt(expected_squared_1);

  if (dim ==3)
    {
      expected_error[6] = expected_error[2];
      expected_error[8] = expected_error[4];
    }

  dealii::deallog << "expected:"<< std::endl;
  for (unsigned int i = 0; i < expected_error.size(); i++)
    deallog<<" " << expected_error[i];
  deallog <<std::endl;

  for (unsigned int i = 0; i < expected_error.size(); i++)
    AssertThrow (std::fabs(expected_error[i] - error[i] ) < 1e-6, dealii::ExcInternalError());

  dof_handler.clear();
}

template<int dim>
class MySecondFunction : public dealii::Function<dim>
{
public:
  MySecondFunction();

  virtual double value(const dealii::Point<dim> &point,
                       const unsigned int component = 0 ) const;
};

template<int dim>
MySecondFunction<dim>::MySecondFunction()
  :
  dealii::Function<dim>(1)
{

}

template<int dim>
double MySecondFunction<dim>::value(const dealii::Point<dim> &point,
                                    const unsigned int ) const
{
  double f = 0.0;
  const double &x = point[0];
  Assert (dim>1, dealii::ExcNotImplemented());
  const double &y = point[1];

  return (1.-x)*(1.-y)*(1.-y)+std::pow(1.0-y,4)*std::exp(-x);
}

template<int dim>
void test(const MySecondFunction<dim> &func)
{
  deallog << std::endl;
  deallog << "More complicated mesh:"<<std::endl;
  deallog << "----------------------"<<std::endl;

  dealii::Triangulation<dim> triangulation;
  dealii::hp::DoFHandler<dim> dof_handler(triangulation);
  dealii::hp::FECollection<dim> fe_collection;
  dealii::hp::QCollection<dim> quadrature_formula;
  dealii::hp::QCollection<dim-1> face_quadrature_formula;
  dealii::ConstraintMatrix constraints;
  for (unsigned int p = 1; p <=3; p++)
    {
      fe_collection.push_back(dealii::FE_Q<dim>(QIterated<1>(QTrapez<1>(),p)));
      quadrature_formula.push_back(dealii::QGauss<dim>(p+1));
      face_quadrature_formula.push_back(dealii::QGauss<dim-1>(p+1));
    }
  dealii::GridGenerator::hyper_cube (triangulation,0.0,1.0); // reference cell

  // refine
  {
    triangulation.refine_global(1);

    // otherwise the next set_active_fe_index
    // will not carry to the child cells.
    dof_handler.distribute_dofs (fe_collection);

    typename dealii::hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell != endc; cell++)
      {
        bool in_top_left = true;
        for (unsigned int d=0; d< dim; d++)
          in_top_left = in_top_left && (cell->center()[d] < 0.5);

        if (in_top_left)
          {
            cell->set_active_fe_index(1);
            cell->set_refine_flag();
            break;
          }
      }

    triangulation.prepare_coarsening_and_refinement();

    triangulation.execute_coarsening_and_refinement();

    cell = dof_handler.begin_active();

    for (; cell != endc; cell++)
      {
        if (cell->center()[0] < 0.25)
          {
            cell->set_active_fe_index(2);
          }
      }
  }

  dof_handler.distribute_dofs (fe_collection);

  // constraints
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();

  // interpolate some function
  dealii::Vector<double> values(dof_handler.n_dofs());

  dealii::VectorTools::interpolate (dof_handler,
                                    func,
                                    values);

  dealii::deallog << "dof values:"<<std::endl;
  for (unsigned int i = 0; i < values.size(); i++)
    dealii::deallog << " " << values[i];
  dealii::deallog << std::endl;

  // call Kelly
  dealii::Vector<float> error(dof_handler.n_dofs());
  dealii::KellyErrorEstimator<dim>::estimate(dof_handler,
                                             face_quadrature_formula,
                                             typename dealii::FunctionMap<dim>::type (),
                                             values,
                                             error,
                                             dealii::ComponentMask(),
                                             0,
                                             dealii::numbers::invalid_unsigned_int,
                                             dealii::numbers::invalid_subdomain_id,
                                             dealii::numbers::invalid_material_id,
                                             dealii::KellyErrorEstimator<dim>::face_diameter_over_twice_max_degree);

  dealii::deallog <<"error:"<<std::endl;
  for (unsigned int i = 0; i < error.size(); i++)
    dealii::deallog << " " << error[i];
  dealii::deallog << std::endl;

  dof_handler.clear();
}


int main ()
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile);
  dealii::deallog.threshold_double(1e-8);

  {
    NeumanBC<2> func(1.25);
    test_neumann(func);
  }

  {
    MyFunction<2> func(0.25);
    test_regular(func);
  }

  {
    MyFunction<2> func(0.75);
    test_irregular(func);
  }

  deallog << "===3d==="<<std::endl;
  {
    NeumanBC<3> func(1.25);
    test_neumann(func);
  }

  {
    MyFunction<3> func(0.25);
    test_regular(func);
  }

  {
    MyFunction<3> func(0.75);
    test_irregular(func);
  }

  {
    MySecondFunction<2> function;
    test(function);
  }



  dealii::deallog << "Ok"<<std::endl;

}
