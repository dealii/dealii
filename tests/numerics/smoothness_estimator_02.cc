// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Check Fourier coefficient for simple functions in 1D/2D/3D.
// Further, estimate regularity in 1D/2D/3D.
// Similar to tests fe/fe_series_01.cc and numerics/smoothness_estimator_01.cc.

// Test functions: 1D: x^p, 2D: (x*y)^p, 3D: (x*y*z)^p
// Below is the MWE in Maxima for x*y*z:
/*********************************************************
integrate2(F,xx,aa,bb,yy,cc,dd):=block(integrate(integrate(F,xx,aa,bb),yy,cc,dd));
integrate3(F,xx,aa,bb,yy,cc,dd,zz,ee,ff):=block(integrate(integrate(integrate(F,xx,aa,bb),yy,cc,dd),zz,ee,ff));
a:0;
b:1;
nmax:3;
Phi(xx,nn):=exp(((-2)*%i*%pi*nn*xx)/(b-a));
Phi2(xx,nn,yy,mm):=Phi(xx,nn)*Phi(yy,mm);
Phi3(xx,nn,yy,mm,zz,ll):=Phi(xx,nn)*Phi(yy,mm)*Phi(zz,ll);
f3:x*y*z;
C3(n,m,l):=integrate3(f3*conjugate(Phi3(x,n,y,m,z,l)),x,a,b,y,a,b,z,a,b)/(b-a)^3;
load(functs);
for i:0 thru nmax do (for j:0 thru nmax do (for k:0 thru nmax do
print([i,j,k],fullratsimp(C3(i,j,k)))));
*********************************************************/


#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <iostream>

#include "../tests.h"


using namespace dealii;


template <int dim>
class PolyFunction : public Function<dim>
{
public:
  PolyFunction(const unsigned int poly_degree)
    : Function<dim>(1)
    , poly_degree(poly_degree)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const;

private:
  const unsigned int poly_degree;
};



template <int dim>
double
PolyFunction<dim>::value(const Point<dim> &point, const unsigned int) const
{
  double f = 1.;
  for (unsigned int d = 0; d < dim; ++d)
    f *= std::pow(point[d], poly_degree);
  return f;
}



template <typename CoefficientType>
void
  prepare_symmetric_coefficients(Table<1, CoefficientType> &         coeff,
                                 const std::vector<CoefficientType> &coeff_1d)
{
  Assert(coeff.size(0) == coeff_1d.size(), ExcInternalError());

  for (unsigned int i = 0; i < coeff.size(0); ++i)
    coeff(i) = coeff_1d[i];
}

template <typename CoefficientType>
void
  prepare_symmetric_coefficients(Table<2, CoefficientType> &         coeff,
                                 const std::vector<CoefficientType> &coeff_1d)
{
  for (unsigned int d = 0; d < 2; ++d)
    Assert(coeff.size(d) == coeff_1d.size(), ExcInternalError());

  for (unsigned int i = 0; i < coeff.size(0); ++i)
    for (unsigned int j = 0; j < coeff.size(1); ++j)
      coeff(i, j) = coeff_1d[i] * coeff_1d[j];
}

template <typename CoefficientType>
void
  prepare_symmetric_coefficients(Table<3, CoefficientType> &         coeff,
                                 const std::vector<CoefficientType> &coeff_1d)
{
  for (unsigned int d = 0; d < 3; ++d)
    Assert(coeff.size(d) == coeff_1d.size(), ExcInternalError());

  for (unsigned int i = 0; i < coeff.size(0); ++i)
    for (unsigned int j = 0; j < coeff.size(1); ++j)
      for (unsigned int k = 0; k < coeff.size(2); ++k)
        coeff(i, j, k) = coeff_1d[i] * coeff_1d[j] * coeff_1d[k];
}



template <typename CoefficientType>
typename CoefficientType::value_type
compare(const Table<1, CoefficientType> &coeff1,
        const Table<1, CoefficientType> &coeff2)
{
  Assert(coeff1.size(0) == coeff2.size(0), ExcInternalError());

  typename CoefficientType::value_type linf = 0.;
  for (unsigned int i = 0; i < coeff1.size(0); i++)
    linf = std::max(linf, std::abs(coeff1(i) - coeff2(i)));

  return linf;
}

template <typename CoefficientType>
typename CoefficientType::value_type
compare(const Table<2, CoefficientType> &coeff1,
        const Table<2, CoefficientType> &coeff2)
{
  for (unsigned int d = 0; d < 2; ++d)
    Assert(coeff1.size(d) == coeff2.size(d), ExcInternalError());

  typename CoefficientType::value_type linf = 0.;
  for (unsigned int i = 0; i < coeff1.size(0); i++)
    for (unsigned int j = 0; j < coeff1.size(1); j++)
      linf = std::max(linf, std::abs(coeff1(i, j) - coeff2(i, j)));

  return linf;
}

template <typename CoefficientType>
typename CoefficientType::value_type
compare(const Table<3, CoefficientType> &coeff1,
        const Table<3, CoefficientType> &coeff2)
{
  for (unsigned int d = 0; d < 3; ++d)
    Assert(coeff1.size(d) == coeff2.size(d), ExcInternalError());

  typename CoefficientType::value_type linf = 0.;
  for (unsigned int i = 0; i < coeff1.size(0); i++)
    for (unsigned int j = 0; j < coeff1.size(1); j++)
      for (unsigned int k = 0; k < coeff1.size(2); k++)
        linf = std::max(linf, std::abs(coeff1(i, j, k) - coeff2(i, j, k)));

  return linf;
}



template <int dim>
void
test(const unsigned int poly_degree)
{
  const unsigned int    max_poly = 3;
  hp::FECollection<dim> fe_collection;
  for (unsigned int p = 1; p <= max_poly; ++p)
    fe_collection.push_back(FE_Q<dim>(p));

  FESeries::Fourier<dim> fourier =
    SmoothnessEstimator::Fourier::default_fe_series(fe_collection);

  const unsigned int fe_index = poly_degree - 1;
  const unsigned int n_modes =
    fourier.get_n_coefficients_per_direction(fe_index);

  Assert((poly_degree >= 1) && (poly_degree <= max_poly), ExcInternalError());
  Assert((n_modes >= 3) && (n_modes <= max_poly + 1), ExcInternalError());

  deallog << "-----------------------------------" << std::endl;
  deallog << dim << "d, p=" << poly_degree << ", max_p=" << max_poly
          << ", n_modes=" << n_modes << std::endl;
  deallog << "-----------------------------------" << std::endl;

  // --- prepare test function ---
  PolyFunction<dim> test_function(poly_degree);

  // exact coefficients in 1D case
  const double &pi  = numbers::PI;
  const double  pi2 = std::pow(pi, 2);
  const double  pi3 = std::pow(pi, 3);

  std::vector<std::complex<double>> exact(n_modes);
  switch (poly_degree)
    {
      case 1:
        exact[0] = std::complex<double>(1., 0.) / 2.;
        if (n_modes > 1)
          exact[1] = std::complex<double>(0., -1.) / (2. * pi);
        if (n_modes > 2)
          exact[2] = std::complex<double>(0., -1.) / (4. * pi);
        if (n_modes > 3)
          exact[3] = std::complex<double>(0., -1.) / (6. * pi);
        break;
      case 2:
        exact[0] = std::complex<double>(1., 0.) / 3.;
        if (n_modes > 1)
          exact[1] = std::complex<double>(1., -pi) / (2. * pi2);
        if (n_modes > 2)
          exact[2] = std::complex<double>(1., -2. * pi) / (8. * pi2);
        if (n_modes > 3)
          exact[3] = std::complex<double>(1., -3. * pi) / (18. * pi2);
        break;
      case 3:
        exact[0] = std::complex<double>(1., 0.) / 4.;
        if (n_modes > 1)
          exact[1] = std::complex<double>(3. * pi, 3. - 2. * pi2) / (4. * pi3);
        if (n_modes > 2)
          exact[2] = std::complex<double>(6. * pi, 3. - 8. * pi2) / (32. * pi3);
        if (n_modes > 3)
          exact[3] = std::complex<double>(3. * pi, 1. - 6. * pi2) / (36. * pi3);
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }

  // coefficient table for multi-dimensional case
  TableIndices<dim> size;
  for (unsigned int d = 0; d < dim; ++d)
    size[d] = n_modes;
  Table<dim, std::complex<double>> coeff_in;
  coeff_in.reinit(size);
  prepare_symmetric_coefficients(coeff_in, exact);


  // --- prepare data structures ---
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0.0, 1.0); // reference cell

  hp::DoFHandler<dim> dof_handler(tria);
  dof_handler.begin_active()->set_active_fe_index(fe_index);
  dof_handler.distribute_dofs(fe_collection);


  // --- calculate coefficients from exact solution ---
  Vector<double> values(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, test_function, values);

  Vector<double> local_dof_values;
  auto           cell = dof_handler.begin_active();
  local_dof_values.reinit(cell->get_fe().dofs_per_cell);
  cell->get_dof_values(values, local_dof_values);

  Table<dim, std::complex<double>> coeff_out;
  coeff_out.reinit(size);
  fourier.calculate(local_dof_values, cell->active_fe_index(), coeff_out);

  // verify results
  const double linf = compare(coeff_in, coeff_out);
  deallog << "Linf norm in exact and calculate Fourier coefficients:"
          << std::endl
          << linf << std::endl;

  // finally test smoothness estimator:
  Vector<float> regularity(1);
  SmoothnessEstimator::Fourier::coefficient_decay(
    fourier,
    dof_handler,
    values,
    regularity,
    /*regression_strategy=*/VectorTools::Linfty_norm,
    /*smallest_abs_coefficient=*/1e-10,
    /*only_flagged_cells=*/false);

  deallog << "estimated regularity:" << std::endl << regularity[0] << std::endl;

  dof_handler.clear();
}



int
main()
{
  std::ofstream logfile("output");
  dealii::deallog.attach(logfile, /*do not print job id*/ false);
  dealii::deallog.depth_console(0);

  for (unsigned int poly_degree = 1; poly_degree <= 3; ++poly_degree)
    {
      test<1>(poly_degree);
      test<2>(poly_degree);
      test<3>(poly_degree);
    }

  dealii::deallog << "Ok" << std::endl;
}
