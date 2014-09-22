// ---------------------------------------------------------------------
//
// Copyright (C) @YEAR@ by the deal.II authors
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

// Generic routines to check consistency of function classes


// Check, whether the various implementations of function values are
// consistent. Arguments include the function, the amount of
// quadrature points in each direction, and the threshold above which
// we consider values unequal

template<int dim>
void
check_function_value_consistency(
  const Function<dim> &f,
  unsigned int sub,
  double threshold = 1.e-15)
{
  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);

  std::vector<double> f1(quadrature.size());
  std::vector<Vector<double> > f2(quadrature.size(), Vector<double>(f.n_components));

  f.vector_value_list(quadrature.get_points(), f2);

  deallog << "value vs vector value list";
  for (unsigned int d=0; d<f.n_components; ++d)
    for (unsigned int i=0; i<f1.size(); ++i)
      {
        const double v = f.value(quadrature.point(i), d);
        if (std::fabs(v-f2[i](d)) > threshold)
          deallog << "v-vl " << d << ':' << i << ':' << v-f2[i](d);
      }
  deallog << std::endl << "value list vs vector value list";
  for (unsigned int d=0; d<f.n_components; ++d)
    {
      f.value_list(quadrature.get_points(), f1, d);
      for (unsigned int i=0; i<f1.size(); ++i)
        {
          if (std::fabs(f1[i]-f2[i](d)) > threshold)
            deallog << ' ' << d << ':' << i << ':' << f1[i]-f2[i](d);
        }
    }
  deallog << std::endl;
}

// Same for gradients
template<int dim>
void
check_function_gradient_consistency(
  const Function<dim> &f,
  unsigned int sub,
  double threshold = 1.e-15)
{
  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);

  std::vector<Tensor<1,dim> > f1(quadrature.size());
  std::vector<std::vector<Tensor<1,dim> > > f2(quadrature.size(),
                                               std::vector<Tensor<1,dim> >(f.n_components));

  f.vector_gradient_list(quadrature.get_points(), f2);

  deallog << "gradient vs vector gradient list";
  for (unsigned int d=0; d<f.n_components; ++d)
    for (unsigned int i=0; i<f1.size(); ++i)
      {
        const Tensor<1,dim> v = f.gradient(quadrature.point(i), d)-f2[i][d];

        if (std::sqrt(v*v) > threshold)
          deallog << "v-vl " << d << ':' << i << ':' << v;
      }
  deallog << std::endl << "gradient list vs vector gradient list";
  for (unsigned int d=0; d<f.n_components; ++d)
    {
      f.gradient_list(quadrature.get_points(), f1, d);
      for (unsigned int i=0; i<f1.size(); ++i)
        {
          const Tensor<1,dim> v = f1[i]-f2[i][d];
          if (std::sqrt(v*v) > threshold)
            deallog << ' ' << d << ':' << i << ':' << v;
        }
    }
  deallog << std::endl;
}

/**
 * A class replacing the implemented derivatives of a function with
 * difference quotients. This way, the correctness of the
 * implementation can be tested.
 */
template <int dim>
class DerivativeTestFunction :
  public AutoDerivativeFunction<dim>
{
public:
  DerivativeTestFunction(const Function<dim> &, const double h);
  ~DerivativeTestFunction();

  virtual void vector_value (const Point<dim> &points, Vector<double> &value) const;
  virtual double value (const Point<dim> &points, const unsigned int component) const;
  virtual void vector_value_list (const std::vector< Point< dim > > &points,
                                  std::vector< Vector< double > > &values) const;

private:
  const Function<dim> &func;
};


template <int dim>
DerivativeTestFunction<dim>::DerivativeTestFunction(const Function<dim> &f,
                                                    const double h)
  :
  AutoDerivativeFunction<dim>(h, f.n_components),
  func(f)
{
  this->set_formula(AutoDerivativeFunction<dim>::FourthOrder);
}


template <int dim>
DerivativeTestFunction<dim>::~DerivativeTestFunction()
{}


template <int dim>
void
DerivativeTestFunction<dim>::vector_value_list (
  const std::vector< Point< dim > > &points,
  std::vector< Vector< double > > &values) const
{
  func.vector_value_list(points, values);
}


template<int dim>
void DerivativeTestFunction<dim>::vector_value (
  const Point<dim> &point,
  Vector<double> &value) const
{
  func.vector_value(point, value);
}


template<int dim>
double DerivativeTestFunction<dim>::value (
  const Point<dim> &point,
  const unsigned int comp) const
{
//  std::cerr << '[' << point << '!' << func.value(point, comp) << ']';

  return func.value(point, comp);
}


// Check whether the difference quotients converge to the gradient
template<int dim>
void
check_gradient(
  const Function<dim> &f,
  unsigned int sub,
  double threshold = 1./14.)
{
  DerivativeTestFunction<dim> dtest1(f, 1.e-2);
  DerivativeTestFunction<dim> dtest2(f, 2.e-2);

  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);
  const std::vector<Point<dim> > &points = quadrature.get_points();

  std::vector<std::vector<Tensor<1,dim> > >
  gradients(f.n_components, std::vector<Tensor<1,dim> >(points.size()));
  std::vector<std::vector<Tensor<1,dim> > >
  gradients1(f.n_components, std::vector<Tensor<1,dim> >(points.size()));
  std::vector<std::vector<Tensor<1,dim> > >
  gradients2(f.n_components, std::vector<Tensor<1,dim> >(points.size()));

  deallog << "gradients vs difference quotients";

  f.vector_gradients(points, gradients);
  dtest1.vector_gradients(points, gradients1);
  dtest2.vector_gradients(points, gradients2);

  // Compare gradients and difference quotients
  for (unsigned int k=0; k<gradients.size(); ++k)
    for (unsigned int i=0; i<gradients[k].size(); ++i)
      {
        // Compute difference
        Tensor<1,dim> d1 = gradients1[k][i] - gradients[k][i];
        Tensor<1,dim> d2 = gradients2[k][i] - gradients[k][i];

        // If the difference is
        // already small, we are fine
        if (d1.norm() > 1.e-13)
          {
            // Check for
            // convergence. For full
            // 4th order, gradients2
            // should be 16 times as
            // large, so let's be a
            // bit generous
            if (threshold * d2.norm() < d1.norm())
              {
                deallog << "Gradient error: point " << i
                        << " (" << points[i] << " )"
                        << " comp " << k
//      << " norms " << d1.norm() << " " << d2.norm()
                        << std::endl;
                for (unsigned int d=0; d<dim; ++d)
                  deallog
                      << " " << gradients[k][i][d]
                      << " " << gradients1[i][k][d]
                      << std::endl;
              }
          }
      }
  deallog << std::endl;
}

