// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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


#include <deal.II/base/polynomials_bernardi_raugel.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN


template <int dim>
PolynomialsBernardiRaugel<dim>::PolynomialsBernardiRaugel(const unsigned int k)
  : TensorPolynomialsBase<dim>(k + 1, n_polynomials(k))
  , polynomial_space_Q(create_polynomials_Q())
  , polynomial_space_bubble(create_polynomials_bubble())
{}


template <int dim>
std::vector<std::vector<Polynomials::Polynomial<double>>>
PolynomialsBernardiRaugel<dim>::create_polynomials_bubble()
{
  std::vector<std::vector<Polynomials::Polynomial<double>>> pols;
  std::vector<Polynomials::Polynomial<double>>              bubble_shapes;
  bubble_shapes.push_back(Polynomials::LagrangeEquidistant(1, 0));
  bubble_shapes.push_back(Polynomials::LagrangeEquidistant(1, 1));
  bubble_shapes.push_back(Polynomials::LagrangeEquidistant(2, 1));

  for (unsigned int d = 0; d < dim; ++d)
    pols.push_back(bubble_shapes);
  // In 2D, the only q_ij polynomials we will use are 31,32,13,23
  // where ij corresponds to index (i-1)+3*(j-1) (2,5,6,7)

  // In 3D, the only q_ijk polynomials we will use are 331,332,313,323,133,233
  // where ijk corresponds to index (i-1)+3*(j-1)+9*(k-1)  (8,17,20,23,24,25)
  return pols;
}



template <int dim>
std::vector<std::vector<Polynomials::Polynomial<double>>>
PolynomialsBernardiRaugel<dim>::create_polynomials_Q()
{
  std::vector<std::vector<Polynomials::Polynomial<double>>> pols;
  std::vector<Polynomials::Polynomial<double>>              Q_shapes;
  Q_shapes.push_back(Polynomials::LagrangeEquidistant(1, 0));
  Q_shapes.push_back(Polynomials::LagrangeEquidistant(1, 1));
  for (unsigned int d = 0; d < dim; ++d)
    pols.push_back(Q_shapes);

  return pols;
}


template <int dim>
void
PolynomialsBernardiRaugel<dim>::evaluate(
  const Point<dim> &           unit_point,
  std::vector<Tensor<1, dim>> &values,
  std::vector<Tensor<2, dim>> &grads,
  std::vector<Tensor<3, dim>> &grad_grads,
  std::vector<Tensor<4, dim>> &third_derivatives,
  std::vector<Tensor<5, dim>> &fourth_derivatives) const
{
  Assert(values.size() == this->n() || values.size() == 0,
         ExcDimensionMismatch(values.size(), this->n()));
  Assert(grads.size() == this->n() || grads.size() == 0,
         ExcDimensionMismatch(grads.size(), this->n()));
  Assert(grad_grads.size() == this->n() || grad_grads.size() == 0,
         ExcDimensionMismatch(grad_grads.size(), this->n()));
  Assert(third_derivatives.size() == this->n() || third_derivatives.size() == 0,
         ExcDimensionMismatch(third_derivatives.size(), this->n()));
  Assert(fourth_derivatives.size() == this->n() ||
           fourth_derivatives.size() == 0,
         ExcDimensionMismatch(fourth_derivatives.size(), this->n()));

  std::vector<double>         Q_values;
  std::vector<Tensor<1, dim>> Q_grads;
  std::vector<Tensor<2, dim>> Q_grad_grads;
  std::vector<Tensor<3, dim>> Q_third_derivatives;
  std::vector<Tensor<4, dim>> Q_fourth_derivatives;
  std::vector<double>         bubble_values;
  std::vector<Tensor<1, dim>> bubble_grads;
  std::vector<Tensor<2, dim>> bubble_grad_grads;
  std::vector<Tensor<3, dim>> bubble_third_derivatives;
  std::vector<Tensor<4, dim>> bubble_fourth_derivatives;

  constexpr int n_bubbles =
    Utilities::pow(3, dim);     // size for create_polynomials_bubble
  constexpr int n_q = 1 << dim; // size for create_polynomials_q

  // don't resize if the provided vector has 0 length
  Q_values.resize((values.size() == 0) ? 0 : n_q);
  Q_grads.resize((grads.size() == 0) ? 0 : n_q);
  Q_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_q);
  Q_third_derivatives.resize((third_derivatives.size() == 0) ? 0 : n_q);
  Q_fourth_derivatives.resize((fourth_derivatives.size() == 0) ? 0 : n_q);
  bubble_values.resize((values.size() == 0) ? 0 : n_bubbles);
  bubble_grads.resize((grads.size() == 0) ? 0 : n_bubbles);
  bubble_grad_grads.resize((grad_grads.size() == 0) ? 0 : n_bubbles);
  bubble_third_derivatives.resize((third_derivatives.size() == 0) ? 0 :
                                                                    n_bubbles);
  bubble_fourth_derivatives.resize(
    (fourth_derivatives.size() == 0) ? 0 : n_bubbles);

  // 1 normal vector per face, ordering consistent with GeometryInfo
  // Normal vectors point in the +x, +y, and +z directions for
  // consistent orientation across edges
  std::vector<Tensor<1, dim>> normals;
  for (unsigned int i : GeometryInfo<dim>::face_indices())
    {
      Tensor<1, dim> normal;
      normal[i / 2] = 1;
      normals.push_back(normal);
    }

  // dim standard basis vectors for R^dim, usual ordering
  std::vector<Tensor<1, dim>> units;
  for (unsigned int i = 0; i < dim; ++i)
    {
      Tensor<1, dim> unit;
      unit[i] = 1;
      units.push_back(unit);
    }

  // set indices for the anisotropic polynomials to find
  // them after polynomial_space_bubble.evaluate is called
  std::vector<int> aniso_indices;
  if (dim == 2)
    {
      aniso_indices.push_back(6);
      aniso_indices.push_back(7);
      aniso_indices.push_back(2);
      aniso_indices.push_back(5);
    }
  else if (dim == 3)
    {
      aniso_indices.push_back(24);
      aniso_indices.push_back(25);
      aniso_indices.push_back(20);
      aniso_indices.push_back(23);
      aniso_indices.push_back(8);
      aniso_indices.push_back(17);
    }

  polynomial_space_bubble.evaluate(unit_point,
                                   bubble_values,
                                   bubble_grads,
                                   bubble_grad_grads,
                                   bubble_third_derivatives,
                                   bubble_fourth_derivatives);
  polynomial_space_Q.evaluate(unit_point,
                              Q_values,
                              Q_grads,
                              Q_grad_grads,
                              Q_third_derivatives,
                              Q_fourth_derivatives);

  // first dim*vertices_per_cell functions are Q_1^2 functions
  for (unsigned int i = 0; i < dim * GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      if (values.size() != 0)
        {
          values[i] = units[i % dim] * Q_values[i / dim];
        }
      if (grads.size() != 0)
        {
          grads[i] = outer_product(units[i % dim], Q_grads[i / dim]);
        }
      if (grad_grads.size() != 0)
        {
          grad_grads[i] = outer_product(units[i % dim], Q_grad_grads[i / dim]);
        }
      if (third_derivatives.size() != 0)
        {
          third_derivatives[i] =
            outer_product(units[i % dim], Q_third_derivatives[i / dim]);
        }
      if (fourth_derivatives.size() != 0)
        {
          fourth_derivatives[i] =
            outer_product(units[i % dim], Q_fourth_derivatives[i / dim]);
        }
    }

  // last faces_per_cell functions are bubble functions
  for (unsigned int i = dim * GeometryInfo<dim>::vertices_per_cell;
       i < dim * GeometryInfo<dim>::vertices_per_cell +
             GeometryInfo<dim>::faces_per_cell;
       ++i)
    {
      unsigned int j =
        i -
        dim *
          GeometryInfo<dim>::vertices_per_cell; // ranges 0 to faces_per_cell-1
      if (values.size() != 0)
        {
          values[i] = normals[j] * bubble_values[aniso_indices[j]];
        }
      if (grads.size() != 0)
        {
          grads[i] = outer_product(normals[j], bubble_grads[aniso_indices[j]]);
        }
      if (grad_grads.size() != 0)
        {
          grad_grads[i] =
            outer_product(normals[j], bubble_grad_grads[aniso_indices[j]]);
        }
      if (third_derivatives.size() != 0)
        {
          third_derivatives[i] =
            outer_product(normals[j],
                          bubble_third_derivatives[aniso_indices[j]]);
        }
      if (fourth_derivatives.size() != 0)
        {
          fourth_derivatives[i] =
            outer_product(normals[j],
                          bubble_fourth_derivatives[aniso_indices[j]]);
        }
    }
}

template <int dim>
unsigned int
PolynomialsBernardiRaugel<dim>::n_polynomials(const unsigned int k)
{
  (void)k;
  Assert(k == 1, ExcNotImplemented());
  if (dim == 2 || dim == 3)
    return dim * GeometryInfo<dim>::vertices_per_cell +
           GeometryInfo<dim>::faces_per_cell;
  // 2*4+4=12 polynomials in 2D and 3*8+6=30 polynomials in 3D

  Assert(false, ExcNotImplemented());
  return 0;
}


template <int dim>
std::unique_ptr<TensorPolynomialsBase<dim>>
PolynomialsBernardiRaugel<dim>::clone() const
{
  return std::make_unique<PolynomialsBernardiRaugel<dim>>(*this);
}

template class PolynomialsBernardiRaugel<1>; // to prevent errors
template class PolynomialsBernardiRaugel<2>;
template class PolynomialsBernardiRaugel<3>;


DEAL_II_NAMESPACE_CLOSE
