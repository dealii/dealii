// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/fe/fe_poly.h>

#include <deal.II/matrix_free/shape_info.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {

    // ----------------- actual ShapeInfo functions --------------------

    template <typename Number>
    ShapeInfo<Number>::ShapeInfo ()
      :
      n_q_points (0),
      dofs_per_cell (0)
    {}



    template <typename Number>
    template <int dim>
    void
    ShapeInfo<Number>::reinit (const Quadrature<1> &quad,
                               const FiniteElement<dim> &fe)
    {
      Assert (fe.n_components() == 1,
              ExcMessage("FEEvaluation only works for scalar finite elements."));

      const unsigned int n_dofs_1d = fe.degree+1,
                         n_q_points_1d = quad.size();
      AssertDimension(fe.dofs_per_cell, Utilities::fixed_power<dim>(n_dofs_1d));
      std::vector<unsigned int> lexicographic (fe.dofs_per_cell);

      // renumber (this is necessary for FE_Q, for example, since there the
      // vertex DoFs come first, which is incompatible with the lexicographic
      // ordering necessary to apply tensor products efficiently)
      {
        const FE_Poly<TensorProductPolynomials<dim>,dim,dim> *fe_poly =
          dynamic_cast<const FE_Poly<TensorProductPolynomials<dim>,dim,dim>*>(&fe);
        const FE_Poly<TensorProductPolynomials<dim,Polynomials::
        PiecewisePolynomial<double> >,dim,dim> *fe_poly_piece =
          dynamic_cast<const FE_Poly<TensorProductPolynomials<dim,
          Polynomials::PiecewisePolynomial<double> >,dim,dim>*> (&fe);
        Assert (fe_poly != 0 || fe_poly_piece, ExcNotImplemented());
        lexicographic = fe_poly != 0 ?
                        fe_poly->get_poly_space_numbering_inverse() :
                        fe_poly_piece->get_poly_space_numbering_inverse();

        // to evaluate 1D polynomials, evaluate along the line where y=z=0,
        // assuming that shape_value(0,Point<dim>()) == 1. otherwise, need
        // other entry point (e.g. generating a 1D element by reading the
        // name, as done before r29356)
        Assert(std::fabs(fe.shape_value(lexicographic[0], Point<dim>())-1) < 1e-13,
               ExcInternalError());
      }

      n_q_points      = Utilities::fixed_power<dim>(n_q_points_1d);
      dofs_per_cell   = Utilities::fixed_power<dim>(n_dofs_1d);
      n_q_points_face = dim>1?Utilities::fixed_power<dim-1>(n_q_points_1d):1;
      dofs_per_face   = dim>1?Utilities::fixed_power<dim-1>(n_dofs_1d):1;

      const unsigned int array_size = n_dofs_1d*n_q_points_1d;
      this->shape_gradients.resize_fast (array_size);
      this->shape_values.resize_fast (array_size);
      this->shape_hessians.resize_fast (array_size);

      this->face_value[0].resize(n_dofs_1d);
      this->face_gradient[0].resize(n_dofs_1d);
      this->subface_value[0].resize(array_size);
      this->face_value[1].resize(n_dofs_1d);
      this->face_gradient[1].resize(n_dofs_1d);
      this->subface_value[1].resize(array_size);
      this->shape_values_number.resize (array_size);
      this->shape_gradient_number.resize (array_size);

      for (unsigned int i=0; i<n_dofs_1d; ++i)
        {
          // need to reorder from hierarchical to lexicographic to get the
          // DoFs correct
          const unsigned int my_i = lexicographic[i];
          for (unsigned int q=0; q<n_q_points_1d; ++q)
            {
              // fill both vectors with
              // VectorizedArray<Number>::n_array_elements
              // copies for the shape information and
              // non-vectorized fields
              Point<dim> q_point;
              q_point[0] = quad.get_points()[q][0];
              shape_values_number[i*n_q_points_1d+q]   = fe.shape_value(my_i,q_point);
              shape_gradient_number[i*n_q_points_1d+q] = fe.shape_grad (my_i,q_point)[0];
              shape_values   [i*n_q_points_1d+q] =
                shape_values_number  [i*n_q_points_1d+q];
              shape_gradients[i*n_q_points_1d+q] =
                shape_gradient_number[i*n_q_points_1d+q];
              shape_hessians[i*n_q_points_1d+q] =
                fe.shape_grad_grad(my_i,q_point)[0][0];
              q_point[0] *= 0.5;
              subface_value[0][i*n_q_points_1d+q] = fe.shape_value(my_i,q_point);
              q_point[0] += 0.5;
              subface_value[1][i*n_q_points_1d+q] = fe.shape_value(my_i,q_point);
            }
          Point<dim> q_point;
          this->face_value[0][i] = fe.shape_value(my_i,q_point);
          this->face_gradient[0][i] = fe.shape_grad(my_i,q_point)[0];
          q_point[0] = 1;
          this->face_value[1][i] = fe.shape_value(my_i,q_point);
          this->face_gradient[1][i] = fe.shape_grad(my_i,q_point)[0];
        }

      // face information
      unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
      this->face_indices.reinit(n_faces, this->dofs_per_face);
      switch (dim)
        {
        case 3:
        {
          for (unsigned int i=0; i<this->dofs_per_face; i++)
            {
              const unsigned int jump_term =
                this->dofs_per_face*((i*n_dofs_1d)/this->dofs_per_face);
              this->face_indices(0,i) = i*n_dofs_1d;
              this->face_indices(1,i) = i*n_dofs_1d + n_dofs_1d-1;
              this->face_indices(2,i) = i%n_dofs_1d + jump_term;
              this->face_indices(3,i) = (i%n_dofs_1d + jump_term +
                                         (n_dofs_1d-1)*n_dofs_1d);
              this->face_indices(4,i) = i;
              this->face_indices(5,i) = (n_dofs_1d-1)*this->dofs_per_face+i;
            }
          break;
        }
        case 2:
        {
          for (unsigned int i=0; i<n_dofs_1d; i++)
            {
              this->face_indices(0,i) = n_dofs_1d*i;
              this->face_indices(1,i) = n_dofs_1d*i + n_dofs_1d-1;
              this->face_indices(2,i) = i;
              this->face_indices(3,i) = (n_dofs_1d-1)*n_dofs_1d+i;
            }
          break;
        }
        case 1:
        {
          this->face_indices(0,0) = 0;
          this->face_indices(1,0) = n_dofs_1d-1;
          break;
        }
        default:
          Assert (false, ExcNotImplemented());
        }
    }



    template <typename Number>
    std::size_t
    ShapeInfo<Number>::memory_consumption () const
    {
      std::size_t memory = sizeof(*this);
      memory += MemoryConsumption::memory_consumption(shape_values);
      memory += MemoryConsumption::memory_consumption(shape_gradients);
      memory += MemoryConsumption::memory_consumption(shape_hessians);
      memory += face_indices.memory_consumption();
      for (unsigned int i=0; i<2; ++i)
        {
          memory += MemoryConsumption::memory_consumption(face_value[i]);
          memory += MemoryConsumption::memory_consumption(face_gradient[i]);
        }
      memory += MemoryConsumption::memory_consumption(shape_values_number);
      memory += MemoryConsumption::memory_consumption(shape_gradient_number);
      return memory;
    }

    // end of functions for ShapeInfo

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE
