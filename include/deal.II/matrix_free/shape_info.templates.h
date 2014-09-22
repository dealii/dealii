// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q_dg0.h>

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
      element_type (tensor_general),
      n_q_points (0),
      dofs_per_cell (0)
    {}



    template <typename Number>
    template <int dim>
    void
    ShapeInfo<Number>::reinit (const Quadrature<1> &quad,
                               const FiniteElement<dim> &fe_in,
                               const unsigned int base_element_number)
    {
      const FiniteElement<dim> *fe = &fe_in;
      fe = &fe_in.base_element(base_element_number);

      Assert (fe->n_components() == 1,
              ExcMessage("FEEvaluation only works for scalar finite elements."));

      fe_degree = fe->degree;

      const unsigned int n_dofs_1d = fe_degree+1,
                         n_q_points_1d = quad.size();

      // renumber (this is necessary for FE_Q, for example, since there the
      // vertex DoFs come first, which is incompatible with the lexicographic
      // ordering necessary to apply tensor products efficiently)
      std::vector<unsigned int> scalar_lexicographic;
      {
        // find numbering to lexicographic
        Assert(fe->n_components() == 1,
               ExcMessage("Expected a scalar element"));

        const FE_Poly<TensorProductPolynomials<dim>,dim,dim> *fe_poly =
          dynamic_cast<const FE_Poly<TensorProductPolynomials<dim>,dim,dim>*>(fe);

        const FE_Poly<TensorProductPolynomials<dim,Polynomials::
        PiecewisePolynomial<double> >,dim,dim> *fe_poly_piece =
          dynamic_cast<const FE_Poly<TensorProductPolynomials<dim,
          Polynomials::PiecewisePolynomial<double> >,dim,dim>*> (fe);

        const FE_DGP<dim> *fe_dgp = dynamic_cast<const FE_DGP<dim>*>(fe);

        const FE_Q_DG0<dim> *fe_q_dg0 = dynamic_cast<const FE_Q_DG0<dim>*>(fe);

        element_type = tensor_general;
        if (fe_poly != 0)
          scalar_lexicographic = fe_poly->get_poly_space_numbering_inverse();
        else if (fe_poly_piece != 0)
          scalar_lexicographic = fe_poly_piece->get_poly_space_numbering_inverse();
        else if (fe_dgp != 0)
          {
            scalar_lexicographic.resize(fe_dgp->dofs_per_cell);
            for (unsigned int i=0; i<fe_dgp->dofs_per_cell; ++i)
              scalar_lexicographic[i] = i;
            element_type = truncated_tensor;
          }
        else if (fe_q_dg0 != 0)
          {
            scalar_lexicographic = fe_q_dg0->get_poly_space_numbering_inverse();
            element_type = tensor_symmetric_plus_dg0;
          }
        else
          Assert(false, ExcNotImplemented());

        // Finally store the renumbering into the member variable of this
        // class
        if (fe_in.n_components() == 1)
          lexicographic_numbering = scalar_lexicographic;
        else
          {
            // have more than one component, get the inverse
            // permutation, invert it, sort the components one after one,
            // and invert back
            std::vector<unsigned int> scalar_inv =
              Utilities::invert_permutation(scalar_lexicographic);
            std::vector<unsigned int> lexicographic (fe_in.dofs_per_cell);
            for (unsigned int comp=0; comp<fe_in.n_components(); ++comp)
              for (unsigned int i=0; i<scalar_inv.size(); ++i)
                lexicographic[fe_in.component_to_system_index(comp,i)]
                  = scalar_inv.size () * comp + scalar_inv[i];

            // invert numbering again
            lexicographic_numbering =
              Utilities::invert_permutation(lexicographic);
          }

        // to evaluate 1D polynomials, evaluate along the line where y=z=0,
        // assuming that shape_value(0,Point<dim>()) == 1. otherwise, need
        // other entry point (e.g. generating a 1D element by reading the
        // name, as done before r29356)
        Assert(std::fabs(fe->shape_value(scalar_lexicographic[0],
                                         Point<dim>())-1) < 1e-13,
               ExcInternalError());
      }

      n_q_points      = Utilities::fixed_power<dim>(n_q_points_1d);
      dofs_per_cell   = fe->dofs_per_cell;
      n_q_points_face = dim>1?Utilities::fixed_power<dim-1>(n_q_points_1d):1;
      dofs_per_face   = fe->dofs_per_face;

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
          const unsigned int my_i = scalar_lexicographic[i];
          for (unsigned int q=0; q<n_q_points_1d; ++q)
            {
              // fill both vectors with
              // VectorizedArray<Number>::n_array_elements
              // copies for the shape information and
              // non-vectorized fields
              Point<dim> q_point;
              q_point[0] = quad.get_points()[q][0];
              shape_values_number[i*n_q_points_1d+q]   = fe->shape_value(my_i,q_point);
              shape_gradient_number[i*n_q_points_1d+q] = fe->shape_grad (my_i,q_point)[0];
              shape_values   [i*n_q_points_1d+q] =
                shape_values_number  [i*n_q_points_1d+q];
              shape_gradients[i*n_q_points_1d+q] =
                shape_gradient_number[i*n_q_points_1d+q];
              shape_hessians[i*n_q_points_1d+q] =
                fe->shape_grad_grad(my_i,q_point)[0][0];
              q_point[0] *= 0.5;
              subface_value[0][i*n_q_points_1d+q] = fe->shape_value(my_i,q_point);
              q_point[0] += 0.5;
              subface_value[1][i*n_q_points_1d+q] = fe->shape_value(my_i,q_point);
            }
          Point<dim> q_point;
          this->face_value[0][i] = fe->shape_value(my_i,q_point);
          this->face_gradient[0][i] = fe->shape_grad(my_i,q_point)[0];
          q_point[0] = 1;
          this->face_value[1][i] = fe->shape_value(my_i,q_point);
          this->face_gradient[1][i] = fe->shape_grad(my_i,q_point)[0];
        }

      if (element_type == tensor_general &&
          check_1d_shapes_symmetric(n_q_points_1d))
        {
          if (check_1d_shapes_gausslobatto())
            element_type = tensor_gausslobatto;
          else
            element_type = tensor_symmetric;
        }
      else if (element_type == tensor_symmetric_plus_dg0)
        check_1d_shapes_symmetric(n_q_points_1d);

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
          for (unsigned int i=0; i<this->dofs_per_face; i++)
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
    bool
    ShapeInfo<Number>::check_1d_shapes_symmetric(const unsigned int n_q_points_1d)
    {
      const double zero_tol =
        types_are_equal<Number,double>::value==true?1e-10:1e-7;
      // symmetry for values
      const unsigned int n_dofs_1d = fe_degree + 1;
      for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
        for (unsigned int j=0; j<n_q_points_1d; ++j)
          if (std::fabs(shape_values[i*n_q_points_1d+j][0] -
                        shape_values[(n_dofs_1d-i)*n_q_points_1d
                                     -j-1][0]) > zero_tol)
            return false;

      // shape values should be zero at x=0.5 for all basis functions except
      // for one which is one
      if (n_q_points_1d%2 == 1 && n_dofs_1d%2 == 1)
        {
          for (unsigned int i=0; i<n_dofs_1d/2; ++i)
            if (std::fabs(shape_values[i*n_q_points_1d+
                                       n_q_points_1d/2][0]) > zero_tol)
              return false;
          if (std::fabs(shape_values[(n_dofs_1d/2)*n_q_points_1d+
                                     n_q_points_1d/2][0]-1.)> zero_tol)
            return false;
        }

      // skew-symmetry for gradient, zero of middle basis function in middle
      // quadrature point
      for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
        for (unsigned int j=0; j<n_q_points_1d; ++j)
          if (std::fabs(shape_gradients[i*n_q_points_1d+j][0] +
                        shape_gradients[(n_dofs_1d-i)*n_q_points_1d-
                                        j-1][0]) > zero_tol)
            return false;
      if (n_dofs_1d%2 == 1 && n_q_points_1d%2 == 1)
        if (std::fabs(shape_gradients[(n_dofs_1d/2)*n_q_points_1d+
                                      (n_q_points_1d/2)][0]) > zero_tol)
          return false;

      // symmetry for Laplacian
      for (unsigned int i=0; i<(n_dofs_1d+1)/2; ++i)
        for (unsigned int j=0; j<n_q_points_1d; ++j)
          if (std::fabs(shape_hessians[i*n_q_points_1d+j][0] -
                        shape_hessians[(n_dofs_1d-i)*n_q_points_1d-
                                       j-1][0]) > zero_tol)
            return false;

      const unsigned int stride = (n_q_points_1d+1)/2;
      shape_val_evenodd.resize((fe_degree+1)*stride);
      shape_gra_evenodd.resize((fe_degree+1)*stride);
      shape_hes_evenodd.resize((fe_degree+1)*stride);
      for (unsigned int i=0; i<(fe_degree+1)/2; ++i)
        for (unsigned int q=0; q<stride; ++q)
          {
            shape_val_evenodd[i*stride+q] =
              0.5 * (shape_values[i*n_q_points_1d+q] +
                     shape_values[i*n_q_points_1d+n_q_points_1d-1-q]);
            shape_val_evenodd[(fe_degree-i)*stride+q] =
              0.5 * (shape_values[i*n_q_points_1d+q] -
                     shape_values[i*n_q_points_1d+n_q_points_1d-1-q]);

            shape_gra_evenodd[i*stride+q] =
              0.5 * (shape_gradients[i*n_q_points_1d+q] +
                     shape_gradients[i*n_q_points_1d+n_q_points_1d-1-q]);
            shape_gra_evenodd[(fe_degree-i)*stride+q] =
              0.5 * (shape_gradients[i*n_q_points_1d+q] -
                     shape_gradients[i*n_q_points_1d+n_q_points_1d-1-q]);

            shape_hes_evenodd[i*stride+q] =
              0.5 * (shape_hessians[i*n_q_points_1d+q] +
                     shape_hessians[i*n_q_points_1d+n_q_points_1d-1-q]);
            shape_hes_evenodd[(fe_degree-i)*stride+q] =
              0.5 * (shape_hessians[i*n_q_points_1d+q] -
                     shape_hessians[i*n_q_points_1d+n_q_points_1d-1-q]);
          }
      if (fe_degree % 2 == 0)
        for (unsigned int q=0; q<stride; ++q)
          {
            shape_val_evenodd[fe_degree/2*stride+q] =
              shape_values[(fe_degree/2)*n_q_points_1d+q];
            shape_gra_evenodd[fe_degree/2*stride+q] =
              shape_gradients[(fe_degree/2)*n_q_points_1d+q];
            shape_hes_evenodd[fe_degree/2*stride+q] =
              shape_hessians[(fe_degree/2)*n_q_points_1d+q];
          }

      return true;
    }



    template <typename Number>
    bool
    ShapeInfo<Number>::check_1d_shapes_gausslobatto()
    {
      if (dofs_per_cell != n_q_points)
        return false;

      const double zero_tol =
        types_are_equal<Number,double>::value==true?1e-10:1e-7;
      // check: - identity operation for shape values
      //        - zero diagonal at interior points for gradients
      //        - gradient equal to unity at element boundary
      const unsigned int n_points_1d = fe_degree+1;
      for (unsigned int i=0; i<n_points_1d; ++i)
        for (unsigned int j=0; j<n_points_1d; ++j)
          if (i!=j)
            {
              if (std::fabs(shape_values[i*n_points_1d+j][0])>zero_tol)
                return false;
            }
          else
            {
              if (std::fabs(shape_values[i*n_points_1d+
                                         j][0]-1.)>zero_tol)
                return false;
            }
      for (unsigned int i=1; i<n_points_1d-1; ++i)
        if (std::fabs(shape_gradients[i*n_points_1d+i][0])>zero_tol)
          return false;
      if (std::fabs(shape_gradients[n_points_1d-1][0]-
                    (n_points_1d%2==0 ? -1. : 1.)) > zero_tol)
        return false;

      return true;
    }



    template <typename Number>
    std::size_t
    ShapeInfo<Number>::memory_consumption () const
    {
      std::size_t memory = sizeof(*this);
      memory += MemoryConsumption::memory_consumption(shape_values);
      memory += MemoryConsumption::memory_consumption(shape_gradients);
      memory += MemoryConsumption::memory_consumption(shape_hessians);
      memory += MemoryConsumption::memory_consumption(shape_val_evenodd);
      memory += MemoryConsumption::memory_consumption(shape_gra_evenodd);
      memory += MemoryConsumption::memory_consumption(shape_hes_evenodd);
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
