//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_tools.h>

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
                                    const FiniteElement<dim> &fe_dim)
  {
    Assert (fe_dim.n_components() == 1,
            ExcMessage("FEEvaluation only works for scalar finite elements."));

                                // take the name of the finite element
                                // and generate a 1d element. read the
                                // name, change the template argument
                                // to one and construct an element
    std::string fe_name = fe_dim.get_name();
    const std::size_t template_starts = fe_name.find_first_of('<');
    Assert (fe_name[template_starts+1] == (dim==1?'1':(dim==2?'2':'3')),
            ExcInternalError());
    fe_name[template_starts+1] = '1';
    std_cxx1x::shared_ptr<FiniteElement<1> > fe_1d
      (FETools::get_fe_from_name<1>(fe_name));
    const FiniteElement<1> & fe = *fe_1d;
    do_initialize (quad, fe, dim);
  }


  template <typename Number>
  void
  ShapeInfo<Number>::do_initialize (const Quadrature<1>    &quad,
                                           const FiniteElement<1> &fe,
                                           const unsigned int dim)
  {
    const unsigned int n_dofs_1d = fe.dofs_per_cell,
      n_q_points_1d = quad.size();
    std::vector<unsigned int> lexicographic (n_dofs_1d);

                                // renumber (this is necessary for FE_Q, for
                                // example, since there the vertex DoFs come
                                // first, which is incompatible with the
                                // lexicographic ordering necessary to apply
                                // tensor products efficiently)
    {
      const FE_Poly<TensorProductPolynomials<1>,1,1> *fe_poly =
        dynamic_cast<const FE_Poly<TensorProductPolynomials<1>,1,1>*>(&fe);
      Assert (fe_poly != 0, ExcNotImplemented());
      lexicographic = fe_poly->get_poly_space_numbering();
    }

    n_q_points      = 1;
    dofs_per_cell   = 1;
    n_q_points_face = 1;
    dofs_per_face   = 1;
    for (unsigned int d=0; d<dim; ++d)
      {
        n_q_points *= n_q_points_1d;
        dofs_per_cell *= n_dofs_1d;
      }
    for (int d=0; d<static_cast<int>(dim)-1; ++d)
      {
        n_q_points_face *= n_q_points_1d;
        dofs_per_face *= n_dofs_1d;
      }

    const unsigned int array_size = n_dofs_1d*n_q_points_1d;
    this->shape_gradients.resize_fast (array_size);
    this->shape_values.resize_fast (array_size);
    this->shape_hessians.resize_fast (array_size);

    this->face_gradient[0].resize(n_dofs_1d);
    this->face_value[0].resize(array_size);
    this->face_gradient[1].resize(n_dofs_1d);
    this->face_value[1].resize(array_size);
    this->shape_values_number.resize (array_size);
    this->shape_gradient_number.resize (array_size);

    for (unsigned int i=0; i<n_dofs_1d; ++i)
      {
                                // need to reorder from hierarchical to
                                // lexicographic to get the DoFs correct
        const unsigned int my_i = lexicographic[i];
        for (unsigned int q=0; q<n_q_points_1d; ++q)
          {
                                // fill both vectors with
                                // VectorizedArray<Number>::n_array_elements
                                // copies for the shape information and
                                // non-vectorized fields
            const Point<1> q_point = quad.get_points()[q];
            shape_values_number[my_i*n_q_points_1d+q]   = fe.shape_value(i,q_point);
            shape_gradient_number[my_i*n_q_points_1d+q] = fe.shape_grad (i,q_point)[0];
            shape_values   [my_i*n_q_points_1d+q] =
              shape_values_number  [my_i*n_q_points_1d+q];
            shape_gradients[my_i*n_q_points_1d+q] =
              shape_gradient_number[my_i*n_q_points_1d+q];
            shape_hessians[my_i*n_q_points_1d+q] =
              fe.shape_grad_grad(i,q_point)[0][0];
            face_value[0][my_i*n_q_points_1d+q] = fe.shape_value(i,q_point*0.5);
            face_value[1][my_i*n_q_points_1d+q] = fe.shape_value(i,Point<1>(0.5)+q_point*0.5);
          }
        this->face_gradient[0][my_i] = fe.shape_grad(i,Point<1>(0.))[0];
        this->face_gradient[1][my_i] = fe.shape_grad(i,Point<1>(1.))[0];
      }

                                // face information
    unsigned int n_faces = 1;
    for (unsigned int d=0; d<dim; ++d)
      n_faces *= 2;
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
