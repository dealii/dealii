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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <cmath>
#include <algorithm>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim>
const unsigned int MappingQ1<dim,spacedim>::n_shape_functions;



template<int dim, int spacedim>
MappingQ1<dim,spacedim>::InternalData::InternalData (const unsigned int n_shape_functions)
  :
  is_mapping_q1_data(true),
  n_shape_functions (n_shape_functions)
{}



template<int dim, int spacedim>
std::size_t
MappingQ1<dim,spacedim>::InternalData::memory_consumption () const
{
  return (Mapping<dim,spacedim>::InternalDataBase::memory_consumption() +
          MemoryConsumption::memory_consumption (shape_values) +
          MemoryConsumption::memory_consumption (shape_derivatives) +
          MemoryConsumption::memory_consumption (covariant) +
          MemoryConsumption::memory_consumption (contravariant) +
          MemoryConsumption::memory_consumption (unit_tangentials) +
          MemoryConsumption::memory_consumption (aux) +
          MemoryConsumption::memory_consumption (mapping_support_points) +
          MemoryConsumption::memory_consumption (cell_of_current_support_points) +
          MemoryConsumption::memory_consumption (volume_elements) +
          MemoryConsumption::memory_consumption (is_mapping_q1_data) +
          MemoryConsumption::memory_consumption (n_shape_functions));
}



template<int dim, int spacedim>
MappingQ1<dim,spacedim>::MappingQ1 ()
{}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_shapes (const std::vector<Point<dim> > &unit_points,
                                         InternalData &data) const
{
  // choose either the function implemented
  // in this class, or whatever a virtual
  // function call resolves to
  if (data.is_mapping_q1_data)
    MappingQ1<dim,spacedim>::compute_shapes_virtual(unit_points, data);
  else
    compute_shapes_virtual(unit_points, data);
}


namespace internal
{
  namespace MappingQ1
  {
    // These are left as templates on the spatial dimension (even though dim
    // == spacedim must be true for them to make sense) because templates are
    // expanded before the compiler eliminates code due to the 'if (dim ==
    // spacedim)' statement (see the body of the general
    // transform_real_to_unit_cell).
    template<int spacedim>
    Point<1>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<1>::vertices_per_cell> &vertices,
     const Point<spacedim> &p)
    {
      Assert(spacedim == 1, ExcInternalError());
      return Point<1>((p[0] - vertices[0](0))/(vertices[1](0) - vertices[0](0)));
    }



    template<int spacedim>
    Point<2>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<2>::vertices_per_cell> &vertices,
     const Point<spacedim> &p)
    {
      Assert(spacedim == 2, ExcInternalError());
      const double x = p(0);
      const double y = p(1);

      const double x0 = vertices[0](0);
      const double x1 = vertices[1](0);
      const double x2 = vertices[2](0);
      const double x3 = vertices[3](0);

      const double y0 = vertices[0](1);
      const double y1 = vertices[1](1);
      const double y2 = vertices[2](1);
      const double y3 = vertices[3](1);

      const double a = (x1 - x3)*(y0 - y2) - (x0 - x2)*(y1 - y3);
      const double b = -(x0 - x1 - x2 + x3)*y + (x - 2*x1 + x3)*y0 - (x - 2*x0 + x2)*y1
                       - (x - x1)*y2 + (x - x0)*y3;
      const double c = (x0 - x1)*y - (x - x1)*y0 + (x - x0)*y1;

      const double discriminant = b*b - 4*a*c;
      // fast exit if the point is not in the cell (this is the only case
      // where the discriminant is negative)
      if (discriminant < 0.0)
        {
          return Point<2>(2, 2);
        }

      double eta1;
      double eta2;
      // special case #1: if a is zero, then use the linear formula
      if (a == 0.0 && b != 0.0)
        {
          eta1 = -c/b;
          eta2 = -c/b;
        }
      // special case #2: if c is very small:
      else if (std::abs(c/b) < 1e-12)
        {
          eta1 = (-b - std::sqrt(discriminant)) / (2*a);
          eta2 = (-b + std::sqrt(discriminant)) / (2*a);
        }
      // finally, use the numerically stable version of the quadratic formula:
      else
        {
          eta1 = 2*c / (-b - std::sqrt(discriminant));
          eta2 = 2*c / (-b + std::sqrt(discriminant));
        }
      // pick the one closer to the center of the cell.
      const double eta = (std::abs(eta1 - 0.5) < std::abs(eta2 - 0.5)) ? eta1 : eta2;

      /*
       * There are two ways to compute xi from eta, but either one may have a
       * zero denominator.
       */
      const double subexpr0 = -eta*x2 + x0*(eta - 1);
      const double xi_denominator0 = eta*x3 - x1*(eta - 1) + subexpr0;
      const double max_x = std::max(std::max(std::abs(x0), std::abs(x1)),
                                    std::max(std::abs(x2), std::abs(x3)));

      if (std::abs(xi_denominator0) > 1e-10*max_x)
        {
          const double xi = (x + subexpr0)/xi_denominator0;
          return Point<2>(xi, eta);
        }
      else
        {
          const double max_y = std::max(std::max(std::abs(y0), std::abs(y1)),
                                        std::max(std::abs(y2), std::abs(y3)));
          const double subexpr1 = -eta*y2 + y0*(eta - 1);
          const double xi_denominator1 = eta*y3 - y1*(eta - 1) + subexpr1;
          if (std::abs(xi_denominator1) > 1e-10*max_y)
            {
              const double xi = (subexpr1 + y)/xi_denominator1;
              return Point<2>(xi, eta);
            }
          else // give up and try Newton iteration
            {
              return Point<2>(2, 2);
            }
        }
    }



    template<int spacedim>
    Point<3>
    transform_real_to_unit_cell
    (const std_cxx11::array<Point<spacedim>, GeometryInfo<3>::vertices_per_cell> &/*vertices*/,
     const Point<spacedim> &/*p*/)
    {
      // It should not be possible to get here
      Assert(false, ExcInternalError());
      return Point<3>();
    }



    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
                            const std::vector<Point<1> > &unit_points,
                            typename dealii::MappingQ1<1,spacedim>::InternalData &data)
    {
      (void)n_shape_functions;
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
        {
          double x = unit_points[k](0);

          if (data.shape_values.size()!=0)
            {
              Assert(data.shape_values.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.shape(k,0) = 1.-x;
              data.shape(k,1) = x;
            }
          if (data.shape_derivatives.size()!=0)
            {
              Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.derivative(k,0)[0] = -1.;
              data.derivative(k,1)[0] = 1.;
            }
          if (data.shape_second_derivatives.size()!=0)
            {
              // the following may or may not
              // work if dim != spacedim
              Assert (spacedim == 1, ExcNotImplemented());

              Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.second_derivative(k,0)[0][0] = 0;
              data.second_derivative(k,1)[0][0] = 0;
            }
        }
    }


    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
                            const std::vector<Point<2> > &unit_points,
                            typename dealii::MappingQ1<2,spacedim>::InternalData &data)
    {
      (void)n_shape_functions;
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
        {
          double x = unit_points[k](0);
          double y = unit_points[k](1);

          if (data.shape_values.size()!=0)
            {
              Assert(data.shape_values.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.shape(k,0) = (1.-x)*(1.-y);
              data.shape(k,1) = x*(1.-y);
              data.shape(k,2) = (1.-x)*y;
              data.shape(k,3) = x*y;
            }
          if (data.shape_derivatives.size()!=0)
            {
              Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.derivative(k,0)[0] = (y-1.);
              data.derivative(k,1)[0] = (1.-y);
              data.derivative(k,2)[0] = -y;
              data.derivative(k,3)[0] = y;
              data.derivative(k,0)[1] = (x-1.);
              data.derivative(k,1)[1] = -x;
              data.derivative(k,2)[1] = (1.-x);
              data.derivative(k,3)[1] = x;
            }
          if (data.shape_second_derivatives.size()!=0)
            {
              Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.second_derivative(k,0)[0][0] = 0;
              data.second_derivative(k,1)[0][0] = 0;
              data.second_derivative(k,2)[0][0] = 0;
              data.second_derivative(k,3)[0][0] = 0;
              data.second_derivative(k,0)[0][1] = 1.;
              data.second_derivative(k,1)[0][1] = -1.;
              data.second_derivative(k,2)[0][1] = -1.;
              data.second_derivative(k,3)[0][1] = 1.;
              data.second_derivative(k,0)[1][0] = 1.;
              data.second_derivative(k,1)[1][0] = -1.;
              data.second_derivative(k,2)[1][0] = -1.;
              data.second_derivative(k,3)[1][0] = 1.;
              data.second_derivative(k,0)[1][1] = 0;
              data.second_derivative(k,1)[1][1] = 0;
              data.second_derivative(k,2)[1][1] = 0;
              data.second_derivative(k,3)[1][1] = 0;
            }
        }
    }



    template <int spacedim>
    void
    compute_shapes_virtual (const unsigned int            n_shape_functions,
                            const std::vector<Point<3> > &unit_points,
                            typename dealii::MappingQ1<3,spacedim>::InternalData &data)
    {
      (void)n_shape_functions;
      const unsigned int n_points=unit_points.size();
      for (unsigned int k = 0 ; k < n_points ; ++k)
        {
          double x = unit_points[k](0);
          double y = unit_points[k](1);
          double z = unit_points[k](2);

          if (data.shape_values.size()!=0)
            {
              Assert(data.shape_values.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.shape(k,0) = (1.-x)*(1.-y)*(1.-z);
              data.shape(k,1) = x*(1.-y)*(1.-z);
              data.shape(k,2) = (1.-x)*y*(1.-z);
              data.shape(k,3) = x*y*(1.-z);
              data.shape(k,4) = (1.-x)*(1.-y)*z;
              data.shape(k,5) = x*(1.-y)*z;
              data.shape(k,6) = (1.-x)*y*z;
              data.shape(k,7) = x*y*z;
            }
          if (data.shape_derivatives.size()!=0)
            {
              Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.derivative(k,0)[0] = (y-1.)*(1.-z);
              data.derivative(k,1)[0] = (1.-y)*(1.-z);
              data.derivative(k,2)[0] = -y*(1.-z);
              data.derivative(k,3)[0] = y*(1.-z);
              data.derivative(k,4)[0] = (y-1.)*z;
              data.derivative(k,5)[0] = (1.-y)*z;
              data.derivative(k,6)[0] = -y*z;
              data.derivative(k,7)[0] = y*z;
              data.derivative(k,0)[1] = (x-1.)*(1.-z);
              data.derivative(k,1)[1] = -x*(1.-z);
              data.derivative(k,2)[1] = (1.-x)*(1.-z);
              data.derivative(k,3)[1] = x*(1.-z);
              data.derivative(k,4)[1] = (x-1.)*z;
              data.derivative(k,5)[1] = -x*z;
              data.derivative(k,6)[1] = (1.-x)*z;
              data.derivative(k,7)[1] = x*z;
              data.derivative(k,0)[2] = (x-1)*(1.-y);
              data.derivative(k,1)[2] = x*(y-1.);
              data.derivative(k,2)[2] = (x-1.)*y;
              data.derivative(k,3)[2] = -x*y;
              data.derivative(k,4)[2] = (1.-x)*(1.-y);
              data.derivative(k,5)[2] = x*(1.-y);
              data.derivative(k,6)[2] = (1.-x)*y;
              data.derivative(k,7)[2] = x*y;
            }
          if (data.shape_second_derivatives.size()!=0)
            {
              // the following may or may not
              // work if dim != spacedim
              Assert (spacedim == 3, ExcNotImplemented());

              Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
                     ExcInternalError());
              data.second_derivative(k,0)[0][0] = 0;
              data.second_derivative(k,1)[0][0] = 0;
              data.second_derivative(k,2)[0][0] = 0;
              data.second_derivative(k,3)[0][0] = 0;
              data.second_derivative(k,4)[0][0] = 0;
              data.second_derivative(k,5)[0][0] = 0;
              data.second_derivative(k,6)[0][0] = 0;
              data.second_derivative(k,7)[0][0] = 0;
              data.second_derivative(k,0)[1][1] = 0;
              data.second_derivative(k,1)[1][1] = 0;
              data.second_derivative(k,2)[1][1] = 0;
              data.second_derivative(k,3)[1][1] = 0;
              data.second_derivative(k,4)[1][1] = 0;
              data.second_derivative(k,5)[1][1] = 0;
              data.second_derivative(k,6)[1][1] = 0;
              data.second_derivative(k,7)[1][1] = 0;
              data.second_derivative(k,0)[2][2] = 0;
              data.second_derivative(k,1)[2][2] = 0;
              data.second_derivative(k,2)[2][2] = 0;
              data.second_derivative(k,3)[2][2] = 0;
              data.second_derivative(k,4)[2][2] = 0;
              data.second_derivative(k,5)[2][2] = 0;
              data.second_derivative(k,6)[2][2] = 0;
              data.second_derivative(k,7)[2][2] = 0;

              data.second_derivative(k,0)[0][1] = (1.-z);
              data.second_derivative(k,1)[0][1] = -(1.-z);
              data.second_derivative(k,2)[0][1] = -(1.-z);
              data.second_derivative(k,3)[0][1] = (1.-z);
              data.second_derivative(k,4)[0][1] = z;
              data.second_derivative(k,5)[0][1] = -z;
              data.second_derivative(k,6)[0][1] = -z;
              data.second_derivative(k,7)[0][1] = z;
              data.second_derivative(k,0)[1][0] = (1.-z);
              data.second_derivative(k,1)[1][0] = -(1.-z);
              data.second_derivative(k,2)[1][0] = -(1.-z);
              data.second_derivative(k,3)[1][0] = (1.-z);
              data.second_derivative(k,4)[1][0] = z;
              data.second_derivative(k,5)[1][0] = -z;
              data.second_derivative(k,6)[1][0] = -z;
              data.second_derivative(k,7)[1][0] = z;

              data.second_derivative(k,0)[0][2] = (1.-y);
              data.second_derivative(k,1)[0][2] = -(1.-y);
              data.second_derivative(k,2)[0][2] = y;
              data.second_derivative(k,3)[0][2] = -y;
              data.second_derivative(k,4)[0][2] = -(1.-y);
              data.second_derivative(k,5)[0][2] = (1.-y);
              data.second_derivative(k,6)[0][2] = -y;
              data.second_derivative(k,7)[0][2] = y;
              data.second_derivative(k,0)[2][0] = (1.-y);
              data.second_derivative(k,1)[2][0] = -(1.-y);
              data.second_derivative(k,2)[2][0] = y;
              data.second_derivative(k,3)[2][0] = -y;
              data.second_derivative(k,4)[2][0] = -(1.-y);
              data.second_derivative(k,5)[2][0] = (1.-y);
              data.second_derivative(k,6)[2][0] = -y;
              data.second_derivative(k,7)[2][0] = y;

              data.second_derivative(k,0)[1][2] = (1.-x);
              data.second_derivative(k,1)[1][2] = x;
              data.second_derivative(k,2)[1][2] = -(1.-x);
              data.second_derivative(k,3)[1][2] = -x;
              data.second_derivative(k,4)[1][2] = -(1.-x);
              data.second_derivative(k,5)[1][2] = -x;
              data.second_derivative(k,6)[1][2] = (1.-x);
              data.second_derivative(k,7)[1][2] = x;
              data.second_derivative(k,0)[2][1] = (1.-x);
              data.second_derivative(k,1)[2][1] = x;
              data.second_derivative(k,2)[2][1] = -(1.-x);
              data.second_derivative(k,3)[2][1] = -x;
              data.second_derivative(k,4)[2][1] = -(1.-x);
              data.second_derivative(k,5)[2][1] = -x;
              data.second_derivative(k,6)[2][1] = (1.-x);
              data.second_derivative(k,7)[2][1] = x;
            }
        }
    }
  }
}


template<int dim, int spacedim>
void
MappingQ1<dim, spacedim>::
compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                        InternalData &data) const
{
  internal::MappingQ1::
  compute_shapes_virtual<spacedim> (n_shape_functions,
                                    unit_points, data);
}



template<int dim, int spacedim>
UpdateFlags
MappingQ1<dim,spacedim>::requires_update_flags (const UpdateFlags in) const
{
  // add flags if the respective quantities are necessary to compute
  // what we need. note that some flags appear in both the conditions
  // and in subsequent set operations. this leads to some circular
  // logic. the only way to treat this is to iterate. since there are
  // 5 if-clauses in the loop, it will take at most 5 iterations to
  // converge. do them:
  UpdateFlags out = in;
  for (unsigned int i=0; i<5; ++i)
    {
      // The following is a little incorrect:
      // If not applied on a face,
      // update_boundary_forms does not
      // make sense. On the other hand,
      // it is necessary on a
      // face. Currently,
      // update_boundary_forms is simply
      // ignored for the interior of a
      // cell.
      if (out & (update_JxW_values
                 | update_normal_vectors))
        out |= update_boundary_forms;

      if (out & (update_covariant_transformation
                 | update_JxW_values
                 | update_jacobians
                 | update_jacobian_grads
                 | update_boundary_forms
                 | update_normal_vectors))
        out |= update_contravariant_transformation;

      if (out & (update_inverse_jacobians))
        out |= update_covariant_transformation;

      // The contravariant transformation
      // is a Piola transformation, which
      // requires the determinant of the
      // Jacobi matrix of the transformation.
      // Therefore these values have to be
      // updated for each cell.
      if (out & update_contravariant_transformation)
        out |= update_JxW_values;

      if (out & update_normal_vectors)
        out |= update_JxW_values;
    }

  return out;
}


template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_data (const UpdateFlags      update_flags,
                                       const Quadrature<dim> &q,
                                       const unsigned int     n_original_q_points,
                                       InternalData          &data) const
{
  // store the flags in the internal data object so we can access them
  // in fill_fe_*_values(). use the transitive hull of the required
  // flags
  data.update_each = requires_update_flags(update_flags);

  const unsigned int n_q_points = q.size();

  // see if we need the (transformation) shape function values
  // and/or gradients and resize the necessary arrays
  if (data.update_each & update_quadrature_points)
    data.shape_values.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & (update_covariant_transformation
                          | update_contravariant_transformation
                          | update_JxW_values
                          | update_boundary_forms
                          | update_normal_vectors
                          | update_jacobians
                          | update_jacobian_grads
                          | update_inverse_jacobians))
    data.shape_derivatives.resize(data.n_shape_functions * n_q_points);

  if (data.update_each & update_covariant_transformation)
    data.covariant.resize(n_original_q_points);

  if (data.update_each & update_contravariant_transformation)
    data.contravariant.resize(n_original_q_points);

  if (data.update_each & update_volume_elements)
    data.volume_elements.resize(n_original_q_points);

  if (data.update_each & update_jacobian_grads)
    data.shape_second_derivatives.resize(data.n_shape_functions * n_q_points);

  compute_shapes (q.get_points(), data);
}



template<int dim, int spacedim>
typename MappingQ1<dim,spacedim>::InternalData *
MappingQ1<dim,spacedim>::get_data (const UpdateFlags update_flags,
                                   const Quadrature<dim> &q) const
{
  InternalData *data = new InternalData(n_shape_functions);
  compute_data (update_flags, q, q.size(), *data);
  return data;
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_face_data (const UpdateFlags update_flags,
                                            const Quadrature<dim> &q,
                                            const unsigned int n_original_q_points,
                                            InternalData &data) const
{
  compute_data (update_flags, q, n_original_q_points, data);

  if (dim > 1)
    {
      if (data.update_each & update_boundary_forms)
        {
          data.aux.resize (dim-1, std::vector<Tensor<1,spacedim> > (n_original_q_points));

          // Compute tangentials to the
          // unit cell.
          const unsigned int nfaces = GeometryInfo<dim>::faces_per_cell;
          data.unit_tangentials.resize (nfaces*(dim-1),
                                        std::vector<Tensor<1,dim> > (n_original_q_points));
          if (dim==2)
            {
              // ensure a counterclockwise
              // orientation of tangentials
              static const int tangential_orientation[4]= {-1,1,1,-1};
              for (unsigned int i=0; i<nfaces; ++i)
                {
                  Tensor<1,dim> tang;
                  tang[1-i/2]=tangential_orientation[i];
                  std::fill (data.unit_tangentials[i].begin(),
                             data.unit_tangentials[i].end(), tang);
                }
            }
          else if (dim==3)
            {
              for (unsigned int i=0; i<nfaces; ++i)
                {
                  Tensor<1,dim> tang1, tang2;

                  const unsigned int nd=
                    GeometryInfo<dim>::unit_normal_direction[i];

                  // first tangential
                  // vector in direction
                  // of the (nd+1)%3 axis
                  // and inverted in case
                  // of unit inward normal
                  tang1[(nd+1)%dim]=GeometryInfo<dim>::unit_normal_orientation[i];
                  // second tangential
                  // vector in direction
                  // of the (nd+2)%3 axis
                  tang2[(nd+2)%dim]=1.;

                  // same unit tangents
                  // for all quadrature
                  // points on this face
                  std::fill (data.unit_tangentials[i].begin(),
                             data.unit_tangentials[i].end(), tang1);
                  std::fill (data.unit_tangentials[nfaces+i].begin(),
                             data.unit_tangentials[nfaces+i].end(), tang2);
                }
            }
        }
    }
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ1<dim,spacedim>::get_face_data (const UpdateFlags        update_flags,
                                        const Quadrature<dim-1> &quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  compute_face_data (update_flags,
                     QProjector<dim>::project_to_all_faces(quadrature),
                     quadrature.size(),
                     *data);

  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ1<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
                                           const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  compute_face_data (update_flags,
                     QProjector<dim>::project_to_all_subfaces(quadrature),
                     quadrature.size(),
                     *data);

  return data;
}



namespace internal
{
  namespace
  {
    /**
     * Compute the locations of quadrature points on the object described by
     * the first argument (and the cell for which the mapping support points
     * have already been set), but only if the update_flags of the @p data
     * argument indicate so.
     */
    template <int dim, int spacedim>
    void
    maybe_compute_q_points (const typename dealii::MappingQ1<dim,spacedim>::DataSetDescriptor  data_set,
                            const typename dealii::MappingQ1<dim,spacedim>::InternalData      &data,
                            std::vector<Point<spacedim> >                                     &quadrature_points)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_quadrature_points)
        {
          for (unsigned int point=0; point<quadrature_points.size(); ++point)
            {
              const double *shape = &data.shape(point+data_set,0);
              Point<spacedim> result = (shape[0] *
                                        data.mapping_support_points[0]);
              for (unsigned int k=1; k<data.n_shape_functions; ++k)
                for (unsigned int i=0; i<spacedim; ++i)
                  result[i] += shape[k] * data.mapping_support_points[k][i];
              quadrature_points[point] = result;
            }
        }
    }


    /**
     * Update the co- and contravariant matrices as well as their determinant, for the cell
     * described stored in the data object, but only if the update_flags of the @p data
     * argument indicate so.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_Jacobians (const CellSimilarity::Similarity                                   cell_similarity,
                            const typename dealii::MappingQ1<dim,spacedim>::DataSetDescriptor  data_set,
                            const typename dealii::MappingQ1<dim,spacedim>::InternalData      &data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_contravariant_transformation)
        // if the current cell is just a
        // translation of the previous one, no
        // need to recompute jacobians...
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();

            std::fill(data.contravariant.begin(), data.contravariant.end(),
                      DerivativeForm<1,dim,spacedim>());

            Assert (data.n_shape_functions > 0, ExcInternalError());
            const Tensor<1,spacedim> *supp_pts =
              &data.mapping_support_points[0];

            for (unsigned int point=0; point<n_q_points; ++point)
              {
                const Tensor<1,dim> *data_derv =
                  &data.derivative(point+data_set, 0);

                double result [spacedim][dim];

                // peel away part of sum to avoid zeroing the
                // entries and adding for the first time
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    result[i][j] = data_derv[0][j] * supp_pts[0][i];
                for (unsigned int k=1; k<data.n_shape_functions; ++k)
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      result[i][j] += data_derv[k][j] * supp_pts[k][i];

                // write result into contravariant data. for
                // j=dim in the case dim<spacedim, there will
                // never be any nonzero data that arrives in
                // here, so it is ok anyway because it was
                // initialized to zero at the initialization
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    data.contravariant[point][i][j] = result[i][j];
              }
          }

      if (update_flags & update_covariant_transformation)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point=0; point<n_q_points; ++point)
              {
                data.covariant[point] = (data.contravariant[point]).covariant_form();
              }
          }

      if (update_flags & update_volume_elements)
        if (cell_similarity != CellSimilarity::translation)
          {
            const unsigned int n_q_points = data.contravariant.size();
            for (unsigned int point=0; point<n_q_points; ++point)
              data.volume_elements[point] = data.contravariant[point].determinant();
          }

    }

    /**
     * Update the Hessian of the transformation from unit to real cell, the
     * Jacobian gradients.
     *
     * Skip the computation if possible as indicated by the first argument.
     */
    template <int dim, int spacedim>
    void
    maybe_update_jacobian_grads (const CellSimilarity::Similarity                                   cell_similarity,
                                 const typename dealii::MappingQ1<dim,spacedim>::DataSetDescriptor  data_set,
                                 const typename dealii::MappingQ1<dim,spacedim>::InternalData      &data,
                                 std::vector<DerivativeForm<2,dim,spacedim> >                      &jacobian_grads)
    {
      const UpdateFlags update_flags = data.update_each;
      if (update_flags & update_jacobian_grads)
        {
          const unsigned int n_q_points = jacobian_grads.size();

          if (cell_similarity != CellSimilarity::translation)
            {
              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  const Tensor<2,dim> *second =
                    &data.second_derivative(point+data_set, 0);
                  double result [spacedim][dim][dim];
                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        result[i][j][l] = (second[0][j][l] *
                                           data.mapping_support_points[0][i]);
                  for (unsigned int k=1; k<data.n_shape_functions; ++k)
                    for (unsigned int i=0; i<spacedim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int l=0; l<dim; ++l)
                          result[i][j][l]
                          += (second[k][j][l]
                              *
                              data.mapping_support_points[k][i]);

                  for (unsigned int i=0; i<spacedim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int l=0; l<dim; ++l)
                        jacobian_grads[point][i][j][l] = result[i][j][l];
                }
            }
        }
    }
  }
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  std::vector<Point<spacedim> > &a) const
{
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  vertices = this->get_vertices(cell);
  a.resize(GeometryInfo<dim>::vertices_per_cell);

  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = vertices[i];
}



template<int dim, int spacedim>
CellSimilarity::Similarity
MappingQ1<dim,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           cell_similarity,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following static_cast is really correct:
  Assert (dynamic_cast<const InternalData *>(&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(internal_data);

  const unsigned int n_q_points=quadrature.size();

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      compute_mapping_support_points(cell, data.mapping_support_points);
      data.cell_of_current_support_points = cell;
    }

  internal::maybe_compute_q_points<dim,spacedim> (DataSetDescriptor::cell (),
                                                  data,
                                                  output_data.quadrature_points);
  internal::maybe_update_Jacobians<dim,spacedim> (cell_similarity,
                                                  DataSetDescriptor::cell (),
                                                  data);

  const UpdateFlags update_flags = data.update_each;
  const std::vector<double> &weights=quadrature.get_weights();

  // Multiply quadrature weights by absolute value of Jacobian determinants or
  // the area element g=sqrt(DX^t DX) in case of codim > 0

  if (update_flags & (update_normal_vectors
                      | update_JxW_values))
    {
      AssertDimension (output_data.JxW_values.size(), n_q_points);

      Assert( !(update_flags & update_normal_vectors ) ||
              (output_data.normal_vectors.size() == n_q_points),
              ExcDimensionMismatch(output_data.normal_vectors.size(), n_q_points));


      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          {

            if (dim == spacedim)
              {
                const double det = data.contravariant[point].determinant();

                // check for distorted cells.

                // TODO: this allows for anisotropies of up to 1e6 in 3D and
                // 1e12 in 2D. might want to find a finer
                // (dimension-independent) criterion
                Assert (det > 1e-12*Utilities::fixed_power<dim>(cell->diameter()/
                                                                std::sqrt(double(dim))),
                        (typename Mapping<dim,spacedim>::ExcDistortedMappedCell(cell->center(), det, point)));

                output_data.JxW_values[point] = weights[point] * det;
              }
            // if dim==spacedim, then there is no cell normal to
            // compute. since this is for FEValues (and not FEFaceValues),
            // there are also no face normals to compute
            else //codim>0 case
              {
                Tensor<1, spacedim> DX_t [dim];
                for (unsigned int i=0; i<spacedim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    DX_t[j][i] = data.contravariant[point][i][j];

                Tensor<2, dim> G; //First fundamental form
                for (unsigned int i=0; i<dim; ++i)
                  for (unsigned int j=0; j<dim; ++j)
                    G[i][j] = DX_t[i] * DX_t[j];

                output_data.JxW_values[point]
                  = sqrt(determinant(G)) * weights[point];

                if (cell_similarity == CellSimilarity::inverted_translation)
                  {
                    // we only need to flip the normal
                    if (update_flags & update_normal_vectors)
                      output_data.normal_vectors[point] *= -1.;
                  }
                else
                  {
                    const unsigned int codim = spacedim-dim;
                    (void)codim;

                    if (update_flags & update_normal_vectors)
                      {
                        Assert( codim==1 , ExcMessage("There is no cell normal in codim 2."));

                        if (dim==1)
                          cross_product(output_data.normal_vectors[point],
                                        -DX_t[0]);
                        else //dim == 2
                          cross_product(output_data.normal_vectors[point],DX_t[0],DX_t[1]);

                        output_data.normal_vectors[point] /= output_data.normal_vectors[point].norm();

                        if (cell->direction_flag() == false)
                          output_data.normal_vectors[point] *= -1.;
                      }

                  }
              } //codim>0 case

          }
    }



  // copy values from InternalData to vector given by reference
  if (update_flags & update_jacobians)
    {
      AssertDimension (output_data.jacobians.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          output_data.jacobians[point] = data.contravariant[point];
    }

  // copy values from InternalData to vector given by reference
  if (update_flags & update_inverse_jacobians)
    {
      AssertDimension (output_data.inverse_jacobians.size(), n_q_points);
      if (cell_similarity != CellSimilarity::translation)
        for (unsigned int point=0; point<n_q_points; ++point)
          output_data.inverse_jacobians[point] = data.covariant[point].transpose();
    }

  internal::maybe_update_jacobian_grads<dim,spacedim> (cell_similarity,
                                                       DataSetDescriptor::cell (),
                                                       data,
                                                       output_data.jacobian_grads);

  return cell_similarity;
}






namespace internal
{
  namespace
  {
    /**
     * Depending on what information is called for in the update flags of the
     * @p data object, compute the various pieces of information that is required
     * by the fill_fe_face_values() and fill_fe_subface_values() functions.
     * This function simply unifies the work that would be done by
     * those two functions.
     *
     * The resulting data is put into the @p output_data argument.
     */
    template <int dim, int spacedim>
    void
    maybe_compute_face_data (const dealii::MappingQ1<dim,spacedim> &mapping,
                             const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                             const unsigned int               face_no,
                             const unsigned int               subface_no,
                             const unsigned int               n_q_points,
                             const std::vector<double>        &weights,
                             const typename dealii::MappingQ1<dim,spacedim>::InternalData &data,
                             internal::FEValues::MappingRelatedData<dim,spacedim>         &output_data)
    {
      const UpdateFlags update_flags = data.update_each;

      if (update_flags & update_boundary_forms)
        {
          AssertDimension (output_data.boundary_forms.size(), n_q_points);
          if (update_flags & update_normal_vectors)
            AssertDimension (output_data.normal_vectors.size(), n_q_points);
          if (update_flags & update_JxW_values)
            AssertDimension (output_data.JxW_values.size(), n_q_points);

          // map the unit tangentials to the real cell. checking for d!=dim-1
          // eliminates compiler warnings regarding unsigned int expressions <
          // 0.
          for (unsigned int d=0; d!=dim-1; ++d)
            {
              Assert (face_no+GeometryInfo<dim>::faces_per_cell*d <
                      data.unit_tangentials.size(),
                      ExcInternalError());
              Assert (data.aux[d].size() <=
                      data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d].size(),
                      ExcInternalError());

              mapping.transform (data.unit_tangentials[face_no+GeometryInfo<dim>::faces_per_cell*d],
                                 data.aux[d],
                                 data,
                                 mapping_contravariant);
            }

          // if dim==spacedim, we can use the unit tangentials to compute the
          // boundary form by simply taking the cross product
          if (dim == spacedim)
            {
              for (unsigned int i=0; i<n_q_points; ++i)
                switch (dim)
                  {
                  case 1:
                    // in 1d, we don't have access to any of the data.aux
                    // fields (because it has only dim-1 components), but we
                    // can still compute the boundary form by simply
                    // looking at the number of the face
                    output_data.boundary_forms[i][0] = (face_no == 0 ?
                                                        -1 : +1);
                    break;
                  case 2:
                    cross_product (output_data.boundary_forms[i], data.aux[0][i]);
                    break;
                  case 3:
                    cross_product (output_data.boundary_forms[i], data.aux[0][i], data.aux[1][i]);
                    break;
                  default:
                    Assert(false, ExcNotImplemented());
                  }
            }
          else //(dim < spacedim)
            {
              // in the codim-one case, the boundary form results from the
              // cross product of all the face tangential vectors and the cell
              // normal vector
              //
              // to compute the cell normal, use the same method used in
              // fill_fe_values for cells above
              AssertDimension (data.contravariant.size(), n_q_points);

              for (unsigned int point=0; point<n_q_points; ++point)
                {
                  if (dim==1)
                    {
                      // J is a tangent vector
                      output_data.boundary_forms[point] = data.contravariant[point].transpose()[0];
                      output_data.boundary_forms[point] /=
                        (face_no == 0 ? -1. : +1.) * output_data.boundary_forms[point].norm();
                    }

                  if (dim==2)
                    {
                      Tensor<1,spacedim> cell_normal;
                      const DerivativeForm<1,spacedim,dim> DX_t =
                        data.contravariant[point].transpose();
                      cross_product(cell_normal,DX_t[0],DX_t[1]);
                      cell_normal /= cell_normal.norm();

                      // then compute the face normal from the face tangent
                      // and the cell normal:
                      cross_product (output_data.boundary_forms[point],
                                     data.aux[0][point], cell_normal);
                    }
                }
            }

          if (update_flags & (update_normal_vectors
                              | update_JxW_values))
            for (unsigned int i=0; i<output_data.boundary_forms.size(); ++i)
              {
                if (update_flags & update_JxW_values)
                  {
                    output_data.JxW_values[i] = output_data.boundary_forms[i].norm() * weights[i];

                    if (subface_no!=numbers::invalid_unsigned_int)
                      {
                        const double area_ratio=GeometryInfo<dim>::subface_ratio(
                                                  cell->subface_case(face_no), subface_no);
                        output_data.JxW_values[i] *= area_ratio;
                      }
                  }

                if (update_flags & update_normal_vectors)
                  output_data.normal_vectors[i] = Point<spacedim>(output_data.boundary_forms[i] /
                                                                  output_data.boundary_forms[i].norm());
              }

          if (update_flags & update_jacobians)
            for (unsigned int point=0; point<n_q_points; ++point)
              output_data.jacobians[point] = data.contravariant[point];

          if (update_flags & update_inverse_jacobians)
            for (unsigned int point=0; point<n_q_points; ++point)
              output_data.inverse_jacobians[point] = data.covariant[point].transpose();
        }
    }


    /**
     * Do the work of MappingQ1::fill_fe_face_values() and
     * MappingQ1::fill_fe_subface_values() in a generic way,
     * using the 'data_set' to differentiate whether we will
     * work on a face (and if so, which one) or subface.
     */
    template<int dim, int spacedim>
    void
    do_fill_fe_face_values (const dealii::MappingQ1<dim,spacedim>                             &mapping,
                            const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                            const unsigned int                                                 face_no,
                            const unsigned int                                                 subface_no,
                            const typename dealii::MappingQ1<dim,spacedim>::DataSetDescriptor  data_set,
                            const Quadrature<dim-1>                                           &quadrature,
                            const typename dealii::MappingQ1<dim,spacedim>::InternalData      &data,
                            internal::FEValues::MappingRelatedData<dim,spacedim>              &output_data)
    {
      maybe_compute_q_points<dim,spacedim> (data_set,
                                            data,
                                            output_data.quadrature_points);
      maybe_update_Jacobians<dim,spacedim> (CellSimilarity::none,
                                            data_set,
                                            data);
      maybe_update_jacobian_grads<dim,spacedim> (CellSimilarity::none,
                                                 data_set,
                                                 data,
                                                 output_data.jacobian_grads);
      maybe_compute_face_data (mapping,
                               cell, face_no, subface_no, quadrature.size(),
                               quadrature.get_weights(), data,
                               output_data);
    }
  }
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int                                         face_no,
                     const Quadrature<dim-1>                                   &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                     internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following cast is really correct:
  Assert ((dynamic_cast<const InternalData *>(&internal_data) != 0),
          ExcInternalError());
  const InternalData &data
    = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      compute_mapping_support_points(cell, data.mapping_support_points);
      data.cell_of_current_support_points = cell;
    }

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, numbers::invalid_unsigned_int,
                                    DataSetDescriptor::face (face_no,
                                                             cell->face_orientation(face_no),
                                                             cell->face_flip(face_no),
                                                             cell->face_rotation(face_no),
                                                             quadrature.size()),
                                    quadrature,
                                    data,
                                    output_data);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const unsigned int                                         face_no,
                        const unsigned int                                         subface_no,
                        const Quadrature<dim-1>                                   &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                        internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // ensure that the following cast is really correct:
  Assert ((dynamic_cast<const InternalData *>(&internal_data) != 0),
          ExcInternalError());
  const InternalData &data
    = static_cast<const InternalData &>(internal_data);

  // if necessary, recompute the support points of the transformation of this cell
  // (note that we need to first check the triangulation pointer, since otherwise
  // the second test might trigger an exception if the triangulations are not the
  // same)
  if ((data.mapping_support_points.size() == 0)
      ||
      (&cell->get_triangulation() !=
       &data.cell_of_current_support_points->get_triangulation())
      ||
      (cell != data.cell_of_current_support_points))
    {
      compute_mapping_support_points(cell, data.mapping_support_points);
      data.cell_of_current_support_points = cell;
    }

  internal::do_fill_fe_face_values (*this,
                                    cell, face_no, subface_no,
                                    DataSetDescriptor::subface (face_no, subface_no,
                                        cell->face_orientation(face_no),
                                        cell->face_flip(face_no),
                                        cell->face_rotation(face_no),
                                        quadrature.size(),
                                        cell->subface_case(face_no)),
                                    quadrature,
                                    data,
                                    output_data);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<1, dim> > >  input,
  VectorSlice<std::vector<Tensor<1, spacedim> > >         output,
  const typename Mapping<dim,spacedim>::InternalDataBase  &mapping_data,
  const MappingType                                       mapping_type) const
{
  transform_fields(input, output, mapping_data, mapping_type);
}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector<DerivativeForm<1, dim,spacedim> > > input,
  VectorSlice<std::vector<Tensor<2, spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  transform_differential_forms(input, output, mapping_data, mapping_type);

}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<2, dim> > >    input,
  VectorSlice<std::vector<Tensor<2, spacedim> > >     output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  switch (mapping_type)
    {
    case mapping_contravariant:
      transform_fields(input, output, mapping_data, mapping_type);
      return;

    case mapping_piola_gradient:
    case mapping_contravariant_gradient:
    case mapping_covariant_gradient:
      transform_gradients(input, output, mapping_data, mapping_type);
      return;
    default:
      Assert(false, ExcNotImplemented());
    }

}



template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector< DerivativeForm<2, dim, spacedim> > > input,
  VectorSlice<std::vector<Tensor<3,spacedim> > >             output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_covariant_gradient:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int q=0; q<output.size(); ++q)
        for (unsigned int i=0; i<spacedim; ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            for (unsigned int k=0; k<spacedim; ++k)
              {
                output[q][i][j][k] = data.covariant[q][j][0]
                                     * data.covariant[q][k][0]
                                     * input[q][i][0][0];
                for (unsigned int J=0; J<dim; ++J)
                  {
                    const unsigned int K0 = (0==J)? 1 : 0;
                    for (unsigned int K=K0; K<dim; ++K)
                      output[q][i][j][k] += data.covariant[q][j][J]
                                            * data.covariant[q][k][K]
                                            * input[q][i][J][K];
                  }

              }
    }
    default:
      Assert(false, ExcNotImplemented());
    }

}

template<int dim, int spacedim>
void
MappingQ1<dim,spacedim>::transform (
  const VectorSlice<const std::vector< Tensor<3,dim> > > input,
  VectorSlice<std::vector<Tensor<3,spacedim> > >             output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  switch (mapping_type)
    {
    case mapping_piola_hessian:
    case mapping_contravariant_hessian:
    case mapping_covariant_hessian:
      transform_hessians(input, output, mapping_data, mapping_type);
      return;
    default:
      Assert(false, ExcNotImplemented());
    }

}

template<int dim, int spacedim>
template < int rank >
void MappingQ1<dim,spacedim>::transform_fields(
  const VectorSlice<const std::vector<Tensor<rank,dim> > > input,
  VectorSlice<std::vector<Tensor<rank,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_contravariant:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.contravariant[i], input[i]);

      return;
    }

    case mapping_piola:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_each & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));
      Assert (rank==1, ExcMessage("Only for rank 1"));
      if (rank!=1)
        return;

      for (unsigned int i=0; i<output.size(); ++i)
        {
          output[i] = apply_transformation(data.contravariant[i], input[i]);
          output[i] /= data.volume_elements[i];
        }
      return;
    }
    //We still allow this operation as in the
    //reference cell Derivatives are Tensor
    //rather than DerivativeForm
    case mapping_covariant:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.covariant[i], input[i]);

      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}


template<int dim, int spacedim>
template < int rank >
void MappingQ1<dim,spacedim>::transform_gradients(
  const VectorSlice<const std::vector<Tensor<rank,dim> > > input,
  VectorSlice<std::vector<Tensor<rank,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_contravariant_gradient:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (rank==2, ExcMessage("Only for rank 2"));

      for (unsigned int i=0; i<output.size(); ++i)
        {
          DerivativeForm<1,spacedim,dim> A =
            apply_transformation(data.contravariant[i], transpose(input[i]) );
          output[i] = apply_transformation(data.covariant[i], A.transpose() );
        }

      return;
    }

    case mapping_covariant_gradient:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
      Assert (rank==2, ExcMessage("Only for rank 2"));

      for (unsigned int i=0; i<output.size(); ++i)
        {
          DerivativeForm<1,spacedim,dim> A =
            apply_transformation(data.covariant[i], transpose(input[i]) );
          output[i] = apply_transformation(data.covariant[i], A.transpose() );
        }

      return;
    }

    case mapping_piola_gradient:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_each & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));
      Assert (rank==2, ExcMessage("Only for rank 2"));

      for (unsigned int i=0; i<output.size(); ++i)
        {
          DerivativeForm<1,spacedim,dim> A =
            apply_transformation(data.covariant[i], input[i] );
          Tensor<2,spacedim> T =
            apply_transformation(data.contravariant[i], A.transpose() );

          output[i] = transpose(T);
          output[i] /= data.volume_elements[i];
        }

      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}




template<int dim, int spacedim>
void MappingQ1<dim,spacedim>::transform_hessians(
  const VectorSlice<const std::vector<Tensor<3,dim> > > input,
  VectorSlice<std::vector<Tensor<3,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_contravariant_hessian:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));

      for (unsigned int q=0; q<output.size(); ++q)
        for (unsigned int i=0; i<spacedim; ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            for (unsigned int k=0; k<spacedim; ++k)
              {
                output[q][i][j][k] =    data.contravariant[q][i][0]
                                        * data.covariant[q][j][0]
                                        * data.covariant[q][k][0]
                                        * input[q][0][0][0];
                for (unsigned int I=0; I<dim; ++I)
                  for (unsigned int J=0; J<dim; ++J)
                    {
                      const unsigned int K0 = (0==(I+J))? 1 : 0;
                      for (unsigned int K=K0; K<dim; ++K)
                        output[q][i][j][k] +=    data.contravariant[q][i][I]
                                                 * data.covariant[q][j][J]
                                                 * data.covariant[q][k][K]
                                                 * input[q][I][J][K];
                    }

              }
      return;
    }

    case mapping_covariant_hessian:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int q=0; q<output.size(); ++q)
        for (unsigned int i=0; i<spacedim; ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            for (unsigned int k=0; k<spacedim; ++k)
              {
                output[q][i][j][k] =    data.covariant[q][i][0]
                                        * data.covariant[q][j][0]
                                        * data.covariant[q][k][0]
                                        * input[q][0][0][0];
                for (unsigned int I=0; I<dim; ++I)
                  for (unsigned int J=0; J<dim; ++J)
                    {
                      const unsigned int K0 = (0==(I+J))? 1 : 0;
                      for (unsigned int K=K0; K<dim; ++K)
                        output[q][i][j][k] +=   data.covariant[q][i][I]
                                                * data.covariant[q][j][J]
                                                * data.covariant[q][k][K]
                                                * input[q][I][J][K];
                    }

              }

      return;
    }

    case mapping_piola_hessian:
    {
      Assert (data.update_each & update_covariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_contravariant_transformation"));
      Assert (data.update_each & update_volume_elements,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_volume_elements"));

      for (unsigned int q=0; q<output.size(); ++q)
        for (unsigned int i=0; i<spacedim; ++i)
          for (unsigned int j=0; j<spacedim; ++j)
            for (unsigned int k=0; k<spacedim; ++k)
              {
                output[q][i][j][k] =    data.contravariant[q][i][0]
                                        / data.volume_elements[q]
                                        * data.covariant[q][j][0]
                                        * data.covariant[q][k][0]
                                        * input[q][0][0][0];
                for (unsigned int I=0; I<dim; ++I)
                  for (unsigned int J=0; J<dim; ++J)
                    {
                      const unsigned int K0 = (0==(I+J))? 1 : 0;
                      for (unsigned int K=K0; K<dim; ++K)
                        output[q][i][j][k] +=    data.contravariant[q][i][I]
                                                 / data.volume_elements[q]
                                                 * data.covariant[q][j][J]
                                                 * data.covariant[q][k][K]
                                                 * input[q][I][J][K];
                    }

              }

      return;
    }

    default:
      Assert(false, ExcNotImplemented());
    }
}




template<int dim, int spacedim>
template < int rank >
void MappingQ1<dim,spacedim>::transform_differential_forms(
  const VectorSlice<const std::vector<DerivativeForm<rank, dim,spacedim> > >    input,
  VectorSlice<std::vector<Tensor<rank+1, spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{

  AssertDimension (input.size(), output.size());
  Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &>(mapping_data);

  switch (mapping_type)
    {
    case mapping_covariant:
    {
      Assert (data.update_each & update_contravariant_transformation,
              typename FEValuesBase<dim>::ExcAccessToUninitializedField("update_covariant_transformation"));

      for (unsigned int i=0; i<output.size(); ++i)
        output[i] = apply_transformation(data.covariant[i], input[i]);

      return;
    }
    default:
      Assert(false, ExcNotImplemented());
    }

}




template<int dim, int spacedim>
Point<spacedim>
MappingQ1<dim,spacedim>::transform_unit_to_real_cell (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Point<dim> &p) const
{
  // Use the get_data function to create an InternalData with data
  // vectors of the right size and transformation shape values already
  // computed at point p.
  const Quadrature<dim> point_quadrature(p);

  std_cxx11::unique_ptr<InternalData> mdata (get_data(update_quadrature_points | update_jacobians,
                                                      point_quadrature));

  // compute the mapping support
  // points
  compute_mapping_support_points(cell, mdata->mapping_support_points);

  // Mapping support points can be
  // bigger than necessary. If this
  // is the case, force them to be
  // Q1.
  if (mdata->mapping_support_points.size() > mdata->shape_values.size())
    mdata->mapping_support_points.resize
    (GeometryInfo<dim>::vertices_per_cell);


  return transform_unit_to_real_cell_internal(*mdata);
}



template<int dim, int spacedim>
Point<spacedim>
MappingQ1<dim,spacedim>::
transform_unit_to_real_cell_internal (const InternalData &data) const
{
  const unsigned int n_mapping_points=data.mapping_support_points.size();
  (void)n_mapping_points;
  AssertDimension (data.shape_values.size(), n_mapping_points);

  // use now the InternalData to
  // compute the point in real space.
  Point<spacedim> p_real;
  for (unsigned int i=0; i<data.mapping_support_points.size(); ++i)
    p_real += data.mapping_support_points[i] * data.shape(0,i);

  return p_real;
}



/* For an explanation of the  KA and Kb
   arrays see the comments in the declaration of
   transform_real_to_unit_cell_initial_guess */
namespace
{
  template <int dim>
  struct TransformR2UInitialGuess
  {
    static const double KA[GeometryInfo<dim>::vertices_per_cell][dim];
    static const double Kb[GeometryInfo<dim>::vertices_per_cell];
  };


  /*
    Octave code:
    M=[0 1; 1 1];
    K1 = transpose(M) * inverse (M*transpose(M));
    printf ("{%f, %f},\n", K1' );
  */
  template <>
  const double
  TransformR2UInitialGuess<1>::
  KA[GeometryInfo<1>::vertices_per_cell][1] =
  {
    {-1.000000},
    {1.000000}
  };

  template <>
  const double
  TransformR2UInitialGuess<1>::
  Kb[GeometryInfo<1>::vertices_per_cell] = {1.000000, 0.000000};


  /*
    Octave code:
    M=[0 1 0 1;0 0 1 1;1 1 1 1];
    K2 = transpose(M) * inverse (M*transpose(M));
    printf ("{%f, %f, %f},\n", K2' );
  */
  template <>
  const double
  TransformR2UInitialGuess<2>::
  KA[GeometryInfo<2>::vertices_per_cell][2] =
  {
    {-0.500000, -0.500000},
    { 0.500000, -0.500000},
    {-0.500000,  0.500000},
    { 0.500000,  0.500000}
  };

  /*
    Octave code:
    M=[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1; 1 1 1 1 1 1 1 1];
    K3 = transpose(M) * inverse (M*transpose(M))
    printf ("{%f, %f, %f, %f},\n", K3' );
  */
  template <>
  const double
  TransformR2UInitialGuess<2>::
  Kb[GeometryInfo<2>::vertices_per_cell] =
  {0.750000,0.250000,0.250000,-0.250000 };


  template <>
  const double
  TransformR2UInitialGuess<3>::
  KA[GeometryInfo<3>::vertices_per_cell][3] =
  {
    {-0.250000, -0.250000, -0.250000},
    { 0.250000, -0.250000, -0.250000},
    {-0.250000,  0.250000, -0.250000},
    { 0.250000,  0.250000, -0.250000},
    {-0.250000, -0.250000,  0.250000},
    { 0.250000, -0.250000,  0.250000},
    {-0.250000,  0.250000,  0.250000},
    { 0.250000,  0.250000,  0.250000}

  };


  template <>
  const double
  TransformR2UInitialGuess<3>::
  Kb[GeometryInfo<3>::vertices_per_cell] =
  {0.500000,0.250000,0.250000,0.000000,0.250000,0.000000,0.000000,-0.250000};

}


template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_initial_guess (const std::vector<Point<spacedim> > &vertex,
                                           const Point<spacedim>               &p) const
{
  Point<dim> p_unit;

  FullMatrix<double>  KA(GeometryInfo<dim>::vertices_per_cell, dim);
  Vector <double>  Kb(GeometryInfo<dim>::vertices_per_cell);

  KA.fill( (double *)(TransformR2UInitialGuess<dim>::KA) );
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    Kb(i)=(TransformR2UInitialGuess<dim>::Kb)[i];

  FullMatrix<double> Y(spacedim, GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
    for (unsigned int i=0; i<spacedim; ++i)
      Y(i,v) = vertex[v][i];

  FullMatrix<double> A(spacedim,dim);
  Y.mmult(A,KA); // A = Y*KA
  Vector< double > b(spacedim);
  Y.vmult(b,Kb); // b = Y*Kb

  for (unsigned int i=0; i<spacedim; ++i)
    b(i) -= p[i];
  b*=-1;

  Vector< double > dest(dim);

  FullMatrix<double> A_1(dim,spacedim);
  if (dim<spacedim)
    A_1.left_invert(A);
  else
    A_1.invert(A);

  A_1.vmult(dest,b); //A^{-1}*b

  for (unsigned int i=0; i<dim; ++i)
    p_unit[i]=dest(i);

  return p_unit;
}



template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
  // Use the exact formula if available
  if (dim == spacedim && (dim == 1 || dim == 2))
    {
      // The dimension-dependent algorithms are much faster (about 25-45x in 2D)
      // but fail most of the time when the given point (p) is not in the
      // cell. The dimension-independent Newton algorithm given below is
      // slower, but more robust (though it still sometimes fails).
      const std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
      vertices = this->get_vertices(cell);
      // These internal routines do not throw exceptions when the point does
      // not lie in the cell because manually checking (see the tolerances
      // before the return statements) is about 10% faster than a throw/catch.
      Point<dim> point = internal::MappingQ1::transform_real_to_unit_cell(vertices, p);

      if (dim == 1)
        {
          // formula not subject to any issues
          return point;
        }
      else if (dim == 2)
        {
          // formula not guaranteed to work for points outside of the cell
          const double eps = 1e-15;
          if (-eps <= point(1) && point(1) <= 1 + eps &&
              -eps <= point(0) && point(0) <= 1 + eps)
            {
              return point;
            }
        }
      else
        {
          Assert(false, ExcInternalError());
        }
    }
  // Find the initial value for the
  // Newton iteration by a normal projection
  // to the least square plane determined by
  // the vertices of the cell
  std::vector<Point<spacedim> > a;
  compute_mapping_support_points (cell,a);
  Point<dim> initial_p_unit =
    transform_real_to_unit_cell_initial_guess(a,p);

  // if dim==1 there is nothing
  // else to do to the initial
  // value, and it is the answer
  if (dim == 1)
    return initial_p_unit;
  else
    {
      // use the full mapping. in case the
      // function above should have given us
      // something back that lies outside the
      // unit cell (that might happen because
      // either the function computing an
      // initial guess gave us a poor initial
      // guess or for the following reason:
      // we call this function here in the Q1
      // mapping to produce an initial guess
      // for a higher order mapping, but
      // we may have given a point 'p' that
      // lies inside the cell with the higher
      // order mapping, but outside the
      // Q1-mapped reference cell), then
      // project it back into the reference
      // cell in hopes that this gives a
      // better starting point to the
      // following iteration
//TODO: the following line was added in r25581 but it leads to
// changes in the test results. investigate why this is so --
// it shouldn't really make any difference...
//      initial_p_unit = GeometryInfo<dim>::project_to_unit_cell(initial_p_unit);

      // Use the get_data function to
      // create an InternalData with data
      // vectors of the right size and
      // transformation shape values and
      // derivatives already computed at
      // point initial_p_unit.
      const Quadrature<dim> point_quadrature(initial_p_unit);

      UpdateFlags update_flags = update_quadrature_points | update_jacobians;
      if (spacedim>dim)
        update_flags |= update_jacobian_grads;

      std_cxx11::unique_ptr<InternalData> mdata(MappingQ1<dim,spacedim>::get_data(update_flags,
                                                point_quadrature));

      compute_mapping_support_points (cell, mdata->mapping_support_points);
      // The support points have to be at
      // least as many as there are
      // vertices.
      Assert(mdata->mapping_support_points.size() >=
             GeometryInfo<dim>::vertices_per_cell,
             ExcInternalError());
      // Ignore non vertex support points.
      mdata->mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);

      // perform the Newton iteration and
      // return the result. note that this
      // statement may throw an exception, which
      // we simply pass up to the caller
      return transform_real_to_unit_cell_internal(cell, p, initial_p_unit,
                                                  *mdata);
    }
}



template<int dim, int spacedim>
Point<dim>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_internal
(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
 const Point<spacedim>                            &p,
 const Point<dim>                                 &initial_p_unit,
 InternalData                                     &mdata) const
{
  const unsigned int n_shapes=mdata.shape_values.size();
  (void)n_shapes;
  Assert(n_shapes!=0, ExcInternalError());
  AssertDimension (mdata.shape_derivatives.size(), n_shapes);

  std::vector<Point<spacedim> > &points=mdata.mapping_support_points;
  AssertDimension (points.size(), n_shapes);


  // Newton iteration to solve
  //    f(x)=p(x)-p=0
  // where we are looking for 'x' and p(x) is the forward transformation
  // from unit to real cell. We solve this using a Newton iteration
  //    x_{n+1}=x_n-[f'(x)]^{-1}f(x)
  // The start value is set to be the linear approximation to the cell

  // The shape values and derivatives of the mapping at this point are
  // previously computed.

  Point<dim> p_unit = initial_p_unit;

  compute_shapes(std::vector<Point<dim> > (1, p_unit), mdata);
  Point<spacedim> p_real = transform_unit_to_real_cell_internal(mdata);
  Tensor<1,spacedim> f = p_real-p;

  // early out if we already have our point
  if (f.norm_square() < 1e-24 * cell->diameter() * cell->diameter())
    return p_unit;

  // we need to compare the position of the computed p(x) against the given
  // point 'p'. We will terminate the iteration and return 'x' if they are
  // less than eps apart. The question is how to choose eps -- or, put maybe
  // more generally: in which norm we want these 'p' and 'p(x)' to be eps
  // apart.
  //
  // the question is difficult since we may have to deal with very elongated
  // cells where we may achieve 1e-12*h for the distance of these two points
  // in the 'long' direction, but achieving this tolerance in the 'short'
  // direction of the cell may not be possible
  //
  // what we do instead is then to terminate iterations if
  //    \| p(x) - p \|_A < eps
  // where the A-norm is somehow induced by the transformation of the cell.
  // in particular, we want to measure distances relative to the sizes of
  // the cell in its principal directions.
  //
  // to define what exactly A should be, note that to first order we have
  // the following (assuming that x* is the solution of the problem, i.e.,
  // p(x*)=p):
  //    p(x) - p = p(x) - p(x*)
  //             = -grad p(x) * (x*-x) + higher order terms
  // This suggest to measure with a norm that corresponds to
  //    A = {[grad p(x]^T [grad p(x)]}^{-1}
  // because then
  //    \| p(x) - p \|_A  \approx  \| x - x* \|
  // Consequently, we will try to enforce that
  //    \| p(x) - p \|_A  =  \| f \|  <=  eps
  //
  // Note that using this norm is a bit dangerous since the norm changes
  // in every iteration (A isn't fixed by depends on xk). However, if the
  // cell is not too deformed (it may be stretched, but not twisted) then
  // the mapping is almost linear and A is indeed constant or nearly so.
  const double eps = 1.e-11;
  const unsigned int newton_iteration_limit = 20;

  unsigned int newton_iteration = 0;
  double last_f_weighted_norm;
  do
    {
#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
      std::cout << "Newton iteration " << newton_iteration << std::endl;
#endif

      // f'(x)
      Tensor<2,spacedim> df;
      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
        {
          const Tensor<1,dim> &grad_transform=mdata.derivative(0,k);
          const Point<spacedim> &point=points[k];

          for (unsigned int i=0; i<spacedim; ++i)
            for (unsigned int j=0; j<dim; ++j)
              df[i][j]+=point[i]*grad_transform[j];
        }

      // Solve  [f'(x)]d=f(x)
      Tensor<1,spacedim> delta;
      Tensor<2,spacedim> df_inverse = invert(df);
      contract (delta, df_inverse, static_cast<const Tensor<1,spacedim>&>(f));

#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
      std::cout << "   delta=" << delta  << std::endl;
#endif

      // do a line search
      double step_length = 1;
      do
        {
          // update of p_unit. The spacedim-th component of transformed point
          // is simply ignored in codimension one case. When this component is
          // not zero, then we are projecting the point to the surface or
          // curve identified by the cell.
          Point<dim> p_unit_trial = p_unit;
          for (unsigned int i=0; i<dim; ++i)
            p_unit_trial[i] -= step_length * delta[i];

          // shape values and derivatives
          // at new p_unit point
          compute_shapes(std::vector<Point<dim> > (1, p_unit_trial), mdata);

          // f(x)
          Point<spacedim> p_real_trial = transform_unit_to_real_cell_internal(mdata);
          const Tensor<1,spacedim> f_trial = p_real_trial-p;

#ifdef DEBUG_TRANSFORM_REAL_TO_UNIT_CELL
          std::cout << "     step_length=" << step_length << std::endl
                    << "       ||f ||   =" << f.norm() << std::endl
                    << "       ||f*||   =" << f_trial.norm() << std::endl
                    << "       ||f*||_A =" << (df_inverse * f_trial).norm() << std::endl;
#endif

          // see if we are making progress with the current step length
          // and if not, reduce it by a factor of two and try again
          //
          // strictly speaking, we should probably use the same norm as we use
          // for the outer algorithm. in practice, line search is just a
          // crutch to find a "reasonable" step length, and so using the l2
          // norm is probably just fine
          if (f_trial.norm() < f.norm())
            {
              p_real = p_real_trial;
              p_unit = p_unit_trial;
              f = f_trial;
              break;
            }
          else if (step_length > 0.05)
            step_length /= 2;
          else
            AssertThrow (false,
                         (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
        }
      while (true);

      ++newton_iteration;
      if (newton_iteration > newton_iteration_limit)
        AssertThrow (false,
                     (typename Mapping<dim,spacedim>::ExcTransformationFailed()));
      last_f_weighted_norm = (df_inverse * f).norm();
    }
  while (last_f_weighted_norm > eps);

  return p_unit;
}



/*
  This function becomes a little tricky in dimension <2,3>.
  There is a surface embedded in R^3 and we pass a point p in R^3, that
  is most likely not lying on the surface.
  We then ask,
  what point in R^2 (hopefully in the unit cell) satisfies that
  map(x) = p.

  An appropriate modification of this question is:
  Find x in R^2 and alpha in R such that

  map(x) + alpha * normal(x) = p


 */

template<>
Point<2>
MappingQ1<2,3>::
transform_real_to_unit_cell_internal (const Triangulation<2,3>::cell_iterator &cell,
                                      const Point<3> &p,
                                      const Point<2> &initial_p_unit,
                                      InternalData   &mdata) const
{
  return
    transform_real_to_unit_cell_internal_codim1(cell, p, initial_p_unit,
                                                mdata);
}




template<>
Point<1>
MappingQ1<1,2>::
transform_real_to_unit_cell_internal (const Triangulation<1,2>::cell_iterator &cell,
                                      const Point<2> &p,
                                      const Point<1> &initial_p_unit,
                                      InternalData   &mdata) const
{
  return
    transform_real_to_unit_cell_internal_codim1(cell, p, initial_p_unit,
                                                mdata);
}


template<>
Point<1>
MappingQ1<1,3>::
transform_real_to_unit_cell_internal (const Triangulation<1,3>::cell_iterator &/*cell*/,
                                      const Point<3> &/*p*/,
                                      const Point<1> &/*initial_p_unit*/,
                                      InternalData   &/*mdata*/) const
{
  Assert(false, ExcNotImplemented());
  return Point<1>();
}




template<int dim, int spacedim>
template<int dim_>
Point<dim_>
MappingQ1<dim,spacedim>::
transform_real_to_unit_cell_internal_codim1
(const typename Triangulation<dim_,dim_ + 1>::cell_iterator &cell,
 const Point<dim_ + 1>                                      &p,
 const Point<dim_ >                                         &initial_p_unit,
 typename MappingQ1<dim,spacedim>::InternalData                      &mdata) const
{
  const unsigned int spacedim1 = dim_+1;
  const unsigned int dim1 = dim_;


  const unsigned int n_shapes=mdata.shape_values.size();
  (void)n_shapes;
  Assert(n_shapes!=0, ExcInternalError());
  Assert(mdata.shape_derivatives.size()==n_shapes, ExcInternalError());
  Assert(mdata.shape_second_derivatives.size()==n_shapes, ExcInternalError());

  std::vector<Point<spacedim1> > &points=mdata.mapping_support_points;
  Assert(points.size()==n_shapes, ExcInternalError());

  Point<spacedim1> p_minus_F;

  Tensor<1,spacedim1>  DF[dim1];
  Tensor<1,spacedim1>  D2F[dim1][dim1];

  Point<dim1> p_unit = initial_p_unit;
  Point<dim1> f;
  Tensor<2,dim1>  df;

  //Evaluate first and second derivatives
  compute_shapes(std::vector<Point<dim1> > (1, p_unit), mdata);
  for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
    {
      const Tensor<1,dim1>   &grad_phi_k = mdata.derivative(0,k);
      const Tensor<2,dim1>   &hessian_k  = mdata.second_derivative(0,k);
      const Point<spacedim1> &point_k = points[k];

      for (unsigned int j=0; j<dim1; ++j)
        {
          DF[j] += grad_phi_k[j] * point_k;
          for (unsigned int l=0; l<dim1; ++l)
            D2F[j][l] += hessian_k[j][l] * point_k;
        }
    }

  p_minus_F = p;
  p_minus_F -= transform_unit_to_real_cell_internal(mdata);


  for (unsigned int j=0; j<dim1; ++j)
    f[j] = DF[j] * p_minus_F;

  for (unsigned int j=0; j<dim1; ++j)
    {
      f[j] = DF[j] * p_minus_F;
      for (unsigned int l=0; l<dim1; ++l)
        df[j][l] = -DF[j]*DF[l] + D2F[j][l] * p_minus_F;
    }


  const double eps = 1.e-12*cell->diameter();
  const unsigned int loop_limit = 10;

  unsigned int loop=0;

  while (f.norm()>eps && loop++<loop_limit)
    {
      // Solve  [df(x)]d=f(x)
      Tensor<1,dim1> d;
      Tensor<2,dim1> df_1;

      df_1 = invert(df);
      contract (d, df_1, static_cast<const Tensor<1,dim1>&>(f));
      p_unit -= d;

      for (unsigned int j=0; j<dim1; ++j)
        {
          DF[j].clear();
          for (unsigned int l=0; l<dim1; ++l)
            D2F[j][l].clear();
        }

      compute_shapes(std::vector<Point<dim1> > (1, p_unit), mdata);
      for (unsigned int k=0; k<mdata.n_shape_functions; ++k)
        {
          const Tensor<1,dim1>   &grad_phi_k = mdata.derivative(0,k);
          const Tensor<2,dim1>   &hessian_k  = mdata.second_derivative(0,k);
          const Point<spacedim1> &point_k = points[k];

          for (unsigned int j=0; j<dim1; ++j)
            {
              DF[j] += grad_phi_k[j] * point_k;
              for (unsigned int l=0; l<dim1; ++l)
                D2F[j][l] += hessian_k[j][l] * point_k;
            }
        }

//TODO: implement a line search here in much the same way as for
// the corresponding function above that does so for codim==0.
      p_minus_F = p;
      p_minus_F -= transform_unit_to_real_cell_internal(mdata);

      for (unsigned int j=0; j<dim1; ++j)
        {
          f[j] = DF[j] * p_minus_F;
          for (unsigned int l=0; l<dim1; ++l)
            df[j][l] = -DF[j]*DF[l] + D2F[j][l] * p_minus_F;
        }

    }


  // Here we check that in the last
  // execution of while the first
  // condition was already wrong,
  // meaning the residual was below
  // eps. Only if the first condition
  // failed, loop will have been
  // increased and tested, and thus
  // have reached the limit.
  AssertThrow (loop<loop_limit, (typename Mapping<dim,spacedim>::ExcTransformationFailed()));

  return p_unit;
}





template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ1<dim,spacedim>::clone () const
{
  return new MappingQ1<dim,spacedim>(*this);
}

//---------------------------------------------------------------------------


template<int dim, int spacedim>
MappingQ1<dim,spacedim> StaticMappingQ1<dim,spacedim>::mapping;



//--------------------------- Explicit instantiations -----------------------
#include "mapping_q1.inst"


DEAL_II_NAMESPACE_CLOSE
