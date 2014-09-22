// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <numeric>
#include <memory>

DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim>
MappingQ<dim,spacedim>::InternalData::InternalData (const unsigned int n_shape_functions)
  :
  MappingQ1<dim,spacedim>::InternalData(n_shape_functions),
  use_mapping_q1_on_current_cell(false),
  mapping_q1_data(1 << dim)
{
  this->is_mapping_q1_data=false;
}



template<int dim, int spacedim>
std::size_t
MappingQ<dim,spacedim>::InternalData::memory_consumption () const
{
  return (MappingQ1<dim,spacedim>::InternalData::memory_consumption () +
          MemoryConsumption::memory_consumption (unit_normals) +
          MemoryConsumption::memory_consumption (use_mapping_q1_on_current_cell) +
          MemoryConsumption::memory_consumption (mapping_q1_data));
}



// in 1d, it is irrelevant which polynomial degree to use, since all
// cells are scaled linearly. Unless codimension is equal to two
template<>
MappingQ<1>::MappingQ (const unsigned int,
                       const bool /*use_mapping_q_on_all_cells*/)
  :
  degree(1),
  n_inner(0),
  n_outer(0),
  tensor_pols(0),
  n_shape_functions(2),
  renumber(0),
  use_mapping_q_on_all_cells (false),
  feq(degree),
  line_support_points(degree+1)
{}


template<>
MappingQ<1>::MappingQ (const MappingQ<1> &m):
  MappingQ1<1> (),
  degree(1),
  n_inner(0),
  n_outer(0),
  tensor_pols(0),
  n_shape_functions(2),
  renumber(0),
  use_mapping_q_on_all_cells (m.use_mapping_q_on_all_cells),
  feq(degree),
  line_support_points(degree+1)
{}

template<>
MappingQ<1>::~MappingQ ()
{}



namespace
{
  template <int dim>
  std::vector<unsigned int>
  get_dpo_vector (const unsigned int degree)
  {
    std::vector<unsigned int> dpo(dim+1, 1U);
    for (unsigned int i=1; i<dpo.size(); ++i)
      dpo[i]=dpo[i-1]*(degree-1);
    return dpo;
  }
}




template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const unsigned int p,
                                  const bool use_mapping_q_on_all_cells)
  :
  degree(p),
  n_inner(Utilities::fixed_power<dim>(degree-1)),
  n_outer((dim==1) ? 2 :
          ((dim==2) ?
           4+4*(degree-1) :
           8+12*(degree-1)+6*(degree-1)*(degree-1))),
  tensor_pols(0),
  n_shape_functions(Utilities::fixed_power<dim>(degree+1)),
  renumber(FETools::
           lexicographic_to_hierarchic_numbering (
             FiniteElementData<dim> (get_dpo_vector<dim>(degree), 1,
                                     degree))),
  use_mapping_q_on_all_cells (use_mapping_q_on_all_cells
                              || (dim != spacedim)),
  feq(degree),
  line_support_points(degree+1)
{
  // Construct the tensor product polynomials used as shape functions for the
  // Qp mapping of cells at the boundary.
  tensor_pols = new TensorProductPolynomials<dim>
  (Polynomials::generate_complete_Lagrange_basis(line_support_points.get_points()));
  Assert (n_shape_functions==tensor_pols->n(),
          ExcInternalError());
  Assert(n_inner+n_outer==n_shape_functions, ExcInternalError());

  // build laplace_on_quad_vector
  if (degree>1)
    {
      if (dim >= 2)
        set_laplace_on_quad_vector(laplace_on_quad_vector);
      if (dim >= 3)
        set_laplace_on_hex_vector(laplace_on_hex_vector);
    }
}


template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const MappingQ<dim,spacedim> &mapping)
  :
  MappingQ1<dim,spacedim>(),
  degree(mapping.degree),
  n_inner(mapping.n_inner),
  n_outer(mapping.n_outer),
  tensor_pols(0),
  n_shape_functions(mapping.n_shape_functions),
  renumber(mapping.renumber),
  use_mapping_q_on_all_cells (mapping.use_mapping_q_on_all_cells),
  feq(degree),
  line_support_points(degree+1)
{
  tensor_pols=new TensorProductPolynomials<dim> (*mapping.tensor_pols);
  laplace_on_quad_vector=mapping.laplace_on_quad_vector;
  laplace_on_hex_vector=mapping.laplace_on_hex_vector;
}


template<int dim, int spacedim>
MappingQ<dim,spacedim>::~MappingQ ()
{
  delete tensor_pols;
}



template<>
void
MappingQ<1>::compute_shapes_virtual (const std::vector<Point<1> > &unit_points,
                                     MappingQ1<1>::InternalData   &data) const
{
  MappingQ1<1>::compute_shapes_virtual(unit_points, data);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
                                                typename MappingQ1<dim,spacedim>::InternalData &data) const
{
  const unsigned int n_points=unit_points.size();
  std::vector<double> values;
  std::vector<Tensor<1,dim> > grads;
  if (data.shape_values.size()!=0)
    {
      Assert(data.shape_values.size()==n_shape_functions*n_points,
             ExcInternalError());
      values.resize(n_shape_functions);
    }
  if (data.shape_derivatives.size()!=0)
    {
      Assert(data.shape_derivatives.size()==n_shape_functions*n_points,
             ExcInternalError());
      grads.resize(n_shape_functions);
    }

//                                 // dummy variable of size 0
  std::vector<Tensor<2,dim> > grad2;
  if (data.shape_second_derivatives.size()!=0)
    {
      Assert(data.shape_second_derivatives.size()==n_shape_functions*n_points,
             ExcInternalError());
      grad2.resize(n_shape_functions);
    }


  if (data.shape_values.size()!=0 || data.shape_derivatives.size()!=0)
    for (unsigned int point=0; point<n_points; ++point)
      {
        tensor_pols->compute(unit_points[point], values, grads, grad2);

        if (data.shape_values.size()!=0)
          for (unsigned int i=0; i<n_shape_functions; ++i)
            data.shape(point,renumber[i]) = values[i];

        if (data.shape_derivatives.size()!=0)
          for (unsigned int i=0; i<n_shape_functions; ++i)
            data.derivative(point,renumber[i]) = grads[i];

        if (data.shape_second_derivatives.size()!=0)
          for (unsigned int i=0; i<n_shape_functions; ++i)
            data.second_derivative(point,renumber[i]) = grad2[i];
      }
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_data (const UpdateFlags update_flags,
                                  const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  this->compute_data (update_flags, quadrature,
                      quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_data (update_flags, quadrature,
                        quadrature.size(), data->mapping_q1_data);
  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_face_data (const UpdateFlags update_flags,
                                       const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.size(),
                             data->mapping_q1_data);
  return data;
}



template<int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
MappingQ<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
                                          const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(n_shape_functions);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));
  this->compute_face_data (update_flags, q,
                           quadrature.size(), *data);
  if (!use_mapping_q_on_all_cells)
    this->compute_face_data (update_flags, q,
                             quadrature.size(),
                             data->mapping_q1_data);
  return data;
}


// Note that the CellSimilarity flag is modifyable, since MappingQ can need to
// recalculate data even when cells are similar.
template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const Quadrature<dim>                                     &q,
  typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
  std::vector<Point<spacedim> >                             &quadrature_points,
  std::vector<double>                                       &JxW_values,
  std::vector<DerivativeForm<1,dim,spacedim>   >    &jacobians,
  std::vector<DerivativeForm<2,dim,spacedim>    >   &jacobian_grads,
  std::vector<DerivativeForm<1,spacedim,dim>    >   &inverse_jacobians,
  std::vector<Point<spacedim> >                             &normal_vectors,
  CellSimilarity::Similarity                           &cell_similarity) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0, ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is in the interior of the domain
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());

  // depending on this result, use this or the other data object for the
  // mapping. furthermore, we need to ensure that the flag indicating whether
  // we can use some similarity has to be modified - for a general MappingQ,
  // the data needs to be recomputed anyway since then the mapping changes the
  // data. this needs to be known also for later operations, so modify the
  // variable here. this also affects the calculation of the next cell -- if
  // we use Q1 data on the next cell, the data will still be invalid.
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    {
      p_data=&data;
      if (get_degree() > 1)
        cell_similarity = CellSimilarity::invalid_next_cell;
    }

  MappingQ1<dim,spacedim>::fill_fe_values(cell, q, *p_data,
                                          quadrature_points, JxW_values,
                                          jacobians, jacobian_grads, inverse_jacobians,
                                          normal_vectors, cell_similarity);

}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_face_values (
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  const unsigned int       face_no,
  const Quadrature<dim-1> &q,
  typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  std::vector<Point<spacedim> >     &quadrature_points,
  std::vector<double>          &JxW_values,
  std::vector<Tensor<1,spacedim> >  &exterior_forms,
  std::vector<Point<spacedim> >     &normal_vectors) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0,
          ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
                                        || cell->has_boundary_lines());

  // depending on this result, use this or the other data object for the
  // mapping
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.size();
  this->compute_fill_face (cell, face_no, deal_II_numbers::invalid_unsigned_int,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           face (face_no,
                                 cell->face_orientation(face_no),
                                 cell->face_flip(face_no),
                                 cell->face_rotation(face_no),
                                 n_q_points),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors);
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                                const unsigned int       face_no,
                                                const unsigned int       sub_no,
                                                const Quadrature<dim-1> &q,
                                                typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
                                                std::vector<Point<spacedim> >     &quadrature_points,
                                                std::vector<double>          &JxW_values,
                                                std::vector<Tensor<1,spacedim> >  &exterior_forms,
                                                std::vector<Point<spacedim> >     &normal_vectors) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<InternalData *> (&mapping_data) != 0,
          ExcInternalError());
  InternalData &data = static_cast<InternalData &> (mapping_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell=!(use_mapping_q_on_all_cells
                                        || cell->has_boundary_lines());

  // depending on this result, use this or the other data object for the
  // mapping
  typename MappingQ1<dim,spacedim>::InternalData *p_data=0;
  if (data.use_mapping_q1_on_current_cell)
    p_data=&data.mapping_q1_data;
  else
    p_data=&data;

  const unsigned int n_q_points=q.size();
  this->compute_fill_face (cell, face_no, sub_no,
                           n_q_points,
                           QProjector<dim>::DataSetDescriptor::
                           subface (face_no, sub_no,
                                    cell->face_orientation(face_no),
                                    cell->face_flip(face_no),
                                    cell->face_rotation(face_no),
                                    n_q_points,
                                    cell->subface_case(face_no)),
                           q.get_weights(),
                           *p_data,
                           quadrature_points, JxW_values,
                           exterior_forms, normal_vectors);
}



template <>
void
MappingQ<1>::set_laplace_on_quad_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::set_laplace_on_quad_vector(Table<2,double> &loqvs) const
{
  Assert(degree>1, ExcInternalError());
  const unsigned int n_inner_2d=(degree-1)*(degree-1);
  const unsigned int n_outer_2d=4+4*(degree-1);

  // first check whether we have precomputed the values for some polynomial
  // degree; the sizes of arrays is n_inner_2d*n_outer_2d
  double const *loqv_ptr=0;
  switch (degree)
    {
    // for degree==1, we shouldn't have to compute any support points, since
    // all of them are on the vertices

    case 2:
    {
      // (checked these values against the output of compute_laplace_vector
      // again, and found they're indeed right -- just in case someone wonders
      // where they come from -- WB)
      static const double loqv2[1*8]
        = {1/16., 1/16., 1/16., 1/16., 3/16., 3/16., 3/16., 3/16.};
      loqv_ptr=&loqv2[0];
      Assert (sizeof(loqv2)/sizeof(loqv2[0]) ==
              n_inner_2d * n_outer_2d,
              ExcInternalError());

      break;
    }

    // no other cases implemented, so simply fall through
    default:
      break;
    }

  if (loqv_ptr!=0)
    {
      // precomputed. copy values to the loqvs array
      loqvs.reinit(n_inner_2d, n_outer_2d);
      for (unsigned int unit_point=0; unit_point<n_inner_2d; ++unit_point)
        for (unsigned int k=0; k<n_outer_2d; ++k)
          loqvs[unit_point][k]=loqv_ptr[unit_point*n_outer_2d+k];
    }
  else
    {
      // not precomputed, then do so now
      if (dim == 2)
        compute_laplace_vector(loqvs);
      else if (dim == 3)
        {
          MappingQ<2,2> mapping_2d(this->degree);
          loqvs = mapping_2d.laplace_on_quad_vector;
        }
    }

  // the sum of weights of the points at the outer rim should be one. check
  // this
  for (unsigned int unit_point=0; unit_point<loqvs.n_rows(); ++unit_point)
    Assert(std::fabs(std::accumulate(loqvs[unit_point].begin(),
                                     loqvs[unit_point].end(),0.)-1)<1e-13*this->degree,
           ExcInternalError());
}



template <>
void
MappingQ<3>::set_laplace_on_hex_vector(Table<2,double> &lohvs) const
{
  Assert(degree>1, ExcInternalError());

  // first check whether we have precomputed the values for some polynomial
  // degree
  double const *lohv_ptr=0;
  if (degree==2)
    {
      static const double loqv2[26]
        = {1/128., 1/128., 1/128., 1/128., 1/128., 1/128., 1/128., 1/128.,
           7/192., 7/192., 7/192., 7/192., 7/192., 7/192., 7/192., 7/192.,
           7/192., 7/192., 7/192., 7/192.,
           1/12., 1/12., 1/12., 1/12., 1/12., 1/12.
          };

      lohv_ptr=&loqv2[0];
    }

  if (lohv_ptr!=0)
    {
      // precomputed. copy values to the lohvs array
      lohvs.reinit(n_inner, n_outer);
      for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
        for (unsigned int k=0; k<n_outer; ++k)
          lohvs[unit_point][k]=lohv_ptr[unit_point*n_outer+k];
    }
  else
    // not precomputed, then do so now
    compute_laplace_vector(lohvs);

  // the sum of weights of the points at the outer rim should be one. check
  // this
  for (unsigned int unit_point=0; unit_point<n_inner; ++unit_point)
    Assert(std::fabs(std::accumulate(lohvs[unit_point].begin(),
                                     lohvs[unit_point].end(),0.) - 1)<1e-13*this->degree,
           ExcInternalError());
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::set_laplace_on_hex_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}



template <>
void
MappingQ<1>::compute_laplace_vector(Table<2,double> &) const
{
  Assert(false, ExcInternalError());
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_laplace_vector(Table<2,double> &lvs) const
{
  Assert(lvs.n_rows()==0, ExcInternalError());
  Assert(dim==2 || dim==3, ExcNotImplemented());

  // for degree==1, we shouldn't have to compute any support points, since all
  // of them are on the vertices
  Assert(degree>1, ExcInternalError());

  // compute the shape gradients at the quadrature points on the unit cell
  const QGauss<dim> quadrature(degree+1);
  const unsigned int n_q_points=quadrature.size();

  InternalData quadrature_data(n_shape_functions);
  quadrature_data.shape_derivatives.resize(n_shape_functions * n_q_points);
  this->compute_shapes(quadrature.get_points(), quadrature_data);

  // Compute the stiffness matrix of the inner dofs
  FullMatrix<long double> S(n_inner);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int j=0; j<n_inner; ++j)
        {
          long double res = 0.;
          for (unsigned int l=0; l<dim; ++l)
            res += (long double)quadrature_data.derivative(point, n_outer+i)[l] *
                   (long double)quadrature_data.derivative(point, n_outer+j)[l];

          S(i,j) += res * (long double)quadrature.weight(point);
        }

  // Compute the components of T to be the product of gradients of inner and
  // outer shape functions.
  FullMatrix<long double> T(n_inner, n_outer);
  for (unsigned int point=0; point<n_q_points; ++point)
    for (unsigned int i=0; i<n_inner; ++i)
      for (unsigned int k=0; k<n_outer; ++k)
        {
          long double res = 0.;
          for (unsigned int l=0; l<dim; ++l)
            res += (long double)quadrature_data.derivative(point, n_outer+i)[l] *
                   (long double)quadrature_data.derivative(point, k)[l];

          T(i,k) += res *(long double)quadrature.weight(point);
        }

  FullMatrix<long double> S_1(n_inner);
  S_1.invert(S);

  FullMatrix<long double> S_1_T(n_inner, n_outer);

  // S:=S_1*T
  S_1.mmult(S_1_T,T);

  // Resize and initialize the lvs
  lvs.reinit (n_inner, n_outer);
  for (unsigned int i=0; i<n_inner; ++i)
    for (unsigned int k=0; k<n_outer; ++k)
      lvs(i,k) = -S_1_T(i,k);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::apply_laplace_vector(const Table<2,double> &lvs,
                                             std::vector<Point<spacedim> > &a) const
{
  // check whether the data we need is really available. if you fail here and
  // if lvs==laplace_on_quad_vector in the calling function, then we didn't
  // compute the quad laplace vector. this is mentioned in the constructor of
  // this class, although I don't understand the reason for not aborting there
  // any more [WB]
  Assert(lvs.n_rows()!=0, ExcLaplaceVectorNotSet(degree));

  const unsigned int n_inner_apply=lvs.n_rows();
  Assert(n_inner_apply==n_inner || n_inner_apply==(degree-1)*(degree-1),
         ExcInternalError());
  const unsigned int n_outer_apply=lvs.n_cols();
  Assert(a.size()==n_outer_apply,
         ExcDimensionMismatch(a.size(), n_outer_apply));

  // compute each inner point as linear combination of the outer points. the
  // weights are given by the lvs entries, the outer points are the first
  // (existing) elements of a
  for (unsigned int unit_point=0; unit_point<n_inner_apply; ++unit_point)
    {
      Assert(lvs.n_cols()==n_outer_apply, ExcInternalError());
      Point<spacedim> p;
      for (unsigned int k=0; k<n_outer_apply; ++k)
        p+=lvs[unit_point][k]*a[k];

      a.push_back(p);
    }
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim,spacedim>::cell_iterator &cell,
  std::vector<Point<spacedim> > &a) const
{
  // if this is a cell for which we want to compute the full mapping, then get
  // them from the following function
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    compute_support_points_laplace(cell, a);
  else
    // otherwise: use a Q1 mapping for which the mapping shape function
    // support points are simply the vertices of the cell
    {
      a.resize(GeometryInfo<dim>::vertices_per_cell);

      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        a[i] = cell->vertex(i);
    }
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::compute_support_points_laplace(const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                                       std::vector<Point<spacedim> > &a) const
{
  // in any case, we need the vertices first
  a.resize(GeometryInfo<dim>::vertices_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    a[i] = cell->vertex(i);

  if (degree>1)
    switch (dim)
      {
      case 1:
        add_line_support_points(cell, a);
        break;
      case 2:
        // in 2d, add the points on the four bounding lines to the exterior
        // (outer) points
        add_line_support_points (cell, a);
        if (dim != spacedim)
          add_quad_support_points(cell, a);
        else
          apply_laplace_vector (laplace_on_quad_vector,a);
        break;

      case 3:
      {
        // in 3d also add the points located on the boundary faces
        add_line_support_points (cell, a);
        add_quad_support_points (cell, a);
        apply_laplace_vector (laplace_on_hex_vector, a);
        break;
      }
      default:
        Assert(false, ExcNotImplemented());
        break;
      }
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::add_line_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                                                 std::vector<Point<spacedim> > &a) const
{
  // if we only need the midpoint, then ask for it.
  if (degree==2)
    {
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
        {
          const typename Triangulation<dim,spacedim>::line_iterator line =
            (dim == 1  ?
             static_cast<typename Triangulation<dim,spacedim>::line_iterator>(cell) :
             cell->line(line_no));

          const Manifold<dim,spacedim> &manifold =
            ( ( line->manifold_id() == numbers::invalid_manifold_id ) &&
              ( dim < spacedim ) ? cell->get_manifold() :
              line->get_manifold() );
          a.push_back(manifold.get_new_point_on_line(line));
        };
    }
  else
    // otherwise call the more complicated functions and ask for inner points
    // from the boundary description
    {
      std::vector<Point<spacedim> > line_points (degree-1);
      // loop over each of the lines, and if it is at the boundary, then first
      // get the boundary description and second compute the points on it
      for (unsigned int line_no=0; line_no<GeometryInfo<dim>::lines_per_cell; ++line_no)
        {
          const typename Triangulation<dim,spacedim>::line_iterator
          line = (dim == 1
                  ?
                  static_cast<typename Triangulation<dim,spacedim>::line_iterator>(cell)
                  :
                  cell->line(line_no));

          const Manifold<dim,spacedim> &manifold =
            ( ( line->manifold_id() == numbers::invalid_manifold_id ) &&
              ( dim < spacedim ) ? cell->get_manifold() :
              line->get_manifold() );

          get_intermediate_points_on_object (manifold, line, line_points);

          if (dim==3)
            {
              // in 3D, lines might be in wrong orientation. if so, reverse
              // the vector
              if (cell->line_orientation(line_no))
                a.insert (a.end(), line_points.begin(), line_points.end());
              else
                a.insert (a.end(), line_points.rbegin(), line_points.rend());
            }
          else
            // in 2D, lines always have the correct orientation. simply append
            // all points
            a.insert (a.end(), line_points.begin(), line_points.end());

        }
    }
}



template<>
void
MappingQ<3>::
add_quad_support_points(const Triangulation<3>::cell_iterator &cell,
                        std::vector<Point<3> >                &a) const
{
  const unsigned int faces_per_cell    = GeometryInfo<3>::faces_per_cell,
                     vertices_per_face = GeometryInfo<3>::vertices_per_face,
                     lines_per_face    = GeometryInfo<3>::lines_per_face,
                     vertices_per_cell = GeometryInfo<3>::vertices_per_cell;

  static const StraightBoundary<3> straight_boundary;
  // used if face quad at boundary or entirely in the interior of the domain
  std::vector<Point<3> > quad_points ((degree-1)*(degree-1));
  // used if only one line of face quad is at boundary
  std::vector<Point<3> > b(4*degree);

  // Used by the new Manifold interface. This vector collects the
  // vertices used to compute the intermediate points.
  std::vector<Point<3> > vertices(4);

  // loop over all faces and collect points on them
  for (unsigned int face_no=0; face_no<faces_per_cell; ++face_no)
    {
      const Triangulation<3>::face_iterator face = cell->face(face_no);

      // select the correct mappings for the present face
      const bool face_orientation = cell->face_orientation(face_no),
                 face_flip        = cell->face_flip       (face_no),
                 face_rotation    = cell->face_rotation   (face_no);

#ifdef DEBUG
      // some sanity checks up front
      for (unsigned int i=0; i<vertices_per_face; ++i)
        Assert(face->vertex_index(i)==cell->vertex_index(
                 GeometryInfo<3>::face_to_cell_vertices(face_no, i,
                                                        face_orientation,
                                                        face_flip,
                                                        face_rotation)),
               ExcInternalError());

      // indices of the lines that bound a face are given by GeometryInfo<3>::
      // face_to_cell_lines
      for (unsigned int i=0; i<lines_per_face; ++i)
        Assert(face->line(i)==cell->line(GeometryInfo<3>::face_to_cell_lines(
                                           face_no, i, face_orientation, face_flip, face_rotation)),
               ExcInternalError());
#endif

      // if face at boundary, then ask boundary object to return intermediate
      // points on it
      if (face->at_boundary())
        {
          get_intermediate_points_on_object(face->get_manifold(), face, quad_points);

          // in 3D, the orientation, flip and rotation of the face might not
          // match what we expect here, namely the standard orientation. thus
          // reorder points accordingly. since a Mapping uses the same shape
          // function as an FEQ, we can ask a FEQ to do the reordering for us.
          for (unsigned int i=0; i<quad_points.size(); ++i)
            a.push_back(quad_points[feq.adjust_quad_dof_index_for_face_orientation(i,
                                    face_orientation,
                                    face_flip,
                                    face_rotation)]);
        }
      else
        {
          // face is not at boundary, but maybe some of its lines are. count
          // them
          unsigned int lines_at_boundary=0;
          for (unsigned int i=0; i<lines_per_face; ++i)
            if (face->line(i)->at_boundary())
              ++lines_at_boundary;

          Assert(lines_at_boundary<=lines_per_face, ExcInternalError());

          // if at least one of the lines bounding this quad is at the
          // boundary, then collect points separately
          if (lines_at_boundary>0)
            {
              // call of function apply_laplace_vector increases size of b
              // about 1. There resize b for the case the mentioned function
              // was already called.
              b.resize(4*degree);

              // b is of size 4*degree, make sure that this is the right size
              Assert(b.size()==vertices_per_face+lines_per_face*(degree-1),
                     ExcDimensionMismatch(b.size(),
                                          vertices_per_face+lines_per_face*(degree-1)));

              // sort the points into b. We used access from the cell (not
              // from the face) to fill b, so we can assume a standard face
              // orientation. Doing so, the calculated points will be in
              // standard orientation as well.
              for (unsigned int i=0; i<vertices_per_face; ++i)
                b[i]=a[GeometryInfo<3>::face_to_cell_vertices(face_no, i)];

              for (unsigned int i=0; i<lines_per_face; ++i)
                for (unsigned int j=0; j<degree-1; ++j)
                  b[vertices_per_face+i*(degree-1)+j]=
                    a[vertices_per_cell + GeometryInfo<3>::face_to_cell_lines(
                        face_no, i)*(degree-1)+j];

              // Now b includes the support points on the quad and we can
              // apply the laplace vector
              apply_laplace_vector(laplace_on_quad_vector, b);
              Assert(b.size()==4*degree+(degree-1)*(degree-1),
                     ExcDimensionMismatch(b.size(), 4*degree+(degree-1)*(degree-1)));

              for (unsigned int i=0; i<(degree-1)*(degree-1); ++i)
                a.push_back(b[4*degree+i]);
            }
          else
            {
              // face is entirely in the interior. get intermediate
              // points from the relevant manifold object.
              vertices.resize(4);
              for (unsigned int i=0; i<4; ++i)
                vertices[i] = face->vertex(i);
              get_intermediate_points (face->get_manifold(), vertices, quad_points);
              // in 3D, the orientation, flip and rotation of the face might
              // not match what we expect here, namely the standard
              // orientation. thus reorder points accordingly. since a Mapping
              // uses the same shape function as an FEQ, we can ask a FEQ to
              // do the reordering for us.
              for (unsigned int i=0; i<quad_points.size(); ++i)
                a.push_back(quad_points[feq.adjust_quad_dof_index_for_face_orientation(i,
                                        face_orientation,
                                        face_flip,
                                        face_rotation)]);
            }
        }
    }
}



template<>
void
MappingQ<2,3>::
add_quad_support_points(const Triangulation<2,3>::cell_iterator &cell,
                        std::vector<Point<3> >                &a) const
{
  std::vector<Point<3> > quad_points ((degree-1)*(degree-1));
  get_intermediate_points_on_object (cell->get_manifold(), cell, quad_points);
  for (unsigned int i=0; i<quad_points.size(); ++i)
    a.push_back(quad_points[i]);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
add_quad_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &,
                        std::vector<Point<spacedim> > &) const
{
  Assert (dim > 2, ExcImpossibleInDim(dim));
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform (
  const VectorSlice<const std::vector<Tensor<1,dim> > > input,
  VectorSlice<std::vector<Tensor<1,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  // The data object may be just a MappingQ1::InternalData, so we have to test
  // for this first.
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  if (!q1_data->is_mapping_q1_data)
    {
      Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
              ExcInternalError());
      const InternalData &data = static_cast<const InternalData &>(mapping_data);
      // If we only use the Q1-portion, we have to extract that data object
      if (data.use_mapping_q1_on_current_cell)
        q1_data = &data.mapping_q1_data;
    }
  // Now, q1_data should have the right tensors in it and we call the base
  // classes transform function
  MappingQ1<dim,spacedim>::transform(input, output, *q1_data, mapping_type);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::transform (
  const VectorSlice<const std::vector<DerivativeForm<1, dim ,spacedim>  > >  input,
  VectorSlice<std::vector<Tensor<2,spacedim> > > output,
  const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
  const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  // The data object may be just a MappingQ1::InternalData, so we have to test
  // for this first.
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  if (!q1_data->is_mapping_q1_data)
    {
      Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
              ExcInternalError());
      const InternalData &data = static_cast<const InternalData &>(mapping_data);
      // If we only use the Q1-portion, we have to extract that data object
      if (data.use_mapping_q1_on_current_cell)
        q1_data = &data.mapping_q1_data;
    }
  // Now, q1_data should have the right tensors in it and we call the base
  // classes transform function
  MappingQ1<dim,spacedim>::transform(input, output, *q1_data, mapping_type);
}


template<int dim, int spacedim>
void MappingQ<dim,spacedim>::transform
(const VectorSlice<const std::vector<Tensor<2, dim> > >     input,
 VectorSlice<std::vector<Tensor<2,spacedim> > >             output,
 const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
 const MappingType mapping_type) const
{
  AssertDimension (input.size(), output.size());
  // The data object may be just a MappingQ1::InternalData, so we have to test
  // for this first.
  const typename MappingQ1<dim,spacedim>::InternalData *q1_data =
    dynamic_cast<const typename MappingQ1<dim,spacedim>::InternalData *> (&mapping_data);
  Assert(q1_data!=0, ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  if (!q1_data->is_mapping_q1_data)
    {
      Assert (dynamic_cast<const InternalData *>(&mapping_data) != 0,
              ExcInternalError());
      const InternalData &data = static_cast<const InternalData &>(mapping_data);
      // If we only use the Q1-portion, we have to extract that data object
      if (data.use_mapping_q1_on_current_cell)
        q1_data = &data.mapping_q1_data;
    }
  // Now, q1_data should have the right tensors in it and we call the base
  // classes transform function
  MappingQ1<dim,spacedim>::transform(input, output, *q1_data, mapping_type);
}


template<int dim, int spacedim>
Point<spacedim>
MappingQ<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
{
  // Use the get_data function to create an InternalData with data vectors of
  // the right size and transformation shape values already computed at point
  // p.
  const Quadrature<dim> point_quadrature(p);
  std::auto_ptr<InternalData>
  mdata (dynamic_cast<InternalData *> (
           get_data(update_transformation_values, point_quadrature)));

  mdata->use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells ||
                                            cell->has_boundary_lines());

  typename MappingQ1<dim,spacedim>::InternalData
  *p_data = (mdata->use_mapping_q1_on_current_cell ?
             &mdata->mapping_q1_data :
             &*mdata);

  compute_mapping_support_points(cell, p_data->mapping_support_points);
  // If this should be Q1, ignore all other support points.
  if (p_data->shape_values.size()<p_data->mapping_support_points.size())
    p_data->mapping_support_points.resize
    (GeometryInfo<dim>::vertices_per_cell);


  return this->transform_unit_to_real_cell_internal(*p_data);
}



template<int dim, int spacedim>
Point<dim>
MappingQ<dim,spacedim>::
transform_real_to_unit_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<spacedim>                            &p) const
{
  // first a Newton iteration based on a Q1 mapping to get a good starting
  // point, the idea being that this is cheaper than trying to start with the
  // real mapping and likely also more robust.
  //
  // that said, this doesn't always work: there are cases where the point is
  // outside the cell and the inverse mapping doesn't converge. in that case,
  // use the center point of the cell as a starting point if we are to go on
  // using the full mapping, or just propagate up the exception if we had no
  // intention of continuing with the full mapping
  Point<dim> initial_p_unit;
  try
    {
      initial_p_unit
        = MappingQ1<dim,spacedim>::transform_real_to_unit_cell(cell, p);
    }
  catch (const typename Mapping<dim,spacedim>::ExcTransformationFailed &)
    {
      // mirror the conditions of the code below to determine if we need to
      // use an arbitrary starting point or if we just need to rethrow the
      // exception
      if (cell->has_boundary_lines()
          ||
          use_mapping_q_on_all_cells
          ||
          (dim!=spacedim) )
        {
          for (unsigned int d=0; d<dim; ++d)
            initial_p_unit[d] = 0.5;
        }
      else
        throw;
    }

  // then a Newton iteration based on the full MappingQ if we need this. note
  // that for interior cells with dim==spacedim, the mapping used is in fact a
  // Q1 mapping, so there is nothing we need to do unless the iteration above
  // failed
  if (cell->has_boundary_lines()
      ||
      use_mapping_q_on_all_cells
      ||
      (dim!=spacedim) )
    {
      // use the full mapping. in case the function above should have given us
      // something back that lies outside the unit cell (that might happen
      // because we may have given a point 'p' that lies inside the cell with
      // the higher order mapping, but outside the Q1-mapped reference cell),
      // then project it back into the reference cell in hopes that this gives
      // a better starting point to the following iteration
      initial_p_unit = GeometryInfo<dim>::project_to_unit_cell(initial_p_unit);

      const Quadrature<dim> point_quadrature(initial_p_unit);

      UpdateFlags update_flags = update_transformation_values|update_transformation_gradients;
      if (spacedim>dim)
        update_flags |= update_jacobian_grads;
      std::auto_ptr<InternalData>
      mdata (dynamic_cast<InternalData *> (
               get_data(update_flags,point_quadrature)));

      mdata->use_mapping_q1_on_current_cell = false;

      compute_mapping_support_points (cell, mdata->mapping_support_points);

      // If this is a q1 mapping, then only use the support points on the
      // vertices.
      if (mdata->shape_values.size() < mdata->mapping_support_points.size())
        mdata->mapping_support_points.resize(GeometryInfo<dim>::vertices_per_cell);

      return this->transform_real_to_unit_cell_internal(cell, p, initial_p_unit, *mdata);
    }
  else
    return initial_p_unit;
}



template<int dim, int spacedim>
unsigned int
MappingQ<dim,spacedim>::get_degree() const
{
  return degree;
}



template<int dim, int spacedim>
Mapping<dim,spacedim> *
MappingQ<dim,spacedim>::clone () const
{
  return new MappingQ<dim,spacedim>(*this);
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::get_intermediate_points (const Manifold<dim, spacedim> &manifold,
                                                 const std::vector<Point<spacedim> > &surrounding_points,
                                                 std::vector<Point<spacedim> > &points) const
{
  Assert(surrounding_points.size() >= 2, ExcMessage("At least 2 surrounding points are required"));
  const unsigned int n=points.size();
  Assert(n>0, ExcMessage("You can't ask for 0 intermediate points."));
  std::vector<double> w(surrounding_points.size());

  switch (surrounding_points.size())
    {
    case 2:
    {
      // If two points are passed, these are the two vertices, and
      // we can only compute degree-1 intermediate points.
      AssertDimension(n, degree-1);
      for (unsigned int i=0; i<n; ++i)
        {
          const double x = line_support_points.point(i+1)[0];
          w[1] = x;
          w[0] = (1-x);
          Quadrature<spacedim> quadrature(surrounding_points, w);
          points[i] = manifold.get_new_point(quadrature);
        }
      break;
    }

    case 4:
    {
      Assert(spacedim >= 2, ExcImpossibleInDim(spacedim));
      const unsigned m=
        static_cast<unsigned int>(std::sqrt(static_cast<double>(n)));
      // is n a square number
      Assert(m*m==n, ExcInternalError());

      // If four points are passed, these are the two vertices, and
      // we can only compute (degree-1)*(degree-1) intermediate
      // points.
      AssertDimension(m, degree-1);

      for (unsigned int i=0; i<m; ++i)
        {
          const double y=line_support_points.point(1+i)[0];
          for (unsigned int j=0; j<m; ++j)
            {
              const double x=line_support_points.point(1+j)[0];

              w[0] = (1-x)*(1-y);
              w[1] =     x*(1-y);
              w[2] = (1-x)*y    ;
              w[3] =     x*y    ;
              Quadrature<spacedim> quadrature(surrounding_points, w);
              points[i*m+j]=manifold.get_new_point(quadrature);
            }
        }
      break;
    }

    case 8:
      Assert(false, ExcNotImplemented());
      break;
    default:
      Assert(false, ExcInternalError());
      break;
    }
}



// explicit instantiations
#include "mapping_q.inst"


DEAL_II_NAMESPACE_CLOSE
