// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2015 by the deal.II authors
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
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <numeric>

DEAL_II_NAMESPACE_OPEN


template<int dim, int spacedim>
MappingQ<dim,spacedim>::InternalData::InternalData (const unsigned int polynomial_degree)
  :
  MappingQGeneric<dim,spacedim>::InternalData(polynomial_degree),
  use_mapping_q1_on_current_cell(false)
{}



template<int dim, int spacedim>
std::size_t
MappingQ<dim,spacedim>::InternalData::memory_consumption () const
{
  return (MappingQGeneric<dim,spacedim>::InternalData::memory_consumption () +
          MemoryConsumption::memory_consumption (use_mapping_q1_on_current_cell) +
          MemoryConsumption::memory_consumption (mapping_q1_data));
}



template<int dim, int spacedim>
MappingQ<dim,spacedim>::MappingQ (const unsigned int degree,
                                  const bool use_mapping_q_on_all_cells)
  :
  MappingQ1<dim,spacedim>(degree),
  n_inner(Utilities::fixed_power<dim>(degree-1)),
  n_outer((dim==1) ? 2 :
          ((dim==2) ?
           4+4*(degree-1) :
           8+12*(degree-1)+6*(degree-1)*(degree-1))),

  // see whether we want to use *this* mapping objects on *all* cells,
  // or defer to an explicit Q1 mapping on interior cells. if
  // degree==1, then we are already that Q1 mapping, so we don't need
  // it; if dim!=spacedim, there is also no need for anything because
  // we're most likely on a curved manifold
  use_mapping_q_on_all_cells (degree==1
                              ||
                              use_mapping_q_on_all_cells
                              ||
                              (dim != spacedim)),
  feq(degree),
  // create a Q1 mapping for use on interior cells (if necessary)
  // or to create a good initial guess in transform_real_to_unit_cell()
  q1_mapping (new MappingQ1<dim,spacedim>()),
  line_support_points(degree+1)
{
  Assert(n_inner+n_outer==Utilities::fixed_power<dim>(degree+1),
         ExcInternalError());

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
  MappingQ1<dim,spacedim>(mapping.get_degree()),
  n_inner(mapping.n_inner),
  n_outer(mapping.n_outer),
  use_mapping_q_on_all_cells (mapping.use_mapping_q_on_all_cells),
  feq(mapping.get_degree()),
  // clone the Q1 mapping for use on interior cells (if necessary)
  // or to create a good initial guess in transform_real_to_unit_cell()
  q1_mapping (dynamic_cast<MappingQ1<dim,spacedim>*>(mapping.q1_mapping->clone())),
  line_support_points(mapping.line_support_points)
{
  Assert(n_inner+n_outer==Utilities::fixed_power<dim>(this->polynomial_degree+1),
         ExcInternalError());

  // build laplace_on_quad_vector
  if (this->polynomial_degree>1)
    {
      if (dim >= 2)
        set_laplace_on_quad_vector(laplace_on_quad_vector);
      if (dim >= 3)
        set_laplace_on_hex_vector(laplace_on_hex_vector);
    }
}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_data (const UpdateFlags update_flags,
                                  const Quadrature<dim> &quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(quadrature),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_data (update_flags, quadrature));

  task.join ();
  return data;
}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_face_data (const UpdateFlags update_flags,
                                       const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_faces(quadrature));

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize_face,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(q),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_face_data (update_flags, quadrature));

  task.join ();
  return data;
}



template<int dim, int spacedim>
typename MappingQ<dim,spacedim>::InternalData *
MappingQ<dim,spacedim>::get_subface_data (const UpdateFlags update_flags,
                                          const Quadrature<dim-1>& quadrature) const
{
  InternalData *data = new InternalData(this->polynomial_degree);
  const Quadrature<dim> q (QProjector<dim>::project_to_all_subfaces(quadrature));

  // fill the data of both the Q_p and the Q_1 objects in parallel
  Threads::Task<>
  task = Threads::new_task (std_cxx11::function<void()>
                            (std_cxx11::bind(&MappingQGeneric<dim,spacedim>::InternalData::initialize_face,
                                             data,
                                             update_flags,
                                             std_cxx11::cref(q),
                                             quadrature.size())));
  if (!use_mapping_q_on_all_cells)
    data->mapping_q1_data.reset (q1_mapping->get_subface_data (update_flags, quadrature));

  task.join ();
  return data;
}


// Note that the CellSimilarity flag is modifiable, since MappingQ can need to
// recalculate data even when cells are similar.
template<int dim, int spacedim>
CellSimilarity::Similarity
MappingQ<dim,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           cell_similarity,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0, ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is in the interior of the domain
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());


  // call the base class. we need to ensure that the flag indicating whether
  // we can use some similarity has to be modified - for a general MappingQ,
  // the data needs to be recomputed anyway since then the mapping changes the
  // data. this needs to be known also for later operations, so modify the
  // variable here. this also affects the calculation of the next cell -- if
  // we use Q1 data on the next cell, the data will still be invalid.
  const CellSimilarity::Similarity updated_cell_similarity
    = ((data.use_mapping_q1_on_current_cell == false)
       &&
       (this->polynomial_degree > 1)
       ?
       CellSimilarity::invalid_next_cell
       :
       cell_similarity);

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_values (cell,
                                updated_cell_similarity,
                                quadrature,
                                *data.mapping_q1_data,
                                output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_values(cell,
                                                  updated_cell_similarity,
                                                  quadrature,
                                                  data,
                                                  output_data);

  return updated_cell_similarity;
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                     const unsigned int                                         face_no,
                     const Quadrature<dim-1>                                   &quadrature,
                     const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                     internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_face_values (cell,
                                     face_no,
                                     quadrature,
                                     *data.mapping_q1_data,
                                     output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_face_values(cell,
                                                       face_no,
                                                       quadrature,
                                                       data,
                                                       output_data);
}


template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                        const unsigned int                                         face_no,
                        const unsigned int                                         subface_no,
                        const Quadrature<dim-1>                                   &quadrature,
                        const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                        internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // convert data object to internal data for this class. fails with an
  // exception if that is not possible
  Assert (dynamic_cast<const InternalData *> (&internal_data) != 0,
          ExcInternalError());
  const InternalData &data = static_cast<const InternalData &> (internal_data);

  // check whether this cell needs the full mapping or can be treated by a
  // reduced Q1 mapping, e.g. if the cell is entirely in the interior of the
  // domain. note that it is not sufficient to ask whether the present _face_
  // is in the interior, as the mapping on the face depends on the mapping of
  // the cell, which in turn depends on the fact whether _any_ of the faces of
  // this cell is at the boundary, not only the present face
  data.use_mapping_q1_on_current_cell = !(use_mapping_q_on_all_cells
                                          || cell->has_boundary_lines());

  // depending on the results above, decide whether the Q1 mapping or
  // the Qp mapping needs to handle this cell
  if (data.use_mapping_q1_on_current_cell)
    q1_mapping->fill_fe_subface_values (cell,
                                        face_no,
                                        subface_no,
                                        quadrature,
                                        *data.mapping_q1_data,
                                        output_data);
  else
    MappingQGeneric<dim,spacedim>::fill_fe_subface_values(cell,
                                                          face_no,
                                                          subface_no,
                                                          quadrature,
                                                          data,
                                                          output_data);
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
  Assert(this->polynomial_degree>1, ExcInternalError());
  const unsigned int n_inner_2d=(this->polynomial_degree-1)*(this->polynomial_degree-1);
  const unsigned int n_outer_2d=4+4*(this->polynomial_degree-1);

  // first check whether we have precomputed the values for some polynomial
  // degree; the sizes of arrays is n_inner_2d*n_outer_2d
  double const *loqv_ptr=0;
  switch (this->polynomial_degree)
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
          MappingQ<2,2> mapping_2d(this->polynomial_degree);
          loqvs = mapping_2d.laplace_on_quad_vector;
        }
    }

  // the sum of weights of the points at the outer rim should be one. check
  // this
  for (unsigned int unit_point=0; unit_point<loqvs.n_rows(); ++unit_point)
    Assert(std::fabs(std::accumulate(loqvs[unit_point].begin(),
                                     loqvs[unit_point].end(),0.)-1)<1e-13*this->polynomial_degree,
           ExcInternalError());
}



template <>
void
MappingQ<3>::set_laplace_on_hex_vector(Table<2,double> &lohvs) const
{
  Assert(this->polynomial_degree>1, ExcInternalError());

  // first check whether we have precomputed the values for some polynomial
  // degree
  double const *lohv_ptr=0;
  if (this->polynomial_degree==2)
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
                                     lohvs[unit_point].end(),0.) - 1)<1e-13*this->polynomial_degree,
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
  Assert(this->polynomial_degree>1, ExcInternalError());

  // compute the shape gradients at the quadrature points on the unit cell
  const QGauss<dim> quadrature(this->polynomial_degree+1);
  const unsigned int n_q_points=quadrature.size();

  InternalData quadrature_data(this->polynomial_degree);
  quadrature_data.shape_derivatives.resize(quadrature_data.n_shape_functions *
                                           n_q_points);
  quadrature_data.compute_shape_function_values(quadrature.get_points());

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
  Assert(lvs.n_rows()!=0, ExcLaplaceVectorNotSet(this->polynomial_degree));

  const unsigned int n_inner_apply=lvs.n_rows();
  Assert(n_inner_apply==n_inner || n_inner_apply==(this->polynomial_degree-1)*(this->polynomial_degree-1),
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

  if (this->polynomial_degree>1)
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
  if (this->polynomial_degree==2)
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
      std::vector<Point<spacedim> > line_points (this->polynomial_degree-1);
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
              ( dim < spacedim )
              ?
              cell->get_manifold() :
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
  std::vector<Point<3> > quad_points ((this->polynomial_degree-1)*(this->polynomial_degree-1));
  // used if only one line of face quad is at boundary
  std::vector<Point<3> > b(4*this->polynomial_degree);

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
              b.resize(4*this->polynomial_degree);

              // b is of size 4*degree, make sure that this is the right size
              Assert(b.size()==vertices_per_face+lines_per_face*(this->polynomial_degree-1),
                     ExcDimensionMismatch(b.size(),
                                          vertices_per_face+lines_per_face*(this->polynomial_degree-1)));

              // sort the points into b. We used access from the cell (not
              // from the face) to fill b, so we can assume a standard face
              // orientation. Doing so, the calculated points will be in
              // standard orientation as well.
              for (unsigned int i=0; i<vertices_per_face; ++i)
                b[i]=a[GeometryInfo<3>::face_to_cell_vertices(face_no, i)];

              for (unsigned int i=0; i<lines_per_face; ++i)
                for (unsigned int j=0; j<this->polynomial_degree-1; ++j)
                  b[vertices_per_face+i*(this->polynomial_degree-1)+j]=
                    a[vertices_per_cell + GeometryInfo<3>::face_to_cell_lines(
                        face_no, i)*(this->polynomial_degree-1)+j];

              // Now b includes the support points on the quad and we can
              // apply the laplace vector
              apply_laplace_vector(laplace_on_quad_vector, b);
              Assert(b.size()==4*this->polynomial_degree+(this->polynomial_degree-1)*(this->polynomial_degree-1),
                     ExcDimensionMismatch(b.size(), 4*this->polynomial_degree+(this->polynomial_degree-1)*(this->polynomial_degree-1)));

              for (unsigned int i=0; i<(this->polynomial_degree-1)*(this->polynomial_degree-1); ++i)
                a.push_back(b[4*this->polynomial_degree+i]);
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
  std::vector<Point<3> > quad_points ((this->polynomial_degree-1)*(this->polynomial_degree-1));
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
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<1,dim> > >   input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<1,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQGeneric<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  // whether we should in fact work on the Q1 portion of it
  if (const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data))
    if (data->use_mapping_q1_on_current_cell)
      return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);

  // otherwise just stick with what we already had
  MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<DerivativeForm<1, dim ,spacedim>  > >  input,
           const MappingType                                                          mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase                    &mapping_data,
           VectorSlice<std::vector<Tensor<2,spacedim> > >                             output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQGeneric<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  // whether we should in fact work on the Q1 portion of it
  if (const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data))
    if (data->use_mapping_q1_on_current_cell)
      return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);

  // otherwise just stick with what we already had
  MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
void MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<2, dim> > >  input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<2,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQGeneric<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  // whether we should in fact work on the Q1 portion of it
  if (const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data))
    if (data->use_mapping_q1_on_current_cell)
      return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);

  // otherwise just stick with what we already had
  MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}



template<int dim, int spacedim>
void
MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<DerivativeForm<2, dim ,spacedim>  > >  input,
           const MappingType                                                          mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase                    &mapping_data,
           VectorSlice<std::vector<Tensor<3,spacedim> > >                             output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQGeneric<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  // whether we should in fact work on the Q1 portion of it
  if (const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data))
    if (data->use_mapping_q1_on_current_cell)
      return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);

  // otherwise just stick with what we already had
  MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
void MappingQ<dim,spacedim>::
transform (const VectorSlice<const std::vector<Tensor<3, dim> > >  input,
           const MappingType                                       mapping_type,
           const typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
           VectorSlice<std::vector<Tensor<3,spacedim> > >          output) const
{
  AssertDimension (input.size(), output.size());
  Assert ((dynamic_cast<const typename MappingQGeneric<dim,spacedim>::InternalData *> (&mapping_data)
           != 0),
          ExcInternalError());

  // If it is a genuine MappingQ::InternalData, we have to test further
  // whether we should in fact work on the Q1 portion of it
  if (const InternalData *data = dynamic_cast<const InternalData *>(&mapping_data))
    if (data->use_mapping_q1_on_current_cell)
      return q1_mapping->transform (input, mapping_type, *data->mapping_q1_data, output);

  // otherwise just stick with what we already had
  MappingQGeneric<dim,spacedim>::transform(input, mapping_type, mapping_data, output);
}


template<int dim, int spacedim>
Point<spacedim>
MappingQ<dim,spacedim>::
transform_unit_to_real_cell (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                             const Point<dim>                                 &p) const
{
  // first see, whether we want to use a linear or a higher order
  // mapping, then either use our own facilities or that of the Q1
  // mapping we store
  if (use_mapping_q_on_all_cells || cell->has_boundary_lines())
    return this->MappingQGeneric<dim,spacedim>::transform_unit_to_real_cell (cell, p);
  else
    return q1_mapping->transform_unit_to_real_cell (cell, p);
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
      initial_p_unit = q1_mapping->transform_real_to_unit_cell(cell, p);
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

      UpdateFlags update_flags = update_quadrature_points | update_jacobians;
      if (spacedim>dim)
        update_flags |= update_jacobian_grads;
      std_cxx11::unique_ptr<InternalData> mdata (get_data(update_flags,
                                                          point_quadrature));

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
Mapping<dim,spacedim> *
MappingQ<dim,spacedim>::clone () const
{
  return new MappingQ<dim,spacedim>(this->polynomial_degree);
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
      AssertDimension(n, this->polynomial_degree-1);
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
      AssertDimension(m, this->polynomial_degree-1);

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
