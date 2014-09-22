// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_poly_face.templates.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/householder.h>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace
{
  std::vector<Point<1> >
  get_QGaussLobatto_points (const unsigned int degree)
  {
    if (degree > 0)
      return QGaussLobatto<1>(degree+1).get_points();
    else
      return std::vector<Point<1> >(1, Point<1>(0.5));
  }
}

template <int dim, int spacedim>
FE_FaceQ<dim,spacedim>::FE_FaceQ (const unsigned int degree)
  :
  FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim> (
    TensorProductPolynomials<dim-1>(Polynomials::generate_complete_Lagrange_basis(get_QGaussLobatto_points(degree))),
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
    std::vector<bool>(1,true))
{
  // initialize unit face support points
  const unsigned int codim = dim-1;
  this->unit_face_support_points.resize(Utilities::fixed_power<codim>(this->degree+1));

  if (this->degree == 0)
    for (unsigned int d=0; d<codim; ++d)
      this->unit_face_support_points[0][d] = 0.5;
  else
    {
      std::vector<Point<1> > points = get_QGaussLobatto_points(degree);

      unsigned int k=0;
      for (unsigned int iz=0; iz <= ((codim>2) ? this->degree : 0) ; ++iz)
        for (unsigned int iy=0; iy <= ((codim>1) ? this->degree : 0) ; ++iy)
          for (unsigned int ix=0; ix<=this->degree; ++ix)
            {
              Point<codim> p;

              p(0) = points[ix][0];
              if (codim>1)
                p(1) = points[iy][0];
              if (codim>2)
                p(2) = points[iz][0];

              this->unit_face_support_points[k++] = p;
            }
      AssertDimension (k, this->unit_face_support_points.size());
    }

  // initialize unit support points (this makes it possible to assign initial
  // values to FE_FaceQ)
  this->unit_support_points.resize(GeometryInfo<dim>::faces_per_cell*
                                   this->unit_face_support_points.size());
  const unsigned int n_face_dofs = this->unit_face_support_points.size();
  for (unsigned int i=0; i<n_face_dofs; ++i)
    for (unsigned int d=0; d<dim; ++d)
      {
        for (unsigned int e=0, c=0; e<dim; ++e)
          if (d!=e)
            {
              // faces in y-direction are oriented differently
              unsigned int renumber = i;
              if (dim == 3 && d == 1)
                renumber = i/(degree+1)+(degree+1)*(i%(degree+1));
              this->unit_support_points[n_face_dofs*2*d+i][e] =
                this->unit_face_support_points[renumber][c];
              this->unit_support_points[n_face_dofs*(2*d+1)+i][e] =
                this->unit_face_support_points[renumber][c];
              this->unit_support_points[n_face_dofs*(2*d+1)+i][d] = 1;
              ++c;
            }
      }
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_FaceQ<dim,spacedim>::clone() const
{
  return new FE_FaceQ<dim,spacedim>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_FaceQ<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceQ<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_FaceQ<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}



template <int dim, int spacedim>
void
FE_FaceQ<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  // this function is similar to the respective method in FE_Q

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // see if source is a FaceQ element
  if (const FE_FaceQ<dim,spacedim> *source_fe
      = dynamic_cast<const FE_FaceQ<dim,spacedim> *>(&x_source_fe))
    {

      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp procedures,
      // which use this method.
      Assert (this->dofs_per_face <= source_fe->dofs_per_face,
              (typename FiniteElement<dim,spacedim>::
               ExcInterpolationNotImplemented ()));

      // generate a quadrature with the unit face support points.
      const Quadrature<dim-1> face_quadrature (source_fe->get_unit_face_support_points ());

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      const double eps = 2e-13*(this->degree+1)*(dim-1);

      // compute the interpolation matrix by simply taking the value at the
      // support points.
      for (unsigned int i=0; i<source_fe->dofs_per_face; ++i)
        {
          const Point<dim-1> p =
            subface == numbers::invalid_unsigned_int
            ?
            face_quadrature.point(i)
            :
            GeometryInfo<dim-1>::child_to_cell_coordinates (face_quadrature.point(i),
                                                            subface);

          for (unsigned int j=0; j<this->dofs_per_face; ++j)
            {
              double matrix_entry = this->poly_space.compute_value (j, p);

              // Correct the interpolated value. I.e. if it is close to 1 or 0,
              // make it exactly 1 or 0. Unfortunately, this is required to avoid
              // problems with higher order elements.
              if (std::fabs (matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs (matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(i,j) = matrix_entry;
            }
        }

      // make sure that the row sum of each of the matrices is 1 at this
      // point. this must be so since the shape functions sum up to 1
      for (unsigned int j=0; j<source_fe->dofs_per_face; ++j)
        {
          double sum = 0.;

          for (unsigned int i=0; i<this->dofs_per_face; ++i)
            sum += interpolation_matrix(j,i);

          Assert (std::fabs(sum-1) < eps, ExcInternalError());
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != 0)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow (false,(typename FiniteElement<dim,spacedim>::
                        ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
bool
FE_FaceQ<dim,spacedim>::has_support_on_face (
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (face_index == (shape_index/this->dofs_per_face));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_FaceQ<dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim-1] = deg+1;
  for (unsigned int i=1; i<dim-1; ++i)
    dpo[dim-1] *= deg+1;
  return dpo;
}



template <int dim, int spacedim>
bool
FE_FaceQ<dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_FaceQ<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_FaceQ<dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_FaceQ<dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no degrees of freedom and it is typically used in
      // a context where we don't require any continuity along the interface
      return FiniteElementDomination::no_requirements;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_FaceQ<dim,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    constant_modes(0,i) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}



// ----------------------------- FE_FaceQ<1,spacedim> ------------------------

template <int spacedim>
FE_FaceQ<1,spacedim>::FE_FaceQ (const unsigned int degree)
  :
  FiniteElement<1,spacedim> (FiniteElementData<1>(get_dpo_vector(degree), 1, degree, FiniteElementData<1>::L2),
                             std::vector<bool>(1,true),
                             std::vector<ComponentMask> (1, ComponentMask(1,true)))
{
  this->unit_face_support_points.resize(1);

  // initialize unit support points (this makes it possible to assign initial
  // values to FE_FaceQ)
  this->unit_support_points.resize(GeometryInfo<1>::faces_per_cell);
  this->unit_support_points[1] = Point<1>(1.);
}



template <int spacedim>
FiniteElement<1,spacedim> *
FE_FaceQ<1,spacedim>::clone() const
{
  return new FE_FaceQ<1,spacedim>(this->degree);
}



template <int spacedim>
std::string
FE_FaceQ<1,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceQ<"
          << Utilities::dim_string(1,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int spacedim>
void
FE_FaceQ<1,spacedim>::
get_face_interpolation_matrix (const FiniteElement<1,spacedim> &source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}



template <int spacedim>
void
FE_FaceQ<1,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<1,spacedim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));
  interpolation_matrix(0,0) = 1.;
}



template <int spacedim>
bool
FE_FaceQ<1,spacedim>::has_support_on_face (
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  AssertIndexRange(shape_index, 2);
  return (face_index == shape_index);
}



template <int spacedim>
std::vector<unsigned int>
FE_FaceQ<1,spacedim>::get_dpo_vector (const unsigned int)
{
  std::vector<unsigned int> dpo(2, 0U);
  dpo[0] = 1;
  return dpo;
}



template <int spacedim>
bool
FE_FaceQ<1,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int spacedim>
FiniteElementDomination::Domination
FE_FaceQ<1,spacedim>::
compare_for_face_domination (const FiniteElement<1,spacedim> &fe_other) const
{
  return FiniteElementDomination::no_requirements;
}



template <int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_FaceQ<1,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    constant_modes(0,i) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1,0));
}



template <int spacedim>
UpdateFlags
FE_FaceQ<1,spacedim>::update_once (const UpdateFlags) const
{
  return update_default;
}



template <int spacedim>
UpdateFlags
FE_FaceQ<1,spacedim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = flags & update_values;
  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation;
  if (flags & update_cell_normal_vectors)
    out |= update_cell_normal_vectors | update_JxW_values;

  return out;
}



template <int spacedim>
typename Mapping<1,spacedim>::InternalDataBase *
FE_FaceQ<1,spacedim>::get_data (
  const UpdateFlags,
  const Mapping<1,spacedim> &,
  const Quadrature<1> &) const
{
  return new typename Mapping<1,spacedim>::InternalDataBase;
}


template <int spacedim>
typename Mapping<1,spacedim>::InternalDataBase *
FE_FaceQ<1,spacedim>::get_face_data (
  const UpdateFlags update_flags,
  const Mapping<1,spacedim> &,
  const Quadrature<0> &quadrature) const
{
  // generate a new data object and initialize some fields
  typename Mapping<1,spacedim>::InternalDataBase *data =
    new typename Mapping<1,spacedim>::InternalDataBase;

  // check what needs to be initialized only once and what on every
  // cell/face/subface we visit
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.size();
  AssertDimension(n_q_points, 1);

  // No derivatives of this element are implemented.
  if (flags & update_gradients || flags & update_hessians)
    {
      Assert(false, ExcNotImplemented());
    }

  return data;
}



template <int spacedim>
typename Mapping<1,spacedim>::InternalDataBase *
FE_FaceQ<1,spacedim>::get_subface_data (
  const UpdateFlags flags,
  const Mapping<1,spacedim> &mapping,
  const Quadrature<0> &quadrature) const
{
  return get_face_data (flags, mapping, quadrature);
}



template <int spacedim>
void
FE_FaceQ<1,spacedim>::fill_fe_values
(const Mapping<1,spacedim> &,
 const typename Triangulation<1,spacedim>::cell_iterator &,
 const Quadrature<1> &,
 typename Mapping<1,spacedim>::InternalDataBase &,
 typename Mapping<1,spacedim>::InternalDataBase &,
 FEValuesData<1,spacedim> &,
 CellSimilarity::Similarity &) const
{
  // Do nothing, since we do not have values in the interior
}



template <int spacedim>
void
FE_FaceQ<1,spacedim>::fill_fe_face_values (
  const Mapping<1,spacedim> &,
  const typename Triangulation<1,spacedim>::cell_iterator &,
  const unsigned int face,
  const Quadrature<0> &quadrature,
  typename Mapping<1,spacedim>::InternalDataBase &,
  typename Mapping<1,spacedim>::InternalDataBase &fedata,
  FEValuesData<1,spacedim> &data) const
{
  const UpdateFlags flags(fedata.update_once | fedata.update_each);

  const unsigned int foffset = face;
  if (flags & update_values)
    {
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
        data.shape_values(k,0) = 0.;
      data.shape_values(foffset,0) = 1;
    }
}


template <int spacedim>
void
FE_FaceQ<1,spacedim>::fill_fe_subface_values (
  const Mapping<1,spacedim> &,
  const typename Triangulation<1,spacedim>::cell_iterator &,
  const unsigned int ,
  const unsigned int ,
  const Quadrature<0> &,
  typename Mapping<1,spacedim>::InternalDataBase &,
  typename Mapping<1,spacedim>::InternalDataBase &,
  FEValuesData<1,spacedim> &) const
{
  Assert(false, ExcMessage("Should not fille subface values in 1D"));
}



// --------------------------------------- FE_FaceP --------------------------

template <int dim, int spacedim>
FE_FaceP<dim,spacedim>::FE_FaceP (const unsigned int degree)
  :
  FE_PolyFace<PolynomialSpace<dim-1>, dim, spacedim>
  (PolynomialSpace<dim-1>(Polynomials::Legendre::generate_complete_basis(degree)),
   FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
   std::vector<bool>(1,true))
{}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_FaceP<dim,spacedim>::clone() const
{
  return new FE_FaceP<dim,spacedim>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_FaceP<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceP<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
bool
FE_FaceP<dim,spacedim>::has_support_on_face (
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (face_index == (shape_index/this->dofs_per_face));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_FaceP<dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim-1] = deg+1;
  for (unsigned int i=1; i<dim-1; ++i)
    {
      dpo[dim-1] *= deg+1+i;
      dpo[dim-1] /= i+1;
    }
  return dpo;
}




template <int dim, int spacedim>
bool
FE_FaceP<dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_FaceP<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_FaceP<dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_FaceP<dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no degrees of freedom and it is typically used in
      // a context where we don't require any continuity along the interface
      return FiniteElementDomination::no_requirements;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}




template <int dim, int spacedim>
void
FE_FaceP<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}



template <int dim, int spacedim>
void
FE_FaceP<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  // this function is similar to the respective method in FE_Q

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // see if source is a FaceP element
  if (const FE_FaceP<dim,spacedim> *source_fe
      = dynamic_cast<const FE_FaceP<dim,spacedim> *>(&x_source_fe))
    {
      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp procedures,
      // which use this method.
      Assert (this->dofs_per_face <= source_fe->dofs_per_face,
              (typename FiniteElement<dim,spacedim>::
               ExcInterpolationNotImplemented ()));

      // do this as in FETools by solving a least squares problem where we
      // force the source FE polynomial to be equal the given FE on all
      // quadrature points
      const QGauss<dim-1> face_quadrature (source_fe->degree+1);

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      const double eps = 2e-13*(this->degree+1)*(dim-1);

      FullMatrix<double> mass (face_quadrature.size(), source_fe->dofs_per_face);

      for (unsigned int k = 0; k < face_quadrature.size(); ++k)
        {
          const Point<dim-1> p =
            subface == numbers::invalid_unsigned_int ?
            face_quadrature.point(k) :
            GeometryInfo<dim-1>::child_to_cell_coordinates (face_quadrature.point(k),
                                                            subface);

          for (unsigned int j = 0; j < source_fe->dofs_per_face; ++j)
            mass (k , j) = source_fe->poly_space.compute_value(j, p);
        }

      Householder<double> H(mass);
      Vector<double> v_in(face_quadrature.size());
      Vector<double> v_out(source_fe->dofs_per_face);


      // compute the interpolation matrix by evaluating on the fine side and
      // then solving the least squares problem
      for (unsigned int i=0; i<this->dofs_per_face; ++i)
        {
          for (unsigned int k = 0; k < face_quadrature.size(); ++k)
            {
              const Point<dim-1> p = numbers::invalid_unsigned_int ?
                                     face_quadrature.point(k) :
                                     GeometryInfo<dim-1>::child_to_cell_coordinates (face_quadrature.point(k),
                                         subface);
              v_in(k) = this->poly_space.compute_value(i, p);
            }
          const double result = H.least_squares(v_out, v_in);
          Assert(result < 1e-12, FETools::ExcLeastSquaresError (result));

          for (unsigned int j = 0; j < source_fe->dofs_per_face; ++j)
            {
              double matrix_entry = v_out(j);

              // Correct the interpolated value. I.e. if it is close to 1 or 0,
              // make it exactly 1 or 0. Unfortunately, this is required to avoid
              // problems with higher order elements.
              if (std::fabs (matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs (matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(j,i) = matrix_entry;
            }
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != 0)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow (false,(typename FiniteElement<dim,spacedim>::
                        ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_FaceP<dim,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    constant_modes(0, face*this->dofs_per_face) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}



template <int spacedim>
FE_FaceP<1,spacedim>::FE_FaceP (const unsigned int degree)
  :
  FE_FaceQ<1,spacedim> (degree)
{}



template <int spacedim>
std::string
FE_FaceP<1,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceP<"
          << Utilities::dim_string(1,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



// explicit instantiations
#include "fe_face.inst"


DEAL_II_NAMESPACE_CLOSE
