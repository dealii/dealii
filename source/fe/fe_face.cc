// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_poly_face.templates.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/householder.h>

#include <memory>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace FE_FaceQImplementation
  {
    namespace
    {
      std::vector<Point<1>>
      get_QGaussLobatto_points(const unsigned int degree)
      {
        if (degree > 0)
          return QGaussLobatto<1>(degree + 1).get_points();
        else
          return std::vector<Point<1>>(1, Point<1>(0.5));
      }
    } // namespace
  }   // namespace FE_FaceQImplementation
} // namespace internal

template <int dim, int spacedim>
FE_FaceQ<dim, spacedim>::FE_FaceQ(const unsigned int degree)
  : FE_PolyFace<TensorProductPolynomials<dim - 1>, dim, spacedim>(
      TensorProductPolynomials<dim - 1>(
        Polynomials::generate_complete_Lagrange_basis(
          internal::FE_FaceQImplementation::get_QGaussLobatto_points(degree))),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(1, true))
{
  // initialize unit face support points
  const unsigned int codim = dim - 1;
  this->unit_face_support_points[0].resize(
    Utilities::fixed_power<codim>(this->degree + 1));

  if (this->degree == 0)
    for (unsigned int d = 0; d < codim; ++d)
      this->unit_face_support_points[0][0][d] = 0.5;
  else
    {
      std::vector<Point<1>> points =
        internal::FE_FaceQImplementation::get_QGaussLobatto_points(degree);

      unsigned int k = 0;
      for (unsigned int iz = 0; iz <= ((codim > 2) ? this->degree : 0); ++iz)
        for (unsigned int iy = 0; iy <= ((codim > 1) ? this->degree : 0); ++iy)
          for (unsigned int ix = 0; ix <= this->degree; ++ix)
            {
              Point<codim> p;

              p[0] = points[ix][0];
              if (codim > 1)
                p[1] = points[iy][0];
              if (codim > 2)
                p[2] = points[iz][0];

              this->unit_face_support_points[0][k++] = p;
            }
      AssertDimension(k, this->unit_face_support_points[0].size());
    }

  // initialize unit support points (this makes it possible to assign initial
  // values to FE_FaceQ)
  this->unit_support_points.resize(GeometryInfo<dim>::faces_per_cell *
                                   this->unit_face_support_points[0].size());
  const unsigned int n_face_dofs = this->unit_face_support_points[0].size();
  for (unsigned int i = 0; i < n_face_dofs; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      {
        for (unsigned int e = 0, c = 0; e < dim; ++e)
          if (d != e)
            {
              // faces in y-direction are oriented differently
              unsigned int renumber = i;
              if (dim == 3 && d == 1)
                renumber = i / (degree + 1) + (degree + 1) * (i % (degree + 1));
              this->unit_support_points[n_face_dofs * 2 * d + i][e] =
                this->unit_face_support_points[0][renumber][c];
              this->unit_support_points[n_face_dofs * (2 * d + 1) + i][e] =
                this->unit_face_support_points[0][renumber][c];
              this->unit_support_points[n_face_dofs * (2 * d + 1) + i][d] = 1;
              ++c;
            }
      }
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_FaceQ<dim, spacedim>::clone() const
{
  return std::make_unique<FE_FaceQ<dim, spacedim>>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_FaceQ<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceQ<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
void
FE_FaceQ<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  get_subface_interpolation_matrix(source_fe,
                                   numbers::invalid_unsigned_int,
                                   interpolation_matrix,
                                   face_no);
}



template <int dim, int spacedim>
void
FE_FaceQ<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int                  subface,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  // this function is similar to the respective method in FE_Q

  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // see if source is a FaceQ element
  if (const FE_FaceQ<dim, spacedim> *source_fe =
        dynamic_cast<const FE_FaceQ<dim, spacedim> *>(&x_source_fe))
    {
      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp-procedures,
      // which use this method.
      Assert(
        this->n_dofs_per_face(face_no) <= source_fe->n_dofs_per_face(face_no),
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));

      // generate a quadrature with the unit face support points.
      const Quadrature<dim - 1> face_quadrature(
        source_fe->get_unit_face_support_points(face_no));

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      const double eps = 2e-13 * (this->degree + 1) * (dim - 1);

      // compute the interpolation matrix by simply taking the value at the
      // support points.
      for (unsigned int i = 0; i < source_fe->n_dofs_per_face(face_no); ++i)
        {
          const Point<dim - 1> p =
            subface == numbers::invalid_unsigned_int ?
              face_quadrature.point(i) :
              GeometryInfo<dim - 1>::child_to_cell_coordinates(
                face_quadrature.point(i), subface);

          for (unsigned int j = 0; j < this->n_dofs_per_face(face_no); ++j)
            {
              double matrix_entry = this->poly_space.compute_value(j, p);

              // Correct the interpolated value. I.e. if it is close to 1 or 0,
              // make it exactly 1 or 0. Unfortunately, this is required to
              // avoid problems with higher order elements.
              if (std::fabs(matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs(matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(i, j) = matrix_entry;
            }
        }

      if constexpr (running_in_debug_mode())
        {
          // make sure that the row sum of each of the matrices is 1 at this
          // point. this must be so since the shape functions sum up to 1
          for (unsigned int j = 0; j < source_fe->n_dofs_per_face(face_no); ++j)
            {
              double sum = 0.;

              for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
                sum += interpolation_matrix(j, i);

              Assert(std::fabs(sum - 1) < eps, ExcInternalError());
            }
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != nullptr)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow(
      false,
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
bool
FE_FaceQ<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (face_index == (shape_index / this->n_dofs_per_face(face_index)));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_FaceQ<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 0U);
  dpo[dim - 1] = deg + 1;
  for (unsigned int i = 1; i < dim - 1; ++i)
    dpo[dim - 1] *= deg + 1;
  return dpo;
}



template <int dim, int spacedim>
bool
FE_FaceQ<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> & /*fe_other*/) const
{
  // this element is always discontinuous at vertices
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  Assert(dim >= 2, ExcInternalError());

  // this element is continuous only for the highest dimensional bounding object
  if (dim > 2)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      // this is similar to the FE_Q_Base class
      if (const FE_FaceQ<dim, spacedim> *fe_q_other =
            dynamic_cast<const FE_FaceQ<dim, spacedim> *>(&fe_other))
        {
          // dofs are located along lines, so two dofs are identical if they are
          // located at identical positions.
          // Therefore, read the points in unit_support_points for the
          // first coordinate direction. We take the lexicographic ordering of
          // the points in the second direction (i.e., y-direction) since we
          // know that the first p+1 dofs are located at the left (x=0) face.
          const unsigned int p = this->degree;
          const unsigned int q = fe_q_other->degree;

          std::vector<std::pair<unsigned int, unsigned int>> identities;

          const std::vector<unsigned int> &index_map_inverse =
            this->poly_space.get_numbering_inverse();
          const std::vector<unsigned int> &index_map_inverse_other =
            fe_q_other->poly_space.get_numbering_inverse();

          for (unsigned int i = 0; i < p + 1; ++i)
            for (unsigned int j = 0; j < q + 1; ++j)
              if (std::fabs(
                    this->unit_support_points[index_map_inverse[i]][dim - 1] -
                    fe_q_other->unit_support_points[index_map_inverse_other[j]]
                                                   [dim - 1]) < 1e-14)
                identities.emplace_back(i, j);

          return identities;
        }
      else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
        {
          // the FE_Nothing has no degrees of freedom, so there are no
          // equivalencies to be recorded
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
      else if (fe_other.n_unique_faces() == 1 &&
               fe_other.n_dofs_per_face(0) == 0)
        {
          // if the other element has no elements on faces at all,
          // then it would be impossible to enforce any kind of
          // continuity even if we knew exactly what kind of element
          // we have -- simply because the other element declares
          // that it is discontinuous because it has no DoFs on
          // its faces. in that case, just state that we have no
          // constraints to declare
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
      else
        {
          DEAL_II_NOT_IMPLEMENTED();
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int) const
{
  Assert(dim >= 3, ExcInternalError());

  // this element is continuous only for the highest dimensional bounding object
  if (dim > 3)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      // this is similar to the FE_Q_Base class
      if (const FE_FaceQ<dim, spacedim> *fe_q_other =
            dynamic_cast<const FE_FaceQ<dim, spacedim> *>(&fe_other))
        {
          // this works exactly like the line case above, except that now we
          // have to have two indices i1, i2 and j1, j2 to characterize the dofs
          // on the face of each of the finite elements. since they are ordered
          // lexicographically along the first line and we have a tensor
          // product, the rest is rather straightforward
          const unsigned int p = this->degree;
          const unsigned int q = fe_q_other->degree;

          std::vector<std::pair<unsigned int, unsigned int>> identities;

          const std::vector<unsigned int> &index_map_inverse =
            this->poly_space.get_numbering_inverse();
          const std::vector<unsigned int> &index_map_inverse_other =
            fe_q_other->poly_space.get_numbering_inverse();

          std::vector<std::pair<unsigned int, unsigned int>> identities_1d;

          for (unsigned int i = 0; i < p + 1; ++i)
            for (unsigned int j = 0; j < q + 1; ++j)
              if (std::fabs(
                    this->unit_support_points[index_map_inverse[i]][dim - 2] -
                    fe_q_other->unit_support_points[index_map_inverse_other[j]]
                                                   [dim - 2]) < 1e-14)
                identities_1d.emplace_back(i, j);

          for (unsigned int n1 = 0; n1 < identities_1d.size(); ++n1)
            for (unsigned int n2 = 0; n2 < identities_1d.size(); ++n2)
              identities.emplace_back(identities_1d[n1].first * (p + 1) +
                                        identities_1d[n2].first,
                                      identities_1d[n1].second * (q + 1) +
                                        identities_1d[n2].second);

          return identities;
        }
      else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
        {
          // the FE_Nothing has no degrees of freedom, so there are no
          // equivalencies to be recorded
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
      else if (fe_other.n_unique_faces() == 1 &&
               fe_other.n_dofs_per_face(0) == 0)
        {
          // if the other element has no elements on faces at all,
          // then it would be impossible to enforce any kind of
          // continuity even if we knew exactly what kind of element
          // we have -- simply because the other element declares
          // that it is discontinuous because it has no DoFs on
          // its faces. in that case, just state that we have no
          // constraints to declare
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
      else
        {
          DEAL_II_NOT_IMPLEMENTED();
          return std::vector<std::pair<unsigned int, unsigned int>>();
        }
    }
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_FaceQ<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face/cell domination
  // --------------------------------
  if (const FE_FaceQ<dim, spacedim> *fe_faceq_other =
        dynamic_cast<const FE_FaceQ<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_faceq_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_faceq_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}

template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_FaceQ<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    constant_modes(0, i) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}

template <int dim, int spacedim>
void
FE_FaceQ<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_values) const
{
  AssertDimension(support_point_values.size(),
                  this->get_unit_support_points().size());
  AssertDimension(support_point_values.size(), nodal_values.size());
  AssertDimension(this->n_dofs_per_cell(), nodal_values.size());

  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    {
      AssertDimension(support_point_values[i].size(), 1);

      nodal_values[i] = support_point_values[i](0);
    }
}

// ----------------------------- FE_FaceQ<1,spacedim> ------------------------

template <int spacedim>
FE_FaceQ<1, spacedim>::FE_FaceQ(const unsigned int degree)
  : FiniteElement<1, spacedim>(
      FiniteElementData<1>(get_dpo_vector(degree),
                           1,
                           degree,
                           FiniteElementData<1>::L2),
      std::vector<bool>(1, true),
      std::vector<ComponentMask>(1, ComponentMask(1, true)))
{
  this->unit_face_support_points[0].resize(1);

  // initialize unit support points (this makes it possible to assign initial
  // values to FE_FaceQ)
  this->unit_support_points.resize(GeometryInfo<1>::faces_per_cell);
  this->unit_support_points[1] = Point<1>(1.);
}



template <int spacedim>
std::unique_ptr<FiniteElement<1, spacedim>>
FE_FaceQ<1, spacedim>::clone() const
{
  return std::make_unique<FE_FaceQ<1, spacedim>>(this->degree);
}



template <int spacedim>
std::string
FE_FaceQ<1, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceQ<" << Utilities::dim_string(1, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int spacedim>
void
FE_FaceQ<1, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<1, spacedim> &source_fe,
  FullMatrix<double>               &interpolation_matrix,
  const unsigned int                face_no) const
{
  get_subface_interpolation_matrix(source_fe,
                                   numbers::invalid_unsigned_int,
                                   interpolation_matrix,
                                   face_no);
}



template <int spacedim>
void
FE_FaceQ<1, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<1, spacedim> &x_source_fe,
  const unsigned int /*subface*/,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int  face_no) const
{
  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));
  interpolation_matrix(0, 0) = 1.;
}



template <int spacedim>
bool
FE_FaceQ<1, spacedim>::has_support_on_face(const unsigned int shape_index,
                                           const unsigned int face_index) const
{
  AssertIndexRange(shape_index, 2);
  return (face_index == shape_index);
}



template <int spacedim>
std::vector<unsigned int>
FE_FaceQ<1, spacedim>::get_dpo_vector(const unsigned int)
{
  std::vector<unsigned int> dpo(2, 0U);
  dpo[0] = 1;
  return dpo;
}



template <int spacedim>
bool
FE_FaceQ<1, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}

template <int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<1, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<1, spacedim> & /*fe_other*/) const
{
  // this element is always discontinuous at vertices
  return std::vector<std::pair<unsigned int, unsigned int>>(1,
                                                            std::make_pair(0U,
                                                                           0U));
}



template <int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<1, spacedim>::hp_line_dof_identities(
  const FiniteElement<1, spacedim> &) const
{
  // this element is continuous only for the highest dimensional bounding object
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_FaceQ<1, spacedim>::hp_quad_dof_identities(
  const FiniteElement<1, spacedim> &,
  const unsigned int) const
{
  // this element is continuous only for the highest dimensional bounding object
  return std::vector<std::pair<unsigned int, unsigned int>>();
}



template <int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_FaceQ<1, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    constant_modes(0, i) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int spacedim>
UpdateFlags
FE_FaceQ<1, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = flags & update_values;
  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation;
  if (flags & update_normal_vectors)
    out |= update_normal_vectors | update_JxW_values;

  return out;
}


template <int spacedim>
void
FE_FaceQ<1, spacedim>::fill_fe_values(
  const typename Triangulation<1, spacedim>::cell_iterator &,
  const CellSimilarity::Similarity,
  const Quadrature<1> &,
  const Mapping<1, spacedim> &,
  const typename Mapping<1, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &,
  const typename FiniteElement<1, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                     spacedim>
    &) const
{
  // Do nothing, since we do not have values in the interior
}



template <int spacedim>
void
FE_FaceQ<1, spacedim>::fill_fe_face_values(
  const typename Triangulation<1, spacedim>::cell_iterator &,
  const unsigned int face,
  const hp::QCollection<0> &,
  const Mapping<1, spacedim> &,
  const typename Mapping<1, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &,
  const typename FiniteElement<1, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                     spacedim>
    &output_data) const
{
  const unsigned int foffset = face;
  if (fe_internal.update_each & update_values)
    {
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        output_data.shape_values(k, 0) = 0.;
      output_data.shape_values(foffset, 0) = 1;
    }
}


template <int spacedim>
void
FE_FaceQ<1, spacedim>::fill_fe_subface_values(
  const typename Triangulation<1, spacedim>::cell_iterator &,
  const unsigned int,
  const unsigned int,
  const Quadrature<0> &,
  const Mapping<1, spacedim> &,
  const typename Mapping<1, spacedim>::InternalDataBase &,
  const internal::FEValuesImplementation::MappingRelatedData<1, spacedim> &,
  const typename FiniteElement<1, spacedim>::InternalDataBase &,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<1,
                                                                     spacedim>
    &) const
{
  Assert(false, ExcMessage("There are no sub-face values to fill in 1d!"));
}



// --------------------------------------- FE_FaceP --------------------------

template <int dim, int spacedim>
FE_FaceP<dim, spacedim>::FE_FaceP(const unsigned int degree)
  : FE_PolyFace<PolynomialSpace<dim - 1>, dim, spacedim>(
      PolynomialSpace<dim - 1>(
        Polynomials::Legendre::generate_complete_basis(degree)),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(1, true))
{}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_FaceP<dim, spacedim>::clone() const
{
  return std::make_unique<FE_FaceP<dim, spacedim>>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_FaceP<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceP<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
bool
FE_FaceP<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (face_index == (shape_index / this->n_dofs_per_face(face_index)));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_FaceP<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 0U);
  dpo[dim - 1] = deg + 1;
  for (unsigned int i = 1; i < dim - 1; ++i)
    {
      dpo[dim - 1] *= deg + 1 + i;
      dpo[dim - 1] /= i + 1;
    }
  return dpo;
}



template <int dim, int spacedim>
bool
FE_FaceP<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_FaceP<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face/cell domination
  // --------------------------------
  if (const FE_FaceP<dim, spacedim> *fe_facep_other =
        dynamic_cast<const FE_FaceP<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_facep_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_facep_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
void
FE_FaceP<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source_fe,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  get_subface_interpolation_matrix(source_fe,
                                   numbers::invalid_unsigned_int,
                                   interpolation_matrix,
                                   face_no);
}



template <int dim, int spacedim>
void
FE_FaceP<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int                  subface,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  // this function is similar to the respective method in FE_Q

  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // see if source is a FaceP element
  if (const FE_FaceP<dim, spacedim> *source_fe =
        dynamic_cast<const FE_FaceP<dim, spacedim> *>(&x_source_fe))
    {
      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp-procedures,
      // which use this method.
      Assert(
        this->n_dofs_per_face(face_no) <= source_fe->n_dofs_per_face(face_no),
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));

      // do this as in FETools by solving a least squares problem where we
      // force the source FE polynomial to be equal the given FE on all
      // quadrature points
      const QGauss<dim - 1> face_quadrature(source_fe->degree + 1);

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      const double eps = 2e-13 * (this->degree + 1) * (dim - 1);

      FullMatrix<double> mass(face_quadrature.size(),
                              source_fe->n_dofs_per_face(face_no));

      for (unsigned int k = 0; k < face_quadrature.size(); ++k)
        {
          const Point<dim - 1> p =
            subface == numbers::invalid_unsigned_int ?
              face_quadrature.point(k) :
              GeometryInfo<dim - 1>::child_to_cell_coordinates(
                face_quadrature.point(k), subface);

          for (unsigned int j = 0; j < source_fe->n_dofs_per_face(face_no); ++j)
            mass(k, j) = source_fe->poly_space.compute_value(j, p);
        }

      Householder<double> H(mass);
      Vector<double>      v_in(face_quadrature.size());
      Vector<double>      v_out(source_fe->n_dofs_per_face(face_no));


      // compute the interpolation matrix by evaluating on the fine side and
      // then solving the least squares problem
      for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
        {
          for (unsigned int k = 0; k < face_quadrature.size(); ++k)
            {
              const Point<dim - 1> p =
                subface == numbers::invalid_unsigned_int ?
                  face_quadrature.point(k) :
                  GeometryInfo<dim - 1>::child_to_cell_coordinates(
                    face_quadrature.point(k), subface);
              v_in(k) = this->poly_space.compute_value(i, p);
            }
          [[maybe_unused]] const double result = H.least_squares(v_out, v_in);
          Assert(result < 1e-12, FETools::ExcLeastSquaresError(result));

          for (unsigned int j = 0; j < source_fe->n_dofs_per_face(face_no); ++j)
            {
              double matrix_entry = v_out(j);

              // Correct the interpolated value. I.e. if it is close to 1 or 0,
              // make it exactly 1 or 0. Unfortunately, this is required to
              // avoid problems with higher order elements.
              if (std::fabs(matrix_entry - 1.0) < eps)
                matrix_entry = 1.0;
              if (std::fabs(matrix_entry) < eps)
                matrix_entry = 0.0;

              interpolation_matrix(j, i) = matrix_entry;
            }
        }
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != nullptr)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow(
      false,
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_FaceP<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  for (const unsigned int face : GeometryInfo<dim>::face_indices())
    constant_modes(0, face * this->n_dofs_per_face(face)) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}



template <int spacedim>
FE_FaceP<1, spacedim>::FE_FaceP(const unsigned int degree)
  : FE_FaceQ<1, spacedim>(degree)
{}



template <int spacedim>
std::string
FE_FaceP<1, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_FaceP<" << Utilities::dim_string(1, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



// explicit instantiations
#include "fe/fe_face.inst"


DEAL_II_NAMESPACE_CLOSE
