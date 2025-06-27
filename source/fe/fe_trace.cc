// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_trace.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_TraceQ<dim, spacedim>::FE_TraceQ(const unsigned int degree)
  : FE_PolyFace<TensorProductPolynomials<dim - 1>, dim, spacedim>(
      TensorProductPolynomials<dim - 1>(
        Polynomials::generate_complete_Lagrange_basis(
          QGaussLobatto<1>(degree + 1).get_points())),
      FiniteElementData<dim>(get_dpo_vector(degree),
                             1,
                             degree,
                             FiniteElementData<dim>::L2),
      std::vector<bool>(1, true))
  , fe_q(degree)
{
  Assert(degree > 0,
         ExcMessage("FE_Trace can only be used for polynomial degrees "
                    "greater than zero"));
  this->poly_space.set_numbering(
    FETools::hierarchic_to_lexicographic_numbering<dim - 1>(degree));

  // Initialize face support points
  AssertDimension(this->n_unique_faces(), fe_q.n_unique_faces());
  for (unsigned int face_no = 0; face_no < this->n_unique_faces(); ++face_no)
    this->unit_face_support_points[face_no] =
      fe_q.get_unit_face_support_points(face_no);

  // initialize unit support points (this makes it possible to assign initial
  // values to FE_TraceQ). Note that we simply take the points of fe_q but
  // skip the last ones which are associated with the interior of FE_Q.
  this->unit_support_points.resize(this->n_dofs_per_cell());
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    this->unit_support_points[i] = fe_q.get_unit_support_points()[i];

  // Initialize constraint matrices
  this->interface_constraints = fe_q.constraints();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_TraceQ<dim, spacedim>::clone() const
{
  return std::make_unique<FE_TraceQ<dim, spacedim>>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_TraceQ<dim, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_TraceQ<" << Utilities::dim_string(dim, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
bool
FE_TraceQ<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  AssertIndexRange(shape_index, this->n_dofs_per_cell());
  AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);

  // FE_TraceQ shares the numbering of elemental degrees of freedom with FE_Q
  // except for the missing interior ones (quad dofs in 2d and hex dofs in
  // 3d). Therefore, it is safe to ask fe_q for the corresponding
  // information. The assertion 'shape_index < this->n_dofs_per_cell()' will
  // make sure that we only access the trace dofs.
  return fe_q.has_support_on_face(shape_index, face_index);
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FE_TraceQ<dim, spacedim>::get_constant_modes() const
{
  Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    constant_modes(0, i) = true;
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(
    constant_modes, std::vector<unsigned int>(1, 0));
}

template <int dim, int spacedim>
void
FE_TraceQ<dim, spacedim>::
  convert_generalized_support_point_values_to_dof_values(
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


template <int dim, int spacedim>
std::vector<unsigned int>
FE_TraceQ<dim, spacedim>::get_dpo_vector(const unsigned int deg)
{
  // This constructs FE_TraceQ in exactly the same way as FE_Q except for the
  // interior degrees of freedom that are not present here (line in 1d, quad
  // in 2d, hex in 3d).
  AssertThrow(deg > 0, ExcMessage("FE_TraceQ needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim + 1, 1U);
  dpo[dim] = 0;
  dpo[0]   = 1;
  for (unsigned int i = 1; i < dim; ++i)
    dpo[i] = dpo[i - 1] * (deg - 1);
  return dpo;
}



template <int dim, int spacedim>
bool
FE_TraceQ<dim, spacedim>::hp_constraints_are_implemented() const
{
  return fe_q.hp_constraints_are_implemented();
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_TraceQ<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));
  (void)codim;

  // vertex/line/face/cell domination
  // --------------------------------
  if (const FE_TraceQ<dim, spacedim> *fe_traceq_other =
        dynamic_cast<const FE_TraceQ<dim, spacedim> *>(&fe_other))
    {
      if (this->degree < fe_traceq_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_traceq_other->degree)
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
FE_TraceQ<dim, spacedim>::get_face_interpolation_matrix(
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
FE_TraceQ<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int                  subface,
  FullMatrix<double>                 &interpolation_matrix,
  const unsigned int                  face_no) const
{
  // this is the code from FE_FaceQ
  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // see if source is a FaceQ element
  if (const FE_TraceQ<dim, spacedim> *source_fe =
        dynamic_cast<const FE_TraceQ<dim, spacedim> *>(&x_source_fe))
    {
      fe_q.get_subface_interpolation_matrix(source_fe->fe_q,
                                            subface,
                                            interpolation_matrix,
                                            face_no);
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



template <int spacedim>
FE_TraceQ<1, spacedim>::FE_TraceQ(const unsigned int degree)
  : FE_FaceQ<1, spacedim>(degree)
{}



template <int spacedim>
std::string
FE_TraceQ<1, spacedim>::get_name() const
{
  // note that the FETools::get_fe_by_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_TraceQ<" << Utilities::dim_string(1, spacedim) << ">("
          << this->degree << ")";

  return namebuf.str();
}



// explicit instantiations
#include "fe/fe_trace.inst"


DEAL_II_NAMESPACE_CLOSE
