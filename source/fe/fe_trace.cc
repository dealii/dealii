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

#include <deal.II/base/config.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_trace.h>

#include <sstream>
#include <fstream>
#include <iostream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_TraceQ<dim,spacedim>::FE_TraceQ (const unsigned int degree)
  :
  FE_PolyFace<TensorProductPolynomials<dim-1>, dim, spacedim> (
    TensorProductPolynomials<dim-1>(Polynomials::generate_complete_Lagrange_basis(QGaussLobatto<1>(degree+1).get_points())),
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
    std::vector<bool>(1,true)),
  fe_q (degree)
{
  Assert (degree > 0,
          ExcMessage ("FE_Trace can only be used for polynomial degrees "
                      "greater than zero"));
  std::vector<unsigned int> renumber (this->dofs_per_face);
  FETools::hierarchic_to_lexicographic_numbering<dim-1> (degree, renumber);
  this->poly_space.set_numbering(renumber);

  // Initialize face support points
  this->unit_face_support_points = fe_q.get_unit_face_support_points();
  // Initialize constraint matrices
  this->interface_constraints = fe_q.constraints();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_TraceQ<dim,spacedim>::clone() const
{
  return new FE_TraceQ<dim,spacedim>(this->degree);
}



template <int dim, int spacedim>
std::string
FE_TraceQ<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_TraceQ<"
          << Utilities::dim_string(dim,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
bool
FE_TraceQ<dim,spacedim>::has_support_on_face (const unsigned int shape_index,
                                              const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  // FE_TraceQ shares the numbering of elemental degrees of freedom with FE_Q
  // except for the missing interior ones (quad dofs in 2D and hex dofs in
  // 3D). Therefore, it is safe to ask fe_q for the corresponding
  // information. The assertion 'shape_index < this->dofs_per_cell' will make
  // sure that we only access the trace dofs.
  return fe_q.has_support_on_face (shape_index, face_index);
}



template <int dim, int spacedim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_TraceQ<dim,spacedim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    constant_modes(0,i) = true;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_TraceQ<dim,spacedim>::get_dpo_vector (const unsigned int deg)
{
  // This constructs FE_TraceQ in exactly the same way as FE_Q except for the
  // interior degrees of freedom that are not present here (line in 1D, quad
  // in 2D, hex in 3D).
  AssertThrow(deg>0,ExcMessage("FE_TraceQ needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim+1, 1U);
  dpo[dim]=0;
  dpo[0]=1;
  for (unsigned int i=1; i<dim; ++i)
    dpo[i] = dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim, int spacedim>
bool
FE_TraceQ<dim,spacedim>::hp_constraints_are_implemented () const
{
  return fe_q.hp_constraints_are_implemented ();
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_TraceQ<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_TraceQ<dim,spacedim> *fe_q_other
      = dynamic_cast<const FE_TraceQ<dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_q_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_q_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing = dynamic_cast<const FE_Nothing<dim>*>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        {
          return FiniteElementDomination::other_element_dominates;
        }
      else
        {
          // the FE_Nothing has no degrees of freedom and it is typically used in
          // a context where we don't require any continuity along the interface
          return FiniteElementDomination::no_requirements;
        }
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
void
FE_TraceQ<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}



template <int dim, int spacedim>
void
FE_TraceQ<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  // this is the code from FE_FaceQ
  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // see if source is a FaceQ element
  if (const FE_TraceQ<dim,spacedim> *source_fe
      = dynamic_cast<const FE_TraceQ<dim,spacedim> *>(&x_source_fe))
    {
      fe_q.get_subface_interpolation_matrix (source_fe->fe_q, subface, interpolation_matrix);
    }
  else if (dynamic_cast<const FE_Nothing<dim> *>(&x_source_fe) != 0)
    {
      // nothing to do here, the FE_Nothing has no degrees of freedom anyway
    }
  else
    AssertThrow (false,(typename FiniteElement<dim,spacedim>::
                        ExcInterpolationNotImplemented()));
}



template <int spacedim>
FE_TraceQ<1,spacedim>::FE_TraceQ (const unsigned int degree)
  :
  FE_FaceQ<1,spacedim> (degree)
{}



template <int spacedim>
std::string
FE_TraceQ<1,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch
  std::ostringstream namebuf;
  namebuf << "FE_TraceQ<"
          << Utilities::dim_string(1,spacedim)
          << ">(" << this->degree << ")";

  return namebuf.str();
}



// explicit instantiations
#include "fe_trace.inst"


DEAL_II_NAMESPACE_CLOSE
