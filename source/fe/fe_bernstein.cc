// ---------------------------------------------------------------------
// $Id: fe_q.cc 30037 2013-07-18 16:55:40Z maier $
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


#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/base/polynomials_bernstein.h>

#include <vector>
#include <sstream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
FE_Bernstein<dim,spacedim>::FE_Bernstein (const unsigned int degree)
  :
  FE_Q_Base<TensorProductPolynomials<dim>, dim, spacedim> (
    this->renumber_bases(degree),
    FiniteElementData<dim>(this->get_dpo_vector(degree),
                           1, degree,
                           FiniteElementData<dim>::H1),
    std::vector<bool> (1, false))
{}


template <int dim, int spacedim>
void
FE_Bernstein<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  Assert (dim > 1, ExcImpossibleInDim(1));
  get_subface_interpolation_matrix (source_fe, numbers::invalid_unsigned_int,
                                    interpolation_matrix);
}


template <int dim, int spacedim>
void
FE_Bernstein<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // see if source is a Bernstein element
  if (const FE_Bernstein<dim,spacedim> *source_fe
      = dynamic_cast<const FE_Bernstein<dim,spacedim> *>(&x_source_fe))
    {
      // have this test in here since a table of size 2x0 reports its size as
      // 0x0
      Assert (interpolation_matrix.n() == this->dofs_per_face,
              ExcDimensionMismatch (interpolation_matrix.n(),
                                    this->dofs_per_face));

      // Make sure that the element for which the DoFs should be constrained
      // is the one with the higher polynomial degree.  Actually the procedure
      // will work also if this assertion is not satisfied. But the matrices
      // produced in that case might lead to problems in the hp procedures,
      // which use this method.
      Assert (this->dofs_per_face <= source_fe->dofs_per_face,
              (typename FiniteElement<dim,spacedim>::
               ExcInterpolationNotImplemented ()));

      const Quadrature<dim-1>
      quad_face_support(FE_Q<dim,spacedim>(QIterated<1>(QTrapez<1>(),source_fe->degree)).get_unit_face_support_points ());

      // Rule of thumb for FP accuracy, that can be expected for a given
      // polynomial degree.  This value is used to cut off values close to
      // zero.
      double eps = 2e-13 * std::max(this->degree, source_fe->degree) * (dim-1);

      // compute the interpolation matrix by simply taking the value at the
      // support points.
//TODO: Verify that all faces are the same with respect to
// these support points. Furthermore, check if something has to
// be done for the face orientation flag in 3D.
      const Quadrature<dim> subface_quadrature
        = subface == numbers::invalid_unsigned_int
          ?
          QProjector<dim>::project_to_face (quad_face_support, 0)
          :
          QProjector<dim>::project_to_subface (quad_face_support, 0, subface);

      for (unsigned int i=0; i<source_fe->dofs_per_face; ++i)
        {
          const Point<dim> &p = subface_quadrature.point (i);
          for (unsigned int j=0; j<this->dofs_per_face; ++j)
            {
              double matrix_entry = this->shape_value (this->face_to_cell_index(j, 0), p);

              // Correct the interpolated value. I.e. if it is close to 1 or
              // 0, make it exactly 1 or 0. Unfortunately, this is required to
              // avoid problems with higher order elements.
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
FE_Bernstein<dim,spacedim>::hp_constraints_are_implemented () const
{
  return true;
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Bernstein<dim,spacedim>::hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  // we can presently only compute these identities if both FEs are FE_Bernsteins
  // or if the other one is an FE_Nothing. in the first case, there should be
  // exactly one single DoF of each FE at a vertex, and they should have
  // identical value
  if (dynamic_cast<const FE_Bernstein<dim,spacedim>*>(&fe_other) != 0)
    {
      return
        std::vector<std::pair<unsigned int, unsigned int> >
        (1, std::make_pair (0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no degrees of freedom, so there are no
      // equivalencies to be recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else if (fe_other.dofs_per_face == 0)
    {
      // if the other element has no elements on faces at all,
      // then it would be impossible to enforce any kind of
      // continuity even if we knew exactly what kind of element
      // we have -- simply because the other element declares
      // that it is discontinuous because it has no DoFs on
      // its faces. in that case, just state that we have no
      // constraints to declare
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Bernstein<dim,spacedim>::hp_line_dof_identities (const FiniteElement<dim,spacedim> &) const
{
  // Since this fe is not interpolatory but on the vertices, we can
  // not identify dofs on lines and on quads even if there are dofs
  // on lines and on quads.
  //
  // we also have nothing to say about interpolation to other finite
  // elements. consequently, we never have anything to say at all
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Bernstein<dim,spacedim>::hp_quad_dof_identities (const FiniteElement<dim,spacedim> &) const
{
  // Since this fe is not interpolatory but on the vertices, we can
  // not identify dofs on lines and on quads even if there are dofs
  // on lines and on quads.
  //
  // we also have nothing to say about interpolation to other finite
  // elements. consequently, we never have anything to say at all
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Bernstein<dim,spacedim>::compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Bernstein<dim,spacedim> *fe_b_other
      = dynamic_cast<const FE_Bernstein<dim,spacedim>*>(&fe_other))
    {
      if (this->degree < fe_b_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_b_other->degree)
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
std::string
FE_Bernstein<dim,spacedim>::get_name () const
{
  // note that the FETools::get_fe_from_name function depends on the
  // particular format of the string this function returns, so they have to be
  // kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Bernstein<" << dim << ">(" << this->degree << ")";
  return namebuf.str();
}


template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Bernstein<dim,spacedim>::clone() const
{
  return new FE_Bernstein<dim,spacedim>(*this);
}


/**
 * Only the assertion differs from the same function in FE_Q_Base!!
 */
template <int dim, int spacedim>
std::vector<unsigned int>
FE_Bernstein<dim,spacedim>::get_dpo_vector(const unsigned int deg)
{
  AssertThrow(deg>0,ExcMessage("FE_Bernstein needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}


template <int dim, int spacedim>
TensorProductPolynomials<dim>
FE_Bernstein<dim, spacedim>::renumber_bases(const unsigned int deg)
{
  TensorProductPolynomials<dim> tpp(dealii::generate_complete_bernstein_basis<double>(deg));
  std::vector<unsigned int> renumber(Utilities::fixed_power<dim>(deg+1));
  const FiniteElementData<dim> fe(this->get_dpo_vector(deg),1,
                                  deg);
  FETools::hierarchic_to_lexicographic_numbering (fe, renumber);
  tpp.set_numbering(renumber);
  return tpp;
}


// explicit instantiations
#include "fe_bernstein.inst"

DEAL_II_NAMESPACE_CLOSE
