// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_nothing.h>

#include <cmath>
#include <sstream>

//TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
//adjust_line_dof_index_for_line_orientation_table fields, and write tests
//similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN

namespace
{
  inline
  std::vector<unsigned int>
  invert_numbering (const std::vector<unsigned int> &in)
  {
    std::vector<unsigned int> out (in.size());
    for (unsigned int i=0; i<in.size(); ++i)
      out[in[i]]=i;
    return out;
  }
}



template <int dim>
FE_Q_Hierarchical<dim>::FE_Q_Hierarchical (const unsigned int degree)
  :
  FE_Poly<TensorProductPolynomials<dim>, dim> (
    Polynomials::Hierarchical::generate_complete_basis(degree),
    FiniteElementData<dim>(get_dpo_vector(degree),1, degree,
                           FiniteElementData<dim>::H1),
    std::vector<bool> (FiniteElementData<dim>(
                         get_dpo_vector(degree),1, degree).dofs_per_cell, false),
    std::vector<ComponentMask>(FiniteElementData<dim>(
                                 get_dpo_vector(degree),1, degree).dofs_per_cell, std::vector<bool>(1,true))),
  face_renumber(face_fe_q_hierarchical_to_hierarchic_numbering (degree))
{
  this->poly_space.set_numbering(
    hierarchic_to_fe_q_hierarchical_numbering(*this));

  // The matrix @p{dofs_cell} contains the
  // values of the linear functionals of
  // the master 1d cell applied to the
  // shape functions of the two 1d subcells.
  // The matrix @p{dofs_subcell} constains
  // the values of the linear functionals
  // on each 1d subcell applied to the
  // shape functions on the master 1d
  // subcell.
  // We use @p{dofs_cell} and
  // @p{dofs_subcell} to compute the
  // @p{prolongation}, @p{restriction} and
  // @p{interface_constraints} matrices
  // for all dimensions.
  std::vector<FullMatrix<double> >
  dofs_cell (GeometryInfo<1>::max_children_per_cell,
             FullMatrix<double> (2*this->dofs_per_vertex + this->dofs_per_line,
                                 2*this->dofs_per_vertex + this->dofs_per_line));
  std::vector<FullMatrix<double> >
  dofs_subcell (GeometryInfo<1>::max_children_per_cell,
                FullMatrix<double> (2*this->dofs_per_vertex + this->dofs_per_line,
                                    2*this->dofs_per_vertex + this->dofs_per_line));
  // build these fields, as they are
  // needed as auxiliary fields later
  // on
  build_dofs_cell (dofs_cell, dofs_subcell);

  // then use them to initialize
  // other fields
  initialize_constraints (dofs_subcell);
  initialize_embedding_and_restriction (dofs_cell, dofs_subcell);

  // finally fill in support points
  // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
}



template <int dim>
std::string
FE_Q_Hierarchical<dim>::get_name () const
{
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Q_Hierarchical<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_Q_Hierarchical<dim>::clone() const
{
  return new FE_Q_Hierarchical<dim>(*this);
}



template <int dim>
bool
FE_Q_Hierarchical<dim>::hp_constraints_are_implemented () const
{
  return true;
}


template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Q_Hierarchical<dim>::
hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const
{
  // we can presently only compute
  // these identities if both FEs are
  // FE_Q_Hierarchicals or if the other
  // one is an FE_Nothing. in the first
  // case, there should be exactly one
  // single DoF of each FE at a
  // vertex, and they should have
  // identical value
  if (dynamic_cast<const FE_Q_Hierarchical<dim>*>(&fe_other) != 0)
    {
      return
        std::vector<std::pair<unsigned int, unsigned int> >
        (1, std::make_pair (0U, 0U));
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
      // the FE_Nothing has no
      // degrees of freedom, so there
      // are no equivalencies to be
      // recorded
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> > ();
    }
}


template <int dim>
FiniteElementDomination::Domination
FE_Q_Hierarchical<dim>::
compare_for_face_domination (const FiniteElement<dim> &fe_other) const
{
  if (const FE_Q_Hierarchical<dim> *fe_q_other
      = dynamic_cast<const FE_Q_Hierarchical<dim>*>(&fe_other))
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
      // the FE_Nothing has no
      // degrees of freedom and it is
      // typically used in a context
      // where we don't require any
      // continuity along the
      // interface
      return FiniteElementDomination::no_requirements;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
void
FE_Q_Hierarchical<dim>::build_dofs_cell (std::vector<FullMatrix<double> > &dofs_cell,
                                         std::vector<FullMatrix<double> > &dofs_subcell) const
{
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;

  for (unsigned int c=0; c<GeometryInfo<1>::max_children_per_cell; ++c)
    for (unsigned int j=0; j<dofs_1d; ++j)
      for (unsigned int k=0; k<dofs_1d; ++k)
        {
          // upper diagonal block
          if ((j<=1) && (k<=1))
            {
              if (((c==0) && (j==0) && (k==0)) ||
                  ((c==1) && (j==1) && (k==1)))
                dofs_cell[c](j,k) = 1.;
              else
                dofs_cell[c](j,k) = 0.;

              if      (((c==0) && (j==1)) || ((c==1) && (j==0)))
                dofs_subcell[c](j,k) = .5;
              else if (((c==0) && (k==0)) || ((c==1) && (k==1)))
                dofs_subcell[c](j,k) = 1.;
              else
                dofs_subcell[c](j,k) = 0.;
            }
          // upper right block
          else if ((j<=1) && (k>=2))
            {
              if (((c==0) && (j==1) && ((k % 2)==0)) ||
                  ((c==1) && (j==0) && ((k % 2)==0)))
                dofs_subcell[c](j,k) = -1.;
            }
          // lower diagonal block
          else if ((j>=2) && (k>=2) && (j<=k))
            {
              double factor = 1.;
              for (unsigned int i=1; i<=j; ++i)
                factor *= ((double) (k-i+1))/((double) i);
              if (c==0)
                {
                  dofs_subcell[c](j,k) = ((k+j) % 2 == 0) ?
                                         std::pow(.5,static_cast<double>(k))*factor :
                                         -std::pow(.5,static_cast<double>(k))*factor;
                  dofs_cell[c](j,k) = std::pow(2.,static_cast<double>(j))*factor;
                }
              else
                {
                  dofs_subcell[c](j,k) = std::pow(.5,static_cast<double>(k))*factor;
                  dofs_cell[c](j,k) = ((k+j) % 2 == 0) ?
                                      std::pow(2.,static_cast<double>(j))*factor :
                                      -std::pow(2.,static_cast<double>(j))*factor;
                }
            }
        }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::
initialize_constraints (const std::vector<FullMatrix<double> > &dofs_subcell)
{
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;
  const unsigned int degree=this->degree;

  this->interface_constraints
  .TableBase<2,double>::reinit (this->interface_constraints_size());

  switch (dim)
    {
    case 1:
    {
      // no constraints in 1d
      break;
    }

    case 2:
    {
      // vertex node
      for (unsigned int i=0; i<dofs_1d; ++i)
        this->interface_constraints(0,i) = dofs_subcell[0](1,i);
      // edge nodes
      for (unsigned int c=0; c<GeometryInfo<1>::max_children_per_cell; ++c)
        for (unsigned int i=0; i<dofs_1d; ++i)
          for (unsigned int j=2; j<dofs_1d; ++j)
            this->interface_constraints(1 + c*(degree-1) + j - 2,i) =
              dofs_subcell[c](j,i);
      break;
    }

    case 3:
    {
      for (unsigned int i=0; i<dofs_1d * dofs_1d; i++)
        {
          // center vertex node
          this->interface_constraints(0,face_renumber[i]) =
            dofs_subcell[0](1,i % dofs_1d) *
            dofs_subcell[0](1,(i - (i % dofs_1d)) / dofs_1d);

          // boundary vertex nodes
          this->interface_constraints(1,face_renumber[i]) =
            dofs_subcell[0](0, i % dofs_1d) *
            dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
          this->interface_constraints(2,face_renumber[i]) =
            dofs_subcell[1](1, i % dofs_1d) *
            dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
          this->interface_constraints(3,face_renumber[i]) =
            dofs_subcell[0](1, i % dofs_1d) *
            dofs_subcell[0](0, (i - (i % dofs_1d)) / dofs_1d);
          this->interface_constraints(4,face_renumber[i]) =
            dofs_subcell[1](0, i % dofs_1d) *
            dofs_subcell[1](1, (i - (i % dofs_1d)) / dofs_1d);

          // interior edges
          for (unsigned int j=0; j<(degree-1); j++)
            {
              this->interface_constraints(5 + j,face_renumber[i]) =
                dofs_subcell[0](1, i % dofs_1d) *
                dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + (degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](1,i % dofs_1d) *
                dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 2*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](2 + j,i % dofs_1d) *
                dofs_subcell[1](0, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 3*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[1](2 + j, i % dofs_1d) *
                dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
            }

          // boundary edges
          for (unsigned int j=0; j<(degree-1); j++)
            {
              // left edge
              this->interface_constraints(5 + 4*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](0,     i % dofs_1d) *
                dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 4*(degree-1) + (degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](0,     i % dofs_1d) *
                dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              // right edge
              this->interface_constraints(5 + 4*(degree-1) + 2*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[1](1,     i % dofs_1d) *
                dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 4*(degree-1) + 3*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[1](1,     i % dofs_1d) *
                dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
              // bottom edge
              this->interface_constraints(5 + 4*(degree-1) + 4*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](2 + j, i % dofs_1d) *
                dofs_subcell[0](0,     (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 4*(degree-1) + 5*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[1](2 + j, i % dofs_1d) *
                dofs_subcell[0](0,     (i - (i % dofs_1d)) / dofs_1d);
              // top edge
              this->interface_constraints(5 + 4*(degree-1) + 6*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[0](2 + j, i % dofs_1d) *
                dofs_subcell[1](1,     (i - (i % dofs_1d)) / dofs_1d);
              this->interface_constraints(5 + 4*(degree-1) + 7*(degree-1) + j,face_renumber[i]) =
                dofs_subcell[1](2 + j, i % dofs_1d) *
                dofs_subcell[1](1,     (i - (i % dofs_1d)) / dofs_1d);
            }

          // interior faces
          for (unsigned int j=0; j<(degree-1); j++)
            for (unsigned int k=0; k<(degree-1); k++)
              {
                // subcell 0
                this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1),face_renumber[i]) =
                  dofs_subcell[0](2 + j, i % dofs_1d) *
                  dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                // subcell 1
                this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + (degree-1)*(degree-1),face_renumber[i]) =
                  dofs_subcell[1](2 + j, i % dofs_1d) *
                  dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                // subcell 2
                this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + 2*(degree-1)*(degree-1),face_renumber[i]) =
                  dofs_subcell[0](2 + j, i % dofs_1d) *
                  dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
                // subcell 3
                this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + 3*(degree-1)*(degree-1),face_renumber[i]) =
                  dofs_subcell[1](2 + j, i % dofs_1d) *
                  dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
              }
        }
      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::
initialize_embedding_and_restriction (const std::vector<FullMatrix<double> > &dofs_cell,
                                      const std::vector<FullMatrix<double> > &dofs_subcell)
{
  unsigned int iso=RefinementCase<dim>::isotropic_refinement-1;

  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;
  const std::vector<unsigned int> &renumber=
    this->poly_space.get_numbering();

  for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
    {
      this->prolongation[iso][c].reinit (this->dofs_per_cell, this->dofs_per_cell);
      this->restriction[iso][c].reinit (this->dofs_per_cell, this->dofs_per_cell);
    }

  // the 1d case is particularly
  // simple, so special case it:
  if (dim==1)
    {
      for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
        {
          this->prolongation[iso][c].fill (dofs_subcell[c]);
          this->restriction[iso][c].fill (dofs_cell[c]);
        }
      return;
    }

  // for higher dimensions, things
  // are a little more tricky:

  // j loops over dofs in the
  // subcell.  These are the rows in
  // the embedding matrix.
  //
  // i loops over the dofs in the
  // master cell. These are the
  // columns in the embedding matrix.
  for (unsigned int j=0; j<this->dofs_per_cell; ++j)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      switch (dim)
        {
        case 2:
        {
          for (unsigned int c=0; c<GeometryInfo<2>::max_children_per_cell; ++c)
            {
              // left/right line: 0/1
              const unsigned int c0 = c%2;
              // bottom/top line: 0/1
              const unsigned int c1 = c/2;

              this->prolongation[iso][c](j,i) =
                dofs_subcell[c0](renumber[j] % dofs_1d,
                                 renumber[i] % dofs_1d) *
                dofs_subcell[c1]((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d,
                                 (renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d);

              this->restriction[iso][c](j,i) =
                dofs_cell[c0](renumber[j] % dofs_1d,
                              renumber[i] % dofs_1d) *
                dofs_cell[c1]((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d,
                              (renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d);
            }
          break;
        }

        case 3:
        {
          for (unsigned int c=0; c<GeometryInfo<3>::max_children_per_cell; ++c)
            {
              // left/right face: 0/1
              const unsigned int c0 = c%2;
              // front/back face: 0/1
              const unsigned int c1 = (c%4)/2;
              // bottom/top face: 0/1
              const unsigned int c2 = c/4;

              this->prolongation[iso][c](j,i) =
                dofs_subcell[c0](renumber[j] % dofs_1d,
                                 renumber[i] % dofs_1d) *
                dofs_subcell[c1](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) % dofs_1d,
                                 ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
                dofs_subcell[c2](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d - (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
                                 ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d - (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);

              this->restriction[iso][c](j,i) =
                dofs_cell[c0](renumber[j] % dofs_1d,
                              renumber[i] % dofs_1d) *
                dofs_cell[c1](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) % dofs_1d,
                              ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
                dofs_cell[c2](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d - (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
                              ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d - (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);
            }
          break;
        }

        default:
          Assert (false, ExcNotImplemented());
        }
}



template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_support_points ()
{
  // number of points: (degree+1)^dim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= this->degree+1;

  this->unit_support_points.resize(n);

  const std::vector<unsigned int> &index_map_inverse=
    this->poly_space.get_numbering_inverse();

  Point<dim> p;
  // the method of numbering allows
  // each dof to be associated with a
  // support point. There is
  // only one support point per
  // vertex, line, quad, hex, etc.
  //
  // note, however, that the support
  // points thus associated with
  // shape functions are not unique:
  // the linear shape functions are
  // associated with the vertices,
  // but all others are associated
  // with either line, quad, or hex
  // midpoints, and there may be
  // multiple shape functions
  // associated with them. there
  // really is no other useful
  // numbering, since the
  // hierarchical shape functions do
  // not vanish at all-but-one
  // interpolation points (like the
  // Lagrange functions used in
  // FE_Q), so there's not much we
  // can do here.
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((dim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
        {
          if (ix==0)
            p(0) =  0.;
          else if (ix==1)
            p(0) =  1.;
          else
            p(0) = .5;
          if (dim>1)
            {
              if (iy==0)
                p(1) =  0.;
              else if (iy==1)
                p(1) =  1.;
              else
                p(1) = .5;
            }
          if (dim>2)
            {
              if (iz==0)
                p(2) =  0.;
              else if (iz==1)
                p(2) =  1.;
              else
                p(2) = .5;
            }
          this->unit_support_points[index_map_inverse[k++]] = p;
        };
}



template <>
void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ()
{
  // no faces in 1d, so nothing to do
}


template <>
void FE_Q_Hierarchical<1>::
get_face_interpolation_matrix (const FiniteElement<1,1> &/*x_source_fe*/,
                               FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false, ExcImpossibleInDim(1));
}


template <>
void
FE_Q_Hierarchical<1>::
get_subface_interpolation_matrix (const FiniteElement<1,1> &/*x_source_fe*/,
                                  const unsigned int      /*subface*/,
                                  FullMatrix<double>     &/*interpolation_matrix*/) const
{
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim>
void
FE_Q_Hierarchical<dim>::
get_face_interpolation_matrix (const FiniteElement<dim> &x_source_fe,
                               FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // Q_Hierarchical element
  typedef FE_Q_Hierarchical<dim> FEQHierarchical;
  typedef FiniteElement<dim> FEL;
  AssertThrow ((x_source_fe.get_name().find ("FE_Q_Hierarchical<") == 0)
               ||
               (dynamic_cast<const FEQHierarchical *>(&x_source_fe) != 0),
               typename FEL::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // ok, source is a Q_Hierarchical element, so
  // we will be able to do the work
  const FE_Q_Hierarchical<dim> &source_fe
    = dynamic_cast<const FE_Q_Hierarchical<dim>&>(x_source_fe);

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp procedures, which use this
  // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FEL::
          ExcInterpolationNotImplemented ());
  interpolation_matrix = 0;

  switch (dim)
    {
    case 2:
    {
      for (unsigned int i = 0; i < this->dofs_per_face; ++i)
        interpolation_matrix (i, i) = 1;

      break;
    }

    case 3:
    {
      for (unsigned int i = 0; i < GeometryInfo<3>::vertices_per_face; ++i)
        interpolation_matrix (i, i) = 1;

      for (unsigned int i = 0; i < this->degree - 1; ++i)
        {
          for (unsigned int j = 0; j < GeometryInfo<3>::lines_per_face; ++j)
            interpolation_matrix (
              i + j * (x_source_fe.degree - 1) + GeometryInfo<3>::vertices_per_face,
              i + j * (this->degree - 1) + GeometryInfo<3>::vertices_per_face) = 1;

          for (unsigned int j = 0; j < this->degree - 1; ++j)
            interpolation_matrix (
              (i + GeometryInfo<3>::lines_per_face) * (x_source_fe.degree - 1) + j
              + GeometryInfo<3>::vertices_per_face,
              (i + GeometryInfo<3>::lines_per_face) * (this->degree - 1) + j
              + GeometryInfo<3>::vertices_per_face) = 1;
        }
    }
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::
get_subface_interpolation_matrix (const FiniteElement<dim> &x_source_fe,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &interpolation_matrix) const
{
  // this is only implemented, if the
  // source FE is also a
  // Q_Hierarchical element
  typedef FE_Q_Hierarchical<dim> FEQHierarchical;
  typedef FiniteElement<dim> FEL;
  AssertThrow ((x_source_fe.get_name().find ("FE_Q_Hierarchical<") == 0)
               ||
               (dynamic_cast<const FEQHierarchical *>(&x_source_fe) != 0),
               typename FEL::
               ExcInterpolationNotImplemented());

  Assert (interpolation_matrix.n() == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n(),
                                this->dofs_per_face));
  Assert (interpolation_matrix.m() == x_source_fe.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m(),
                                x_source_fe.dofs_per_face));

  // ok, source is a Q_Hierarchical element, so
  // we will be able to do the work
  const FE_Q_Hierarchical<dim> &source_fe
    = dynamic_cast<const FE_Q_Hierarchical<dim>&>(x_source_fe);

  // Make sure, that the element,
  // for which the DoFs should be
  // constrained is the one with
  // the higher polynomial degree.
  // Actually the procedure will work
  // also if this assertion is not
  // satisfied. But the matrices
  // produced in that case might
  // lead to problems in the
  // hp procedures, which use this
  // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FEL::
          ExcInterpolationNotImplemented ());

  switch (dim)
    {
    case 2:
    {
      switch (subface)
        {
        case 0:
        {
          interpolation_matrix (0, 0) = 1.0;
          interpolation_matrix (1, 0) = 0.5;
          interpolation_matrix (1, 1) = 0.5;

          for (unsigned int dof = 2; dof < this->dofs_per_face;)
            {
              interpolation_matrix (1, dof) = -1.0;
              dof = dof + 2;
            }

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;

          for (int i = 2; i < (int) this->dofs_per_face; ++i)
            {
              interpolation_matrix (i, i) = std::pow (0.5, i);
              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j < (int) this->dofs_per_face; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;

                  if ((i + j) & 1)
                    interpolation_matrix (i, j)
                      = -1.0 * std::pow (0.5, j) *
                        factorial_j / (factorial_i * factorial_ij);

                  else
                    interpolation_matrix (i, j)
                      = std::pow (0.5, j) *
                        factorial_j / (factorial_i * factorial_ij);
                }
            }

          break;
        }

        case 1:
        {
          interpolation_matrix (0, 0) = 0.5;
          interpolation_matrix (0, 1) = 0.5;

          for (unsigned int dof = 2; dof < this->dofs_per_face;)
            {
              interpolation_matrix (0, dof) = -1.0;
              dof = dof + 2;
            }

          interpolation_matrix (1, 1) = 1.0;

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;

          for (int i = 2; i < (int) this->dofs_per_face; ++i)
            {
              interpolation_matrix (i, i) = std::pow (0.5, i);
              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j < (int) this->dofs_per_face; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;
                  interpolation_matrix (i, j)
                    = std::pow (0.5, j) * factorial_j / (factorial_i * factorial_ij);
                }
            }
        }
        }

      break;
    }

    case 3:
    {
      switch (subface)
        {
        case 0:
        {
          interpolation_matrix (0, 0) = 1.0;
          interpolation_matrix (1, 0) = 0.5;
          interpolation_matrix (1, 1) = 0.5;
          interpolation_matrix (2, 0) = 0.5;
          interpolation_matrix (2, 2) = 0.5;

          for (unsigned int i = 0; i < this->degree - 1;)
            {
              for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_face; ++line)
                interpolation_matrix (3, i + line * (this->degree - 1) + 4) = -0.5;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (3, i + (j + 4) * this->degree - j) = 1.0;
                  j = j + 2;
                }

              interpolation_matrix (1, i + 2 * (this->degree + 1)) = -1.0;
              interpolation_matrix (2, i + 4) = -1.0;
              i = i + 2;
            }

          for (unsigned int vertex = 0; vertex < GeometryInfo<3>::vertices_per_face; ++vertex)
            interpolation_matrix (3, vertex) = 0.25;

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;
          int factorial_k;
          int factorial_kl;
          int factorial_l;

          for (int i = 2; i <= (int) this->degree; ++i)
            {
              double tmp = std::pow (0.5, i);
              interpolation_matrix (i + 2, i + 2) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 2 * this->degree) = tmp;
              tmp *= 0.5;
              interpolation_matrix (i + source_fe.degree + 1, i + 2) = tmp;
              interpolation_matrix (i + source_fe.degree + 1, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 2 * this->degree) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 3 * this->degree - 1) = tmp;
              tmp *= -2.0;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (i + source_fe.degree + 1, (i + 2) * this->degree + j + 2 - i) = tmp;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, i + (j + 4) * this->degree - j - 2) = tmp;
                  j = j + 2;
                }

              factorial_k = 1;

              for (int j = 2; j <= (int) this->degree; ++j)
                {
                  interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (j + 2) * this->degree - j) = std::pow (0.5, i + j);
                  factorial_k *= j;
                  factorial_kl = 1;
                  factorial_l = factorial_k;

                  for (int k = j + 1; k < (int) this->degree; ++k)
                    {
                      factorial_kl *= k - j;
                      factorial_l *= k;

                      if ((j + k) & 1)
                        interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = -1.0 * std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);

                      else
                        interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);
                    }
                }

              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j <= (int) this->degree; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;

                  if ((i + j) & 1)
                    {
                      tmp = -1.0 * std::pow (0.5, j) * factorial_j / (factorial_i * factorial_ij);
                      interpolation_matrix (i + 2, j + 2) = tmp;
                      interpolation_matrix (i + 2 * source_fe.degree, j + 2 * this->degree) = tmp;
                      factorial_k = 1;

                      for (int k = 2; k <= (int) this->degree; ++k)
                        {
                          interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (k + 2) * this->degree - k) = tmp * std::pow (0.5, k);
                          factorial_k *= k;
                          factorial_l = factorial_k;
                          factorial_kl = 1;

                          for (int l = k + 1; l <= (int) this->degree; ++l)
                            {
                              factorial_kl *= l - k;
                              factorial_l *= l;

                              if ((k + l) & 1)
                                interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = -1.0 * tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);

                              else
                                interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);
                            }
                        }

                      tmp *= 0.5;
                      interpolation_matrix (i + source_fe.degree + 1, j + 2) = tmp;
                      interpolation_matrix (i + source_fe.degree + 1, j + this->degree + 1) = tmp;
                      interpolation_matrix (i + 3 * source_fe.degree - 1, j + 2 * this->degree) = tmp;
                      interpolation_matrix (i + 3 * source_fe.degree - 1, j + 3 * this->degree - 1) = tmp;
                      tmp *= -2.0;

                      for (unsigned int k = 0; k < this->degree - 1;)
                        {
                          interpolation_matrix (i + source_fe.degree + 1, (j + 2) * this->degree + k + 2 - j) = tmp;
                          interpolation_matrix (i + 3 * source_fe.degree - 1, j + (k + 4) * this->degree - k - 2) = tmp;
                          k = k + 2;
                        }
                    }
                  else
                    {
                      tmp = std::pow (0.5, j) * factorial_j / (factorial_i * factorial_ij);
                      interpolation_matrix (i + 2, j + 2) = tmp;
                      interpolation_matrix (i + 2 * source_fe.degree, j + 2 * this->degree) = tmp;
                      factorial_k = 1;

                      for (int k = 2; k <= (int) this->degree; ++k)
                        {
                          interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (k + 2) * this->degree - k) = tmp * std::pow (0.5, k);
                          factorial_k *= k;
                          factorial_l = factorial_k;
                          factorial_kl = 1;

                          for (int l = k + 1; l <= (int) this->degree; ++l)
                            {
                              factorial_kl *= l - k;
                              factorial_l *= l;

                              if ((k + l) & 1)
                                interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = -1.0 * tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);

                              else
                                interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);
                            }
                        }

                      tmp *= 0.5;
                      interpolation_matrix (i + source_fe.degree + 1, j + 2) = tmp;
                      interpolation_matrix (i + source_fe.degree + 1, j + this->degree + 1) = tmp;
                      interpolation_matrix (i + 3 * source_fe.degree - 1, j + 2 * this->degree) = tmp;
                      interpolation_matrix (i + 3 * source_fe.degree - 1, j + 3 * this->degree - 1) = tmp;
                      tmp *= -2.0;

                      for (unsigned int k = 0; k < this->degree - 1;)
                        {
                          interpolation_matrix (i + source_fe.degree + 1, (j + 2) * this->degree + k + 2 - j) = tmp;
                          interpolation_matrix (i + 3 * source_fe.degree - 1, j + (k + 4) * this->degree - k - 2) = tmp;
                          k = k + 2;
                        }
                    }
                }
            }

          break;
        }

        case 1:
        {
          interpolation_matrix (0, 0) = 0.5;
          interpolation_matrix (0, 1) = 0.5;
          interpolation_matrix (1, 1) = 1.0;
          interpolation_matrix (3, 1) = 0.5;
          interpolation_matrix (3, 3) = 0.5;

          for (unsigned int i = 0; i < this->degree - 1;)
            {
              for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_face; ++line)
                interpolation_matrix (2, i + line * (this->degree - 1) + 4) = -0.5;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (2, i + (j + 4) * this->degree - j) = 1.0;
                  j = j + 2;
                }

              interpolation_matrix (0, i + 2 * (this->degree + 1)) = -1.0;
              interpolation_matrix (3, i + this->degree + 3) = -1.0;
              i = i + 2;
            }

          for (unsigned int vertex = 0; vertex < GeometryInfo<3>::vertices_per_face; ++vertex)
            interpolation_matrix (2, vertex) = 0.25;

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;
          int factorial_k;
          int factorial_kl;
          int factorial_l;

          for (int i = 2; i <= (int) this->degree; ++i)
            {
              double tmp = std::pow (0.5, i + 1);
              interpolation_matrix (i + 2, i + 2) = tmp;
              interpolation_matrix (i + 2, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 2 * this->degree) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 3 * this->degree - 1) = tmp;
              tmp *= -2.0;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (i + 2, j + (i + 2) * this->degree + 2 - i) = tmp;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, i + (j + 4) * this->degree - j - 2) = tmp;
                  j = j + 2;
                }

              tmp *= - 1.0;
              interpolation_matrix (i + source_fe.degree + 1, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 2 * this->degree) = tmp;
              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j <= (int) this->degree; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;
                  tmp = std::pow (0.5, j) * factorial_j / (factorial_i * factorial_ij);
                  interpolation_matrix (i + 2 * source_fe.degree, j + 2 * this->degree) = tmp;
                  factorial_k = 1;

                  for (int k = 2; k <= (int) this->degree; ++k)
                    {
                      interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (k + 2) * this->degree - k) = tmp * std::pow (0.5, k);
                      factorial_k *= k;
                      factorial_l = factorial_k;
                      factorial_kl = 1;

                      for (int l = k + 1; l <= (int) this->degree; ++l)
                        {
                          factorial_kl *= l - k;
                          factorial_l *= l;

                          if ((k + l) & 1)
                            interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = -1.0 * tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);

                          else
                            interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);
                        }
                    }

                  tmp *= -1.0;

                  for (unsigned int k = 0; k < this->degree - 1;)
                    {
                      interpolation_matrix (i + 3 * source_fe.degree - 1, j + (k + 4) * this->degree - k - 2) = tmp;
                      k = k + 2;
                    }

                  tmp *= -0.5;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, j + 2 * this->degree) = tmp;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, j + 3 * this->degree - 1) = tmp;

                  if ((i + j) & 1)
                    tmp *= -1.0;

                  interpolation_matrix (i + 2, j + 2) = tmp;
                  interpolation_matrix (i + 2, j + this->degree + 1) = tmp;
                  interpolation_matrix (i + source_fe.degree + 1, j + this->degree + 1) = 2.0 * tmp;
                  tmp *= -2.0;

                  for (unsigned int k = 0; k < this->degree - 1;)
                    {
                      interpolation_matrix (i + 2, k + (j + 2) * this->degree + 2 - j) = tmp;
                      k = k + 2;
                    }
                }

              factorial_k = 1;

              for (int j = 2; j <= (int) this->degree; ++j)
                {
                  interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (j + 2) * this->degree - j) = std::pow (0.5, i + j);
                  factorial_k *= j;
                  factorial_l = factorial_k;
                  factorial_kl = 1;

                  for (int k = j + 1; k <= (int) this->degree; ++k)
                    {
                      factorial_kl *= k - j;
                      factorial_l *= k;

                      if ((j + k) & 1)
                        interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = -1.0 * std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);

                      else
                        interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);
                    }
                }
            }

          break;
        }

        case 2:
        {
          interpolation_matrix (0, 0) = 0.5;
          interpolation_matrix (0, 2) = 0.5;
          interpolation_matrix (2, 2) = 1.0;
          interpolation_matrix (3, 2) = 0.5;
          interpolation_matrix (3, 3) = 0.5;

          for (unsigned int i = 0; i < this->degree - 1;)
            {
              for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_face; ++line)
                interpolation_matrix (1, i + line * (this->degree - 1) + 4) = -0.5;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (1, i + (j + 4) * this->degree - j) = 1.0;
                  j = j + 2;
                }

              interpolation_matrix (0, i + 4) = -1.0;
              interpolation_matrix (3, i + 3 * this->degree + 1) = -1.0;
              i = i + 2;
            }

          for (unsigned int vertex = 0; vertex < GeometryInfo<3>::vertices_per_face; ++vertex)
            interpolation_matrix (1, vertex) = 0.25;

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;
          int factorial_k;
          int factorial_kl;
          int factorial_l;

          for (int i = 2; i <= (int) this->degree; ++i)
            {
              double tmp = std::pow (0.5, i);
              interpolation_matrix (i + 2, i + 2) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 3 * this->degree - 1) = tmp;
              tmp *= 0.5;
              interpolation_matrix (i + source_fe.degree + 1, i + 2) = tmp;
              interpolation_matrix (i + source_fe.degree + 1, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 2 * this->degree) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 3 * this->degree - 1) = tmp;
              tmp *= -2.0;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (i + source_fe.degree + 1, j + (i + 2) * this->degree + 2 - i) = tmp;
                  interpolation_matrix (i + 2 * source_fe.degree, i + (j + 4) * this->degree - j - 2) = tmp;
                  j = j + 2;
                }

              factorial_k = 1;

              for (int j = 2; j <= (int) this->degree; ++j)
                {
                  interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (j + 2) * this->degree - j) = std::pow (0.5, i + j);
                  factorial_k *= j;
                  factorial_l = factorial_k;
                  factorial_kl = 1;

                  for (int k = j + 1; k <= (int) this->degree; ++k)
                    {
                      factorial_kl *= k - j;
                      factorial_l *= k;
                      interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);
                    }
                }

              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j <= (int) this->degree; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;
                  tmp = std::pow (0.5, j) * factorial_j / (factorial_i * factorial_ij);
                  interpolation_matrix (i + 2, j + 2) = tmp;
                  tmp *= -1.0;

                  for (unsigned int k = 0; k < this->degree - 1;)
                    {
                      interpolation_matrix (i + source_fe.degree + 1, k + (j + 2) * this->degree + 2 - j) = tmp;
                      k = k + 2;
                    }

                  tmp *= -0.5;
                  interpolation_matrix (i + source_fe.degree + 1, j + 2) = tmp;
                  interpolation_matrix (i + source_fe.degree + 1, j + this->degree + 1) = tmp;

                  if ((i + j) & 1)
                    tmp *= -1.0;

                  interpolation_matrix (i + 2 * source_fe.degree, j + 2 * this->degree) = tmp;
                  interpolation_matrix (i + 2 * source_fe.degree, j + 3 * this->degree - 1) = tmp;
                  tmp *= 2.0;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, j + 3 * this->degree - 1) = tmp;
                  factorial_k = 1;

                  for (int k = 2; k <= (int) this->degree; ++k)
                    {
                      interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (k + 2) * this->degree - k) = tmp * std::pow (0.5, k);
                      factorial_k *= k;
                      factorial_l = factorial_k;
                      factorial_kl = 1;

                      for (int l = k + 1; l <= (int) this->degree; ++l)
                        {
                          factorial_kl *= l - k;
                          factorial_l *= l;
                          interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);
                        }
                    }

                  tmp *= -1.0;

                  for (unsigned int k = 0; k < this->degree - 1;)
                    {
                      interpolation_matrix (i + 2 * source_fe.degree, j + (k + 4) * this->degree - k - 2) = tmp;
                      k = k + 2;
                    }
                }
            }

          break;
        }

        case 3:
        {
          for (unsigned int vertex = 0; vertex < GeometryInfo<3>::vertices_per_face; ++vertex)
            interpolation_matrix (0, vertex) = 0.25;

          for (unsigned int i = 0; i < this->degree - 1;)
            {
              for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_face; ++line)
                interpolation_matrix (0, i + line * (this->degree - 1) + 4) = -0.5;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (0, i + (j + 4) * this->degree - j) = 1.0;
                  j = j + 2;
                }

              interpolation_matrix (1, i + 4) = -1.0;
              interpolation_matrix (2, i + 3 * this->degree + 1) = -1.0;
              i = i + 2;
            }

          interpolation_matrix (1, 0) = 0.5;
          interpolation_matrix (1, 1) = 0.5;
          interpolation_matrix (2, 2) = 0.5;
          interpolation_matrix (2, 3) = 0.5;
          interpolation_matrix (3, 3) = 1.0;

          int factorial_i = 1;
          int factorial_ij;
          int factorial_j;
          int factorial_k;
          int factorial_kl;
          int factorial_l;

          for (int i = 2; i <= (int) this->degree; ++i)
            {
              double tmp = std::pow (0.5, i + 1);
              interpolation_matrix (i + 2, i + 2) = tmp;
              interpolation_matrix (i + 2, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 2 * this->degree) = tmp;
              interpolation_matrix (i + 2 * source_fe.degree, i + 3 * this->degree - 1) = tmp;
              tmp *= -2.0;

              for (unsigned int j = 0; j < this->degree - 1;)
                {
                  interpolation_matrix (i + 2, j + (i + 2) * this->degree + 2 - i) = tmp;
                  interpolation_matrix (i + 2 * source_fe.degree, i + (j + 4) * this->degree - 2) = tmp;
                  j = j + 2;
                }

              tmp *= -1.0;
              interpolation_matrix (i + source_fe.degree + 1, i + this->degree + 1) = tmp;
              interpolation_matrix (i + 3 * source_fe.degree - 1, i + 3 * this->degree - 1) = tmp;
              factorial_k = 1;

              for (int j = 2; j <= (int) this->degree; ++j)
                {
                  interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (j + 2) * this->degree - j) = std::pow (0.5, i + j);
                  factorial_k *= j;
                  factorial_l = factorial_k;
                  factorial_kl = 1;

                  for (int k = j + 1; k <= (int) this->degree; ++k)
                    {
                      factorial_kl *= k - j;
                      factorial_l *= k;
                      interpolation_matrix (i + (j + 2) * source_fe.degree - j, i + (k + 2) * this->degree - k) = std::pow (0.5, i + k) * factorial_l / (factorial_k * factorial_kl);
                    }
                }

              factorial_i *= i;
              factorial_j = factorial_i;
              factorial_ij = 1;

              for (int j = i + 1; j <= (int) this->degree; ++j)
                {
                  factorial_ij *= j - i;
                  factorial_j *= j;
                  tmp = std::pow (0.5, j + 1) * factorial_j / (factorial_i * factorial_ij);
                  interpolation_matrix (i + 2, j + 2) = tmp;
                  interpolation_matrix (i + 2, j + this->degree + 1) = tmp;
                  interpolation_matrix (i + 2 * source_fe.degree, j + 2 * this->degree) = tmp;
                  interpolation_matrix (i + 2 * source_fe.degree, j + 3 * this->degree - 1) = tmp;
                  tmp *= 2.0;
                  interpolation_matrix (i + source_fe.degree + 1, j + this->degree + 1) = tmp;
                  interpolation_matrix (i + 3 * source_fe.degree - 1, j + 3 * this->degree - 1) = tmp;
                  factorial_k = 1;

                  for (int k = 2; k <= (int) this->degree; ++k)
                    {
                      interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (k + 2) * this->degree - k) = tmp * std::pow (0.5, k);
                      factorial_k *= k;
                      factorial_l = factorial_k;
                      factorial_kl = 1;

                      for (int l = k + 1; l <= (int) this->degree; ++l)
                        {
                          factorial_kl *= l - k;
                          factorial_l *= l;
                          interpolation_matrix (i + (k + 2) * source_fe.degree - k, j + (l + 2) * this->degree - l) = tmp * std::pow (0.5, l) * factorial_l / (factorial_k * factorial_kl);
                        }
                    }

                  tmp *= -1.0;

                  for (unsigned int k = 0; k < this->degree - 1;)
                    {
                      interpolation_matrix (i + 2, k + (j + 2) * this->degree + 2 - j) = tmp;
                      interpolation_matrix (i + 2 * source_fe.degree, j + (k + 4) * this->degree - 2) = tmp;
                      k = k + 2;
                    }
                }
            }
        }
        }
    }
    }
}



template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_face_support_points ()
{
  const unsigned int codim = dim-1;

  // number of points: (degree+1)^codim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= this->degree+1;

  this->unit_face_support_points.resize(n);

  Point<codim> p;

  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
        {
          if (ix==0)
            p(0) =  0.;
          else if (ix==1)
            p(0) =  1.;
          else
            p(0) = .5;
          if (codim>1)
            {
              if (iy==0)
                p(1) =  0.;
              else if (iy==1)
                p(1) =  1.;
              else
                p(1) = .5;
            }
          if (codim>2)
            {
              if (iz==0)
                p(2) =  0.;
              else if (iz==1)
                p(2) =  1.;
              else
                p(2) = .5;
            }
          this->unit_face_support_points[face_renumber[k++]] = p;
        };
}


// we use same dpo_vector as FE_Q
template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::
hierarchic_to_fe_q_hierarchical_numbering (const FiniteElementData<dim> &fe)
{
  Assert (fe.n_components() == 1, ExcInternalError());
  std::vector<unsigned int> h2l(fe.dofs_per_cell);

  // polynomial degree
  const unsigned int degree = fe.dofs_per_line+1;
  // number of grid points in each
  // direction
  const unsigned int n = degree+1;

  // the following lines of code are
  // somewhat odd, due to the way the
  // hierarchic numbering is
  // organized. if someone would
  // really want to understand these
  // lines, you better draw some
  // pictures where you indicate the
  // indices and orders of vertices,
  // lines, etc, along with the
  // numbers of the degrees of
  // freedom in hierarchical and
  // lexicographical order
  switch (dim)
    {
    case 1:
    {
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        h2l[i] = i;

      break;
    }

    case 2:
    {
      // Example: degree=3
      //
      // hierarchical numbering:
      //  2 10 11  3
      //  5 14 15  7
      //  4 12 13  6
      //  0  8  9  1
      //
      // fe_q_hierarchical numbering:
      //  4  6  7  5
      // 12 14 15 13
      //  8 10 11  9
      //  0  2  3  1
      unsigned int next_index = 0;
      // first the four vertices
      h2l[next_index++] = 0;
      h2l[next_index++] = 1;
      h2l[next_index++] = n;
      h2l[next_index++] = n+1;
      // left line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n;
      // right line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n+1;
      // bottom line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = 2+i;
      // top line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n+2+i;
      // inside quad
      Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
              ExcInternalError());
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n+2+j;

      Assert (next_index == fe.dofs_per_cell, ExcInternalError());

      break;
    }

    case 3:
    {
      unsigned int next_index = 0;
      const unsigned int n2=n*n;
      // first the eight vertices
      // bottom face, lexicographic
      h2l[next_index++] = 0;
      h2l[next_index++] = 1;
      h2l[next_index++] = n;
      h2l[next_index++] = n+1;
      // top face, lexicographic
      h2l[next_index++] = n2;
      h2l[next_index++] = n2+1;
      h2l[next_index++] = n2+n;
      h2l[next_index++] = n2+n+1;

      // now the lines
      // bottom face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n+1;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = 2+i;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n+2+i;
      // top face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n2+(2+i)*n;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n2+(2+i)*n+1;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n2+2+i;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = n2+n+2+i;
      // lines in z-direction
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n2;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n2+1;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n2+n;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        h2l[next_index++] = (2+i)*n2+n+1;

      // inside quads
      Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
              ExcInternalError());
      // left face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n2+(2+j)*n;
      // right face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n2+(2+j)*n+1;
      // front face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n2+2+j;
      // back face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n2+n+2+j;
      // bottom face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = (2+i)*n+2+j;
      // top face
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          h2l[next_index++] = n2+(2+i)*n+2+j;

      // inside hex
      Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
              ExcInternalError());
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          for (unsigned int k=0; k<fe.dofs_per_line; ++k)
            h2l[next_index++] = (2+i)*n2+(2+j)*n+2+k;

      Assert (next_index == fe.dofs_per_cell, ExcInternalError());

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }
  return h2l;
}


template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::
face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int degree)
{
  FiniteElementData<dim-1> fe_data(FE_Q_Hierarchical<dim-1>::get_dpo_vector(degree),1,degree);
  return invert_numbering(FE_Q_Hierarchical<dim-1>::
                          hierarchic_to_fe_q_hierarchical_numbering (fe_data));
}



template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int)
{
  return std::vector<unsigned int> ();
}


template <>
bool
FE_Q_Hierarchical<1>::has_support_on_face (const unsigned int shape_index,
                                           const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<1>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<1>::faces_per_cell));


  // in 1d, things are simple. since
  // there is only one degree of
  // freedom per vertex in this
  // class, the first is on vertex 0
  // (==face 0 in some sense), the
  // second on face 1:
  return (((shape_index == 0) && (face_index == 0)) ||
          ((shape_index == 1) && (face_index == 1)));
}




template <int dim>
bool
FE_Q_Hierarchical<dim>::has_support_on_face (const unsigned int shape_index,
                                             const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  // first, special-case interior
  // shape functions, since they
  // have no support no-where on
  // the boundary
  if (((dim==2) && (shape_index>=this->first_quad_index))
      ||
      ((dim==3) && (shape_index>=this->first_hex_index)))
    return false;

  // let's see whether this is a
  // vertex
  if (shape_index < this->first_line_index)
    {
      // for Q elements, there is
      // one dof per vertex, so
      // shape_index==vertex_number. check
      // whether this vertex is
      // on the given face.
      const unsigned int vertex_no = shape_index;
      Assert (vertex_no < GeometryInfo<dim>::vertices_per_cell,
              ExcInternalError());
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
        if (GeometryInfo<dim>::face_to_cell_vertices(face_index,i) == vertex_no)
          return true;
      return false;
    }
  else if (shape_index < this->first_quad_index)
    // ok, dof is on a line
    {
      const unsigned int line_index
        = (shape_index - this->first_line_index) / this->dofs_per_line;
      Assert (line_index < GeometryInfo<dim>::lines_per_cell,
              ExcInternalError());

      for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_face; ++i)
        if (GeometryInfo<dim>::face_to_cell_lines(face_index,i) == line_index)
          return true;
      return false;
    }
  else if (shape_index < this->first_hex_index)
    // dof is on a quad
    {
      const unsigned int quad_index
        = (shape_index - this->first_quad_index) / this->dofs_per_quad;
      Assert (static_cast<signed int>(quad_index) <
              static_cast<signed int>(GeometryInfo<dim>::quads_per_cell),
              ExcInternalError());

      // in 2d, cell bubble are
      // zero on all faces. but
      // we have treated this
      // case above already
      Assert (dim != 2, ExcInternalError());

      // in 3d,
      // quad_index=face_index
      if (dim == 3)
        return (quad_index == face_index);
      else
        Assert (false, ExcNotImplemented());
    }
  else
    // dof on hex
    {
      // can only happen in 3d, but
      // this case has already been
      // covered above
      Assert (false, ExcNotImplemented());
      return false;
    }

  // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;
}



template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_embedding_dofs (const unsigned int sub_degree) const
{
  Assert ((sub_degree>0) && (sub_degree<=this->degree),
          ExcIndexRange(sub_degree, 1, this->degree));

  if (dim==1)
    {
      std::vector<unsigned int> embedding_dofs (sub_degree+1);
      for (unsigned int i=0; i<(sub_degree+1); ++i)
        embedding_dofs[i] = i;

      return embedding_dofs;
    }

  if (sub_degree==1)
    {
      std::vector<unsigned int> embedding_dofs (GeometryInfo<dim>::vertices_per_cell);
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        embedding_dofs[i] = i;

      return embedding_dofs;
    }
  else if (sub_degree==this->degree)
    {
      std::vector<unsigned int> embedding_dofs (this->dofs_per_cell);
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        embedding_dofs[i] = i;

      return embedding_dofs;
    }

  if ((dim==2) || (dim==3))
    {
      std::vector<unsigned int> embedding_dofs ( (dim==2) ?
                                                 (sub_degree+1) * (sub_degree+1) :
                                                 (sub_degree+1) * (sub_degree+1) * (sub_degree+1) );

      for (unsigned int i=0; i<( (dim==2) ?
                                 (sub_degree+1) * (sub_degree+1) :
                                 (sub_degree+1) * (sub_degree+1) * (sub_degree+1) ); ++i)
        {
          // vertex
          if (i<GeometryInfo<dim>::vertices_per_cell)
            embedding_dofs[i] = i;
          // line
          else if (i<(GeometryInfo<dim>::vertices_per_cell +
                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1)))
            {
              const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell) %
                                     (sub_degree-1);
              const unsigned int line = (i - GeometryInfo<dim>::vertices_per_cell - j) / (sub_degree-1);

              embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell +
                                  line * (this->degree-1) + j;
            }
          // quad
          else if (i<(GeometryInfo<dim>::vertices_per_cell +
                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) +
                   GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1))
            {
              const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell -
                                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) % (sub_degree-1);
              const unsigned int k = ( (i - GeometryInfo<dim>::vertices_per_cell -
                                        GeometryInfo<dim>::lines_per_cell * (sub_degree-1) - j) / (sub_degree-1) ) % (sub_degree-1);
              const unsigned int face = (i - GeometryInfo<dim>::vertices_per_cell -
                                         GeometryInfo<dim>::lines_per_cell * (sub_degree-1) - k * (sub_degree-1) - j) / ( (sub_degree-1) * (sub_degree-1) );

              embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell +
                                  GeometryInfo<dim>::lines_per_cell * (this->degree-1) +
                                  face * (this->degree-1) * (this->degree-1) +
                                  k * (this->degree-1) + j;
            }
          // hex
          else if (i<(GeometryInfo<dim>::vertices_per_cell +
                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) +
                   GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) +
                   GeometryInfo<dim>::hexes_per_cell * (sub_degree-1) * (sub_degree-1) * (sub_degree-1))
            {
              const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell -
                                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
                                      GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) ) % (sub_degree-1);
              const unsigned int k = ( (i - GeometryInfo<dim>::vertices_per_cell -
                                        GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
                                        GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) - j) / (sub_degree-1) ) % (sub_degree-1);
              const unsigned int l = (i - GeometryInfo<dim>::vertices_per_cell -
                                      GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
                                      GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) - j - k * (sub_degree-1)) / ( (sub_degree-1) * (sub_degree-1) );

              embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell +
                                  GeometryInfo<dim>::lines_per_cell * (this->degree-1) +
                                  GeometryInfo<dim>::quads_per_cell * (this->degree-1) * (this->degree-1) +
                                  l * (this->degree-1) * (this->degree-1) + k * (this->degree-1) + j;
            }
        }

      return embedding_dofs;
    }
  else
    {
      Assert(false, ExcNotImplemented ());
      return std::vector<unsigned int> ();
    }
}



template <int dim>
std::pair<Table<2,bool>, std::vector<unsigned int> >
FE_Q_Hierarchical<dim>::get_constant_modes () const
{
  Table<2,bool> constant_modes(1, this->dofs_per_cell);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    constant_modes(0,i) = true;
  for (unsigned int i=GeometryInfo<dim>::vertices_per_cell; i<this->dofs_per_cell; ++i)
    constant_modes(0,i) = false;
  return std::pair<Table<2,bool>, std::vector<unsigned int> >
         (constant_modes, std::vector<unsigned int>(1, 0));
}



template <int dim>
std::size_t
FE_Q_Hierarchical<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



// explicit instantiations
#include "fe_q_hierarchical.inst"


DEAL_II_NAMESPACE_CLOSE
