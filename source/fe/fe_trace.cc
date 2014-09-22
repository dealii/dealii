// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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
#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_q.h>
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
    TensorProductPolynomials<dim-1>(Polynomials::LagrangeEquidistant::generate_complete_basis(degree)),
    FiniteElementData<dim>(get_dpo_vector(degree), 1, degree, FiniteElementData<dim>::L2),
    std::vector<bool>(1,true))
{
  Assert (degree > 0,
          ExcMessage ("FE_Trace can only be used for polynomial degrees "
                      "greater than zero"));
  std::vector<unsigned int> renumber (this->dofs_per_face);
  FETools::hierarchic_to_lexicographic_numbering<dim-1> (degree, renumber);
  this->poly_space.set_numbering(renumber);

  // Initialize face support points
  FE_Q<dim,spacedim> fe_q(degree);
  this->unit_face_support_points = fe_q.get_unit_face_support_points();
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
  // note that the
  // FETools::get_fe_from_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

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

  // let's see whether this is a
  // vertex
  if (shape_index < this->first_line_index)
    {
      // for Q elements, there is
      // one dof per vertex, so
      // shape_index==vertex_number. check
      // whether this vertex is
      // on the given face. thus,
      // for each face, give a
      // list of vertices
      const unsigned int vertex_no = shape_index;
      Assert (vertex_no < GeometryInfo<dim>::vertices_per_cell,
              ExcInternalError());

      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
        if (GeometryInfo<dim>::face_to_cell_vertices(face_index, v) == vertex_no)
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

      // in 2d, the line is the
      // face, so get the line
      // index
      if (dim == 2)
        return (line_index == face_index);
      else if (dim == 3)
        {
          // see whether the
          // given line is on the
          // given face.
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
            if (GeometryInfo<3>::face_to_cell_lines(face_index, l) == line_index)
              return true;

          return false;
        }
      else
        Assert (false, ExcNotImplemented());
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
      // can only happen in 3d,
      // but this case has
      // already been covered
      // above
      Assert (false, ExcNotImplemented());
      return false;
    }

  // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;
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
  AssertThrow(deg>0,ExcMessage("FE_TraceQ needs to be of degree > 0."));
  std::vector<unsigned int> dpo(dim+1, 1U);
  dpo[dim]=0;
  dpo[0]=1;
  for (unsigned int i=1; i<dim; ++i)
    dpo[i] = dpo[i-1]*(deg-1);
  return dpo;
}

// explicit instantiations
#include "fe_trace.inst"


DEAL_II_NAMESPACE_CLOSE
