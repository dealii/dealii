//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_accessor.templates.h>
#include <grid/tria_iterator.templates.h>
#include <base/geometry_info.h>
#include <grid/grid_tools.h>
#include <fe/mapping_q1.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


template <int structdim,int dim>
const unsigned int TriaAccessor<structdim,dim>::objectdim;


// anonymous namespace with two little helper functions
namespace{
				   // given the number of face's child
				   // (subface_no), return the number of the
				   // subface concerning the FaceRefineCase of
				   // the face
  inline
  unsigned int translate_subface_no(const TriaIterator<3,TriaObjectAccessor<2, 3> > &face,
				    const unsigned int                               subface_no)
  {
    Assert(face->has_children(), ExcInternalError());
    Assert(subface_no<face->n_children(), ExcInternalError());

    if(face->child(subface_no)->has_children())
				       // although the subface is refine, it
				       // still matches the face of the cell
				       // invoking the
				       // neighbor_of_coarser_neighbor
				       // function. this means that we are
				       // looking from one cell (anisotropic
				       // child) to a coarser neighbor which is
				       // refined stronger than we are
				       // (isotropically). So we won't be able
				       // to use the neighbor_child_on_subface
				       // function anyway, as the neighbor is
				       // not active. In this case, simply
				       // return the subface_no.
      return subface_no;
    
    const bool first_child_has_children=face->child(0)->has_children();
				     // if the first child has children
				     // (FaceRefineCase case_x1y or case_y1x),
				     // then the current subface_no needs to be
				     // 1 and the result of this function is 2,
				     // else simply return the given number,
				     // which is 0 or 1 in an anisotropic case
				     // (case_x, case_y, casex2y or casey2x) or
				     // 0...3 in an isotropic case (case_xy)
    return subface_no + first_child_has_children;
  }


  
				   // given the number of face's child
				   // (subface_no) and grandchild
				   // (subsubface_no), return the number of the
				   // subface concerning the FaceRefineCase of
				   // the face
  inline
  unsigned int translate_subface_no(const TriaIterator<3,TriaObjectAccessor<2, 3> > &face,
				    const unsigned int                               subface_no,
				    const unsigned int                               subsubface_no)
  {
    Assert(face->has_children(), ExcInternalError());
				     // the subface must be refined, otherwise
				     // we would have ended up in the second
				     // function of this name...
    Assert(face->child(subface_no)->has_children(), ExcInternalError());
    Assert(subsubface_no<face->child(subface_no)->n_children(), ExcInternalError());
				     // This can only be an anisotropic refinement case
    Assert(face->refinement_case() < RefinementCase<2>::isotropic_refinement,
	   ExcInternalError());
    
    const bool first_child_has_children=face->child(0)->has_children();

    const unsigned int e = deal_II_numbers::invalid_unsigned_int;
    
				     // array containing the translation of the
				     // numbers,
				     //
				     // first index: subface_no
				     // second index: subsubface_no
				     // third index: does the first subface have children? -> no and yes
    unsigned int translated_subface_no[2][2][2]
      =
      {{{e,0},       // first  subface, first  subsubface, first_child_has_children==no and yes
	{e,1}},      // first  subface, second subsubface, first_child_has_children==no and yes
       {{1,2},       // second subface, first  subsubface, first_child_has_children==no and yes
	{2,3}}};     // second subface, second subsubface, first_child_has_children==no and yes
    
    Assert(translated_subface_no[subface_no][subsubface_no][first_child_has_children]!=e,
	   ExcInternalError());

    return translated_subface_no[subface_no][subsubface_no][first_child_has_children];
  }
}





/*------------------------ Functions: LineAccessor ---------------------------*/

template <int dim>
const unsigned int TriaObjectAccessor<1, dim>::objectdim;


template <int dim>
void
TriaObjectAccessor<1, dim>::set (const internal::Triangulation::TriaObject<1>  &line) const
{
  objects().cells[this->present_index] = line;
}



template <int dim>
unsigned int TriaObjectAccessor<1, dim>::vertex_index (const unsigned int i) const
{
  Assert (i<2, ExcIndexRange(i,0,2));
  return objects().cells[this->present_index].face (i);
}



template <int dim>
Point<dim> &
TriaObjectAccessor<1, dim>::vertex (const unsigned int i) const
{
  return const_cast<Point<dim> &> (this->tria->vertices[vertex_index(i)]);
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  objects().used[this->present_index] = true;
}



template <int dim>
void TriaObjectAccessor<1, dim>::clear_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  objects().used[this->present_index] = false;
}



template <int dim>
void TriaObjectAccessor<1, dim>::recursively_set_user_flag () const
{
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_set_user_flag ();
}



template <int dim>
void TriaObjectAccessor<1, dim>::recursively_clear_user_flag () const
{
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_clear_user_flag ();
}



template <int dim>
void TriaObjectAccessor<1, dim>::clear_user_data () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().clear_user_data(this->present_index);
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_user_pointer (void *p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_pointer(this->present_index) = p;
}



template <int dim>
void TriaObjectAccessor<1, dim>::clear_user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_pointer(this->present_index) = 0;
}



template <int dim>
void * TriaObjectAccessor<1, dim>::user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_pointer(this->present_index);
}



template <int dim>
void
TriaObjectAccessor<1, dim>::recursively_set_user_pointer (void *p) const
{
  set_user_pointer (p);

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_set_user_pointer (p);
}



template <int dim>
void
TriaObjectAccessor<1, dim>::recursively_clear_user_pointer () const
{
  clear_user_pointer ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_clear_user_pointer ();
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_user_index (unsigned int p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_index(this->present_index) = p;
}



template <int dim>
void TriaObjectAccessor<1, dim>::clear_user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_index(this->present_index) = 0;
}



template <int dim>
unsigned int TriaObjectAccessor<1, dim>::user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_index(this->present_index);
}



template <int dim>
void
TriaObjectAccessor<1, dim>::recursively_set_user_index (unsigned int p) const
{
  set_user_index (p);

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_set_user_index (p);
}



template <int dim>
void
TriaObjectAccessor<1, dim>::recursively_clear_user_index () const
{
  clear_user_index ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_clear_user_index ();
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_children (const unsigned int i, const int index) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i%2==0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));
  Assert ((index==-1) ||
	  (!has_children() && (index>=0)),
	  TriaAccessorExceptions::ExcCantSetChildren(index));
  
  objects().children[this->present_index] = index;
}



template <int dim>
void TriaObjectAccessor<1, dim>::clear_children () const
{
  set_children (0,-1);
}



template <int dim>
unsigned char TriaObjectAccessor<1, dim>::boundary_indicator () const
{
  Assert (dim>=2, ExcImpossibleInDim(dim));
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return objects().material_id[this->present_index];
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_boundary_indicator (const unsigned char boundary_ind) const
{
  Assert (dim>=2, ExcImpossibleInDim(dim));
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  objects().material_id[this->present_index] = boundary_ind;
}



template <int dim>
void TriaObjectAccessor<1, dim>::set_all_boundary_indicators (const unsigned char boundary_ind) const
{
  set_boundary_indicator (boundary_ind);
}



template <int dim>
bool TriaObjectAccessor<1, dim>::at_boundary () const
{
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
}



template <int dim>
double TriaObjectAccessor<1, dim>::diameter () const
{
  return std::sqrt((vertex(1)-vertex(0)).square());
}



template <int dim>
Point<dim> TriaObjectAccessor<1, dim>::center () const
{
  return (vertex(1)+vertex(0))/2.;
}



template <int dim>
Point<dim> TriaObjectAccessor<1, dim>::barycenter () const
{
  return (vertex(1)+vertex(0))/2.;
}



template <int dim>
double TriaObjectAccessor<1, dim>::measure () const
{
  return std::sqrt((vertex(1)-vertex(0)).square());
}



template <int dim>
unsigned int TriaObjectAccessor<1, dim>::number_of_children () const
{
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<2; ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
}



/*------------------------ Functions: QuadAccessor ---------------------------*/

template <int dim>
const unsigned int TriaObjectAccessor<2, dim>::objectdim;


template <int dim>
void
TriaObjectAccessor<2, dim>::set (const internal::Triangulation::TriaObject<2> &quad) const
{
  objects().cells[this->present_index] = quad;
}



template <int dim>
unsigned int TriaObjectAccessor<2, dim>::vertex_index (const unsigned int corner) const
{
  Assert (corner<4, ExcIndexRange(corner,0,4));

				   // table used to switch the vertices, if the
				   // line orientation is wrong,
				   //
				   // first index: line orientation 0: false or
				   // 1: true=standard
				   //
				   // second index: vertex index to be switched
				   // (or not)
  
  static const unsigned int switch_table[2][2]={{1,0},{0,1}};
  
  return line(corner%2)->vertex_index(switch_table[line_orientation(corner%2)][corner/2]);
}



template <int dim>
Point<dim> &
TriaObjectAccessor<2, dim>::vertex (const unsigned int i) const
{
  return const_cast<Point<dim> &> (this->tria->vertices[vertex_index(i)]);
}



template <int dim>
void
TriaObjectAccessor<2, dim>::set_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  objects().used[this->present_index] = true;
}



template <int dim>
void TriaObjectAccessor<2, dim>::clear_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  objects().used[this->present_index] = false;
}



template <int dim>
void TriaObjectAccessor<2, dim>::recursively_set_user_flag () const
{
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_flag ();
}



template <int dim>
void TriaObjectAccessor<2, dim>::recursively_clear_user_flag () const
{
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_flag ();
}



template <int dim>
void TriaObjectAccessor<2, dim>::clear_user_data () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().clear_user_data(this->present_index);
}



template <int dim>
void TriaObjectAccessor<2, dim>::set_user_pointer (void *p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_pointer(this->present_index) = p;
}



template <int dim>
void TriaObjectAccessor<2, dim>::clear_user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_pointer(this->present_index) = 0;
}



template <int dim>
void * TriaObjectAccessor<2, dim>::user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_pointer(this->present_index);
}



template <int dim>
void
TriaObjectAccessor<2, dim>::recursively_set_user_pointer (void *p) const
{
  set_user_pointer (p);

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_pointer (p);
}



template <int dim>
void
TriaObjectAccessor<2, dim>::recursively_clear_user_pointer () const
{
  clear_user_pointer ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_pointer ();
}



template <int dim>
void TriaObjectAccessor<2, dim>::set_user_index (unsigned int p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_index(this->present_index) = p;
}



template <int dim>
void TriaObjectAccessor<2, dim>::clear_user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  objects().user_index(this->present_index) = 0;
}



template <int dim>
unsigned int  TriaObjectAccessor<2, dim>::user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return objects().user_index(this->present_index);
}



template <int dim>
void
TriaObjectAccessor<2, dim>::recursively_set_user_index (unsigned int p) const
{
  set_user_index (p);

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_index (p);
}



template <int dim>
void
TriaObjectAccessor<2, dim>::recursively_clear_user_index () const
{
  clear_user_index ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_index ();
}



template <int dim>
void TriaObjectAccessor<2, dim>::set_children (const unsigned int i, const int index) const
{
  Assert (used(),
	  TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i%2==0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));
  Assert ((index==-1) ||
	  (i==0 && !has_children() && (index>=0)) ||
	  (i>0  &&  has_children() && (index>=0) &&
	   objects().children[2*this->present_index+i/2] == -1),
	  TriaAccessorExceptions::ExcCantSetChildren(index));

  objects().children[2*this->present_index+i/2] = index;
}



template <int dim>
void TriaObjectAccessor<2, dim>::clear_children () const
{
  set_children (0, -1);
  set_children (2, -1);
}



template <int dim>
unsigned char TriaObjectAccessor<2, dim>::boundary_indicator () const
{
  Assert (dim>2, ExcImpossibleInDim(dim));
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  return objects().material_id[this->present_index];
}



template <int dim>
void TriaObjectAccessor<2, dim>::set_boundary_indicator (const unsigned char boundary_ind) const
{
  Assert (dim>2, ExcImpossibleInDim(dim));
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());

  objects().material_id[this->present_index] = boundary_ind;
}



template <int dim>
void TriaObjectAccessor<2, dim>::set_all_boundary_indicators (const unsigned char boundary_ind) const
{
  set_boundary_indicator (boundary_ind);

  for (unsigned int i=0; i<4; ++i)
    this->line(i)->set_boundary_indicator (boundary_ind);
}



template <int dim>
bool TriaObjectAccessor<2, dim>::at_boundary () const
{
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
}



template <int dim>
double TriaObjectAccessor<2, dim>::diameter () const
{
  return std::sqrt(std::max((vertex(3)-vertex(0)).square(),
			    (vertex(2)-vertex(1)).square()));
}



template <int dim>
Point<dim> TriaObjectAccessor<2, dim>::center () const
{
  return (vertex(0)+vertex(1)+vertex(2)+vertex(3))/4.;
}


#if deal_II_dimension == 2

template <>
Point<2> TriaObjectAccessor<2, 2>::barycenter () const
{
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independently, so we write this function
				   // for 2D and 3D separately
/*
  Get the computation of the barycenter by this little Maple script. We
  use the bilinear mapping of the unit quad to the real quad. However,
  every transformation mapping the unit faces to strait lines should
  do.

  Remember that the area of the quad is given by
  |K| = \int_K 1 dx dy  = \int_{\hat K} |det J| d(xi) d(eta)
  and that the barycenter is given by
  \vec x_s = 1/|K| \int_K \vec x dx dy
  = 1/|K| \int_{\hat K} \vec x(xi,eta) |det J| d(xi) d(eta)

  # x and y are arrays holding the x- and y-values of the four vertices
  # of this cell in real space. 
  x := array(0..3);
  y := array(0..3);
  tphi[0] := (1-xi)*(1-eta):
  tphi[1] :=     xi*(1-eta):
  tphi[2] := (1-xi)*eta:
  tphi[3] :=     xi*eta:
  x_real := sum(x[s]*tphi[s], s=0..3):
  y_real := sum(y[s]*tphi[s], s=0..3):
  detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi):

  measure := simplify ( int ( int (detJ, xi=0..1), eta=0..1)):

  xs := simplify (1/measure * int ( int (x_real * detJ, xi=0..1), eta=0..1)):
  ys := simplify (1/measure * int ( int (y_real * detJ, xi=0..1), eta=0..1)):
  readlib(C):

  C(array(1..2, [xs, ys]), optimized);
*/

  const double x[4] = { vertex(0)(0),
			vertex(1)(0),
			vertex(2)(0),
			vertex(3)(0)  };
  const double y[4] = { vertex(0)(1),
			vertex(1)(1),
			vertex(2)(1),
			vertex(3)(1)  };
  const double t1 = x[0]*x[1];
  const double t3 = x[0]*x[0];
  const double t5 = x[1]*x[1];
  const double t9 = y[0]*x[0];
  const double t11 = y[1]*x[1];
  const double t14 = x[2]*x[2];
  const double t16 = x[3]*x[3];
  const double t20 = x[2]*x[3];
  const double t27 = t1*y[1]+t3*y[1]-t5*y[0]-t3*y[2]+t5*y[3]+t9*x[2]-t11*x[3]-t1*y[0]-t14*y[3]+t16*y[2]-t16*y[1]+t14*y[0]-t20*y[3]-x[0]*x[2]*y[2]+x[1]*x[3]*y[3]+t20*y[2];
  const double t37 = 1/(-x[1]*y[0]+x[1]*y[3]+y[0]*x[2]+x[0]*y[1]-x[0]*y[2]-y[1]*x[3]-x[2]*y[3]+x[3]*y[2]);
  const double t39 = y[2]*y[2];
  const double t51 = y[0]*y[0];
  const double t53 = y[1]*y[1];
  const double t59 = y[3]*y[3];
  const double t63 = t39*x[3]+y[2]*y[0]*x[2]+y[3]*x[3]*y[2]-y[2]*x[2]*y[3]-y[3]*y[1]*x[3]-t9*y[2]+t11*y[3]+t51*x[2]-t53*x[3]-x[1]*t51+t9*y[1]-t11*y[0]+x[0]*t53-t59*x[2]+t59*x[1]-t39*x[0];

  return Point<2> (t27*t37/3, t63*t37/3);
}



template <>
double TriaObjectAccessor<2, 2>::measure () const
{
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independently, so we write this function
				   // for 2D and 3D separately
/*
  Get the computation of the measure by this little Maple script. We
  use the blinear mapping of the unit quad to the real quad. However,
  every transformation mapping the unit faces to straight lines should
  do.

  Remember that the area of the quad is given by
  \int_K 1 dx dy  = \int_{\hat K} |det J| d(xi) d(eta)

  # x and y are arrays holding the x- and y-values of the four vertices
  # of this cell in real space. 
  x := array(0..3);
  y := array(0..3);
  tphi[0] := (1-xi)*(1-eta):
  tphi[1] :=     xi*(1-eta):
  tphi[2] := (1-xi)*eta:
  tphi[3] :=     xi*eta:
  x_real := sum(x[s]*tphi[s], s=0..3):
  y_real := sum(y[s]*tphi[s], s=0..3):
  detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi):

  measure := simplify ( int ( int (detJ, xi=0..1), eta=0..1)):
  readlib(C):

  C(measure, optimized);

  additional optimizaton: divide by 2 only one time
*/

  const double x[4] = { vertex(0)(0),
			vertex(1)(0),
			vertex(2)(0),
			vertex(3)(0)  };
  const double y[4] = { vertex(0)(1),
			vertex(1)(1),
			vertex(2)(1),
			vertex(3)(1)  };

  return (-x[1]*y[0]+x[1]*y[3]+y[0]*x[2]+x[0]*y[1]-x[0]*y[2]-y[1]*x[3]-x[2]*y[3]+x[3]*y[2])/2;
}


#endif


#if deal_II_dimension == 3

template <>
Point<3> TriaObjectAccessor<2, 3>::barycenter () const
{
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independently, so we write this function
				   // for 2D and 3D separately
/*
  To compute the barycenter, we first have to find out the size of
  an area element in real space; this equals the Jacobian determinant
  at this point, then. To do so, find out the points in real space
  belonging to
  xi,eta
  xi+dxi, eta
  xi, eta+deta
  To this end, remember that the mapping is
  x(xi,eta) = \sum_i=0^3 x_i phi_i(xi,eta)
  y(xi,eta) = \sum_i=0^3 y_i phi_i(xi,eta)
  z(xi,eta) = \sum_i=0^3 z_i phi_i(xi,eta)
  with x_i, y_i being the four vertices and the phi_i the shape functions
  corresponding to these four vertices of this face. Now the real space
  points belonging to the above points on the unit face are:
  x, y, z
  x+sum x_i \partial_xi phi_i(xi,eta) dxi,
  y+sum y_i \partial_xi phi_i(xi,eta) dxi,
  z+sum z_i \partial_xi phi_i(xi,eta) dxi
  x+sum x_i \partial_eta phi_i(xi,eta) deta,
  y+sum y_i \partial_eta phi_i(xi,eta) deta,
  z+sum z_i \partial_eta phi_i(xi,eta) deta
  The unit infinitesimal vectors at the point xi,eta have the direction
  dxi, 0
  0, deta
  and are therefore in real space
  sum x_i \partial_xi phi_i(xi,eta) dxi,
  sum y_i \partial_xi phi_i(xi,eta) dxi,
  sum z_i \partial_xi phi_i(xi,eta) dxi
  sum x_i \partial_eta phi_i(xi,eta) deta,
  sum y_i \partial_eta phi_i(xi,eta) deta
  sum z_i \partial_eta phi_i(xi,eta) deta
  or in other form:
  \partial_xi (x,y,z) dxi
  \partial_eta (x,y,z) deta
  Then the area element is the length of the cross-product of these two vectors and
  the Jacobian determinant is this expression divided by dxi deta:
  |J| = |(\partial_xi (x,y,z) \times (\partial_eta (x,y,z)|

  There is a script in the deal.II/source/fe/scripts/3d directory, which does
  these computations in Maple.
*/
  Assert (false, ExcNotImplemented());
  
  return Point<3> ();
}



template <>
double TriaObjectAccessor<2, 3>::measure () const
{
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independently, so we write this function
				   // for 2D and 3D separately
				   //
				   // for documentation, see the barycenter
				   // function above.
  Assert (false, ExcNotImplemented());
  return 0;
}

#endif


template <int dim>
unsigned int TriaObjectAccessor<2, dim>::number_of_children () const
{
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<n_children(); ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
}


template <int dim>
void
TriaObjectAccessor<2, dim>::
set_line_orientation (const unsigned int ,
                      const bool         ) const
{
  Assert(false,
	 ExcMessage("The line orientation of Quads can only be set in 3D. In 2D the lines are alway in standard orientation."));}



#if deal_II_dimension == 3

template <>
void
TriaObjectAccessor<2, 3>::
set_line_orientation (const unsigned int line,
                      const bool         orientation) const
{
  const int dim=3;
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (line<GeometryInfo<dim>::lines_per_face,
          ExcIndexRange (line, 0, GeometryInfo<3>::lines_per_face));
  Assert (this->present_index * GeometryInfo<dim>::lines_per_face + line
          < this->tria->faces->quads.line_orientations.size(),
          ExcInternalError());
          
  this->tria->faces
    ->quads.line_orientations[this->present_index *
                              GeometryInfo<dim>::lines_per_face
                              +
                              line]  = orientation;
}

#endif


/*------------------------ Functions: TriaObjectAccessor ---------------------------*/


template <int dim>
const unsigned int TriaObjectAccessor<3, dim>::objectdim;


#if deal_II_dimension == 3

template <>
void
TriaObjectAccessor<3, 3>::set (const internal::Triangulation::TriaObject<3> &hex) const
{
  this->tria->levels[this->present_level]
    ->cells.cells[this->present_index] = hex;
}

#endif


template <int dim>
unsigned int TriaObjectAccessor<3, dim>::vertex_index (const unsigned int corner) const
{
  Assert (corner<8, ExcIndexRange(corner,0,8));

				   // get the corner indices by asking either
				   // the bottom or the top face for its
				   // vertices. handle non-standard faces by
				   // calling the vertex reordering function
				   // from GeometryInfo

  				   // bottom face (4) for first four vertices,
				   // top face (5) for the rest
  const unsigned int face_index=4+corner/4;

  return quad(face_index)->vertex_index(GeometryInfo<dim>
					::standard_to_real_face_vertex(corner%4,
								       face_orientation(face_index),
								       face_flip(face_index),
								       face_rotation(face_index)));
}



template <int dim>
Point<dim> &
TriaObjectAccessor<3, dim>::vertex (const unsigned int i) const
{
  return const_cast<Point<dim> &> (this->tria->vertices[vertex_index(i)]);
}


#if deal_II_dimension == 3

template <>
void TriaObjectAccessor<3, 3>::set_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  this->tria->levels[this->present_level]->cells.used[this->present_index] = true;
}


template <>
void TriaObjectAccessor<3, 3>::clear_used_flag () const
{
  Assert (this->state() == IteratorState::valid,
	  TriaAccessorExceptions::ExcDereferenceInvalidObject());
  this->tria->levels[this->present_level]->cells.used[this->present_index] = false;
}

#endif


template <int dim>
void TriaObjectAccessor<3, dim>::recursively_set_user_flag () const
{
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_flag ();
}



template <int dim>
void TriaObjectAccessor<3, dim>::recursively_clear_user_flag () const
{
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_flag ();
}


#if deal_II_dimension == 3

template <>
void TriaObjectAccessor<3, 3>::clear_user_data () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.clear_user_data(this->present_index);
}



template <>
void TriaObjectAccessor<3, 3>::set_user_pointer (void *p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_pointer(this->present_index) = p;
}


template <>
void TriaObjectAccessor<3, 3>::clear_user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_pointer(this->present_index) = 0;
}


template <>
void * TriaObjectAccessor<3, 3>::user_pointer () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->cells.user_pointer(this->present_index);
}


template <>
void TriaObjectAccessor<3, 3>::set_user_index (unsigned int p) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_index(this->present_index) = p;
}


template <>
void TriaObjectAccessor<3, 3>::clear_user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.user_index(this->present_index) = 0;
}


template <>
unsigned int  TriaObjectAccessor<3, 3>::user_index () const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->cells.user_index(this->present_index);
}

#endif


template <int dim>
void
TriaObjectAccessor<3, dim>::recursively_set_user_pointer (void *p) const
{
  set_user_pointer (p);

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_pointer (p);
}



template <int dim>
void
TriaObjectAccessor<3, dim>::recursively_clear_user_pointer () const
{
  clear_user_pointer ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_pointer ();
}


template <int dim>
void
TriaObjectAccessor<3, dim>::recursively_set_user_index (unsigned int p) const
{
  set_user_index (p);

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_set_user_index (p);
}



template <int dim>
void
TriaObjectAccessor<3, dim>::recursively_clear_user_index () const
{
  clear_user_index ();

  if (has_children())
    for (unsigned int c=0; c<n_children(); ++c)
      child(c)->recursively_clear_user_index ();
}


#if deal_II_dimension == 3

template <>
void TriaObjectAccessor<3,3>::set_children (const unsigned int i, const int index) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i%2==0, TriaAccessorExceptions::ExcSetOnlyEvenChildren(i));
  Assert ((index==-1) ||
	  (i==0 && !has_children() && (index>=0)) ||
	  (i>0  &&  has_children() && (index>=0) &&
	   this->tria->levels[this->present_level]->
	   cells.children[4*this->present_index+i/2] == -1),
	  TriaAccessorExceptions::ExcCantSetChildren(index));

  this->tria->levels[this->present_level]->cells.children[4*this->present_index+i/2] = index;
}

#endif


template <int dim>
void TriaObjectAccessor<3, dim>::clear_children () const
{
  set_children (0,-1);
  set_children (2,-1);
  set_children (4,-1);
  set_children (6,-1);
}



template <int dim>
unsigned char TriaObjectAccessor<3, dim>::boundary_indicator () const
{
  Assert (false, ExcImpossibleInDim(dim));
  return 255;
}



template <int dim>
void TriaObjectAccessor<3, dim>::set_boundary_indicator (const unsigned char) const
{
  Assert (false, ExcImpossibleInDim(dim));
}



template <int dim>
void TriaObjectAccessor<3, dim>::set_all_boundary_indicators (const unsigned char) const
{
  Assert (false, ExcImpossibleInDim(dim));
}



template <int dim>
bool TriaObjectAccessor<3, dim>::at_boundary () const
{
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
}



template <int dim>
double TriaObjectAccessor<3, dim>::diameter () const
{
  return std::sqrt(std::max( std::max((vertex(7)-vertex(0)).square(),
				      (vertex(6)-vertex(1)).square()),
			     std::max((vertex(2)-vertex(5)).square(),
				      (vertex(3)-vertex(4)).square()) ));
}



template <int dim>
Point<dim> TriaObjectAccessor<3, dim>::center () const
{
  return (vertex(0)+vertex(1)+vertex(2)+vertex(3)+
	  vertex(4)+vertex(5)+vertex(6)+vertex(7))/8.;
}


#if deal_II_dimension == 3

template <>
Point<3> TriaObjectAccessor<3, 3>::barycenter () const
{
/*
  Get the computation of the barycenter by this little Maple script. We
  use the trilinear mapping of the unit hex to the real hex.

  Remember that the area of the hex is given by
  |K| = \int_K 1 dx dy dz = \int_{\hat K} |det J| d(xi) d(eta) d(zeta)
  and that the barycenter is given by
  \vec x_s = 1/|K| \int_K \vec x dx dy dz
  = 1/|K| \int_{\hat K} \vec x(xi,eta,zeta) |det J| d(xi) d(eta) d(zeta)

  Note, that in the ordering of the shape functions tphi[0]-tphi[7]
  below, eta and zeta have been exchanged (zeta belongs to the y, and
  eta to the z direction). However, the resulting Jacobian determinant
  detJ should be the same, as a matrix and the matrix created from it
  by exchanging two consecutive lines and two neighboring columns have
  the same determinant.
  
  # x, y and z are arrays holding the x-, y- and z-values of the four vertices
  # of this cell in real space. 
  x := array(0..7):
  y := array(0..7):
  z := array(0..7):
  tphi[0] := (1-xi)*(1-eta)*(1-zeta):
  tphi[1] := xi*(1-eta)*(1-zeta):
  tphi[2] := xi*eta*(1-zeta):
  tphi[3] := (1-xi)*eta*(1-zeta):
  tphi[4] := (1-xi)*(1-eta)*zeta:
  tphi[5] := xi*(1-eta)*zeta:
  tphi[6] := xi*eta*zeta:
  tphi[7] := (1-xi)*eta*zeta:
  x_real := sum(x[s]*tphi[s], s=0..7):
  y_real := sum(y[s]*tphi[s], s=0..7):
  z_real := sum(z[s]*tphi[s], s=0..7):
  with (linalg):
  J := matrix(3,3, [[diff(x_real, xi), diff(x_real, eta), diff(x_real, zeta)],
  [diff(y_real, xi), diff(y_real, eta), diff(y_real, zeta)],
  [diff(z_real, xi), diff(z_real, eta), diff(z_real, zeta)]]): 
  detJ := det (J):

  measure := simplify ( int ( int ( int (detJ, xi=0..1), eta=0..1), zeta=0..1)):

  xs := simplify (1/measure * int ( int ( int (x_real * detJ, xi=0..1), eta=0..1), zeta=0..1)):
  ys := simplify (1/measure * int ( int ( int (y_real * detJ, xi=0..1), eta=0..1), zeta=0..1)):
  zs := simplify (1/measure * int ( int ( int (z_real * detJ, xi=0..1), eta=0..1), zeta=0..1)):

  readlib(C):

  C(array(1..3, [xs, ys, zs]));


  This script takes more than several hours when using an old version
  of maple on an old and slow computer. Therefore, when changing to
  the new deal.II numbering scheme (lexicographic numbering) the code
  lines below have not been reproduced with maple but only the
  ordering of points in the definitions of x[], y[] and z[] have been
  changed.

  For the case, someone is willing to rerun the maple script, he/she
  should use following ordering of shape functions:
  
  tphi[0] := (1-xi)*(1-eta)*(1-zeta):
  tphi[1] :=     xi*(1-eta)*(1-zeta):
  tphi[2] := (1-xi)*    eta*(1-zeta):
  tphi[3] :=     xi*    eta*(1-zeta):
  tphi[4] := (1-xi)*(1-eta)*zeta:
  tphi[5] :=     xi*(1-eta)*zeta:
  tphi[6] := (1-xi)*    eta*zeta:
  tphi[7] :=     xi*    eta*zeta:

  and change the ordering of points in the definitions of x[], y[] and
  z[] back to the standard ordering.
*/

  const double x[8] = { vertex(0)(0),
			vertex(1)(0),
			vertex(5)(0),
			vertex(4)(0),
			vertex(2)(0),
			vertex(3)(0),
			vertex(7)(0),
			vertex(6)(0)   };
  const double y[8] = { vertex(0)(1),
			vertex(1)(1),
			vertex(5)(1),
			vertex(4)(1),
			vertex(2)(1),
			vertex(3)(1),
			vertex(7)(1),
			vertex(6)(1)   };
  const double z[8] = { vertex(0)(2),
			vertex(1)(2),
			vertex(5)(2),
			vertex(4)(2),
			vertex(2)(2),
			vertex(3)(2),
			vertex(7)(2),
			vertex(6)(2)   };

  double s1, s2, s3, s4, s5, s6, s7, s8;
  
  s1 = 1.0/6.0;
  s8 = -x[2]*x[2]*y[0]*z[3]-2.0*z[6]*x[7]*x[7]*y[4]-z[5]*x[7]*x[7]*y[4]-z
       [6]*x[7]*x[7]*y[5]+2.0*y[6]*x[7]*x[7]*z[4]-z[5]*x[6]*x[6]*y[4]+x[6]*x[6]*y[4]*z
       [7]-z[1]*x[0]*x[0]*y[2]-x[6]*x[6]*y[7]*z[4]+2.0*x[6]*x[6]*y[5]*z[7]-2.0*x[6]*x
       [6]*y[7]*z[5]+y[5]*x[6]*x[6]*z[4]+2.0*x[5]*x[5]*y[4]*z[6]+x[0]*x[0]*y[7]*z[4]
       -2.0*x[5]*x[5]*y[6]*z[4];
  s7 = s8-y[6]*x[5]*x[5]*z[7]+z[6]*x[5]*x[5]*y[7]-y[1]*x[0]*x[0]*z[5]+x[7]*
       z[5]*x[4]*y[7]-x[7]*y[6]*x[5]*z[7]-2.0*x[7]*x[6]*y[7]*z[4]+2.0*x[7]*x[6]*y[4]*z
       [7]-x[7]*x[5]*y[7]*z[4]-2.0*x[7]*y[6]*x[4]*z[7]-x[7]*y[5]*x[4]*z[7]+x[2]*x[2]*y
       [3]*z[0]-x[7]*x[6]*y[7]*z[5]+x[7]*x[6]*y[5]*z[7]+2.0*x[1]*x[1]*y[0]*z[5]+x[7]*z
       [6]*x[5]*y[7];
  s8 = -2.0*x[1]*x[1]*y[5]*z[0]+z[1]*x[0]*x[0]*y[5]+2.0*x[2]*x[2]*y[3]*z[1]
       -z[5]*x[4]*x[4]*y[1]+y[5]*x[4]*x[4]*z[1]-2.0*x[5]*x[5]*y[4]*z[1]+2.0*x[5]*x[5]*
       y[1]*z[4]-2.0*x[2]*x[2]*y[1]*z[3]-y[1]*x[2]*x[2]*z[0]+x[7]*y[2]*x[3]*z[7]+x[7]*
       z[2]*x[6]*y[3]+2.0*x[7]*z[6]*x[4]*y[7]+z[5]*x[1]*x[1]*y[4]+z[1]*x[2]*x[2]*y[0]
       -2.0*y[0]*x[3]*x[3]*z[7];
  s6 = s8+2.0*z[0]*x[3]*x[3]*y[7]-x[7]*x[2]*y[3]*z[7]-x[7]*z[2]*x[3]*y[7]+x
       [7]*x[2]*y[7]*z[3]-x[7]*y[2]*x[6]*z[3]+x[4]*x[5]*y[1]*z[4]-x[4]*x[5]*y[4]*z[1]+
       x[4]*z[5]*x[1]*y[4]-x[4]*y[5]*x[1]*z[4]-2.0*x[5]*z[5]*x[4]*y[1]-2.0*x[5]*y[5]*x
       [1]*z[4]+2.0*x[5]*z[5]*x[1]*y[4]+2.0*x[5]*y[5]*x[4]*z[1]-x[6]*z[5]*x[7]*y[4]-z
       [2]*x[3]*x[3]*y[6]+s7;
  s8 = -2.0*x[6]*z[6]*x[7]*y[5]-x[6]*y[6]*x[4]*z[7]+y[2]*x[3]*x[3]*z[6]+x
       [6]*y[6]*x[7]*z[4]+2.0*y[2]*x[3]*x[3]*z[7]+x[0]*x[1]*y[0]*z[5]+x[0]*y[1]*x[5]*z
       [0]-x[0]*z[1]*x[5]*y[0]-2.0*z[2]*x[3]*x[3]*y[7]+2.0*x[6]*z[6]*x[5]*y[7]-x[0]*x
       [1]*y[5]*z[0]-x[6]*y[5]*x[4]*z[6]-2.0*x[3]*z[0]*x[7]*y[3]-x[6]*z[6]*x[7]*y[4]
       -2.0*x[1]*z[1]*x[5]*y[0];
  s7 = s8+2.0*x[1]*y[1]*x[5]*z[0]+2.0*x[1]*z[1]*x[0]*y[5]+2.0*x[3]*y[0]*x
       [7]*z[3]+2.0*x[3]*x[0]*y[3]*z[7]-2.0*x[3]*x[0]*y[7]*z[3]-2.0*x[1]*y[1]*x[0]*z
       [5]-2.0*x[6]*y[6]*x[5]*z[7]+s6-y[5]*x[1]*x[1]*z[4]+x[6]*z[6]*x[4]*y[7]-2.0*x[2]
       *y[2]*x[3]*z[1]+x[6]*z[5]*x[4]*y[6]+x[6]*x[5]*y[4]*z[6]-y[6]*x[7]*x[7]*z[2]-x
       [6]*x[5]*y[6]*z[4];
  s8 = x[3]*x[3]*y[7]*z[4]-2.0*y[6]*x[7]*x[7]*z[3]+z[6]*x[7]*x[7]*y[2]+2.0*
       z[6]*x[7]*x[7]*y[3]+2.0*y[1]*x[0]*x[0]*z[3]+2.0*x[0]*x[1]*y[3]*z[0]-2.0*x[0]*y
       [0]*x[3]*z[4]-2.0*x[0]*z[1]*x[4]*y[0]-2.0*x[0]*y[1]*x[3]*z[0]+2.0*x[0]*y[0]*x
       [4]*z[3]-2.0*x[0]*z[0]*x[4]*y[3]+2.0*x[0]*x[1]*y[0]*z[4]+2.0*x[0]*z[1]*x[3]*y
       [0]-2.0*x[0]*x[1]*y[0]*z[3]-2.0*x[0]*x[1]*y[4]*z[0]+2.0*x[0]*y[1]*x[4]*z[0];
  s5 = s8+2.0*x[0]*z[0]*x[3]*y[4]+x[1]*y[1]*x[0]*z[3]-x[1]*z[1]*x[4]*y[0]-x
       [1]*y[1]*x[0]*z[4]+x[1]*z[1]*x[0]*y[4]-x[1]*y[1]*x[3]*z[0]-x[1]*z[1]*x[0]*y[3]-
       x[0]*z[5]*x[4]*y[1]+x[0]*y[5]*x[4]*z[1]-2.0*x[4]*x[0]*y[4]*z[7]-2.0*x[4]*y[5]*x
       [0]*z[4]+2.0*x[4]*z[5]*x[0]*y[4]-2.0*x[4]*x[5]*y[4]*z[0]-2.0*x[4]*y[0]*x[7]*z
       [4]-x[5]*y[5]*x[0]*z[4]+s7;
  s8 = x[5]*z[5]*x[0]*y[4]-x[5]*z[5]*x[4]*y[0]+x[1]*z[5]*x[0]*y[4]+x[5]*y
       [5]*x[4]*z[0]-x[0]*y[0]*x[7]*z[4]-x[0]*z[5]*x[4]*y[0]-x[1]*y[5]*x[0]*z[4]+x[0]*
       z[0]*x[7]*y[4]+x[0]*y[5]*x[4]*z[0]-x[0]*z[0]*x[4]*y[7]+x[0]*x[5]*y[0]*z[4]+x[0]
       *y[0]*x[4]*z[7]-x[0]*x[5]*y[4]*z[0]-x[3]*x[3]*y[4]*z[7]+2.0*x[2]*z[2]*x[3]*y[1]
       ;
  s7 = s8-x[5]*x[5]*y[4]*z[0]+2.0*y[5]*x[4]*x[4]*z[0]-2.0*z[0]*x[4]*x[4]*y
       [7]+2.0*y[0]*x[4]*x[4]*z[7]-2.0*z[5]*x[4]*x[4]*y[0]+x[5]*x[5]*y[4]*z[7]-x[5]*x
       [5]*y[7]*z[4]-2.0*y[5]*x[4]*x[4]*z[7]+2.0*z[5]*x[4]*x[4]*y[7]-x[0]*x[0]*y[7]*z
       [3]+y[2]*x[0]*x[0]*z[3]+x[0]*x[0]*y[3]*z[7]-x[5]*x[1]*y[4]*z[0]+x[5]*y[1]*x[4]*
       z[0]-x[4]*y[0]*x[3]*z[4];
  s8 = -x[4]*y[1]*x[0]*z[4]+x[4]*z[1]*x[0]*y[4]+x[4]*x[0]*y[3]*z[4]-x[4]*x
       [0]*y[4]*z[3]+x[4]*x[1]*y[0]*z[4]-x[4]*x[1]*y[4]*z[0]+x[4]*z[0]*x[3]*y[4]+x[5]*
       x[1]*y[0]*z[4]+x[1]*z[1]*x[3]*y[0]+x[1]*y[1]*x[4]*z[0]-x[5]*z[1]*x[4]*y[0]-2.0*
       y[1]*x[0]*x[0]*z[4]+2.0*z[1]*x[0]*x[0]*y[4]+2.0*x[0]*x[0]*y[3]*z[4]-2.0*z[1]*x
       [0]*x[0]*y[3];
  s6 = s8-2.0*x[0]*x[0]*y[4]*z[3]+x[1]*x[1]*y[3]*z[0]+x[1]*x[1]*y[0]*z[4]-x
       [1]*x[1]*y[0]*z[3]-x[1]*x[1]*y[4]*z[0]-z[1]*x[4]*x[4]*y[0]+y[0]*x[4]*x[4]*z[3]-
       z[0]*x[4]*x[4]*y[3]+y[1]*x[4]*x[4]*z[0]-x[0]*x[0]*y[4]*z[7]-y[5]*x[0]*x[0]*z[4]
       +z[5]*x[0]*x[0]*y[4]+x[5]*x[5]*y[0]*z[4]-x[0]*y[0]*x[3]*z[7]+x[0]*z[0]*x[3]*y
       [7]+s7;
  s8 = s6+x[0]*x[2]*y[3]*z[0]-x[0]*x[2]*y[0]*z[3]+x[0]*y[0]*x[7]*z[3]-x[0]*
       y[2]*x[3]*z[0]+x[0]*z[2]*x[3]*y[0]-x[0]*z[0]*x[7]*y[3]+x[1]*x[2]*y[3]*z[0]-z[2]
       *x[0]*x[0]*y[3]+x[3]*z[2]*x[6]*y[3]-x[3]*x[2]*y[3]*z[6]+x[3]*x[2]*y[6]*z[3]-x
       [3]*y[2]*x[6]*z[3]-2.0*x[3]*y[2]*x[7]*z[3]+2.0*x[3]*z[2]*x[7]*y[3];
  s7 = s8+2.0*x[4]*y[5]*x[7]*z[4]+2.0*x[4]*x[5]*y[4]*z[7]-2.0*x[4]*z[5]*x
       [7]*y[4]-2.0*x[4]*x[5]*y[7]*z[4]+x[5]*y[5]*x[7]*z[4]-x[5]*z[5]*x[7]*y[4]-x[5]*y
       [5]*x[4]*z[7]+x[5]*z[5]*x[4]*y[7]+2.0*x[3]*x[2]*y[7]*z[3]-2.0*x[2]*z[2]*x[1]*y
       [3]+2.0*x[4]*z[0]*x[7]*y[4]+2.0*x[4]*x[0]*y[7]*z[4]+2.0*x[4]*x[5]*y[0]*z[4]-x
       [7]*x[6]*y[2]*z[7]-2.0*x[3]*x[2]*y[3]*z[7]-x[0]*x[4]*y[7]*z[3];
  s8 = x[0]*x[3]*y[7]*z[4]-x[0]*x[3]*y[4]*z[7]+x[0]*x[4]*y[3]*z[7]-2.0*x[7]
       *z[6]*x[3]*y[7]+x[3]*x[7]*y[4]*z[3]-x[3]*x[4]*y[7]*z[3]-x[3]*x[7]*y[3]*z[4]+x
       [3]*x[4]*y[3]*z[7]+2.0*x[2]*y[2]*x[1]*z[3]+y[6]*x[3]*x[3]*z[7]-z[6]*x[3]*x[3]*y
       [7]-x[1]*z[5]*x[4]*y[1]-x[1]*x[5]*y[4]*z[1]-x[1]*z[2]*x[0]*y[3]-x[1]*x[2]*y[0]*
       z[3]+x[1]*y[2]*x[0]*z[3];
  s4 = s8+x[1]*x[5]*y[1]*z[4]+x[1]*y[5]*x[4]*z[1]+x[4]*y[0]*x[7]*z[3]-x[4]*
       z[0]*x[7]*y[3]-x[4]*x[4]*y[7]*z[3]+x[4]*x[4]*y[3]*z[7]+x[3]*z[6]*x[7]*y[3]-x[3]
       *x[6]*y[3]*z[7]+x[3]*x[6]*y[7]*z[3]-x[3]*z[6]*x[2]*y[7]-x[3]*y[6]*x[7]*z[3]+x
       [3]*z[6]*x[7]*y[2]+x[3]*y[6]*x[2]*z[7]+2.0*x[5]*z[5]*x[4]*y[6]+s5+s7;
  s8 = s4-2.0*x[5]*z[5]*x[6]*y[4]-x[5]*z[6]*x[7]*y[5]+x[5]*x[6]*y[5]*z[7]-x
       [5]*x[6]*y[7]*z[5]-2.0*x[5]*y[5]*x[4]*z[6]+2.0*x[5]*y[5]*x[6]*z[4]-x[3]*y[6]*x
       [7]*z[2]+x[4]*x[7]*y[4]*z[3]+x[4]*x[3]*y[7]*z[4]-x[4]*x[7]*y[3]*z[4]-x[4]*x[3]*
       y[4]*z[7]-z[1]*x[5]*x[5]*y[0]+y[1]*x[5]*x[5]*z[0]+x[4]*y[6]*x[7]*z[4];
  s7 = s8-x[4]*x[6]*y[7]*z[4]+x[4]*x[6]*y[4]*z[7]-x[4]*z[6]*x[7]*y[4]-x[5]*
       y[6]*x[4]*z[7]-x[5]*x[6]*y[7]*z[4]+x[5]*x[6]*y[4]*z[7]+x[5]*z[6]*x[4]*y[7]-y[6]
       *x[4]*x[4]*z[7]+z[6]*x[4]*x[4]*y[7]+x[7]*x[5]*y[4]*z[7]-y[2]*x[7]*x[7]*z[3]+z
       [2]*x[7]*x[7]*y[3]-y[0]*x[3]*x[3]*z[4]-y[1]*x[3]*x[3]*z[0]+z[1]*x[3]*x[3]*y[0];
  s8 = z[0]*x[3]*x[3]*y[4]-x[2]*y[1]*x[3]*z[0]+x[2]*z[1]*x[3]*y[0]+x[3]*y
       [1]*x[0]*z[3]+x[3]*x[1]*y[3]*z[0]+x[3]*x[0]*y[3]*z[4]-x[3]*z[1]*x[0]*y[3]-x[3]*
       x[0]*y[4]*z[3]+x[3]*y[0]*x[4]*z[3]-x[3]*z[0]*x[4]*y[3]-x[3]*x[1]*y[0]*z[3]+x[3]
       *z[0]*x[7]*y[4]-x[3]*y[0]*x[7]*z[4]+z[0]*x[7]*x[7]*y[4]-y[0]*x[7]*x[7]*z[4];
  s6 = s8+y[1]*x[0]*x[0]*z[2]-2.0*y[2]*x[3]*x[3]*z[0]+2.0*z[2]*x[3]*x[3]*y
       [0]-2.0*x[1]*x[1]*y[0]*z[2]+2.0*x[1]*x[1]*y[2]*z[0]-y[2]*x[3]*x[3]*z[1]+z[2]*x
       [3]*x[3]*y[1]-y[5]*x[4]*x[4]*z[6]+z[5]*x[4]*x[4]*y[6]+x[7]*x[0]*y[7]*z[4]-x[7]*
       z[0]*x[4]*y[7]-x[7]*x[0]*y[4]*z[7]+x[7]*y[0]*x[4]*z[7]-x[0]*x[1]*y[0]*z[2]+x[0]
       *z[1]*x[2]*y[0]+s7;
  s8 = s6+x[0]*x[1]*y[2]*z[0]-x[0]*y[1]*x[2]*z[0]-x[3]*z[1]*x[0]*y[2]+2.0*x
       [3]*x[2]*y[3]*z[0]+y[0]*x[7]*x[7]*z[3]-z[0]*x[7]*x[7]*y[3]-2.0*x[3]*z[2]*x[0]*y
       [3]-2.0*x[3]*x[2]*y[0]*z[3]+2.0*x[3]*y[2]*x[0]*z[3]+x[3]*x[2]*y[3]*z[1]-x[3]*x
       [2]*y[1]*z[3]-x[5]*y[1]*x[0]*z[5]+x[3]*y[1]*x[0]*z[2]+x[4]*y[6]*x[7]*z[5];
  s7 = s8-x[5]*x[1]*y[5]*z[0]+2.0*x[1]*z[1]*x[2]*y[0]-2.0*x[1]*z[1]*x[0]*y
       [2]+x[1]*x[2]*y[3]*z[1]-x[1]*x[2]*y[1]*z[3]+2.0*x[1]*y[1]*x[0]*z[2]-2.0*x[1]*y
       [1]*x[2]*z[0]-z[2]*x[1]*x[1]*y[3]+y[2]*x[1]*x[1]*z[3]+y[5]*x[7]*x[7]*z[4]+y[6]*
       x[7]*x[7]*z[5]+x[7]*x[6]*y[7]*z[2]+x[7]*y[6]*x[2]*z[7]-x[7]*z[6]*x[2]*y[7]-2.0*
       x[7]*x[6]*y[3]*z[7];
  s8 = s7+2.0*x[7]*x[6]*y[7]*z[3]+2.0*x[7]*y[6]*x[3]*z[7]-x[3]*z[2]*x[1]*y
       [3]+x[3]*y[2]*x[1]*z[3]+x[5]*x[1]*y[0]*z[5]+x[4]*y[5]*x[6]*z[4]+x[5]*z[1]*x[0]*
       y[5]-x[4]*z[6]*x[7]*y[5]-x[4]*x[5]*y[6]*z[4]+x[4]*x[5]*y[4]*z[6]-x[4]*z[5]*x[6]
       *y[4]-x[1]*y[2]*x[3]*z[1]+x[1]*z[2]*x[3]*y[1]-x[2]*x[1]*y[0]*z[2]-x[2]*z[1]*x
       [0]*y[2];
  s5 = s8+x[2]*x[1]*y[2]*z[0]-x[2]*z[2]*x[0]*y[3]+x[2]*y[2]*x[0]*z[3]-x[2]*
       y[2]*x[3]*z[0]+x[2]*z[2]*x[3]*y[0]+x[2]*y[1]*x[0]*z[2]+x[5]*y[6]*x[7]*z[5]+x[6]
       *y[5]*x[7]*z[4]+2.0*x[6]*y[6]*x[7]*z[5]-x[7]*y[0]*x[3]*z[7]+x[7]*z[0]*x[3]*y[7]
       -x[7]*x[0]*y[7]*z[3]+x[7]*x[0]*y[3]*z[7]+2.0*x[7]*x[7]*y[4]*z[3]-2.0*x[7]*x[7]*
       y[3]*z[4]-2.0*x[1]*x[1]*y[2]*z[5];
  s8 = s5-2.0*x[7]*x[4]*y[7]*z[3]+2.0*x[7]*x[3]*y[7]*z[4]-2.0*x[7]*x[3]*y
       [4]*z[7]+2.0*x[7]*x[4]*y[3]*z[7]+2.0*x[1]*x[1]*y[5]*z[2]-x[1]*x[1]*y[2]*z[6]+x
       [1]*x[1]*y[6]*z[2]+z[1]*x[5]*x[5]*y[2]-y[1]*x[5]*x[5]*z[2]-x[1]*x[1]*y[6]*z[5]+
       x[1]*x[1]*y[5]*z[6]+x[5]*x[5]*y[6]*z[2]-x[5]*x[5]*y[2]*z[6]-2.0*y[1]*x[5]*x[5]*
       z[6];
  s7 = s8+2.0*z[1]*x[5]*x[5]*y[6]+2.0*x[1]*z[1]*x[5]*y[2]+2.0*x[1]*y[1]*x
       [2]*z[5]-2.0*x[1]*z[1]*x[2]*y[5]-2.0*x[1]*y[1]*x[5]*z[2]-x[1]*y[1]*x[6]*z[2]-x
       [1]*z[1]*x[2]*y[6]+x[1]*z[1]*x[6]*y[2]+x[1]*y[1]*x[2]*z[6]-x[5]*x[1]*y[2]*z[5]+
       x[5]*y[1]*x[2]*z[5]-x[5]*z[1]*x[2]*y[5]+x[5]*x[1]*y[5]*z[2]-x[5]*y[1]*x[6]*z[2]
       -x[5]*x[1]*y[2]*z[6];
  s8 = s7+x[5]*x[1]*y[6]*z[2]+x[5]*z[1]*x[6]*y[2]+x[1]*x[2]*y[5]*z[6]-x[1]*
       x[2]*y[6]*z[5]-x[1]*z[1]*x[6]*y[5]-x[1]*y[1]*x[5]*z[6]+x[1]*z[1]*x[5]*y[6]+x[1]
       *y[1]*x[6]*z[5]-x[5]*x[6]*y[5]*z[2]+x[5]*x[2]*y[5]*z[6]-x[5]*x[2]*y[6]*z[5]+x
       [5]*x[6]*y[2]*z[5]-2.0*x[5]*z[1]*x[6]*y[5]-2.0*x[5]*x[1]*y[6]*z[5]+2.0*x[5]*x
       [1]*y[5]*z[6];
  s6 = s8+2.0*x[5]*y[1]*x[6]*z[5]+2.0*x[2]*x[1]*y[6]*z[2]+2.0*x[2]*z[1]*x
       [6]*y[2]-2.0*x[2]*x[1]*y[2]*z[6]+x[2]*x[5]*y[6]*z[2]+x[2]*x[6]*y[2]*z[5]-x[2]*x
       [5]*y[2]*z[6]+y[1]*x[2]*x[2]*z[5]-z[1]*x[2]*x[2]*y[5]-2.0*x[2]*y[1]*x[6]*z[2]-x
       [2]*x[6]*y[5]*z[2]-2.0*z[1]*x[2]*x[2]*y[6]+x[2]*x[2]*y[5]*z[6]-x[2]*x[2]*y[6]*z
       [5]+2.0*y[1]*x[2]*x[2]*z[6]+x[2]*z[1]*x[5]*y[2];
  s8 = s6-x[2]*x[1]*y[2]*z[5]+x[2]*x[1]*y[5]*z[2]-x[2]*y[1]*x[5]*z[2]+x[6]*
       y[1]*x[2]*z[5]-x[6]*z[1]*x[2]*y[5]-z[1]*x[6]*x[6]*y[5]+y[1]*x[6]*x[6]*z[5]-y[1]
       *x[6]*x[6]*z[2]-2.0*x[6]*x[6]*y[5]*z[2]+2.0*x[6]*x[6]*y[2]*z[5]+z[1]*x[6]*x[6]*
       y[2]-x[6]*x[1]*y[6]*z[5]-x[6]*y[1]*x[5]*z[6]+x[6]*x[1]*y[5]*z[6];
  s7 = s8+x[6]*z[1]*x[5]*y[6]-x[6]*z[1]*x[2]*y[6]-x[6]*x[1]*y[2]*z[6]+2.0*x
       [6]*x[5]*y[6]*z[2]+2.0*x[6]*x[2]*y[5]*z[6]-2.0*x[6]*x[2]*y[6]*z[5]-2.0*x[6]*x
       [5]*y[2]*z[6]+x[6]*x[1]*y[6]*z[2]+x[6]*y[1]*x[2]*z[6]-x[2]*x[2]*y[3]*z[7]+x[2]*
       x[2]*y[7]*z[3]-x[2]*z[2]*x[3]*y[7]-x[2]*y[2]*x[7]*z[3]+x[2]*z[2]*x[7]*y[3]+x[2]
       *y[2]*x[3]*z[7]-x[6]*x[6]*y[3]*z[7];
  s8 = s7+x[6]*x[6]*y[7]*z[3]-x[6]*x[2]*y[3]*z[7]+x[6]*x[2]*y[7]*z[3]-x[6]*
       y[6]*x[7]*z[3]+x[6]*y[6]*x[3]*z[7]-x[6]*z[6]*x[3]*y[7]+x[6]*z[6]*x[7]*y[3]+y[6]
       *x[2]*x[2]*z[7]-z[6]*x[2]*x[2]*y[7]+2.0*x[2]*x[2]*y[6]*z[3]-x[2]*y[6]*x[7]*z[2]
       -2.0*x[2]*y[2]*x[6]*z[3]-2.0*x[2]*x[2]*y[3]*z[6]+2.0*x[2]*y[2]*x[3]*z[6]-x[2]*x
       [6]*y[2]*z[7];
  s3 = s8+x[2]*x[6]*y[7]*z[2]+x[2]*z[6]*x[7]*y[2]+2.0*x[2]*z[2]*x[6]*y[3]
       -2.0*x[2]*z[2]*x[3]*y[6]-y[2]*x[6]*x[6]*z[3]-2.0*x[6]*x[6]*y[2]*z[7]+2.0*x[6]*x
       [6]*y[7]*z[2]+z[2]*x[6]*x[6]*y[3]-2.0*x[6]*y[6]*x[7]*z[2]+x[6]*y[2]*x[3]*z[6]-x
       [6]*x[2]*y[3]*z[6]+2.0*x[6]*z[6]*x[7]*y[2]+2.0*x[6]*y[6]*x[2]*z[7]-2.0*x[6]*z
       [6]*x[2]*y[7]+x[6]*x[2]*y[6]*z[3]-x[6]*z[2]*x[3]*y[6];
  s8 = y[1]*x[0]*z[3]+x[1]*y[3]*z[0]-y[0]*x[3]*z[7]-x[1]*y[5]*z[0]-y[0]*x
       [3]*z[4]-x[1]*y[0]*z[2]+z[1]*x[2]*y[0]-y[1]*x[0]*z[5]-z[1]*x[0]*y[2]-y[1]*x[0]*
       z[4]+z[1]*x[5]*y[2]+z[0]*x[7]*y[4]+z[0]*x[3]*y[7]+z[1]*x[0]*y[4]-x[1]*y[2]*z[5]
       +x[2]*y[3]*z[0]+y[1]*x[2]*z[5]-x[2]*y[3]*z[7];
  s7 = s8-z[1]*x[2]*y[5]-y[1]*x[3]*z[0]-x[0]*y[7]*z[3]-z[1]*x[0]*y[3]+y[5]*
       x[4]*z[0]-x[0]*y[4]*z[3]+y[5]*x[7]*z[4]-z[0]*x[4]*y[3]+x[1]*y[0]*z[4]-z[2]*x[3]
       *y[7]-y[6]*x[7]*z[2]+x[1]*y[5]*z[2]+y[6]*x[7]*z[5]+x[0]*y[7]*z[4]+x[1]*y[2]*z
       [0]-z[1]*x[4]*y[0]-z[0]*x[4]*y[7]-z[2]*x[0]*y[3];
  s8 = x[5]*y[0]*z[4]+z[1]*x[0]*y[5]-x[2]*y[0]*z[3]-z[1]*x[5]*y[0]+y[1]*x
       [5]*z[0]-x[1]*y[0]*z[3]-x[1]*y[4]*z[0]-y[1]*x[5]*z[2]+x[2]*y[7]*z[3]+y[0]*x[4]*
       z[3]-x[0]*y[4]*z[7]+x[1]*y[0]*z[5]-y[1]*x[6]*z[2]-y[2]*x[6]*z[3]+y[0]*x[7]*z[3]
       -y[2]*x[7]*z[3]+z[2]*x[7]*y[3]+y[2]*x[0]*z[3];
  s6 = s8+y[2]*x[3]*z[7]-y[2]*x[3]*z[0]-x[6]*y[5]*z[2]-y[5]*x[0]*z[4]+z[2]*
       x[3]*y[0]+x[2]*y[3]*z[1]+x[0]*y[3]*z[7]-x[2]*y[1]*z[3]+y[1]*x[4]*z[0]+y[1]*x[0]
       *z[2]-z[1]*x[2]*y[6]+y[2]*x[3]*z[6]-y[1]*x[2]*z[0]+z[1]*x[3]*y[0]-x[1]*y[2]*z
       [6]-x[2]*y[3]*z[6]+x[0]*y[3]*z[4]+z[0]*x[3]*y[4]+s7;
  s8 = x[5]*y[4]*z[7]+s6+y[5]*x[6]*z[4]-y[5]*x[4]*z[6]+z[6]*x[5]*y[7]-x[6]*
       y[2]*z[7]-x[6]*y[7]*z[5]+x[5]*y[6]*z[2]+x[6]*y[5]*z[7]+x[6]*y[7]*z[2]+y[6]*x[7]
       *z[4]-y[6]*x[4]*z[7]-y[6]*x[7]*z[3]+z[6]*x[7]*y[2]+x[2]*y[5]*z[6]-x[2]*y[6]*z
       [5]+y[6]*x[2]*z[7]+x[6]*y[2]*z[5];
  s7 = s8-x[5]*y[2]*z[6]-z[6]*x[7]*y[5]-z[5]*x[7]*y[4]+z[5]*x[0]*y[4]-y[5]*
       x[4]*z[7]+y[0]*x[4]*z[7]-z[6]*x[2]*y[7]-x[5]*y[4]*z[0]-x[5]*y[7]*z[4]-y[0]*x[7]
       *z[4]+y[5]*x[4]*z[1]-x[6]*y[7]*z[4]+x[7]*y[4]*z[3]-x[4]*y[7]*z[3]+x[3]*y[7]*z
       [4]-x[7]*y[3]*z[4]-x[6]*y[3]*z[7]+x[6]*y[4]*z[7];
  s8 = -x[3]*y[4]*z[7]+x[4]*y[3]*z[7]-z[6]*x[7]*y[4]-z[1]*x[6]*y[5]+x[6]*y
       [7]*z[3]-x[1]*y[6]*z[5]-y[1]*x[5]*z[6]+z[5]*x[4]*y[7]-z[5]*x[4]*y[0]+x[1]*y[5]*
       z[6]-y[6]*x[5]*z[7]-y[2]*x[3]*z[1]+z[1]*x[5]*y[6]-y[5]*x[1]*z[4]+z[6]*x[4]*y[7]
       +x[5]*y[1]*z[4]-x[5]*y[6]*z[4]+y[6]*x[3]*z[7]-x[5]*y[4]*z[1];
  s5 = s8+x[5]*y[4]*z[6]+z[5]*x[1]*y[4]+y[1]*x[6]*z[5]-z[6]*x[3]*y[7]+z[6]*
       x[7]*y[3]-z[5]*x[6]*y[4]-z[5]*x[4]*y[1]+z[5]*x[4]*y[6]+x[1]*y[6]*z[2]+x[2]*y[6]
       *z[3]+z[2]*x[6]*y[3]+z[1]*x[6]*y[2]+z[2]*x[3]*y[1]-z[2]*x[1]*y[3]-z[2]*x[3]*y
       [6]+y[2]*x[1]*z[3]+y[1]*x[2]*z[6]-z[0]*x[7]*y[3]+s7;
  s4 = 1/s5;
  s2 = s3*s4;
  const double unknown0 = s1*s2;
  s1 = 1.0/6.0;
  s8 = 2.0*x[1]*y[0]*y[0]*z[4]+x[5]*y[0]*y[0]*z[4]-x[1]*y[4]*y[4]*z[0]+z[1]
       *x[0]*y[4]*y[4]+x[1]*y[0]*y[0]*z[5]-z[1]*x[5]*y[0]*y[0]-2.0*z[1]*x[4]*y[0]*y[0]
       +2.0*z[1]*x[3]*y[0]*y[0]+z[2]*x[3]*y[0]*y[0]+y[0]*y[0]*x[7]*z[3]+2.0*y[0]*y[0]*
       x[4]*z[3]-2.0*x[1]*y[0]*y[0]*z[3]-2.0*x[5]*y[4]*y[4]*z[0]+2.0*z[5]*x[0]*y[4]*y
       [4]+2.0*y[4]*y[5]*x[7]*z[4];
  s7 = s8-x[3]*y[4]*y[4]*z[7]+x[7]*y[4]*y[4]*z[3]+z[0]*x[3]*y[4]*y[4]-2.0*x
       [0]*y[4]*y[4]*z[7]-y[1]*x[1]*y[4]*z[0]-x[0]*y[4]*y[4]*z[3]+2.0*z[0]*x[7]*y[4]*y
       [4]+y[4]*z[6]*x[4]*y[7]-y[0]*y[0]*x[7]*z[4]+y[0]*y[0]*x[4]*z[7]+2.0*y[4]*z[5]*x
       [4]*y[7]-2.0*y[4]*x[5]*y[7]*z[4]-y[4]*x[6]*y[7]*z[4]-y[4]*y[6]*x[4]*z[7]-2.0*y
       [4]*y[5]*x[4]*z[7];
  s8 = y[4]*y[6]*x[7]*z[4]-y[7]*y[2]*x[7]*z[3]+y[7]*z[2]*x[7]*y[3]+y[7]*y
       [2]*x[3]*z[7]+2.0*x[5]*y[4]*y[4]*z[7]-y[7]*x[2]*y[3]*z[7]-y[0]*z[0]*x[4]*y[7]+z
       [6]*x[7]*y[3]*y[3]-y[0]*x[0]*y[4]*z[7]+y[0]*x[0]*y[7]*z[4]-2.0*x[2]*y[3]*y[3]*z
       [7]-z[5]*x[4]*y[0]*y[0]+y[0]*z[0]*x[7]*y[4]-2.0*z[6]*x[3]*y[7]*y[7]+z[1]*x[2]*y
       [0]*y[0];
  s6 = s8+y[4]*y[0]*x[4]*z[3]-2.0*y[4]*z[0]*x[4]*y[7]+2.0*y[4]*x[0]*y[7]*z
       [4]-y[4]*z[0]*x[4]*y[3]-y[4]*x[0]*y[7]*z[3]+y[4]*z[0]*x[3]*y[7]-y[4]*y[0]*x[3]*
       z[4]+y[0]*x[4]*y[3]*z[7]-y[0]*x[7]*y[3]*z[4]-y[0]*x[3]*y[4]*z[7]+y[0]*x[7]*y[4]
       *z[3]+x[2]*y[7]*y[7]*z[3]-z[2]*x[3]*y[7]*y[7]-2.0*z[2]*x[0]*y[3]*y[3]+2.0*y[0]*
       z[1]*x[0]*y[4]+s7;
  s8 = -2.0*y[0]*y[1]*x[0]*z[4]-y[0]*y[1]*x[0]*z[5]-y[0]*y[0]*x[3]*z[7]-z
       [1]*x[0]*y[3]*y[3]-y[0]*x[1]*y[5]*z[0]-2.0*z[0]*x[7]*y[3]*y[3]+x[0]*y[3]*y[3]*z
       [4]+2.0*x[0]*y[3]*y[3]*z[7]-z[0]*x[4]*y[3]*y[3]+2.0*x[2]*y[3]*y[3]*z[0]+x[1]*y
       [3]*y[3]*z[0]+2.0*y[7]*z[6]*x[7]*y[3]+2.0*y[7]*y[6]*x[3]*z[7]-2.0*y[7]*y[6]*x
       [7]*z[3]-2.0*y[7]*x[6]*y[3]*z[7];
  s7 = s8+y[4]*x[4]*y[3]*z[7]-y[4]*x[4]*y[7]*z[3]+y[4]*x[3]*y[7]*z[4]-y[4]*
       x[7]*y[3]*z[4]+2.0*y[4]*y[0]*x[4]*z[7]-2.0*y[4]*y[0]*x[7]*z[4]+2.0*x[6]*y[7]*y
       [7]*z[3]+y[4]*x[0]*y[3]*z[4]+y[0]*y[1]*x[5]*z[0]+y[0]*z[1]*x[0]*y[5]-x[2]*y[0]*
       y[0]*z[3]+x[4]*y[3]*y[3]*z[7]-x[7]*y[3]*y[3]*z[4]-x[5]*y[4]*y[4]*z[1]+y[3]*z[0]
       *x[3]*y[4];
  s8 = y[3]*y[0]*x[4]*z[3]+2.0*y[3]*y[0]*x[7]*z[3]+2.0*y[3]*y[2]*x[0]*z[3]
       -2.0*y[3]*y[2]*x[3]*z[0]+2.0*y[3]*z[2]*x[3]*y[0]+y[3]*z[1]*x[3]*y[0]-2.0*y[3]*x
       [2]*y[0]*z[3]-y[3]*x[1]*y[0]*z[3]-y[3]*y[1]*x[3]*z[0]-2.0*y[3]*x[0]*y[7]*z[3]-y
       [3]*x[0]*y[4]*z[3]-2.0*y[3]*y[0]*x[3]*z[7]-y[3]*y[0]*x[3]*z[4]+2.0*y[3]*z[0]*x
       [3]*y[7]+y[3]*y[1]*x[0]*z[3]+z[5]*x[1]*y[4]*y[4];
  s5 = s8-2.0*y[0]*y[0]*x[3]*z[4]-2.0*y[0]*x[1]*y[4]*z[0]+y[3]*x[7]*y[4]*z
       [3]-y[3]*x[4]*y[7]*z[3]+y[3]*x[3]*y[7]*z[4]-y[3]*x[3]*y[4]*z[7]+y[3]*x[0]*y[7]*
       z[4]-y[3]*z[0]*x[4]*y[7]-2.0*y[4]*y[5]*x[0]*z[4]+s6+y[7]*x[0]*y[3]*z[7]-y[7]*z
       [0]*x[7]*y[3]+y[7]*y[0]*x[7]*z[3]-y[7]*y[0]*x[3]*z[7]+2.0*y[0]*y[1]*x[4]*z[0]+
       s7;
  s8 = -2.0*y[7]*x[7]*y[3]*z[4]-2.0*y[7]*x[3]*y[4]*z[7]+2.0*y[7]*x[4]*y[3]*
       z[7]+y[7]*y[0]*x[4]*z[7]-y[7]*y[0]*x[7]*z[4]+2.0*y[7]*x[7]*y[4]*z[3]-y[7]*x[0]*
       y[4]*z[7]+y[7]*z[0]*x[7]*y[4]+z[5]*x[4]*y[7]*y[7]+2.0*z[6]*x[4]*y[7]*y[7]-x[5]*
       y[7]*y[7]*z[4]-2.0*x[6]*y[7]*y[7]*z[4]+2.0*y[7]*x[6]*y[4]*z[7]-2.0*y[7]*z[6]*x
       [7]*y[4]+2.0*y[7]*y[6]*x[7]*z[4];
  s7 = s8-2.0*y[7]*y[6]*x[4]*z[7]-y[7]*z[5]*x[7]*y[4]-y[7]*y[5]*x[4]*z[7]-x
       [0]*y[7]*y[7]*z[3]+z[0]*x[3]*y[7]*y[7]+y[7]*x[5]*y[4]*z[7]+y[7]*y[5]*x[7]*z[4]-
       y[4]*x[1]*y[5]*z[0]-x[1]*y[0]*y[0]*z[2]-y[4]*y[5]*x[1]*z[4]-2.0*y[4]*z[5]*x[4]*
       y[0]-y[4]*y[1]*x[0]*z[4]+y[4]*y[5]*x[4]*z[1]+y[0]*z[0]*x[3]*y[7]-y[0]*z[1]*x[0]
       *y[2];
  s8 = 2.0*y[0]*x[1]*y[3]*z[0]+y[4]*y[1]*x[4]*z[0]+2.0*y[0]*y[1]*x[0]*z[3]+
       y[4]*x[1]*y[0]*z[5]-y[4]*z[1]*x[5]*y[0]+y[4]*z[1]*x[0]*y[5]-y[4]*z[1]*x[4]*y[0]
       +y[4]*x[1]*y[0]*z[4]-y[4]*z[5]*x[4]*y[1]+x[5]*y[4]*y[4]*z[6]-z[5]*x[6]*y[4]*y
       [4]+y[4]*x[5]*y[1]*z[4]-y[0]*z[2]*x[0]*y[3]+y[0]*y[5]*x[4]*z[0]+y[0]*x[1]*y[2]*
       z[0];
  s6 = s8-2.0*y[0]*z[0]*x[4]*y[3]-2.0*y[0]*x[0]*y[4]*z[3]-2.0*y[0]*z[1]*x
       [0]*y[3]-y[0]*x[0]*y[7]*z[3]-2.0*y[0]*y[1]*x[3]*z[0]+y[0]*x[2]*y[3]*z[0]-y[0]*y
       [1]*x[2]*z[0]+y[0]*y[1]*x[0]*z[2]-y[0]*x[2]*y[1]*z[3]+y[0]*x[0]*y[3]*z[7]+y[0]*
       x[2]*y[3]*z[1]-y[0]*y[2]*x[3]*z[0]+y[0]*y[2]*x[0]*z[3]-y[0]*y[5]*x[0]*z[4]-y[4]
       *y[5]*x[4]*z[6]+s7;
  s8 = s6+y[4]*z[6]*x[5]*y[7]-y[4]*x[6]*y[7]*z[5]+y[4]*x[6]*y[5]*z[7]-y[4]*
       z[6]*x[7]*y[5]-y[4]*x[5]*y[6]*z[4]+y[4]*z[5]*x[4]*y[6]+y[4]*y[5]*x[6]*z[4]-2.0*
       y[1]*y[1]*x[0]*z[5]+2.0*y[1]*y[1]*x[5]*z[0]-2.0*y[2]*y[2]*x[6]*z[3]+x[5]*y[1]*y
       [1]*z[4]-z[5]*x[4]*y[1]*y[1]-x[6]*y[2]*y[2]*z[7]+z[6]*x[7]*y[2]*y[2];
  s7 = s8-x[1]*y[5]*y[5]*z[0]+z[1]*x[0]*y[5]*y[5]+y[1]*y[5]*x[4]*z[1]-y[1]*
       y[5]*x[1]*z[4]-2.0*y[2]*z[2]*x[3]*y[6]+2.0*y[1]*z[1]*x[0]*y[5]-2.0*y[1]*z[1]*x
       [5]*y[0]+2.0*y[1]*x[1]*y[0]*z[5]-y[2]*x[2]*y[3]*z[7]-y[2]*z[2]*x[3]*y[7]+y[2]*x
       [2]*y[7]*z[3]+y[2]*z[2]*x[7]*y[3]-2.0*y[2]*x[2]*y[3]*z[6]+2.0*y[2]*x[2]*y[6]*z
       [3]+2.0*y[2]*z[2]*x[6]*y[3]-y[3]*y[2]*x[6]*z[3];
  s8 = y[3]*y[2]*x[3]*z[6]+y[3]*x[2]*y[6]*z[3]-y[3]*z[2]*x[3]*y[6]-y[2]*y
       [2]*x[7]*z[3]+2.0*y[2]*y[2]*x[3]*z[6]+y[2]*y[2]*x[3]*z[7]-2.0*y[1]*x[1]*y[5]*z
       [0]-x[2]*y[3]*y[3]*z[6]+z[2]*x[6]*y[3]*y[3]+2.0*y[6]*x[2]*y[5]*z[6]+2.0*y[6]*x
       [6]*y[2]*z[5]-2.0*y[6]*x[5]*y[2]*z[6]+2.0*y[3]*x[2]*y[7]*z[3]-2.0*y[3]*z[2]*x
       [3]*y[7]-y[0]*z[0]*x[7]*y[3]-y[0]*z[2]*x[1]*y[3];
  s4 = s8-y[2]*y[6]*x[7]*z[2]+y[0]*z[2]*x[3]*y[1]+y[1]*z[5]*x[1]*y[4]-y[1]*
       x[5]*y[4]*z[1]+2.0*y[0]*z[0]*x[3]*y[4]+2.0*y[0]*x[0]*y[3]*z[4]+2.0*z[2]*x[7]*y
       [3]*y[3]-2.0*z[5]*x[7]*y[4]*y[4]+x[6]*y[4]*y[4]*z[7]-z[6]*x[7]*y[4]*y[4]+y[1]*y
       [1]*x[0]*z[3]+y[3]*x[6]*y[7]*z[2]-y[3]*z[6]*x[2]*y[7]+2.0*y[3]*y[2]*x[3]*z[7]+
       s5+s7;
  s8 = s4+y[2]*x[6]*y[7]*z[2]-y[2]*y[6]*x[7]*z[3]+y[2]*y[6]*x[2]*z[7]-y[2]*
       z[6]*x[2]*y[7]-y[2]*x[6]*y[3]*z[7]+y[2]*y[6]*x[3]*z[7]+y[2]*z[6]*x[7]*y[3]-2.0*
       y[3]*y[2]*x[7]*z[3]-x[6]*y[3]*y[3]*z[7]+y[1]*y[1]*x[4]*z[0]-y[1]*y[1]*x[3]*z[0]
       +x[2]*y[6]*y[6]*z[3]-z[2]*x[3]*y[6]*y[6]-y[1]*y[1]*x[0]*z[4];
  s7 = s8+y[5]*x[1]*y[0]*z[5]+y[6]*x[2]*y[7]*z[3]-y[6]*y[2]*x[6]*z[3]+y[6]*
       y[2]*x[3]*z[6]-y[6]*x[2]*y[3]*z[6]+y[6]*z[2]*x[6]*y[3]-y[5]*y[1]*x[0]*z[5]-y[5]
       *z[1]*x[5]*y[0]+y[5]*y[1]*x[5]*z[0]-y[6]*z[2]*x[3]*y[7]-y[7]*y[6]*x[7]*z[2]+2.0
       *y[6]*y[6]*x[2]*z[7]+y[6]*y[6]*x[3]*z[7]+x[6]*y[7]*y[7]*z[2]-z[6]*x[2]*y[7]*y
       [7];
  s8 = -x[2]*y[1]*y[1]*z[3]+2.0*y[1]*y[1]*x[0]*z[2]-2.0*y[1]*y[1]*x[2]*z[0]
       +z[2]*x[3]*y[1]*y[1]-z[1]*x[0]*y[2]*y[2]+x[1]*y[2]*y[2]*z[0]+y[2]*y[2]*x[0]*z
       [3]-y[2]*y[2]*x[3]*z[0]-2.0*y[2]*y[2]*x[3]*z[1]+y[1]*x[1]*y[3]*z[0]-2.0*y[6]*y
       [6]*x[7]*z[2]+2.0*y[5]*y[5]*x[4]*z[1]-2.0*y[5]*y[5]*x[1]*z[4]-y[6]*y[6]*x[7]*z
       [3]-2.0*y[1]*x[1]*y[0]*z[2];
  s6 = s8+2.0*y[1]*z[1]*x[2]*y[0]-2.0*y[1]*z[1]*x[0]*y[2]+2.0*y[1]*x[1]*y
       [2]*z[0]+y[1]*x[2]*y[3]*z[1]-y[1]*y[2]*x[3]*z[1]-y[1]*z[2]*x[1]*y[3]+y[1]*y[2]*
       x[1]*z[3]-y[2]*x[1]*y[0]*z[2]+y[2]*z[1]*x[2]*y[0]+y[2]*x[2]*y[3]*z[0]-y[7]*x[6]
       *y[2]*z[7]+y[7]*z[6]*x[7]*y[2]+y[7]*y[6]*x[2]*z[7]-y[6]*x[6]*y[3]*z[7]+y[6]*x
       [6]*y[7]*z[3]+s7;
  s8 = s6-y[6]*z[6]*x[3]*y[7]+y[6]*z[6]*x[7]*y[3]+2.0*y[2]*y[2]*x[1]*z[3]+x
       [2]*y[3]*y[3]*z[1]-z[2]*x[1]*y[3]*y[3]+y[1]*x[1]*y[0]*z[4]+y[1]*z[1]*x[3]*y[0]-
       y[1]*x[1]*y[0]*z[3]+2.0*y[5]*x[5]*y[1]*z[4]-2.0*y[5]*x[5]*y[4]*z[1]+2.0*y[5]*z
       [5]*x[1]*y[4]-2.0*y[5]*z[5]*x[4]*y[1]-2.0*y[6]*x[6]*y[2]*z[7]+2.0*y[6]*x[6]*y
       [7]*z[2];
  s7 = s8+2.0*y[6]*z[6]*x[7]*y[2]-2.0*y[6]*z[6]*x[2]*y[7]-y[1]*z[1]*x[4]*y
       [0]+y[1]*z[1]*x[0]*y[4]-y[1]*z[1]*x[0]*y[3]+2.0*y[6]*y[6]*x[7]*z[5]+2.0*y[5]*y
       [5]*x[6]*z[4]-2.0*y[5]*y[5]*x[4]*z[6]+x[6]*y[5]*y[5]*z[7]-y[3]*x[2]*y[1]*z[3]-y
       [3]*y[2]*x[3]*z[1]+y[3]*z[2]*x[3]*y[1]+y[3]*y[2]*x[1]*z[3]-y[2]*x[2]*y[0]*z[3]+
       y[2]*z[2]*x[3]*y[0];
  s8 = s7+2.0*y[2]*x[2]*y[3]*z[1]-2.0*y[2]*x[2]*y[1]*z[3]+y[2]*y[1]*x[0]*z
       [2]-y[2]*y[1]*x[2]*z[0]+2.0*y[2]*z[2]*x[3]*y[1]-2.0*y[2]*z[2]*x[1]*y[3]-y[2]*z
       [2]*x[0]*y[3]+y[5]*z[6]*x[5]*y[7]-y[5]*x[6]*y[7]*z[5]-y[5]*y[6]*x[4]*z[7]-y[5]*
       y[6]*x[5]*z[7]-2.0*y[5]*x[5]*y[6]*z[4]+2.0*y[5]*x[5]*y[4]*z[6]-2.0*y[5]*z[5]*x
       [6]*y[4]+2.0*y[5]*z[5]*x[4]*y[6];
  s5 = s8-y[1]*y[5]*x[0]*z[4]-z[6]*x[7]*y[5]*y[5]+y[6]*y[6]*x[7]*z[4]-y[6]*
       y[6]*x[4]*z[7]-2.0*y[6]*y[6]*x[5]*z[7]-x[5]*y[6]*y[6]*z[4]+z[5]*x[4]*y[6]*y[6]+
       z[6]*x[5]*y[7]*y[7]-x[6]*y[7]*y[7]*z[5]+y[1]*y[5]*x[4]*z[0]+y[7]*y[6]*x[7]*z[5]
       +y[6]*y[5]*x[7]*z[4]+y[5]*y[6]*x[7]*z[5]+y[6]*y[5]*x[6]*z[4]-y[6]*y[5]*x[4]*z
       [6]+2.0*y[6]*z[6]*x[5]*y[7];
  s8 = s5-2.0*y[6]*x[6]*y[7]*z[5]+2.0*y[6]*x[6]*y[5]*z[7]-2.0*y[6]*z[6]*x
       [7]*y[5]-y[6]*x[5]*y[7]*z[4]-y[6]*x[6]*y[7]*z[4]+y[6]*x[6]*y[4]*z[7]-y[6]*z[6]*
       x[7]*y[4]+y[6]*z[5]*x[4]*y[7]+y[6]*z[6]*x[4]*y[7]+y[6]*x[5]*y[4]*z[6]-y[6]*z[5]
       *x[6]*y[4]+y[7]*x[6]*y[5]*z[7]-y[7]*z[6]*x[7]*y[5]-2.0*y[6]*x[6]*y[5]*z[2];
  s7 = s8-y[7]*y[6]*x[5]*z[7]+2.0*y[4]*y[5]*x[4]*z[0]+2.0*x[3]*y[7]*y[7]*z
       [4]-2.0*x[4]*y[7]*y[7]*z[3]-z[0]*x[4]*y[7]*y[7]+x[0]*y[7]*y[7]*z[4]-y[0]*z[5]*x
       [4]*y[1]+y[0]*x[5]*y[1]*z[4]-y[0]*x[5]*y[4]*z[0]+y[0]*z[5]*x[0]*y[4]-y[5]*y[5]*
       x[0]*z[4]+y[5]*y[5]*x[4]*z[0]+2.0*y[1]*y[1]*x[2]*z[5]-2.0*y[1]*y[1]*x[5]*z[2]+z
       [1]*x[5]*y[2]*y[2];
  s8 = s7-x[1]*y[2]*y[2]*z[5]-y[5]*z[5]*x[4]*y[0]+y[5]*z[5]*x[0]*y[4]-y[5]*
       x[5]*y[4]*z[0]-y[2]*x[1]*y[6]*z[5]-y[2]*y[1]*x[5]*z[6]+y[2]*z[1]*x[5]*y[6]+y[2]
       *y[1]*x[6]*z[5]-y[1]*z[1]*x[6]*y[5]-y[1]*x[1]*y[6]*z[5]+y[1]*x[1]*y[5]*z[6]+y
       [1]*z[1]*x[5]*y[6]+y[5]*x[5]*y[0]*z[4]+y[2]*y[1]*x[2]*z[5]-y[2]*z[1]*x[2]*y[5];
  s6 = s8+y[2]*x[1]*y[5]*z[2]-y[2]*y[1]*x[5]*z[2]-y[1]*y[1]*x[5]*z[6]+y[1]*
       y[1]*x[6]*z[5]-z[1]*x[2]*y[5]*y[5]+x[1]*y[5]*y[5]*z[2]+2.0*y[1]*z[1]*x[5]*y[2]
       -2.0*y[1]*x[1]*y[2]*z[5]-2.0*y[1]*z[1]*x[2]*y[5]+2.0*y[1]*x[1]*y[5]*z[2]-y[1]*y
       [1]*x[6]*z[2]+y[1]*y[1]*x[2]*z[6]-2.0*y[5]*x[1]*y[6]*z[5]-2.0*y[5]*y[1]*x[5]*z
       [6]+2.0*y[5]*z[1]*x[5]*y[6]+2.0*y[5]*y[1]*x[6]*z[5];
  s8 = s6-y[6]*z[1]*x[6]*y[5]-y[6]*y[1]*x[5]*z[6]+y[6]*x[1]*y[5]*z[6]+y[6]*
       y[1]*x[6]*z[5]-2.0*z[1]*x[6]*y[5]*y[5]+2.0*x[1]*y[5]*y[5]*z[6]-x[1]*y[6]*y[6]*z
       [5]+z[1]*x[5]*y[6]*y[6]+y[5]*z[1]*x[5]*y[2]-y[5]*x[1]*y[2]*z[5]+y[5]*y[1]*x[2]*
       z[5]-y[5]*y[1]*x[5]*z[2]-y[6]*z[1]*x[2]*y[5]+y[6]*x[1]*y[5]*z[2];
  s7 = s8-y[1]*z[1]*x[2]*y[6]-y[1]*x[1]*y[2]*z[6]+y[1]*x[1]*y[6]*z[2]+y[1]*
       z[1]*x[6]*y[2]+y[5]*x[5]*y[6]*z[2]-y[5]*x[2]*y[6]*z[5]+y[5]*x[6]*y[2]*z[5]-y[5]
       *x[5]*y[2]*z[6]-x[6]*y[5]*y[5]*z[2]+x[2]*y[5]*y[5]*z[6]-y[5]*y[5]*x[4]*z[7]+y
       [5]*y[5]*x[7]*z[4]-y[1]*x[6]*y[5]*z[2]+y[1]*x[2]*y[5]*z[6]-y[2]*x[6]*y[5]*z[2]
       -2.0*y[2]*y[1]*x[6]*z[2];
  s8 = s7-2.0*y[2]*z[1]*x[2]*y[6]+2.0*y[2]*x[1]*y[6]*z[2]+2.0*y[2]*y[1]*x
       [2]*z[6]-2.0*x[1]*y[2]*y[2]*z[6]+2.0*z[1]*x[6]*y[2]*y[2]+x[6]*y[2]*y[2]*z[5]-x
       [5]*y[2]*y[2]*z[6]+2.0*x[5]*y[6]*y[6]*z[2]-2.0*x[2]*y[6]*y[6]*z[5]-z[1]*x[2]*y
       [6]*y[6]-y[6]*y[1]*x[6]*z[2]-y[6]*x[1]*y[2]*z[6]+y[6]*z[1]*x[6]*y[2]+y[6]*y[1]*
       x[2]*z[6]+x[1]*y[6]*y[6]*z[2];
  s3 = s8+y[2]*x[5]*y[6]*z[2]+y[2]*x[2]*y[5]*z[6]-y[2]*x[2]*y[6]*z[5]+y[5]*
       z[5]*x[4]*y[7]+y[5]*x[5]*y[4]*z[7]-y[5]*z[5]*x[7]*y[4]-y[5]*x[5]*y[7]*z[4]+2.0*
       y[4]*x[5]*y[0]*z[4]-y[3]*z[6]*x[3]*y[7]+y[3]*y[6]*x[3]*z[7]+y[3]*x[6]*y[7]*z[3]
       -y[3]*y[6]*x[7]*z[3]-y[2]*y[1]*x[3]*z[0]-y[2]*z[1]*x[0]*y[3]+y[2]*y[1]*x[0]*z
       [3]+y[2]*x[1]*y[3]*z[0];
  s8 = y[1]*x[0]*z[3]+x[1]*y[3]*z[0]-y[0]*x[3]*z[7]-x[1]*y[5]*z[0]-y[0]*x
       [3]*z[4]-x[1]*y[0]*z[2]+z[1]*x[2]*y[0]-y[1]*x[0]*z[5]-z[1]*x[0]*y[2]-y[1]*x[0]*
       z[4]+z[1]*x[5]*y[2]+z[0]*x[7]*y[4]+z[0]*x[3]*y[7]+z[1]*x[0]*y[4]-x[1]*y[2]*z[5]
       +x[2]*y[3]*z[0]+y[1]*x[2]*z[5]-x[2]*y[3]*z[7];
  s7 = s8-z[1]*x[2]*y[5]-y[1]*x[3]*z[0]-x[0]*y[7]*z[3]-z[1]*x[0]*y[3]+y[5]*
       x[4]*z[0]-x[0]*y[4]*z[3]+y[5]*x[7]*z[4]-z[0]*x[4]*y[3]+x[1]*y[0]*z[4]-z[2]*x[3]
       *y[7]-y[6]*x[7]*z[2]+x[1]*y[5]*z[2]+y[6]*x[7]*z[5]+x[0]*y[7]*z[4]+x[1]*y[2]*z
       [0]-z[1]*x[4]*y[0]-z[0]*x[4]*y[7]-z[2]*x[0]*y[3];
  s8 = x[5]*y[0]*z[4]+z[1]*x[0]*y[5]-x[2]*y[0]*z[3]-z[1]*x[5]*y[0]+y[1]*x
       [5]*z[0]-x[1]*y[0]*z[3]-x[1]*y[4]*z[0]-y[1]*x[5]*z[2]+x[2]*y[7]*z[3]+y[0]*x[4]*
       z[3]-x[0]*y[4]*z[7]+x[1]*y[0]*z[5]-y[1]*x[6]*z[2]-y[2]*x[6]*z[3]+y[0]*x[7]*z[3]
       -y[2]*x[7]*z[3]+z[2]*x[7]*y[3]+y[2]*x[0]*z[3];
  s6 = s8+y[2]*x[3]*z[7]-y[2]*x[3]*z[0]-x[6]*y[5]*z[2]-y[5]*x[0]*z[4]+z[2]*
       x[3]*y[0]+x[2]*y[3]*z[1]+x[0]*y[3]*z[7]-x[2]*y[1]*z[3]+y[1]*x[4]*z[0]+y[1]*x[0]
       *z[2]-z[1]*x[2]*y[6]+y[2]*x[3]*z[6]-y[1]*x[2]*z[0]+z[1]*x[3]*y[0]-x[1]*y[2]*z
       [6]-x[2]*y[3]*z[6]+x[0]*y[3]*z[4]+z[0]*x[3]*y[4]+s7;
  s8 = x[5]*y[4]*z[7]+s6+y[5]*x[6]*z[4]-y[5]*x[4]*z[6]+z[6]*x[5]*y[7]-x[6]*
       y[2]*z[7]-x[6]*y[7]*z[5]+x[5]*y[6]*z[2]+x[6]*y[5]*z[7]+x[6]*y[7]*z[2]+y[6]*x[7]
       *z[4]-y[6]*x[4]*z[7]-y[6]*x[7]*z[3]+z[6]*x[7]*y[2]+x[2]*y[5]*z[6]-x[2]*y[6]*z
       [5]+y[6]*x[2]*z[7]+x[6]*y[2]*z[5];
  s7 = s8-x[5]*y[2]*z[6]-z[6]*x[7]*y[5]-z[5]*x[7]*y[4]+z[5]*x[0]*y[4]-y[5]*
       x[4]*z[7]+y[0]*x[4]*z[7]-z[6]*x[2]*y[7]-x[5]*y[4]*z[0]-x[5]*y[7]*z[4]-y[0]*x[7]
       *z[4]+y[5]*x[4]*z[1]-x[6]*y[7]*z[4]+x[7]*y[4]*z[3]-x[4]*y[7]*z[3]+x[3]*y[7]*z
       [4]-x[7]*y[3]*z[4]-x[6]*y[3]*z[7]+x[6]*y[4]*z[7];
  s8 = -x[3]*y[4]*z[7]+x[4]*y[3]*z[7]-z[6]*x[7]*y[4]-z[1]*x[6]*y[5]+x[6]*y
       [7]*z[3]-x[1]*y[6]*z[5]-y[1]*x[5]*z[6]+z[5]*x[4]*y[7]-z[5]*x[4]*y[0]+x[1]*y[5]*
       z[6]-y[6]*x[5]*z[7]-y[2]*x[3]*z[1]+z[1]*x[5]*y[6]-y[5]*x[1]*z[4]+z[6]*x[4]*y[7]
       +x[5]*y[1]*z[4]-x[5]*y[6]*z[4]+y[6]*x[3]*z[7]-x[5]*y[4]*z[1];
  s5 = s8+x[5]*y[4]*z[6]+z[5]*x[1]*y[4]+y[1]*x[6]*z[5]-z[6]*x[3]*y[7]+z[6]*
       x[7]*y[3]-z[5]*x[6]*y[4]-z[5]*x[4]*y[1]+z[5]*x[4]*y[6]+x[1]*y[6]*z[2]+x[2]*y[6]
       *z[3]+z[2]*x[6]*y[3]+z[1]*x[6]*y[2]+z[2]*x[3]*y[1]-z[2]*x[1]*y[3]-z[2]*x[3]*y
       [6]+y[2]*x[1]*z[3]+y[1]*x[2]*z[6]-z[0]*x[7]*y[3]+s7;
  s4 = 1/s5;
  s2 = s3*s4;
  const double unknown1 = s1*s2;
  s1 = 1.0/6.0;
  s8 = -z[2]*x[1]*y[2]*z[5]+z[2]*y[1]*x[2]*z[5]-z[2]*z[1]*x[2]*y[5]+z[2]*z
       [1]*x[5]*y[2]+2.0*y[5]*x[7]*z[4]*z[4]-y[1]*x[2]*z[0]*z[0]+x[0]*y[3]*z[7]*z[7]
       -2.0*z[5]*z[5]*x[4]*y[1]+2.0*z[5]*z[5]*x[1]*y[4]+z[5]*z[5]*x[0]*y[4]-2.0*z[2]*z
       [2]*x[1]*y[3]+2.0*z[2]*z[2]*x[3]*y[1]-x[0]*y[4]*z[7]*z[7]-y[0]*x[3]*z[7]*z[7]+x
       [1]*y[0]*z[5]*z[5];
  s7 = s8-y[1]*x[0]*z[5]*z[5]+z[1]*y[1]*x[2]*z[6]+y[1]*x[0]*z[2]*z[2]+z[2]*
       z[2]*x[3]*y[0]-z[2]*z[2]*x[0]*y[3]-x[1]*y[0]*z[2]*z[2]+2.0*z[5]*z[5]*x[4]*y[6]
       -2.0*z[5]*z[5]*x[6]*y[4]-z[5]*z[5]*x[7]*y[4]-x[6]*y[7]*z[5]*z[5]+2.0*z[2]*y[1]*
       x[2]*z[6]-2.0*z[2]*x[1]*y[2]*z[6]+2.0*z[2]*z[1]*x[6]*y[2]-y[6]*x[5]*z[7]*z[7]+
       2.0*x[6]*y[4]*z[7]*z[7];
  s8 = -2.0*y[6]*x[4]*z[7]*z[7]+x[6]*y[5]*z[7]*z[7]-2.0*z[2]*z[1]*x[2]*y[6]
       +z[4]*y[6]*x[7]*z[5]+x[5]*y[4]*z[6]*z[6]+z[6]*z[6]*x[4]*y[7]-z[6]*z[6]*x[7]*y
       [4]-2.0*z[6]*z[6]*x[7]*y[5]+2.0*z[6]*z[6]*x[5]*y[7]-y[5]*x[4]*z[6]*z[6]+2.0*z
       [0]*z[0]*x[3]*y[4]-x[6]*y[5]*z[2]*z[2]+z[1]*z[1]*x[5]*y[6]-z[1]*z[1]*x[6]*y[5]-
       z[5]*z[5]*x[4]*y[0];
  s6 = s8+2.0*x[1]*y[3]*z[0]*z[0]+2.0*x[1]*y[6]*z[2]*z[2]-2.0*y[1]*x[6]*z
       [2]*z[2]-y[1]*x[5]*z[2]*z[2]-z[1]*z[1]*x[2]*y[6]-2.0*z[1]*z[1]*x[2]*y[5]+2.0*z
       [1]*z[1]*x[5]*y[2]+z[1]*y[1]*x[6]*z[5]+y[1]*x[2]*z[5]*z[5]+z[2]*z[1]*x[2]*y[0]+
       z[1]*x[1]*y[5]*z[6]-z[1]*x[1]*y[6]*z[5]-z[1]*y[1]*x[5]*z[6]-z[1]*x[2]*y[6]*z[5]
       +z[1]*x[6]*y[2]*z[5]+s7;
  s8 = -x[1]*y[2]*z[5]*z[5]+z[1]*x[5]*y[6]*z[2]-2.0*z[2]*z[2]*x[3]*y[6]+2.0
       *z[2]*z[2]*x[6]*y[3]+z[2]*z[2]*x[7]*y[3]-z[2]*z[2]*x[3]*y[7]-z[1]*x[6]*y[5]*z
       [2]+2.0*z[1]*x[1]*y[5]*z[2]-2.0*x[3]*y[4]*z[7]*z[7]+2.0*x[4]*y[3]*z[7]*z[7]+x
       [5]*y[6]*z[2]*z[2]+y[1]*x[2]*z[6]*z[6]+y[0]*x[4]*z[7]*z[7]+z[2]*x[2]*y[3]*z[0]-
       x[1]*y[2]*z[6]*z[6];
  s7 = s8-z[7]*z[2]*x[3]*y[7]+x[2]*y[6]*z[3]*z[3]-y[2]*x[6]*z[3]*z[3]-z[6]*
       x[2]*y[3]*z[7]-z[2]*z[1]*x[0]*y[2]+z[6]*z[2]*x[6]*y[3]-z[6]*z[2]*x[3]*y[6]+z[6]
       *x[2]*y[6]*z[3]+z[2]*x[1]*y[2]*z[0]+z[6]*y[2]*x[3]*z[7]-z[4]*z[5]*x[6]*y[4]+z
       [4]*z[5]*x[4]*y[6]-z[4]*y[6]*x[5]*z[7]+z[4]*z[6]*x[4]*y[7]+z[4]*x[5]*y[4]*z[6];
  s8 = -z[6]*y[2]*x[6]*z[3]-z[4]*y[5]*x[4]*z[6]-z[2]*y[1]*x[5]*z[6]+z[2]*x
       [1]*y[5]*z[6]+z[4]*x[6]*y[4]*z[7]+2.0*z[4]*z[5]*x[4]*y[7]-z[4]*z[6]*x[7]*y[4]+x
       [6]*y[7]*z[3]*z[3]-2.0*z[4]*z[5]*x[7]*y[4]-2.0*z[4]*y[5]*x[4]*z[7]-z[4]*y[6]*x
       [4]*z[7]+z[4]*x[6]*y[5]*z[7]-z[4]*x[6]*y[7]*z[5]+2.0*z[4]*x[5]*y[4]*z[7]+z[2]*x
       [2]*y[5]*z[6]-z[2]*x[2]*y[6]*z[5];
  s5 = s8+z[2]*x[6]*y[2]*z[5]-z[2]*x[5]*y[2]*z[6]-z[2]*x[2]*y[3]*z[7]-x[2]*
       y[3]*z[7]*z[7]+2.0*z[2]*x[2]*y[3]*z[1]-z[2]*y[2]*x[3]*z[0]+z[2]*y[2]*x[0]*z[3]-
       z[2]*x[2]*y[0]*z[3]-z[7]*y[2]*x[7]*z[3]+z[7]*z[2]*x[7]*y[3]+z[7]*x[2]*y[7]*z[3]
       +z[6]*y[1]*x[2]*z[5]-z[6]*x[1]*y[2]*z[5]+z[5]*x[1]*y[5]*z[2]+s6+s7;
  s8 = z[5]*z[1]*x[5]*y[2]-z[5]*z[1]*x[2]*y[5]-y[6]*x[7]*z[2]*z[2]+2.0*z[2]
       *x[2]*y[6]*z[3]-2.0*z[2]*x[2]*y[3]*z[6]+2.0*z[2]*y[2]*x[3]*z[6]+y[2]*x[3]*z[6]*
       z[6]+y[6]*x[7]*z[5]*z[5]+z[2]*y[2]*x[3]*z[7]-z[2]*y[2]*x[7]*z[3]-2.0*z[2]*y[2]*
       x[6]*z[3]+z[2]*x[2]*y[7]*z[3]+x[6]*y[2]*z[5]*z[5]-2.0*z[2]*x[2]*y[1]*z[3]-x[2]*
       y[6]*z[5]*z[5];
  s7 = s8-y[1]*x[5]*z[6]*z[6]+z[6]*x[1]*y[6]*z[2]-z[3]*z[2]*x[3]*y[6]+z[6]*
       z[1]*x[6]*y[2]-z[6]*z[1]*x[2]*y[6]-z[6]*y[1]*x[6]*z[2]-2.0*x[5]*y[2]*z[6]*z[6]+
       z[4]*z[1]*x[0]*y[4]-z[3]*x[2]*y[3]*z[6]-z[5]*y[1]*x[5]*z[2]+z[3]*y[2]*x[3]*z[6]
       +2.0*x[2]*y[5]*z[6]*z[6]-z[5]*x[1]*y[5]*z[0]+y[2]*x[3]*z[7]*z[7]-x[2]*y[3]*z[6]
       *z[6];
  s8 = z[5]*y[5]*x[4]*z[0]+z[3]*z[2]*x[6]*y[3]+x[1]*y[5]*z[6]*z[6]+z[5]*y
       [5]*x[7]*z[4]-z[1]*x[1]*y[2]*z[6]+z[1]*x[1]*y[6]*z[2]+2.0*z[6]*y[6]*x[7]*z[5]-z
       [7]*y[6]*x[7]*z[2]-z[3]*y[6]*x[7]*z[2]+x[6]*y[7]*z[2]*z[2]-2.0*z[6]*y[6]*x[7]*z
       [2]-2.0*x[6]*y[3]*z[7]*z[7]-x[6]*y[2]*z[7]*z[7]-z[5]*x[6]*y[5]*z[2]+y[6]*x[2]*z
       [7]*z[7];
  s6 = s8+2.0*y[6]*x[3]*z[7]*z[7]+z[6]*z[6]*x[7]*y[3]-y[6]*x[7]*z[3]*z[3]+z
       [5]*x[5]*y[0]*z[4]+2.0*z[6]*z[6]*x[7]*y[2]-2.0*z[6]*z[6]*x[2]*y[7]-z[6]*z[6]*x
       [3]*y[7]+z[7]*y[6]*x[7]*z[5]+z[7]*y[5]*x[7]*z[4]-2.0*z[7]*x[7]*y[3]*z[4]+2.0*z
       [7]*x[3]*y[7]*z[4]-2.0*z[7]*x[4]*y[7]*z[3]+2.0*z[7]*x[7]*y[4]*z[3]-z[7]*y[0]*x
       [7]*z[4]-2.0*z[7]*z[6]*x[3]*y[7]+s7;
  s8 = s6+2.0*z[7]*z[6]*x[7]*y[3]+2.0*z[7]*x[6]*y[7]*z[3]+z[7]*x[6]*y[7]*z
       [2]-2.0*z[7]*y[6]*x[7]*z[3]+z[7]*z[6]*x[7]*y[2]-z[7]*z[6]*x[2]*y[7]+z[5]*y[1]*x
       [5]*z[0]-z[5]*z[1]*x[5]*y[0]+2.0*y[1]*x[6]*z[5]*z[5]-2.0*x[1]*y[6]*z[5]*z[5]+z
       [5]*z[1]*x[0]*y[5]+z[6]*y[6]*x[3]*z[7]+2.0*z[6]*x[6]*y[7]*z[2]-z[6]*y[6]*x[7]*z
       [3];
  s7 = s8+2.0*z[6]*y[6]*x[2]*z[7]-z[6]*x[6]*y[3]*z[7]+z[6]*x[6]*y[7]*z[3]
       -2.0*z[6]*x[6]*y[2]*z[7]-2.0*z[1]*y[1]*x[5]*z[2]-z[1]*y[1]*x[6]*z[2]-z[7]*z[0]*
       x[7]*y[3]-2.0*z[6]*x[6]*y[5]*z[2]-z[2]*z[6]*x[3]*y[7]+z[2]*x[6]*y[7]*z[3]-z[2]*
       z[6]*x[2]*y[7]+y[5]*x[6]*z[4]*z[4]+z[2]*y[6]*x[2]*z[7]+y[6]*x[7]*z[4]*z[4]+z[2]
       *z[6]*x[7]*y[2]-2.0*x[5]*y[7]*z[4]*z[4];
  s8 = -x[6]*y[7]*z[4]*z[4]-z[5]*y[5]*x[0]*z[4]-z[2]*x[6]*y[2]*z[7]-x[5]*y
       [6]*z[4]*z[4]-2.0*z[5]*y[1]*x[5]*z[6]+2.0*z[5]*z[1]*x[5]*y[6]+2.0*z[5]*x[1]*y
       [5]*z[6]-2.0*z[5]*z[1]*x[6]*y[5]-z[5]*x[5]*y[2]*z[6]+z[5]*x[5]*y[6]*z[2]+z[5]*x
       [2]*y[5]*z[6]+z[5]*z[5]*x[4]*y[7]-y[5]*x[4]*z[7]*z[7]+x[5]*y[4]*z[7]*z[7]+z[6]*
       z[1]*x[5]*y[6]+z[6]*y[1]*x[6]*z[5];
  s4 = s8-z[6]*z[1]*x[6]*y[5]-z[6]*x[1]*y[6]*z[5]+z[2]*z[6]*x[7]*y[3]+2.0*z
       [6]*x[6]*y[2]*z[5]+2.0*z[6]*x[5]*y[6]*z[2]-2.0*z[6]*x[2]*y[6]*z[5]+z[7]*z[0]*x
       [3]*y[7]+z[7]*z[0]*x[7]*y[4]+z[3]*z[6]*x[7]*y[3]-z[3]*z[6]*x[3]*y[7]-z[3]*x[6]*
       y[3]*z[7]+z[3]*y[6]*x[2]*z[7]-z[3]*x[6]*y[2]*z[7]+z[5]*x[5]*y[4]*z[7]+s5+s7;
  s8 = s4+z[3]*y[6]*x[3]*z[7]-z[7]*x[0]*y[7]*z[3]+z[6]*x[5]*y[4]*z[7]+z[7]*
       y[0]*x[7]*z[3]+z[5]*z[6]*x[4]*y[7]-2.0*z[5]*x[5]*y[6]*z[4]+2.0*z[5]*x[5]*y[4]*z
       [6]-z[5]*x[5]*y[7]*z[4]-z[5]*y[6]*x[5]*z[7]-z[5]*z[6]*x[7]*y[4]-z[7]*z[0]*x[4]*
       y[7]-z[5]*z[6]*x[7]*y[5]-z[5]*y[5]*x[4]*z[7]+z[7]*x[0]*y[7]*z[4];
  s7 = s8-2.0*z[5]*y[5]*x[4]*z[6]+z[5]*z[6]*x[5]*y[7]+z[5]*x[6]*y[5]*z[7]+
       2.0*z[5]*y[5]*x[6]*z[4]+z[6]*z[5]*x[4]*y[6]-z[6]*x[5]*y[6]*z[4]-z[6]*z[5]*x[6]*
       y[4]-z[6]*x[6]*y[7]*z[4]-2.0*z[6]*y[6]*x[5]*z[7]+z[6]*x[6]*y[4]*z[7]-z[6]*y[5]*
       x[4]*z[7]-z[6]*y[6]*x[4]*z[7]+z[6]*y[6]*x[7]*z[4]+z[6]*y[5]*x[6]*z[4]+2.0*z[6]*
       x[6]*y[5]*z[7];
  s8 = -2.0*z[6]*x[6]*y[7]*z[5]-z[2]*y[1]*x[2]*z[0]+2.0*z[7]*z[6]*x[4]*y[7]
       -2.0*z[7]*x[6]*y[7]*z[4]-2.0*z[7]*z[6]*x[7]*y[4]+z[7]*z[5]*x[4]*y[7]-z[7]*z[5]*
       x[7]*y[4]-z[7]*x[5]*y[7]*z[4]+2.0*z[7]*y[6]*x[7]*z[4]-z[7]*z[6]*x[7]*y[5]+z[7]*
       z[6]*x[5]*y[7]-z[7]*x[6]*y[7]*z[5]+z[1]*z[1]*x[6]*y[2]+s7+x[1]*y[5]*z[2]*z[2];
  s6 = s8+2.0*z[2]*y[2]*x[1]*z[3]-2.0*z[2]*y[2]*x[3]*z[1]-2.0*x[1]*y[4]*z
       [0]*z[0]+2.0*y[1]*x[4]*z[0]*z[0]+2.0*x[2]*y[7]*z[3]*z[3]-2.0*y[2]*x[7]*z[3]*z
       [3]-x[1]*y[5]*z[0]*z[0]+z[0]*z[0]*x[7]*y[4]+z[0]*z[0]*x[3]*y[7]+x[2]*y[3]*z[0]*
       z[0]-2.0*y[1]*x[3]*z[0]*z[0]+y[5]*x[4]*z[0]*z[0]-2.0*z[0]*z[0]*x[4]*y[3]+x[1]*y
       [2]*z[0]*z[0]-z[0]*z[0]*x[4]*y[7]+y[1]*x[5]*z[0]*z[0];
  s8 = s6-y[2]*x[3]*z[0]*z[0]+y[1]*x[0]*z[3]*z[3]-2.0*x[0]*y[7]*z[3]*z[3]-x
       [0]*y[4]*z[3]*z[3]-2.0*x[2]*y[0]*z[3]*z[3]-x[1]*y[0]*z[3]*z[3]+y[0]*x[4]*z[3]*z
       [3]-2.0*z[0]*y[1]*x[0]*z[4]+2.0*z[0]*z[1]*x[0]*y[4]+2.0*z[0]*x[1]*y[0]*z[4]-2.0
       *z[0]*z[1]*x[4]*y[0]-2.0*z[3]*x[2]*y[3]*z[7]-2.0*z[3]*z[2]*x[3]*y[7]+2.0*z[3]*z
       [2]*x[7]*y[3];
  s7 = s8+2.0*z[3]*y[2]*x[3]*z[7]+2.0*z[5]*y[5]*x[4]*z[1]+2.0*z[0]*y[1]*x
       [0]*z[3]-z[0]*y[0]*x[3]*z[7]-2.0*z[0]*y[0]*x[3]*z[4]-z[0]*x[1]*y[0]*z[2]+z[0]*z
       [1]*x[2]*y[0]-z[0]*y[1]*x[0]*z[5]-z[0]*z[1]*x[0]*y[2]-z[0]*x[0]*y[7]*z[3]-2.0*z
       [0]*z[1]*x[0]*y[3]-z[5]*x[5]*y[4]*z[0]-2.0*z[0]*x[0]*y[4]*z[3]+z[0]*x[0]*y[7]*z
       [4]-z[0]*z[2]*x[0]*y[3];
  s8 = s7+z[0]*x[5]*y[0]*z[4]+z[0]*z[1]*x[0]*y[5]-z[0]*x[2]*y[0]*z[3]-z[0]*
       z[1]*x[5]*y[0]-2.0*z[0]*x[1]*y[0]*z[3]+2.0*z[0]*y[0]*x[4]*z[3]-z[0]*x[0]*y[4]*z
       [7]+z[0]*x[1]*y[0]*z[5]+z[0]*y[0]*x[7]*z[3]+z[0]*y[2]*x[0]*z[3]-z[0]*y[5]*x[0]*
       z[4]+z[0]*z[2]*x[3]*y[0]+z[0]*x[2]*y[3]*z[1]+z[0]*x[0]*y[3]*z[7]-z[0]*x[2]*y[1]
       *z[3];
  s5 = s8+z[0]*y[1]*x[0]*z[2]+z[3]*x[1]*y[3]*z[0]-2.0*z[3]*y[0]*x[3]*z[7]-z
       [3]*y[0]*x[3]*z[4]-z[3]*x[1]*y[0]*z[2]+z[3]*z[0]*x[7]*y[4]+2.0*z[3]*z[0]*x[3]*y
       [7]+2.0*z[3]*x[2]*y[3]*z[0]-z[3]*y[1]*x[3]*z[0]-z[3]*z[1]*x[0]*y[3]-z[3]*z[0]*x
       [4]*y[3]+z[3]*x[1]*y[2]*z[0]-z[3]*z[0]*x[4]*y[7]-2.0*z[3]*z[2]*x[0]*y[3]-z[3]*x
       [0]*y[4]*z[7]-2.0*z[3]*y[2]*x[3]*z[0];
  s8 = s5+2.0*z[3]*z[2]*x[3]*y[0]+z[3]*x[2]*y[3]*z[1]+2.0*z[3]*x[0]*y[3]*z
       [7]+z[3]*y[1]*x[0]*z[2]-z[4]*y[0]*x[3]*z[7]-z[4]*x[1]*y[5]*z[0]-z[4]*y[1]*x[0]*
       z[5]+2.0*z[4]*z[0]*x[7]*y[4]+z[4]*z[0]*x[3]*y[7]+2.0*z[4]*y[5]*x[4]*z[0]+2.0*y
       [0]*x[7]*z[3]*z[3]+2.0*y[2]*x[0]*z[3]*z[3]-x[2]*y[1]*z[3]*z[3]-y[0]*x[3]*z[4]*z
       [4];
  s7 = s8-y[1]*x[0]*z[4]*z[4]+x[1]*y[0]*z[4]*z[4]+2.0*x[0]*y[7]*z[4]*z[4]+
       2.0*x[5]*y[0]*z[4]*z[4]-2.0*y[5]*x[0]*z[4]*z[4]+2.0*z[1]*z[1]*x[2]*y[0]-2.0*z
       [1]*z[1]*x[0]*y[2]+z[1]*z[1]*x[0]*y[4]-z[1]*z[1]*x[0]*y[3]-z[1]*z[1]*x[4]*y[0]+
       2.0*z[1]*z[1]*x[0]*y[5]-2.0*z[1]*z[1]*x[5]*y[0]+x[2]*y[3]*z[1]*z[1]-x[5]*y[4]*z
       [0]*z[0]-z[0]*z[0]*x[7]*y[3];
  s8 = s7+x[7]*y[4]*z[3]*z[3]-x[4]*y[7]*z[3]*z[3]+y[2]*x[1]*z[3]*z[3]+x[0]*
       y[3]*z[4]*z[4]-2.0*y[0]*x[7]*z[4]*z[4]+x[3]*y[7]*z[4]*z[4]-x[7]*y[3]*z[4]*z[4]-
       y[5]*x[1]*z[4]*z[4]+x[5]*y[1]*z[4]*z[4]+z[1]*z[1]*x[3]*y[0]+y[5]*x[4]*z[1]*z[1]
       -y[2]*x[3]*z[1]*z[1]-x[5]*y[4]*z[1]*z[1]-z[4]*x[0]*y[4]*z[3]-z[4]*z[0]*x[4]*y
       [3];
  s6 = s8-z[4]*z[1]*x[4]*y[0]-2.0*z[4]*z[0]*x[4]*y[7]+z[4]*y[1]*x[5]*z[0]
       -2.0*z[5]*x[5]*y[4]*z[1]-z[4]*x[1]*y[4]*z[0]+z[4]*y[0]*x[4]*z[3]-2.0*z[4]*x[0]*
       y[4]*z[7]+z[4]*x[1]*y[0]*z[5]-2.0*z[1]*x[1]*y[2]*z[5]+z[4]*x[0]*y[3]*z[7]+2.0*z
       [5]*x[5]*y[1]*z[4]+z[4]*y[1]*x[4]*z[0]+z[1]*y[1]*x[0]*z[3]+z[1]*x[1]*y[3]*z[0]
       -2.0*z[1]*x[1]*y[5]*z[0]-2.0*z[1]*x[1]*y[0]*z[2];
  s8 = s6-2.0*z[1]*y[1]*x[0]*z[5]-z[1]*y[1]*x[0]*z[4]+2.0*z[1]*y[1]*x[2]*z
       [5]-z[1]*y[1]*x[3]*z[0]-2.0*z[5]*y[5]*x[1]*z[4]+z[1]*y[5]*x[4]*z[0]+z[1]*x[1]*y
       [0]*z[4]+2.0*z[1]*x[1]*y[2]*z[0]-z[1]*z[2]*x[0]*y[3]+2.0*z[1]*y[1]*x[5]*z[0]-z
       [1]*x[1]*y[0]*z[3]-z[1]*x[1]*y[4]*z[0]+2.0*z[1]*x[1]*y[0]*z[5]-z[1]*y[2]*x[3]*z
       [0];
  s7 = s8+z[1]*z[2]*x[3]*y[0]-z[1]*x[2]*y[1]*z[3]+z[1]*y[1]*x[4]*z[0]+2.0*z
       [1]*y[1]*x[0]*z[2]+2.0*z[0]*z[1]*x[3]*y[0]+2.0*z[0]*x[0]*y[3]*z[4]+z[0]*z[5]*x
       [0]*y[4]+z[0]*y[0]*x[4]*z[7]-z[0]*y[0]*x[7]*z[4]-z[0]*x[7]*y[3]*z[4]-z[0]*z[5]*
       x[4]*y[0]-z[0]*x[5]*y[4]*z[1]+z[3]*z[1]*x[3]*y[0]+z[3]*x[0]*y[3]*z[4]+z[3]*z[0]
       *x[3]*y[4]+z[3]*y[0]*x[4]*z[7];
  s8 = s7+z[3]*x[3]*y[7]*z[4]-z[3]*x[7]*y[3]*z[4]-z[3]*x[3]*y[4]*z[7]+z[3]*
       x[4]*y[3]*z[7]-z[3]*y[2]*x[3]*z[1]+z[3]*z[2]*x[3]*y[1]-z[3]*z[2]*x[1]*y[3]-2.0*
       z[3]*z[0]*x[7]*y[3]+z[4]*z[0]*x[3]*y[4]+2.0*z[4]*z[5]*x[0]*y[4]+2.0*z[4]*y[0]*x
       [4]*z[7]-2.0*z[4]*x[5]*y[4]*z[0]+z[4]*y[5]*x[4]*z[1]+z[4]*x[7]*y[4]*z[3]-z[4]*x
       [4]*y[7]*z[3];
  s3 = s8-z[4]*x[3]*y[4]*z[7]+z[4]*x[4]*y[3]*z[7]-2.0*z[4]*z[5]*x[4]*y[0]-z
       [4]*x[5]*y[4]*z[1]+z[4]*z[5]*x[1]*y[4]-z[4]*z[5]*x[4]*y[1]-2.0*z[1]*y[1]*x[2]*z
       [0]+z[1]*z[5]*x[0]*y[4]-z[1]*z[5]*x[4]*y[0]-z[1]*y[5]*x[1]*z[4]+z[1]*x[5]*y[1]*
       z[4]+z[1]*z[5]*x[1]*y[4]-z[1]*z[5]*x[4]*y[1]+z[1]*z[2]*x[3]*y[1]-z[1]*z[2]*x[1]
       *y[3]+z[1]*y[2]*x[1]*z[3];
  s8 = y[1]*x[0]*z[3]+x[1]*y[3]*z[0]-y[0]*x[3]*z[7]-x[1]*y[5]*z[0]-y[0]*x
       [3]*z[4]-x[1]*y[0]*z[2]+z[1]*x[2]*y[0]-y[1]*x[0]*z[5]-z[1]*x[0]*y[2]-y[1]*x[0]*
       z[4]+z[1]*x[5]*y[2]+z[0]*x[7]*y[4]+z[0]*x[3]*y[7]+z[1]*x[0]*y[4]-x[1]*y[2]*z[5]
       +x[2]*y[3]*z[0]+y[1]*x[2]*z[5]-x[2]*y[3]*z[7];
  s7 = s8-z[1]*x[2]*y[5]-y[1]*x[3]*z[0]-x[0]*y[7]*z[3]-z[1]*x[0]*y[3]+y[5]*
       x[4]*z[0]-x[0]*y[4]*z[3]+y[5]*x[7]*z[4]-z[0]*x[4]*y[3]+x[1]*y[0]*z[4]-z[2]*x[3]
       *y[7]-y[6]*x[7]*z[2]+x[1]*y[5]*z[2]+y[6]*x[7]*z[5]+x[0]*y[7]*z[4]+x[1]*y[2]*z
       [0]-z[1]*x[4]*y[0]-z[0]*x[4]*y[7]-z[2]*x[0]*y[3];
  s8 = x[5]*y[0]*z[4]+z[1]*x[0]*y[5]-x[2]*y[0]*z[3]-z[1]*x[5]*y[0]+y[1]*x
       [5]*z[0]-x[1]*y[0]*z[3]-x[1]*y[4]*z[0]-y[1]*x[5]*z[2]+x[2]*y[7]*z[3]+y[0]*x[4]*
       z[3]-x[0]*y[4]*z[7]+x[1]*y[0]*z[5]-y[1]*x[6]*z[2]-y[2]*x[6]*z[3]+y[0]*x[7]*z[3]
       -y[2]*x[7]*z[3]+z[2]*x[7]*y[3]+y[2]*x[0]*z[3];
  s6 = s8+y[2]*x[3]*z[7]-y[2]*x[3]*z[0]-x[6]*y[5]*z[2]-y[5]*x[0]*z[4]+z[2]*
       x[3]*y[0]+x[2]*y[3]*z[1]+x[0]*y[3]*z[7]-x[2]*y[1]*z[3]+y[1]*x[4]*z[0]+y[1]*x[0]
       *z[2]-z[1]*x[2]*y[6]+y[2]*x[3]*z[6]-y[1]*x[2]*z[0]+z[1]*x[3]*y[0]-x[1]*y[2]*z
       [6]-x[2]*y[3]*z[6]+x[0]*y[3]*z[4]+z[0]*x[3]*y[4]+s7;
  s8 = x[5]*y[4]*z[7]+s6+y[5]*x[6]*z[4]-y[5]*x[4]*z[6]+z[6]*x[5]*y[7]-x[6]*
       y[2]*z[7]-x[6]*y[7]*z[5]+x[5]*y[6]*z[2]+x[6]*y[5]*z[7]+x[6]*y[7]*z[2]+y[6]*x[7]
       *z[4]-y[6]*x[4]*z[7]-y[6]*x[7]*z[3]+z[6]*x[7]*y[2]+x[2]*y[5]*z[6]-x[2]*y[6]*z
       [5]+y[6]*x[2]*z[7]+x[6]*y[2]*z[5];
  s7 = s8-x[5]*y[2]*z[6]-z[6]*x[7]*y[5]-z[5]*x[7]*y[4]+z[5]*x[0]*y[4]-y[5]*
       x[4]*z[7]+y[0]*x[4]*z[7]-z[6]*x[2]*y[7]-x[5]*y[4]*z[0]-x[5]*y[7]*z[4]-y[0]*x[7]
       *z[4]+y[5]*x[4]*z[1]-x[6]*y[7]*z[4]+x[7]*y[4]*z[3]-x[4]*y[7]*z[3]+x[3]*y[7]*z
       [4]-x[7]*y[3]*z[4]-x[6]*y[3]*z[7]+x[6]*y[4]*z[7];
  s8 = -x[3]*y[4]*z[7]+x[4]*y[3]*z[7]-z[6]*x[7]*y[4]-z[1]*x[6]*y[5]+x[6]*y
       [7]*z[3]-x[1]*y[6]*z[5]-y[1]*x[5]*z[6]+z[5]*x[4]*y[7]-z[5]*x[4]*y[0]+x[1]*y[5]*
       z[6]-y[6]*x[5]*z[7]-y[2]*x[3]*z[1]+z[1]*x[5]*y[6]-y[5]*x[1]*z[4]+z[6]*x[4]*y[7]
       +x[5]*y[1]*z[4]-x[5]*y[6]*z[4]+y[6]*x[3]*z[7]-x[5]*y[4]*z[1];
  s5 = s8+x[5]*y[4]*z[6]+z[5]*x[1]*y[4]+y[1]*x[6]*z[5]-z[6]*x[3]*y[7]+z[6]*
       x[7]*y[3]-z[5]*x[6]*y[4]-z[5]*x[4]*y[1]+z[5]*x[4]*y[6]+x[1]*y[6]*z[2]+x[2]*y[6]
       *z[3]+z[2]*x[6]*y[3]+z[1]*x[6]*y[2]+z[2]*x[3]*y[1]-z[2]*x[1]*y[3]-z[2]*x[3]*y
       [6]+y[2]*x[1]*z[3]+y[1]*x[2]*z[6]-z[0]*x[7]*y[3]+s7;
  s4 = 1/s5;
  s2 = s3*s4;
  const double unknown2 = s1*s2;

  return Point<3> (unknown0, unknown1, unknown2);
}



template <>
double TriaObjectAccessor<3, 3>::measure () const
{
  unsigned int vertex_indices[GeometryInfo<3>::vertices_per_cell];
  for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
    vertex_indices[i]=vertex_index(i);

  return GridTools::cell_measure<3>(this->tria->vertices, vertex_indices);
}

#endif



template <int dim>
unsigned int TriaObjectAccessor<3, dim>::number_of_children () const
{
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<n_children(); ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
}


#if deal_II_dimension == 3 

template <>
void
TriaObjectAccessor<3, 3>::
set_face_orientation (const unsigned int face,
                      const bool         orientation) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (face<GeometryInfo<3>::faces_per_cell,
          ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
  Assert (this->present_index * GeometryInfo<3>::faces_per_cell + face
          < this->tria->levels[this->present_level]
          ->cells.face_orientations.size(),
          ExcInternalError());
          
  this->tria->levels[this->present_level]
    ->cells.face_orientations[this->present_index *
                              GeometryInfo<3>::faces_per_cell
                              +
                              face]  = orientation;
}


template <>
void
TriaObjectAccessor<3, 3>::
set_face_flip (const unsigned int face,
	       const bool         flip) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (face<GeometryInfo<3>::faces_per_cell,
          ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
  Assert (this->present_index * GeometryInfo<3>::faces_per_cell + face
          < this->tria->levels[this->present_level]
          ->cells.face_flips.size(),
          ExcInternalError());
          
  this->tria->levels[this->present_level]
    ->cells.face_flips[this->present_index *
		       GeometryInfo<3>::faces_per_cell
		       +
		       face]  = flip;
}


template <>
void
TriaObjectAccessor<3, 3>::
set_face_rotation (const unsigned int face,
		   const bool         rotation) const
{
  Assert (used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (face<GeometryInfo<3>::faces_per_cell,
          ExcIndexRange (face, 0, GeometryInfo<3>::faces_per_cell));
  Assert (this->present_index * GeometryInfo<3>::faces_per_cell + face
          < this->tria->levels[this->present_level]
          ->cells.face_rotations.size(),
          ExcInternalError());
          
  this->tria->levels[this->present_level]
    ->cells.face_rotations[this->present_index *
			   GeometryInfo<3>::faces_per_cell
			   +
			   face]  = rotation;
}

#endif
// Remark: The following explicit instantiations needed to be moved to
// this place here to work around a problem with gcc3.3 on Apple MacOSX.
// The reason is that some of the functions instantiated here are used
// further down; if they are not explicitly instantiated here, then the
// compiler will do an implicit instantiation and give it internal linkage
// (despite the later explicit instantiation that should make sure it
// gets external linkage). To make sure the functions have external
// linkage, we need to place the explicit instantiation before the first
// use.
//
// For more information, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331

template class TriaObjectAccessor<1, deal_II_dimension>;

#if deal_II_dimension >= 2
template class TriaObjectAccessor<2, deal_II_dimension>;
#endif

#if deal_II_dimension >= 3
template class TriaObjectAccessor<3, deal_II_dimension>;
#endif



/*------------------------ Functions: CellAccessor<1> -----------------------*/


#if deal_II_dimension == 1

template <>
bool CellAccessor<1>::at_boundary () const
{
  return at_boundary(0) || at_boundary(1);
}



template <>
bool CellAccessor<1>::point_inside (const Point<1> &p) const
{
  return (this->vertex(0)[0] <= p[0]) && (p[0] <= this->vertex(1)[0]);
}



template <>
std::pair<unsigned int, unsigned int>
CellAccessor<1>::neighbor_of_coarser_neighbor (const unsigned int) const
{
  Assert(false, ExcImpossibleInDim(1));
  return std::make_pair (numbers::invalid_unsigned_int,
			 numbers::invalid_unsigned_int);
}


#endif


/*------------------------ Functions: CellAccessor<2> -----------------------*/


#if deal_II_dimension == 2

template <>
bool CellAccessor<2>::at_boundary () const
{
  return at_boundary(0) || at_boundary(1) || at_boundary(2) || at_boundary(3);
}



template <>
bool CellAccessor<2>::point_inside (const Point<2> &p) const
{
				   // we check whether the point is
				   // inside the cell by making sure
				   // that it on the inner side of
				   // each line defined by the faces,
				   // i.e. for each of the four faces
				   // we take the line that connects
				   // the two vertices and subdivide
				   // the whole domain by that in two
				   // and check whether the point is
				   // on the `cell-side' (rather than
				   // the `out-side') of this line. if
				   // the point is on the `cell-side'
				   // for all four faces, it must be
				   // inside the cell.

				   // we want the faces in counter
				   // clockwise orientation
  static const int direction[4]={-1,1,1,-1};
  for (unsigned int f=0; f<4; ++f)
    {
				       // vector from the first vertex
				       // of the line to the point
      const Point<2> to_p = p-this->vertex(
	GeometryInfo<2>::face_to_cell_vertices(f,0));
				       // vector describing the line
      const Point<2> face = direction[f]*(
	this->vertex(GeometryInfo<2>::face_to_cell_vertices(f,1)) -
	this->vertex(GeometryInfo<2>::face_to_cell_vertices(f,0)));

				       // if we rotate the face vector
				       // by 90 degrees to the left
				       // (i.e. it points to the
				       // inside) and take the scalar
				       // product with the vector from
				       // the vertex to the point,
				       // then the point is in the
				       // `cell-side' if the scalar
				       // product is positive. if this
				       // is not the case, we can be
				       // sure that the point is
				       // outside
      if ((-face(1)*to_p(0)+face(0)*to_p(1))<0)
	return false;
    };

				   // if we arrived here, then the
				   // point is inside for all four
				   // faces, and thus inside
  return true;
}


#endif


/*------------------------ Functions: CellAccessor<3> -----------------------*/


#if deal_II_dimension == 3

template <>
bool CellAccessor<3>::at_boundary () const
{
  return (at_boundary(0) || at_boundary(1) ||
	  at_boundary(2) || at_boundary(3) ||
	  at_boundary(4) || at_boundary(5));
}



template <>
bool CellAccessor<3>::point_inside (const Point<3> &p) const
{
				   // original implementation by Joerg
				   // Weimar
  
                                   // we first eliminate points based
                                   // on the maximum and minumum of
                                   // the corner coordinates, then
                                   // transform to the unit cell, and
                                   // check there.
  const unsigned int dim = 3;
  Point<dim> maxp = this->vertex(0);
  Point<dim> minp = this->vertex(0);

  for (unsigned int v=1; v<GeometryInfo<dim>::vertices_per_cell; ++v)
    for (unsigned int d=0; d<dim; ++d)
      {
	maxp[d] = std::max (maxp[d],this->vertex(v)[d]);
	minp[d] = std::min (minp[d],this->vertex(v)[d]);
      }

				   // rule out points outside the
				   // bounding box of this cell
  for (unsigned int d=0; d<dim; d++)
    if ((p[d] < minp[d]) || (p[d] > maxp[d]))
      return false;

				   // now we need to check more
				   // carefully: transform to the
				   // unit cube
				   // and check there.
  const TriaRawIterator<dim, CellAccessor<dim> > cell_iterator (*this);
  return (GeometryInfo<dim>::is_inside_unit_cell (
    StaticMappingQ1<dim>::mapping.transform_real_to_unit_cell(cell_iterator, p)));
}


#endif


/*------------------------ Functions: CellAccessor<dim> -----------------------*/


template <int dim>
unsigned char CellAccessor<dim>::material_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->cells.material_id[this->present_index];
}



template <int dim>
void CellAccessor<dim>::set_material_id (const unsigned char mat_id) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->cells.material_id[this->present_index] = mat_id;
}



template <int dim>
void CellAccessor<dim>::recursively_set_material_id (const unsigned char mat_id) const
{
  set_material_id (mat_id);

  if (this->has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_set_material_id (mat_id);
}



template <int dim>
unsigned int CellAccessor<dim>::subdomain_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->subdomain_ids[this->present_index];
}



template <int dim>
void
CellAccessor<dim>::set_subdomain_id (const unsigned int new_subdomain_id) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->subdomain_ids[this->present_index]
    = new_subdomain_id;
}



template <int dim>
void CellAccessor<dim>::set_neighbor (const unsigned int i,
				      const TriaIterator<dim,CellAccessor<dim> > &pointer) const
{
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(i));

  if (pointer.state() == IteratorState::valid)
    {
      this->tria->levels[this->present_level]->
	neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].first
	= pointer->present_level;
      this->tria->levels[this->present_level]->
	neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].second
	= pointer->present_index;
    }
  else
    {
      this->tria->levels[this->present_level]->
	neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].first
	= -1;
      this->tria->levels[this->present_level]->
	neighbors[this->present_index*GeometryInfo<dim>::faces_per_cell+i].second
	= -1;
    };
}



template <int dim>
unsigned int CellAccessor<dim>::neighbor_of_neighbor_internal (const unsigned int neighbor) const
{
  Assert (neighbor < GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(neighbor));

  if (dim==1)
    return GeometryInfo<dim>::opposite_face[neighbor];
  
  const TriaIterator<dim,CellAccessor<dim> > neighbor_cell = this->neighbor(neighbor);
  
				   // usually, on regular patches of
				   // the grid, this cell is just on
				   // the opposite side of the
				   // neighbor that the neighbor is of
				   // this cell. for example in 2d, if
				   // we want to know the
				   // neighbor_of_neighbor if
				   // neighbor==1 (the right
				   // neighbor), then we will get 3
				   // (the left neighbor) in most
				   // cases. look up this relationship
				   // in the table provided by
				   // GeometryInfo and try it

  const unsigned int this_face_index=face_index(neighbor);
  
  const unsigned int neighbor_guess
    = GeometryInfo<dim>::opposite_face[neighbor];
  
  if (neighbor_cell->face_index (neighbor_guess) == this_face_index)
    return neighbor_guess;
  else
				     // if the guess was false, then
				     // we need to loop over all
				     // neighbors and find the number
				     // the hard way
    {
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	if (neighbor_cell->face_index (face_no) == this_face_index)
	  return face_no;

				       // running over all neighbors
				       // faces we did not find the
				       // present face. Thereby the
				       // neighbor must be coarser
				       // than the present
				       // cell. Return an invalid
				       // unsigned int in this case.
      return numbers::invalid_unsigned_int;
    }
}



template <int dim>
unsigned int CellAccessor<dim>::neighbor_of_neighbor (const unsigned int neighbor) const
{
  const unsigned int n2=neighbor_of_neighbor_internal(neighbor);
  Assert (n2!=numbers::invalid_unsigned_int, TriaAccessorExceptions::ExcNeighborIsCoarser());

  return n2;
}



template <int dim>
bool
CellAccessor<dim>::neighbor_is_coarser (const unsigned int neighbor) const
{
  return neighbor_of_neighbor_internal(neighbor)==numbers::invalid_unsigned_int;
}


# if deal_II_dimension == 2

template <int dim>
std::pair<unsigned int, unsigned int>
CellAccessor<dim>::neighbor_of_coarser_neighbor (const unsigned int neighbor) const
{
  Assert (neighbor < GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(neighbor));
				   // make sure that the neighbor is
				   // on a coarser level
  Assert (neighbor_is_coarser(neighbor),
	  TriaAccessorExceptions::ExcNeighborIsNotCoarser());
  
  const int this_face_index=face_index(neighbor);
  const TriaIterator<dim,CellAccessor<dim> > neighbor_cell = this->neighbor(neighbor);
  
				   // usually, on regular patches of
				   // the grid, this cell is just on
				   // the opposite side of the
				   // neighbor that the neighbor is of
				   // this cell. for example in 2d, if
				   // we want to know the
				   // neighbor_of_neighbor if
				   // neighbor==1 (the right
				   // neighbor), then we will get 0
				   // (the left neighbor) in most
				   // cases. look up this relationship
				   // in the table provided by
				   // GeometryInfo and try it
  const unsigned int face_no_guess
    = GeometryInfo<dim>::opposite_face[neighbor];

  const TriaIterator<dim,TriaObjectAccessor<dim-1, dim> > face_guess
    =neighbor_cell->face(face_no_guess);
  
  if (face_guess->has_children())
    for (unsigned int subface_no=0; subface_no<face_guess->n_children(); ++subface_no)
      if (face_guess->child_index(subface_no)==this_face_index)
	return std::make_pair (face_no_guess, subface_no);

				     // if the guess was false, then
				     // we need to loop over all faces
				     // and subfaces and find the
				     // number the hard way
  for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
      if (face_no!=face_no_guess)
	{
	  const TriaIterator<dim,TriaObjectAccessor<dim-1, dim> > face
	    =neighbor_cell->face(face_no);
	  if (face->has_children())
	    for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
	      if (face->child_index(subface_no)==this_face_index)
		return std::make_pair (face_no, subface_no);
	}
    }
  
				   // we should never get here,
				   // since then we did not find
				   // our way back...
  Assert (false, ExcInternalError());
  return std::make_pair (deal_II_numbers::invalid_unsigned_int,
			 deal_II_numbers::invalid_unsigned_int);
}

#endif

#if deal_II_dimension == 3

template <int dim>
std::pair<unsigned int, unsigned int>
CellAccessor<dim>::neighbor_of_coarser_neighbor (const unsigned int neighbor) const
{
  Assert (neighbor < GeometryInfo<dim>::faces_per_cell,
	  TriaAccessorExceptions::ExcInvalidNeighbor(neighbor));
				   // make sure that the neighbor is
				   // on a coarser level
  Assert (neighbor_is_coarser(neighbor),
	  TriaAccessorExceptions::ExcNeighborIsNotCoarser());

  const int this_face_index=face_index(neighbor);
  const TriaIterator<dim,CellAccessor<dim> > neighbor_cell = this->neighbor(neighbor);
  
				   // usually, on regular patches of
				   // the grid, this cell is just on
				   // the opposite side of the
				   // neighbor that the neighbor is of
				   // this cell. for example in 2d, if
				   // we want to know the
				   // neighbor_of_neighbor if
				   // neighbor==1 (the right
				   // neighbor), then we will get 0
				   // (the left neighbor) in most
				   // cases. look up this relationship
				   // in the table provided by
				   // GeometryInfo and try it
  const unsigned int face_no_guess
    = GeometryInfo<dim>::opposite_face[neighbor];

  const TriaIterator<dim,TriaObjectAccessor<dim-1, dim> > face_guess
    =neighbor_cell->face(face_no_guess);
  
  if (face_guess->has_children())
    for (unsigned int subface_no=0; subface_no<face_guess->n_children(); ++subface_no)
      if (face_guess->child_index(subface_no)==this_face_index)
					     // call a helper function, that
					     // translates the current subface
					     // number to a subface number for
					     // the current FaceRefineCase
	return std::make_pair (face_no_guess, translate_subface_no(face_guess, subface_no));
      else if (face_guess->child(subface_no)->has_children())
	for (unsigned int subsub_no=0; subsub_no<face_guess->child(subface_no)->n_children(); ++subsub_no)
	  if (face_guess->child(subface_no)->child_index(subsub_no)==this_face_index)
					     // call a helper function, that
					     // translates the current subface
					     // number and subsubface number to
					     // a subface number for the current
					     // FaceRefineCase
	    return std::make_pair (face_no_guess, translate_subface_no(face_guess, subface_no, subsub_no));
	  
  

				     // if the guess was false, then
				     // we need to loop over all faces
				     // and subfaces and find the
				     // number the hard way
  for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
      if (face_no!=face_no_guess)
	{
	  const TriaIterator<dim,TriaObjectAccessor<dim-1, dim> > face
	    =neighbor_cell->face(face_no);
	  if (face->has_children())
	    for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
	      if (face->child_index(subface_no)==this_face_index)
					     // call a helper function, that
					     // translates the current subface
					     // number to a subface number for
					     // the current FaceRefineCase
		return std::make_pair (face_no, translate_subface_no(face, subface_no));
	      else if (face->child(subface_no)->has_children())
		for (unsigned int subsub_no=0; subsub_no<face->child(subface_no)->n_children(); ++subsub_no)
		  if (face->child(subface_no)->child_index(subsub_no)==this_face_index)
					     // call a helper function, that
					     // translates the current subface
					     // number and subsubface number to
					     // a subface number for the current
					     // FaceRefineCase
		    return std::make_pair (face_no, translate_subface_no(face, subface_no, subsub_no));
	}
    }
  
				   // we should never get here,
				   // since then we did not find
				   // our way back...
  Assert (false, ExcInternalError());
  return std::make_pair (numbers::invalid_unsigned_int,
			 numbers::invalid_unsigned_int);
}

#endif



template <int dim>
bool CellAccessor<dim>::at_boundary (const unsigned int i) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (i,0,GeometryInfo<dim>::faces_per_cell));
  
  return (neighbor_index(i) == -1);
}



#if deal_II_dimension == 1

template <>
bool CellAccessor<1>::has_boundary_lines () const
{
  return at_boundary();
}

#else

template <int dim>
bool CellAccessor<dim>::has_boundary_lines () const
{
  for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
    if (this->line(l)->at_boundary())
      return true;
  
  return false;
}

#endif

#if deal_II_dimension == 1

template <>
TriaIterator<1,CellAccessor<1> >
CellAccessor<1>::
neighbor_child_on_subface (const unsigned int,
                           const unsigned int) const
{
  Assert (false, ExcNotImplemented());
  return TriaIterator<1,CellAccessor<1> >();
}

#endif

#if deal_II_dimension == 2

template <>
TriaIterator<2,CellAccessor<2> >
CellAccessor<2>::
neighbor_child_on_subface (const unsigned int face,
                           const unsigned int subface) const
{
  Assert (!this->has_children(),
          ExcMessage ("The present cell must not have children!"));
  Assert (!this->at_boundary(face),
          ExcMessage ("The present cell must have a valid neighbor!"));
  Assert (this->neighbor(face)->has_children() == true,
          ExcMessage ("The neighbor must have children!"));

  const unsigned int neighbor_neighbor
    = this->neighbor_of_neighbor (face);
  const unsigned int neighbor_child_index
    = GeometryInfo<2>::child_cell_on_face(
      this->neighbor(face)->refinement_case(),neighbor_neighbor,subface);

  TriaIterator<2,CellAccessor<2> > sub_neighbor = this->neighbor(face)->child(neighbor_child_index);
				   // the neighbors child can have children,
				   // which are not further refined along the
				   // face under consideration. as we are
				   // normally interested in one of this
				   // child's child, search for the right one.
  while(sub_neighbor->has_children())
    {
      Assert ((GeometryInfo<2>::face_refinement_case(sub_neighbor->refinement_case(),
						     neighbor_neighbor) ==
	       RefinementCase<2>::no_refinement),
	      ExcInternalError());
      sub_neighbor = sub_neighbor->child(GeometryInfo<2>::child_cell_on_face(
					   sub_neighbor->refinement_case(),neighbor_neighbor,0));
      
    }
  
  return sub_neighbor;
}

#endif

#if deal_II_dimension == 3

template <>
TriaIterator<3,CellAccessor<3> >
CellAccessor<3>::
neighbor_child_on_subface (const unsigned int face,
                           const unsigned int subface) const
{
  Assert (!this->has_children(),
          ExcMessage ("The present cell must not have children!"));
  Assert (!this->at_boundary(face),
          ExcMessage ("The present cell must have a valid neighbor!"));
  Assert (this->neighbor(face)->has_children() == true,
          ExcMessage ("The neighbor must have children!"));
				   // this function returns the neighbor's
				   // child on a given face and
				   // subface. 

				   // we have to consider one other aspect here:
				   // The face might be refined
				   // anisotropically. In this case, the subface
				   // number refers to the following, where we
				   // look at the face from the current cell,
				   // thus the subfaces are in standard
				   // orientation concerning the cell
				   //
				   // for isotropic refinement
				   //
				   // *---*---*
				   // | 2 | 3 |
				   // *---*---*
				   // | 0 | 1 |
				   // *---*---*
				   //
				   // for 2*anisotropic refinement
				   // (first cut_y, then cut_x)
				   //
				   // *---*---*
				   // | 2 | 3 |
				   // *---*---*
				   // | 0 | 1 |
				   // *---*---*
				   //
				   // for 2*anisotropic refinement
				   // (first cut_x, then cut_y)
				   //
				   // *---*---*
				   // | 1 | 3 |
				   // *---*---*
				   // | 0 | 2 |
				   // *---*---*
				   //
				   // for purely anisotropic refinement:
				   //
				   // *---*---*      *-------*
				   // |   |   |	     |   1   |
				   // | 0 | 1 |  or  *-------* 
				   // |   |   |	     |   0   |
				   // *---*---*	     *-------*
				   //
				   // for "mixed" refinement:
				   //
				   // *---*---*      *---*---*      *---*---*      *-------*
				   // |   | 2 |      | 1 |   |      | 1 | 2 |      |   2   |
				   // | 0 *---*  or  *---* 2 |  or  *---*---*  or  *---*---*
				   // |   | 1 |      | 0 |   |      |   0   |      | 0 | 1 |
				   // *---*---*      *---*---*      *-------*      *---*---*

  const Triangulation<3>::face_iterator mother_face=this->face(face);
  const unsigned int total_children=mother_face->number_of_children();
  Assert (subface<total_children,ExcIndexRange(subface,0,total_children));
  Assert(total_children<=GeometryInfo<3>::max_children_per_face, ExcInternalError());
  
  unsigned int neighbor_neighbor;
  TriaIterator<3,CellAccessor<3> > neighbor_child;
  const TriaIterator<3,CellAccessor<3> > neighbor=this->neighbor(face);

  
  const RefinementCase<2> mother_face_ref_case=mother_face->refinement_case();
  if (mother_face_ref_case==RefinementCase<2>::cut_xy) // total_children==4
    {
				       // this case is quite easy. we are sure,
				       // that the neighbor is not coarser.

				       // get the neighbor's number for the given
				       // face and the neighbor
      neighbor_neighbor
	= this->neighbor_of_neighbor (face);
      
				       // now use the info provided by GeometryInfo
				       // to extract the neighbors child number
      const unsigned int neighbor_child_index
	= GeometryInfo<3>::child_cell_on_face(neighbor->refinement_case(),
					      neighbor_neighbor, subface,
					      neighbor->face_orientation(neighbor_neighbor),
					      neighbor->face_flip(neighbor_neighbor),
					      neighbor->face_rotation(neighbor_neighbor));
      neighbor_child=
	neighbor->child(neighbor_child_index);
      
				       // make sure that the neighbor child cell we
				       // have found shares the desired subface.
      Assert((this->face(face)->child(subface) ==
	      neighbor_child->face(neighbor_neighbor)),
	     ExcInternalError());
    }
  else //-> the face is refined anisotropically
    {
					   // first of all, we have to find the
					   // neighbor at one of the anisotropic
					   // children of the
					   // mother_face. determine, which of
					   // these we need.
	  unsigned int first_child_to_find;
	  unsigned int neighbor_child_index;
	  if (total_children==2)
	    first_child_to_find=subface;
	  else 
	    {
	      first_child_to_find=subface/2;
	      if (total_children==3 &&
		  subface==1 &&
		  !mother_face->child(0)->has_children())
		first_child_to_find=1;
	    }
      if (neighbor_is_coarser(face))
	{
	  std::pair<unsigned int, unsigned int> indices=neighbor_of_coarser_neighbor(face);
	  neighbor_neighbor=indices.first;

	  
					   // we have to translate our
					   // subface_index according to the
					   // RefineCase and subface index of
					   // the coarser face (our face is an
					   // anisotropic child of the coarser
					   // face), 'a' denotes our
					   // subface_index 0 and 'b' denotes
					   // our subface_index 1, whereas 0...3
					   // denote isotropic subfaces of the
					   // coarser face
					   //
					   // cut_x and coarser_subface_index=0
					   //
					   // *---*---*
					   // |b=2|   |
					   // |   |   |
					   // |a=0|   |
					   // *---*---*
					   //
					   // cut_x and coarser_subface_index=1
					   //
					   // *---*---*
					   // |   |b=3|
					   // |   |   |
					   // |   |a=1|
					   // *---*---*
					   //
					   // cut_y and coarser_subface_index=0
					   //
					   // *-------*
					   // |       |
					   // *-------*
					   // |a=0 b=1|
					   // *-------*
					   //
					   // cut_y and coarser_subface_index=1
					   //
					   // *-------*
					   // |a=2 b=3|
					   // *-------*
					   // |       |
					   // *-------*
	  unsigned int iso_subface;
 	  if (neighbor->face(neighbor_neighbor)->refinement_case()==RefinementCase<2>::cut_x)
 	    iso_subface=2*first_child_to_find + indices.second;
 	  else
   	    {
  	      Assert(neighbor->face(neighbor_neighbor)->refinement_case()==RefinementCase<2>::cut_y,
  		     ExcInternalError());
  	      iso_subface=first_child_to_find + 2*indices.second;
  	    }
	  neighbor_child_index
	    =GeometryInfo<3>::child_cell_on_face(neighbor->refinement_case(),
						 neighbor_neighbor,
						 iso_subface,
						 neighbor->face_orientation(neighbor_neighbor),
						 neighbor->face_flip(neighbor_neighbor),
						 neighbor->face_rotation(neighbor_neighbor));
	}
      else //neighbor is not coarser
	{
	  neighbor_neighbor=neighbor_of_neighbor(face);
	  neighbor_child_index
	    =GeometryInfo<3>::child_cell_on_face(neighbor->refinement_case(),
						 neighbor_neighbor,
						 first_child_to_find,
						 neighbor->face_orientation(neighbor_neighbor),
						 neighbor->face_flip(neighbor_neighbor),
						 neighbor->face_rotation(neighbor_neighbor),
						 mother_face_ref_case);
	}
      
      neighbor_child=neighbor->child(neighbor_child_index);
				       // it might be, that the neighbor_child
				       // has children, which are not refined
				       // along the given subface. go down that
				       // list and deliver the last of those.
      while (neighbor_child->has_children() &&
	     GeometryInfo<3>::face_refinement_case(neighbor_child->refinement_case(),neighbor_neighbor)==RefinementCase<2>::no_refinement)
	neighbor_child=neighbor_child->child(GeometryInfo<3>::child_cell_on_face(neighbor_child->refinement_case(),
										 neighbor_neighbor,0));

				       // if there are two total subfaces, we
				       // are finished. if there are four we
				       // have to get a child of our current
				       // neighbor_child. If there are three,
				       // we have to check which of the two
				       // possibilities applies.
      if (total_children==3)
	{
	  if (mother_face->child(0)->has_children())
	    {
	      if (subface<2)
		neighbor_child=neighbor_child->child(GeometryInfo<3>::child_cell_on_face(neighbor_child->refinement_case(),
											 neighbor_neighbor,subface,
											 neighbor_child->face_orientation(neighbor_neighbor),
											 neighbor_child->face_flip(neighbor_neighbor),
											 neighbor_child->face_rotation(neighbor_neighbor),
											 mother_face->child(0)->refinement_case()));
	    }
	  else
	    {
	      Assert(mother_face->child(1)->has_children(), ExcInternalError());
	      if (subface>0)
		neighbor_child=neighbor_child->child(GeometryInfo<3>::child_cell_on_face(neighbor_child->refinement_case(),
											 neighbor_neighbor,subface-1,
											 neighbor_child->face_orientation(neighbor_neighbor),
											 neighbor_child->face_flip(neighbor_neighbor),
											 neighbor_child->face_rotation(neighbor_neighbor),
											 mother_face->child(1)->refinement_case()));
	    }
	}
      else if (total_children==4)
	{
	  neighbor_child=neighbor_child->child(GeometryInfo<3>::child_cell_on_face(neighbor_child->refinement_case(),
										   neighbor_neighbor,subface%2,
										   neighbor_child->face_orientation(neighbor_neighbor),
										   neighbor_child->face_flip(neighbor_neighbor),
										   neighbor_child->face_rotation(neighbor_neighbor),
										   mother_face->child(subface/2)->refinement_case()));
	}
    }

				   // it might be, that the neighbor_child has
				   // children, which are not refined along the
				   // given subface. go down that list and
				   // deliver the last of those.
  while (neighbor_child->has_children())
    neighbor_child=neighbor_child->child(GeometryInfo<3>::child_cell_on_face(neighbor_child->refinement_case(),
									     neighbor_neighbor,0));

#ifdef DEBUG
				   // check, whether the face neighbor_child
				   // matches the requested subface
  Triangulation<3>::face_iterator requested;
  switch (this->subface_case(face))
    {
      case internal::SubfaceCase<3>::case_x:
      case internal::SubfaceCase<3>::case_y:
      case internal::SubfaceCase<3>::case_xy:
	    requested=mother_face->child(subface);
	    break;
      case internal::SubfaceCase<3>::case_x1y2y:
      case internal::SubfaceCase<3>::case_y1x2x:
	    requested=mother_face->child(subface/2)->child(subface%2);
	    break;

      case internal::SubfaceCase<3>::case_x1y:
      case internal::SubfaceCase<3>::case_y1x:
	    switch (subface)
	      {
		case 0:
		case 1:
		      requested=mother_face->child(0)->child(subface);
		      break;
		case 2:
		      requested=mother_face->child(1);
		      break;
		default:
		      Assert(false, ExcInternalError());
	      }
	    break;
      case internal::SubfaceCase<3>::case_x2y:
      case internal::SubfaceCase<3>::case_y2x:
	    switch (subface)
	      {
		case 0:
		      requested=mother_face->child(0);
		      break;
		case 1:
		case 2:
		      requested=mother_face->child(1)->child(subface-1);
		      break;
		default:
		      Assert(false, ExcInternalError());
	      }
	    break;
      default:
	    Assert(false, ExcInternalError());
	    break;
    }
  Assert(requested==neighbor_child->face(neighbor_neighbor),
	 ExcInternalError());
#endif

  return neighbor_child;
  
}

#endif


// Remark: The explicit instantiations for "TriaObjectAccessor" were moved
// to the top of this source file. The reason is a slightly buggy version
// of the Apple gcc v.3.3.
// For more information, see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24331

// explicit instantiations
template class TriaAccessor<1,deal_II_dimension>;
#if deal_II_dimension >= 2
template class TriaAccessor<2,deal_II_dimension>;
#endif
#if deal_II_dimension >= 3
template class TriaAccessor<3,deal_II_dimension>;
#endif

template class CellAccessor<deal_II_dimension>;
template class TriaRawIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaRawIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;

#if deal_II_dimension >= 2
template class TriaRawIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
#endif

#if deal_II_dimension >= 3
template class TriaRawIterator<deal_II_dimension,TriaObjectAccessor<3, deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,TriaObjectAccessor<3, deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,TriaObjectAccessor<3, deal_II_dimension> >;
#endif

DEAL_II_NAMESPACE_CLOSE

