/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_accessor.templates.h>
#include <grid/tria_iterator.templates.h>
#include <grid/geometry_info.h>

#include <cmath>






/*------------------------ Functions: LineAccessor ---------------------------*/

template <int dim>
void TriaObjectAccessor<1, dim>::set (const Line &line) const {
  tria->levels[present_level]->lines.lines[present_index] = line;
};



template <int dim>
int TriaObjectAccessor<1, dim>::vertex_index (const unsigned int i) const {
  Assert (i<2, ExcIndexRange(i,0,2));
  return tria->levels[present_level]->lines.lines[present_index].vertex (i);
};



template <int dim>
Point<dim> &
TriaObjectAccessor<1, dim>::vertex (const unsigned int i) const {
  return tria->vertices[vertex_index(i)];
};



template <int dim>
void TriaObjectAccessor<1, dim>::set_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->lines.used[present_index] = true;
};
  


template <int dim>
void TriaObjectAccessor<1, dim>::clear_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->lines.used[present_index] = false;
};



template <int dim>
void TriaObjectAccessor<1, dim>::recursively_set_user_flag () const {
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_set_user_flag ();
};



template <int dim>
void TriaObjectAccessor<1, dim>::recursively_clear_user_flag () const {
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<2; ++c)
      child(c)->recursively_clear_user_flag ();
};



template <int dim>
void TriaObjectAccessor<1, dim>::set_user_pointer (void *p) const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_pointers[present_index] = p;
};



template <int dim>
void TriaObjectAccessor<1, dim>::clear_user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->lines.user_pointers[present_index] = 0;
};



template <int dim>
void * TriaObjectAccessor<1, dim>::user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->lines.user_pointers[present_index];
};
  

template <int dim>
void TriaObjectAccessor<1, dim>::set_children (const int index) const {
  Assert (used(), ExcCellNotUsed());
  Assert ((index==-1) ||
	  (!has_children() && (index>=0)), ExcCantSetChildren(index));
  
  tria->levels[present_level]->lines.children[present_index] = index;
};



template <int dim>
void TriaObjectAccessor<1, dim>::clear_children () const {
  set_children (-1);
};



template <int dim>
unsigned char TriaObjectAccessor<1, dim>::boundary_indicator () const {
  Assert (dim>=2, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  return tria->levels[present_level]->lines.material_id[present_index];
};



template <int dim>
void TriaObjectAccessor<1, dim>::set_boundary_indicator (unsigned char boundary_ind) const {
  Assert (dim>=2, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  tria->levels[present_level]->lines.material_id[present_index] = boundary_ind;
};



template <int dim>
bool TriaObjectAccessor<1, dim>::at_boundary () const {
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
};



template <int dim>
double TriaObjectAccessor<1, dim>::diameter () const {
  return sqrt((vertex(1)-vertex(0)).square());
};



template <int dim>
Point<dim> TriaObjectAccessor<1, dim>::center () const {
  return (vertex(1)+vertex(0))/2.;
};



template <int dim>
Point<dim> TriaObjectAccessor<1, dim>::barycenter () const {
  return (vertex(1)+vertex(0))/2.;
};



template <int dim>
double TriaObjectAccessor<1, dim>::measure () const {
  return sqrt((vertex(1)-vertex(0)).square());
};



template <int dim>
unsigned int TriaObjectAccessor<1, dim>::number_of_children () const {
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<2; ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
};




/*------------------------ Functions: QuadAccessor ---------------------------*/

template <int dim>
void TriaObjectAccessor<2, dim>::set (const Quad &quad) const {
  tria->levels[present_level]->quads.quads[present_index] = quad;
};



template <int dim>
int TriaObjectAccessor<2, dim>::vertex_index (const unsigned int corner) const {
  Assert (corner<4, ExcIndexRange(corner,0,4));

  const int corner_convention[4] = { 0,0,1,1 };
  return line(corner)->vertex_index(corner_convention[corner]);
};



template <int dim>
Point<dim> &
TriaObjectAccessor<2, dim>::vertex (const unsigned int i) const {
  return tria->vertices[vertex_index(i)];
};




template <int dim>
void
TriaObjectAccessor<2, dim>::set_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->quads.used[present_index] = true;
};
  


template <int dim>
void TriaObjectAccessor<2, dim>::clear_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->quads.used[present_index] = false;
};



template <int dim>
void TriaObjectAccessor<2, dim>::recursively_set_user_flag () const {
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<4; ++c)
      child(c)->recursively_set_user_flag ();
};



template <int dim>
void TriaObjectAccessor<2, dim>::recursively_clear_user_flag () const {
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<4; ++c)
      child(c)->recursively_clear_user_flag ();
};



template <int dim>
void TriaObjectAccessor<2, dim>::set_user_pointer (void *p) const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_pointers[present_index] = p;
};



template <int dim>
void TriaObjectAccessor<2, dim>::clear_user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->quads.user_pointers[present_index] = 0;
};



template <int dim>
void * TriaObjectAccessor<2, dim>::user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->quads.user_pointers[present_index];
};



template <int dim>
void TriaObjectAccessor<2, dim>::set_children (const int index) const {
  Assert (used(), ExcCellNotUsed());
  Assert ((index==-1) ||
	  (!has_children() && (index>=0)), ExcCantSetChildren(index));

  tria->levels[present_level]->quads.children[present_index] = index;
};



template <int dim>
void TriaObjectAccessor<2, dim>::clear_children () const {
  set_children (-1);
};



template <int dim>
unsigned char TriaObjectAccessor<2, dim>::boundary_indicator () const {
  Assert (dim>=3, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  return tria->levels[present_level]->quads.material_id[present_index];
};



template <int dim>
void TriaObjectAccessor<2, dim>::set_boundary_indicator (unsigned char boundary_ind) const {
  Assert (dim>=3, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  tria->levels[present_level]->quads.material_id[present_index] = boundary_ind;
};



template <int dim>
bool TriaObjectAccessor<2, dim>::at_boundary () const {
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
};



template <int dim>
double TriaObjectAccessor<2, dim>::diameter () const {
  return sqrt(max((vertex(2)-vertex(0)).square(),
		  (vertex(3)-vertex(1)).square()));
};



template <int dim>
Point<dim> TriaObjectAccessor<2, dim>::center () const {
  return (vertex(0)+vertex(1)+vertex(2)+vertex(3))/4.;
};



#if deal_II_dimension == 2

template <>
Point<2> TriaObjectAccessor<2, 2>::barycenter () const {
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independant, so we write this function
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
  tphi[1] := xi*(1-eta):
  tphi[2] := xi*eta:
  tphi[3] := (1-xi)*eta:
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
  const double t4 = x[2]*x[3];
  const double t7 = x[0]*x[3];
  const double t9 = x[1]*x[2];
  const double t11 = x[0]*x[0];
  const double t13 = x[1]*x[1];
  const double t15 = x[2]*x[2];
  const double t17 = x[3]*x[3];
  const double t25 = -t1*y[0]+t1*y[1]-t4*y[2]+t4*y[3]+t7*y[0]-t9*y[1]
		     +t11*y[1]-t13*y[0]+t15*y[3]-t17*y[2]-t11*y[3]
		     +t13*y[2]-t15*y[1]+t17*y[0]+t9*y[2]-t7*y[3];
  const double t35 = 1/(-x[0]*y[3]-x[1]*y[0]+x[2]*y[3]+x[3]*y[0]
			+x[0]*y[1]-x[3]*y[2]+x[1]*y[2]-x[2]*y[1]);
  const double t38 = y[3]*x[3];
  const double t40 = y[0]*x[0];
  const double t42 = y[1]*x[1];
  const double t44 = y[2]*x[2];
  const double t50 = y[0]*y[0];
  const double t52 = y[1]*y[1];
  const double t54 = y[2]*y[2];
  const double t56 = y[3]*y[3];
  const double t62 = -t38*y[2]+t40*y[1]-t42*y[0]+t44*y[3]-t40*y[3]
		     +t42*y[2]-t44*y[1]+t38*y[0]-x[1]*t50+x[0]*t52
		     -x[3]*t54+x[2]*t56+x[3]*t50-x[2]*t52+x[1]*t54-x[0]*t56;
  return Point<2> (t25*t35/3.0, t62*t35/3.0);
};




template <>
double TriaObjectAccessor<2, 2>::measure () const {
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independant, so we write this function
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
  tphi[1] := xi*(1-eta):
  tphi[2] := xi*eta:
  tphi[3] := (1-xi)*eta:
  x_real := sum(x[s]*tphi[s], s=0..3):
  y_real := sum(y[s]*tphi[s], s=0..3):
  detJ := diff(x_real,xi)*diff(y_real,eta) - diff(x_real,eta)*diff(y_real,xi):

  measure := simplify ( int ( int (detJ, xi=0..1), eta=0..1)):
  readlib(C):

  C(measure, optimized);
*/

  const double x[4] = { vertex(0)(0),
			vertex(1)(0),
			vertex(2)(0),
			vertex(3)(0)  };
  const double y[4] = { vertex(0)(1),
			vertex(1)(1),
			vertex(2)(1),
			vertex(3)(1)  };

  return (-x[0]*y[3]/2.0-x[1]*y[0]/2.0+x[2]*y[3]/2.0+x[3]*y[0]/2.0+
	  x[0]*y[1]/2.0-x[3]*y[2]/2.0+x[1]*y[2]/2.0-x[2]*y[1]/2.0);
};

#endif



#if deal_II_dimension == 3

template <>
Point<3> TriaObjectAccessor<2, 3>::barycenter () const {
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independant, so we write this function
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
};




template <>
double TriaObjectAccessor<2, 3>::measure () const {
				   // the evaluation of the formulae
				   // is a bit tricky when done dimension
				   // independant, so we write this function
				   // for 2D and 3D separately
				   //
				   // for documentation, see the barycenter
				   // function above.
  Assert (false, ExcNotImplemented());
  return 0;
};

#endif


template <int dim>
unsigned int TriaObjectAccessor<2, dim>::number_of_children () const {
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<4; ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
};





/*------------------------ Functions: HexAccessor ---------------------------*/

template <int dim>
void HexAccessor<dim>::set (const Hexahedron &hex) const {
  tria->levels[present_level]->hexes.hexes[present_index] = hex;
};



template <int dim>
int HexAccessor<dim>::vertex_index (const unsigned int corner) const {
  Assert (corner<8, ExcIndexRange(corner,0,8));

				   // get the corner indices by asking
				   // either the front or the back face
				   // for its vertices.
  if (corner<4)
    return quad(0)->vertex_index(corner);
  else
    return quad(1)->vertex_index(corner-4);
};



template <int dim>
Point<dim> &
HexAccessor<dim>::vertex (const unsigned int i) const {
  return tria->vertices[vertex_index(i)];
};



template <int dim>
void HexAccessor<dim>::set_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->hexes.used[present_index] = true;
};
  


template <int dim>
void HexAccessor<dim>::clear_used_flag () const {
  Assert (state() == valid, ExcDereferenceInvalidObject());
  tria->levels[present_level]->hexes.used[present_index] = false;
};




template <int dim>
void HexAccessor<dim>::recursively_set_user_flag () const {
  set_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<8; ++c)
      child(c)->recursively_set_user_flag ();
};



template <int dim>
void HexAccessor<dim>::recursively_clear_user_flag () const {
  clear_user_flag ();

  if (has_children())
    for (unsigned int c=0; c<8; ++c)
      child(c)->recursively_clear_user_flag ();
};



template <int dim>
void HexAccessor<dim>::set_user_pointer (void *p) const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_pointers[present_index] = p;
};



template <int dim>
void HexAccessor<dim>::clear_user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  tria->levels[present_level]->hexes.user_pointers[present_index] = 0;
};



template <int dim>
void * HexAccessor<dim>::user_pointer () const {
  Assert (used(), ExcCellNotUsed());
  return tria->levels[present_level]->hexes.user_pointers[present_index];
};



template <int dim>
void HexAccessor<dim>::set_children (const int index) const {
  Assert (used(), ExcCellNotUsed());
  Assert ((index==-1) ||
	  (!has_children() && (index>=0)), ExcCantSetChildren(index));

  tria->levels[present_level]->hexes.children[present_index] = index;
};



template <int dim>
void HexAccessor<dim>::clear_children () const {
  set_children (-1);
};



template <int dim>
unsigned char HexAccessor<dim>::boundary_indicator () const {
  Assert (dim<4, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  return tria->levels[present_level]->hexes.material_id[present_index];
};



template <int dim>
void HexAccessor<dim>::set_boundary_indicator (unsigned char boundary_ind) const {
  Assert (dim<4, ExcNotUsefulForThisDimension());
  Assert (used(), ExcCellNotUsed());

  tria->levels[present_level]->hexes.material_id[present_index] = boundary_ind;
};



template <int dim>
bool HexAccessor<dim>::at_boundary () const {
				   // error checking is done
				   // in boundary_indicator()
  return (boundary_indicator() != 255);
};



template <int dim>
double HexAccessor<dim>::diameter () const {
  return sqrt(max( max((vertex(6)-vertex(0)).square(),
		       (vertex(7)-vertex(1)).square()),
		   max((vertex(4)-vertex(2)).square(),
		       (vertex(5)-vertex(3)).square()) ));
};



template <int dim>
Point<dim> HexAccessor<dim>::center () const {
  return (vertex(0)+vertex(1)+vertex(2)+vertex(3)+
	  vertex(4)+vertex(5)+vertex(6)+vertex(7))/8.;
};



#if deal_II_dimension == 3

template <>
Point<3> HexAccessor<3>::barycenter () const {
/*
  Get the computation of the barycenter by this little Maple script. We
  use the trilinear mapping of the unit hex to the real hex.

  Remember that the area of the hex is given by
  |K| = \int_K 1 dx dy dz = \int_{\hat K} |det J| d(xi) d(eta) d(zeta)
  and that the barycenter is given by
  \vec x_s = 1/|K| \int_K \vec x dx dy dz
           = 1/|K| \int_{\hat K} \vec x(xi,eta,zeta) |det J| d(xi) d(eta) d(zeta)

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
*/

  const double x[8] = { vertex(0)(0),
			vertex(1)(0),
			vertex(2)(0),
			vertex(3)(0),
			vertex(4)(0),
			vertex(5)(0),
			vertex(6)(0),
			vertex(7)(0)   };
  const double y[8] = { vertex(0)(1),
			vertex(1)(1),
			vertex(2)(1),
			vertex(3)(1),
			vertex(4)(1),
			vertex(5)(1),
			vertex(6)(1),
			vertex(7)(1)  };
  const double z[8] = { vertex(0)(2),
			vertex(1)(2),
			vertex(2)(2),
			vertex(3)(2),
			vertex(4)(2),
			vertex(5)(2),
			vertex(6)(2),
			vertex(7)(2)  };

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
};




template <>
double HexAccessor<3>::measure () const {
/*
  Get the computation of the measure by the same little Maple script as above.
*/

  const double x[8] = { vertex(0)(0),
			vertex(1)(0),
			vertex(2)(0),
			vertex(3)(0),
			vertex(4)(0),
			vertex(5)(0),
			vertex(6)(0),
			vertex(7)(0)   };
  const double y[8] = { vertex(0)(1),
			vertex(1)(1),
			vertex(2)(1),
			vertex(3)(1),
			vertex(4)(1),
			vertex(5)(1),
			vertex(6)(1),
			vertex(7)(1)  };
  const double z[8] = { vertex(0)(2),
			vertex(1)(2),
			vertex(2)(2),
			vertex(3)(2),
			vertex(4)(2),
			vertex(5)(2),
			vertex(6)(2),
			vertex(7)(2)  };

  double s1, s2, s3;
  
  s3 = -x[2]*y[7]*z[3]/12.0+x[1]*y[4]*z[0]/12.0-x[0]*y[3]*z[7]/12.0-z[1]*x
[2]*y[0]/12.0+x[5]*y[4]*z[0]/12.0+x[1]*y[5]*z[0]/12.0-x[1]*y[2]*z[0]/12.0-y[2]*
x[3]*z[7]/12.0-x[1]*y[0]*z[4]/12.0-y[1]*x[0]*z[2]/12.0-y[0]*x[7]*z[3]/12.0+y[5]
*x[4]*z[7]/12.0-x[0]*y[3]*z[4]/12.0-y[1]*x[0]*z[3]/12.0-z[0]*x[7]*y[4]/12.0+x
[1]*y[0]*z[3]/12.0-x[2]*y[3]*z[0]/12.0-y[1]*x[2]*z[5]/12.0;
  s2 = s3+y[0]*x[3]*z[4]/12.0+y[1]*x[2]*z[0]/12.0+y[2]*x[7]*z[3]/12.0+z[1]*
x[4]*y[0]/12.0-y[5]*x[4]*z[0]/12.0-y[1]*x[5]*z[0]/12.0+x[5]*y[7]*z[4]/12.0+y[1]
*x[3]*z[0]/12.0+x[0]*y[4]*z[7]/12.0+z[0]*x[7]*y[3]/12.0-z[2]*x[3]*y[0]/12.0+y
[1]*x[0]*z[5]/12.0-x[1]*y[5]*z[2]/12.0+y[1]*x[0]*z[4]/12.0+z[1]*x[2]*y[5]/12.0-
x[5]*y[0]*z[4]/12.0+y[2]*x[3]*z[1]/12.0+x[1]*y[0]*z[2]/12.0;
  s3 = -z[1]*x[0]*y[4]/12.0+z[1]*x[0]*y[2]/12.0-x[0]*y[7]*z[4]/12.0+z[0]*x
[4]*y[7]/12.0+x[2]*y[3]*z[7]/12.0+y[2]*x[6]*z[3]/12.0+y[1]*x[6]*z[2]/12.0-y[2]*
x[1]*z[3]/12.0+x[1]*y[2]*z[5]/12.0+y[0]*x[3]*z[7]/12.0+y[2]*x[3]*z[0]/12.0-y[0]
*x[4]*z[3]/12.0-y[2]*x[0]*z[3]/12.0-y[1]*x[4]*z[0]/12.0+z[1]*x[5]*y[0]/12.0+x
[1]*y[2]*z[6]/12.0-x[1]*y[0]*z[5]/12.0-x[2]*y[3]*z[1]/12.0;
  s1 = s3-z[0]*x[3]*y[4]/12.0+z[1]*x[0]*y[3]/12.0-z[1]*x[5]*y[2]/12.0+z[1]*
x[2]*y[6]/12.0-z[2]*x[3]*y[1]/12.0+z[2]*x[3]*y[6]/12.0+x[2]*y[1]*z[3]/12.0-z[1]
*x[3]*y[0]/12.0+x[2]*y[3]*z[6]/12.0-z[0]*x[3]*y[7]/12.0+y[5]*x[0]*z[4]/12.0-z
[1]*x[6]*y[2]/12.0-z[2]*x[6]*y[3]/12.0+z[2]*x[0]*y[3]/12.0+z[2]*x[3]*y[7]/12.0+
z[2]*x[1]*y[3]/12.0+y[1]*x[5]*z[2]/12.0-z[2]*x[7]*y[3]/12.0+s2;
  s3 = x[0]*y[7]*z[3]/12.0-z[1]*x[0]*y[5]/12.0-x[1]*y[3]*z[0]/12.0-x[1]*y
[6]*z[2]/12.0-x[2]*y[6]*z[3]/12.0-x[5]*y[4]*z[7]/12.0+z[0]*x[4]*y[3]/12.0+x[0]*
y[4]*z[3]/12.0-x[6]*y[5]*z[7]/12.0-x[6]*y[7]*z[2]/12.0+z[1]*x[6]*y[5]/12.0+y[6]
*x[4]*z[7]/12.0+y[1]*x[5]*z[6]/12.0+x[1]*y[6]*z[5]/12.0-y[5]*x[7]*z[4]/12.0+y
[0]*x[7]*z[4]/12.0-y[0]*x[4]*z[7]/12.0-z[5]*x[0]*y[4]/12.0;
  s2 = s3+z[6]*x[3]*y[7]/12.0+z[6]*x[7]*y[4]/12.0-z[6]*x[4]*y[7]/12.0-z[6]*
x[7]*y[3]/12.0+y[6]*x[5]*z[7]/12.0-y[6]*x[7]*z[5]/12.0-x[2]*y[5]*z[6]/12.0+z[6]
*x[2]*y[7]/12.0+z[6]*x[7]*y[5]/12.0-z[6]*x[5]*y[7]/12.0-z[6]*x[7]*y[2]/12.0+x
[6]*y[7]*z[5]/12.0+z[5]*x[7]*y[4]/12.0-z[5]*x[4]*y[7]/12.0+y[5]*x[1]*z[4]/12.0-
y[5]*x[6]*z[4]/12.0-y[5]*x[4]*z[1]/12.0+y[5]*x[4]*z[6]/12.0;
  s3 = y[6]*x[7]*z[3]/12.0+x[6]*y[5]*z[2]/12.0-x[5]*y[6]*z[2]/12.0-x[6]*y
[2]*z[5]/12.0+x[5]*y[2]*z[6]/12.0+x[6]*y[2]*z[7]/12.0+z[5]*x[4]*y[0]/12.0-y[6]*
x[2]*z[7]/12.0+x[2]*y[6]*z[5]/12.0+y[6]*x[7]*z[2]/12.0-y[6]*x[3]*z[7]/12.0-x[7]
*y[4]*z[3]/12.0-x[6]*y[4]*z[7]/12.0-x[5]*y[1]*z[4]/12.0+x[6]*y[7]*z[4]/12.0+x
[7]*y[3]*z[4]/12.0-z[1]*x[5]*y[6]/12.0+x[3]*y[4]*z[7]/12.0+x[5]*y[6]*z[4]/12.0;
  
  return s3+x[5]*y[4]*z[1]/12.0-x[5]*y[4]*z[6]/12.0-x[3]*y[7]*z[4]/12.0+x[6]*
y[3]*z[7]/12.0-y[1]*x[6]*z[5]/12.0-z[5]*x[1]*y[4]/12.0+z[5]*x[6]*y[4]/12.0+z[5]
*x[4]*y[1]/12.0-z[5]*x[4]*y[6]/12.0-x[6]*y[7]*z[3]/12.0-x[4]*y[3]*z[7]/12.0+x
[4]*y[7]*z[3]/12.0-x[1]*y[5]*z[6]/12.0-y[6]*x[7]*z[4]/12.0-y[1]*x[2]*z[6]/12.0-
y[2]*x[3]*z[6]/12.0+s2+s1+x[2]*y[0]*z[3]/12.0;
};

#endif


template <int dim>
unsigned int HexAccessor<dim>::number_of_children () const {
  if (!has_children())
    return 1;
  else
    {
      unsigned int sum = 0;
      for (unsigned int c=0; c<8; ++c)
	sum += child(c)->number_of_children();
      return sum;
    };
};





/*------------------------ Functions: CellAccessor<1> -----------------------*/


#if deal_II_dimension == 1

template <>
bool CellAccessor<1>::at_boundary () const {
  return at_boundary(0) || at_boundary(1);
};



template <>
unsigned char CellAccessor<1>::material_id () const {
  Assert (used(), TriaSubstructAccessor<1>::ExcCellNotUsed());
  return tria->levels[present_level]->lines.material_id[present_index];
};



template <>
void CellAccessor<1>::set_material_id (const unsigned char mat_id) const {
  Assert (used(), TriaSubstructAccessor<1>::ExcCellNotUsed());
  tria->levels[present_level]->lines.material_id[present_index]
    = mat_id;						 
};


#endif


/*------------------------ Functions: CellAccessor<2> -----------------------*/


#if deal_II_dimension == 2

template <>
bool CellAccessor<2>::at_boundary () const {
  return at_boundary(0) || at_boundary(1) || at_boundary(2) || at_boundary(3);
};



template <>
unsigned char CellAccessor<2>::material_id () const {
  Assert (used(), TriaSubstructAccessor<2>::ExcCellNotUsed());
  return tria->levels[present_level]->quads.material_id[present_index];
};



template <>
void CellAccessor<2>::set_material_id (const unsigned char mat_id) const {
  Assert (used(), TriaSubstructAccessor<2>::ExcCellNotUsed());
  tria->levels[present_level]->quads.material_id[present_index]
    = mat_id;						 
};

#endif



/*------------------------ Functions: CellAccessor<3> -----------------------*/


#if deal_II_dimension == 3

template <>
bool CellAccessor<3>::at_boundary () const {
  return (at_boundary(0) || at_boundary(1) ||
	  at_boundary(2) || at_boundary(3) ||
	  at_boundary(4) || at_boundary(5));
};



template <>
unsigned char CellAccessor<3>::material_id () const {
  Assert (used(), TriaSubstructAccessor<3>::ExcCellNotUsed());
  return tria->levels[present_level]->hexes.material_id[present_index];
};



template <>
void CellAccessor<3>::set_material_id (const unsigned char mat_id) const {
  Assert (used(), TriaSubstructAccessor<3>::ExcCellNotUsed());
  tria->levels[present_level]->hexes.material_id[present_index]
    = mat_id;						 
};

#endif


/*------------------------ Functions: CellAccessor<dim> -----------------------*/


template <int dim>
void CellAccessor<dim>::set_neighbor (const unsigned int i,
				      const TriaIterator<dim,CellAccessor<dim> > &pointer) const {
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  typename TriaSubstructAccessor<dim>::ExcInvalidNeighbor(i));

  if (pointer.state() == valid)
    {
      tria->levels[present_level]->
	neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].first
	= pointer->present_level;
      tria->levels[present_level]->
	neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].second
	= pointer->present_index;
    }
  else
    {
      tria->levels[present_level]->
	neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].first
	= -1;
      tria->levels[present_level]->
	neighbors[present_index*GeometryInfo<dim>::faces_per_cell+i].second
	= -1;
    };
};



template <int dim>
bool CellAccessor<dim>::at_boundary (const unsigned int i) const {
  Assert (used(), typename TriaSubstructAccessor<dim>::ExcCellNotUsed());
  Assert (i<GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (i,0,GeometryInfo<dim>::faces_per_cell));
  
  return (neighbor_index(i) == -1);
};





// explicit instantiations
template class TriaAccessor<deal_II_dimension>;
template class TriaObjectAccessor<1, deal_II_dimension>;
template class CellAccessor<deal_II_dimension>;
template class TriaRawIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaRawIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,TriaObjectAccessor<1, deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,CellAccessor<deal_II_dimension> >;

#if deal_II_dimension >= 2
template class TriaObjectAccessor<2, deal_II_dimension>;
template class TriaRawIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,TriaObjectAccessor<2, deal_II_dimension> >;
#endif

#if deal_II_dimension >= 3
template class HexAccessor<deal_II_dimension>;
template class TriaRawIterator<deal_II_dimension,HexAccessor<deal_II_dimension> >;
template class TriaIterator<deal_II_dimension,HexAccessor<deal_II_dimension> >;
template class TriaActiveIterator<deal_II_dimension,HexAccessor<deal_II_dimension> >;
#endif






