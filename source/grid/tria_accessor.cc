// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_levels.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_accessor.templates.h>
#include <deal.II/grid/tria_iterator.templates.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

// anonymous namespace for helper functions
namespace
{
  // given the number of face's child
  // (subface_no), return the number of the
  // subface concerning the FaceRefineCase of
  // the face
  inline
  unsigned int translate_subface_no(const TriaIterator<TriaAccessor<2, 3, 3> > &face,
                                    const unsigned int                           subface_no)
  {
    Assert(face->has_children(), ExcInternalError());
    Assert(subface_no<face->n_children(), ExcInternalError());

    if (face->child(subface_no)->has_children())
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
  unsigned int translate_subface_no(const TriaIterator<TriaAccessor<2, 3, 3> > &face,
                                    const unsigned int                           subface_no,
                                    const unsigned int                           subsubface_no)
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

    static const unsigned int e = numbers::invalid_unsigned_int;

    // array containing the translation of the
    // numbers,
    //
    // first index: subface_no
    // second index: subsubface_no
    // third index: does the first subface have children? -> no and yes
    static const unsigned int translated_subface_no[2][2][2]
    =
    {
      { {e,0},       // first  subface, first  subsubface, first_child_has_children==no and yes
        {e,1}
      },      // first  subface, second subsubface, first_child_has_children==no and yes
      { {1,2},       // second subface, first  subsubface, first_child_has_children==no and yes
        {2,3}
      }
    };     // second subface, second subsubface, first_child_has_children==no and yes

    Assert(translated_subface_no[subface_no][subsubface_no][first_child_has_children]!=e,
           ExcInternalError());

    return translated_subface_no[subface_no][subsubface_no][first_child_has_children];
  }


  template <int dim, int spacedim>
  Point<spacedim>
  barycenter (const TriaAccessor<1, dim, spacedim> &accessor)
  {
    return (accessor.vertex(1)+accessor.vertex(0))/2.;
  }


  Point<2>
  barycenter (const TriaAccessor<2, 2, 2> &accessor)
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

    const double x[4] = { accessor.vertex(0)(0),
                          accessor.vertex(1)(0),
                          accessor.vertex(2)(0),
                          accessor.vertex(3)(0)
                        };
    const double y[4] = { accessor.vertex(0)(1),
                          accessor.vertex(1)(1),
                          accessor.vertex(2)(1),
                          accessor.vertex(3)(1)
                        };
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



  Point<3>
  barycenter (const TriaAccessor<3,3,3> &accessor)
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

    const double x[8] = { accessor.vertex(0)(0),
                          accessor.vertex(1)(0),
                          accessor.vertex(5)(0),
                          accessor.vertex(4)(0),
                          accessor.vertex(2)(0),
                          accessor.vertex(3)(0),
                          accessor.vertex(7)(0),
                          accessor.vertex(6)(0)
                        };
    const double y[8] = { accessor.vertex(0)(1),
                          accessor.vertex(1)(1),
                          accessor.vertex(5)(1),
                          accessor.vertex(4)(1),
                          accessor.vertex(2)(1),
                          accessor.vertex(3)(1),
                          accessor.vertex(7)(1),
                          accessor.vertex(6)(1)
                        };
    const double z[8] = { accessor.vertex(0)(2),
                          accessor.vertex(1)(2),
                          accessor.vertex(5)(2),
                          accessor.vertex(4)(2),
                          accessor.vertex(2)(2),
                          accessor.vertex(3)(2),
                          accessor.vertex(7)(2),
                          accessor.vertex(6)(2)
                        };

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



  template <int structdim, int dim, int spacedim>
  Point<spacedim>
  barycenter (const TriaAccessor<structdim, dim, spacedim> &)
  {
    // this function catches all the cases not
    // explicitly handled above
    Assert (false, ExcNotImplemented());
    return Point<spacedim>();
  }




  template <int dim, int spacedim>
  double
  measure (const TriaAccessor<1, dim, spacedim> &accessor)
  {
    // remember that we use (dim-)linear
    // mappings
    return (accessor.vertex(1)-accessor.vertex(0)).norm();
  }



  double
  measure (const TriaAccessor<2,2,2> &accessor)
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

    const double x[4] = { accessor.vertex(0)(0),
                          accessor.vertex(1)(0),
                          accessor.vertex(2)(0),
                          accessor.vertex(3)(0)
                        };
    const double y[4] = { accessor.vertex(0)(1),
                          accessor.vertex(1)(1),
                          accessor.vertex(2)(1),
                          accessor.vertex(3)(1)
                        };

    return (-x[1]*y[0]+x[1]*y[3]+y[0]*x[2]+x[0]*y[1]-x[0]*y[2]-y[1]*x[3]-x[2]*y[3]+x[3]*y[2])/2;
  }


  double
  measure (const TriaAccessor<3, 3, 3> &accessor)
  {
    unsigned int vertex_indices[GeometryInfo<3>::vertices_per_cell];
    for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
      vertex_indices[i] = accessor.vertex_index(i);

    return GridTools::cell_measure<3>(accessor.get_triangulation().get_vertices(),
                                      vertex_indices);
  }


  // a 2d face in 3d space
  double measure (const dealii::TriaAccessor<2,3,3> &accessor)
  {
    // If the face is planar, the diagonal from vertex 0 to vertex 3,
    // v_03, should be in the plane P_012 of vertices 0, 1 and 2.  Get
    // the normal vector of P_012 and test if v_03 is orthogonal to
    // that. If so, the face is planar and computing its area is simple.
    const Tensor<1,3> v01 = accessor.vertex(1) - accessor.vertex(0);
    const Tensor<1,3> v02 = accessor.vertex(2) - accessor.vertex(0);

    Tensor<1,3> normal = cross_product_3d(v01, v02);

    const Tensor<1,3> v03 = accessor.vertex(3) - accessor.vertex(0);

    // check whether v03 does not lie in the plane of v01 and v02
    // (i.e., whether the face is not planar). we do so by checking
    // whether the triple product (v01 x v02) * v03 forms a positive
    // volume relative to |v01|*|v02|*|v03|. the test checks the
    // squares of these to avoid taking norms/square roots:
    if (std::abs((v03 * normal) * (v03 * normal) /
                 ((v03 * v03) * (v01 * v01) * (v02 * v02)))
        >=
        1e-24)
      {
        Assert (false,
                ExcMessage("Computing the measure of a nonplanar face is not implemented!"));
        return std::numeric_limits<double>::quiet_NaN();
      }

    // the face is planar. then its area is 1/2 of the norm of the
    // cross product of the two diagonals
    const Tensor<1,3> v12 = accessor.vertex(2) - accessor.vertex(1);
    Tensor<1,3> twice_area = cross_product_3d(v03, v12);
    return 0.5 * twice_area.norm();
  }



  template <int structdim, int dim, int spacedim>
  double
  measure (const TriaAccessor<structdim, dim, spacedim> &)
  {
    // catch-all for all cases not explicitly
    // listed above
    Assert (false, ExcNotImplemented());
    return std::numeric_limits<double>::quiet_NaN();
  }


  template <int dim, int spacedim>
  Point<spacedim> get_new_point_on_object(const TriaAccessor<1, dim, spacedim> &obj)
  {
    TriaIterator<TriaAccessor<1,dim,spacedim> > it(obj);
    return obj.get_manifold().get_new_point_on_line(it);
  }

  template <int dim, int spacedim>
  Point<spacedim> get_new_point_on_object(const TriaAccessor<2, dim, spacedim> &obj)
  {
    TriaIterator<TriaAccessor<2,dim,spacedim> > it(obj);
    return obj.get_manifold().get_new_point_on_quad(it);
  }

  template <int dim, int spacedim>
  Point<spacedim> get_new_point_on_object(const TriaAccessor<3, dim, spacedim> &obj)
  {
    TriaIterator<TriaAccessor<3,dim,spacedim> > it(obj);
    return obj.get_manifold().get_new_point_on_hex(it);
  }

  template <int structdim, int dim, int spacedim>
  Point<spacedim> get_new_point_on_object(const TriaAccessor<structdim, dim, spacedim> &obj,
                                          const bool use_laplace)
  {
    if (use_laplace == false)
      return get_new_point_on_object(obj);
    else
      {
        TriaRawIterator<TriaAccessor<structdim, dim, spacedim> > it(obj);
        Quadrature<spacedim> quadrature = Manifolds::get_default_quadrature(it, use_laplace);
        return obj.get_manifold().get_new_point(quadrature);
      }
  }
}



/*------------------------ Static variables: TriaAccessorBase ---------------------------*/

template <int structdim, int dim, int spacedim>
const unsigned int TriaAccessorBase<structdim, dim, spacedim>::dimension;

template <int structdim, int dim, int spacedim>
const unsigned int TriaAccessorBase<structdim, dim, spacedim>::space_dimension;

template <int structdim, int dim, int spacedim>
const unsigned int TriaAccessorBase<structdim, dim, spacedim>::structure_dimension;


/*------------------------ Functions: TriaAccessor ---------------------------*/

template <int structdim, int dim, int spacedim>
void
TriaAccessor<structdim, dim, spacedim>::
set (const internal::Triangulation::TriaObject<structdim> &object) const
{
  this->objects().cells[this->present_index] = object;
}



template <int structdim, int dim, int spacedim>
Point<spacedim>
TriaAccessor<structdim, dim, spacedim>::
barycenter () const
{
  // call the function in the anonymous
  // namespace above
  return dealii::barycenter (*this);
}



template <int structdim, int dim, int spacedim>
double
TriaAccessor<structdim, dim, spacedim>::
measure () const
{
  // call the function in the anonymous
  // namespace above
  return dealii::measure (*this);
}



template <>
double TriaAccessor<1,1,1>::extent_in_direction(const unsigned int axis) const
{
  (void)axis;
  Assert (axis == 0, ExcIndexRange (axis, 0, 1));

  return this->diameter();
}


template <>
double TriaAccessor<1,1,2>::extent_in_direction(const unsigned int axis) const
{
  (void)axis;
  Assert (axis == 0, ExcIndexRange (axis, 0, 1));

  return this->diameter();
}


template <>
double TriaAccessor<2,2,2>::extent_in_direction(const unsigned int axis) const
{
  const unsigned int lines[2][2] = {{2,3}, /// Lines along x-axis, see GeometryInfo
    {0,1}
  };/// Lines along y-axis

  Assert (axis < 2, ExcIndexRange (axis, 0, 2));

  return std::max(this->line(lines[axis][0])->diameter(),
                  this->line(lines[axis][1])->diameter());
}

template <>
double TriaAccessor<2,2,3>::extent_in_direction(const unsigned int axis) const
{
  const unsigned int lines[2][2] = {{2,3}, /// Lines along x-axis, see GeometryInfo
    {0,1}
  };/// Lines along y-axis

  Assert (axis < 2, ExcIndexRange (axis, 0, 2));

  return std::max(this->line(lines[axis][0])->diameter(),
                  this->line(lines[axis][1])->diameter());
}


template <>
double TriaAccessor<3,3,3>::extent_in_direction(const unsigned int axis) const
{
  const unsigned int lines[3][4] = {{2,3,6,7},     /// Lines along x-axis, see GeometryInfo
    {0,1,4,5},    /// Lines along y-axis
    {8,9,10,11}
  }; /// Lines along z-axis

  Assert (axis < 3, ExcIndexRange (axis, 0, 3));

  double lengths[4] = { this->line(lines[axis][0])->diameter(),
                        this->line(lines[axis][1])->diameter(),
                        this->line(lines[axis][2])->diameter(),
                        this->line(lines[axis][3])->diameter()
                      };

  return std::max(std::max(lengths[0], lengths[1]),
                  std::max(lengths[2], lengths[3]));
}


// Recursively set manifold ids on hex iterators.
template <>
void
TriaAccessor<3,3,3>::
set_all_manifold_ids (const types::manifold_id manifold_ind) const
{
  set_manifold_id (manifold_ind);

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->set_all_manifold_ids (manifold_ind);

  // for hexes also set manifold_id
  // of bounding quads and lines

  // Six bonding quads
  for (unsigned int i=0; i<6; ++i)
    this->quad(i)->set_manifold_id(manifold_ind);
  // Twelve bounding lines
  for (unsigned int i=0; i<12; ++i)
    this->line(i)->set_manifold_id (manifold_ind);
}


template <int structdim, int dim, int spacedim>
Point<spacedim>
TriaAccessor<structdim, dim, spacedim>::intermediate_point (const Point<structdim> &coordinates) const
{
  // We use an FE_Q<structdim>(1) to extract the "weights" of each
  // vertex, used to get a point from the manifold.
  static FE_Q<structdim> fe(1);

  // Surrounding points and weights.
  std::vector<Point<spacedim> > p(GeometryInfo<structdim>::vertices_per_cell);
  std::vector<double>   w(GeometryInfo<structdim>::vertices_per_cell);

  for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
    {
      p[i] = this->vertex(i);
      w[i] = fe.shape_value(i, coordinates);
    }

  Quadrature<spacedim> quadrature(p, w);
  return this->get_manifold().get_new_point(quadrature);
}


template <int structdim, int dim, int spacedim>
Point<spacedim>
TriaAccessor<structdim, dim, spacedim>::center (const bool respect_manifold,
                                                const bool use_laplace) const
{
  if (respect_manifold == false)
    {
      Assert(use_laplace == false, ExcNotImplemented());
      Point<spacedim> p;
      for (unsigned int v=0; v<GeometryInfo<structdim>::vertices_per_cell; ++v)
        p += vertex(v);
      return p/GeometryInfo<structdim>::vertices_per_cell;
    }
  else
    return get_new_point_on_object(*this, use_laplace);
}


/*------------------------ Functions: CellAccessor<1> -----------------------*/



template <>
bool CellAccessor<1>::point_inside (const Point<1> &p) const
{
  return (this->vertex(0)[0] <= p[0]) && (p[0] <= this->vertex(1)[0]);
}






/*------------------------ Functions: CellAccessor<2> -----------------------*/



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
  static const int direction[4]= {-1,1,1,-1};
  for (unsigned int f=0; f<4; ++f)
    {
      // vector from the first vertex
      // of the line to the point
      const Tensor<1,2> to_p = p-this->vertex(
                                 GeometryInfo<2>::face_to_cell_vertices(f,0));
      // vector describing the line
      const Tensor<1,2> face = direction[f]*(
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
      if ((-face[1]*to_p[0]+face[0]*to_p[1])<0)
        return false;
    };

  // if we arrived here, then the
  // point is inside for all four
  // faces, and thus inside
  return true;
}







/*------------------------ Functions: CellAccessor<3> -----------------------*/



template <>
bool CellAccessor<3>::point_inside (const Point<3> &p) const
{
  // original implementation by Joerg
  // Weimar

  // we first eliminate points based
  // on the maximum and minimum of
  // the corner coordinates, then
  // transform to the unit cell, and
  // check there.
  const unsigned int dim = 3;
  const unsigned int spacedim = 3;
  Point<spacedim> maxp = this->vertex(0);
  Point<spacedim> minp = this->vertex(0);

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

  // now we need to check more carefully: transform to the
  // unit cube and check there. unfortunately, this isn't
  // completely trivial since the transform_real_to_unit_cell
  // function may throw an exception that indicates that the
  // point given could not be inverted. we take this as a sign
  // that the point actually lies outside, as also documented
  // for that function
  try
    {
      const TriaRawIterator< CellAccessor<dim,spacedim> > cell_iterator (*this);
      return (GeometryInfo<dim>::is_inside_unit_cell
              (StaticMappingQ1<dim,spacedim>::mapping.transform_real_to_unit_cell(cell_iterator, p)));
    }
  catch (const Mapping<dim,spacedim>::ExcTransformationFailed &)
    {
      return false;
    }
}





/*------------------------ Functions: CellAccessor<dim,spacedim> -----------------------*/

// For codim>0 we proceed as follows:
// 1) project point onto manifold and
// 2) transform to the unit cell with a Q1 mapping
// 3) then check if inside unit cell
template<int dim, int spacedim>
template<int dim_,int spacedim_ >
bool CellAccessor<dim,spacedim>::
point_inside_codim(const Point<spacedim_> &p) const
{
  const TriaRawIterator< CellAccessor<dim_,spacedim_> > cell_iterator (*this);
  const Point< dim_ > p_unit =
    StaticMappingQ1<dim_,spacedim_>::mapping.transform_real_to_unit_cell(cell_iterator, p);

  return GeometryInfo< dim_ >::is_inside_unit_cell(p_unit);

}



template <>
bool CellAccessor<1,2>::point_inside (const Point<2> &p) const
{
  return point_inside_codim<1,2>(p);
}


template <>
bool CellAccessor<1,3>::point_inside (const Point<3> &p) const
{
  return point_inside_codim<1,3>(p);
}


template <>
bool CellAccessor<2,3>::point_inside (const Point<3> &p) const
{
  return point_inside_codim<2,3>(p);
}



template <int dim, int spacedim>
bool CellAccessor<dim, spacedim>::at_boundary () const
{
  switch (dim)
    {
    case 1:
      return at_boundary(0) || at_boundary(1);
    case 2:
      return (at_boundary(0) || at_boundary(1) ||
              at_boundary(2) || at_boundary(3));
    case 3:
      return (at_boundary(0) || at_boundary(1) ||
              at_boundary(2) || at_boundary(3) ||
              at_boundary(4) || at_boundary(5));
    default:
      Assert (false, ExcNotImplemented());
      return false;
    }
}



template <int dim, int spacedim>
types::material_id CellAccessor<dim, spacedim>::material_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->cells.boundary_or_material_id[this->present_index].material_id;
}



template <int dim, int spacedim>
void CellAccessor<dim, spacedim>::set_material_id (const types::material_id mat_id) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert ( mat_id < numbers::invalid_material_id, ExcIndexRange(mat_id,0,numbers::invalid_material_id));
  this->tria->levels[this->present_level]->cells.boundary_or_material_id[this->present_index].material_id = mat_id;
}



template <int dim, int spacedim>
void CellAccessor<dim, spacedim>::recursively_set_material_id (const types::material_id mat_id) const
{
  set_material_id (mat_id);

  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_set_material_id (mat_id);
}



template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::set_subdomain_id (const types::subdomain_id new_subdomain_id) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (this->active(),
          ExcMessage("set_subdomain_id() can only be called on active cells!"));
  this->tria->levels[this->present_level]->subdomain_ids[this->present_index]
    = new_subdomain_id;
}


template <int dim, int spacedim>
types::subdomain_id CellAccessor<dim, spacedim>::level_subdomain_id () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  return this->tria->levels[this->present_level]->level_subdomain_ids[this->present_index];
}



template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::set_level_subdomain_id (const types::subdomain_id new_level_subdomain_id) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  this->tria->levels[this->present_level]->level_subdomain_ids[this->present_index]
    = new_level_subdomain_id;
}


template <int dim, int spacedim>
bool CellAccessor<dim, spacedim>::direction_flag () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  if (dim==spacedim)
    return true;
  else
    return this->tria->levels[this->present_level]->direction_flags[this->present_index];
}



template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::set_direction_flag (const bool new_direction_flag) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  if (dim<spacedim)
    this->tria->levels[this->present_level]->direction_flags[this->present_index]
      = new_direction_flag;
  else
    Assert (new_direction_flag == true,
            ExcMessage ("If dim==spacedim, direction flags are always true and "
                        "can not be set to anything else."));
}



template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::set_active_cell_index (const unsigned int active_cell_index)
{
  // set the active cell index. allow setting it also for non-active (and unused)
  // cells to allow resetting the index after refinement
  this->tria->levels[this->present_level]->active_cell_indices[this->present_index]
    = active_cell_index;
}



template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::set_parent (const unsigned int parent_index)
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (this->present_level > 0, TriaAccessorExceptions::ExcCellHasNoParent ());
  this->tria->levels[this->present_level]->parents[this->present_index / 2]
    = parent_index;
}



template <int dim, int spacedim>
int
CellAccessor<dim, spacedim>::
parent_index () const
{
  Assert (this->present_level > 0, TriaAccessorExceptions::ExcCellHasNoParent ());

  // the parent of two consecutive cells
  // is stored only once, since it is
  // the same
  return this->tria->levels[this->present_level]->parents[this->present_index / 2];
}



template <int dim, int spacedim>
unsigned int
CellAccessor<dim, spacedim>::
active_cell_index () const
{
  Assert (this->has_children()==false, TriaAccessorExceptions::ExcCellNotActive());
  return this->tria->levels[this->present_level]->active_cell_indices[this->present_index];
}



template <int dim, int spacedim>
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim, spacedim>::parent () const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (this->present_level > 0, TriaAccessorExceptions::ExcCellHasNoParent ());
  TriaIterator<CellAccessor<dim,spacedim> >
  q (this->tria, this->present_level-1, parent_index ());

  return q;
}


template <int dim, int spacedim>
void
CellAccessor<dim, spacedim>::
recursively_set_subdomain_id (const types::subdomain_id new_subdomain_id) const
{
  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); ++c)
      this->child(c)->recursively_set_subdomain_id (new_subdomain_id);
  else
    set_subdomain_id (new_subdomain_id);
}



template <int dim, int spacedim>
void CellAccessor<dim, spacedim>::set_neighbor (const unsigned int i,
                                                const TriaIterator<CellAccessor<dim, spacedim> > &pointer) const
{
  AssertIndexRange (i, GeometryInfo<dim>::faces_per_cell);

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



template <int dim, int spacedim>
unsigned int CellAccessor<dim,spacedim>::neighbor_of_neighbor_internal (const unsigned int neighbor) const
{
  AssertIndexRange (neighbor, GeometryInfo<dim>::faces_per_cell);

  // if we have a 1d mesh in 1d, we
  // can assume that the left
  // neighbor of the right neighbor is
  // the current cell. but that is an
  // invariant that isn't true if the
  // mesh is embedded in a higher
  // dimensional space, so we have to
  // fall back onto the generic code
  // below
  if ((dim==1) && (spacedim==dim))
    return GeometryInfo<dim>::opposite_face[neighbor];

  const TriaIterator<CellAccessor<dim, spacedim> > neighbor_cell = this->neighbor(neighbor);

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



template <int dim, int spacedim>
unsigned int CellAccessor<dim,spacedim>::neighbor_of_neighbor (const unsigned int neighbor) const
{
  const unsigned int n2=neighbor_of_neighbor_internal(neighbor);
  Assert (n2!=numbers::invalid_unsigned_int,
          TriaAccessorExceptions::ExcNeighborIsCoarser());

  return n2;
}



template <int dim, int spacedim>
bool
CellAccessor<dim,spacedim>::neighbor_is_coarser (const unsigned int neighbor) const
{
  return neighbor_of_neighbor_internal(neighbor)==numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
std::pair<unsigned int, unsigned int>
CellAccessor<dim, spacedim>::neighbor_of_coarser_neighbor (const unsigned int neighbor) const
{
  AssertIndexRange (neighbor, GeometryInfo<dim>::faces_per_cell);
  // make sure that the neighbor is
  // on a coarser level
  Assert (neighbor_is_coarser(neighbor),
          TriaAccessorExceptions::ExcNeighborIsNotCoarser());

  switch (dim)
    {
    case 2:
    {
      const int this_face_index=face_index(neighbor);
      const TriaIterator<CellAccessor<2,spacedim> > neighbor_cell = this->neighbor(neighbor);

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
        = GeometryInfo<2>::opposite_face[neighbor];

      const TriaIterator<TriaAccessor<1, 2, spacedim> > face_guess
        =neighbor_cell->face(face_no_guess);

      if (face_guess->has_children())
        for (unsigned int subface_no=0; subface_no<face_guess->n_children(); ++subface_no)
          if (face_guess->child_index(subface_no)==this_face_index)
            return std::make_pair (face_no_guess, subface_no);

      // if the guess was false, then
      // we need to loop over all faces
      // and subfaces and find the
      // number the hard way
      for (unsigned int face_no=0; face_no<GeometryInfo<2>::faces_per_cell; ++face_no)
        {
          if (face_no!=face_no_guess)
            {
              const TriaIterator<TriaAccessor<1, 2, spacedim> > face
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
      return std::make_pair (numbers::invalid_unsigned_int,
                             numbers::invalid_unsigned_int);
    }

    case 3:
    {
      const int this_face_index=face_index(neighbor);
      const TriaIterator<CellAccessor<3, spacedim> >
      neighbor_cell = this->neighbor(neighbor);

      // usually, on regular patches of the grid, this cell is just on the
      // opposite side of the neighbor that the neighbor is of this cell.
      // for example in 2d, if we want to know the neighbor_of_neighbor if
      // neighbor==1 (the right neighbor), then we will get 0 (the left
      // neighbor) in most cases. look up this relationship in the table
      // provided by GeometryInfo and try it
      const unsigned int face_no_guess
        = GeometryInfo<3>::opposite_face[neighbor];

      const TriaIterator<TriaAccessor<3-1, 3, spacedim> > face_guess
        =neighbor_cell->face(face_no_guess);

      if (face_guess->has_children())
        for (unsigned int subface_no=0; subface_no<face_guess->n_children(); ++subface_no)
          {
            if (face_guess->child_index(subface_no)==this_face_index)
              // call a helper function, that translates the current
              // subface number to a subface number for the current
              // FaceRefineCase
              return std::make_pair (face_no_guess, translate_subface_no(face_guess, subface_no));

            if (face_guess->child(subface_no)->has_children())
              for (unsigned int subsub_no=0; subsub_no<face_guess->child(subface_no)->n_children(); ++subsub_no)
                if (face_guess->child(subface_no)->child_index(subsub_no)==this_face_index)
                  // call a helper function, that translates the current
                  // subface number and subsubface number to a subface
                  // number for the current FaceRefineCase
                  return std::make_pair (face_no_guess, translate_subface_no(face_guess, subface_no, subsub_no));
          }

      // if the guess was false, then we need to loop over all faces and
      // subfaces and find the number the hard way
      for (unsigned int face_no=0; face_no<GeometryInfo<3>::faces_per_cell; ++face_no)
        {
          if (face_no==face_no_guess)
            continue;

          const TriaIterator<TriaAccessor<3-1, 3, spacedim> > face
            =neighbor_cell->face(face_no);

          if (!face->has_children())
            continue;

          for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
            {
              if (face->child_index(subface_no)==this_face_index)
                // call a helper function, that translates the current
                // subface number to a subface number for the current
                // FaceRefineCase
                return std::make_pair (face_no, translate_subface_no(face, subface_no));

              if (face->child(subface_no)->has_children())
                for (unsigned int subsub_no=0; subsub_no<face->child(subface_no)->n_children(); ++subsub_no)
                  if (face->child(subface_no)->child_index(subsub_no)==this_face_index)
                    // call a helper function, that translates the current
                    // subface number and subsubface number to a subface
                    // number for the current FaceRefineCase
                    return std::make_pair (face_no, translate_subface_no(face, subface_no, subsub_no));
            }
        }

      // we should never get here, since then we did not find our way
      // back...
      Assert (false, ExcInternalError());
      return std::make_pair (numbers::invalid_unsigned_int,
                             numbers::invalid_unsigned_int);
    }

    default:
    {
      Assert(false, ExcImpossibleInDim(1));
      return std::make_pair (numbers::invalid_unsigned_int,
                             numbers::invalid_unsigned_int);
    }
    }
}



template <int dim, int spacedim>
bool
CellAccessor<dim, spacedim>::has_periodic_neighbor (const unsigned int i_face) const
{
  /*
   * Implementation note: In all of the functions corresponding to periodic faces
   * we mainly use the Triangulation::periodic_face_map to find the
   * information about periodically connected faces. So, we actually search in
   * this std::map and return the cell_face on the other side of the periodic
   * boundary. For this search process, we have two options:
   *
   * 1- Using the [] operator of std::map: This option results in a more readalbe
   *    code, but requires an extra iteration in the map. Because when we call [] on
   *    std::map, with a key which does not exist in the std::map, that key will be
   *    created and the default value will be returned by []. This is not desirable.
   *    So, one has to first check if the key exists in the std::map and if it exists,
   *    then use the [] operator. The existence check is possible using std::map::find()
   *    or std::map::count(). Using this option will result in two iteration cycles
   *    through the map. First, existence check, then returning the value.
   *
   * 2- Using std::map::find(): This option is less readble, but theoretically
   *    faster, because it results in one iteration through std::map object.
   *
   * We decided to use the 2nd option.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  // my_it : is the iterator to the current cell.
  cell_iterator my_it(*this);
  if (this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face))
      != this->tria->periodic_face_map.end())
    return true;
  return false;
}



template <int dim, int spacedim>
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim, spacedim>::
periodic_neighbor (const unsigned int i_face) const
{
  /*
   * To know, why we are using std::map::find() instead of [] operator, refer
   * to the implementation note in has_periodic_neighbor() function.
   *
   * my_it        : the iterator to the current cell.
   * my_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the current cell_face.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  cell_iterator my_it(*this);

  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        my_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face));
  // Assertion is required to check that we are actually on a periodic boundary.
  Assert (my_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  return my_face_pair->second.first.first;
}



template <int dim, int spacedim>
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim, spacedim>::
neighbor_or_periodic_neighbor (const unsigned int i_face) const
{
  if (!(this->face(i_face)->at_boundary()))
    return this->neighbor(i_face);
  else if (this->has_periodic_neighbor(i_face))
    return this->periodic_neighbor(i_face);
  else
    AssertThrow (false, TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  // we can't come here
  return this->neighbor(i_face);
}



template <int dim, int spacedim>
TriaIterator<CellAccessor<dim, spacedim> >
CellAccessor<dim, spacedim>::
periodic_neighbor_child_on_subface (const unsigned int i_face,
                                    const unsigned int i_subface) const
{
  /*
   * To know, why we are using std::map::find() instead of [] operator, refer
   * to the implementation note in has_periodic_neighbor() function.
   *
   * my_it            : the iterator to the current cell.
   * my_face_pair     : the pair reported by periodic_face_map as its first pair being
   *                    the current cell_face.
   * nb_it            : the iterator to the neighbor of current cell at i_face.
   * face_num_of_nb   : the face number of the periodically neighboring face in the
   *                    relevant element.
   * nb_parent_face_it: the iterator to the parent face of the periodically
   *                    neighboring face.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  cell_iterator my_it(*this);
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        my_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face));
  /*
   * There should be an assertion, which tells the user that this function should not be
   * used for a cell which is not located at a periodic boundary.
   */
  Assert (my_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  cell_iterator parent_nb_it = my_face_pair->second.first.first;
  unsigned int nb_face_num = my_face_pair->second.first.second;
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim> > nb_parent_face_it = parent_nb_it->face(nb_face_num);
  /*
   * We should check if the parent face of the neighbor has at least the same number of
   * children as i_subface.
   */
  AssertIndexRange (i_subface, nb_parent_face_it->n_children());
  unsigned int sub_neighbor_num =
    GeometryInfo<dim>::child_cell_on_face(parent_nb_it->refinement_case(),
                                          nb_face_num,
                                          i_subface,
                                          my_face_pair->second.second[0],
                                          my_face_pair->second.second[1],
                                          my_face_pair->second.second[2],
                                          nb_parent_face_it->refinement_case());
  return parent_nb_it->child(sub_neighbor_num);
}



template <int dim, int spacedim>
std::pair<unsigned int, unsigned int>
CellAccessor<dim, spacedim>::
periodic_neighbor_of_coarser_periodic_neighbor (const unsigned int i_face) const
{
  /*
   * To know, why we are using std::map::find() instead of [] operator, refer
   * to the implementation note in has_periodic_neighbor() function.
   *
   * my_it        : the iterator to the current cell.
   * my_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the current cell_face.
   * nb_it        : the iterator to the periodic neighbor.
   * nb_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the periodic neighbor cell_face.
   * p_nb_of_p_nb : the iterator of the periodic neighbor of the periodic neighbor
   *                of the current cell.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  const int my_face_index = this->face_index(i_face);
  cell_iterator my_it(*this);
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        my_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face));
  /*
   * There should be an assertion, which tells the user that this function should not be
   * used for a cell which is not located at a periodic boundary.
   */
  Assert (my_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  cell_iterator nb_it = my_face_pair->second.first.first;
  unsigned int face_num_of_nb = my_face_pair->second.first.second;
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        nb_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(nb_it, face_num_of_nb));
  /*
   * Since, we store periodic neighbors for every cell (either active or artificial or inactive)
   * the nb_face_pair should also be mapped to some cell_face pair. We assert this here.
   */
  Assert (nb_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  cell_iterator p_nb_of_p_nb = nb_face_pair->second.first.first;
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim> > parent_face_it = p_nb_of_p_nb->face(nb_face_pair->second.first.second);
  for (unsigned int i_subface = 0; i_subface < parent_face_it->n_children(); ++i_subface)
    if (parent_face_it->child_index(i_subface) == my_face_index)
      return (std::pair<unsigned int, unsigned int>(face_num_of_nb, i_subface));
  /*
   * Obviously, if the execution reaches to this point, some of our assumptions should have
   * been false. The most important one is, the user has called this function on a face
   * which does not have a coarser periodic neighbor.
   */
  Assert (false, TriaAccessorExceptions::ExcNeighborIsNotCoarser());
  return std::pair<unsigned int, unsigned int>(numbers::invalid_unsigned_int,
                                               numbers::invalid_unsigned_int);
}



template <int dim, int spacedim>
int CellAccessor<dim, spacedim>::
periodic_neighbor_index(const unsigned int i_face) const
{
  return periodic_neighbor(i_face)->index();
}



template <int dim, int spacedim>
int CellAccessor<dim, spacedim>::
periodic_neighbor_level(const unsigned int i_face) const
{
  return periodic_neighbor(i_face)->level();
}



template <int dim, int spacedim>
unsigned int CellAccessor<dim, spacedim>::
periodic_neighbor_of_periodic_neighbor(const unsigned int i_face) const
{
  return periodic_neighbor_face_no(i_face);
}



template <int dim, int spacedim>
unsigned int
CellAccessor<dim, spacedim>::periodic_neighbor_face_no (const unsigned int i_face) const
{
  /*
   * To know, why we are using std::map::find() instead of [] operator, refer
   * to the implementation note in has_periodic_neighbor() function.
   *
   * my_it        : the iterator to the current cell.
   * my_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the current cell_face.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  cell_iterator my_it(*this);
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        my_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face));
  /*
   * There should be an assertion, which tells the user that this function should not be
   * called for a cell which is not located at a periodic boundary !
   */
  Assert (my_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  return my_face_pair->second.first.second;
}



template <int dim, int spacedim>
bool
CellAccessor<dim, spacedim>::periodic_neighbor_is_coarser (const unsigned int i_face) const
{
  /*
   * To know, why we are using std::map::find() instead of [] operator, refer
   * to the implementation note in has_periodic_neighbor() function.
   *
   * Implementation note: Let p_nb_of_p_nb be the periodic neighbor of the periodic
   * neighbor of the current cell. Also, let p_face_of_p_nb_of_p_nb be the periodic
   * face of the p_nb_of_p_nb. If p_face_of_p_nb_of_p_nb has children , then the
   * periodic neighbor of the current cell is coarser than itself. Although not tested,
   * this implementation should work for anisotropic refinement as well.
   *
   * my_it        : the iterator to the current cell.
   * my_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the current cell_face.
   * nb_it        : the iterator to the periodic neighbor.
   * nb_face_pair : the pair reported by periodic_face_map as its first pair being
   *                the periodic neighbor cell_face.
   */
  AssertIndexRange (i_face, GeometryInfo<dim>::faces_per_cell);
  typedef TriaIterator<CellAccessor<dim, spacedim> > cell_iterator;
  cell_iterator my_it(*this);
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        my_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(my_it, i_face));
  /*
   * There should be an assertion, which tells the user that this function should not be
   * used for a cell which is not located at a periodic boundary.
   */
  Assert (my_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  cell_iterator nb_it = my_face_pair->second.first.first;
  unsigned int face_num_of_nb = my_face_pair->second.first.second;
  const typename std::map<std::pair<cell_iterator, unsigned int>,
        std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3> > >::const_iterator
        nb_face_pair =
          this->tria->periodic_face_map.find(std::pair<cell_iterator, unsigned int>(nb_it, face_num_of_nb));
  /*
   * Since, we store periodic neighbors for every cell (either active or artificial or inactive)
   * the nb_face_pair should also be mapped to some cell_face pair. We assert this here.
   */
  Assert (nb_face_pair != this->tria->periodic_face_map.end(),
          TriaAccessorExceptions::ExcNoPeriodicNeighbor());
  const unsigned int my_level = this->level();
  const unsigned int neighbor_level = nb_face_pair->second.first.first->level();
  Assert (my_level >= neighbor_level, ExcInternalError());
  return my_level > neighbor_level;
}



template <int dim, int spacedim>
bool CellAccessor<dim, spacedim>::at_boundary (const unsigned int i) const
{
  Assert (this->used(), TriaAccessorExceptions::ExcCellNotUsed());
  Assert (i<GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (i,0,GeometryInfo<dim>::faces_per_cell));

  return (neighbor_index(i) == -1);
}



template <int dim, int spacedim>
bool CellAccessor<dim, spacedim>::has_boundary_lines () const
{
  if (dim == 1)
    return at_boundary ();
  else
    {
      for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_cell; ++l)
        if (this->line(l)->at_boundary())
          return true;

      return false;
    }
}



template <int dim, int spacedim>
TriaIterator<CellAccessor<dim,spacedim> >
CellAccessor<dim, spacedim>::
neighbor_child_on_subface (const unsigned int face,
                           const unsigned int subface) const
{
  Assert (!this->has_children(),
          ExcMessage ("The present cell must not have children!"));
  Assert (!this->at_boundary(face),
          ExcMessage ("The present cell must have a valid neighbor!"));
  Assert (this->neighbor(face)->has_children() == true,
          ExcMessage ("The neighbor must have children!"));

  switch (dim)
    {
    case 2:
    {
      const unsigned int neighbor_neighbor
        = this->neighbor_of_neighbor (face);
      const unsigned int neighbor_child_index
        = GeometryInfo<dim>::child_cell_on_face(
            this->neighbor(face)->refinement_case(),neighbor_neighbor,subface);

      TriaIterator<CellAccessor<dim,spacedim> > sub_neighbor
        = this->neighbor(face)->child(neighbor_child_index);
      // the neighbors child can have children,
      // which are not further refined along the
      // face under consideration. as we are
      // normally interested in one of this
      // child's child, search for the right one.
      while (sub_neighbor->has_children())
        {
          Assert ((GeometryInfo<dim>::face_refinement_case(sub_neighbor->refinement_case(),
                                                           neighbor_neighbor) ==
                   RefinementCase<dim>::no_refinement),
                  ExcInternalError());
          sub_neighbor = sub_neighbor->child(GeometryInfo<dim>::child_cell_on_face(
                                               sub_neighbor->refinement_case(),neighbor_neighbor,0));

        }

      return sub_neighbor;
    }


    case 3:
    {
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
      // |   |   |        |   1   |
      // | 0 | 1 |  or  *-------*
      // |   |   |        |   0   |
      // *---*---*        *-------*
      //
      // for "mixed" refinement:
      //
      // *---*---*      *---*---*      *---*---*      *-------*
      // |   | 2 |      | 1 |   |      | 1 | 2 |      |   2   |
      // | 0 *---*  or  *---* 2 |  or  *---*---*  or  *---*---*
      // |   | 1 |      | 0 |   |      |   0   |      | 0 | 1 |
      // *---*---*      *---*---*      *-------*      *---*---*

      const typename Triangulation<3,spacedim>::face_iterator
      mother_face = this->face(face);
      const unsigned int total_children=mother_face->number_of_children();
      Assert (subface<total_children,ExcIndexRange(subface,0,total_children));
      Assert (total_children<=GeometryInfo<3>::max_children_per_face, ExcInternalError());

      unsigned int neighbor_neighbor;
      TriaIterator<CellAccessor<3,spacedim> > neighbor_child;
      const TriaIterator<CellAccessor<3,spacedim> > neighbor
        = this->neighbor(face);


      const RefinementCase<2> mother_face_ref_case
        = mother_face->refinement_case();
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
          neighbor_child = neighbor->child(neighbor_child_index);

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
                = GeometryInfo<3>::child_cell_on_face(neighbor->refinement_case(),
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
                = GeometryInfo<3>::child_cell_on_face(neighbor->refinement_case(),
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
                 GeometryInfo<3>::face_refinement_case(neighbor_child->refinement_case(),
                                                       neighbor_neighbor)
                 == RefinementCase<2>::no_refinement)
            neighbor_child =
              neighbor_child->child(GeometryInfo<3>::
                                    child_cell_on_face(neighbor_child->refinement_case(),
                                                       neighbor_neighbor,
                                                       0));

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
                    neighbor_child =
                      neighbor_child->child(GeometryInfo<3>::
                                            child_cell_on_face(neighbor_child->refinement_case(),
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
                    neighbor_child =
                      neighbor_child->child(GeometryInfo<3>::
                                            child_cell_on_face(neighbor_child->refinement_case(),
                                                               neighbor_neighbor,subface-1,
                                                               neighbor_child->face_orientation(neighbor_neighbor),
                                                               neighbor_child->face_flip(neighbor_neighbor),
                                                               neighbor_child->face_rotation(neighbor_neighbor),
                                                               mother_face->child(1)->refinement_case()));
                }
            }
          else if (total_children==4)
            {
              neighbor_child =
                neighbor_child->child(GeometryInfo<3>::
                                      child_cell_on_face(neighbor_child->refinement_case(),
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
        neighbor_child
          = neighbor_child->child(GeometryInfo<3>::
                                  child_cell_on_face(neighbor_child->refinement_case(),
                                                     neighbor_neighbor,
                                                     0));

#ifdef DEBUG
      // check, whether the face neighbor_child
      // matches the requested subface
      typename Triangulation<3,spacedim>::face_iterator requested;
      switch (this->subface_case(face))
        {
        case internal::SubfaceCase<3>::case_x:
        case internal::SubfaceCase<3>::case_y:
        case internal::SubfaceCase<3>::case_xy:
          requested = mother_face->child(subface);
          break;
        case internal::SubfaceCase<3>::case_x1y2y:
        case internal::SubfaceCase<3>::case_y1x2x:
          requested = mother_face->child(subface/2)->child(subface%2);
          break;

        case internal::SubfaceCase<3>::case_x1y:
        case internal::SubfaceCase<3>::case_y1x:
          switch (subface)
            {
            case 0:
            case 1:
              requested = mother_face->child(0)->child(subface);
              break;
            case 2:
              requested = mother_face->child(1);
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
      Assert (requested==neighbor_child->face(neighbor_neighbor),
              ExcInternalError());
#endif

      return neighbor_child;

    }

    default:
      // 1d or more than 3d
      Assert (false, ExcNotImplemented());
      return TriaIterator<CellAccessor<dim,spacedim> >();
    }
}



// explicit instantiations
#include "tria_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
