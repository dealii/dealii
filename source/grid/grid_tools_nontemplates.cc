// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_tools.h>

#include <vector>

// GridTools functions that are template specializations (i.e., only compiled
// once without expand_instantiations)

DEAL_II_NAMESPACE_OPEN


namespace GridTools
{
  template <>
  double
  cell_measure<1>(
    const std::vector<Point<1>> &all_vertices,
    const unsigned int (&vertex_indices)[GeometryInfo<1>::vertices_per_cell])
  {
    return all_vertices[vertex_indices[1]][0] -
           all_vertices[vertex_indices[0]][0];
  }



  template <>
  double
  cell_measure<2>(
    const std::vector<Point<2>> &all_vertices,
    const unsigned int (&vertex_indices)[GeometryInfo<2>::vertices_per_cell])
  {
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
      z := array(0..3);
      tphi[0] := (1-xi)*(1-eta):
      tphi[1] :=     xi*(1-eta):
      tphi[2] := (1-xi)*eta:
      tphi[3] :=     xi*eta:
      x_real := sum(x[s]*tphi[s], s=0..3):
      y_real := sum(y[s]*tphi[s], s=0..3):
      z_real := sum(z[s]*tphi[s], s=0..3):

      Jxi := <diff(x_real,xi)  | diff(y_real,xi) | diff(z_real,xi)>;
      Jeta := <diff(x_real,eta)| diff(y_real,eta)| diff(z_real,eta)>;
      with(VectorCalculus):
      J := CrossProduct(Jxi, Jeta);
      detJ := sqrt(J[1]^2 + J[2]^2 +J[3]^2);

      # measure := evalf (Int (Int (detJ, xi=0..1, method = _NCrule ) ,
      eta=0..1, method = _NCrule  ) ): # readlib(C):

      # C(measure, optimized);

      additional optimizaton: divide by 2 only one time
    */

    const double x[4] = {all_vertices[vertex_indices[0]](0),
                         all_vertices[vertex_indices[1]](0),
                         all_vertices[vertex_indices[2]](0),
                         all_vertices[vertex_indices[3]](0)};

    const double y[4] = {all_vertices[vertex_indices[0]](1),
                         all_vertices[vertex_indices[1]](1),
                         all_vertices[vertex_indices[2]](1),
                         all_vertices[vertex_indices[3]](1)};

    return (-x[1] * y[0] + x[1] * y[3] + y[0] * x[2] + x[0] * y[1] -
            x[0] * y[2] - y[1] * x[3] - x[2] * y[3] + x[3] * y[2]) /
           2;
  }



  template <>
  double
  cell_measure<3>(
    const std::vector<Point<3>> &all_vertices,
    const unsigned int (&vertex_indices)[GeometryInfo<3>::vertices_per_cell])
  {
    // note that this is the
    // cell_measure based on the new
    // deal.II numbering. When called
    // from inside GridReordering make
    // sure that you reorder the
    // vertex_indices before
    const double x[8] = {all_vertices[vertex_indices[0]](0),
                         all_vertices[vertex_indices[1]](0),
                         all_vertices[vertex_indices[2]](0),
                         all_vertices[vertex_indices[3]](0),
                         all_vertices[vertex_indices[4]](0),
                         all_vertices[vertex_indices[5]](0),
                         all_vertices[vertex_indices[6]](0),
                         all_vertices[vertex_indices[7]](0)};
    const double y[8] = {all_vertices[vertex_indices[0]](1),
                         all_vertices[vertex_indices[1]](1),
                         all_vertices[vertex_indices[2]](1),
                         all_vertices[vertex_indices[3]](1),
                         all_vertices[vertex_indices[4]](1),
                         all_vertices[vertex_indices[5]](1),
                         all_vertices[vertex_indices[6]](1),
                         all_vertices[vertex_indices[7]](1)};
    const double z[8] = {all_vertices[vertex_indices[0]](2),
                         all_vertices[vertex_indices[1]](2),
                         all_vertices[vertex_indices[2]](2),
                         all_vertices[vertex_indices[3]](2),
                         all_vertices[vertex_indices[4]](2),
                         all_vertices[vertex_indices[5]](2),
                         all_vertices[vertex_indices[6]](2),
                         all_vertices[vertex_indices[7]](2)};

    /*
      This is the same Maple script as in the barycenter method above
      except of that here the shape functions tphi[0]-tphi[7] are ordered
      according to the lexicographic numbering.

      x := array(0..7):
      y := array(0..7):
      z := array(0..7):
      tphi[0] := (1-xi)*(1-eta)*(1-zeta):
      tphi[1] :=     xi*(1-eta)*(1-zeta):
      tphi[2] := (1-xi)*    eta*(1-zeta):
      tphi[3] :=     xi*    eta*(1-zeta):
      tphi[4] := (1-xi)*(1-eta)*zeta:
      tphi[5] :=     xi*(1-eta)*zeta:
      tphi[6] := (1-xi)*    eta*zeta:
      tphi[7] :=     xi*    eta*zeta:
      x_real := sum(x[s]*tphi[s], s=0..7):
      y_real := sum(y[s]*tphi[s], s=0..7):
      z_real := sum(z[s]*tphi[s], s=0..7):
      with (linalg):
      J := matrix(3,3, [[diff(x_real, xi), diff(x_real, eta), diff(x_real,
      zeta)], [diff(y_real, xi), diff(y_real, eta), diff(y_real, zeta)],
      [diff(z_real, xi), diff(z_real, eta), diff(z_real, zeta)]]):
      detJ := det (J):

      measure := simplify ( int ( int ( int (detJ, xi=0..1), eta=0..1),
      zeta=0..1)):

      readlib(C):

      C(measure, optimized);

      The C code produced by this maple script is further optimized by
      hand. In particular, division by 12 is performed only once, not
      hundred of times.
    */

    const double t3  = y[3] * x[2];
    const double t5  = z[1] * x[5];
    const double t9  = z[3] * x[2];
    const double t11 = x[1] * y[0];
    const double t14 = x[4] * y[0];
    const double t18 = x[5] * y[7];
    const double t20 = y[1] * x[3];
    const double t22 = y[5] * x[4];
    const double t26 = z[7] * x[6];
    const double t28 = x[0] * y[4];
    const double t34 =
      z[3] * x[1] * y[2] + t3 * z[1] - t5 * y[7] + y[7] * x[4] * z[6] +
      t9 * y[6] - t11 * z[4] - t5 * y[3] - t14 * z[2] + z[1] * x[4] * y[0] -
      t18 * z[3] + t20 * z[0] - t22 * z[0] - y[0] * x[5] * z[4] - t26 * y[3] +
      t28 * z[2] - t9 * y[1] - y[1] * x[4] * z[0] - t11 * z[5];
    const double t37 = y[1] * x[0];
    const double t44 = x[1] * y[5];
    const double t46 = z[1] * x[0];
    const double t49 = x[0] * y[2];
    const double t52 = y[5] * x[7];
    const double t54 = x[3] * y[7];
    const double t56 = x[2] * z[0];
    const double t58 = x[3] * y[2];
    const double t64 = -x[6] * y[4] * z[2] - t37 * z[2] + t18 * z[6] -
                       x[3] * y[6] * z[2] + t11 * z[2] + t5 * y[0] +
                       t44 * z[4] - t46 * y[4] - t20 * z[7] - t49 * z[6] -
                       t22 * z[1] + t52 * z[3] - t54 * z[2] - t56 * y[4] -
                       t58 * z[0] + y[1] * x[2] * z[0] + t9 * y[7] + t37 * z[4];
    const double t66 = x[1] * y[7];
    const double t68 = y[0] * x[6];
    const double t70 = x[7] * y[6];
    const double t73 = z[5] * x[4];
    const double t76 = x[6] * y[7];
    const double t90 = x[4] * z[0];
    const double t92 = x[1] * y[3];
    const double t95 = -t66 * z[3] - t68 * z[2] - t70 * z[2] + t26 * y[5] -
                       t73 * y[6] - t14 * z[6] + t76 * z[2] - t3 * z[6] +
                       x[6] * y[2] * z[4] - z[3] * x[6] * y[2] + t26 * y[4] -
                       t44 * z[3] - x[1] * y[2] * z[0] + x[5] * y[6] * z[4] +
                       t54 * z[5] + t90 * y[2] - t92 * z[2] + t46 * y[2];
    const double t102 = x[2] * y[0];
    const double t107 = y[3] * x[7];
    const double t114 = x[0] * y[6];
    const double t125 =
      y[0] * x[3] * z[2] - z[7] * x[5] * y[6] - x[2] * y[6] * z[4] +
      t102 * z[6] - t52 * z[6] + x[2] * y[4] * z[6] - t107 * z[5] - t54 * z[6] +
      t58 * z[6] - x[7] * y[4] * z[6] + t37 * z[5] - t114 * z[4] + t102 * z[4] -
      z[1] * x[2] * y[0] + t28 * z[6] - y[5] * x[6] * z[4] -
      z[5] * x[1] * y[4] - t73 * y[7];
    const double t129 = z[0] * x[6];
    const double t133 = y[1] * x[7];
    const double t145 = y[1] * x[5];
    const double t156 = t90 * y[6] - t129 * y[4] + z[7] * x[2] * y[6] -
                        t133 * z[5] + x[5] * y[3] * z[7] - t26 * y[2] -
                        t70 * z[3] + t46 * y[3] + z[5] * x[7] * y[4] +
                        z[7] * x[3] * y[6] - t49 * z[4] + t145 * z[7] -
                        x[2] * y[7] * z[6] + t70 * z[5] + t66 * z[5] -
                        z[7] * x[4] * y[6] + t18 * z[4] + x[1] * y[4] * z[0];
    const double t160 = x[5] * y[4];
    const double t165 = z[1] * x[7];
    const double t178 = z[1] * x[3];
    const double t181 =
      t107 * z[6] + t22 * z[7] + t76 * z[3] + t160 * z[1] - x[4] * y[2] * z[6] +
      t70 * z[4] + t165 * y[5] + x[7] * y[2] * z[6] - t76 * z[5] - t76 * z[4] +
      t133 * z[3] - t58 * z[1] + y[5] * x[0] * z[4] + t114 * z[2] - t3 * z[7] +
      t20 * z[2] + t178 * y[7] + t129 * y[2];
    const double t207 = t92 * z[7] + t22 * z[6] + z[3] * x[0] * y[2] -
                        x[0] * y[3] * z[2] - z[3] * x[7] * y[2] - t165 * y[3] -
                        t9 * y[0] + t58 * z[7] + y[3] * x[6] * z[2] +
                        t107 * z[2] + t73 * y[0] - x[3] * y[5] * z[7] +
                        t3 * z[0] - t56 * y[6] - z[5] * x[0] * y[4] +
                        t73 * y[1] - t160 * z[6] + t160 * z[0];
    const double t228 = -t44 * z[7] + z[5] * x[6] * y[4] - t52 * z[4] -
                        t145 * z[4] + t68 * z[4] + t92 * z[5] - t92 * z[0] +
                        t11 * z[3] + t44 * z[0] + t178 * y[5] - t46 * y[5] -
                        t178 * y[0] - t145 * z[0] - t20 * z[5] - t37 * z[3] -
                        t160 * z[7] + t145 * z[3] + x[4] * y[6] * z[2];

    return (t34 + t64 + t95 + t125 + t156 + t181 + t207 + t228) / 12.;
  }



  namespace
  {
    // the following class is only
    // needed in 2d, so avoid trouble
    // with compilers warning otherwise
    class Rotate2d
    {
    public:
      explicit Rotate2d(const double angle)
        : angle(angle)
      {}
      Point<2>
      operator()(const Point<2> &p) const
      {
        return {std::cos(angle) * p(0) - std::sin(angle) * p(1),
                std::sin(angle) * p(0) + std::cos(angle) * p(1)};
      }

    private:
      const double angle;
    };
  } // namespace



  template <>
  void
  rotate(const double angle, Triangulation<2> &triangulation)
  {
    transform(Rotate2d(angle), triangulation);
  }



  template <>
  void
  rotate(const double angle, Triangulation<3> &triangulation)
  {
    (void)angle;
    (void)triangulation;

    AssertThrow(
      false, ExcMessage("GridTools::rotate() is not available for dim = 3."));
  }
} /* namespace GridTools */

DEAL_II_NAMESPACE_CLOSE
