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
  cell_measure<1>(const std::vector<Point<1>> &        all_vertices,
                  const ArrayView<const unsigned int> &vertex_indices)
  {
    AssertDimension(vertex_indices.size(), GeometryInfo<1>::vertices_per_cell);

    return all_vertices[vertex_indices[1]][0] -
           all_vertices[vertex_indices[0]][0];
  }



  template <>
  double
  cell_measure<2>(const std::vector<Point<2>> &        all_vertices,
                  const ArrayView<const unsigned int> &vertex_indices)
  {
    if (vertex_indices.size() == 3) // triangle
      {
        const double x[3] = {all_vertices[vertex_indices[0]](0),
                             all_vertices[vertex_indices[1]](0),
                             all_vertices[vertex_indices[2]](0)};

        const double y[3] = {all_vertices[vertex_indices[0]](1),
                             all_vertices[vertex_indices[1]](1),
                             all_vertices[vertex_indices[2]](1)};

        return 0.5 * std::abs((x[0] - x[2]) * (y[1] - y[0]) -
                              (x[0] - x[1]) * (y[2] - y[0]));
      }

    AssertDimension(vertex_indices.size(), GeometryInfo<2>::vertices_per_cell);

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
  cell_measure<3>(const std::vector<Point<3>> &        all_vertices,
                  const ArrayView<const unsigned int> &vertex_indices)
  {
    if (vertex_indices.size() == 4) // tetrahedron
      {
        const auto &a = all_vertices[vertex_indices[0]];
        const auto &b = all_vertices[vertex_indices[1]];
        const auto &c = all_vertices[vertex_indices[2]];
        const auto &d = all_vertices[vertex_indices[3]];

        return (1.0 / 6.0) * std::abs((a - d) * cross_product_3d(b - d, c - d));
      }
    else if (vertex_indices.size() == 6) // wedge
      {
        /*
         The following python/sympy script was used:

         #!/usr/bin/env python
         # coding: utf-8
         import sympy as sp
         from sympy.simplify.cse_main import cse
         xs = list(sp.symbols(" ".join(["x{}".format(i) for i in range(6)])))
         ys = list(sp.symbols(" ".join(["y{}".format(i) for i in range(6)])))
         zs = list(sp.symbols(" ".join(["z{}".format(i) for i in range(6)])))
         xi, eta, zeta = sp.symbols("xi eta zeta")
         tphi = [(1 - xi - eta)*(1 - zeta),
                 (xi)*(1 - zeta),
                 (eta)*(1 - zeta),
                 (1 - xi - zeta)*(zeta),
                 (xi)*(zeta),
                 (eta)*(zeta)]
         x_real = sum(xs[i]*tphi[i] for i in range(len(xs)))
         y_real = sum(ys[i]*tphi[i] for i in range(len(xs)))
         z_real = sum(zs[i]*tphi[i] for i in range(len(xs)))
         J = sp.Matrix([[var.diff(v) for v in [xi, eta, zeta]]
                        for var in [x_real, y_real, z_real]])
         detJ = J.det()
         detJ2 = detJ.expand().collect(zeta).collect(eta).collect(xi)
         for x in xs:
             detJ2 = detJ2.collect(x)
         for y in ys:
             detJ2 = detJ2.collect(y)
         for z in zs:
             detJ2 = detJ2.collect(z)
         measure = sp.integrate(sp.integrate(sp.integrate(detJ2, (xi, 0, 1)),
                           (eta, 0, 1)), (zeta, 0, 1))
         measure2 = measure
         for vs in [xs, ys, zs]:
             for v in vs:
                 measure2 = measure2.collect(v)
         pairs, expression = cse(measure2)
         for pair in pairs:
             print("const double " + sp.ccode(pair[0]) + " = "
                                   + sp.ccode(pair[1]) + ";")
         print("const double result = " + sp.ccode(expression[0]) + ";")
         print("return result;")
         */
        const double x0 = all_vertices[vertex_indices[0]](0);
        const double y0 = all_vertices[vertex_indices[0]](1);
        const double z0 = all_vertices[vertex_indices[0]](2);
        const double x1 = all_vertices[vertex_indices[1]](0);
        const double y1 = all_vertices[vertex_indices[1]](1);
        const double z1 = all_vertices[vertex_indices[1]](2);
        const double x2 = all_vertices[vertex_indices[2]](0);
        const double y2 = all_vertices[vertex_indices[2]](1);
        const double z2 = all_vertices[vertex_indices[2]](2);
        const double x3 = all_vertices[vertex_indices[3]](0);
        const double y3 = all_vertices[vertex_indices[3]](1);
        const double z3 = all_vertices[vertex_indices[3]](2);
        const double x4 = all_vertices[vertex_indices[4]](0);
        const double y4 = all_vertices[vertex_indices[4]](1);
        const double z4 = all_vertices[vertex_indices[4]](2);
        const double x5 = all_vertices[vertex_indices[5]](0);
        const double y5 = all_vertices[vertex_indices[5]](1);
        const double z5 = all_vertices[vertex_indices[5]](2);

        const double x6  = (1.0 / 6.0) * z5;
        const double x7  = (1.0 / 12.0) * z1;
        const double x8  = -x7;
        const double x9  = (1.0 / 12.0) * z2;
        const double x10 = -x9;
        const double x11 = (1.0 / 4.0) * z5;
        const double x12 = -x11;
        const double x13 = (1.0 / 12.0) * z0;
        const double x14 = x12 + x13;
        const double x15 = (1.0 / 4.0) * z2;
        const double x16 = (1.0 / 6.0) * z4;
        const double x17 = (1.0 / 4.0) * z1;
        const double x18 = (1.0 / 6.0) * z0;
        const double x19 = x17 - x18;
        const double x20 = -x6;
        const double x21 = (1.0 / 4.0) * z0;
        const double x22 = -x21;
        const double x23 = -x17;
        const double x24 = -x15;
        const double x25 = (1.0 / 6.0) * z3;
        const double x26 = x24 - x25;
        const double x27 = x18 + x23;
        const double x28 = (1.0 / 3.0) * z2;
        const double x29 = (1.0 / 12.0) * z5;
        const double x30 = (1.0 / 12.0) * z3;
        const double x31 = -x30;
        const double x32 = (1.0 / 4.0) * z4;
        const double x33 = x31 + x32;
        const double x34 = (1.0 / 3.0) * z1;
        const double x35 = (1.0 / 12.0) * z4;
        const double x36 = -x16;
        const double x37 = x15 + x25;
        const double x38 = -x13;
        const double x39 = x11 + x38;
        const double x40 = -x32;
        const double x41 = x30 + x40;
        const double x42 = (1.0 / 3.0) * z0;
        const double x43 = (1.0 / 4.0) * z3;
        const double x44 = x32 - x43;
        const double x45 = x40 + x43;
        return x0 * (y1 * (-x28 + x29 + x33) + y2 * (x12 + x31 + x34 - x35) +
                     y3 * (x20 + x7 + x9) + y4 * (x23 + x6 + x9) +
                     y5 * (x36 + x37 + x8)) +
               x1 * (y0 * (x28 - x29 + x41) + y2 * (x11 + x33 - x42) +
                     y3 * (x39 + x9) + y4 * (x12 + x21 + x24) +
                     y5 * (x13 + x24 + x44)) +
               x2 * (y0 * (x11 + x30 - x34 + x35) + y1 * (x12 + x41 + x42) +
                     y3 * (x39 + x8) + y4 * (x12 + x17 + x38) +
                     y5 * (x17 + x22 + x44)) +
               x3 * (-x6 * y4 + y0 * (x10 + x6 + x8) + y1 * (x10 + x14) +
                     y2 * (x14 + x7) + y5 * (x15 + x16 + x19)) +
               x4 * (x6 * y3 + y0 * (x10 + x17 + x20) + y1 * (x11 + x15 + x22) +
                     y2 * (x11 + x13 + x23) + y5 * (x26 + x27)) +
               x5 * (y0 * (x16 + x26 + x7) + y1 * (x15 + x38 + x45) +
                     y2 * (x21 + x23 + x45) + y3 * (x24 + x27 + x36) +
                     y4 * (x19 + x37));
      }

    AssertDimension(vertex_indices.size(), GeometryInfo<3>::vertices_per_cell);

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
