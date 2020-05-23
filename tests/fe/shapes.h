// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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


// Show the shape functions implemented.

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"

char fname[50];

////////////////////////////////////////////////////////////////////////////
// Plot shape function values at quadrature points inside the cell [0,1]^d
//
// Output values in each line are
//
// x (y) (z) value[0]+1 value[1]+1 ...
////////////////////////////////////////////////////////////////////////////
template <int dim>
inline void
plot_shape_functions(Mapping<dim> &      mapping,
                     FiniteElement<dim> &finel,
                     const char *        name)
{
  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  dof.distribute_dofs(finel);
  typename DoFHandler<dim>::cell_iterator c = dof.begin();

  const unsigned int div = 4;

  QTrapez<1>     q_trapez;
  QIterated<dim> q(q_trapez, div);
  FEValues<dim>  fe(mapping,
                   finel,
                   q,
                   UpdateFlags(update_values | update_gradients |
                               update_hessians));

  sprintf(fname, "Cell%dd-%s", dim, name);
  //  cerr << "\n" << fname << "\n";
  deallog.push(fname);

  fe.reinit(c);

  unsigned int k = 0;
  for (unsigned int mz = 0; mz <= ((dim > 2) ? div : 0); ++mz)
    {
      for (unsigned int my = 0; my <= ((dim > 1) ? div : 0); ++my)
        {
          for (unsigned int mx = 0; mx <= div; ++mx)
            {
              deallog << q.point(k);

              for (unsigned int i = 0; i < finel.dofs_per_cell; ++i)
                {
                  deallog << " " << fe.shape_value(i, k) + 1.;

                  // some additional
                  // checks
                  for (unsigned int c = 0; c < fe.get_fe().n_components(); ++c)
                    {
                      if (fe.get_fe().system_to_component_index(i).first == c)
                        {
                          AssertThrow((fe.shape_value(i, k) ==
                                         fe.shape_value_component(i, k, c) &&
                                       fe.shape_grad(i, k) ==
                                         fe.shape_grad_component(i, k, c) &&
                                       fe.shape_hessian(i, k) ==
                                         fe.shape_hessian_component(i, k, c)),
                                      ExcInternalError());
                        }
                      else
                        {
                          AssertThrow((fe.shape_value_component(i, k, c) == 0 &&
                                       fe.shape_grad_component(i, k, c) ==
                                         Tensor<1, dim>() &&
                                       fe.shape_hessian_component(i, k, c) ==
                                         Tensor<2, dim>()),
                                      ExcInternalError());
                        }
                    };
                }
              deallog << std::endl;
              k++;
            }
          deallog << std::endl;
        }
      deallog << std::endl;
    }
  deallog.pop();
}



template <int dim>
inline void
plot_face_shape_functions(Mapping<dim> &      mapping,
                          FiniteElement<dim> &finel,
                          const char *        name,
                          UpdateFlags uflags = UpdateFlags(update_values |
                                                           update_gradients |
                                                           update_hessians))
{
  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  tr.refine_global(1);
  typename DoFHandler<dim>::active_cell_iterator c = dof.begin_active();
  ++c;
  c->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  c = dof.begin_active();

  dof.distribute_dofs(finel);

  const unsigned int div = 4;

  QTrapez<1>           q_trapez;
  QIterated<dim - 1>   q(q_trapez, div);
  FEFaceValues<dim>    fe(mapping,
                       finel,
                       q,
                       UpdateFlags(uflags | update_quadrature_points));
  FESubfaceValues<dim> sub(mapping,
                           finel,
                           q,
                           UpdateFlags(uflags | update_quadrature_points));

  // Test the iterator-based initialization routines too:
  FEFaceValues<dim>    fe2(mapping,
                        finel,
                        q,
                        UpdateFlags(uflags | update_quadrature_points));
  FESubfaceValues<dim> sub2(mapping,
                            finel,
                            q,
                            UpdateFlags(uflags | update_quadrature_points));


  sprintf(fname, "Face%dd-%s", dim, name);
  deallog.push(fname);

  for (const unsigned int f : GeometryInfo<dim>::face_indices())
    {
      if (!c->face(f)->has_children())
        {
          fe.reinit(c, f);
          fe2.reinit(c, c->face(f));

          unsigned int k = 0;
          for (unsigned int my = 0; my <= ((dim > 2) ? div : 0); ++my)
            {
              for (unsigned int mx = 0; mx <= div; ++mx)
                {
                  deallog << fe.quadrature_point(k);

                  for (unsigned int i = 0; i < finel.dofs_per_cell; ++i)
                    {
                      if (uflags & update_values)
                        deallog << " " << fe.shape_value(i, k) + 1.;

                      // some additional
                      // checks
                      for (unsigned int c = 0; c < fe.get_fe().n_components();
                           ++c)
                        {
                          if (fe.get_fe().system_to_component_index(i).first ==
                              c)
                            {
                              if (uflags & update_values)
                                {
                                  AssertThrow(
                                    (fe.shape_value(i, k) ==
                                     fe.shape_value_component(i, k, c)),
                                    ExcInternalError());
                                  AssertThrow((fe.shape_value(i, k) ==
                                               fe2.shape_value(i, k)),
                                              ExcInternalError());
                                }
                              if (uflags & update_gradients)
                                {
                                  AssertThrow(
                                    (fe.shape_grad(i, k) ==
                                     fe.shape_grad_component(i, k, c)),
                                    ExcInternalError());
                                  AssertThrow((fe.shape_grad(i, k) ==
                                               fe2.shape_grad(i, k)),
                                              ExcInternalError());
                                }
                              if (uflags & update_hessians)
                                {
                                  AssertThrow(
                                    (fe.shape_hessian(i, k) ==
                                     fe.shape_hessian_component(i, k, c)),
                                    ExcInternalError());
                                  AssertThrow((fe.shape_hessian(i, k) ==
                                               fe2.shape_hessian(i, k)),
                                              ExcInternalError());
                                }
                              if (uflags & update_3rd_derivatives)
                                {
                                  AssertThrow(
                                    (fe.shape_3rd_derivative(i, k) ==
                                     fe.shape_3rd_derivative_component(i,
                                                                       k,
                                                                       c)),
                                    ExcInternalError());
                                  AssertThrow((fe.shape_3rd_derivative(i, k) ==
                                               fe2.shape_3rd_derivative(i, k)),
                                              ExcInternalError());
                                }
                            }
                          else
                            {
                              if (uflags & update_values)
                                AssertThrow(
                                  (fe.shape_value_component(i, k, c) == 0),
                                  ExcInternalError());
                              if (uflags & update_gradients)
                                AssertThrow((fe.shape_grad_component(i, k, c) ==
                                             Tensor<1, dim>()),
                                            ExcInternalError());
                              if (uflags & update_hessians)
                                AssertThrow((fe.shape_hessian_component(
                                               i, k, c) == Tensor<2, dim>()),
                                            ExcInternalError());
                              if (uflags & update_3rd_derivatives)
                                AssertThrow((fe.shape_3rd_derivative_component(
                                               i, k, c) == Tensor<3, dim>()),
                                            ExcInternalError());
                            }
                        }
                    }
                  deallog << std::endl;
                  k++;
                }
              deallog << std::endl;
            }
          deallog << std::endl;
        }
      else
        {
          for (unsigned int s = 0; s < GeometryInfo<dim>::max_children_per_face;
               ++s)
            {
              sub.reinit(c, f, s);
              sub2.reinit(c, c->face(f), c->face(f)->child(s));

              unsigned int k = 0;
              for (unsigned int my = 0; my <= ((dim > 2) ? div : 0); ++my)
                {
                  for (unsigned int mx = 0; mx <= div; ++mx)
                    {
                      deallog << sub.quadrature_point(k);

                      for (unsigned int i = 0; i < finel.dofs_per_cell; ++i)
                        {
                          if (uflags & update_values)
                            deallog << " " << sub.shape_value(i, k) + 1.;

                          // some additional
                          // checks
                          for (unsigned int c = 0;
                               c < fe.get_fe().n_components();
                               ++c)
                            {
                              if (fe.get_fe()
                                    .system_to_component_index(i)
                                    .first == c)
                                {
                                  if (uflags & update_values)
                                    {
                                      const double v1 = sub.shape_value(i, k),
                                                   v2 =
                                                     sub.shape_value_component(
                                                       i, k, c),
                                                   v3 = sub2.shape_value(i, k);
                                      Assert(v1 == v2, ExcInternalError());
                                      Assert(v1 == v3, ExcInternalError());
                                    }
                                  if (uflags & update_gradients)
                                    {
                                      const Tensor<1, dim>
                                        g1 = sub.shape_grad(i, k),
                                        g2 = sub.shape_grad_component(i, k, c),
                                        g3 = sub2.shape_grad(i, k);
                                      Assert(g1 == g2, ExcInternalError());
                                      Assert(g1 == g3, ExcInternalError());
                                    }
                                  if (uflags & update_hessians)
                                    {
                                      const Tensor<2, dim>
                                        s1 = sub.shape_hessian(i, k),
                                        s2 =
                                          sub.shape_hessian_component(i, k, c),
                                        s3 = sub2.shape_hessian(i, k);
                                      Assert(s1 == s2, ExcInternalError());
                                      Assert(s1 == s3, ExcInternalError());
                                    }
                                  if (uflags & update_3rd_derivatives)
                                    {
                                      const Tensor<3, dim>
                                        t1 = sub.shape_3rd_derivative(i, k),
                                        t2 =
                                          sub.shape_3rd_derivative_component(i,
                                                                             k,
                                                                             c),
                                        t3 = sub2.shape_3rd_derivative(i, k);

                                      Assert(t1 == t2, ExcInternalError());
                                      Assert(t1 == t3, ExcInternalError());
                                    }
                                }
                              else
                                {
                                  if (uflags & update_values)
                                    Assert(
                                      (sub.shape_value_component(i, k, c) == 0),
                                      ExcInternalError());
                                  if (uflags & update_gradients)
                                    Assert((sub.shape_grad_component(i, k, c) ==
                                            Tensor<1, dim>()),
                                           ExcInternalError());
                                  if (uflags & update_hessians)
                                    Assert((sub.shape_hessian_component(
                                              i, k, c) == Tensor<2, dim>()),
                                           ExcInternalError());
                                  if (uflags & update_3rd_derivatives)
                                    Assert((sub.shape_3rd_derivative_component(
                                              i, k, c) == Tensor<3, dim>()),
                                           ExcInternalError());
                                }
                            };
                        }
                      deallog << std::endl;
                      k++;
                    }
                  deallog << std::endl;
                }
              deallog << std::endl;
            }
        }
    }
  deallog.pop();
}


template <>
void plot_face_shape_functions(Mapping<1> &,
                               FiniteElement<1> &,
                               const char *,
                               UpdateFlags)
{}



// given an FEValues object for a cell that is equal to the unit cell,
// check that the values and gradients that the FEValues object
// generated are equivalent to the values the finite element returns
// for the unit cell
template <int dim>
void
check_values_and_derivatives(const FiniteElement<dim> &fe,
                             const FEValuesBase<dim> & fe_values,
                             const Quadrature<dim> &   q)
{
  // check values
  for (unsigned int x = 0; x < q.size(); ++x)
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        for (unsigned int c = 0; c < fe.n_components(); ++c)
          {
            const double val1 = fe_values.shape_value_component(i, x, c);
            const double val2 = fe.shape_value_component(i, q.point(x), c);
            const double diff = std::fabs(val1 - val2);
            if (diff > 1e-13)
              deallog << " values differ v" << i << "(x" << x << ") diff "
                      << diff << std::endl;
          };

        // test something about the
        // correctness of indices
        // etc, except for the more
        // complicated case of
        // non-primitive elements
        if (fe.is_primitive(i))
          for (unsigned int c = 0; c < fe.n_components(); ++c)
            Assert(((c == fe.system_to_component_index(i).first) &&
                    (fe_values.shape_value(i, x) ==
                     fe_values.shape_value_component(i, x, c))) ||
                     ((c != fe.system_to_component_index(i).first) &&
                      (fe_values.shape_value_component(i, x, c) == 0)),
                   ExcInternalError());
      };

  // check gradients
  for (unsigned int x = 0; x < q.size(); ++x)
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        for (unsigned int c = 0; c < fe.n_components(); ++c)
          {
            Tensor<1, dim> tmp = fe_values.shape_grad_component(i, x, c);
            tmp -= fe.shape_grad_component(i, q.point(x), c);
            Assert(std::sqrt(tmp * tmp) < 1e-14, ExcInternalError());
          };

        if (fe.is_primitive(i))
          for (unsigned int c = 0; c < fe.n_components(); ++c)
            Assert(((c == fe.system_to_component_index(i).first) &&
                    (fe_values.shape_grad(i, x) ==
                     fe_values.shape_grad_component(i, x, c))) ||
                     ((c != fe.system_to_component_index(i).first) &&
                      (fe_values.shape_grad_component(i, x, c) ==
                       Tensor<1, dim>())),
                   ExcInternalError());
      }

  // check second derivatives
  double max_diff = 0.;
  for (unsigned int x = 0; x < q.size(); ++x)
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      {
        for (unsigned int c = 0; c < fe.n_components(); ++c)
          {
            Tensor<2, dim> tmp = fe_values.shape_hessian_component(i, x, c);
            tmp -= fe.shape_grad_grad_component(i, q.point(x), c);
            for (unsigned int j = 0; j < dim; ++j)
              for (unsigned int k = 0; k < dim; ++k)
                {
                  const double diff = std::fabs(tmp[j][k]);
                  if (diff > max_diff)
                    max_diff = diff;
                  const double tmpabs = std::fabs(tmp[j][k]);
                  if (tmpabs > 1.e-6)
                    deallog << "Second derivatives differ " << tmpabs
                            << std::endl;
                }
          };

        if (fe.is_primitive(i))
          for (unsigned int c = 0; c < fe.n_components(); ++c)
            Assert(((c == fe.system_to_component_index(i).first) &&
                    (fe_values.shape_hessian(i, x) ==
                     fe_values.shape_hessian_component(i, x, c))) ||
                     ((c != fe.system_to_component_index(i).first) &&
                      (fe_values.shape_hessian_component(i, x, c) ==
                       Tensor<2, dim>())),
                   ExcInternalError());
      }
}



template <int dim>
void
test_compute_functions(const Mapping<dim> &      mapping,
                       const FiniteElement<dim> &fe,
                       const char *)
{
  // generate a grid with only one
  // cell, which furthermore has the
  // shape of the unit cell. then the
  // values/gradients/... we get from
  // the FEValues object on this cell
  // should really be equal to what
  // we get from the finite element
  // itself on the unit cell:
  Triangulation<dim> tr;
  DoFHandler<dim>    dof(tr);
  GridGenerator::hyper_cube(tr, 0., 1.);
  dof.distribute_dofs(fe);

  const UpdateFlags update_all =
    (update_values | update_gradients | update_hessians);

  // first check this for FEValues
  // objects
  if (true)
    {
      const QGauss<dim> q(6);
      FEValues<dim>     fe_values(mapping, fe, q, update_all);
      fe_values.reinit(dof.begin_active());
      check_values_and_derivatives(fe, fe_values, q);
    };
}
