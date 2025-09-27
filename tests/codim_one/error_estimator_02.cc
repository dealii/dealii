// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check that the kelly error estimator uses correct normals on surfaces
// with kinks by interpolating a function that is inside the FE space
// and should produce zero errors.



#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"



// function is x^2+2xy
template <int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    return p[0] * p[0] + 2.0 * p[0] * p[1];
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    values(0) = value(p, 0);
  }
};


// normal derivative of function above is 0 except for x=1, y=1 with n=(0,1,0)
template <int dim>
class MyNormalDerivative : public Function<dim>
{
public:
  MyNormalDerivative()
    : Function<dim>(1)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const
  {
    double val = 0.0;
    if (std::abs(p[1] - 1.0) < 1e-5)
      val = 2.0;

    deallog << "evaluate normal derivative at " << p << " with value " << val
            << std::endl;
    return val;
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const
  {
    values(0) = value(p, 0);
  }
};



template <int dim>
Quadrature<dim - 1> &
get_q_face()
{
  static QGauss<dim - 1> q(4);
  return q;
}



template <int dim, int spacedim>
void
make_mesh(Triangulation<dim, spacedim> &tria)
{
  // two faces of a hyper_cube
  Triangulation<spacedim, spacedim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh, 0, 1, true);
  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert(1);
  boundary_ids.insert(2);
  GridGenerator::extract_boundary_mesh(volume_mesh, tria, boundary_ids);
  tria.refine_global(2);
  typename Triangulation<dim, spacedim>::active_cell_iterator cell =
    tria.begin_active();
  for (; cell != tria.end(); ++cell)
    if (cell->center()[1] < 1e-10 && (spacedim != 3 || cell->center()[2] < 0.5))
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
}



template <int dim, int spacedim>
void
check()
{
  MyFunction<spacedim> function;

  Triangulation<dim, spacedim> tria;
  make_mesh(tria);

  FE_Q<dim, spacedim>       element(2);
  DoFHandler<dim, spacedim> dof(tria);
  dof.distribute_dofs(element);

  MappingQ<dim, spacedim> mapping(3);
  Quadrature<dim - 1>    &q_face = get_q_face<dim>();

  Vector<double> v(dof.n_dofs());
  VectorTools::interpolate(mapping, dof, function, v);

  Vector<float> error(tria.n_active_cells());

  std::map<types::boundary_id, const Function<spacedim> *> neumann_bc;
  MyNormalDerivative<spacedim>                             function_normal;
  neumann_bc[0] = &function_normal;
  neumann_bc[1] = &function_normal;

  deallog << "estimating..." << std::endl;
  KellyErrorEstimator<dim, spacedim>::estimate(
    mapping, dof, q_face, neumann_bc, v, error);
  deallog << "Estimated error indicators:" << std::endl;
  for (unsigned int i = 0; i < error.size(); ++i)
    deallog << error(i) << std::endl;

  {
    DataOut<dim, spacedim> data_out;
    data_out.attach_dof_handler(dof);
    data_out.add_data_vector(v,
                             "solution",
                             DataOut<dim, spacedim>::type_dof_data);
    data_out.add_data_vector(error, "error");
    data_out.build_patches();
    std::string filename = spacedim == 2 ? "solution-2d-" : "solution-3d-";
    filename += Utilities::int_to_string(0, 2) + ".vtk";
    std::ofstream output(filename);
    data_out.write_vtk(output);
  }

  deallog << "OK" << std::endl;
}


int
main()
{
  // do not run multithreaded estimate() because that would break the order of
  // the function evaluations that we log above
  MultithreadInfo::set_thread_limit(1);

  initlog(true);

  deallog.push("<1,2>");
  check<1, 2>();
  deallog.pop();
  deallog.push("<2,3>");
  check<2, 3>();
  deallog.pop();
}
