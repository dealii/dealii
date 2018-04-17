// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2017 by the deal.II authors
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



// add test for ExcNotImplemented() error in compute_no_normal_flux_constraints.
// reported by Keith Galvin, mailing list, 2013/10/13. Simplified.
// note that the cylinder boundary is not even the problem here!

/*
5: An error occurred in line <4586> of file </scratch/deal-trunk/deal.II/include/deal.II/numerics/vector_tools.templates.h> in function
5:     void dealii::VectorTools::compute_no_normal_flux_constraints(const DoFHandlerType<dim, spacedim>&, unsigned int, const std::set<types::boundary_id>&, dealii::ConstraintMatrix&, const dealii::Mapping<dim, spacedim>&) [with int dim = 3; DoFHandlerType = dealii::DoFHandler; int spacedim = 3]
5: The violated condition was:
5:     contribution->second.size() == dim-1
5: The name and call sequence of the exception was:
5:     ExcNotImplemented()

 */

// The problem comes down to this: The domain (a box with a cylinder
// through it) is discretized with eight cells that go around the
// perimeter of the cylinder and are extruded to full height in
// z-direction. In particular, for the top and bottom surfaces, there
// is a cell edge from the circle to the corners of the top and bottom
// outer squares (as well as an edge from the circle to the midpoints
// of the sides of the squares).
//
// now, we try to compute the boundary conditions for no normal flux
// with a set of boundary indicators that includes the top part of the
// domain plus one of the four outer sides of the box. consider what
// happens at one of the vertices of the box that are part of both the
// top and the selected side: there, one of the eight cells of the
// mesh contributes two normal vectors (one from its top face and from
// its "side" face) but since at each of the eight vertices of the box
// two cells come together, there is also a cell that contributes only
// one normal vector (because it has one face at the top surface, but
// its other outer face does not have a selected boundary
// indicator). we used to be unable to deal with this situation where
// one cell contributed twice and another only once. this has been
// fixed now by simply ignoring the one that contributes only once.

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>


template <int dim>
void run()
{
  Triangulation<dim> tria;

  // indicator 6 = cylinder
  GridGenerator::hyper_cube_with_cylindrical_hole (tria, 0.25, 0.5, 0.5, 1, true);
  tria.reset_manifold(0);

  /*  std::string filename = "Mesh.eps";
  std::ofstream output (filename.c_str());
  GridOut grid_out;
  grid_out.write_eps (tria, output);
  */

  FESystem<dim> fe(FE_Q<dim>(1), dim);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  ConstraintMatrix constraints;
  std::set<types::boundary_id> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0); // x=0
  no_normal_flux_boundaries.insert (5); // z=1

  VectorTools::compute_no_normal_flux_constraints
  (dof_handler, 0,
   no_normal_flux_boundaries,
   constraints);

  constraints.print(deallog.get_file_stream());

  deallog.get_file_stream() << std::flush;
  constraints.close();

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile ("output");
  logfile.precision (7);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);

  run<3> ();
}
