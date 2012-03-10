//----------------------------  vectors_rhs_hp_02.cc  ---------------------------
//    $Id: vectors_rhs_hp_02.cc 23710 2011-05-17 04:50:10Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2006, 2007, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vectors_rhs_hp_02.cc  ---------------------------

// we get this crash:
//
// An error occurred in line <4677> of file </scratch/deal-trunk/deal.II/include/deal.II/numerics/vectors.templates.h> in function
//     void dealii::VectorTools::compute_no_normal_flux_constraints(const DH<dim, spacedim>&, unsigned int, const std::set<unsigned char>&, dealii::ConstraintMatrix&, const dealii::Mapping<dim, spacedim>&) [with int dim = 2, DH = dealii::DoFHandler, int spacedim = 2]
// The violated condition was: 
//     std::fabs(determinant (t)) > 1e-3
// The name and call sequence of the exception was:
//     ExcMessage("Found a set of normal vectors that are nearly collinear.")
// Additional Information: 
// Found a set of normal vectors that are nearly collinear.

// the crash does not happen without a boundary object. Same error occures with HyperShellBoundary instead of HalfHyperShellBoundary.

#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1.h>
#include <numerics/vectors.h>

using namespace dealii;

template <int dim>
void run()
{
 Triangulation<dim> triangulation;
 FESystem<dim> fe(FE_Q<dim>(1), dim);
 DoFHandler<dim> dof_handler (triangulation);
 ConstraintMatrix constraints;
 MappingQ<dim> mapping(1);

 GridGenerator::half_hyper_shell (triangulation,
                                  Point<dim>(),
                                  0.5,
                                  1);

 static HalfHyperShellBoundary<dim> boundary;
 triangulation.set_boundary (0, boundary);

 dof_handler.distribute_dofs (fe);

 std::set<unsigned char> no_normal_flux_boundaries;
 no_normal_flux_boundaries.insert (0);
 VectorTools::compute_no_normal_flux_constraints
   (dof_handler, 0,
    no_normal_flux_boundaries,
    constraints, mapping);
}

int main (int argc, char *argv[])
{
  std::ofstream logfile ("no_flux_10/output");
  logfile.precision (4);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console (0);

  run<2> ();
}
