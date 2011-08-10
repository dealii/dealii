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


// we used to get this crash:
//
// An error occurred in line <4646> of file </w/heister/deal-trunk/deal.II/include/deal.II/numerics/vectors.templates.h> in function
//     static void dealii::VectorTools::compute_no_normal_flux_constraints(const DH<dim, spacedim>&, unsigned int, const std::set<unsigned char>&, dealii::ConstraintMatrix&, const dealii::Mapping<dim, spacedim>&) [with int dim = 3, DH = dealii::DoFHandler, int spacedim = 3]
// The violated condition was:
//     std::fabs (tangent.norm()-1) < 1e-12
// The name and call sequence of the exception was:
//     ExcInternalError()
//
// quarter_hyper_shell works, even though it is a very similar mesh.
//
// this was fixed with r24044 together with the no_flux_07 test that
// reduces it to its essence

#include "../tests.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vectors.h>


template <int dim>
void
check ()
{
  Triangulation<dim> tr;
   GridGenerator::half_hyper_shell (tr,
			      Point<dim>(),
			      0.5,
			      1,0);

  ConstraintMatrix cm;
  MappingQ<dim> mapping(1);

  FESystem<dim> fe(FE_Q<dim>(1),dim);
  DoFHandler<dim> dofh(tr);

  dofh.distribute_dofs (fe);

  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  VectorTools::compute_no_normal_flux_constraints (dofh, 0, no_normal_flux_boundaries, cm, mapping);

  cm.print (deallog.get_file_stream ());
}



int main ()
{
  std::ofstream logfile ("no_flux_08/output");
  logfile.precision (4);
  logfile.setf(std::ios::fixed);
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
