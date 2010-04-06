//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Check downstream numbering


#include "../tests.h"
#include <base/logstream.h>
#include <base/function_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <fstream>



template <int dim>
void
print_dofs (const DoFHandler<dim> &dof)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  std::vector<unsigned int> v (fe.dofs_per_cell);
  std_cxx1x::shared_ptr<FEValues<dim> > fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std_cxx1x::shared_ptr<FEValues<dim> >(new FEValues<dim>(fe, quad, update_q_points));
    }

  for (typename DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != 0)
	fevalues->reinit(cell);

      cell->get_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
	if (fevalues.get() != 0)
	  deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
	else
	  deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs (const MGDoFHandler<dim> &dof, unsigned int level)
{
  const FiniteElement<dim>& fe = dof.get_fe();
  std::vector<unsigned int> v (fe.dofs_per_cell);
  std_cxx1x::shared_ptr<FEValues<dim> > fevalues;

  if (fe.has_support_points())
    {
      Quadrature<dim> quad(fe.get_unit_support_points());
      fevalues = std_cxx1x::shared_ptr<FEValues<dim> >(new FEValues<dim>(fe, quad, update_q_points));
    }

  for (typename MGDoFHandler<dim>::cell_iterator cell=dof.begin(level);
       cell != dof.end(level); ++cell)
    {
      Point<dim> p = cell->center();
      if (fevalues.get() != 0)
	fevalues->reinit(cell);

      cell->get_mg_dof_indices (v);
      for (unsigned int i=0; i<v.size(); ++i)
	if (fevalues.get() != 0)
	  deallog << fevalues->quadrature_point(i) << '\t' << v[i] << std::endl;
	else
	  deallog << p << '\t' << v[i] << std::endl;
      deallog << std::endl;
    }
}


template <int dim>
void
check_renumbering(MGDoFHandler<dim>& mgdof)
{
  const FiniteElement<dim>& element = mgdof.get_fe();
  DoFHandler<dim>& dof = mgdof;
  deallog << element.get_name() << std::endl;

				   // Prepare a reordering of
				   // components for later use
  std::vector<unsigned int> order(element.n_components());
  for (unsigned int i=0; i<order.size(); ++i) order[i] = order.size()-i-1;

  Point<dim> direction;
  for (unsigned int i=0;i<dim;++i)
    direction(i) = -5.+i;

  deallog << std::endl << "Downstream numbering cell-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction);
  print_dofs (dof);
				   // Check level ordering
  for (unsigned int level=0;level<dof.get_tria().n_levels();++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction);
      print_dofs (mgdof, level);
    }

  deallog << std::endl << "Downstream numbering dof-wise" << std::endl;
  DoFRenumbering::downstream(dof, direction, true);
  print_dofs (dof);
				   // Check level ordering
  for (unsigned int level=0;level<dof.get_tria().n_levels();++level)
    {
      deallog << "Level " << level << std::endl;
      DoFRenumbering::downstream(mgdof, level, direction, true);
      print_dofs (mgdof, level);
    }

}


template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1., 1.);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  tr.refine_global(3-dim);

  MGDoFHandler<dim> mgdof(tr);

  FE_Q<dim> q2(2);
  FE_DGQ<dim> dgq1(1);
  FESystem<dim> system (q2, 2, dgq1, 1);

  mgdof.distribute_dofs(q2);
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(dgq1);
  check_renumbering(mgdof);
  mgdof.clear();

  mgdof.distribute_dofs(system);
  check_renumbering(mgdof);
  mgdof.clear();
}


int main ()
{
  std::ofstream logfile ("dof_renumbering_02/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("1d");
  check<1> ();
  deallog.pop ();
  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
