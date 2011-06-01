//----------------------------  refinement_listener_01.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2008, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  refinement_listener_01.cc  ---------------------------


// test the various refinement listener functions


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("refinement_listener_01/output");


template <int dim, int spacedim = dim>
class RefinementListener :
  public Triangulation<dim,spacedim>::RefinementListener
{
  public:
    RefinementListener (const std::string &prefix)
		    :
		    prefix(prefix)
      {}

    virtual
    void
    pre_refinement_notification (const Triangulation<dim, spacedim> &tria)
      {
	deallog << prefix << ' ' << "Pre-refinement: " << tria.n_active_cells() << std::endl;
      }

    virtual
    void
    post_refinement_notification (const Triangulation<dim, spacedim> &tria)
      {
	deallog << prefix << ' ' << "Post-refinement: " << tria.n_active_cells() << std::endl;
      }

    virtual
    void
    copy_notification (const Triangulation<dim, spacedim> &old_tria,
		       const Triangulation<dim, spacedim> &new_tria)
      {
	deallog << prefix << ' ' << "Copy: "
		<< old_tria.n_active_cells() << ' '
		<< new_tria.n_active_cells() << std::endl;
      }

    virtual
    void
    create_notification (const Triangulation<dim, spacedim> &tria)
      {
	deallog << prefix << ' ' << "Create: " << tria.n_active_cells() << std::endl;
      }

  private:
    std::string prefix;
};


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;

  Triangulation<dim> tria_1, tria_2;

  GridGenerator::hyper_cube(tria_2);

  RefinementListener<dim> rl_1 ("tria_1");
  RefinementListener<dim> rl_2 ("tria_2");
  tria_1.add_refinement_listener (rl_1);
  tria_2.add_refinement_listener (rl_2);

				   // this should print the create note
  GridGenerator::hyper_cube(tria_1);

				   // this should print the pre- and
				   // post-refinement note
  tria_1.refine_global (1);

				   // this should print the copy note
  tria_1.clear ();
  tria_1.copy_triangulation (tria_2);

				   // no longer print anything
  tria_1.remove_refinement_listener (rl_1);
  tria_1.refine_global (2);
}


int main ()
{
  deallog << std::setprecision(2);
  logfile << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();
  test<3> ();

  return 0;
}
