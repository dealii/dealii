//----------------------------  refinement_listener_02.cc  ---------------------------
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
//----------------------------  refinement_listener_02.cc  ---------------------------


// test the various refinement listener functions. do this new-style,
// i.e. through the signals mechanism


#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <cstdio>

std::ofstream logfile("refinement_listener_02/output");



template <int dim, int spacedim>
void
pre_refinement_notification (const std::string &prefix,
			     const Triangulation<dim, spacedim> &tria)
{
  deallog << prefix << ' ' << "Pre-refinement: " << tria.n_active_cells() << std::endl;
}


template <int dim, int spacedim>
void
post_refinement_notification (const std::string &prefix,
			      const Triangulation<dim, spacedim> &tria)
{
  deallog << prefix << ' ' << "Post-refinement: " << tria.n_active_cells() << std::endl;
}


template <int dim, int spacedim>
void
copy_notification (const std::string &prefix,
		   const Triangulation<dim, spacedim> &old_tria,
		   const Triangulation<dim, spacedim> &new_tria)
{
  deallog << prefix << ' ' << "Copy: "
	  << old_tria.n_active_cells() << ' '
	  << new_tria.n_active_cells() << std::endl;
}


template <int dim, int spacedim>
void
create_notification (const std::string &prefix,
		     const Triangulation<dim, spacedim> &tria)
{
  deallog << prefix << ' ' << "Create: " << tria.n_active_cells() << std::endl;
}


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;

  Triangulation<dim> tria_1, tria_2;

  GridGenerator::hyper_cube(tria_2);

  boost::signals2::connection connections_1[4]
    = {tria_1.signals.pre_refinement
       .connect (std_cxx1x::bind (&pre_refinement_notification<dim,dim>,
				  "tria_1",
				  std_cxx1x::cref(tria_1))),
       tria_1.signals.post_refinement
       .connect (std_cxx1x::bind (&post_refinement_notification<dim,dim>,
				  "tria_1",
				  std_cxx1x::cref(tria_1))),
       tria_1.signals.create
       .connect (std_cxx1x::bind (&create_notification<dim,dim>,
				  "tria_1",
				  std_cxx1x::cref(tria_1))),
       tria_1.signals.copy
       .connect (std_cxx1x::bind (&copy_notification<dim,dim>,
				  "tria_1",
				  std_cxx1x::_1,
				  std_cxx1x::cref(tria_1)))
  };
  boost::signals2::connection connections_2[4]
    = {tria_2.signals.pre_refinement
       .connect (std_cxx1x::bind (&pre_refinement_notification<dim,dim>,
				  "tria_2",
				  std_cxx1x::cref(tria_2))),
       tria_2.signals.post_refinement
       .connect (std_cxx1x::bind (&post_refinement_notification<dim,dim>,
				  "tria_2",
				  std_cxx1x::cref(tria_2))),
       tria_2.signals.create
       .connect (std_cxx1x::bind (&create_notification<dim,dim>,
				  "tria_2",
				  std_cxx1x::cref(tria_2))),
       tria_2.signals.copy
       .connect (std_cxx1x::bind (&copy_notification<dim,dim>,
				  "tria_2",
				  std_cxx1x::_1,
				  std_cxx1x::cref(tria_2)))
  };



				   // this should print the create note
  GridGenerator::hyper_cube(tria_1);

				   // this should print the pre- and
				   // post-refinement note
  tria_1.refine_global (1);

				   // this should print the copy note
  tria_1.clear ();
  tria_1.copy_triangulation (tria_2);

				   // no longer print anything
  for (unsigned int i=0; i<4; ++i)
    connections_1[i].disconnect ();

  tria_1.refine_global (2);

  (void)connections_2;
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
