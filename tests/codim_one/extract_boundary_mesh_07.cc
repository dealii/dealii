//----------------------------  extract_boundary_mesh_07.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  extract_boundary_mesh_07.cc  ---------------------------

// verify that we can indeed call extract_boundary_mesh with DoFHandler
// arguments. based on a test by Korosh Taebi

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <fstream>
#include <iostream>

namespace Step38
{
  using namespace dealii;

  template <int spacedim>
  class Extract_Mesh_Test
  {
  public:
    Extract_Mesh_Test ();
    void run ();

  private:
    static const unsigned int boundary_dim = spacedim-1;

    Triangulation<spacedim> 				 volume_mesh_triangulation;
    Triangulation<boundary_dim,spacedim>   boundary_triangulation;

    FE_Q<spacedim>					 	 space_fe;
    FE_Q<boundary_dim,spacedim>            boundary_fe;

    DoFHandler<spacedim>					 space_dof_handler;
    DoFHandler<boundary_dim,spacedim>      contact_dof_handler;

  };


  template <int spacedim>
  Extract_Mesh_Test<spacedim>::
  Extract_Mesh_Test ()
  :
  space_fe (spacedim),
  boundary_fe (1),
  space_dof_handler(volume_mesh_triangulation),
  contact_dof_handler(boundary_triangulation)
  {}

  template <int spacedim>
  void Extract_Mesh_Test<spacedim>::run ()
  {
    GridGenerator::hyper_cube(volume_mesh_triangulation);
    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);

    space_dof_handler.distribute_dofs(space_fe);

    std::map<typename DoFHandler<boundary_dim, spacedim>::cell_iterator,
	     typename DoFHandler<spacedim>::face_iterator>
      element_assignment =
      GridTools::extract_boundary_mesh(space_dof_handler,
				       contact_dof_handler,
				       boundary_ids);

    contact_dof_handler.distribute_dofs(boundary_fe);

    typename std::map<typename DoFHandler<boundary_dim, spacedim>::cell_iterator,typename DoFHandler<spacedim>::face_iterator>::iterator Iterator;

    for (Iterator = element_assignment.begin();Iterator != element_assignment.end(); ++Iterator)
      {
	deallog << "element_assignment maps " << Iterator->first << " onto " << Iterator->second << std::endl;
      }

  }
}


int main ()
{
  std::ofstream logfile("extract_boundary_mesh_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {
    using namespace dealii;
    using namespace Step38;

    Extract_Mesh_Test<2> Test;
    Test.run();
  }
 }
