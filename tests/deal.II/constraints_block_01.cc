// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2013 by the deal.II authors
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



// test that we fixed a case where we tried to deal with some constraints in a
// block matrix where the last few rows of one of the blocks were empty and we
// ran into an unrelated assertion because we were accessing something beyond
// the end of the array
//
// testcase by Jason Sheldon

#include "../tests.h"
#include <base/logstream.h>
#include <base/utilities.h>
#include <base/quadrature_lib.h>

#include <numerics/vector_tools.h>
#include <numerics/matrix_tools.h>
#include <numerics/data_out.h>

#include <lac/vector.h>
#include <lac/sparse_direct.h>
#include <lac/full_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/constraint_matrix.h>

#include <dofs/dof_accessor.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>

#include <hp/dof_handler.h>
#include <hp/q_collection.h>
#include <hp/fe_collection.h>
#include <hp/fe_values.h>

#include <fe/fe.h>
#include <fe/fe_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_nothing.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <grid/grid_generator.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria.h>

#include <fstream>
#include <cmath>
#include <cstdlib>


std::ofstream logfile("output");


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // setting up some constants
  const unsigned int dim = 2;
  const unsigned int solid_dim = 2*dim;
  const unsigned int fluid_dim = dim+1;
  const unsigned int mesh_dim  = dim;
  const unsigned int total_dim = solid_dim + fluid_dim + mesh_dim;

  // make the tria and domain
  Triangulation<dim>   tria;

  GridGenerator::hyper_cube (tria);
  tria.refine_global(1);

  Point<dim> real_cell_center;

  // create a solid domain (0) in the
  // lower left hand corner of a 2x2 grid
  // the rest is fluid/mesh (1)

  for (Triangulation<dim>::active_cell_iterator
       cell = tria.begin_active();
       cell != tria.end(); ++cell)
    {
      real_cell_center = cell->center();
      if (real_cell_center(0) < 0.5
          && real_cell_center(1) < 0.5 )
        cell->set_material_id (0); // solid
      else
        cell->set_material_id (1); //fluid
    }//cell

  //create the FE spaces for the solid and the fluid/mesh
  //each are padded with FE_Nothing to be equal length

  std::string solid_fe_name = "FESystem[FE_Q(2)^2-FE_Q(2)^2-FE_Nothing()^2-FE_Nothing()-FE_Nothing()^2]";
  std::string fluid_fe_name = "FESystem[FE_Nothing()^2-FE_Nothing()^2-FE_Q(2)^2-FE_Q(1)-FE_Q(2)^2]";

  hp::FECollection<dim> fe_collection;
  FiniteElement<dim> *solid_fe = FETools::get_fe_from_name<dim>(solid_fe_name);
  FiniteElement<dim> *fluid_fe = FETools::get_fe_from_name<dim>(fluid_fe_name);

  deallog << "Solid FE Space: " << solid_fe->get_name() << std::endl;
  deallog << "Fluid/Mesh FE Space: " << fluid_fe->get_name() << std::endl;

  fe_collection.push_back(*solid_fe);
  fe_collection.push_back(*fluid_fe);

  hp::DoFHandler<dim> dh(tria);

  for (hp::DoFHandler<dim>::active_cell_iterator
       cell = dh.begin_active();
       cell != dh.end(); ++cell)
    {
      if (int(cell->material_id()) == 0)
        cell->set_active_fe_index (0);
      else
        cell->set_active_fe_index (1);
    }//cell
  dh.distribute_dofs(fe_collection);

  std::vector<unsigned int> block_component (total_dim,0);

  for (unsigned int comp=0; comp < total_dim; comp++)
    {
      if ( comp < solid_dim )
        block_component[comp] = 0;
      else if ( comp < solid_dim+fluid_dim )
        block_component[comp] = 1;
      else
        block_component[comp] = 2;
    }//comp

  std::vector<types::global_dof_index> dofs_per_block(3,0);//3 blocks, count dofs:
  DoFTools::count_dofs_per_component (dh, dofs_per_block, false, block_component);

  DoFRenumbering::component_wise(dh, block_component);

  //build the sparsitypattern

  BlockSparsityPattern   block_sparsity_pattern;
  {
    BlockCompressedSimpleSparsityPattern  csp(3,3);
    csp.block(0,0).reinit (dofs_per_block[0],dofs_per_block[0]);//solid-solid
    csp.block(0,1).reinit (dofs_per_block[0],dofs_per_block[1]);//solid-fluid
    csp.block(0,2).reinit (dofs_per_block[0],dofs_per_block[2]);//solid-mesh

    csp.block(1,0).reinit (dofs_per_block[1],dofs_per_block[0]);//fluid-solid
    csp.block(1,1).reinit (dofs_per_block[1],dofs_per_block[1]);//fluid-fluid
    csp.block(1,2).reinit (dofs_per_block[1],dofs_per_block[2]);//fluid-mesh

    csp.block(2,0).reinit (dofs_per_block[2],dofs_per_block[0]);//mesh-solid
    csp.block(2,1).reinit (dofs_per_block[2],dofs_per_block[1]);//mesh-fluid
    csp.block(2,2).reinit (dofs_per_block[2],dofs_per_block[2]);//mesh-mesh

    csp.collect_sizes();

    //enforce coupling across cells and interface

    Table<2,DoFTools::Coupling> cell_coupling (fe_collection.n_components(), fe_collection.n_components());
    Table<2,DoFTools::Coupling> face_coupling (fe_collection.n_components(), fe_collection.n_components());

    for (unsigned int c=0; c<fe_collection.n_components(); ++c)
      {
        for (unsigned int d=0; d<fe_collection.n_components(); ++d)
          {
            if (((c<solid_dim) && (d<solid_dim)) //couples solid dims with solid dims
                ||
                (((c>=solid_dim) && (d>=solid_dim)) //couples fluid dims and mesh dims
                 && !((c==solid_dim+dim) && (d==solid_dim+dim)))) //fluid pressure does not couple with itself
              {
                cell_coupling[c][d] = DoFTools::always;
                cell_coupling[d][c] = DoFTools::always;
              }//cell_coupling

            if ((c<solid_dim) && (d>=solid_dim)) //couples solid dims with fluid and mesh dims on the interface
              {
                face_coupling[c][d] = DoFTools::always;
                face_coupling[d][c] = DoFTools::always;
              }//face_coupling
          }//d
      }//c
    DoFTools::make_flux_sparsity_pattern (dh, csp, cell_coupling, face_coupling);
    block_sparsity_pattern.copy_from (csp);
  }

  // build matrices and vectors

  BlockSparseMatrix<double> system_matrix(block_sparsity_pattern);
  BlockVector<double> system_rhs(dofs_per_block);

  // set up constraints for solution and solution update
  /**
   * These constraints enforce that the solid displacement
   * and mesh displacement are the same on the interface.
   *
   * They also constrain that the solid velocity
   * and fluid velocity are the same on the interface
   *
   * The interface check is simplified for this 2x2 case
   */

  ConstraintMatrix constraints;

  const unsigned int   dofs_per_fl_msh_face = fluid_fe->dofs_per_face;
  const unsigned int   dofs_per_solid_face  = solid_fe->dofs_per_face;
  std::vector<types::global_dof_index> fl_msh_face_dof_indices (dofs_per_fl_msh_face);
  std::vector<types::global_dof_index> solid_face_dof_indices  (dofs_per_solid_face );

  std::vector<std::pair<unsigned int, unsigned int> > solid_fluid_pairs;
  std::vector<std::pair<unsigned int, unsigned int> > solid_mesh_pairs;

  for (hp::DoFHandler<dim>::active_cell_iterator
       cell = dh.begin_active(); cell!=dh.end(); ++cell) //loops over the cells
    {
      if (int(cell->material_id()) == 1)
        {
          //Only loop over cells in the fluid region
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            {
              if (!cell->at_boundary(f))
                {
                  bool face_is_on_interface = false;
                  //checks to see if cell face neighbor is in the solid domain
                  if (int(cell->neighbor(f)->material_id()) == 0 )
                    face_is_on_interface = true;

                  if (face_is_on_interface)
                    {
                      std::vector<unsigned int> solid_disp_dof, solid_vel_dof, fluid_vel_dof, mesh_disp_dof;

                      cell->face(f)->get_dof_indices (solid_face_dof_indices, 0);
                      cell->face(f)->get_dof_indices (fl_msh_face_dof_indices, 1);

                      for (unsigned int i=0; i<dofs_per_solid_face; ++i)
                        {
                          unsigned int comp = solid_fe->face_system_to_component_index(i).first;
                          if (comp <dim)
                            solid_disp_dof.push_back(i);
                          if (comp >= dim && comp < solid_dim)
                            solid_vel_dof.push_back(i);
                        }//i
                      for (unsigned int i=0; i<dofs_per_fl_msh_face; ++i)
                        {
                          unsigned int comp = fluid_fe->face_system_to_component_index(i).first;
                          if (comp >= solid_dim && comp < solid_dim+dim)
                            fluid_vel_dof.push_back(i);
                          if (comp>=solid_dim+fluid_dim)
                            mesh_disp_dof.push_back(i);
                        }//i

                      for (unsigned int i=0; i<solid_vel_dof.size(); ++i)
                        {
                          //in this example solid_vel_dof.size()==fluid_vel_dof.size()
                          constraints.add_line  (fl_msh_face_dof_indices[fluid_vel_dof[i]]);
                          constraints.add_entry (fl_msh_face_dof_indices[fluid_vel_dof[i]],
                                                 solid_face_dof_indices[solid_vel_dof[i]], 1.0);

                          solid_fluid_pairs.push_back(std::pair<unsigned int,unsigned int>
                                                      (solid_face_dof_indices[solid_vel_dof[i]],
                                                       fl_msh_face_dof_indices[fluid_vel_dof[i]]));
                        }//i
                      for (unsigned int i=0; i<solid_disp_dof.size(); ++i)
                        {
                          constraints.add_line  (fl_msh_face_dof_indices[mesh_disp_dof[i]]);
                          constraints.add_entry (fl_msh_face_dof_indices[mesh_disp_dof[i]],
                                                 solid_face_dof_indices[solid_disp_dof[i]], 1.0);

                          solid_mesh_pairs.push_back(std::pair<unsigned int,unsigned int>
                                                     (solid_face_dof_indices[solid_disp_dof[i]],
                                                      fl_msh_face_dof_indices[mesh_disp_dof[i]]));
                        }
                    }//at interface check
                }//not at boundary check
            }//face
        }//is this in the fluid material region?
    }//cell

  constraints.close();

  // prints out which dofs are coupled
  deallog<<"---------------Coupled dofs---------------"<<std::endl;
  for (unsigned int i=0; i<solid_fluid_pairs.size(); i++)
    {
      deallog<<"solid dof: "<< solid_fluid_pairs[i].first
             <<", fluid dof: "<<solid_fluid_pairs[i].second<<std::endl;
    }
  for (unsigned int i=0; i<solid_fluid_pairs.size(); i++)
    {
      deallog<<"solid dof: "<<solid_mesh_pairs[i].first
             <<", mesh dof: "<<solid_mesh_pairs[i].second<<std::endl;
    }
  deallog<<"------------------------------------------"<<std::endl;


  //code crashes in the fluid assembly, mocked up below
  {
    //fluid assembly
    const unsigned int dofs_per_cell   = fluid_fe->dofs_per_cell;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    for (hp::DoFHandler<dim>::active_cell_iterator
         cell = dh.begin_active();
         cell!=dh.end();
         ++cell) //loops over the cells
      {
        if (int(cell->material_id()) == 1)//are we in the fluid region?
          {
            cell->get_dof_indices (local_dof_indices);

            // normally the governing equations would be added
            // to the cell matrix and right hand side here

            deallog<<"I am about to distribute cell: "<<cell<<std::endl;
            deallog<<"Cell: "<<cell<<" has it's center at "<<cell->center()<<std::endl;

            constraints.distribute_local_to_global(
              cell_matrix,
              cell_rhs,
              local_dof_indices,
              system_matrix,
              system_rhs);

          }//is this in the fluid material region?

      }//cell

  }//assemble fluid system

  return 0;
}

