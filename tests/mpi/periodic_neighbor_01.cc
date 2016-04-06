// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
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
//
/*
 * ====================
 * [ Test objectives: ]
 * ====================
 *
 * 1. To control if periodic faces have at most one level of refinement
 *    difference.
 *
 * 2. To test the following functions:
 *    - cell->has_periodic_neighbor(i_face)
 *    - cell->periodic_neighbor(i_face)
 *    - cell->periodic_neighbor_is_coarser(i_face)
 *    - cell->periodic_face_no(i_face)
 *    - cell->periodic_neighbor_of_periodic_neighbor(i_face)
 *    - cell->periodic_neighbor_of_coarser_periodic_neighbor(i_face)
 *    - cell->periodic_neighbor_child_on_subface(i_face, i_subface)
 *
 *
 * ===============
 * [ Test Setup: ]
 * ===============
 *
 * We start with a 4x4 mesh (in 2D) and 4x4x4 mesh (in 3D) and refine on of the
 * top middle elements twice, to get the following mesh:
 *
 *
 *                      Top periodic boundary
 *       -------------------------------------------------
 *       |           |     |  |  |     |     |           |
 *       |           |     |-----|     |     |           |
 *       |           |     |  |  |     |     |           |
 *       |           -------------------------           |
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       -------------------------------------------------
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       -------------------------------------------------
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       |           |           |           |           |
 *       -------------------------------------------------
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       |           -------------------------           |
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       |           |     |     |     |     |           |
 *       -------------------------------------------------
 *                    Bottom periodic boundary
 *
 *
 * =====================
 * [ Passing criteria: ]
 * =====================
 *
 * The test will be passed if the main three assertions are passed.
 */

#include "../tests.h"
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <string>


using namespace dealii;

template <int dim>
struct periodicity_tests
{
  typedef TriaIterator<CellAccessor<dim> > cell_iterator;
  typedef typename Triangulation<dim>::active_cell_iterator active_cell_iterator;
  periodicity_tests();

  unsigned refn_cycle;
  MPI_Comm mpi_comm;
  int comm_rank, comm_size;
  parallel::distributed::Triangulation<dim> the_grid;

  void refine_grid(const unsigned);
  void write_grid();
  void check_periodicity();
};

template <int dim>
periodicity_tests<dim>::periodicity_tests()
  : refn_cycle(0), mpi_comm(MPI_COMM_WORLD), the_grid(mpi_comm)
{
  Assert(dim == 2 || dim == 3, ExcMessage("Only implemented for the 2D and 3D case!"));
  comm_rank = Utilities::MPI::this_mpi_process(mpi_comm);
  comm_size = Utilities::MPI::n_mpi_processes(mpi_comm);
  std::vector<unsigned> repeats(dim, 2);
  Point<dim> p1, p2, periodic_transfer;
  if (dim == 2)
    {
      p2 = Point<dim>( 16., 16. );
      periodic_transfer = Point<dim>( 0.0, 16. );
    }
  if (dim == 3)
    {
      p2 = Point<dim>( 16., 16., 16. );
      periodic_transfer = Point<dim>( 0.0, 16., 0.0 );
    }
  GridGenerator::subdivided_hyper_rectangle(the_grid, repeats, p1, p2, true);
  std::vector<GridTools::PeriodicFacePair<cell_iterator> > periodic_faces;
  GridTools::collect_periodic_faces(the_grid, 2, 3, 0, periodic_faces, periodic_transfer);
  the_grid.add_periodicity(periodic_faces);
}

template <int dim>
void periodicity_tests<dim>::refine_grid(const unsigned n)
{
  if (n != 0 && refn_cycle == 0)
    {
      the_grid.refine_global(n);
      refn_cycle += n;
    }
  else if (n != 0)
    {
      Point<dim> refn_point;
      if (dim == 2)
        refn_point = Point<dim>( 6.5, 15.5 );
      if (dim == 3)
        refn_point = Point<dim>( 6.5, 15.5, 6.5 );
      for (unsigned i_refn = 0; i_refn < n; ++i_refn)
        {
          active_cell_iterator cell_it = the_grid.begin_active();
          for (; cell_it != the_grid.end(); ++cell_it)
            {
              if (cell_it->is_locally_owned() && cell_it->point_inside(refn_point))
                {
                  cell_it->set_refine_flag();
                  break;
                }
            }
          the_grid.execute_coarsening_and_refinement();
          ++refn_cycle;
        }
    }
}

template <int dim>
void periodicity_tests<dim>::write_grid()
{
  GridOut Grid1_Out;
  GridOutFlags::Svg svg_flags(
    1,                                       // line_thickness = 2,
    2,                                       // boundary_line_thickness = 4,
    false,                                   // margin = true,
    dealii::GridOutFlags::Svg::transparent,  // background = white,
    0,                                       // azimuth_angle = 0,
    0,                                       // polar_angle = 0,
    dealii::GridOutFlags::Svg::subdomain_id, // coloring = level_number,
    false, // convert_level_number_to_height = false,
    false, // label_level_number = true,
    true,  // label_cell_index = true,
    false, // label_material_id = false,
    false, // label_subdomain_id = false,
    true,  // draw_colorbar = true,
    true); // draw_legend = true
  Grid1_Out.set_flags(svg_flags);
  if (dim == 2)
    {
      std::ofstream Grid1_OutFile("Grid1" + std::to_string(refn_cycle) +
                                  std::to_string(comm_rank) + ".svg");
      Grid1_Out.write_svg(the_grid, Grid1_OutFile);
    }
  else
    {
      std::ofstream Grid1_OutFile("Grid1" + std::to_string(refn_cycle) +
                                  std::to_string(comm_rank) + ".msh");
      Grid1_Out.write_msh(the_grid, Grid1_OutFile);
    }
}

template <int dim>
void periodicity_tests<dim>::check_periodicity()
{
  typedef std::pair<cell_iterator, unsigned> cell_face_pair;
  typedef typename std::map<cell_face_pair, cell_face_pair>::iterator cell_face_map_it;

  for (int rank_i = 0; rank_i < comm_size; ++rank_i)
    {
      if (comm_rank == rank_i)
        {
          deallog << "All of the cells with periodic neighbors on rank "
                  << comm_rank << std::endl;
          active_cell_iterator cell_it = the_grid.begin_active();
          for (; cell_it != the_grid.end(); ++cell_it)
            {
              for (unsigned i_face = 0; i_face < GeometryInfo<dim>::faces_per_cell; ++i_face)
                {
                  if (cell_it->has_periodic_neighbor(i_face))
                    {
                      const cell_iterator &nb_it = cell_it->periodic_neighbor(i_face);
                      unsigned nb_i_face = cell_it->periodic_neighbor_face_no(i_face);
                      const cell_iterator &nb_of_nb_it = nb_it->periodic_neighbor(nb_i_face);
                      unsigned nb_of_nb_i_face = nb_it->periodic_neighbor_face_no(nb_i_face);
                      deallog << nb_it->face(nb_i_face)->center()
                              << " is a periodic neighbor of "
                              << cell_it->face(i_face)->center();
                      if (cell_it->periodic_neighbor_is_coarser(i_face))
                        {
                          deallog << ". And periodic neighbor is coarser." << std::endl;
                          /*
                           * Now, we want to do the standard test of:
                           *
                           * cell->neighbor(i_face)
                           *     ->periodic_neighbor_child_on_subface(face_no, subface_no)
                           * == cell
                           *
                           * The other important test for coarser neighbors is:
                           * cell->periodic_neighbor(i_face)
                           *     ->face(face_no)->child(subface_no)
                           * == cell->face(i_face)
                           * But, this is not the case for periodic neighbors. Because, we
                           * are not talking about the exact same faces.
                           */
                          std::pair<unsigned, unsigned> face_subface =
                            cell_it->periodic_neighbor_of_coarser_periodic_neighbor(i_face);
                          Assert(nb_it->periodic_neighbor_child_on_subface(
                                   face_subface.first, face_subface.second) == cell_it,
                                 ExcInternalError());
                        }
                      else if (nb_it->face(nb_i_face)->has_children())
                        {
                          deallog << ". And periodic neighbor is refined." << std::endl;
                        }
                      else
                        deallog << ". And periodic neighbor has the same level." << std::endl;
                      deallog << "If the periodic neighbor is not coarser, "
                              << nb_of_nb_it->face(nb_of_nb_i_face)->center()
                              << " should be equal to "
                              << cell_it->face(i_face)->center() << std::endl;
                      deallog
                          << "-------------------------------------------------------"
                          << std::endl;
                      Assert((cell_it->face(i_face)->center() -
                              nb_of_nb_it->face(nb_of_nb_i_face)->center())
                             .norm() < 1e-10 ||
                             cell_it->periodic_neighbor_is_coarser(i_face), ExcInternalError());

                      Assert(cell_it == nb_of_nb_it ||
                             cell_it->periodic_neighbor_is_coarser(i_face), ExcInternalError());
                    }
                }
            }
        }
      MPI_Barrier(mpi_comm);
    }
}

int main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv ,1);
      MPILogInitAll log;

      periodicity_tests<2> test_2D;
      /*
       * Let us first check the periodicity of the simple case of uniformly refined
       * mesh.
       */
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "First, we check the periodic faces on a uniformly refined mesh."
                << std::endl;
      test_2D.refine_grid(1);
      //  test_2D.write_grid();
      test_2D.check_periodicity();

      /* Next, we want to check the case of nonuniform mesh */
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "Next, we check the periodic faces on an adaptively refined mesh."
                << std::endl;
      test_2D.refine_grid(2);
      //  test_2D.write_grid();
      test_2D.check_periodicity();

      periodicity_tests<3> test_3D;
      /*
       * Let us first check the periodicity of the simple case of uniformly refined
       * mesh.
       */
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "First, we check the periodic faces on a uniformly refined mesh."
                << std::endl;
      test_3D.refine_grid(1);
      //  test_3D.write_grid();
      test_3D.check_periodicity();

      /* Next, we want to check the case of nonuniform mesh */
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "Next, we check the periodic faces on an adaptively refined mesh."
                << std::endl;
      test_3D.refine_grid(2);
      //  test_3D.write_grid();
      test_3D.check_periodicity();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}

