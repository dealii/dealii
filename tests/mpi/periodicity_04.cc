/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

//
// Test parallel::distributed::Triangulation::add_periodicity
// for correct behavior on non standard oriented meshes.
// This test case is similar to dof_tools_21_b.cc.
//

#include "../tests.h"
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <int dim>
void set_periodicity(parallel::distributed::Triangulation<dim> &triangulation, bool reverse)
{
  typename Triangulation<dim>::cell_iterator cell_1 = triangulation.begin();
  typename Triangulation<dim>::cell_iterator cell_2 = cell_1++;
  typename Triangulation<dim>::face_iterator face_1;
  typename Triangulation<dim>::face_iterator face_2;

  // Look for the two outermost faces:
  for (unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
    {
      if (cell_1->face(j)->center()(dim-1) > 2.9)
        face_1 = cell_1->face(j);
      if (cell_2->face(j)->center()(dim-1) < -2.9)
        face_2 = cell_2->face(j);
    }
  face_1->set_boundary_id(42);
  face_2->set_boundary_id(43);

  std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
  periodicity_vector;

  if (reverse)
    GridTools::collect_periodic_faces (triangulation, 43, 42, dim-1, periodicity_vector);
  else
    GridTools::collect_periodic_faces (triangulation, 42, 43, dim-1, periodicity_vector);

  triangulation.add_periodicity(periodicity_vector);

  triangulation.refine_global(3);
}

/* The 2D case */
void generate_grid(parallel::distributed::Triangulation<2> &triangulation, int orientation)
{
  Point<2> vertices_1[]
  =
  {
    Point<2> (-1.,-3.),
    Point<2> (+1.,-3.),
    Point<2> (-1.,-1.),
    Point<2> (+1.,-1.),
    Point<2> (-1.,+1.),
    Point<2> (+1.,+1.),
    Point<2> (-1.,+3.),
    Point<2> (+1.,+3.),
  };
  std::vector<Point<2> > vertices (&vertices_1[0], &vertices_1[8]);

  std::vector<CellData<2> > cells (2, CellData<2>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<2>::vertices_per_cell] = {0, 1,  2,  3};

  /* cell 1 */
  int cell_vertices_1[2][GeometryInfo<2>::vertices_per_cell]
  =
  {
    {4,5,6,7},
    {7,6,5,4},
  };

  for (unsigned int j=0; j<GeometryInfo<2>::vertices_per_cell; ++j)
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;

  triangulation.create_triangulation(vertices, cells, SubCellData());
}




/* The 3D case */
void generate_grid(parallel::distributed::Triangulation<3> &triangulation, int orientation)
{
  Point<3> vertices_1[]
  =
  {
    Point<3> (-1.,-1.,-3.),
    Point<3> (+1.,-1.,-3.),
    Point<3> (-1.,+1.,-3.),
    Point<3> (+1.,+1.,-3.),
    Point<3> (-1.,-1.,-1.),
    Point<3> (+1.,-1.,-1.),
    Point<3> (-1.,+1.,-1.),
    Point<3> (+1.,+1.,-1.),
    Point<3> (-1.,-1.,+1.),
    Point<3> (+1.,-1.,+1.),
    Point<3> (-1.,+1.,+1.),
    Point<3> (+1.,+1.,+1.),
    Point<3> (-1.,-1.,+3.),
    Point<3> (+1.,-1.,+3.),
    Point<3> (-1.,+1.,+3.),
    Point<3> (+1.,+1.,+3.)
  };
  std::vector<Point<3> > vertices (&vertices_1[0], &vertices_1[16]);

  std::vector<CellData<3> > cells (2, CellData<3>());

  /* cell 0 */
  int cell_vertices_0[GeometryInfo<3>::vertices_per_cell] = {0, 1,  2,  3,  4,  5,  6,  7};

  /* cell 1 */
  int cell_vertices_1[8][GeometryInfo<3>::vertices_per_cell]
  =
  {
    {8,9,10,11,12,13,14,15},
    {9,11,8,10,13,15,12,14},
    {11,10,9,8,15,14,13,12},
    {10,8,11,9,14,12,15,13},
    {13,12,15,14,9,8,11,10},
    {12,14,13,15,8,10,9,11},
    {14,15,12,13,10,11,8,9},
    {15,13,14,12,11,9,10,8},
  };

  for (unsigned int j=0; j<GeometryInfo<3>::vertices_per_cell; ++j)
    {
      cells[0].vertices[j] = cell_vertices_0[j];
      cells[1].vertices[j] = cell_vertices_1[orientation][j];
    }
  cells[0].material_id = 0;
  cells[1].material_id = 0;


  triangulation.create_triangulation(vertices, cells, SubCellData());
}

template<int dim>
void check
(const unsigned int  orientation,
 bool                reverse)
{
  dealii::parallel::distributed::Triangulation<dim> triangulation (MPI_COMM_WORLD);

  generate_grid(triangulation, orientation);
  set_periodicity(triangulation, reverse);

  //first without refinement
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  ConstraintMatrix      constraints;

  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);

  constraints.reinit(locally_relevant_dofs);
  {
    std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> >
    periodicity_vector;

    GridTools::collect_periodic_faces (dof_handler, 42, 43, dim-1, periodicity_vector);

    DoFTools::make_periodicity_constraints<DoFHandler<dim> >
    (periodicity_vector, constraints);
  }
  constraints.close ();

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  constraints.print(deallog.get_file_stream());

  unsigned int n_local_constraints =0 ;

  std::map<types::global_dof_index, Point<dim> > support_points;
  DoFTools::map_dofs_to_support_points (MappingQGeneric<dim>(1), dof_handler, support_points);
  IndexSet constraints_lines = constraints.get_local_lines();

  for (unsigned int i=0; i<constraints_lines.n_elements(); ++i)
    {
      const unsigned int line = constraints_lines.nth_index_in_set(i);
      if (constraints.is_constrained(line))
        {
          const std::vector<std::pair<types::global_dof_index, double > > *entries
            = constraints.get_constraint_entries(line);
          Assert(entries->size()==1, ExcInternalError());
          const Point<dim> point1 = support_points[line];
          const Point<dim> point2 = support_points[(*entries)[0].first];
          Tensor<1, dim> difference = point1-point2;
          difference[dim-1] = 0.;
          AssertThrow(difference.norm()<1.e-9, ExcInternalError());
          if (locally_owned_dofs.is_element(line))
            ++n_local_constraints;
        }
    }
  const unsigned int n_constraints
    = Utilities::MPI::sum(n_local_constraints, MPI_COMM_WORLD);
  const unsigned int n_expected_constraints = Utilities::fixed_int_power<9,dim-1>::value;
  if (myid==0)
    deallog << "n_constraints: " << n_constraints
            << " n_expected_constraints: " << n_expected_constraints
            << std::endl;
  AssertThrow(n_constraints == n_expected_constraints,
              ExcDimensionMismatch(n_constraints, n_expected_constraints));

  //now refine and check if the neighboring faces are correctly found
  typename Triangulation<dim>::active_cell_iterator cell;
  for (cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned() && cell->center()(dim-1)>0)
      cell->set_refine_flag();

  triangulation.execute_coarsening_and_refinement();

  typedef std::pair<typename Triangulation<dim>::cell_iterator, unsigned int> CellFace;
  const typename std::map<CellFace, std::pair<CellFace, std::bitset<3> > > &face_map = triangulation.get_periodic_face_map();
  typename std::map<CellFace, std::pair<CellFace, std::bitset<3> > >::const_iterator it;
  int sum_of_pairs_local = face_map.size();
  int sum_of_pairs_global;
  MPI_Allreduce(&sum_of_pairs_local, &sum_of_pairs_global, 1, MPI_INT, MPI_SUM, triangulation.get_communicator());
  Assert(sum_of_pairs_global>0, ExcInternalError());
  for (it = face_map.begin(); it != face_map.end(); ++it)
    {
      const typename Triangulation<dim>::cell_iterator cell_1 = it->first.first;
      const unsigned int face_no_1 = it->first.second;
      const typename Triangulation<dim>::cell_iterator cell_2 = it->second.first.first;
      const unsigned int face_no_2 = it->second.first.second;
      const Point<dim> face_center_1 = cell_1->face(face_no_1)->center();
      const Point<dim> face_center_2 = cell_2->face(face_no_2)->center();
      Assert(std::min(std::abs(face_center_1(dim-1)-3.), std::abs(face_center_1(dim-1)+3.))<1.e-8, ExcInternalError());
      Assert(std::min(std::abs(face_center_2(dim-1)-3.), std::abs(face_center_2(dim-1)+3.))<1.e-8, ExcInternalError());
      if (cell_1->level()== cell_2->level())
        for (unsigned int c = 0; c<dim-1; ++c)
          if (std::abs(face_center_1(c)-face_center_2(c))>1.e-8)
            {
              std::cout << "face_center_1: " << face_center_1 << std::endl;
              std::cout << "face_center_2: " << face_center_2 << std::endl;
              typename std::map<CellFace, std::pair<CellFace, std::bitset<3> > >::const_iterator it;
              for (it = triangulation.get_periodic_face_map().begin();
                   it !=  triangulation.get_periodic_face_map().end(); ++it)
                {
                  std::cout << "The cell with center " << it->first.first->center()
                            << " has on face " << it->first.second
                            << " as neighbor the cell with center " << it->second.first.first->center()
                            << " on face " << it->second.first.second
                            << std::endl;
                }
              Assert(false, ExcInternalError());
            }
    }
}

int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

      MPILogInitAll log;

      deallog.threshold_double(1.e-10);
      {
        deallog << "Test for 2D" << std::endl << std::endl;
        for (int i = 0; i < 2; ++i)
          {
            deallog << "Triangulation: " << i << std::endl;
            check<2> (i, false);
            check<2> (i, true);
          }

        deallog << "Test for 3D" << std::endl << std::endl;
        for (int i = 0; i < 8; ++i)
          {
            // Generate a triangulation and match:
            deallog << "Triangulation: " << i << std::endl;
            check<3> (i, false);
            check<3> (i, true);
          }
      }
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
