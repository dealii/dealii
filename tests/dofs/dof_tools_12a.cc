// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------
#include <deal.II/base/logstream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

#include "../tests.h"

// check
//  DoFTools::extract_dofs as in dof_tools_12 in parallel for fewer elements
// The output of this file is for every process the (boolean) vector of
// extracted locally owned DoFs (corresponding to local DoF indices), not the
// extracted vector of the global DoF, which has to be mapped with the locally
// owned IndexSet first.

void
output_bool_vector(std::vector<bool> &v)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << (v[i] ? '1' : '0');
  deallog << std::endl;
}


template <int dim, typename DoFHandlerType>
void
check_this(const DoFHandlerType &dof_handler)
{
  std::vector<bool> selected_dofs(dof_handler.n_locally_owned_dofs());
  std::vector<bool> mask(dof_handler.get_fe().n_components(), false);

  // only select first component
  mask[0] = true;
  DoFTools::extract_dofs<dim, dim, DoFHandlerType>(dof_handler,
                                                   ComponentMask(mask),
                                                   selected_dofs);
  output_bool_vector(selected_dofs);

  // also select last component
  mask.back() = true;
  DoFTools::extract_dofs<dim, dim, DoFHandlerType>(dof_handler,
                                                   ComponentMask(mask),
                                                   selected_dofs);
  output_bool_vector(selected_dofs);
}


template <int dim>
void
check(const FiniteElement<dim> &fe, const std::string &name)
{
  deallog << "Checking " << name << " in " << dim << "d:" << std::endl;

  // create tria and dofhandler
  // objects. set different boundary
  // and sub-domain ids
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  for (unsigned int ref = 0; ref < 2; ++ref)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned() && cell->center()(0) < .5 &&
            cell->center()(1) < .5)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  // setup DoFHandler
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // setup hp DoFHandler
  hp::FECollection<dim> fe_collection(fe);
  hp::DoFHandler<dim>   hp_dof_handler(tria);
  hp_dof_handler.distribute_dofs(fe_collection);

  check_this<dim, DoFHandler<dim>>(dof_handler);
  check_this<dim, hp::DoFHandler<dim>>(hp_dof_handler);
}


#define CHECK(EL, deg, dim) \
  {                         \
    FE_##EL<dim> EL(deg);   \
    check(EL, #EL #deg);    \
  }

#define CHECK_SYS1(sub1, N1, dim) \
  {                               \
    FESystem<dim> q(sub1, N1);    \
    check(q, #sub1 #N1);          \
  }

#define CHECK_SYS2(sub1, N1, sub2, N2, dim) \
  {                                         \
    FESystem<dim> q(sub1, N1, sub2, N2);    \
    check(q, #sub1 #N1 #sub2 #N2);          \
  }

#define CHECK_SYS3(sub1, N1, sub2, N2, sub3, N3, dim) \
  {                                                   \
    FESystem<dim> q(sub1, N1, sub2, N2, sub3, N3);    \
    check(q, #sub1 #N1 #sub2 #N2 #sub3 #N3);          \
  }

#define CHECK_ALL(EL, deg) \
  {                        \
    CHECK(EL, deg, 2);     \
    CHECK(EL, deg, 3);     \
  }


int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    logfile;

      CHECK_ALL(Q, 1)
      CHECK_ALL(Q, 2)
      CHECK_ALL(Q, 3)

      CHECK_ALL(DGQ, 0)
      CHECK_ALL(DGQ, 1)
      CHECK_ALL(DGQ, 3)

      CHECK_ALL(DGP, 0)
      CHECK_ALL(DGP, 1)
      CHECK_ALL(DGP, 3)

      CHECK(Nedelec, 0, 2)
      CHECK(Nedelec, 0, 3)

      CHECK(RaviartThomas, 0, 2)
      CHECK(RaviartThomas, 1, 2)
      CHECK(RaviartThomas, 2, 2)

      CHECK(RaviartThomasNodal, 0, 2)
      CHECK(RaviartThomasNodal, 1, 2)
      CHECK(RaviartThomasNodal, 2, 2)
      CHECK(RaviartThomasNodal, 0, 3)
      CHECK(RaviartThomasNodal, 1, 3)
      CHECK(RaviartThomasNodal, 2, 3)

      CHECK_SYS1(FE_Q<2>(1), 3, 2)
      CHECK_SYS1(FE_DGQ<2>(2), 2, 2)
      CHECK_SYS1(FE_DGP<2>(3), 1, 2)

      CHECK_SYS1(FE_Q<3>(1), 3, 3)
      CHECK_SYS1(FE_DGQ<3>(2), 2, 3)
      CHECK_SYS1(FE_DGP<3>(3), 1, 3)

      CHECK_SYS2(FE_Q<2>(1), 3, FE_DGQ<2>(2), 2, 2)
      CHECK_SYS2(FE_DGQ<2>(2), 2, FE_DGP<2>(3), 1, 2)
      CHECK_SYS2(FE_DGP<2>(3), 1, FE_DGQ<2>(2), 2, 2)

      CHECK_SYS3(FE_Q<2>(1), 3, FE_DGP<2>(3), 1, FE_Q<2>(1), 3, 2)
      CHECK_SYS3(FE_DGQ<2>(2), 2, FE_DGQ<2>(2), 2, FE_Q<2>(3), 3, 2)
      CHECK_SYS3(FE_DGP<2>(3), 1, FE_DGP<2>(3), 1, FE_Q<2>(2), 3, 2)

      CHECK_SYS3(FE_DGQ<3>(1), 3, FE_DGP<3>(3), 1, FE_Q<3>(1), 3, 3)

      // systems of systems
      CHECK_SYS3(
        (FESystem<2>(FE_Q<2>(1), 3)), 3, FE_DGQ<2>(0), 1, FE_Q<2>(1), 3, 2)
      CHECK_SYS3(FE_DGQ<2>(3),
                 1,
                 FESystem<2>(FE_DGQ<2>(0), 3),
                 1,
                 FESystem<2>(FE_Q<2>(2), 1, FE_DGQ<2>(0), 1),
                 2,
                 2)

      // systems with Nedelec elements
      CHECK_SYS2(FE_DGQ<2>(3), 1, FE_Nedelec<2>(0), 2, 2)
      CHECK_SYS3(FE_Nedelec<2>(0),
                 1,
                 FESystem<2>(FE_DGQ<2>(1), 2),
                 1,
                 FESystem<2>(FE_Q<2>(2), 1, FE_Nedelec<2>(0), 2),
                 2,
                 2)

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
