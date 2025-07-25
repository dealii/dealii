// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test whether the various assembler classes put the right data in the
// right place.

// This test is a modification of the test with the same name in ../integrators

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <functional>

#include "../tests.h"



// Define a class that fills all available entries in the info objects
// with recognizable numbers.
template <int dim>
class Local : public EnableObserverPointer
{
public:
  using CellInfo = MeshWorker::IntegrationInfo<dim>;

  void
  cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  bdry(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void
  face(MeshWorker::DoFInfo<dim> &dinfo1,
       MeshWorker::DoFInfo<dim> &dinfo2,
       CellInfo                 &info1,
       CellInfo                 &info2) const;

  bool cells;
  bool faces;
};

template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  if (!cells)
    return;
  for (unsigned int k = 0; k < info.n_matrices(); ++k)
    {
      const unsigned int  block_row = info.matrix(k).row;
      const unsigned int  block_col = info.matrix(k).column;
      FullMatrix<double> &M1        = info.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i = 0; i < M1.m(); ++i)
          for (unsigned int j = 0; j < M1.n(); ++j)
            {
              M1(i, j) = 10.;
            }
    }
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  if (!faces)
    return;
  for (unsigned int k = 0; k < info.n_matrices(); ++k)
    {
      const unsigned int  block_row = info.matrix(k).row;
      const unsigned int  block_col = info.matrix(k).column;
      FullMatrix<double> &M1        = info.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i = 0; i < M1.m(); ++i)
          for (unsigned int j = 0; j < M1.n(); ++j)
            {
              M1(i, j) = 1.;
            }
    }
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim> &info1,
                 MeshWorker::DoFInfo<dim> &info2,
                 CellInfo &,
                 CellInfo &) const
{
  if (!faces)
    return;
  for (unsigned int k = 0; k < info1.n_matrices(); ++k)
    {
      const unsigned int  block_row = info1.matrix(k).row;
      const unsigned int  block_col = info1.matrix(k).column;
      FullMatrix<double> &M1        = info1.matrix(k).matrix;
      if (block_row == block_col)
        for (unsigned int i = 0; i < M1.m(); ++i)
          for (unsigned int j = 0; j < M1.n(); ++j)
            {
              info1.matrix(k, false).matrix(i, j) = 1.;
              info2.matrix(k, false).matrix(i, j) = 1.;
              info1.matrix(k, true).matrix(i, j)  = -1.;
              info2.matrix(k, true).matrix(i, j)  = -1.;
            }
    }
}


template <int dim>
void
test_simple(DoFHandler<dim> &dofs, bool faces)
{
  TrilinosWrappers::SparseMatrix matrix;

  const FiniteElement<dim> &fe = dofs.get_fe();
  DynamicSparsityPattern    csp(dofs.n_dofs(), dofs.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dofs, csp);
  matrix.reinit(dofs.locally_owned_dofs(), csp, MPI_COMM_WORLD, true);

  Local<dim> local;
  local.cells = true;
  local.faces = faces;

  MappingQ<dim> mapping(1);

  MeshWorker::IntegrationInfoBox<dim> info_box;
  info_box.initialize_gauss_quadrature(1, 1, 1);
  info_box.initialize_update_flags();
  info_box.initialize(fe, mapping);

  MeshWorker::DoFInfo<dim> dof_info(dofs.block_info());

  MeshWorker::Assembler::MatrixSimple<TrilinosWrappers::SparseMatrix> assembler;
  assembler.initialize(matrix);

  FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> cell(
    IteratorFilters::LocallyOwnedCell(), dofs.begin_active());
  FilteredIterator<typename DoFHandler<dim>::active_cell_iterator> end(
    IteratorFilters::LocallyOwnedCell(), dofs.end());

  MeshWorker::LoopControl lctrl;
  lctrl.cells_first = true;
  lctrl.own_faces   = MeshWorker::LoopControl::one;
  MeshWorker::loop<dim,
                   dim,
                   MeshWorker::DoFInfo<dim>,
                   MeshWorker::IntegrationInfoBox<dim>>(
    cell,
    end,
    dof_info,
    info_box,
    std::bind(
      &Local<dim>::cell, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(
      &Local<dim>::bdry, local, std::placeholders::_1, std::placeholders::_2),
    std::bind(&Local<dim>::face,
              local,
              std::placeholders::_1,
              std::placeholders::_2,
              std::placeholders::_3,
              std::placeholders::_4),
    assembler,
    lctrl);

  matrix.compress(VectorOperation::add);
  matrix.print(deallog.get_file_stream());
}


template <int dim>
void
test(const FiniteElement<dim> &fe)
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none/*,
                   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy*/);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid == 0)
    tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  deallog << "Triangulation levels";
  for (unsigned int l = 0; l < tr.n_levels(); ++l)
    deallog << ' ' << l << ':' << tr.n_cells(l);
  deallog << std::endl;

  unsigned int cn = 0;
  for (typename Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end();
       ++cell, ++cn)
    cell->set_user_index(cn);

  DoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);



  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  if (numprocs == 1) // renumber DoFs
    {
      std::map<std::string, std::vector<types::global_dof_index>> dofmap;
      std::ifstream f(SOURCE_DIR "/mesh_worker_matrix_01/ordering.2");
      while (!f.eof())
        {
          CellId id;
          f >> id;
          if (f.eof())
            break;
          std::vector<types::global_dof_index> &d = dofmap[id.to_string()];
          d.reserve(fe.dofs_per_cell);
          for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
            {
              unsigned int temp;
              f >> temp;
              d.push_back(temp);
            }
        }

      for (typename DoFHandler<dim>::active_cell_iterator cell =
             dofs.begin_active();
           cell != dofs.end();
           ++cell)
        {
          if (!cell->is_locally_owned())
            continue;

          std::vector<types::global_dof_index> &renumbered =
            dofmap[cell->id().to_string()];
          cell->set_dof_indices(renumbered);
        }
    }

  if (0)
    {
      unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

      for (unsigned int i = 0; i < numprocs; ++i)
        {
          MPI_Barrier(MPI_COMM_WORLD);
          if (myid == i)
            {
              std::ofstream                        f("ordering",
                              (myid > 0) ? std::ofstream::app :
                                                                  std::ofstream::out);
              std::vector<types::global_dof_index> local_dof_indices(
                fe.dofs_per_cell);
              for (typename DoFHandler<dim>::active_cell_iterator cell =
                     dofs.begin_active();
                   cell != dofs.end();
                   ++cell)
                {
                  if (!cell->is_locally_owned())
                    continue;

                  f << cell->id() << ' ';
                  cell->get_dof_indices(local_dof_indices);
                  for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                    f << local_dof_indices[i] << ' ';
                  f << std::endl;
                }
              f.close();
            }
        }
    }


  dofs.initialize_local_block_info();
  deallog << "DoFHandler " << dofs.n_dofs() << std::endl;

  test_simple(dofs, false);
  deallog << "now with jump terms" << std::endl;
  test_simple(dofs, true);
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());
  MPILogInitAll log;

  FE_DGP<2>                       p0(0);
  FE_Q<2>                         q1(1);
  FESystem<2, 2>                  sys1(p0, 1, q1, 1);
  std::vector<FiniteElement<2> *> fe2;
  fe2.push_back(&p0);
  //  fe2.push_back(&q1);
  // fe2.push_back(&sys1);

  for (unsigned int i = 0; i < fe2.size(); ++i)
    test(*fe2[i]);
}
