// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// test meshworker LoopControl

#include "../tests.h"
#include <deal.II/meshworker/loop.h>
#include <deal.II/meshworker/assembler.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/filtered_iterator.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


template <int dim>
class myIntegrator: public dealii::MeshWorker::LocalIntegrator<dim>
{
public:
  typedef MeshWorker::IntegrationInfo<dim> CellInfo;

  void cell(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void boundary(MeshWorker::DoFInfo<dim> &dinfo, CellInfo &info) const;
  void face(MeshWorker::DoFInfo<dim> &dinfo1, MeshWorker::DoFInfo<dim> &dinfo2,
            CellInfo &info1, CellInfo &info2) const;


};

template <int dim>
void
myIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  deallog << "C " << info.cell->id() << std::endl;
}


template <int dim>
void
myIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &info, CellInfo &) const
{
  //deallog << "B cell = " << info.cell->id() << " face = " << info.face_number << std::endl;
}


template <int dim>
void
myIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &info1, MeshWorker::DoFInfo<dim> &info2,
                 CellInfo &, CellInfo &) const
{
  deallog << "F cell1 = " << info1.cell->id()
      << " face = " << info1.face_number
      << " cell2 = " << info2.cell->id()
      << " face2 = " << info2.face_number
      << std::endl;
}


class DoNothingAssembler
{
  public:
    template <class DOFINFO>
    void initialize_info(DOFINFO &info, bool face) const {}
    template<class DOFINFO>
    void assemble(const DOFINFO &info){}
    template<class DOFINFO>
    void assemble(const DOFINFO &info1,
                  const DOFINFO &info2) {}

};

template <int dim>
void
test_simple(DoFHandler<dim> &dofs, MeshWorker::LoopControl &lctrl)
{
  myIntegrator<dim> local;
  DoNothingAssembler assembler;
  MeshWorker::IntegrationInfoBox<dim> info_box;

  MeshWorker::DoFInfo<dim> dof_info(dofs.block_info());

//  integration_loop(ITERATOR begin,
//                typename identity<ITERATOR>::type end,
//                DOFINFO &dinfo,
//                INFOBOX &info,
//                const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &cell_worker,
//                const std_cxx11::function<void (DOFINFO &, typename INFOBOX::CellInfo &)> &boundary_worker,
//                const std_cxx11::function<void (DOFINFO &, DOFINFO &,
//                                                typename INFOBOX::CellInfo &,
//                                                typename INFOBOX::CellInfo &)> &face_worker,
//                ASSEMBLER &assembler,
//                const LoopControl &lctrl)
//


  MeshWorker::integration_loop<dim, dim, typename DoFHandler<dim>::active_cell_iterator, DoNothingAssembler>
    (dofs.begin_active(), dofs.end(),
        dof_info, info_box,
        local,
        assembler,
        lctrl);

//  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
//    (dofs.begin_active(), dofs.end(),
//   dof_info, info_box,
//       std_cxx11::bind (&Integrator<dim>::cell, local, std_cxx11::_1, std_cxx11::_2),
//   std_cxx11::bind (&Integrator<dim>::bdry, local, std_cxx11::_1, std_cxx11::_2),
//   std_cxx11::bind (&Integrator<dim>::face, local, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
//     local,
//     lctrl);
}

std::string id_to_string(const CellId &id)
{
  std::ostringstream ss;
  ss << id;
  return ss.str();
}

template<int dim>
void test_loop(DoFHandler<dim> &dofs, MeshWorker::LoopControl &lctrl)
{
  deallog << "* own_cells=" << lctrl.own_cells
      << " ghost_cells=" << lctrl.ghost_cells
      << " own_faces=" << lctrl.own_faces
      << " faces_to_ghost=" << lctrl.faces_to_ghost
      << std::endl;
  test_simple(dofs, lctrl);
}

template<int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                               Triangulation<dim>::none/*,
                   parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy*/);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  unsigned int myid=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (myid==0)
    tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  FE_DGP<dim> fe(0);
  
  DoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);

  dofs.initialize_local_block_info();
  deallog << "DoFHandler ndofs=" << dofs.n_dofs() << std::endl;

  MeshWorker::LoopControl lctrl;

  lctrl.own_cells = false;
  lctrl.ghost_cells = false;

  lctrl.own_faces = MeshWorker::LoopControl::one;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::never;
  test_loop(dofs, lctrl);

  lctrl.own_faces = MeshWorker::LoopControl::never;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::one;
  test_loop(dofs, lctrl);

  lctrl.own_faces = MeshWorker::LoopControl::never;
  lctrl.faces_to_ghost = MeshWorker::LoopControl::both;
  test_loop(dofs, lctrl);

}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  MPILogInitAll log;

  test<2>();
}
