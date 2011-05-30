//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Test whether Assember::Functional adds up correctly.

#include "../tests.h"
#include "empty_info.h"
#include <deal.II/numerics/mesh_worker.h>
#include <deal.II/numerics/mesh_worker_assembler.h>
#include <deal.II/numerics/mesh_worker_loop.h>

#include <deal.II/base/std_cxx1x/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgp.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


// Define a class that fills all available entries in the info objects
// with recognizable numbers.

// Fills the following functionals:
// 0: number of cells
// 1: number of interior faces
// 2: number of boundary faces

const unsigned int n_functionals = 3;

template <int dim>
class Local : public Subscriptor
{
  public:
    typedef EmptyInfo CellInfo;
    
    void cell(MeshWorker::DoFInfo<dim>& dinfo, CellInfo& info) const;
    void bdry(MeshWorker::DoFInfo<dim>& dinfo, CellInfo& info) const;
    void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
	      CellInfo& info1, CellInfo& info2) const;
};


template <int dim>
void
Local<dim>::cell(MeshWorker::DoFInfo<dim>& info, CellInfo&) const
{
  info.value(0) = 1.;
}


template <int dim>
void
Local<dim>::bdry(MeshWorker::DoFInfo<dim>&  info, CellInfo&) const
{
  info.value(2) = 1.;
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim>&  info1, MeshWorker::DoFInfo<dim>& info2,
		 CellInfo&, CellInfo&) const
{
  info1.value(1) = 1./2.;
  info2.value(1) = 1./2.;
}


template <int dim>
void
test_mesh(MGDoFHandler<dim>& mgdofs)
{
  const DoFHandler<dim>& dofs = mgdofs;
  
  Local<dim> local;
  EmptyInfoBox info_box;
  MeshWorker::DoFInfo<dim> dof_info(dofs);
  
  MeshWorker::Assembler::Functional<double> assembler;
  assembler.initialize(n_functionals);
  
  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>
    (dofs.begin_active(), dofs.end(),
     dof_info, info_box,
     std_cxx1x::bind (&Local<dim>::cell, local, std_cxx1x::_1, std_cxx1x::_2),
     std_cxx1x::bind (&Local<dim>::bdry, local, std_cxx1x::_1, std_cxx1x::_2),
     std_cxx1x::bind (&Local<dim>::face, local, std_cxx1x::_1, std_cxx1x::_2, std_cxx1x::_3, std_cxx1x::_4),
     assembler, true);

  deallog << "  Results";
  for (unsigned int i=0;i<n_functionals;++i)
    deallog << '\t' << assembler(i);
  deallog << std::endl;

  assembler.initialize(n_functionals);
  MeshWorker::DoFInfo<dim> mg_dof_info(mgdofs);
  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, EmptyInfoBox>
    (mgdofs.begin(), mgdofs.end(),
     mg_dof_info, info_box,
     std_cxx1x::bind (&Local<dim>::cell, local, std_cxx1x::_1, std_cxx1x::_2),
     std_cxx1x::bind (&Local<dim>::bdry, local, std_cxx1x::_1, std_cxx1x::_2),
     std_cxx1x::bind (&Local<dim>::face, local, std_cxx1x::_1, std_cxx1x::_2, std_cxx1x::_3, std_cxx1x::_4),
     assembler, true);

  deallog << "MGResults";
  for (unsigned int i=0;i<n_functionals;++i)
    deallog << '\t' << assembler(i);
  deallog << std::endl;
}


template<int dim>
void
test(const FiniteElement<dim>& fe)
{
  Triangulation<dim> tr;
  MGDoFHandler<dim> dofs(tr);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  deallog.push("1");
  dofs.distribute_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
  tr.begin(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("2");
  dofs.distribute_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  deallog.push("3");
  dofs.distribute_dofs(fe);
  test_mesh(dofs);
  deallog.pop();
}


int main ()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);
  
  FE_DGP<2> el2(0);
  FE_DGP<3> el3(0);
  
  deallog.push("2D");
  test(el2);
  deallog.pop();
  deallog.push("3D");
  test(el3);
  deallog.pop();
}
