//----------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

/**
 * @file Test initialization of Assembler::MatrixSimple and
 * DoFInfo
 */

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/simple.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;

template <class DOFINFO, class MATRIX>
void test(DOFINFO& info, MeshWorker::Assembler::MatrixSimple<MATRIX>& ass)
{  
  ass.initialize_info(info, false);
  deallog << "No faces" << std::endl;
  info.print_debug(deallog);
  
  ass.initialize_info(info, true);
  deallog << "With faces" << std::endl;
  info.print_debug(deallog);
}

int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);
  
  Triangulation<2,2> tr;
  GridGenerator::hyper_cube(tr);
  FE_DGP<2,2> fe1(1);
  FE_DGP<2,2> fe2(2);
  FE_DGP<2,2> fe3(3);
  FE_DGP<2,2> fe5(5);
  FESystem<2,2> fes1(fe3,1,fe5,1,fe1,1);
  FESystem<2,2> fes2(fe3,1,fe5,1,fe1,1, fe2,1);

  DoFHandler<2,2> dof1(tr);
  dof1.distribute_dofs(fes1);
  DoFHandler<2,2> dof2(tr);
  dof2.distribute_dofs(fes2);    
  dof1.initialize_local_block_info();
  dof2.initialize_local_block_info();
  MeshWorker::DoFInfo<2,2,double> info1(dof1);
  MeshWorker::DoFInfo<2,2,double> info1b(dof1.block_info());
  MeshWorker::DoFInfo<2,2,double> info2b(dof2.block_info());

  MeshWorker::Assembler::MatrixSimple<FullMatrix<double> > ass1;
  deallog.push("Single block");
  test(info1, ass1);
  deallog.pop();
  deallog.push("Multiple blocks");
  test(info1b, ass1);
  deallog.pop();
  deallog.push("More blocks");
  test(info2b, ass1);  
}
