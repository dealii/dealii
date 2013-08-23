// ---------------------------------------------------------------------
// $Id$
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


/**
 * @file Test whether 
 */

#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;

template <int dim>
class LaplaceMatrix : public MeshWorker::LocalIntegrator<dim>
{
public:
  LaplaceMatrix();
  virtual void cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void boundary(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
		    MeshWorker::IntegrationInfo<dim>& info1, MeshWorker::IntegrationInfo<dim>& info2) const;
};


template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix()
{}


template <int dim>
void LaplaceMatrix<dim>::cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const
{
  const unsigned int n_blocks = dinfo.block_info->local().size();
  AssertDimension (dinfo.n_matrices(), n_blocks*n_blocks);
  for(unsigned int b=0; b<n_blocks; ++b)
    Laplace::cell_matrix(dinfo.matrix(b+b*n_blocks,false).matrix, info.fe_values(b));
}


template <int dim>
void LaplaceMatrix<dim>::boundary(MeshWorker::DoFInfo<dim>& dinfo,
				       typename MeshWorker::IntegrationInfo<dim>& info) const
{
  const unsigned int n_blocks = dinfo.block_info->local().size();
  for(unsigned int b=0; b<n_blocks; ++b)
    {
      const unsigned int deg = info.fe_values(b).get_fe().tensor_degree();
      Laplace::nitsche_matrix(dinfo.matrix(b+b*n_blocks,false).matrix, info.fe_values(b),
			      Laplace::compute_penalty(dinfo, dinfo, deg, deg));
    } 
}


template <int dim>
void LaplaceMatrix<dim>::face(
  MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
  MeshWorker::IntegrationInfo<dim>& info1, MeshWorker::IntegrationInfo<dim>& info2) const
{
  const unsigned int n_blocks = dinfo1.block_info->local().size();
  for (unsigned int b=0;b<n_blocks;++b)
    {
      const unsigned int diag = b+b*n_blocks;
      if (info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::H1))
	continue;
      
      const unsigned int deg = info1.fe_values(b).get_fe().tensor_degree();
      
      if (info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::Hdiv) &&
	  !info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::Hcurl))
	Laplace::ip_tangential_matrix(dinfo1.matrix(diag,false).matrix, dinfo1.matrix(diag,true).matrix, 
				      dinfo2.matrix(diag,true).matrix, dinfo2.matrix(diag,false).matrix,
				      info1.fe_values(b), info2.fe_values(b),
				      Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
      else
	Laplace::ip_matrix(dinfo1.matrix(diag,false).matrix, dinfo1.matrix(diag,true).matrix, 
			   dinfo2.matrix(diag,true).matrix, dinfo2.matrix(diag,false).matrix,
			   info1.fe_values(b), info2.fe_values(b),
			   Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
    }
}

template <int dim>
class LaplaceMatrix : public MeshWorker::LocalIntegrator<dim>
{
public:
  LaplaceMatrix();
  virtual void cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void boundary(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
		    MeshWorker::IntegrationInfo<dim>& info1, MeshWorker::IntegrationInfo<dim>& info2) const;
};


template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix()
{}


template <int dim>
void LaplaceMatrix<dim>::cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const
{
  const unsigned int n_blocks = dinfo.block_info->local().size();
  AssertDimension (dinfo.n_matrices(), n_blocks*n_blocks);
  for(unsigned int b=0; b<n_blocks; ++b)
    Laplace::cell_matrix(dinfo.matrix(b+b*n_blocks,false).matrix, info.fe_values(b));
}


template <int dim>
void LaplaceMatrix<dim>::boundary(MeshWorker::DoFInfo<dim>& dinfo,
				       typename MeshWorker::IntegrationInfo<dim>& info) const
{
  const unsigned int n_blocks = dinfo.block_info->local().size();
  for(unsigned int b=0; b<n_blocks; ++b)
    {
      const unsigned int deg = info.fe_values(b).get_fe().tensor_degree();
      Laplace::nitsche_matrix(dinfo.matrix(b+b*n_blocks,false).matrix, info.fe_values(b),
			      Laplace::compute_penalty(dinfo, dinfo, deg, deg));
    } 
}


template <int dim>
void LaplaceMatrix<dim>::face(
  MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
  MeshWorker::IntegrationInfo<dim>& info1, MeshWorker::IntegrationInfo<dim>& info2) const
{
  const unsigned int n_blocks = dinfo1.block_info->local().size();
  for (unsigned int b=0;b<n_blocks;++b)
    {
      const unsigned int diag = b+b*n_blocks;
      if (info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::H1))
	continue;
      
      const unsigned int deg = info1.fe_values(b).get_fe().tensor_degree();
      
      if (info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::Hdiv) &&
	  !info1.fe_values(b).get_fe().conforms(FiniteElementData<dim>::Hcurl))
	Laplace::ip_tangential_matrix(dinfo1.matrix(diag,false).matrix, dinfo1.matrix(diag,true).matrix, 
				      dinfo2.matrix(diag,true).matrix, dinfo2.matrix(diag,false).matrix,
				      info1.fe_values(b), info2.fe_values(b),
				      Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
      else
	Laplace::ip_matrix(dinfo1.matrix(diag,false).matrix, dinfo1.matrix(diag,true).matrix, 
			   dinfo2.matrix(diag,true).matrix, dinfo2.matrix(diag,false).matrix,
			   info1.fe_values(b), info2.fe_values(b),
			   Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
    }
}

template <int dim>
void assemble_mg_matrix(FiniteElement<dim> &fe, Mapping<dim>& mapping, DoFHandler<dim>& dof_handler,
   MeshWorker::LocalIntegrator<dim>& matrix_integrator, mg::SparseMatrixCollection<double>& mg)
{
  mg.set_zero();
  
  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags update_flags = update_values | update_gradients | update_hessians;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize(fe, mapping, &dof_handler.block_info());

  MeshWorker::DoFInfo<dim> dof_info(dof_handler.block_info());
  
  MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(mg_constraints);
  assembler.initialize(mg.matrix);
  assembler.initialize_interfaces(mg.matrix_in, mg.matrix_out);
  assembler.initialize_fluxes(mg.matrix_up, mg.matrix_down);
  
  MeshWorker::integration_loop<dim, dim> (
    mg_dof_handler.begin(), mg_dof_handler.end(),
    dof_info, info_box, matrix_integrator, assembler);

  const unsigned int nlevels = triangulation.n_levels();
  for (unsigned int level=0;level<nlevels;++level)
  {
    for(unsigned int i=0; i<mg_dof_handler.n_dofs(level); ++i)
      if(mg.matrix[level].diag_element(i)==0)
        mg.matrix[level].set(i,i,1.);
  }
}


template <int dim>
void test(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  std::vector<unsigned int> repititions (dim, 1);
  repititions[0] = 2;
  const Point<dim> p1 = (dim == 1 ? Point<dim> (-1.) : (dim == 2 ? Point<dim>(-1.,-1.) : Point<dim> (-1.,-1.,-1.)));
  const Point<dim> p2 = (dim == 1 ? Point<dim> (1.) : (dim == 2 ? Point<dim>(1.,1.) : Point<dim> (1.,1.,1.)));
  GridGenerator::subdivided_hyper_rectangle(tr, repititions, p1, p2);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

//  MGDoFHandler<dim> dof(tr);
//  dof.distribute_dofs(fe);
//  dof.initialize_local_block_info();
//  DoFRenumbering::component_wise(dof);
//
//  deallog << "DoFs " << dof.n_dofs() << std::endl;
//
//  typename MGDoFHandler<dim>::cell_iterator cell = dof.begin_active();
//  typename MGDoFHandler<dim>::face_iterator face = cell->face(1);
//  typename MGDoFHandler<dim>::cell_iterator neighbor = cell->neighbor(1);
//
//  MGLevelObject<SparsityPattern> sparsity(0, tr.n_levels()-1);
//  MGLevelObject<SparseMatrix<double> > matrix(0, tr.n_levels()-1);
//
//  for (unsigned int level=0; level<tr.n_levels(); ++level)
//    {
//      CompressedSparsityPattern csp(dof.n_dofs(level),dof.n_dofs(level));
//      MGTools::make_flux_sparsity_pattern(dof, csp, level);
//      sparsity[level].copy_from(csp);
//      matrix[level].reinit(sparsity[level]);
//    }
//
//  MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > ass;
//  ass.initialize(matrix);
//  ass.initialize_local_blocks(dof.block_info().local());
//  MeshWorker::DoFInfo<dim> info(dof.block_info());
//  ass.initialize_info(info, false);
//  MeshWorker::DoFInfo<dim> infon(dof.block_info());
//  ass.initialize_info(infon, true);
//
//  deallog << "cell" << std::endl;
//  info.reinit(cell);
//  fill_matrices(info, false);
//  ass.assemble(info);
//  matrix[1].print_formatted(deallog.get_file_stream(), 0, false, 6);
//  matrix[1] = 0.;
//
//  deallog << "face" << std::endl;
//  ass.initialize_info(info, true);
//  info.reinit(cell, face, 1);
//  infon.reinit(neighbor, neighbor->face(0), 0);
//  fill_matrices(info, true);
//  fill_matrices(infon, true);
//  ass.assemble(info, infon);
//  matrix[1].print_formatted(deallog.get_file_stream(), 0, false, 6);
}

int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  FE_DGP<2> p0(0);
  FE_DGP<2> p1(1);
  FE_RaviartThomas<2> rt0(0);
  FE_Q<2> q2(2);

  FESystem<2> sys1(p0, 2, p1, 1);
  FESystem<2> sys2(p0, 2, rt0, 1);
  FESystem<2> sys3(rt0, 1, p0, 2);
  FESystem<2> sys4(p1, 2, q2, 2);
  FESystem<2> sys5(q2, 2, p1, 2);

  test(sys1);
//  test(sys2);
//  test(sys3);
//  test(sys4);
//  test(sys5);
}
