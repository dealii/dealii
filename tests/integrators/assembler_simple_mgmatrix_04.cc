// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2015 by the deal.II authors
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
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/integrators/laplace.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/multigrid/sparse_matrix_collection.h>

using namespace dealii;
using namespace LocalIntegrators;


template <int dim>
class LaplaceMatrix : public MeshWorker::LocalIntegrator<dim>
{
public:
  LaplaceMatrix();
  virtual void cell(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info) const;
  virtual void boundary(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info) const;
  virtual void face(MeshWorker::DoFInfo<dim> &dinfo1, MeshWorker::DoFInfo<dim> &dinfo2,
                    MeshWorker::IntegrationInfo<dim> &info1, MeshWorker::IntegrationInfo<dim> &info2) const;
};


template <int dim>
LaplaceMatrix<dim>::LaplaceMatrix()
{}


template <int dim>
void LaplaceMatrix<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info) const
{
  Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0));
}


template <int dim>
void LaplaceMatrix<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo,
                                  typename MeshWorker::IntegrationInfo<dim> &info) const
{
  const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
  Laplace::nitsche_matrix(dinfo.matrix(0,false).matrix, info.fe_values(0),
                          Laplace::compute_penalty(dinfo, dinfo, deg, deg));
}


template <int dim>
void LaplaceMatrix<dim>::face(
  MeshWorker::DoFInfo<dim> &dinfo1, MeshWorker::DoFInfo<dim> &dinfo2,
  MeshWorker::IntegrationInfo<dim> &info1, MeshWorker::IntegrationInfo<dim> &info2) const
{
  if (info1.fe_values(0).get_fe().conforms(FiniteElementData<dim>::H1))
    return;

  const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();

  if (info1.fe_values(0).get_fe().conforms(FiniteElementData<dim>::Hdiv) &&
      !info1.fe_values(0).get_fe().conforms(FiniteElementData<dim>::Hcurl))
    Laplace::ip_tangential_matrix(dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
                                  dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
                                  info1.fe_values(0), info2.fe_values(0),
                                  Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
  else
    Laplace::ip_matrix(dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
                       dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
                       info1.fe_values(0), info2.fe_values(0),
                       Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
}


template <int dim>
void assemble_mg_matrix(DoFHandler<dim> &dof_handler,
                        MeshWorker::LocalIntegrator<dim> &matrix_integrator, mg::SparseMatrixCollection<double> &mg)
{
  MGConstrainedDoFs mg_constraints;
  mg_constraints.clear();
  mg_constraints.initialize(dof_handler);

  mg.set_zero();

  MappingQGeneric<dim> mapping(1);

  MeshWorker::IntegrationInfoBox<dim> info_box;
  UpdateFlags update_flags = update_values | update_gradients | update_hessians;
  info_box.add_update_flags_all(update_flags);
  info_box.initialize(dof_handler.get_fe(), mapping, &dof_handler.block_info());

  MeshWorker::DoFInfo<dim> dof_info(dof_handler.block_info());

  MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double> > assembler;
  assembler.initialize(mg_constraints);
  assembler.initialize(mg.matrix);
  assembler.initialize_interfaces(mg.matrix_in, mg.matrix_out);
  assembler.initialize_fluxes(mg.matrix_up, mg.matrix_down);

  MeshWorker::integration_loop<dim, dim> (
    dof_handler.begin_mg(), dof_handler.end_mg(),
    dof_info, info_box, matrix_integrator, assembler);

  const unsigned int nlevels = dof_handler.get_triangulation().n_levels();
  for (unsigned int level=0; level<nlevels; ++level)
    {
      for (unsigned int i=0; i<dof_handler.n_dofs(level); ++i)
        if (mg.matrix[level].diag_element(i)==0)
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

  LaplaceMatrix<dim> matrix_integrator;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);
  dof.distribute_mg_dofs(fe);
  dof.initialize_local_block_info();
  mg::SparseMatrixCollection<double> mg;
  mg.resize(0,tr.n_levels()-1);
  mg.reinit(dof);
  assemble_mg_matrix(dof, matrix_integrator, mg);

  for (unsigned int level=0; level<tr.n_levels(); ++level)
    {
      const unsigned int prec = 3;
      const unsigned int wd = 2;

      deallog << "Level " << level << std::endl << "mg" << std::endl;
      mg.matrix[level].print_formatted(deallog.get_file_stream(), prec, false, wd, "0.");
      if (level>0)
        {
          deallog << "in" << std::endl;
          mg.matrix_in[level].print_formatted(deallog.get_file_stream(), prec, false, wd, "0.");
          deallog << "out" << std::endl;
          mg.matrix_out[level].print_formatted(deallog.get_file_stream(), prec, false, wd, "0.");
          deallog << "up" << std::endl;
          mg.matrix_up[level].print_formatted(deallog.get_file_stream(), prec, false, wd, "0.");
          deallog << "down" << std::endl;
          mg.matrix_down[level].print_formatted(deallog.get_file_stream(), prec, false, wd, "0.");
        }
    }
}


int main()
{
  const std::string logname("output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);

  FE_DGP<2> p0(0);
  FE_DGP<2> p1(1);
  FE_RaviartThomas<2> rt0(0);
  FE_RaviartThomas<2> rt1(1);
  FE_Q<2> q2(2);

  FESystem<2> sys1(p0, 2, p0, 1);
  FESystem<2> sys2(rt0, 1, p0, 1);
  FESystem<2> sys3(rt1, 1, p1, 1);
  FESystem<2> sys4(q2, 1, p0, 1);

  test(sys1);
  test(sys2);
  test(sys3);
  test(sys4);
}
