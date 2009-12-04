//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Test whether the various assember classes put the right data in the
// right place.

#include "../tests.h"
#include <numerics/mesh_worker.h>
#include <numerics/mesh_worker_assembler.h>
#include <numerics/mesh_worker_loop.h>

#include <base/logstream.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <grid/grid_generator.h>
#include <dofs/dof_tools.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <fe/fe_dgp.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


// Define a class that fills all available entries in the info objects
// with recognizable numbers.
template <int dim>
class Local : public Subscriptor
{
  public:
    typedef typename MeshWorker::IntegrationWorker<dim>::CellInfo CellInfo;
    typedef typename MeshWorker::IntegrationWorker<dim>::FaceInfo FaceInfo;

    void cell(CellInfo& info) const;
    void bdry(FaceInfo& info) const;
    void face(FaceInfo& info1, FaceInfo& info2) const;

    bool cells;
    bool faces;
};

// Fill local structures with the following format:
// 1. Residuals: CC.DDD
// 2. Matrices:  CC.DDDEEE

//   CC : 2 digits number of cell
//
//   DDD: degree of freedom of test function in cell, 1 digit block, 2
//        digits number in block
//   EEE: degree of freedom of trial function in cell, same format as
//        DDD

template <int dim>
void
Local<dim>::cell(CellInfo& info) const
{
  if (!cells) return;
  const unsigned int cell = info.cell->user_index();

  deallog << "Cell " << std::setw(2) << cell << ':';
  for (unsigned int i=0;i<info.indices.size();++i)
    deallog << ' ' << info.indices[i];
  deallog << std::endl;

				   // Fill local residuals
  for (unsigned int k=0;k<info.R.size();++k)
    for (unsigned int b=0;b<info.R[k].n_blocks();++b)
      for (unsigned int i=0;i<info.R[k].block(b).size();++i)
	{
	  const double x = cell + 0.1 * b + 0.001 * i;
	  info.R[k].block(b)(i) = x;
	}

  for (unsigned int k=0;k<info.M1.size();++k)
    {
      const unsigned int block_row = info.M1[k].row;
      const unsigned int block_col = info.M1[k].column;
      FullMatrix<double>& M1 = info.M1[k].matrix;
      for (unsigned int i=0;i<M1.m();++i)
	for (unsigned int j=0;j<M1.n();++j)
	  {
	    double x = .1 * block_row + .001 * i;
	    x = .1 * block_col + .001 * j + .001 * x;
	    M1(i,j) = cell + x;
	  }
    }
}


template <int dim>
void
Local<dim>::bdry(FaceInfo& info) const
{
  const unsigned int cell = info.cell->user_index();
  deallog << "Bdry " << std::setw(2) << cell;
  deallog << std::endl;
}


template <int dim>
void
Local<dim>::face(FaceInfo& info1, FaceInfo& info2) const
{
  const unsigned int cell1 = info1.cell->user_index();
  const unsigned int cell2 = info2.cell->user_index();
  deallog << "Face " << cell1 << '|' << cell2;
  deallog << std::endl;
}


template <int dim, class WORKER>
void
assemble(const DoFHandler<dim>& dof_handler, WORKER& worker)
{
  const FiniteElement<dim>& fe = dof_handler.get_fe();
  MappingQ1<dim> mapping;

  MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
  info_box.initialize(worker, fe, mapping);

  MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(), info_box, worker);
}


template <int dim>
void
test_simple(MGDoFHandler<dim>& mgdofs)
{
  SparsityPattern pattern;
  SparseMatrix<double> matrix;
  Vector<double> v;

  const DoFHandler<dim>& dofs = mgdofs;
  const FiniteElement<dim>& fe = dofs.get_fe();
  pattern.reinit (dofs.n_dofs(), dofs.n_dofs(),
		  (GeometryInfo<dim>::faces_per_cell
		   *GeometryInfo<dim>::max_children_per_face+1)*fe.dofs_per_cell);
  DoFTools::make_flux_sparsity_pattern (dofs, pattern);
  pattern.compress();
  matrix.reinit (pattern);
  v.reinit (dofs.n_dofs());

  Local<dim> local;
  local.cells = true;
  local.faces = false;

  MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> >, Local<dim> >
    integrator(local);
  integrator.initialize_gauss_quadrature(1, 1, 1);
  integrator.initialize(matrix, v);
  integrator.boundary_fluxes = local.faces;
  integrator.interior_fluxes = local.faces;

  assemble(dofs, integrator);
  for (unsigned int i=0;i<v.size();++i)
    deallog << ' ' << std::setprecision(3) << v(i);
  deallog << std::endl;

  deallog << std::setprecision(6);
  matrix.print(deallog, true);
}


template<int dim>
void
test(const FiniteElement<dim>& fe)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);
  tr.begin(1)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  tr.begin(2)->set_refine_flag();
  tr.execute_coarsening_and_refinement();
//  tr.refine_global(1);
  deallog << "Triangulation levels";
  for (unsigned int l=0;l<tr.n_levels();++l)
    deallog << ' ' << l << ':' << tr.n_cells(l);
  deallog << std::endl;

  unsigned int cn = 0;
  for (typename Triangulation<dim>::cell_iterator cell = tr.begin();
       cell != tr.end(); ++cell, ++cn)
    cell->set_user_index(cn);

  MGDoFHandler<dim> dofs(tr);
  dofs.distribute_dofs(fe);
  deallog << "DoFHandler " << dofs.n_dofs() << " levels";
  for (unsigned int l=0;l<tr.n_levels();++l)
    deallog << ' ' << l << ':' << dofs.n_dofs(l);
  deallog << std::endl;

  test_simple(dofs);
}


int main ()
{
  std::ofstream logfile ("mesh_worker_01/output");
  logfile << std::setprecision (2);
  logfile << std::fixed;
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  std::vector<boost::shared_ptr<FiniteElement<2> > > fe2;
  fe2.push_back(boost::shared_ptr<FiniteElement<2> >(new  FE_DGP<2>(1)));
  fe2.push_back(boost::shared_ptr<FiniteElement<2> >(new  FE_Q<2>(1)));

  for (unsigned int i=0;i<fe2.size();++i)
    test(*fe2[i]);
}
