//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2009, 2010 by the deal.II authors
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

#include <base/std_cxx1x/function.h>
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

    void cell(MeshWorker::DoFInfo<dim>& dinfo, CellInfo& info) const;
    void bdry(MeshWorker::DoFInfo<dim>& dinfo, FaceInfo& info) const;
    void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
	      FaceInfo& info1, FaceInfo& info2) const;

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
Local<dim>::cell(MeshWorker::DoFInfo<dim>& info, CellInfo&) const
{
  if (!cells) return;
  const unsigned int cell = info.cell->user_index();

  deallog << "Cell " << std::setw(2) << cell << ':';
  for (unsigned int i=0;i<info.indices.size();++i)
    deallog << ' ' << info.indices[i];
  deallog << std::endl;

				   // Fill local residuals
  for (unsigned int k=0;k<info.n_vectors();++k)
    for (unsigned int b=0;b<info.vector(k).n_blocks();++b)
      for (unsigned int i=0;i<info.vector(k).block(b).size();++i)
	{
	  const double x = cell + 0.1 * b + 0.001 * i;
	  info.vector(k).block(b)(i) = x;
	}

  for (unsigned int k=0;k<info.n_matrices();++k)
    {
      const unsigned int block_row = info.matrix(k).row;
      const unsigned int block_col = info.matrix(k).column;
      FullMatrix<double>& M1 = info.matrix(k).matrix;
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
Local<dim>::bdry(MeshWorker::DoFInfo<dim>& info, FaceInfo&) const
{
  const unsigned int cell = info.cell->user_index();
  deallog << "Bdry " << std::setw(2) << cell;
  deallog << std::endl;
}


template <int dim>
void
Local<dim>::face(MeshWorker::DoFInfo<dim>& info1, MeshWorker::DoFInfo<dim>& info2,
		 FaceInfo&, FaceInfo&) const
{
  const unsigned int cell1 = info1.cell->user_index();
  const unsigned int cell2 = info2.cell->user_index();
  deallog << "Face " << cell1 << '|' << cell2;
  deallog << std::endl;
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

  MeshWorker::IntegrationWorker<dim> integrator;
  integrator.initialize_gauss_quadrature(1, 1, 1);
  integrator.boundary_fluxes = local.faces;
  integrator.interior_fluxes = local.faces;

  MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double> >
    assembler;
  assembler.initialize(matrix, v);

  {
    MappingQ1<dim> mapping;

    MeshWorker::IntegrationInfoBox<dim> info_box;
    info_box.initialize(integrator, fe, mapping);
    MeshWorker::DoFInfo<dim> dof_info(dofs);
    
    MeshWorker::loop
      <MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
      (dofs.begin_active(), dofs.end(),
       dof_info, info_box,
       std_cxx1x::bind (&Local<dim>::cell, local, _1, _2),
       std_cxx1x::bind (&Local<dim>::bdry, local, _1, _2),
       std_cxx1x::bind (&Local<dim>::face, local, _1, _2, _3, _4),
       assembler, true);
  }

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
