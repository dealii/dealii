// ---------------------------------------------------------------------
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



// a problem uncovered by Baerbel Janssen in that
// DoFTools::make_flux_sparsity_pattern aborted in 1d with adaptively refined
// meshes and hp DoFHandlers. this actually uncovered all sorts of other
// problems that led to a long sequence of assertions triggered everytime one
// of them was fixed. this test cumulatively makes sure everything is ok

// while there, also test the same code in 2d (3d appears to take too long)


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <fstream>

namespace Step
{
  using namespace dealii;

  template <int dim>
  class Problem : public Subscriptor
  {
  public:
    Problem ();
    ~Problem ();

    void refine_mesh ();
    void setup_system ();

  private:
    Triangulation<dim>    triangulation;

    hp::DoFHandler<dim>     dof_handler;
    hp::FECollection<dim>     fe_collection;
    hp::QCollection<dim>      quadrature_collection;
    hp::QCollection<dim-1>    face_quadrature_collection;

    ConstraintMatrix      constraints;

    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  system_matrix;
  };



  template <int dim>
  Problem<dim>::Problem ()
    :
    dof_handler (triangulation)
  {
    GridGenerator::hyper_cube(triangulation);
    for (unsigned int degree=1; degree<=7; ++degree)
      {
        fe_collection.push_back (FE_DGQ<dim>(degree));
        quadrature_collection.push_back (QGauss<dim>(degree+1));
        face_quadrature_collection.push_back (QGauss<dim-1>(degree+1));
      }
  }


  template <int dim>
  Problem<dim>::~Problem ()
  {
    dof_handler.clear ();
  }

  template <int dim>
  void Problem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe_collection);

    constraints.clear ();
    constraints.close ();


    CompressedSetSparsityPattern csp (dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
    DoFTools::make_flux_sparsity_pattern (dof_handler, csp, constraints, false);
    sparsity_pattern.copy_from (csp);
    system_matrix.reinit (sparsity_pattern);

    deallog << "   Number of active cells:       "
            << triangulation.n_active_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl
            << "   Number of constraints       : "
            << constraints.n_constraints()
            << std::endl;
    deallog << "nnz=" << sparsity_pattern.n_nonzero_elements() << std::endl;
  }


  template <int dim>
  void Problem<dim>::refine_mesh ()
  {
    dof_handler.begin_active()->set_refine_flag();
    triangulation.execute_coarsening_and_refinement ();
  }
}


template <int dim>
void test ()
{
  using namespace dealii;
  using namespace Step;

  Problem<dim> problem;

  for (unsigned int cycle=0; cycle<3; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;
      problem.setup_system();
      problem.refine_mesh();
    }
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<1> ();
  test<2> ();

  deallog << "OK" << std::endl;
}
