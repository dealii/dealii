// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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



// this function tests the correctness of the copy operation in the
// matrix-free class on matrix-vector products similar to
// thread_correctness.cc

#include "../tests.h"
#include <deal.II/base/function.h>

std::ofstream logfile("output");

#include "matrix_vector_common.h"


template <int dim, int fe_degree, typename number>
void sub_test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);
  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active (),
  endc = tria.end();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();
  for (; cell!=endc; ++cell)
    if (cell->center().norm()<0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels()-1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active ();

  FE_Q<dim> fe (fe_degree);
  DoFHandler<dim> dof (tria);
  deallog << "Testing " << fe.get_name() << std::endl;

  // run test for several different meshes
  for (unsigned int i=0; i<8-2*dim; ++i)
    {
      cell = tria.begin_active ();
      unsigned int counter = 0;
      for (; cell!=endc; ++cell, ++counter)
        if (counter % (9-i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();

      dof.distribute_dofs(fe);
      ConstraintMatrix constraints;
      DoFTools::make_hanging_node_constraints(dof, constraints);
      VectorTools::interpolate_boundary_values (dof, 0, ZeroFunction<dim>(),
                                                constraints);
      constraints.close();

      //std::cout << "Number of cells: " << dof.get_tria().n_active_cells() << std::endl;
      //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
      //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

      MatrixFree<dim,number> mf_data;
      {
        const QGauss<1> quad (fe_degree+1);
        mf_data.reinit (dof, constraints, quad,
                        typename MatrixFree<dim,number>::AdditionalData(MPI_COMM_SELF,MatrixFree<dim,number>::AdditionalData::none));
      }

      MatrixFreeTest<dim,fe_degree,number> mf_ref (mf_data);

      Vector<number> in_dist (dof.n_dofs());
      Vector<number> out_ref (in_dist), out_copy (in_dist);
      MatrixFree<dim,number> mf_copy;
      mf_copy.copy_from (mf_data);
      MatrixFreeTest<dim,fe_degree,number> copied (mf_copy);

      for (unsigned int i=0; i<dof.n_dofs(); ++i)
        {
          if (constraints.is_constrained(i))
            continue;
          const double entry = Testing::rand()/(double)RAND_MAX;
          in_dist(i) = entry;
        }

      mf_ref.vmult (out_ref, in_dist);
      copied.vmult (out_copy, in_dist);

      out_copy -= out_ref;
      double diff_norm = out_copy.linfty_norm();
      deallog << "Error in copied MF: " << diff_norm
              << std::endl;
    }
  deallog << std::endl;
}


template <int dim, int fe_degree>
void test ()
{
  deallog << "Test doubles" << std::endl;
  sub_test<dim,fe_degree,double>();
  deallog.threshold_double(1.e-6);
  deallog << "Test floats" << std::endl;
  sub_test<dim,fe_degree,float>();
}
