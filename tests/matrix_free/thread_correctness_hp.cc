// ---------------------------------------------------------------------
//
// Copyright (C) 2013-2014 by the deal.II authors
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



// this function tests the correctness of the implementation of parallel
// matrix free matrix-vector products for hp elements by comparing to the
// serial version

#include "../tests.h"

std::ofstream logfile("output");

#include "create_mesh.h"
#include "matrix_vector_common.h"
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/base/function.h>
#include <deal.II/base/template_constraints.h>



template <int dim, typename Number>
class MatrixFreeTestHP
{
public:
  MatrixFreeTestHP(const MatrixFree<dim,Number> &data_in):
    data (data_in)
  {};

  void
  local_apply (const MatrixFree<dim,Number> &data,
               Vector<Number> &dst,
               const Vector<Number> &src,
               const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    // Ask MatrixFree for cell_range for different
    // orders
    std::pair<unsigned int,unsigned int> subrange_deg;
#define CALL_METHOD(degree)                                             \
    subrange_deg = data.create_cell_subrange_hp(cell_range, degree);    \ 
    if (subrange_deg.second > subrange_deg.first)                       \
      helmholtz_operator<dim,degree,Vector<Number> > (data, dst, src, subrange_deg)

    CALL_METHOD(1);
    CALL_METHOD(2);
    CALL_METHOD(3);
    CALL_METHOD(4);
    CALL_METHOD(5);
    CALL_METHOD(6);
    CALL_METHOD(7);

#undef CALL_METHOD
  }

  void vmult (Vector<Number>       &dst,
              const Vector<Number> &src) const
  {
    dst = 0;
    data.cell_loop (&MatrixFreeTestHP<dim,Number>::local_apply, this, dst, src);
  };

private:
  const MatrixFree<dim,Number> &data;
};



template <int dim, typename number>
void do_test (const unsigned int parallel_option)
{
  Triangulation<dim> tria;
  create_mesh (tria);
  // larger mesh in release mode
#ifndef DEBUG
  tria.refine_global(2);
#endif

  // refine a few cells
  for (unsigned int i=0; i<11-3*dim; ++i)
    {
      typename Triangulation<dim>::active_cell_iterator
      cell = tria.begin_active (),
      endc = tria.end();
      for (; cell!=endc; ++cell)
        if (Testing::rand() % (7-i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  const unsigned int max_degree = 9-2*dim;

  hp::FECollection<dim>    fe_collection;
  hp::QCollection<1>       quadrature_collection_mf;

  for (unsigned int deg=1; deg<=max_degree; ++deg)
    {
      fe_collection.push_back (FE_Q<dim>(QGaussLobatto<1>(deg+1)));
      quadrature_collection_mf.push_back (QGauss<1>(deg+1));
    }

  hp::DoFHandler<dim> dof(tria);
  // set the active FE index in a random order
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof.begin_active(),
    endc = dof.end();
    for (; cell!=endc; ++cell)
      {
        const unsigned int fe_index = Testing::rand() % max_degree;
        cell->set_active_fe_index (fe_index);
      }
  }

  // setup DoFs
  dof.distribute_dofs(fe_collection);
  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof,
                                           constraints);
  VectorTools::interpolate_boundary_values (dof,
                                            0,
                                            ZeroFunction<dim>(),
                                            constraints);
  constraints.close ();

  //std::cout << "Number of cells: " << dof.get_tria().n_active_cells() << std::endl;
  //std::cout << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;
  //std::cout << "Number of constraints: " << constraints.n_constraints() << std::endl;

  // set up reference MatrixFree
  MatrixFree<dim,number> mf_data;
  typename MatrixFree<dim,number>::AdditionalData data;
  data.tasks_parallel_scheme =
    MatrixFree<dim,number>::AdditionalData::none;
  mf_data.reinit (dof, constraints, quadrature_collection_mf, data);
  MatrixFreeTestHP<dim,number> mf (mf_data);

  // test different block sizes, starting from
  // auto setting (= 0)
  for (unsigned int block_size = 0; block_size < 5; ++block_size)
    {
      deallog.push ("blk_" + Utilities::int_to_string(block_size,1));
      MatrixFree<dim,number> mf_data_par;
      if (parallel_option == 0)
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::partition_partition;
          deallog << "Parallel option partition/partition" << std::endl;
        }
      else if (parallel_option == 1)
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::partition_color;
          deallog << "Parallel option partition/color" << std::endl;
        }
      else
        {
          data.tasks_parallel_scheme =
            MatrixFree<dim,number>::AdditionalData::color;
          deallog << "Parallel option partition/color" << std::endl;
        }
      data.tasks_block_size = 1;
      mf_data_par.reinit (dof, constraints, quadrature_collection_mf, data);
      MatrixFreeTestHP<dim,number> mf_par(mf_data_par);

      // fill a right hand side vector with random
      // numbers in unconstrained degrees of freedom
      Vector<number> src (dof.n_dofs());
      Vector<number> result_ref(src), result_mf (src);

      for (unsigned int i=0; i<dof.n_dofs(); ++i)
        {
          if (constraints.is_constrained(i) == false)
            src(i) = (double)Testing::rand()/RAND_MAX;
        }

      // now perform 30 matrix-vector products in
      // parallel and check their correctness (take
      // many of them to make sure that we hit an
      // error)
      mf.vmult (result_ref, src);
      deallog << "Norm of difference: ";
      for (unsigned int i=0; i<50; ++i)
        {
          mf_par.vmult (result_mf, src);
          result_mf -= result_ref;
          double diff_norm = result_mf.linfty_norm()/result_ref.linfty_norm();
          deallog << diff_norm << "  ";
        }
      deallog << std::endl << std::endl;
      deallog.pop();
    }
}


template <int dim, int fe_degree>
void test ()
{
  // 'misuse' fe_degree for setting the parallel
  // option here
  unsigned int parallel_option = 0;
  if (fe_degree == 1)
    parallel_option = 0;
  else if (fe_degree == 2)
    parallel_option = 1;
  else
    return;
  deallog.push("double");
  deallog.threshold_double(1.e-12);
  do_test<dim,double>(parallel_option);
  deallog.pop();
  deallog.push("float");
  deallog.threshold_double(1.e-6);
  do_test<dim,float>(parallel_option);
  deallog.pop();
}

