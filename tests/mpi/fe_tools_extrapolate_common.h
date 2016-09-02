// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/std_cxx11/unique_ptr.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <string>


template <int dim>
class  TestFunction : public Function<dim>
{
public:
  TestFunction (unsigned int deg)
    : Function<dim>(), degree(deg)
  {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component = 0) const
  {
    double return_value = 0.;
    for (unsigned int d=0; d<dim; ++d)
      return_value += std::pow(std::abs(.5-p(d)), degree);

    return return_value;
  }
private:
  const unsigned int degree;
};



template <int dim>
parallel::distributed::Triangulation<dim> *make_tria ()
{
  parallel::distributed::Triangulation<dim> *tria = new parallel::distributed::Triangulation<dim>(MPI_COMM_WORLD);
  typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
  GridGenerator::hyper_cube(*tria, 0., 1.);
  tria->refine_global (2);
  for (int i=0; i<2; ++i)
    {
      cell = tria->begin_active();
      cell->set_refine_flag();
      ++cell;
      cell->set_refine_flag();

      tria->execute_coarsening_and_refinement ();
    }
  return tria;
}



template <int dim>
DoFHandler<dim> *make_dof_handler (const parallel::distributed::Triangulation<dim> &tria,
                                   const FiniteElement<dim> &fe)
{
  DoFHandler<dim> *dof_handler = new DoFHandler<dim>(tria);
  dof_handler->distribute_dofs (fe);
  return dof_handler;
}



// output some indicators for a given vector
template <unsigned int dim, typename VectorType>
void
output_vector (const VectorType &v, const std::string &output_name, const DoFHandler<dim> &dof_handler)
{
  //v.print(deallog.get_file_stream());

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (v, output_name, DataOut<dim>::type_dof_data);

  const unsigned int degree = dof_handler.get_fe().degree;
  data_out.build_patches (degree);

  const std::string filename = (output_name +
                                "." +
                                Utilities::int_to_string
                                (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 1));
  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);

}



template <int dim, typename VectorType>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  deallog << std::setprecision (10);

  // only check if both elements have support points. 
  // otherwise, interpolation doesn't really work
  if ((fe1.get_unit_support_points().size() == 0) ||
      (fe2.get_unit_support_points().size() == 0))
    return;
  //  likewise for non-primitive elements
  if (!fe1.is_primitive() || !fe2.is_primitive())
    return;
  // we need to have dof_constraints for this test
  if (!fe2.constraints_are_implemented())
    return;
  // we need prolongation matrices in fe2
  if (!fe2.isotropic_restriction_is_implemented())
    return;

  std_cxx11::unique_ptr<parallel::distributed::Triangulation<dim> > tria(make_tria<dim>());

  std_cxx11::unique_ptr<DoFHandler<dim> >    dof1(make_dof_handler (*tria, fe1));
  std_cxx11::unique_ptr<DoFHandler<dim> >    dof2(make_dof_handler (*tria, fe2));
  ConstraintMatrix cm1;
  DoFTools::make_hanging_node_constraints (*dof1, cm1);
  cm1.close ();
  ConstraintMatrix cm2;
  DoFTools::make_hanging_node_constraints (*dof2, cm2);
  cm2.close ();

  IndexSet locally_owned_dofs1 = dof1->locally_owned_dofs();
  IndexSet locally_relevant_dofs1;
  DoFTools::extract_locally_relevant_dofs (*dof1, locally_relevant_dofs1);
  IndexSet locally_owned_dofs2 = dof2->locally_owned_dofs();
  IndexSet locally_relevant_dofs2;
  DoFTools::extract_locally_relevant_dofs (*dof2, locally_relevant_dofs2);

  VectorType in_ghosted (locally_owned_dofs1, locally_relevant_dofs1, MPI_COMM_WORLD);
  VectorType in_distributed (locally_owned_dofs1, MPI_COMM_WORLD);
  VectorType out_distributed (locally_owned_dofs2, MPI_COMM_WORLD);
  VectorType out_ghosted (locally_owned_dofs2, locally_relevant_dofs2, MPI_COMM_WORLD);
  VectorType out_reference (locally_owned_dofs2, MPI_COMM_WORLD);

  // Choose some reference function of sufficiently high polynomial degree.
  TestFunction<dim> function (std::max(fe1.degree,fe2.degree));
  VectorTools::interpolate (*dof1, function, in_distributed);
  cm1.distribute(in_distributed);
  in_ghosted = in_distributed;

  Vector<double> difference_before;
  VectorTools::integrate_difference (*dof1, in_ghosted, function, difference_before,
                                     QGauss<dim>(fe1.degree+2), VectorTools::L2_norm);
  const double local_error_before = difference_before.l2_norm();
  const double global_error_before
    = std::sqrt(Utilities::MPI::sum(std::pow(local_error_before,2.), MPI_COMM_WORLD));

  /*output_vector<dim, VectorType> (in_ghosted,
                                  Utilities::int_to_string(fe1.degree,1)+Utilities::int_to_string(dim,1)+std::string("in"),
                                  *dof1);*/
  {
    const double l1_norm = in_distributed.l1_norm();
    const double l2_norm = in_distributed.l2_norm();
    const double linfty_norm = in_distributed.linfty_norm();
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
      deallog << l1_norm << ' '
              << l2_norm << ' '
              << linfty_norm << std::endl;
  }

  FETools::extrapolate_parallel (*dof1, in_ghosted, *dof2, cm2, out_distributed);
  out_distributed.compress(VectorOperation::insert);
  out_ghosted = out_distributed;

  /*output_vector<dim, VectorType> (out_ghosted,
                                  Utilities::int_to_string(fe2.degree,1)+Utilities::int_to_string(dim,1)+std::string("out"),
                                  *dof2);*/

  {
    const double l1_norm = out_distributed.l1_norm();
    const double l2_norm = out_distributed.l2_norm();
    const double linfty_norm = out_distributed.linfty_norm();
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
      deallog << l1_norm << ' '
              << l2_norm << ' '
              << linfty_norm << std::endl;
  }

  Vector<double> difference_after;
  VectorTools::integrate_difference (*dof2, out_ghosted, function, difference_after,
                                     QGauss<dim>(fe2.degree+2), VectorTools::L2_norm);
  const double local_error_after = difference_after.l2_norm();
  const double global_error_after
    = std::sqrt(Utilities::MPI::sum(std::pow(local_error_after,2.), MPI_COMM_WORLD));

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    {
      deallog << "global_error_before: "
              << global_error_before << std::endl;
      deallog << "global_error_after: "
              << global_error_after << std::endl;
    }
  if (fe2.degree == fe1.degree)
    {
      out_distributed -= in_distributed;
      AssertThrow(out_distributed.l2_norm()<1.e-8,
                  ExcInternalError());

    }
  if (fe2.degree > fe1.degree)
    {
      AssertThrow(global_error_after < global_error_before,
                  ExcInternalError());
    }
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    deallog << "OK" << std::endl;
}


template <int dim, typename VectorType>
void
check_this_dealii (const FiniteElement<dim> &fe1,
                   const FiniteElement<dim> &fe2)
{
  deallog << std::setprecision (10);

  // only check if both elements have support points. 
  // otherwise, interpolation doesn't really work
  if ((fe1.get_unit_support_points().size() == 0) ||
      (fe2.get_unit_support_points().size() == 0))
    return;
  //  likewise for non-primitive elements
  if (!fe1.is_primitive() || !fe2.is_primitive())
    return;
  // we need to have dof_constraints for this test
  if (!fe2.constraints_are_implemented())
    return;
  // we need prolongation matrices in fe2
  if (!fe2.isotropic_restriction_is_implemented())
    return;

  std_cxx11::unique_ptr<parallel::distributed::Triangulation<dim> > tria(make_tria<dim>());

  std_cxx11::unique_ptr<DoFHandler<dim> >    dof1(make_dof_handler (*tria, fe1));
  std_cxx11::unique_ptr<DoFHandler<dim> >    dof2(make_dof_handler (*tria, fe2));
  ConstraintMatrix cm1;
  DoFTools::make_hanging_node_constraints (*dof1, cm1);
  cm1.close ();
  ConstraintMatrix cm2;
  DoFTools::make_hanging_node_constraints (*dof2, cm2);
  cm2.close ();

  IndexSet locally_owned_dofs1 = dof1->locally_owned_dofs();
  IndexSet locally_relevant_dofs1;
  DoFTools::extract_locally_relevant_dofs (*dof1, locally_relevant_dofs1);
  IndexSet locally_owned_dofs2 = dof2->locally_owned_dofs();
  IndexSet locally_relevant_dofs2;
  DoFTools::extract_locally_relevant_dofs (*dof2, locally_relevant_dofs2);

  VectorType in_ghosted (locally_owned_dofs1, locally_relevant_dofs1, MPI_COMM_WORLD);
  VectorType out_ghosted (locally_owned_dofs2, locally_relevant_dofs2, MPI_COMM_WORLD);
  VectorType out_reference (locally_owned_dofs2, locally_relevant_dofs2, MPI_COMM_WORLD);

  // Choose some reference function of sufficiently high polynomial degree.
  TestFunction<dim> function (std::max(fe1.degree,fe2.degree));
  VectorTools::interpolate (*dof1, function, in_ghosted);
  cm1.distribute(in_ghosted);
  in_ghosted.update_ghost_values();
  VectorTools::interpolate (*dof2, function, out_reference);
  cm2.distribute(out_reference);
  out_reference.update_ghost_values();

  Vector<double> difference_before;
  VectorTools::integrate_difference (*dof1, in_ghosted, function, difference_before,
                                     QGauss<dim>(fe1.degree+2), VectorTools::L2_norm);
  const double local_error_before = difference_before.l2_norm();
  const double global_error_before
    = std::sqrt(Utilities::MPI::sum(std::pow(local_error_before,2.), MPI_COMM_WORLD));

  /*output_vector<dim, VectorType> (in_ghosted,
                                  Utilities::int_to_string(fe1.degree,1)+Utilities::int_to_string(dim,1)+std::string("in"),
                                  *dof1);*/
  {
    const double l1_norm = in_ghosted.l1_norm();
    const double l2_norm = in_ghosted.l2_norm();
    const double linfty_norm = in_ghosted.linfty_norm();
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
      deallog << l1_norm << ' '
              << l2_norm << ' '
              << linfty_norm << std::endl;
  }

  if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)==1)
    FETools::extrapolate (*dof1, in_ghosted, *dof2, cm2, out_ghosted);
  else
    FETools::extrapolate_parallel (*dof1, in_ghosted, *dof2, cm2, out_ghosted);

  out_ghosted.update_ghost_values();

  /*output_vector<dim, VectorType> (out_ghosted,
                                  Utilities::int_to_string(fe2.degree,1)+Utilities::int_to_string(dim,1)+std::string("out"),
                                  *dof2);*/
  {
    const double l1_norm = out_ghosted.l1_norm();
    const double l2_norm = out_ghosted.l2_norm();
    const double linfty_norm = out_ghosted.linfty_norm();
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
      deallog << l1_norm << ' '
              << l2_norm << ' '
              << linfty_norm << std::endl;
  }

  Vector<double> difference_after;
  VectorTools::integrate_difference (*dof2, out_ghosted, function, difference_after,
                                     QGauss<dim>(fe2.degree+2), VectorTools::L2_norm);
  const double local_error_after = difference_after.l2_norm();
  const double global_error_after
    = std::sqrt(Utilities::MPI::sum(std::pow(local_error_after,2.), MPI_COMM_WORLD));

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    {
      deallog << "global_error_before: "
              << global_error_before << std::endl;
      deallog << "global_error_after: "
              << global_error_after << std::endl;
    }
  if (fe2.degree == fe1.degree)
    {
      out_ghosted -= in_ghosted;
      AssertThrow(out_ghosted.l2_norm()<1.e-8,
                  ExcInternalError());

    }
  if (fe2.degree > fe1.degree)
    {
      AssertThrow(global_error_after < global_error_before,
                  ExcInternalError());
    }
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    deallog << "OK" << std::endl;
}
