// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// tests matrix-free partitioners for update_ghost_values and compress(add)


#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "../tests.h"


template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class Test
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  void
  run(unsigned int fe_degree)
  {
    // test FE_DGQ
    run(fe_degree, FE_DGQ<dim>(fe_degree));
    // test FE_Q
    run(fe_degree, FE_DGQHermite<dim>(fe_degree));
  }

  void
  run(unsigned int fe_degree, const FiniteElement<dim> &fe)
  {
    // test different update flags
    run(
      fe_degree,
      fe,
      MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::values);
    run(fe_degree,
        fe,
        MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::
          values_all_faces);
    run(fe_degree,
        fe,
        MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::
          gradients);
    run(fe_degree,
        fe,
        MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::
          gradients_all_faces);
    run(fe_degree,
        fe,
        MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces::
          unspecified);
  }


private:
  void
  run(unsigned int              fe_degree,
      const FiniteElement<dim> &fe,
      typename MatrixFree<dim, Number, VectorizedArrayType>::DataAccessOnFaces
        update_flag)
  {
    this->fe_degree = fe_degree;

    parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);

    {
      std::vector<unsigned int> repetitions(dim, 1);
      repetitions[0] = 2;

      Point<dim> p1, p2;

      for (unsigned int d = 0; d < dim; d++)
        p2[d] = 1.0;

      GridGenerator::subdivided_hyper_rectangle(tria, repetitions, p1, p2);
    }

    DoFHandler<dim> dof_handler(tria);
    dof_handler.distribute_dofs(fe);

    MappingQ<dim> mapping(1);

    QGauss<1> quad(fe_degree + 1);

    AffineConstraints<Number> constraint;
    constraint.close();

    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
      additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData::none;
    additional_data.mapping_update_flags =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.mapping_update_flags_inner_faces =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.mapping_update_flags_boundary_faces =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.mapping_update_flags_faces_by_cells =
      update_gradients | update_JxW_values | update_quadrature_points;
    additional_data.hold_all_faces_to_owned_cells = true;

    MatrixFree<dim, Number, VectorizedArrayType> matrix_free;

    matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

    VectorType src, dst;
    matrix_free.initialize_dof_vector(src);
    matrix_free.initialize_dof_vector(dst);

    for (unsigned int i = 0; i < Utilities::pow(fe_degree + 1, dim); i++)
      src.begin()[i] = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1;

    matrix_free.loop(&Test::dummy_operation_1,
                     &Test::dummy_operation_2,
                     &Test::dummy_operation_2,
                     this,
                     dst,
                     src,
                     false,
                     update_flag,
                     update_flag);

    deallog << "dst:" << std::endl;
    for (unsigned int i = Utilities::pow(fe_degree + 1, dim);
         i < 2 * Utilities::pow(fe_degree + 1, dim);
         i++)
      dst.begin()[i] = 0;
    for (unsigned int i = 0; i < 2 * Utilities::pow(fe_degree + 1, dim); i++)
      deallog << static_cast<int>(dst[i]) << " ";
    deallog << std::endl << std::endl;
  }

  void
  dummy_operation_1(const MatrixFree<dim, Number, VectorizedArrayType> &,
                    VectorType &      dst,
                    const VectorType &src,
                    const std::pair<unsigned int, unsigned int> &) const
  {
    deallog << "src:" << std::endl;
    for (unsigned int i = 0; i < 2 * Utilities::pow(fe_degree + 1, dim); i++)
      deallog << static_cast<int>(src[i]) << " ";
    deallog << std::endl;

    for (unsigned int i = Utilities::pow(fe_degree + 1, dim);
         i < 2 * Utilities::pow(fe_degree + 1, dim);
         i++)
      dst.begin()[i] = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1;
  }

  void
  dummy_operation_2(const MatrixFree<dim, Number, VectorizedArrayType> &,
                    VectorType &,
                    const VectorType &,
                    const std::pair<unsigned int, unsigned int> &) const
  {}

  unsigned int fe_degree;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

  MPILogInitAll all;

  {
    deallog.push("2d");
    Test<2> runner;
    runner.run(3);
    deallog.pop();
  }
}
