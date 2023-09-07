// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

// We consider a triangulation that connects two mesh regions via a
// non-conforming interface.
// |----||----|
// |   A||B   |
// |----||----|
// This test checks the access of values on face B in quadrature points of face
// A (and vice versa).

#include <deal.II/base/mpi.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_block_vector.h>

#include <deal.II/matrix_free/fe_remote_evaluation.h>

#include "../tests.h"

using namespace dealii;

using VectorType = LinearAlgebra::distributed::Vector<double>;

constexpr unsigned int dim = 2;

// the integrator used on cells with inner faces at the non-conforming interface
using NonmatchingIntegratorInner = FEFaceEvaluation<dim, -1, 0, dim>;
// communicator for NonmatchingIntegratorOuter
using NonmatchingCommunicator =
  FERemoteEvaluationCommunicator<NonmatchingIntegratorInner, true>;
// integrator which accesses values at corresponding outer faces
using NonmatchingIntegratorOuter =
  FERemoteEvaluation<NonmatchingCommunicator, dim>;

template <int dim>
struct Coordinates : public Function<dim>
{
  Coordinates()
    : Function<dim>(dim, 0.0)
  {}

  double
  value(const Point<dim> &p, const unsigned int component) const final
  {
    return p[component];
  }
};


struct NonmatchingFaceLoop
{
  NonmatchingFaceLoop(
    const MatrixFree<dim>                             &matrix_free,
    std::vector<std::pair<unsigned int, unsigned int>> face_pairs)
    : matrix_free(matrix_free),
      fe_remote_integrator(nonmatching_comm,
    matrix_free
      .get_dof_handler())
  {
    for (const auto &face_pair : face_pairs)
      nonmatching_faces.insert(face_pair.first);

    NonmatchingIntegratorInner fe_eval(matrix_free);
    nonmatching_comm.initialize_face_pairs(face_pairs, fe_eval);
  }


  void
  evaluate(VectorType &dst, const VectorType &src)
  {
    // cache all relevant values on non-matching faces
    fe_remote_integrator.gather_evaluate(src, EvaluationFlags::values);

    matrix_free.loop({}, {}, &NonmatchingFaceLoop::do_evaluate, this, dst, src);
  }

private:
  void
  do_evaluate(const MatrixFree<dim>                       &mf,
              VectorType                                  &dst,
              const VectorType                            &src,
              const std::pair<unsigned int, unsigned int> &range) const
  {
    NonmatchingIntegratorInner integrator_m(mf);

    for (unsigned int face = range.first; face < range.second; ++face)
      {
        if (is_internal_face(face))
          {
            for (unsigned int q = 0; q < integrator_m.n_q_points; ++q)
              {
                integrator_m.reinit(face);
                integrator_m.gather_evaluate(src, EvaluationFlags::values);

                fe_remote_integrator.reinit(face);
                // fe_remote_integrator.gather_evaluate() already called before
                // loop

                integrator_m.get_value(q);
                fe_remote_integrator.get_value(q);
                deallog << "Boundary ID      : "
                        << matrix_free.get_boundary_id(face) << std::endl;
                deallog << "Inner Integrator : " << integrator_m.get_value(q)
                        << std::endl;
                deallog << "Remote Integrator: "
                        << fe_remote_integrator.get_value(q) << std::endl;
                deallog << std::endl;
              }
          }
      }
  }

  bool
  is_internal_face(const unsigned int face) const
  {
    return nonmatching_faces.find(matrix_free.get_boundary_id(face)) !=
           nonmatching_faces.end();
  }

  const MatrixFree<dim>  &matrix_free;
  NonmatchingCommunicator nonmatching_comm;
  // mutable because reinit sets an unsigned integer
  mutable NonmatchingIntegratorOuter fe_remote_integrator;

  std::set<unsigned int> nonmatching_faces;
};


void
test()
{
  const unsigned int fe_degree = 3;

  // create non-matching grid
  Triangulation<dim> tria_0;
  GridGenerator::subdivided_hyper_rectangle(
    tria_0, {7, 7}, {0.0, 0.0}, {1.0, 1.0}, true);

  Triangulation<dim> tria_1;
  GridGenerator::subdivided_hyper_rectangle(
    tria_1, {6, 3}, {1.0, 0.0}, {3.0, 1.0}, true);
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 2 * dim);

  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 0., false, true);

  AssertDimension(tria_0.n_vertices() + tria_1.n_vertices(), tria.n_vertices());

  // store non-matching face pairs
  std::vector<std::pair<unsigned int, unsigned int>> face_pairs;
  face_pairs.emplace_back(1, 2 * dim);
  face_pairs.emplace_back(2 * dim, 1);

  // create DoFHandler and Mapping
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(FESystem<dim>(FE_DGQ<dim>(fe_degree), dim));
  const MappingQ1<dim> mapping;

  // create MatrixFree
  typename MatrixFree<dim>::AdditionalData data;
  data.mapping_update_flags =
    update_quadrature_points | update_gradients | update_values;
  data.mapping_update_flags_boundary_faces = data.mapping_update_flags;
  data.mapping_update_flags_inner_faces    = data.mapping_update_flags;

  MatrixFree<dim> matrix_free;

  matrix_free.reinit(mapping,
                     dof_handler,
                     AffineConstraints<double>(),
                     QGauss<dim>(fe_degree + 1),
                     data);

  VectorType src, dst;
  matrix_free.initialize_dof_vector(src);
  matrix_free.initialize_dof_vector(dst);

  VectorTools::interpolate(mapping, dof_handler, Coordinates<dim>(), src);


  NonmatchingFaceLoop nonmatching_face_loop(matrix_free, face_pairs);
  nonmatching_face_loop.evaluate(dst, src);
}



int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
  MPILogInitAll                    all;

  test();
}
