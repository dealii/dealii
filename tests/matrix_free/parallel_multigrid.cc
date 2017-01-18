// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2017 by the deal.II authors
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



// test running a multigrid solver on continuous FE_Q finite elements with
// MatrixFree data structures (otherwise similar to step-37)

#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>

std::ofstream logfile("output");


template <int dim, int fe_degree, int n_q_points_1d = fe_degree+1, typename number=double>
class LaplaceOperator : public Subscriptor
{
public:
  LaplaceOperator() {};

  void initialize (const Mapping<dim> &mapping,
                   const DoFHandler<dim> &dof_handler,
                   const std::set<types::boundary_id> &dirichlet_boundaries,
                   const unsigned int level = numbers::invalid_unsigned_int)
  {
    const QGauss<1> quad (n_q_points_1d);
    typename MatrixFree<dim,number>::AdditionalData addit_data;
    addit_data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::none;
    addit_data.level_mg_handler = level;

    // extract the constraints due to Dirichlet boundary conditions
    ConstraintMatrix constraints;
    ZeroFunction<dim> zero;
    typename FunctionMap<dim>::type functions;
    for (std::set<types::boundary_id>::const_iterator it=dirichlet_boundaries.begin();
         it != dirichlet_boundaries.end(); ++it)
      functions[*it] = &zero;
    if (level == numbers::invalid_unsigned_int)
      VectorTools::interpolate_boundary_values(dof_handler, functions, constraints);
    else
      {
        std::vector<types::global_dof_index> local_dofs;
        typename DoFHandler<dim>::cell_iterator
        cell = dof_handler.begin(level),
        endc = dof_handler.end(level);
        for (; cell!=endc; ++cell)
          {
            if (dof_handler.get_triangulation().locally_owned_subdomain()!=numbers::invalid_subdomain_id
                && cell->level_subdomain_id()==numbers::artificial_subdomain_id)
              continue;
            const FiniteElement<dim> &fe = cell->get_fe();
            local_dofs.resize(fe.dofs_per_face);

            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
                 ++face_no)
              if (cell->at_boundary(face_no) == true)
                {
                  const typename DoFHandler<dim>::face_iterator
                  face = cell->face(face_no);
                  const types::boundary_id bi = face->boundary_id();
                  if (functions.find(bi) != functions.end())
                    {
                      face->get_mg_dof_indices(level, local_dofs);
                      for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                        constraints.add_line(local_dofs[i]);
                    }
                }
          }
      }
    constraints.close();

    data.reinit (mapping, dof_handler, constraints, quad, addit_data);

    compute_inverse_diagonal();
  }

  void vmult(LinearAlgebra::distributed::Vector<number> &dst,
             const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void Tvmult(LinearAlgebra::distributed::Vector<number> &dst,
              const LinearAlgebra::distributed::Vector<number> &src) const
  {
    dst = 0;
    vmult_add(dst, src);
  }

  void Tvmult_add(LinearAlgebra::distributed::Vector<number> &dst,
                  const LinearAlgebra::distributed::Vector<number> &src) const
  {
    vmult_add(dst, src);
  }

  void vmult_add(LinearAlgebra::distributed::Vector<number> &dst,
                 const LinearAlgebra::distributed::Vector<number> &src) const
  {
    data.cell_loop (&LaplaceOperator::local_apply,
                    this, dst, src);

    const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
    for (unsigned int i=0; i<constrained_dofs.size(); ++i)
      dst.local_element(constrained_dofs[i]) += src.local_element(constrained_dofs[i]);
  }

  types::global_dof_index m() const
  {
    return data.get_vector_partitioner()->size();
  }

  types::global_dof_index n() const
  {
    return data.get_vector_partitioner()->size();
  }

  number el (const unsigned int row,  const unsigned int col) const
  {
    AssertThrow(false, ExcMessage("Matrix-free does not allow for entry access"));
    return number();
  }

  void
  initialize_dof_vector(LinearAlgebra::distributed::Vector<number> &vector) const
  {
    if (!vector.partitioners_are_compatible(*data.get_dof_info(0).vector_partitioner))
      data.initialize_dof_vector(vector);
    Assert(vector.partitioners_are_globally_compatible(*data.get_dof_info(0).vector_partitioner),
           ExcInternalError());
  }

  const LinearAlgebra::distributed::Vector<number> &
  get_matrix_diagonal_inverse() const
  {
    Assert(inverse_diagonal_entries.size() > 0, ExcNotInitialized());
    return inverse_diagonal_entries;
  }


private:
  void
  local_apply (const MatrixFree<dim,number>                &data,
               LinearAlgebra::distributed::Vector<number>       &dst,
               const LinearAlgebra::distributed::Vector<number> &src,
               const std::pair<unsigned int,unsigned int>  &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,1,number> phi (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        phi.evaluate (false,true,false);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_gradient (phi.get_gradient(q), q);
        phi.integrate (false,true);
        phi.distribute_local_to_global (dst);
      }
  }

  void
  compute_inverse_diagonal ()
  {
    data.initialize_dof_vector(inverse_diagonal_entries);
    unsigned int dummy;
    data.cell_loop (&LaplaceOperator::local_diagonal_cell,
                    this, inverse_diagonal_entries, dummy);

    for (unsigned int i=0; i<inverse_diagonal_entries.local_size(); ++i)
      if (std::abs(inverse_diagonal_entries.local_element(i)) > 1e-10)
        inverse_diagonal_entries.local_element(i) = 1./inverse_diagonal_entries.local_element(i);
      else
        inverse_diagonal_entries.local_element(i) = 1.;
  }

  void
  local_diagonal_cell (const MatrixFree<dim,number>                &data,
                       LinearAlgebra::distributed::Vector<number>       &dst,
                       const unsigned int &,
                       const std::pair<unsigned int,unsigned int>  &cell_range) const
  {
    FEEvaluation<dim,fe_degree,n_q_points_1d,1,number> phi (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);

        VectorizedArray<number> local_diagonal_vector[phi.tensor_dofs_per_cell];
        for (unsigned int i=0; i<phi.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
              phi.begin_dof_values()[j] = VectorizedArray<number>();
            phi.begin_dof_values()[i] = 1.;
            phi.evaluate (false,true,false);
            for (unsigned int q=0; q<phi.n_q_points; ++q)
              phi.submit_gradient (phi.get_gradient(q), q);
            phi.integrate (false,true);
            local_diagonal_vector[i] = phi.begin_dof_values()[i];
          }
        for (unsigned int i=0; i<phi.tensor_dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = local_diagonal_vector[i];
        phi.distribute_local_to_global (dst);
      }
  }

  MatrixFree<dim,number> data;
  LinearAlgebra::distributed::Vector<number> inverse_diagonal_entries;
};



template <typename MatrixType>
class MGTransferPrebuiltMF : public MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double> >
{
public:
  MGTransferPrebuiltMF(const MGLevelObject<MatrixType> &laplace)
    :
    laplace_operator (laplace)
  {};

  /**
   * Overload copy_to_mg from MGTransferPrebuilt to get the vectors compatible
   * with MatrixFree and bypass the crude initialization in MGTransferPrebuilt
   */
  template <int dim, class InVector, int spacedim>
  void
  copy_to_mg (const DoFHandler<dim,spacedim> &mg_dof,
              MGLevelObject<LinearAlgebra::distributed::Vector<double> > &dst,
              const InVector &src) const
  {
    for (unsigned int level=dst.min_level();
         level<=dst.max_level(); ++level)
      laplace_operator[level].initialize_dof_vector(dst[level]);
    MGTransferPrebuilt<LinearAlgebra::distributed::Vector<double> >::copy_to_mg(mg_dof, dst, src);
  }

private:
  const MGLevelObject<MatrixType> &laplace_operator;
};



template<typename MatrixType, typename Number>
class MGCoarseIterative : public MGCoarseGridBase<LinearAlgebra::distributed::Vector<Number> >
{
public:
  MGCoarseIterative() {}

  void initialize(const MatrixType &matrix)
  {
    coarse_matrix = &matrix;
  }

  virtual void operator() (const unsigned int   level,
                           LinearAlgebra::distributed::Vector<double> &dst,
                           const LinearAlgebra::distributed::Vector<double> &src) const
  {
    ReductionControl solver_control (1e4, 1e-50, 1e-10);
    SolverCG<LinearAlgebra::distributed::Vector<double> > solver_coarse (solver_control);
    solver_coarse.solve (*coarse_matrix, dst, src, PreconditionIdentity());
  }

  const MatrixType *coarse_matrix;
};




template <int dim, int fe_degree, int n_q_points_1d, typename number>
void do_test (const DoFHandler<dim>  &dof)
{
  if (types_are_equal<number,float>::value == true)
    {
      deallog.push("float");
      deallog.threshold_double(1e-6);
    }
  else
    {
      deallog.threshold_double(5.e-11);
    }

  deallog << "Testing " << dof.get_fe().get_name();
  deallog << std::endl;
  deallog << "Number of degrees of freedom: " << dof.n_dofs() << std::endl;

  MappingQ<dim> mapping(fe_degree+1);
  LaplaceOperator<dim,fe_degree,n_q_points_1d,number> fine_matrix;
  std::set<types::boundary_id> dirichlet_boundaries;
  dirichlet_boundaries.insert(0);
  fine_matrix.initialize(mapping, dof, dirichlet_boundaries);

  LinearAlgebra::distributed::Vector<number> in, sol;
  fine_matrix.initialize_dof_vector(in);
  fine_matrix.initialize_dof_vector(sol);

  // set constant rhs vector
  in = 1.;

  // set up multigrid in analogy to step-37
  typedef LaplaceOperator<dim,fe_degree,n_q_points_1d,number> LevelMatrixType;

  MGLevelObject<LevelMatrixType> mg_matrices;
  mg_matrices.resize(0, dof.get_triangulation().n_global_levels()-1);
  for (unsigned int level = 0; level<dof.get_triangulation().n_global_levels(); ++level)
    {
      mg_matrices[level].initialize(mapping, dof, dirichlet_boundaries, level);
    }

  MGTransferPrebuiltMF<LevelMatrixType> mg_transfer(mg_matrices);
  mg_transfer.build_matrices(dof);

  MGCoarseIterative<LevelMatrixType,number> mg_coarse;
  mg_coarse.initialize(mg_matrices[0]);

  typedef PreconditionChebyshev<LevelMatrixType,LinearAlgebra::distributed::Vector<number> > SMOOTHER;
  MGSmootherPrecondition<LevelMatrixType, SMOOTHER, LinearAlgebra::distributed::Vector<number> >
  mg_smoother;

  MGLevelObject<typename SMOOTHER::AdditionalData> smoother_data;
  smoother_data.resize(0, dof.get_triangulation().n_global_levels()-1);
  for (unsigned int level = 0; level<dof.get_triangulation().n_global_levels(); ++level)
    {
      smoother_data[level].smoothing_range = 15.;
      smoother_data[level].degree = 5;
      smoother_data[level].eig_cg_n_iterations = 15;
      smoother_data[level].preconditioner.
      reset(new DiagonalMatrix<LinearAlgebra::distributed::Vector<number> >());
      smoother_data[level].preconditioner->get_vector() =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }
  mg_smoother.initialize(mg_matrices, smoother_data);

  mg::Matrix<LinearAlgebra::distributed::Vector<double> >
  mg_matrix(mg_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<double> > mg(dof,
                                                            mg_matrix,
                                                            mg_coarse,
                                                            mg_transfer,
                                                            mg_smoother,
                                                            mg_smoother);
  PreconditionMG<dim, LinearAlgebra::distributed::Vector<double>,
                 MGTransferPrebuiltMF<LevelMatrixType> >
                 preconditioner(dof, mg, mg_transfer);

  {
    ReductionControl control(30, 1e-20, 1e-7);
    SolverCG<LinearAlgebra::distributed::Vector<double> > solver(control);
    solver.solve(fine_matrix, sol, in, preconditioner);
  }

  if (types_are_equal<number,float>::value == true)
    deallog.pop();
}



template <int dim, int fe_degree>
void test ()
{
  for (unsigned int i=5; i<8; ++i)
    {
      parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD,
                                                     Triangulation<dim>::limit_level_difference_at_vertices,
                                                     parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);
      GridGenerator::hyper_cube (tria);
      tria.refine_global(i-dim);

      FE_Q<dim> fe (fe_degree);
      DoFHandler<dim> dof (tria);
      dof.distribute_dofs(fe);
      dof.distribute_mg_dofs(fe);

      do_test<dim, fe_degree, fe_degree+1, double> (dof);
    }
}



int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      deallog.attach(logfile);
      deallog << std::setprecision (4);
    }

  {
    deallog.threshold_double(1.e-10);
    deallog.push("2d");
    test<2,1>();
    test<2,2>();
    deallog.pop();
    deallog.push("3d");
    test<3,1>();
    test<3,2>();
    deallog.pop();
  }
}
