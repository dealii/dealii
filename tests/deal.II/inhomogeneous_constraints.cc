// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// this function tests the correctness of the implementation of
// inhomogeneous constraints. The program is a modification of the step-27
// tutorial program with hp elements and the constraints arising in that
// situation. the idea of the test is to set up a matrix with standard tools
// (i.e., constraints and the boundary value list), and compare that with
// the new function.

#include "../tests.h"

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <fstream>
#include <iostream>
#include <complex>

std::ofstream logfile("output");

using namespace dealii;

template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem ();
  ~LaplaceProblem ();

  void run ();

private:
  void setup_system ();
  void test_equality ();
  void assemble_reference ();
  void assemble_test_1 ();
  void assemble_test_2 ();
  void solve ();
  void create_coarse_grid ();
  void estimate_smoothness (Vector<float> &smoothness_indicators) const;
  void postprocess ();

  Triangulation<dim>   triangulation;

  hp::DoFHandler<dim>      dof_handler;
  hp::FECollection<dim>    fe_collection;
  hp::QCollection<dim>     quadrature_collection;
  hp::QCollection<dim-1>   face_quadrature_collection;

  ConstraintMatrix     hanging_nodes_only;
  ConstraintMatrix     test_all_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> reference_matrix;
  SparseMatrix<double> test_matrix;

  Vector<double>       solution;
  Vector<double>       reference_rhs;
  Vector<double>       test_rhs;

  const unsigned int max_degree;
};



template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const;
};


template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>   &p,
                            const unsigned int  /*component*/) const
{
  double sum = 0;
  for (unsigned int d=0; d<dim; ++d)
    sum += std::sin(deal_II_numbers::PI*p[d]);
  return sum;
}


template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim> () {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
                           const unsigned int  /*component*/) const
{
  double product = 1;
  for (unsigned int d=0; d<dim; ++d)
    product *= (p[d]+1);
  return product;
}


template <int dim>
LaplaceProblem<dim>::LaplaceProblem ()
  :
  dof_handler (triangulation),
  max_degree (5)
{
  if (dim == 2)
    for (unsigned int degree=2; degree<=max_degree; ++degree)
      {
        fe_collection.push_back (FE_Q<dim>(degree));
        quadrature_collection.push_back (QGauss<dim>(degree+1));
        face_quadrature_collection.push_back (QGauss<dim-1>(degree+1));
      }
  else
    for (unsigned int degree=1; degree<max_degree-1; ++degree)
      {
        fe_collection.push_back (FE_Q<dim>(degree));
        quadrature_collection.push_back (QGauss<dim>(degree+1));
        face_quadrature_collection.push_back (QGauss<dim-1>(degree+1));
      }
}


template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
  dof_handler.clear ();
}


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe_collection);

  solution.reinit (dof_handler.n_dofs());
  reference_rhs.reinit (dof_handler.n_dofs());
  test_rhs.reinit (dof_handler.n_dofs());

  hanging_nodes_only.clear ();
  test_all_constraints.clear ();

  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_nodes_only);
  hanging_nodes_only.close ();

  // add boundary conditions as inhomogeneous
  // constraints here, do it after having
  // added the hanging node constraints in
  // order to be consistent and skip dofs
  // that are already constrained (i.e., are
  // hanging nodes on the boundary in 3D). In
  // contrast to step-27, we choose a sine
  // function.
  {
    test_all_constraints.merge (hanging_nodes_only);
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              BoundaryValues<dim>(),
                                              boundary_values);
    std::map<types::global_dof_index,double>::const_iterator boundary_value = boundary_values.begin();
    for ( ; boundary_value !=boundary_values.end(); ++boundary_value)
      {
        if (!test_all_constraints.is_constrained(boundary_value->first))
          {
            test_all_constraints.add_line(boundary_value->first);
            test_all_constraints.set_inhomogeneity (boundary_value->first,
                                                    boundary_value->second);
          }
      }
  }
  test_all_constraints.close ();

  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
                                       dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp,
                                   hanging_nodes_only, true);
  sparsity_pattern.copy_from (csp);

  reference_matrix.reinit (sparsity_pattern);
  test_matrix.reinit (sparsity_pattern);
}



// test whether we are equal with the
// standard matrix and right hand side
template <int dim>
void LaplaceProblem<dim>::test_equality ()
{
  // need to manually go through the
  // matrix, since we can have different
  // entries in constrained lines.
  for (unsigned int i=0; i<reference_matrix.m(); ++i)
    {
      SparseMatrix<double>::const_iterator reference = reference_matrix.begin(i);
      SparseMatrix<double>::iterator test = test_matrix.begin(i);
      if (test_all_constraints.is_constrained(i) == false)
        {
          for ( ; test != test_matrix.end(i); ++test, ++reference)
            test->value() -= reference->value();
        }
      else
        for ( ; test != test_matrix.end(i); ++test)
          test->value() = 0;
    }

  deallog << "  Matrix difference norm: "
          << test_matrix.frobenius_norm() << std::endl;
  Assert (test_matrix.frobenius_norm() < 1e-13, ExcInternalError());

  // same here -- Dirichlet lines will have
  // nonzero rhs, whereas we will have zero
  // rhs when using inhomogeneous
  // constraints.
  for (unsigned int i=0; i<reference_matrix.m(); ++i)
    if (test_all_constraints.is_constrained(i) == false)
      test_rhs(i) -= reference_rhs(i);
    else
      test_rhs(i) = 0;

  deallog << "  RHS difference norm: "
          << test_rhs.l2_norm() << std::endl;

  Assert (test_rhs.l2_norm() < 1e-14, ExcInternalError());
}




template <int dim>
void LaplaceProblem<dim>::assemble_reference ()
{
  reference_matrix = 0;
  reference_rhs = 0;

  hp::FEValues<dim> hp_fe_values (fe_collection,
                                  quadrature_collection,
                                  update_values    |  update_gradients |
                                  update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;

  FullMatrix<double>   cell_matrix;
  Vector<double>       cell_rhs;

  std::vector<types::global_dof_index> local_dof_indices;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      cell_matrix = 0;

      cell_rhs.reinit (dofs_per_cell);
      cell_rhs = 0;

      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

      std::vector<double>  rhs_values (fe_values.n_quadrature_points);
      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      for (unsigned int q_point=0;
           q_point<fe_values.n_quadrature_points;
           ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      reference_matrix.add(local_dof_indices, cell_matrix);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        reference_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  hanging_nodes_only.condense (reference_matrix, reference_rhs);
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            BoundaryValues<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      reference_matrix,
                                      solution,
                                      reference_rhs);

  deallog << "  Reference matrix nonzeros: " << reference_matrix.n_nonzero_elements()
          << ", actually: " << reference_matrix.n_actually_nonzero_elements (1e-10)
          << std::endl;
}



template <int dim>
void LaplaceProblem<dim>::assemble_test_1 ()
{
  test_matrix = 0;
  test_rhs = 0;

  hp::FEValues<dim> hp_fe_values (fe_collection,
                                  quadrature_collection,
                                  update_values    |  update_gradients |
                                  update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;

  FullMatrix<double>   cell_matrix;
  Vector<double>       cell_rhs;

  std::vector<types::global_dof_index> local_dof_indices;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      cell_matrix = 0;

      cell_rhs.reinit (dofs_per_cell);
      cell_rhs = 0;

      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

      std::vector<double>  rhs_values (fe_values.n_quadrature_points);
      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      for (unsigned int q_point=0;
           q_point<fe_values.n_quadrature_points;
           ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      test_matrix.add(local_dof_indices, cell_matrix);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        test_rhs(local_dof_indices[i]) += cell_rhs(i);

    }

  test_all_constraints.condense (test_matrix, test_rhs);
  deallog << "  Test matrix 1 nonzeros: " << test_matrix.n_nonzero_elements()
          << ", actually: " << test_matrix.n_actually_nonzero_elements (1e-10)
          << std::endl;

  test_equality();
}



template <int dim>
void LaplaceProblem<dim>::assemble_test_2 ()
{
  test_matrix = 0;
  test_rhs = 0;

  hp::FEValues<dim> hp_fe_values (fe_collection,
                                  quadrature_collection,
                                  update_values    |  update_gradients |
                                  update_quadrature_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;

  FullMatrix<double>   cell_matrix;
  Vector<double>       cell_rhs;

  std::vector<types::global_dof_index> local_dof_indices;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      cell_matrix = 0;

      cell_rhs.reinit (dofs_per_cell);
      cell_rhs = 0;

      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

      std::vector<double>  rhs_values (fe_values.n_quadrature_points);
      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      for (unsigned int q_point=0;
           q_point<fe_values.n_quadrature_points;
           ++q_point)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                   fe_values.shape_grad(j,q_point) *
                                   fe_values.JxW(q_point));

            cell_rhs(i) += (fe_values.shape_value(i,q_point) *
                            rhs_values[q_point] *
                            fe_values.JxW(q_point));
          }

      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      test_all_constraints.distribute_local_to_global (cell_matrix,
                                                       cell_rhs,
                                                       local_dof_indices,
                                                       test_matrix,
                                                       test_rhs);
    }
  deallog << "  Test matrix 2 nonzeros: " << test_matrix.n_nonzero_elements()
          << ", actually: " << test_matrix.n_actually_nonzero_elements (1e-10)
          << std::endl;
  test_equality();
}



template <int dim>
void LaplaceProblem<dim>::solve ()
{
  SolverControl           solver_control (reference_rhs.size(),
                                          1e-8*reference_rhs.l2_norm());
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(reference_matrix, 1.2);

  cg.solve (reference_matrix, solution, reference_rhs,
            preconditioner);

  Vector<double> solution_test (solution);
  hanging_nodes_only.distribute (solution);

  // test also distribute function
  test_all_constraints.distribute(solution_test);
  solution_test -= solution;
  deallog << "Distribute error: " << solution_test.l2_norm () << std::endl;
  Assert (solution_test.l2_norm() < 1e-8, ExcInternalError());
}


template <int dim>
void LaplaceProblem<dim>::postprocess ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  for (unsigned int i=0; i<estimated_error_per_cell.size(); ++i)
    estimated_error_per_cell(i) = i;

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);
  triangulation.execute_coarsening_and_refinement ();
}


template <>
void LaplaceProblem<2>::create_coarse_grid ()
{
  const unsigned int dim = 2;

  static const Point<2> vertices_1[]
    = {  Point<2> (-1.,   -1.),
         Point<2> (-1./2, -1.),
         Point<2> (0.,    -1.),
         Point<2> (+1./2, -1.),
         Point<2> (+1,    -1.),

         Point<2> (-1.,   -1./2.),
         Point<2> (-1./2, -1./2.),
         Point<2> (0.,    -1./2.),
         Point<2> (+1./2, -1./2.),
         Point<2> (+1,    -1./2.),

         Point<2> (-1.,   0.),
         Point<2> (-1./2, 0.),
         Point<2> (+1./2, 0.),
         Point<2> (+1,    0.),

         Point<2> (-1.,   1./2.),
         Point<2> (-1./2, 1./2.),
         Point<2> (0.,    1./2.),
         Point<2> (+1./2, 1./2.),
         Point<2> (+1,    1./2.),

         Point<2> (-1.,   1.),
         Point<2> (-1./2, 1.),
         Point<2> (0.,    1.),
         Point<2> (+1./2, 1.),
         Point<2> (+1,    1.)
      };
  const unsigned int
  n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
  const std::vector<Point<dim> > vertices (&vertices_1[0],
                                           &vertices_1[n_vertices]);
  static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
  = {{0, 1, 5, 6},
    {1, 2, 6, 7},
    {2, 3, 7, 8},
    {3, 4, 8, 9},
    {5, 6, 10, 11},
    {8, 9, 12, 13},
    {10, 11, 14, 15},
    {12, 13, 17, 18},
    {14, 15, 19, 20},
    {15, 16, 20, 21},
    {16, 17, 21, 22},
    {17, 18, 22, 23}
  };
  const unsigned int
  n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
  for (unsigned int i=0; i<n_cells; ++i)
    {
      for (unsigned int j=0;
           j<GeometryInfo<dim>::vertices_per_cell;
           ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  triangulation.create_triangulation (vertices,
                                      cells,
                                      SubCellData());
  triangulation.refine_global (2);
}


template <>
void LaplaceProblem<3>::create_coarse_grid ()
{
  GridGenerator::hyper_cube (triangulation);
  triangulation.refine_global (1);
}



template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<3; ++cycle)
    {
      if (cycle == 0)
        create_coarse_grid ();

      setup_system ();

      deallog << std::endl << std::endl
              << "   Number of active cells:             "
              << triangulation.n_active_cells()
              << std::endl
              << "   Number of degrees of freedom:       "
              << dof_handler.n_dofs()
              << std::endl
              << "   Number of hanging node constraints: "
              << hanging_nodes_only.n_constraints()
              << std::endl
              << "   Total number of constraints:        "
              << test_all_constraints.n_constraints()
              << std::endl;

      assemble_reference ();
      assemble_test_1 ();
      assemble_test_2 ();

      solve ();
      if (cycle < 2)
        postprocess ();
    }
}


// this function is copied verbatim from
// step-27
template <int dim>
void
LaplaceProblem<dim>::
estimate_smoothness (Vector<float> &smoothness_indicators) const
{
  const unsigned int N = max_degree;

  std::vector<Tensor<1,dim> > k_vectors;
  std::vector<unsigned int>   k_vectors_magnitude;
  switch (dim)
    {
    case 2:
    {
      for (unsigned int i=0; i<N; ++i)
        for (unsigned int j=0; j<N; ++j)
          if (!((i==0) && (j==0))
              &&
              (i*i + j*j < N*N))
            {
              k_vectors.push_back (Point<dim>(numbers::PI * i,
                                              numbers::PI * j));
              k_vectors_magnitude.push_back (i*i+j*j);
            }

      break;
    }

    case 3:
    {
      for (unsigned int i=0; i<N; ++i)
        for (unsigned int j=0; j<N; ++j)
          for (unsigned int k=0; k<N; ++k)
            if (!((i==0) && (j==0) && (k==0))
                &&
                (i*i + j*j + k*k < N*N))
              {
                k_vectors.push_back (Point<dim>(numbers::PI * i,
                                                numbers::PI * j,
                                                numbers::PI * k));
                k_vectors_magnitude.push_back (i*i+j*j+k*k);
              }

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }

  const unsigned n_fourier_modes = k_vectors.size();
  std::vector<double> ln_k (n_fourier_modes);
  for (unsigned int i=0; i<n_fourier_modes; ++i)
    ln_k[i] = std::log (k_vectors[i].norm());

  std::vector<Table<2,std::complex<double> > >
  fourier_transform_matrices (fe_collection.size());
  QGauss<1>      base_quadrature (2);
  QIterated<dim> quadrature (base_quadrature, N);

  for (unsigned int fe=0; fe<fe_collection.size(); ++fe)
    {
      fourier_transform_matrices[fe].reinit (n_fourier_modes,
                                             fe_collection[fe].dofs_per_cell);

      for (unsigned int k=0; k<n_fourier_modes; ++k)
        for (unsigned int j=0; j<fe_collection[fe].dofs_per_cell; ++j)
          {
            std::complex<double> sum = 0;
            for (unsigned int q=0; q<quadrature.size(); ++q)
              {
                const Point<dim> x_q = quadrature.point(q);
                sum += std::exp(std::complex<double>(0,1) *
                                (k_vectors[k] * x_q)) *
                       fe_collection[fe].shape_value(j,x_q) *
                       quadrature.weight(q);
              }
            fourier_transform_matrices[fe](k,j)
              = sum / std::pow(2*numbers::PI, 1.*dim/2);
          }
    }

  std::vector<std::complex<double> > fourier_coefficients (n_fourier_modes);
  Vector<double>                     local_dof_values;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (unsigned int index=0; cell!=endc; ++cell, ++index)
    {
      local_dof_values.reinit (cell->get_fe().dofs_per_cell);
      cell->get_dof_values (solution, local_dof_values);

      for (unsigned int f=0; f<n_fourier_modes; ++f)
        {
          fourier_coefficients[f] = 0;

          for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
            fourier_coefficients[f] +=
              fourier_transform_matrices[cell->active_fe_index()](f,i)
              *
              local_dof_values(i);
        }

      std::map<unsigned int, double> k_to_max_U_map;
      for (unsigned int f=0; f<n_fourier_modes; ++f)
        if ((k_to_max_U_map.find (k_vectors_magnitude[f]) ==
             k_to_max_U_map.end())
            ||
            (k_to_max_U_map[k_vectors_magnitude[f]] <
             std::abs (fourier_coefficients[f])))
          k_to_max_U_map[k_vectors_magnitude[f]]
            = std::abs (fourier_coefficients[f]);
      double  sum_1           = 0,
              sum_ln_k        = 0,
              sum_ln_k_square = 0,
              sum_ln_U        = 0,
              sum_ln_U_ln_k   = 0;
      for (unsigned int f=0; f<n_fourier_modes; ++f)
        if (k_to_max_U_map[k_vectors_magnitude[f]] ==
            std::abs (fourier_coefficients[f]))
          {
            sum_1 += 1;
            sum_ln_k += ln_k[f];
            sum_ln_k_square += ln_k[f]*ln_k[f];
            sum_ln_U += std::log (std::abs (fourier_coefficients[f]));
            sum_ln_U_ln_k += std::log (std::abs (fourier_coefficients[f])) *
                             ln_k[f];
          }

      const double mu
        = (1./(sum_1*sum_ln_k_square - sum_ln_k*sum_ln_k)
           *
           (sum_ln_k*sum_ln_U - sum_1*sum_ln_U_ln_k));

      smoothness_indicators(index) = mu - 1.*dim/2;
    }
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    deallog.push("2d");
    LaplaceProblem<2> laplace_problem;
    laplace_problem.run ();
    deallog.pop();
  }

  {
    deallog.push("3d");
    LaplaceProblem<3> laplace_problem;
    laplace_problem.run ();
    deallog.pop();
  }
}

