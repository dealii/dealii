// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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



// a modified version of step-27 that crashes due to circular constraints

char logname[] = "output";


#include "../tests.h"


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>

#include <complex>

// Finally, this is as in previous
// programs:
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
  void assemble_system ();
  void solve ();
  void create_coarse_grid ();
  void refine_grid ();
  void estimate_smoothness (Vector<float> &smoothness_indicators) const;
  void output_results (const unsigned int cycle) const;

  Triangulation<dim>   triangulation;

  hp::DoFHandler<dim>      dof_handler;
  hp::FECollection<dim>    fe_collection;
  hp::QCollection<dim>     quadrature_collection;
  hp::QCollection<dim-1>   face_quadrature_collection;

  ConstraintMatrix     hanging_node_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;

  Timer distr, condense, hang, assemble, solver;
};



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
LaplaceProblem<dim>::LaplaceProblem () :
  dof_handler (triangulation)
{
  for (unsigned int degree=2; degree<(dim == 2 ? 8 : 5); ++degree)
    {
      fe_collection.push_back (FE_Q<dim>(degree));
      quadrature_collection.push_back (QGauss<dim>(degree+2));
      face_quadrature_collection.push_back (QGauss<dim-1>(degree+2));
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
  distr.reset();
  distr.start();
  dof_handler.distribute_dofs (fe_collection);
  distr.stop();

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  hanging_node_constraints.clear ();
  hang.reset();
  hang.start();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           hanging_node_constraints);

  hanging_node_constraints.close ();
  hang.stop();

  if (dim < 3)
    {
      sparsity_pattern.reinit (dof_handler.n_dofs(),
                               dof_handler.n_dofs(),
                               dof_handler.max_couplings_between_dofs());
      DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

      condense.reset();
      condense.start();
      hanging_node_constraints.condense (sparsity_pattern);
      condense.stop();
      sparsity_pattern.compress();
    }
  else
    {
      CompressedSparsityPattern csp (dof_handler.n_dofs(),
                                     dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern (dof_handler, csp);

      condense.reset();
      condense.start();
      hanging_node_constraints.condense (csp);
      condense.stop();

      sparsity_pattern.copy_from (csp);
    }

  system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  hp::FEValues<dim> hp_fe_values (fe_collection,
                                  quadrature_collection,
                                  update_values    |  update_gradients |
                                  update_q_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
      FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
      Vector<double>       cell_rhs (dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      cell_matrix = 0;
      cell_rhs = 0;

      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

      std::vector<double>  rhs_values (fe_values.n_quadrature_points);
      rhs_function.value_list (fe_values.get_quadrature_points(),
                               rhs_values);

      for (unsigned int q_point=0; q_point<fe_values.n_quadrature_points; ++q_point)
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

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

  condense.start();
  hanging_node_constraints.condense (system_matrix);
  hanging_node_constraints.condense (system_rhs);
  condense.stop();
  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}

template <int dim>
void LaplaceProblem<dim>::solve ()
{
  SolverControl           solver_control (system_rhs.size(),
                                          1e-8*system_rhs.l2_norm());
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
            preconditioner);

  hanging_node_constraints.distribute (solution);
}


unsigned int
int_pow (const unsigned int x,
         const unsigned int n)
{
  unsigned int p=1;
  for (unsigned int i=0; i<n; ++i)
    p *= x;
  return p;
}


template <int dim>
void
LaplaceProblem<dim>::
estimate_smoothness (Vector<float> &smoothness_indicators) const
{
  const unsigned int N = (dim == 2 ? 7 : 4);

  // form all the Fourier vectors
  // that we want to
  // consider. exclude k=0 to avoid
  // problems with |k|^{-mu} and also
  // logarithms of |k|
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


  // assemble the matrices that do
  // the Fourier transforms for each
  // of the finite elements we deal
  // with. note that these matrices
  // are complex-valued, so we can't
  // use FullMatrix
  QGauss<1>      base_quadrature (2);
  QIterated<dim> quadrature (base_quadrature, N);

  std::vector<Table<2,std::complex<double> > >
  fourier_transform_matrices (fe_collection.size());
  for (unsigned int fe=0; fe<fe_collection.size(); ++fe)
    {
      fourier_transform_matrices[fe].reinit (n_fourier_modes,
                                             fe_collection[fe].dofs_per_cell);

      for (unsigned int k=0; k<n_fourier_modes; ++k)
        for (unsigned int i=0; i<fe_collection[fe].dofs_per_cell; ++i)
          {
            std::complex<double> sum = 0;
            for (unsigned int q=0; q<quadrature.size(); ++q)
              {
                const Point<dim> x_q = quadrature.point(q);
                sum += std::exp(std::complex<double>(0,1) *
                                (k_vectors[k] * x_q)) *
                       fe_collection[fe].shape_value(i,x_q) *
                       quadrature.weight(q);
              }
            fourier_transform_matrices[fe](k,i)
              = sum / std::pow(2*numbers::PI, 1.*dim/2);
          }
    }

  // the next thing is to loop over
  // all cells and do our work there,
  // i.e. to locally do the Fourier
  // transform and estimate the decay
  // coefficient
  std::vector<std::complex<double> > fourier_coefficients (n_fourier_modes);
  Vector<double>                     local_dof_values;

  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (unsigned int index=0; cell!=endc; ++cell, ++index)
    {
      local_dof_values.reinit (cell->get_fe().dofs_per_cell);
      cell->get_dof_values (solution, local_dof_values);

      // first compute the Fourier
      // transform of the local
      // solution
      std::fill (fourier_coefficients.begin(), fourier_coefficients.end(), 0);
      for (unsigned int f=0; f<n_fourier_modes; ++f)
        for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
          fourier_coefficients[f] +=
            fourier_transform_matrices[cell->active_fe_index()](f,i)
            *
            local_dof_values(i);

      // enter the Fourier
      // coefficients into a map with
      // the k-magnitudes, to make
      // sure that we get only the
      // largest magnitude for each
      // value of |k|
      std::map<unsigned int, double> k_to_max_U_map;
      for (unsigned int f=0; f<n_fourier_modes; ++f)
        if ((k_to_max_U_map.find (k_vectors_magnitude[f]) ==
             k_to_max_U_map.end())
            ||
            (k_to_max_U_map[k_vectors_magnitude[f]] <
             std::abs (fourier_coefficients[f])))
          k_to_max_U_map[k_vectors_magnitude[f]]
            = std::abs (fourier_coefficients[f]);

      // now we have to calculate the
      // various contributions to the
      // formula for mu. we'll only
      // take those fourier
      // coefficients with the
      // largest value for a given
      // |k|
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
            sum_ln_U_ln_k += std::log (std::abs (fourier_coefficients[f])) * ln_k[f];
          }

      const double mu
        = (1./(sum_1*sum_ln_k_square - sum_ln_k*sum_ln_k)
           *
           (sum_ln_k*sum_ln_U - sum_1*sum_ln_U_ln_k));

      smoothness_indicators(index) = mu - 1.*dim/2;
    }
}




template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      face_quadrature_collection,
                                      typename FunctionMap<dim>::type(),
                                      solution,
                                      estimated_error_per_cell);

  Vector<float> smoothness_indicators (triangulation.n_active_cells());
  estimate_smoothness (smoothness_indicators);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);

  float max_smoothness = 0,
        min_smoothness = smoothness_indicators.linfty_norm();
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      if (cell->refine_flag_set())
        {
          max_smoothness = std::max (max_smoothness,
                                     smoothness_indicators(index));
          min_smoothness = std::min (min_smoothness,
                                     smoothness_indicators(index));
        }
  }
  const float cutoff_smoothness = (max_smoothness + min_smoothness) / 2;
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      if (cell->refine_flag_set()
          &&
          (smoothness_indicators(index) > cutoff_smoothness)
          &&
          !(cell->active_fe_index() == fe_collection.size() - 1))
        {
          cell->clear_refine_flag();
          cell->set_active_fe_index (std::min (cell->active_fe_index() + 1,
                                               fe_collection.size() - 1));
        }
  }

  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());

  {
    const std::string filename = "grid-" +
                                 Utilities::int_to_string (cycle, 2) +
                                 ".eps";
    std::ofstream output (filename.c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, output);
  }

  {
    Vector<float> smoothness_indicators (triangulation.n_active_cells());
    estimate_smoothness (smoothness_indicators);

    Vector<double> smoothness_field (dof_handler.n_dofs());
    DoFTools::distribute_cell_to_dof_vector (dof_handler,
                                             smoothness_indicators,
                                             smoothness_field);

    Vector<float> fe_indices (triangulation.n_active_cells());
    {
      typename hp::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (unsigned int index=0; cell!=endc; ++cell, ++index)
        {

          fe_indices(index) = cell->active_fe_index();
//    smoothness_indicators(index) *= std::sqrt(cell->diameter());
        }
    }

    const std::string filename = "solution-" +
                                 Utilities::int_to_string (cycle, 2) +
                                 ".vtk";
    DataOut<dim,hp::DoFHandler<dim> > data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.add_data_vector (smoothness_indicators, "smoothness1");
    data_out.add_data_vector (smoothness_field, "smoothness2");
    data_out.add_data_vector (fe_indices, "fe_index");
    data_out.build_patches ();

    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }
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
  triangulation.refine_global (3);
}



template <>
void LaplaceProblem<3>::create_coarse_grid ()
{
  const unsigned int dim = 3;

  static const Point<3> vertices_1[]
  =
  {
    // points on the lower surface
    Point<dim>(0,  0, -4),
    Point<dim>(std::cos(0*numbers::PI/6),
    std::sin(0*numbers::PI/6),
    -4),
    Point<dim>(std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    -4),
    Point<dim>(std::cos(4*numbers::PI/6),
    std::sin(4*numbers::PI/6),
    -4),
    Point<dim>(std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    -4),
    Point<dim>(std::cos(8*numbers::PI/6),
    std::sin(8*numbers::PI/6),
    -4),
    Point<dim>(std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    -4),

    // same points on the top
    // of the stem, with
    // indentation in the middle
    Point<dim>(0,  0, 4-std::sqrt(2.)/2),
    Point<dim>(std::cos(0*numbers::PI/6),
    std::sin(0*numbers::PI/6),
    4),
    Point<dim>(std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    4),
    Point<dim>(std::cos(4*numbers::PI/6),
    std::sin(4*numbers::PI/6),
    4),
    Point<dim>(std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    4),
    Point<dim>(std::cos(8*numbers::PI/6),
    std::sin(8*numbers::PI/6),
    4),
    Point<dim>(std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    4),

    // point at top of chevron
    Point<dim>(0,0,4+std::sqrt(2.)/2),

    // points at the top of the
    // first extension
    // points 15-18
    Point<dim>(0,  0, 7) + Point<dim> (std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(0*numbers::PI/6),
    std::sin(0*numbers::PI/6),
    7) + Point<dim> (std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    7) + Point<dim> (std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(4*numbers::PI/6),
    std::sin(4*numbers::PI/6),
    7) + Point<dim> (std::cos(2*numbers::PI/6),
    std::sin(2*numbers::PI/6),
    0) * 4,

    // points at the top of the
    // second extension
    // points 19-22
    Point<dim>(0,  0, 7) + Point<dim> (std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(4*numbers::PI/6),
    std::sin(4*numbers::PI/6),
    7) + Point<dim> (std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    7) + Point<dim> (std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(8*numbers::PI/6),
    std::sin(8*numbers::PI/6),
    7) + Point<dim> (std::cos(6*numbers::PI/6),
    std::sin(6*numbers::PI/6),
    0) * 4,

    // points at the top of the
    // third extension
    // points 23-26
    Point<dim>(0,  0, 7) + Point<dim> (std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(8*numbers::PI/6),
    std::sin(8*numbers::PI/6),
    7) + Point<dim> (std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    7) + Point<dim> (std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    0) * 4,
    Point<dim>(std::cos(0*numbers::PI/6),
    std::sin(0*numbers::PI/6),
    7) + Point<dim> (std::cos(10*numbers::PI/6),
    std::sin(10*numbers::PI/6),
    0) * 4,

  };

  const unsigned int
  n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
  const std::vector<Point<dim> > vertices (&vertices_1[0],
                                           &vertices_1[n_vertices]);
  static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
  =
  {
    // the three cells in the stem
    {0, 2, 4, 3, 7, 9, 11, 10},
    {6, 0, 5, 4, 13, 7, 12, 11},
    {6, 1, 0, 2, 13, 8, 7, 9},
    // the chevron at the center
    {13, 8, 7, 9, 12, 14, 11, 10},
    // first extension
    {14, 8, 10, 9, 15, 16, 18, 17},
    // second extension
    {11, 12, 10, 14, 21, 22, 20, 19},
    // third extension
    {12, 13, 14, 8, 24, 25, 23, 26},
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

  for (Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if ((cell->face(f)->center()[2] != -4)
          &&
          (cell->face(f)->center()[2] != 7)
          &&
          (cell->face(f)->at_boundary()))
        cell->face(f)->set_boundary_indicator (1);

  triangulation.refine_global (1);
}



template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<30; ++cycle)
    {
      deallog << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        create_coarse_grid ();
      else
        refine_grid ();


      deallog << "   Number of active cells:       "
              << triangulation.n_active_cells()
              << std::endl;

      Timer all;
      all.start();
      setup_system ();

      deallog << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;
      deallog << "   Number of constraints       : "
              << hanging_node_constraints.n_constraints()
              << std::endl;

      assemble.reset ();
      assemble.start ();
      assemble_system ();
      assemble.stop();


      solver.reset();
      solver.start();
      solve ();
      solver.stop();

      all.stop();
    }
}

int main ()
{
  std::ofstream logfile(logname);
  logfile.precision (3);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  try
    {
      LaplaceProblem<3> laplace_problem;
      laplace_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
