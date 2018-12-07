// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// test SmoothnessEstimator::legendre_coefficient_decay() on
// problem 4 (peak) in Mitchel 2014.

#include "laplace.h"


template <int dim>
class ForcingFunction : public Function<dim>
{
public:
  ForcingFunction(const double alpha, const Point<dim> center)
    : Function<dim>(1)
    , alpha(alpha)
    , center(center)
  {}

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const;

private:
  const double     alpha;
  const Point<dim> center;
};

template <int dim>
double
ForcingFunction<dim>::value(const Point<dim> &point, const unsigned int) const
{
  const double x = point[0];
  const double y = point[1];

  return -exp(-alpha * (point - center).norm_square()) * 4 * alpha *
         (alpha * (point - center).norm_square() - 1);
}

template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution(const double alpha, const Point<dim> center)
    : Function<dim>(1)
    , alpha(alpha)
    , center(center){};

  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const;

  virtual Tensor<1, dim>
  gradient(const Point<dim> &point, const unsigned int component = 0) const;

private:
  const double     alpha;
  const Point<dim> center;
};

template <int dim>
double
ExactSolution<dim>::value(const Point<dim> &point, const unsigned int) const
{
  return exp(-alpha * ((point - center).norm_square()));
}

template <int dim>
Tensor<1, dim>
ExactSolution<dim>::gradient(const Point<dim> &point, const unsigned int) const
{
  Tensor<1, dim> grad_u = point - center;
  grad_u *= -2 * alpha * exp(-alpha * ((point - center).norm_square()));
  return grad_u;
}


template <int dim>
class Problem4 : public Laplace<dim>
{
public:
  Problem4(const Function<dim> &force_function,
           const Function<dim> &exact_solution,
           const Function<dim> &boundary_conditions,
           const unsigned int   n_cycles,
           const std::string    output_name);


private:
  void
  setup_geometry();
  void
  estimate_error();
  void
  mark_h_cells();

  std::pair<unsigned int, unsigned int>
  substitute_h_for_p(
    std::vector<typename Triangulation<dim>::active_cell_iterator> &p_cells);

  hp::QCollection<dim - 1> quadrature_face;
};

template <int dim>
Problem4<dim>::Problem4(const Function<dim> &force_function,
                        const Function<dim> &exact_solution,
                        const Function<dim> &boundary_conditions,
                        const unsigned int   n_cycles,
                        const std::string    output_name)
  : Laplace<dim>(force_function,
                 exact_solution,
                 boundary_conditions,
                 n_cycles,
                 output_name)
{
  for (unsigned int p = 1; p <= n_cycles; p++)
    {
      // Laplace<dim>::fe.push_back(FE_Q_Hierarchical<dim>(p));
      Laplace<dim>::fe.push_back(FE_Q<dim>(p));
      Laplace<dim>::quadrature.push_back(QSorted<dim>(QGauss<dim>(p + 1)));

      quadrature_face.push_back(QSorted<dim - 1>(QGauss<dim - 1>(p + 1)));

      const QTrapez<1>     q_trapez;
      const QIterated<dim> q_iterated(q_trapez, p + 3);
      Laplace<dim>::quadrature_infty.push_back(QSorted<dim>(q_iterated));
    }
}



template <int dim>
std::pair<unsigned int, unsigned int>
Problem4<dim>::substitute_h_for_p(
  std::vector<typename Triangulation<dim>::active_cell_iterator> &p_cells)
{
  Vector<float> smoothness_indicators(
    Laplace<dim>::triangulation.n_active_cells());
  SmoothnessEstimator::legendre_coefficient_decay(Laplace<dim>::dof_handler,
                                                  Laplace<dim>::solution,
                                                  smoothness_indicators);

  unsigned int num_p_cells = 0;
  unsigned int num_h_cells = 0;
  for (auto &cell : Laplace<dim>::dof_handler.active_cell_iterators())
    if (cell->refine_flag_set())
      {
        typename Triangulation<dim>::active_cell_iterator tria_cell(
          &(Laplace<dim>::triangulation), cell->level(), cell->index());

        const unsigned int cur_fe_index = cell->active_fe_index();
        const bool p_ref = smoothness_indicators(cell->index()) < exp(-1.);

        if (cur_fe_index < Laplace<dim>::fe.size() - 1 && p_ref)
          {
            ++num_p_cells;
            cell->clear_refine_flag();
            p_cells.push_back(tria_cell);
          }
        else
          {
            ++num_h_cells;
          }
      }

  return std::make_pair(num_h_cells, num_p_cells);
}



template <int dim>
void
Problem4<dim>::setup_geometry()
{
  std::vector<unsigned int> number_elements(2);
  number_elements[0] = 16;
  number_elements[1] = 16;

  GridGenerator::subdivided_hyper_rectangle(Laplace<dim>::triangulation,
                                            number_elements,
                                            Point<dim>(0, 0),
                                            Point<dim>(1, 1),
                                            false);
}



template <int dim>
void
Problem4<dim>::estimate_error()
{
  KellyErrorEstimator<dim>::estimate(
    Laplace<dim>::dof_handler,
    quadrature_face,
    std::map<types::boundary_id, const Function<dim> *>(),
    Laplace<dim>::solution,
    Laplace<dim>::estimated_error_per_cell);
}

template <int dim>
void
Problem4<dim>::mark_h_cells()
{
  GridRefinement::refine_and_coarsen_fixed_number(
    Laplace<dim>::triangulation,
    Laplace<dim>::estimated_error_per_cell,
    0.2,
    0.0);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const int dim = 2;

  initlog();

  // peak strength
  const double alpha = 1000;
  // peak position:
  const double xc = 0.5;
  const double yc = 0.5;

  const Point<dim> center(xc, yc);

  ForcingFunction<dim> ff(alpha, center);
  ExactSolution<dim>   ex(alpha, center);

  Problem4<dim> problem(ff, ex, ex, 10, "convergence");
  problem.run();
}
