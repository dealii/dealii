// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Testing the ColorEnriched::Helper class by solving problems
 */

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_cspline.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/slepc_solver.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <math.h>

#include <map>
#include <set>
#include <vector>

#include "../tests.h"

template <int dim>
class SigmaFunction : public Function<dim>
{
  Point<dim>          center;
  FunctionParser<dim> func;

public:
  SigmaFunction()
    : Function<dim>()
    , func(1)
  {}

  // to help with resize function. doesn't copy function parser(func)!
  SigmaFunction(SigmaFunction &&other)
    : center(other.center)
    , func(1)
  {}

  void
  initialize(const Point<dim>  &center,
             const double      &sigma,
             const std::string &func_expr);
  double
  value(const Point<dim> &p, const unsigned int component = 0) const;
  Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &value_list) const;
};

template <int dim>
void
SigmaFunction<dim>::initialize(const Point<dim>  &_center,
                               const double      &sigma,
                               const std::string &func_expr)
{
  center = _center;
  std::string                   variables;
  std::map<std::string, double> constants = {{"sigma", sigma},
                                             {"pi", numbers::PI}};

  AssertThrow(dim == 1 || dim == 2 || dim == 3,
              ExcMessage("Dimension not implemented"));
  switch (dim)
    {
      case 1:
        variables = "x";
        break;
      case 2:
        variables = "x, y";
        break;
      case 3:
        variables = "x, y, z";
        break;
    }

  func.initialize(variables, func_expr, constants);
}

template <int dim>
inline double
SigmaFunction<dim>::value(const Point<dim>  &p,
                          const unsigned int component) const
{
  const Point<dim> d(p - center);
  return func.value(d, component);
}


template <int dim>
inline Tensor<1, dim>
SigmaFunction<dim>::gradient(const Point<dim>  &p,
                             const unsigned int component) const
{
  const Point<dim> d(p - center);
  return func.gradient(d, component);
}

template <int dim>
void
SigmaFunction<dim>::value_list(const std::vector<Point<dim>> &points,
                               std::vector<double>           &value_list) const
{
  const unsigned int n_points = points.size();

  AssertDimension(points.size(), value_list.size());

  for (unsigned int p = 0; p < n_points; ++p)
    value_list[p] = value(points[p]);
}



template <int dim>
struct EnrichmentPredicate
{
  EnrichmentPredicate(const Point<dim> origin, const double radius)
    : origin(origin)
    , radius(radius)
  {}

  template <class Iterator>
  bool
  operator()(const Iterator &i) const
  {
    return ((i->center() - origin).norm_square() < radius * radius);
  }

  const Point<dim> &
  get_origin()
  {
    return origin;
  }
  const double &
  get_radius()
  {
    return radius;
  }

private:
  const Point<dim> origin;
  const double     radius;
};


template <int dim>
class SplineEnrichmentFunction : public Function<dim>
{
public:
  SplineEnrichmentFunction(const Point<dim>          &origin,
                           const std::vector<double> &interpolation_points_1d,
                           const std::vector<double> &interpolation_values_1d)
    : Function<dim>(1)
    , origin(origin)
    , interpolation_points(interpolation_points_1d)
    , interpolation_values(interpolation_values_1d)
    , cspline(interpolation_points, interpolation_values)
  {}

  SplineEnrichmentFunction(SplineEnrichmentFunction &&other)
    : Function<dim>(1)
    , origin(other.origin)
    , interpolation_points(other.interpolation_points)
    , interpolation_values(other.interpolation_values)
    , cspline(interpolation_points, interpolation_values)
  {}

  SplineEnrichmentFunction(const SplineEnrichmentFunction &other)
    : Function<dim>(1)
    , origin(other.origin)
    , interpolation_points(other.interpolation_points)
    , interpolation_values(other.interpolation_values)
    , cspline(interpolation_points, interpolation_values)
  {}



  virtual double
  value(const Point<dim> &point, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = point - origin;
    const double   r    = dist.norm();
    return cspline.value(Point<1>(r), component);
  }

  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const
  {
    Tensor<1, dim> dist = p - origin;
    const double   r    = dist.norm();
    Assert(r > 0., ExcDivideByZero());
    dist /= r;
    Assert(component == 0, ExcMessage("Not implemented"));
    return cspline.gradient(Point<1>(r))[0] * dist;
  }

private:
  /**
   * origin
   */
  const Point<dim>    origin;
  std::vector<double> interpolation_points;
  std::vector<double> interpolation_values;
  // enrichment function as CSpline based on radius
  Functions::CSpline<1> cspline;
};



struct ParameterCollection
{
  ParameterCollection(const std::string &file_name);

  ParameterCollection(const int                 &dim,
                      const double              &size,
                      const unsigned int        &shape,
                      const unsigned int        &global_refinement,
                      const unsigned int        &cycles,
                      const unsigned int        &fe_base_degree,
                      const unsigned int        &fe_enriched_degree,
                      const unsigned int        &max_iterations,
                      const double              &tolerance,
                      const std::string         &rhs_value_expr,
                      const std::string         &boundary_value_expr,
                      const std::string         &rhs_radial_problem,
                      const std::string         &boundary_radial_problem,
                      const std::string         &exact_soln_expr,
                      const unsigned int        &patches,
                      const unsigned int        &debug_level,
                      const unsigned int        &n_enrichments,
                      const std::vector<double> &points_enrichments,
                      const std::vector<double> &radii_predicates,
                      const std::vector<double> &sigmas);

  void
  print();

  void
  set_enrichment_point(Point<2> &p, const unsigned int i)
  {
    AssertDimension(dim, 2);
    p(0) = points_enrichments[2 * i];
    p(1) = points_enrichments[2 * i + 1];
  }
  void
  set_enrichment_point(Point<3> &p, const unsigned int i)
  {
    AssertDimension(dim, 3);
    p(0) = points_enrichments[3 * i];
    p(1) = points_enrichments[3 * i + 1];
    p(2) = points_enrichments[3 * i + 2];
  }

  int          dim;
  double       size;
  unsigned int shape; // 0 = ball, 1 = cube
  unsigned int global_refinement;
  unsigned int cycles;
  unsigned int fe_base_degree;
  unsigned int fe_enriched_degree;
  unsigned int max_iterations;
  double       tolerance;

  // parameters related to exact solution
  std::string rhs_value_expr;
  std::string boundary_value_expr;

  // value = true ==> estimate exact solution from radial problem
  std::string rhs_radial_problem;
  std::string boundary_radial_problem;

  std::string exact_soln_expr;

  unsigned int patches;
  // debug level = 0(output nothing),
  // 1 (print statements)
  // 2 (output solution)
  // 3 (+ output grid data as well)
  // 9 (+ shape functions as well)
  unsigned int        debug_level;
  unsigned int        n_enrichments;
  std::vector<double> points_enrichments;
  std::vector<double> radii_predicates;
  std::vector<double> sigmas;
};



ParameterCollection::ParameterCollection(const std::string &file_name)
{
  // std::cout << "...reading parameters" << std::endl;

  ParameterHandler prm;

  // declare parameters
  prm.enter_subsection("geometry");
  prm.declare_entry("dim", "2", Patterns::Integer());
  prm.declare_entry("size", "1", Patterns::Double(0));
  prm.declare_entry("shape", "1", Patterns::Integer(0));
  prm.declare_entry("Global refinement", "1", Patterns::Integer(1));
  prm.declare_entry("cycles", "0", Patterns::Integer(0));
  prm.leave_subsection();


  prm.enter_subsection("solver");
  prm.declare_entry("fe base degree", "1", Patterns::Integer(1));
  prm.declare_entry("fe enriched degree", "1", Patterns::Integer(1));
  prm.declare_entry("max iterations", "1000", Patterns::Integer(1));
  prm.declare_entry("tolerance", "1e-8", Patterns::Double(0));
  prm.leave_subsection();


  prm.enter_subsection("expressions");
  prm.declare_entry("rhs value", "0", Patterns::Anything());
  prm.declare_entry("boundary value", "0", Patterns::Anything());
  prm.declare_entry("rhs value radial problem", "0", Patterns::Anything());
  prm.declare_entry("boundary value radial problem", "0", Patterns::Anything());
  prm.declare_entry("exact solution expression", "", Patterns::Anything());
  prm.declare_entry("estimate exact solution", "false", Patterns::Bool());
  prm.leave_subsection();


  prm.enter_subsection("output");
  prm.declare_entry("patches", "1", Patterns::Integer(1));
  prm.declare_entry("debug level", "0", Patterns::Integer(0, 9));
  prm.leave_subsection();

  // parse parameter file
  prm.parse_input(file_name, "#end-of-dealii parser");


  // get parameters
  prm.enter_subsection("geometry");
  dim               = prm.get_integer("dim");
  size              = prm.get_double("size");
  shape             = prm.get_integer("shape");
  global_refinement = prm.get_integer("Global refinement");
  cycles            = prm.get_integer("cycles");
  prm.leave_subsection();

  prm.enter_subsection("solver");
  fe_base_degree     = prm.get_integer("fe base degree");
  fe_enriched_degree = prm.get_integer("fe enriched degree");
  max_iterations     = prm.get_integer("max iterations");
  tolerance          = prm.get_double("tolerance");
  prm.leave_subsection();

  prm.enter_subsection("expressions");
  rhs_value_expr          = prm.get("rhs value");
  boundary_value_expr     = prm.get("boundary value");
  rhs_radial_problem      = prm.get("rhs value radial problem");
  boundary_radial_problem = prm.get("boundary value radial problem");
  exact_soln_expr         = prm.get("exact solution expression");
  prm.leave_subsection();

  prm.enter_subsection("output");
  patches     = prm.get_integer("patches");
  debug_level = prm.get_integer("debug level");
  prm.leave_subsection();


  // manual parsing
  // open parameter file
  std::ifstream prm_file(file_name);

  // read lines until "#end-of-dealii parser" is reached
  std::string line;
  while (getline(prm_file, line))
    if (line == "#end-of-dealii parser")
      break;

  AssertThrow(line == "#end-of-dealii parser",
              ExcMessage(
                "line missing in parameter file = \'#end-of-dealii parser\' "));

  // function to read next line not starting with # or empty
  auto read_next_proper_line = [&](std::string &line) {
    while (getline(prm_file, line))
      {
        if (line.size() == 0 || line[0] == '#' || line[0] == ' ')
          continue;
        else
          break;
      }
  };

  std::stringstream s_stream;

  // read num of enrichement points
  read_next_proper_line(line);
  s_stream.str(line);
  s_stream >> n_enrichments;

  // note vector of points
  for (unsigned int i = 0; i != n_enrichments; ++i)
    {
      read_next_proper_line(line);
      s_stream.clear();
      s_stream.str(line);

      points_enrichments.resize(dim * n_enrichments);

      if (dim == 2)
        {
          double x, y;
          s_stream >> x >> y;
          points_enrichments[2 * i]     = x;
          points_enrichments[2 * i + 1] = y;
        }
      else if (dim == 3)
        {
          double x, y, z;
          s_stream >> x >> y >> z;
          points_enrichments[3 * i]     = x;
          points_enrichments[3 * i + 1] = y;
          points_enrichments[3 * i + 2] = z;
        }
      else
        AssertThrow(false, ExcMessage("Dimension not implemented"));
    }

  // note vector of radii for predicates
  for (unsigned int i = 0; i != n_enrichments; ++i)
    {
      read_next_proper_line(line);
      s_stream.clear();
      s_stream.str(line);

      double r;
      s_stream >> r;
      radii_predicates.push_back(r);
    }

  // note vector of sigmas for rhs
  for (unsigned int i = 0; i != n_enrichments; ++i)
    {
      read_next_proper_line(line);
      s_stream.clear();
      s_stream.str(line);

      double r;
      s_stream >> r;
      sigmas.push_back(r);
    }
}



ParameterCollection::ParameterCollection(
  const int                 &dim,
  const double              &size,
  const unsigned int        &shape,
  const unsigned int        &global_refinement,
  const unsigned int        &cycles,
  const unsigned int        &fe_base_degree,
  const unsigned int        &fe_enriched_degree,
  const unsigned int        &max_iterations,
  const double              &tolerance,
  const std::string         &rhs_value_expr,
  const std::string         &boundary_value_expr,
  const std::string         &rhs_radial_problem,
  const std::string         &boundary_radial_problem,
  const std::string         &exact_soln_expr,
  const unsigned int        &patches,
  const unsigned int        &debug_level,
  const unsigned int        &n_enrichments,
  const std::vector<double> &points_enrichments,
  const std::vector<double> &radii_predicates,
  const std::vector<double> &sigmas)
  : dim(dim)
  , size(size)
  , shape(shape)
  , global_refinement(global_refinement)
  , cycles(cycles)
  , fe_base_degree(fe_base_degree)
  , fe_enriched_degree(fe_enriched_degree)
  , max_iterations(max_iterations)
  , tolerance(tolerance)
  , rhs_value_expr(rhs_value_expr)
  , boundary_value_expr(boundary_value_expr)
  , rhs_radial_problem(rhs_radial_problem)
  , boundary_radial_problem(boundary_radial_problem)
  , exact_soln_expr(exact_soln_expr)
  , patches(patches)
  , debug_level(debug_level)
  , n_enrichments(n_enrichments)
  , points_enrichments(points_enrichments)
  , radii_predicates(radii_predicates)
  , sigmas(sigmas)
{}



void
ParameterCollection::print()
{
  std::cout << "Dim : " << dim << std::endl
            << "Size : " << size << std::endl
            << "Shape : " << shape << std::endl
            << "Global refinement : " << global_refinement << std::endl
            << "Cycles : " << cycles << std::endl
            << "FE base degree : " << fe_base_degree << std::endl
            << "FE enriched degree : " << fe_enriched_degree << std::endl
            << "Max Iterations : " << max_iterations << std::endl
            << "Tolerance : " << tolerance << std::endl
            << "rhs - main problem : " << rhs_value_expr << std::endl
            << "boundary value - main problem : " << boundary_value_expr
            << std::endl
            << "rhs of radial problem : " << rhs_radial_problem << std::endl
            << "boundary value of radial problem : " << boundary_radial_problem
            << std::endl
            << "exact solution expr : " << exact_soln_expr << std::endl
            << "estimate exact solution using radial problem : "
            << "Patches used for output: " << patches << std::endl
            << "Debug level: " << debug_level << std::endl
            << "Number of enrichments: " << n_enrichments << std::endl;

  std::cout << "Enrichment points : " << std::endl;
  for (unsigned int i = 0; i < points_enrichments.size(); i = i + dim)
    {
      for (int d = 0; d < dim; ++d)
        std::cout << points_enrichments[i + d] << ' ';

      std::cout << std::endl;
    }

  std::cout << "Enrichment radii : " << std::endl;
  for (auto r : radii_predicates)
    std::cout << r << std::endl;

  std::cout << "Sigma values of different sources : " << std::endl;
  for (auto r : sigmas)
    std::cout << r << std::endl;
}



/*
 * EstimateEnrichmentFunction is used to estimate enrichment function by
 * solving a 1D poisson problem with right hand side and boundary
 * expression provided as a string of single variable 'x', to be
 * interpreted as distance from @par center.
 *
 * Eg: For a 2D poisson problem with right hand side expression R and boundary
 * expression B given by functions R( x^2 + y^2) and B( x^2 + y^2 ), an
 * equivalent radial problem can be solved using this class by setting
 * @par rhs_expr = R( x^2)
 * @par boundary_expr = B (x^2)
 * Note that the original poisson problem is defined by right hand side
 * and boundary expression dependent on square of distance (x^2 + y^2)
 * from center.
 */
class EstimateEnrichmentFunction
{
public:
  EstimateEnrichmentFunction(const Point<1>    &center,
                             const double      &domain_size,
                             const double      &sigma,
                             const std::string &rhs_expr,
                             const std::string &boundary_expr,
                             const unsigned int refinement = 11);
  EstimateEnrichmentFunction(const Point<1>    &center,
                             const double      &left_bound,
                             const double      &right_bound,
                             const double      &sigma,
                             const std::string &rhs_expr,
                             const std::string &boundary_expr,
                             const unsigned int refinement = 11);
  ~EstimateEnrichmentFunction();
  void
  run();
  void
  evaluate_at_x_values(std::vector<double> &interpolation_points,
                       std::vector<double> &interpolation_values);
  double
  value(const Point<1> &p, const unsigned int &component = 0);

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  refine_grid();
  void
                      output_results() const;
  Point<1>            center;
  double              domain_size;
  double              left_bound, right_bound;
  double              sigma;
  std::string         rhs_expr;
  std::string         boundary_expr;
  std::vector<double> rhs_values;

public:
  unsigned int debug_level;

private:
  Triangulation<1>     triangulation;
  unsigned int         refinement;
  FE_Q<1>              fe;
  DoFHandler<1>        dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;
};

EstimateEnrichmentFunction::EstimateEnrichmentFunction(
  const Point<1>    &center,
  const double      &domain_size,
  const double      &sigma,
  const std::string &rhs_expr,
  const std::string &boundary_expr,
  const unsigned int refinement)
  : center(center)
  , domain_size(domain_size)
  , sigma(sigma)
  , rhs_expr(rhs_expr)
  , boundary_expr(boundary_expr)
  , debug_level(0)
  , refinement(refinement)
  , fe(1)
  , dof_handler(triangulation)
{
  left_bound  = center[0] - domain_size / 2;
  right_bound = center[0] + domain_size / 2;
}


EstimateEnrichmentFunction::EstimateEnrichmentFunction(
  const Point<1>    &center,
  const double      &left_bound,
  const double      &right_bound,
  const double      &sigma,
  const std::string &rhs_expr,
  const std::string &boundary_expr,
  const unsigned int refinement)
  : center(center)
  , left_bound(left_bound)
  , right_bound(right_bound)
  , sigma(sigma)
  , rhs_expr(rhs_expr)
  , boundary_expr(boundary_expr)
  , debug_level(0)
  , refinement(refinement)
  , fe(1)
  , dof_handler(triangulation)
{
  domain_size = right_bound - left_bound;
}


void
EstimateEnrichmentFunction::make_grid()
{
  GridGenerator::hyper_cube(triangulation, left_bound, right_bound);
  triangulation.refine_global(refinement);
}


void
EstimateEnrichmentFunction::setup_system()
{
  dof_handler.distribute_dofs(fe);
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


void
EstimateEnrichmentFunction::assemble_system()
{
  QGauss<1>        quadrature_formula(2);
  SigmaFunction<1> rhs;
  rhs.initialize(center, sigma, rhs_expr);
  FEValues<1>        fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      rhs_values.resize(n_q_points);
      rhs.value_list(fe_values.get_quadrature_points(), rhs_values);

      for (const auto q_index : fe_values.quadrature_point_indices())
        {
          double radius = center.distance(fe_values.quadrature_point(q_index));

          //-1/r (r*u_r) = f form converts to
          // r(u_r, v_r) = (r*f,v)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                cell_matrix(i, j) +=
                  (radius * fe_values.shape_grad(i, q_index) *
                   fe_values.shape_grad(j, q_index)) *
                  fe_values.JxW(q_index);
              cell_rhs(i) +=
                radius * (fe_values.shape_value(i, q_index) *
                          rhs_values[q_index] * fe_values.JxW(q_index));
            }
        }
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }
  std::map<types::global_dof_index, double> boundary_values;
  SigmaFunction<1>                          boundary_func;
  boundary_func.initialize(center, sigma, boundary_expr);

  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           boundary_func,
                                           boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           1,
                                           boundary_func,
                                           boundary_values);

  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


void
EstimateEnrichmentFunction::solve()
{
  SolverControl solver_control(50000, 1e-12, false, false);
  SolverCG<>    solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}


void
EstimateEnrichmentFunction::refine_grid()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  KellyErrorEstimator<1>::estimate(
    dof_handler,
    QGauss<1 - 1>(3),
    std::map<types::boundary_id, const Function<1> *>{},
    solution,
    estimated_error_per_cell);
  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.2,
                                                  0.01);
  triangulation.execute_coarsening_and_refinement();
}


void
EstimateEnrichmentFunction::output_results() const
{
  DataOut<1> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::ofstream output("radial_solution.vtk");
  data_out.write_vtk(output);
}


void
EstimateEnrichmentFunction::run()
{
  if (debug_level >= 1)
    std::cout << "Solving problem in 1.: " << 1 << " with center: " << center
              << ", size: " << domain_size << ", sigma: " << sigma << std::endl;

  make_grid();

  double old_value = 0, value = 1, relative_change = 1;
  bool   start = true;
  do
    {
      if (!start)
        {
          refine_grid();
          ++refinement;
        }

      if (debug_level >= 1)
        std::cout << "Refinement level: " << refinement << std::endl;

      setup_system();
      assemble_system();
      solve();

      value = VectorTools::point_value(dof_handler, solution, center);
      if (!start)
        {
          relative_change = fabs((old_value - value) / old_value);
        }
      start     = false;
      old_value = value;
    }
  while (relative_change > 0.005);

  if (debug_level >= 1)
    std::cout << "Radial solution at origin = " << value
              << " after global refinement " << refinement << std::endl;

  if (debug_level >= 1)
    output_results();
}


void
EstimateEnrichmentFunction::evaluate_at_x_values(
  std::vector<double> &interpolation_points,
  std::vector<double> &interpolation_values)
{
  if (interpolation_values.size() != interpolation_points.size())
    interpolation_values.resize(interpolation_points.size());

  // x varies from 0 to 2*sigma.
  // factor 2 because once a cell is decided to be enriched based on its center,
  // its quadrature points can cause x to be twice!
  for (unsigned int i = 0; i != interpolation_values.size(); ++i)
    {
      double value =
        VectorTools::point_value(dof_handler,
                                 solution,
                                 Point<1>(interpolation_points[i]));
      interpolation_values[i] = value;
    }
}


double
EstimateEnrichmentFunction::value(const Point<1>     &p,
                                  const unsigned int &component)
{
  return VectorTools::point_value(dof_handler, solution, p);
}


EstimateEnrichmentFunction::~EstimateEnrichmentFunction()
{
  triangulation.clear();
}



template <int dim>
void
plot_shape_function(DoFHandler<dim> &dof_handler, unsigned int patches = 5)
{
  std::cout << "...start plotting shape function" << std::endl;
  std::cout << "Patches for output: " << patches << std::endl;

  AffineConstraints<double> constraints;
  constraints.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // find set of dofs which belong to enriched cells
  std::set<unsigned int> enriched_cell_dofs;
  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->active_fe_index() != 0)
      {
        unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        enriched_cell_dofs.insert(local_dof_indices.begin(),
                                  local_dof_indices.end());
      }

  // output to check if all is good:
  std::vector<Vector<double>> shape_functions;
  std::vector<std::string>    names;
  for (auto dof : enriched_cell_dofs)
    {
      Vector<double> shape_function;
      shape_function.reinit(dof_handler.n_dofs());
      shape_function[dof] = 1.0;

      // if the dof is constrained, first output unconstrained vector
      names.push_back(std::string("C_") +
                      dealii::Utilities::int_to_string(dof, 2));
      shape_functions.push_back(shape_function);

      //      names.push_back(std::string("UC_") +
      //                      dealii::Utilities::int_to_string(s,2));

      //      // make continuous/constraint:
      //      constraints.distribute(shape_function);
      //      shape_functions.push_back(shape_function);
    }

  if (dof_handler.n_dofs() < 100)
    {
      std::cout << "...start printing support points" << std::endl;

      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      hp::MappingCollection<dim>                    hp_mapping;
      for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
        hp_mapping.push_back(mapping);
      DoFTools::map_dofs_to_support_points(hp_mapping,
                                           dof_handler,
                                           support_points);

      const std::string base_filename =
        "DOFs" + dealii::Utilities::int_to_string(dim) + "_p" +
        dealii::Utilities::int_to_string(0);

      const std::string filename = base_filename + ".gp";
      std::ofstream     f(filename);

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"" << base_filename << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics                                                                                   "
        << std::endl
        << "unset ytics" << std::endl
        << "unset grid" << std::endl
        << "unset border" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
        << std::endl;
      GridOut grid_out;
      grid_out.write_gnuplot(dof_handler.get_triangulation(), f);
      f << 'e' << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

      f << 'e' << std::endl;

      std::cout << "...finished printing support points" << std::endl;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get material ids:
  Vector<float> fe_index(dof_handler.get_triangulation().n_active_cells());
  for (auto &cell : dof_handler.active_cell_iterators())
    {
      fe_index[cell->active_cell_index()] = cell->active_fe_index();
    }
  data_out.add_data_vector(fe_index, "fe_index");

  for (unsigned int i = 0; i < shape_functions.size(); ++i)
    data_out.add_data_vector(shape_functions[i], names[i]);

  data_out.build_patches(patches);

  std::string   filename = "shape_functions.vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);

  std::cout << "...finished plotting shape functions" << std::endl;
}



template <int dim>
using predicate_function =
  std::function<bool(const typename Triangulation<dim>::cell_iterator &)>;



/**
 * Main class
 */
template <int dim>
class LaplaceProblem
{
public:
  LaplaceProblem(const ParameterCollection &prm);
  virtual ~LaplaceProblem();
  void
  run();

protected:
  void
  initialize();
  void
  build_fe_space();
  virtual void
  make_enrichment_functions();
  void
  setup_system();

private:
  void
  build_tables();
  void
  assemble_system();
  unsigned int
  solve();
  void
  refine_grid();
  void
  output_results(const unsigned int cycle);
  void
  process_solution();

protected:
  ParameterCollection prm;
  unsigned int        n_enriched_cells;

  Triangulation<dim> triangulation;
  DoFHandler<dim>    dof_handler;

  std::shared_ptr<const hp::FECollection<dim>> fe_collection;
  hp::QCollection<dim>                         q_collection;

  FE_Q<dim>       fe_base;
  FE_Q<dim>       fe_enriched;
  FE_Nothing<dim> fe_nothing;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  PETScWrappers::MPI::SparseMatrix system_matrix;
  PETScWrappers::MPI::Vector       solution;
  Vector<double>                   localized_solution;
  PETScWrappers::MPI::Vector       system_rhs;

  AffineConstraints<double> constraints;

  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  ConditionalOStream pcout;

  std::vector<SigmaFunction<dim>> vec_rhs;

  using cell_iterator_function = std::function<Function<dim> *(
    const typename DoFHandler<dim>::active_cell_iterator &)>;

  std::vector<std::shared_ptr<Function<dim>>> vec_enrichments;
  std::vector<predicate_function<dim>>        vec_predicates;
  std::vector<cell_iterator_function>         color_enrichments;

  // output vectors. size triangulation.n_active_cells()
  // change to Vector
  std::vector<Vector<float>> predicate_output;
  Vector<float>              color_output;
  Vector<float>              vec_fe_index;
  Vector<float>              mat_id;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem(const ParameterCollection &_par)
  : prm(_par)
  , n_enriched_cells(0)
  , dof_handler(triangulation)
  , fe_base(prm.fe_base_degree)
  , fe_enriched(prm.fe_enriched_degree)
  , fe_nothing(1, true)
  , mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, (this_mpi_process == 0) && (prm.debug_level >= 1))
{
  AssertThrow(prm.dim == dim, ExcMessage("parameter file dim != problem dim"));
  prm.print();
  pcout << "...parameters set" << std::endl;
}


/*
 * Construct basic grid, vector of predicate functions and
 * right hand side function of the problem.
 */
template <int dim>
void
LaplaceProblem<dim>::initialize()
{
  pcout << "...Start initializing" << std::endl;

  /*
   * set up basic grid which is a hyper cube or hyper ball based on
   * parameter file. Refine as per the global refinement value in the
   * parameter file.
   */
  if (prm.shape == 1)
    GridGenerator::hyper_cube(triangulation, -prm.size / 2.0, prm.size / 2.0);
  else if (prm.shape == 0)
    {
      Point<dim> center = Point<dim>();
      GridGenerator::hyper_ball(triangulation, center, prm.size / 2.0);
      triangulation.set_all_manifold_ids_on_boundary(0);
      static SphericalManifold<dim> spherical_manifold(center);
      triangulation.set_manifold(0, spherical_manifold);
    }
  else
    AssertThrow(false, ExcMessage("Shape not implemented."));
  triangulation.refine_global(prm.global_refinement);

  /*
   * Ensure that num of radii, sigma and focal points of enrichment domains
   * provided matches the number of predicates. Since focal points of
   * enrichment domains are stored as vector of numbers, the dimension of
   * points_enrichments is dim times the n_enrichments in prm file.
   */
  Assert(
    prm.points_enrichments.size() / dim == prm.n_enrichments &&
      prm.radii_predicates.size() == prm.n_enrichments &&
      prm.sigmas.size() == prm.n_enrichments,
    ExcMessage(
      "Number of enrichment points, predicate radii and sigmas should be equal"));

  /*
   * Construct vector of predicate functions, where a function at index i
   * returns true for a given cell if it belongs to enrichment domain i.
   * The decision is based on cell center being at a distance within radius
   * of the domain from focal point of the domain.
   */
  for (unsigned int i = 0; i != prm.n_enrichments; ++i)
    {
      Point<dim> p;
      prm.set_enrichment_point(p, i);
      vec_predicates.push_back(
        EnrichmentPredicate<dim>(p, prm.radii_predicates[i]));
    }

  /*
   * Construct a vector of right hand side functions where the actual right hand
   * side of problem is evaluated through summation of contributions from each
   * of the functions in the vector. Each function i at a given point is
   * evaluated with respect to point's position {x_i, y_i, z_i} relative to
   * focal point of enrichment domain i. The value is then a function of x_i,
   * y_i, z_i and sigmas[i] given by the rhs_value_expr[i] in parameter file.
   */
  vec_rhs.resize(prm.n_enrichments);
  for (unsigned int i = 0; i != prm.n_enrichments; ++i)
    {
      Point<dim> p;
      prm.set_enrichment_point(p, i);
      vec_rhs[i].initialize(p, prm.sigmas[i], prm.rhs_value_expr);
    }

  pcout << "...finish initializing" << std::endl;
}



/*
 * Create a enrichment function associated with each enrichment domain.
 * Here the approximate solution is assumed to be only dependent on
 * distance (r) from the focal point of enrichment domain. This is
 * true for poisson problem (a linear PDE) provided the right hand side is a
 * radial function and drops exponentially for points farther from focal
 * point of enrichment domain.
 */
template <int dim>
void
LaplaceProblem<dim>::make_enrichment_functions()
{
  pcout << "!!! Make enrichment function called" << std::endl;

  for (unsigned int i = 0; i < vec_predicates.size(); ++i)
    {
      /*
       * Formulate a 1d/radial problem with center and size appropriate
       * to cover the enrichment domain. The function determining
       * right hand side and boundary value of the problem is provided in
       * parameter file. The sigma governing this function is the same as
       * sigma provided for the corresponding enrichment domain and predicate
       * function for which the enrichment function is to be estimated.
       *
       * The center is 0 since the function is evaluated with respect to
       * the relative position from the focal point anyway.
       *
       * For hexahedral cells, the dimension can extend up to sqrt(3) < 2 times!
       * So take a factor of 4 as size of the problem. This ensures that the
       * enrichment function can be evaluated at all points in the enrichment
       * domain.
       */
      double center = 0;
      double sigma  = prm.sigmas[i];
      double size   = prm.radii_predicates[i] * 4;

      /*
       * Radial problem is solved only when enrichment domain has a positive
       * radius and hence non-empty domain.
       */
      if (prm.radii_predicates[i] != 0)
        {
          EstimateEnrichmentFunction radial_problem(
            Point<1>(center),
            size,
            sigma,
            prm.rhs_radial_problem,
            prm.boundary_radial_problem);
          radial_problem.debug_level = prm.debug_level; // print output
          radial_problem.run();
          pcout << "solved problem with "
                << "x and sigma : " << center << ", " << sigma << std::endl;

          // make points at which solution needs to interpolated
          std::vector<double> interpolation_points, interpolation_values;
          double              cut_point = 3 * sigma;
          unsigned int        n1 = 15, n2 = 15;
          double              radius      = size / 2;
          double              right_bound = center + radius;
          double h1 = cut_point / n1, h2 = (radius - cut_point) / n2;
          for (double p = center; p < center + cut_point; p += h1)
            interpolation_points.push_back(p);
          for (double p = center + cut_point; p < right_bound; p += h2)
            interpolation_points.push_back(p);
          interpolation_points.push_back(right_bound);

          // add enrichment function only when predicate radius is non-zero
          radial_problem.evaluate_at_x_values(interpolation_points,
                                              interpolation_values);


          // construct enrichment function and push
          Point<dim> p;
          prm.set_enrichment_point(p, i);
          SplineEnrichmentFunction<dim> func(p,
                                             interpolation_points,
                                             interpolation_values);
          vec_enrichments.push_back(
            std::make_shared<SplineEnrichmentFunction<dim>>(func));
        }
      else
        {
          pcout << "Dummy function added at " << i << std::endl;
          Functions::ConstantFunction<dim> func(0);
          vec_enrichments.push_back(
            std::make_shared<Functions::ConstantFunction<dim>>(func));
        }
    }
}



/*
 * Since each enrichment domain has different enrichment function
 * associated with it and the cells common to different enrichment
 * domains need to treated differently, we use helper function in
 * ColorEnriched namespace to construct finite element space and necessary
 * data structures required by it. The helper function is also used to
 * set FE indices of dof handler cells. We also set the cell with a unique
 * material id for now, which is used to map cells with a pairs of
 * color and corresponding enrichment function. All this is internally
 * used to figure out the correct set of enrichment functions to be used
 * for a cell.
 *
 * The quadrature point collection, with size equal to FE collection
 * is also constructed here.
 */
template <int dim>
void
LaplaceProblem<dim>::build_fe_space()
{
  pcout << "...building fe space" << std::endl;

  make_enrichment_functions();

  static std::unique_ptr<ColorEnriched::Helper<dim>> fe_space;
  fe_space = std::make_unique<ColorEnriched::Helper<dim>>(fe_base,
                                                          fe_enriched,
                                                          vec_predicates,
                                                          vec_enrichments);

  fe_collection = std::make_shared<const hp::FECollection<dim>>(
    fe_space->build_fe_collection(dof_handler));
  pcout << "size of fe collection: " << fe_collection->size() << std::endl;

  if (prm.debug_level == 9)
    {
      if (triangulation.n_active_cells() < 100)
        {
          pcout << "...start print fe indices" << std::endl;

          // print FE index
          const std::string base_filename =
            "fe_indices" + dealii::Utilities::int_to_string(dim) + "_p" +
            dealii::Utilities::int_to_string(0);
          const std::string filename = base_filename + ".gp";
          std::ofstream     f(filename);

          f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
            << std::endl
            << "set output \"" << base_filename << ".png\"" << std::endl
            << "set size square" << std::endl
            << "set view equal xy" << std::endl
            << "unset xtics" << std::endl
            << "unset ytics" << std::endl
            << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
            << std::endl;
          GridOut().write_gnuplot(triangulation, f);
          f << 'e' << std::endl;

          for (auto it : dof_handler.active_cell_iterators())
            f << it->center() << " \"" << it->active_fe_index() << "\"\n";

          f << std::flush << 'e' << std::endl;
          pcout << "...finished print fe indices" << std::endl;
        }

      if (triangulation.n_active_cells() < 100)
        {
          pcout << "...start print cell indices" << std::endl;

          // print cell ids
          const std::string base_filename =
            "cell_id" + dealii::Utilities::int_to_string(dim) + "_p" +
            dealii::Utilities::int_to_string(0);
          const std::string filename = base_filename + ".gp";
          std::ofstream     f(filename);

          f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
            << std::endl
            << "set output \"" << base_filename << ".png\"" << std::endl
            << "set size square" << std::endl
            << "set view equal xy" << std::endl
            << "unset xtics" << std::endl
            << "unset ytics" << std::endl
            << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
            << std::endl;
          GridOut().write_gnuplot(triangulation, f);
          f << 'e' << std::endl;

          for (auto it : dof_handler.active_cell_iterators())
            f << it->center() << " \"" << it->index() << "\"\n";

          f << std::flush << 'e' << std::endl;

          pcout << "...end print cell indices" << std::endl;
        }
    }

  // q collections the same size as different material identities
  q_collection.push_back(QGauss<dim>(4));
  for (unsigned int i = 1; i < fe_collection->size(); ++i)
    q_collection.push_back(QGauss<dim>(10));

  pcout << "...building fe space" << std::endl;
}



template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  pcout << "...start setup system" << std::endl;

  GridTools::partition_triangulation(n_mpi_processes,
                                     triangulation,
                                     SparsityTools::Partitioner::zoltan);

  dof_handler.distribute_dofs(*fe_collection);

  DoFRenumbering::subdomain_wise(dof_handler);
  std::vector<IndexSet> locally_owned_dofs_per_proc =
    DoFTools::locally_owned_dofs_per_subdomain(dof_handler);
  locally_owned_dofs    = locally_owned_dofs_per_proc[this_mpi_process];
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  SigmaFunction<dim> boundary_value_func;
  Point<dim>         p;
  prm.set_enrichment_point(p, 0);
  boundary_value_func.initialize(p, prm.sigmas[0], prm.boundary_value_expr);


  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           boundary_value_func,
                                           constraints);
  constraints.close();

  // Initialise the stiffness and mass matrices
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  std::vector<types::global_dof_index> n_locally_owned_dofs(n_mpi_processes);
  for (unsigned int i = 0; i < n_mpi_processes; ++i)
    n_locally_owned_dofs[i] = locally_owned_dofs_per_proc[i].n_elements();

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             n_locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);

  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  solution.reinit(locally_owned_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  pcout << "...finished setup system" << std::endl;
}


template <int dim>
void
LaplaceProblem<dim>::assemble_system()
{
  pcout << "...assemble system" << std::endl;

  system_matrix = 0;
  system_rhs    = 0;

  FullMatrix<double> cell_system_matrix;
  Vector<double>     cell_rhs;

  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<double> rhs_value, tmp_rhs_value;

  hp::FEValues<dim> fe_values_hp(*fe_collection,
                                 q_collection,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);


  for (auto &cell : dof_handler.active_cell_iterators())
    if (cell->subdomain_id() == this_mpi_process)
      {
        fe_values_hp.reinit(cell);
        const FEValues<dim> &fe_values = fe_values_hp.get_present_fe_values();

        const unsigned int &dofs_per_cell = cell->get_fe().dofs_per_cell;
        const unsigned int &n_q_points    = fe_values.n_quadrature_points;

        /*
         * Initialize rhs values vector to zero. Add values calculated
         * from each of different rhs functions (vec_rhs).
         */
        rhs_value.assign(n_q_points, 0);
        tmp_rhs_value.assign(n_q_points, 0);
        for (unsigned int i = 0; i < vec_rhs.size(); ++i)
          {
            vec_rhs[i].value_list(fe_values.get_quadrature_points(),
                                  tmp_rhs_value);

            // add tmp to the total one at quadrature points
            for (const auto q_point : fe_values.quadrature_point_indices())
              {
                rhs_value[q_point] += tmp_rhs_value[q_point];
              }
          }

        local_dof_indices.resize(dofs_per_cell);
        cell_system_matrix.reinit(dofs_per_cell, dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);

        cell_system_matrix = 0;
        cell_rhs           = 0;

        for (const auto q_point : fe_values.quadrature_point_indices())
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = i; j < dofs_per_cell; ++j)
                cell_system_matrix(i, j) +=
                  (fe_values.shape_grad(i, q_point) *
                   fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point));

              cell_rhs(i) +=
                (rhs_value[q_point] * fe_values.shape_value(i, q_point) *
                 fe_values.JxW(q_point));
            }

        // exploit symmetry
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i; j < dofs_per_cell; ++j)
            cell_system_matrix(j, i) = cell_system_matrix(i, j);

        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_system_matrix,
                                               cell_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
      }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  pcout << "...finished assemble system" << std::endl;
}

template <int dim>
unsigned int
LaplaceProblem<dim>::solve()
{
  pcout << "...solving" << std::endl;
  SolverControl solver_control(prm.max_iterations, prm.tolerance, false, false);
  PETScWrappers::SolverCG cg(solver_control);

  PETScWrappers::PreconditionSOR preconditioner(system_matrix);

  cg.solve(system_matrix, solution, system_rhs, preconditioner);

  Vector<double> local_soln(solution);

  constraints.distribute(local_soln);
  solution = local_soln;
  pcout << "...finished solving" << std::endl;
  return solver_control.last_step();
}



template <int dim>
void
LaplaceProblem<dim>::refine_grid()
{
  const Vector<double> localized_solution(solution);
  Vector<float>        local_error_per_cell(triangulation.n_active_cells());

  hp::QCollection<dim - 1> q_collection_face;
  for (unsigned int i = 0; i < q_collection.size(); ++i)
    q_collection_face.push_back(QGauss<dim - 1>(1));

  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    q_collection_face,
    std::map<types::boundary_id, const Function<dim> *>{},
    localized_solution,
    local_error_per_cell,
    ComponentMask(),
    nullptr,
    n_mpi_processes,
    this_mpi_process);
  const unsigned int n_local_cells =
    GridTools::count_cells_with_subdomain_association(triangulation,
                                                      this_mpi_process);
  PETScWrappers::MPI::Vector distributed_all_errors(
    mpi_communicator, triangulation.n_active_cells(), n_local_cells);
  for (unsigned int i = 0; i < local_error_per_cell.size(); ++i)
    if (local_error_per_cell(i) != 0)
      distributed_all_errors(i) = local_error_per_cell(i);
  distributed_all_errors.compress(VectorOperation::insert);
  const Vector<float> localized_all_errors(distributed_all_errors);
  GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                    localized_all_errors,
                                                    0.85,
                                                    0);
  triangulation.execute_coarsening_and_refinement();
  ++prm.global_refinement;
}



template <int dim>
void
LaplaceProblem<dim>::output_results(const unsigned int cycle)
{
  pcout << "...output results" << std::endl;
  pcout << "Patches used: " << prm.patches << std::endl;

  Vector<double>     exact_soln_vector, error_vector;
  SigmaFunction<dim> exact_solution;

  if (prm.exact_soln_expr != "")
    {
      // create exact solution vector
      exact_soln_vector.reinit(dof_handler.n_dofs());
      exact_solution.initialize(Point<dim>(),
                                prm.sigmas[0],
                                prm.exact_soln_expr);
      VectorTools::project(dof_handler,
                           constraints,
                           q_collection,
                           exact_solution,
                           exact_soln_vector);

      // create error vector
      error_vector.reinit(dof_handler.n_dofs());
      Vector<double> full_solution(localized_solution);
      error_vector += full_solution;
      error_vector -= exact_soln_vector;
    }

  Assert(cycle < 10, ExcNotImplemented());
  if (this_mpi_process == 0)
    {
      std::string filename = "solution-";
      filename += Utilities::to_string(cycle);
      filename += ".vtk";
      std::ofstream output(filename);

      DataOut<dim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(localized_solution, "solution");
      if (prm.exact_soln_expr != "")
        {
          data_out.add_data_vector(exact_soln_vector, "exact_solution");
          data_out.add_data_vector(error_vector, "error_vector");
        }
      data_out.build_patches(prm.patches);
      data_out.write_vtk(output);
      output.close();
    }
  pcout << "...finished output results" << std::endl;
}



// use this only when exact solution is known
template <int dim>
void
LaplaceProblem<dim>::process_solution()
{
  Vector<float> difference_per_cell(triangulation.n_active_cells());
  double        L2_error, H1_error;

  if (!prm.exact_soln_expr.empty())
    {
      pcout << "...using exact solution for error calculation" << std::endl;

      SigmaFunction<dim> exact_solution;
      exact_solution.initialize(Point<dim>(),
                                prm.sigmas[0],
                                prm.exact_soln_expr);

      VectorTools::integrate_difference(dof_handler,
                                        localized_solution,
                                        exact_solution,
                                        difference_per_cell,
                                        q_collection,
                                        VectorTools::L2_norm);
      L2_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::L2_norm);

      VectorTools::integrate_difference(dof_handler,
                                        localized_solution,
                                        exact_solution,
                                        difference_per_cell,
                                        q_collection,
                                        VectorTools::H1_norm);
      H1_error = VectorTools::compute_global_error(triangulation,
                                                   difference_per_cell,
                                                   VectorTools::H1_norm);
    }

  pcout << "refinement h_smallest Dofs L2_norm H1_norm" << std::endl;
  pcout << prm.global_refinement << ' '
        << prm.size / std::pow(2.0, prm.global_refinement) << ' '
        << dof_handler.n_dofs() << ' ' << L2_error << ' ' << H1_error
        << std::endl;
}



template <int dim>
void
LaplaceProblem<dim>::run()
{
  pcout << "...run problem" << std::endl;
  double norm_soln_old(0), norm_rel_change_old(1);

  // Run making grids and building FE space only once.
  initialize();
  build_fe_space();


  if (this_mpi_process == 0)
    deallog << "Solving problem with number of sources: " << prm.n_enrichments
            << std::endl;

  for (unsigned int cycle = 0; cycle <= prm.cycles; ++cycle)
    {
      pcout << "Cycle " << cycle << std::endl;

      setup_system();

      pcout << "Number of active cells:       "
            << triangulation.n_active_cells() << std::endl
            << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

      if (prm.debug_level == 9 && this_mpi_process == 0)
        plot_shape_function<dim>(dof_handler);

      assemble_system();
      auto n_iterations = solve();
      pcout << "Number of iterations: " << n_iterations << std::endl;
      localized_solution.reinit(dof_handler.n_dofs());
      localized_solution = solution;
      double value =
        VectorTools::point_value(dof_handler, localized_solution, Point<dim>());
      pcout << "Solution at origin:   " << value << std::endl;


      // calculate L2 norm of solution
      if (this_mpi_process == 0)
        {
          pcout << "calculating L2 norm of soln" << std::endl;
          double        norm_soln_new(0), norm_rel_change_new(0);
          Vector<float> difference_per_cell(triangulation.n_active_cells());
          VectorTools::integrate_difference(dof_handler,
                                            localized_solution,
                                            Functions::ZeroFunction<dim>(),
                                            difference_per_cell,
                                            q_collection,
                                            VectorTools::H1_norm);
          norm_soln_new =
            VectorTools::compute_global_error(triangulation,
                                              difference_per_cell,
                                              VectorTools::H1_norm);
          // relative change can only be calculated for cycle > 0
          if (cycle > 0)
            {
              norm_rel_change_new =
                std::abs((norm_soln_new - norm_soln_old) / norm_soln_old);
              pcout << "relative change of solution norm "
                    << norm_rel_change_new << std::endl;
            }

          // monitor relative change of norm in later stages
          if (cycle > 1)
            {
              deallog << (norm_rel_change_new < norm_rel_change_old)
                      << std::endl;
            }

          norm_soln_old = norm_soln_new;

          // first sample of relative change of norm comes only cycle = 1
          if (cycle > 0)
            norm_rel_change_old = norm_rel_change_new;

          pcout << "End of L2 calculation" << std::endl;
        }

      if (prm.debug_level >= 2 && this_mpi_process == 0)
        output_results(cycle);

      // Do not refine if loop is at the end
      if (cycle != prm.cycles)
        refine_grid();

      pcout << "...step run complete" << std::endl;
    }
  pcout << "...finished run problem" << std::endl;
}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem()
{
  dof_handler.clear();
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  /*
   * Case 1: single source with known solution
   */
  {
    ParameterCollection prm(
      2,     // dimension
      2,     // domain size
      1,     // cube shape
      3,     // global refinement
      2,     // num of cycles grid is refined and solved again
      1,     // fe base degree
      1,     // fe enriched degree
      50000, // max iterations
      1e-9,  // tolerance
      // rhs value
      "-(exp(-(x*x + y*y)/(2*sigma*sigma))*(- 2*sigma*sigma + x*x + y*y))/(2*sigma*sigma*sigma*sigma*sigma*sigma*pi)",
      // boundary value
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x + y*y)/(2*sigma*sigma))",
      // rhs value for radial problem solved to find enrichment function
      "-(exp(-(x*x)/(2*sigma*sigma))*(- 2*sigma*sigma + x*x))/(2*sigma*sigma*sigma*sigma*sigma*sigma*pi)",
      // boundary value for radial problem
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x)/(2*sigma*sigma))",
      // exact solution expression. If null nothing is done
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x + y*y)/(2*sigma*sigma))",
      1, // patches
      1, // debug level
      1, // num enrichments
      // enrichment points interpreted 2 at a time if dimension is 2
      {0, 0},
      // radii defining different predicates
      {0.4},
      // sigmas defining different predicates
      {0.1});

    LaplaceProblem<2> problem(prm);
    problem.run();
  }

  /*
   * Case 2: 3 sources
   */
  {
    ParameterCollection prm(
      2,     // dimension
      4,     // domain size
      1,     // cube shape
      3,     // global refinement
      1,     // num of cycles grid is refined and solved again
      1,     // fe base degree
      1,     // fe enriched degree
      50000, // max iterations
      1e-9,  // tolerance
      // rhs value
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x + y*y)/(2*sigma*sigma))",
      // boundary value
      "0",
      // rhs value for radial problem solved to find enrichment function
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x)/(2*sigma*sigma))",
      // boundary value for radial problem
      "0",
      // exact solution expression. If null nothing is done
      "",
      1, // patches
      1, // debug level
      3, // num enrichments
      // enrichment points interpreted 2 at a time if dimension is 2
      {0.5, 0.5, 0, 0, -1, -1},
      // radii defining different predicates
      {0.4, 0.4, 0.4},
      // sigmas defining different predicates
      {0.1, 0.1, 0.1});


    LaplaceProblem<2> problem(prm);
    problem.run();
  }

  /*
   * Case 3: five sources 3d
   */
  {
    ParameterCollection prm(
      3,     // dimension
      8,     // domain size
      1,     // cube shape
      4,     // global refinement
      0,     // num of cycles grid is refined and solved again
      1,     // fe base degree
      1,     // fe enriched degree
      50000, // max iterations
      1e-9,  // tolerance
      // rhs value
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x + y*y + z*z)/(2*sigma*sigma))",
      // boundary value
      "0",
      // rhs value for radial problem solved to find enrichment function
      "1.0/(2*pi*sigma*sigma)*exp(-(x*x)/(2*sigma*sigma))",
      // boundary value for radial problem
      "0",
      // exact solution expression. If null nothing is done
      "",
      1, // patches
      1, // debug level
      5, // num enrichments
      // enrichment points interpreted 3 at a time if dimension is 3
      {1.5, 1.5, 1.5, 1, 1, 1, 0, 0, 0, -1, -1, -1, 1, -1, -1},
      // radii defining different predicates
      {0.45, 0.45, 0.45, 0.45, 0.45},
      // sigmas defining different predicates
      {0.1, 0.1, 0.1, 0.1, 0.1});

    LaplaceProblem<3> problem(prm);
    problem.run();
  }
}
