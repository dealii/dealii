/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2000 - 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */


// @sect3{Include files}

// First we include header files that we have already used in previous
// example programs:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

// This again is C++:
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

// The last step is as in all previous programs:
namespace Step101
{
  using namespace dealii;


  // @sect3{Parameter struct for CMS}

  // We start with the parameter struct for CMS. We need to have it at the beginning
  // of the code because we call it in other parts of the code and else, we would
  // get an error message.
  struct CmsParameters
  {
    // Material parameters (Lamé parameters lambda and mu + density)
    double lambda = 1.0;
    double mu = 1.0;
    double rho = 1.0;

    // We use the parameter thickness to compute the mass for the assessment of the
    // impact of the artificially added mass to the model. We still have a 2D problem,
    // this is a non-existing thickness only used for the assessment, therefore, we
    // set it to 1.0.
    double thickness = 1.0; 

    // To select a bigger critical step, we define the parameter dt_factor.
    // The newly selected critical step, the target critical step,
    // <code>dt_target</code> will then be calculated through
    // <code>dt_target = dt_factor * dt_crit_global</code>, which is based on the selected
    // critical cells.
    double dt_factor = 1.1;

    // Acceptability criteria

    // Safety cap used in the assessment of the impact of the artificially added
    // mass to the system's response. If the mass increase of the model is greater
    // than <code>max_added_mass_fraction</code>, the target critical step
    // <code>dt_target</code> will be reduced iteratively until the mass of the model
    // equals or is smaller than <code>max_added_mass_fraction</code> or until
    // the maximum number of iterations is reached.
    double max_added_mass_fraction = 0.10;
    double max_scaled_cell_fraction = 0.10;

    bool auto_reduce_dt_target = true;

    // When the mass increase is too large, shrink <code>dt_target</code> by
    // the factor <code>dt_shrink_factor</code> and retry.
    double dt_shrink_factor = 0.95;

    // Maximum number of iterations in the loop that reduces the target critical
    // step so that the artificially added mass to the system does not have
    // a major impact on the system's response. If this number of iterations
    // is reached and the artificially added mass is bigger than
    // <code>max_added_mass_fraction</code> multiplied by the original mass,
    // there are two possible options. The first option is that you increase
    // the factor <code>max_added_mass_fraction</code>. However, you should do
    // this carefully considering your specific system. Will an increase of this
    // factor result in too big inertia forces? Then you should not do this.
    // The second option is that CMS is not suitable for your model. You may consider
    // remeshing it with a coarser mesh.
    unsigned int max_dt_adjustment_steps = 25;

    // Output a VTK field with the resulting density scaling factor.
    bool write_density_scaling_output = true;
  };

  // @sect3{The <code>ElasticProblem</code> class template}

  // The first part of the program is the declaration of the main class. All of
  // this has been copied verbatim from step-8.

  // TODO: If something is changed/added compared to step-8, comment it
  template <int dim>
  // TODO: Check if this class is needed
   class ElasticProblem
   {
   public:
     //ElasticProblem();
     ElasticProblem(const CmsParameters &parameters); // TODO: Check if this is correct
     void run();

   private:
   struct CmsResult
     {
      double dt_crit_global_before = 0.0;
      double dt_target = 0.0;
      double dt_crit_global_after = 0.0;

      double total_mass = 0.0;
      double added_mass = 0.0;
      double added_mass_ratio = 0.0;

      unsigned int n_scaled = 0;
      double scaled_fraction = 0.0;

      std::vector<double> h_char;

      // Per-cell density scaling alpha = rho_new / rho_old
      std::vector<double> density_scaling;
     };

     // Mesh and geometry
     void make_conforming_graded_mesh();

     static double compute_min_edge_length_of_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell
     );

     void compute_characteristic_length(std::vector<double> &h_char) const;

     // CFL / dt computations
     double compute_dt_cell(const double h, const double rho_local) const;
     double compute_dt_crit_global(const std::vector<double> &h_char, const std::vector<double> &density_scaling) const;

     // CMS
     CmsResult apply_cms_for_target_dt(const std::vector<double> &h_char, const double dt_target) const;

     CmsResult choose_dt_and_apply_cms(const std::vector<double> &h_char) const;

     void verify_mesh_is_quad_and_non_degenerate() const;

     void print_cell_debug_table(const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target) const;

     // Output
     void output_density_scaling(const std::vector<double> &density_scaling) const;

     const CmsParameters parameters;

     Triangulation<dim> triangulation;
     DoFHandler<dim> dof_handler;
     const FESystem<dim> fe;

    };

    template <int dim>
    ElasticProblem<dim>::ElasticProblem(const CmsParameters &parameters)
      : parameters(parameters)
      , dof_handler(triangulation)
      , fe(FE_Q<dim>(1) ^ dim)
    {}

    // Mesh generation without hanging nodes
    // Build a structured quad mesh on [-1,1]² but non-uniform spacing in x/y (clustering).
    // This yields different cell sizes AND is conforming (no hanging modes / no T-junctions).

    template <int dim>
    void ElasticProblem<dim>::make_conforming_graded_mesh()
    {
      Assert(dim == 2, ExcMessage("This demo is implemented for dim=2 only."));

      // Coordinates with clustering in several regions:
      //  - coarse on the left/bottom
      //  - finer near x ~ 0.6 and y ~ 0.6
      //  - also some refinement near x ~ -0.4, y ~ 0.2 (second region)
      const std::vector<double> x = {
        -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4,
        0.52, 0.58, 0.62, 0.8, 0.9, 1.0
      };

      const std::vector<double> y = {
        -1.0, -0.7, -0.4, -0.2, 0.0, 0.15, 0.22, 0.3,
        0.45, 0.55, 0.6, 0.63, 0.7, 0.85, 1.0
      };

      std::vector<Point<2>> vertices;
      vertices.reserve(x.size() * y.size());

      auto vid = [&](unsigned int i, unsigned int j){
        return j * x.size() + i;
      };

      for (unsigned int j = 0; j < y.size(); ++j)
        for (unsigned int i = 0; i < x.size(); ++i)
          vertices.emplace_back(x[i],y[j]);

      std::vector<CellData<2>> cells;
      cells.reserve((x.size()-1) * (y.size()-1));

      for (unsigned int j = 0; j < y.size()-1; ++j)
        for (unsigned int i = 0; i < x.size()-1; ++i)
          {
            CellData<2> c;
            c.vertices[0] = vid(i,j); // bottom-left
            c.vertices[1] = vid(i+1,j); // bottom-right
            c.vertices[2] = vid(i,j+1); // top-left
            c.vertices[3] = vid(i+1,j+1); //top-right
            c.material_id = 0;
            cells.push_back(c);
          }

      SubCellData subcell_data;
      triangulation.create_triangulation(vertices,cells,subcell_data);

      // No global / local refine here -> mesh stays conforming.
      // No hanging nodes can appear because we never do adaptive refinement

    }


  // @sect3{The <code>ElasticProblem</code> class implementation}

  // @sect4{Characteristic element length: minimum edge length}

  template <int dim>
  double ElasticProblem<dim>::compute_min_edge_length_of_cell(
    //const typename DoFHandler<dim>::active_cell_iterator &cell
    const typename Triangulation<dim>::active_cell_iterator &cell
  )
  {
    double min_length = std::numeric_limits<double>::max();

    for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell; ++line)
    {
      const unsigned int v0_index = GeometryInfo<dim>::line_to_cell_vertices(line,0);
      const unsigned int v1_index = GeometryInfo<dim>::line_to_cell_vertices(line,1);

      const Point<dim> &v0 = cell->vertex(v0_index);
      const Point<dim> &v1 = cell->vertex(v1_index);

      min_length = std::min(min_length, v0.distance(v1));
    }

    Assert(min_length > 0.0, ExcInternalError());
    return min_length;
  }

  template <int dim>
  void ElasticProblem<dim>::compute_characteristic_length(
    std::vector<double> &h_char
  ) const
  {
    h_char.clear();
    h_char.reserve(triangulation.n_active_cells());

    for (const auto &cell : dof_handler.active_cell_iterators())
      h_char.push_back(compute_min_edge_length_of_cell(cell));

    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
  }


  // @sect{Compute global critical time step}
  template <int dim>
  double ElasticProblem<dim>::compute_dt_cell(
    const double h,
    const double rho_local

  ) const
  {

    const double denom = parameters.lambda + 2.0 * parameters.mu;
    Assert(denom > 0.0, ExcMessage("lambda + 2 must be positive."));
    Assert(rho_local > 0.0, ExcMessage("rho_local must be positive."));
    Assert(h > 0.0, ExcMessage("h must be positive."));

    return h * std::sqrt(rho_local / denom);
  }

  template <int dim>
  double ElasticProblem<dim>::compute_dt_crit_global(
    const std::vector<double> &h_char,
    const std::vector<double> &density_scaling
  ) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    double dt_min = std::numeric_limits<double>::max();

    for (unsigned int i = 0; i < h_char.size(); ++i)
    {
      const double rho_i = parameters.rho * density_scaling[i];
      const double dt_i = compute_dt_cell(h_char[i],rho_i);
      dt_min = std::min(dt_min, dt_i);
    }

    return dt_min;
  }

  // @sect4{Apply classical mass scaling (CMS)}

  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::apply_cms_for_target_dt(
    const std::vector<double> &h_char, const double dt_target
  ) const
  {
    Assert(dt_target > 0.0, ExcMessage("dt_target must be positive."));
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());

    CmsResult result;
    result.h_char = h_char;
    result.dt_target = dt_target;

    const double rho0 = parameters.rho;
    const double denom = parameters.lambda + 2.0 * parameters.mu;

    result.density_scaling.assign(h_char.size(), 1.0);

    // First compute dt_crit_global_before (all scaling=1)
    std::vector<double> ones(h_char.size(),1.0);
    result.dt_crit_global_before = compute_dt_crit_global(h_char, ones);

    // Mass accounting
    double total_mass = 0.0;
    double added_mass = 0.0;

    unsigned int cell_index = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const double cell_volume = cell->measure() * parameters.thickness;
      total_mass += rho0 * cell_volume;

      // Per-cell dt before scaling
      const double dt_i = compute_dt_cell(h_char[cell_index],rho0);

      // Scale exactly those cells that violate dt_target:
      if (dt_i < dt_target)
      {
        // Want: dt_target = h * sqrt(rho_required / denom)
        // => rho_required = denom * (dt_target / h)^2
        const double h = h_char[cell_index];
        const double rho_required = denom * (dt_target * dt_target) / (h * h);

        if (rho_required > rho0)
        {
          const double alpha = rho_required / rho0;
          result.density_scaling[cell_index] = alpha;
          added_mass += (rho_required - rho0) * cell_volume;
        }
      }

      ++cell_index;
    }

    result.total_mass = total_mass;
    result.added_mass = added_mass;
    result.added_mass_ratio = (total_mass > 0.0 ? added_mass / total_mass : 0.0);

    // Count scaled cells
    for (const double a : result.density_scaling)
      if (a > 1.0 + 1e-12)
        ++result.n_scaled;

    result.scaled_fraction = (
      result.density_scaling.empty() ? 0.0 : static_cast<double>(result.n_scaled) / static_cast<double>(result.density_scaling.size()));

    // dt_crit_global_after
    result.dt_crit_global_after = compute_dt_crit_global(h_char, result.density_scaling);

    return result;

  }


  //@sect4{Choose dt_target and ensure mass increase stays acceptable}

  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::choose_dt_and_apply_cms(const std::vector<double> &h_char) const
  {
    // global dt before CMS
    std::vector<double> ones(h_char.size(), 1.0);
    const double dt_global_before = compute_dt_crit_global(h_char, ones);

    double dt_target = parameters.dt_factor * dt_global_before;

    CmsResult best = apply_cms_for_target_dt(h_char, dt_target);

    if (!parameters.auto_reduce_dt_target)
      return best;

    // If constraints violated, reduce dt_target iteratively.
    for (unsigned int k = 0; k < parameters.max_dt_adjustment_steps; ++k)
    {
      const bool too_much_mass = (best.added_mass_ratio > parameters.max_added_mass_fraction);

      const bool too_many_cells = (best.scaled_fraction > parameters.max_scaled_cell_fraction);

      if(!too_much_mass && !too_many_cells)
        return best;

      dt_target *= parameters.dt_shrink_factor;
      best = apply_cms_for_target_dt(h_char, dt_target);
    }

    // Return best effort after max iterations
    return best;
  }

  // @sect4{}

  template <int dim>
  void ElasticProblem<dim>::verify_mesh_is_quad_and_non_degenerate() const
  {
    Assert(dim == 2, ExcNotImplemented());

    // 1) Check number of vertices per cell (quad should have 4)
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      Assert(cell->n_vertices() == 4, ExcInternalError());
      Assert(cell->measure() > 0.0, ExcMessage("Found a cell with non-positive area."));
    }

    // 2) Check total area
    double area_sum = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators())
      area_sum += cell->measure();

    const double expected_area = 4.0; // [-1,1]x[-1,1]
    std::cout << "Mesh check: sum(cell areas) = " << area_sum << " (expected ~ " << expected_area << ")\n";

    AssertThrow(std::abs(area_sum - expected_area) < 1e-12, ExcMessage("Total mesh area does not match expected domain area."));
  }

  // @sect4{}
  template <int dim>
  void ElasticProblem<dim>::print_cell_debug_table
  (
    const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target
  ) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    const double denom = parameters.lambda + 2.0 * parameters.mu;
    const double rho0 = parameters.rho;

    std::cout << "\n--- Per-cell debug (first all cells) ---\n";
    std::cout << "idx   area     e0       e1       e2     e3      h_min      dt_cell      alpha       scaled\n";

    unsigned int idx = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
    {
      // edge lengths for quad
      double e[4];
      for (unsigned int l = 0; l < 4; ++l)
      {
        const unsigned int v0 = GeometryInfo<dim>::line_to_cell_vertices(l, 0);
        const unsigned int v1 = GeometryInfo<dim>::line_to_cell_vertices(l, 1);
        e[l] = cell->vertex(v0).distance(cell->vertex(v1));
      }

      const double hmin = *std::min_element(std::begin(e), std::end(e));
      const double dt_i = hmin * std::sqrt(rho0 / denom);

      const double alpha = density_scaling[idx];
      const bool scaled = (alpha > 1.0 + 1e-12);

      std::cout << idx << "   "
                << cell->measure() << "    "
                << e[0] << "    " << e[1] << "    " << e[2] << "    " << e[3] << "   "
                << hmin << "    "
                << dt_i << "    "
                << alpha << "    "
                << (scaled ? "YES" : "NO")
                << "\n";

      // check h_char matches what we just computed
      AssertThrow(std::abs(h_char[idx] - hmin) < 1e-14, ExcMessage("h_char mismatch: stored h_char != recomputed min edge length"));

      ++idx;
    }

    std::cout << "--- End per-cell debug ---\n";
  }

  // @sect4{Output density scaling to VTK}

  template <int dim>
  void ElasticProblem<dim>::output_density_scaling(
    const std::vector<double> &density_scaling
  ) const
  {
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());

    Vector<float> scaling_field(triangulation.n_active_cells());

    for (unsigned int i = 0; i < scaling_field.size(); ++i)
      scaling_field[i] = static_cast<float>(density_scaling[i]);

    DataOut<dim> data_out;
    data_out.attach_triangulation(triangulation);
    data_out.add_data_vector(scaling_field, "rho_scaling");
    data_out.build_patches();

    std::ofstream output("cms_rho_scaling.vtk");
    data_out.write_vtk(output);
  }

  // @sect4{ElasticProblem::run}

  // The <code>run</code> function does the same things as in step-6, for
  // example. This time, we use the square [-1,1]^d as domain, and we refine
  // it globally four times before starting the first iteration.
  //
  // The reason for refining is a bit accidental: we use the QGauss
  // quadrature formula with two points in each direction for integration of the
  // right hand side; that means that there are four quadrature points on each
  // cell (in 2d). If we only refine the initial grid once globally, then there
  // will be only four quadrature points in each direction on the
  // domain. However, the right hand side function was chosen to be rather
  // localized and in that case, by pure chance, it happens that all quadrature
  // points lie at points where the right hand side function is zero (in
  // mathematical terms, the quadrature points happen to be at points outside
  // the <i>support</i> of the right hand side function). The right hand side
  // vector computed with quadrature will then contain only zeroes (even though
  // it would of course be nonzero if we had computed the right hand side vector
  // exactly using the integral) and the solution of the system of
  // equations is the zero vector, i.e., a finite element function that is zero
  // everywhere. In a sense, we
  // should not be surprised that this is happening since we have chosen
  // an initial grid that is totally unsuitable for the problem at hand.
  //
  // The unfortunate thing is that if the discrete solution is constant, then
  // the error indicators computed by the KellyErrorEstimator class are zero
  // for each cell as well, and the call to
  // Triangulation::refine_and_coarsen_fixed_number() will not flag any cells
  // for refinement (why should it if the indicated error is zero for each
  // cell?). The grid in the next iteration will therefore consist of four
  // cells only as well, and the same problem occurs again.
  //
  // The conclusion needs to be: while of course we will not choose the
  // initial grid to be well-suited for the accurate solution of the problem,
  // we must at least choose it such that it has the chance to capture the
  // important features of the solution. In this case, it needs to be able to
  // see the right hand side. Thus, we refine globally four times. (Any larger
  // number of global refinement steps would of course also work.)
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    Assert(dim == 2, ExcMessage("This step-101 is written for dim=2."));

    make_conforming_graded_mesh();

    dof_handler.distribute_dofs(fe);

    verify_mesh_is_quad_and_non_degenerate();

    std::cout << "Step-101 (CMS) (CMS timestep targetting)\n";
    std::cout << "Active cells: " << triangulation.n_active_cells() << '\n';
    std::cout << "DoFs: " << dof_handler.n_dofs() << '\n';

    // assemble_stiffness_and_rhs();
    // assemble_consistent_mass_matrix();

    // Characteristic sizes
    std::vector<double> h_char;
    compute_characteristic_length(h_char);

    const auto cms = choose_dt_and_apply_cms(h_char);

    print_cell_debug_table(h_char, cms.density_scaling, cms.dt_target);

    double total_area = 0.0;
    for (const auto &cell : triangulation.active_cell_iterators())
      total_area += cell->measure();

    std::cout << "\nTotal model area = " << total_area << "\n";

    std::cout << "\nTargetting:\n";
    std::cout << "dt_factor = " << parameters.dt_factor << '\n';
    std::cout << "dt_crit_global (before) = " << cms.dt_crit_global_before << '\n';
    std::cout << "dt_target = " << cms.dt_target << '\n';
    std::cout << "dt_crit_global (after) = " << cms.dt_crit_global_after << "\n";
    std::cout << "improvement factor = " << (cms.dt_crit_global_after / cms.dt_crit_global_before) << "\n";

    std::cout << "\nCMS impact:\n";
    std::cout << "Scaled cells: " << cms.n_scaled << " / " << cms.density_scaling.size() << " (" << 100.0 * cms.scaled_fraction << "%)\n";

    const bool too_many_cells = (cms.scaled_fraction > parameters.max_scaled_cell_fraction);

    const bool too_much_mass = (cms.added_mass_ratio > parameters.max_added_mass_fraction);

    if (too_many_cells || too_much_mass)
    {
      std::cout << "\nWARNING: CMS may not be reasonable for this setup.\n";
      
      if (too_many_cells)
      {
        std::cout << " - Scaled-cell fraction exceeds limit: " << 100.0 * cms.scaled_fraction << " %\n";
      }

      if (too_much_mass)
      {
        std::cout << " - Added-mass ratio exceeds limit: " << 100.0 * cms.added_mass_ratio << " % > " << 100.0 * parameters.max_added_mass_fraction << " %\n";
      }

      std::cout << "Recommendations:\n";
      std::cout << " - Reduce dt_factor (smaller timestep increase), or \n";
      std::cout << " - Remesh with fewer extremely small cells / smoother grading.\n";
    } 
    
    if (parameters.write_density_scaling_output)
    {
      output_density_scaling(cms.density_scaling);
      std::cout << "\nWrote: cms_rho_scaling.vtk (cell data: rho_scaling)\n";
    }
  }

} // namespace Step101

// @sect3{The <code>main</code> function}

// After closing the <code>Step101</code> namespace in the last line above, the
// following is the main function of the program and is again exactly like in
// step-6 (apart from the changed class names, of course).
int main()
{
  try
    {
      //Step101::ElasticProblem<2> elastic_problem_2d;
      Step101::CmsParameters prm;

      prm.lambda = 1.0;
      prm.mu = 1.0;
      prm.rho = 1.0;

      prm.thickness = 1.0;

      prm.dt_factor = 1.05;

      prm.max_added_mass_fraction = 0.10;
      prm.max_scaled_cell_fraction = 0.10;

      prm.auto_reduce_dt_target = true;
      prm.dt_shrink_factor = 0.95;
      prm.max_dt_adjustment_steps = 25;

      prm.write_density_scaling_output = true;

      Step101::ElasticProblem<2> problem(prm);
      problem.run();
      //elastic_problem_2d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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