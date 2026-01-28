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

// This is only necessary for the output of the debugging.
#include <iomanip>

// The last step is as in all previous programs: We put everything into a namespace
// that corresponds to the step number.
namespace Step101
{
  using namespace dealii;


  // @sect3{Parameter struct for CMS}

  // In this tutorial, we isolate input parameters in a struct so the main class can
  // store them and remain easily configurable. We start with the parameter struct for
  // CMS.
  struct CmsParameters
  {
    CmsParameters();
    // @sect4{Material parameters}

    // For isotropic linear elasticity, stress is related to strain via Lamé parameters,
    // which depend on the Young's modulus and the Poisson's ratio. Additionally, in dynamics,
    // the relevant quantity for explicit stability is the wave speed, which depends also
    // on the density of our isotropic medium.
    //
    // We set the following parameters all to 1.0 to keep this tutorial numerically simple. 
    // However, these parameters can be adapted according to your chosen material.
    // Changing these parameters affects the wave speed and thus, the stable explicit time
    // step, so the CMS acceptability limits set in this struct should be reconsidered
    // when switching to a real material.
    double lambda;
    double mu;
    double rho;

    // @sect4{Thickness}

    // We use the parameter thickness to compute the mass for the assessment of the
    // impact of the artificially added mass to the model. We still have a 2D problem,
    // this is a non-existing thickness only used for the assessment, therefore, we
    // set it to 1.0.
    double thickness;

    // @sect4{Selection of a bigger critical time step}

    // To select a bigger critical step, we define the parameter dt_factor.
    // The newly selected critical step, the target critical step,
    // <code>dt_target</code> will then be calculated through multiplying
    // the parameter <code>dt_factor</code> by the actual critical time step
    // of the model <code>dt_crit_global</code>.
    //
    // If the value of <code>dt_factor</code> is bigger than 1, we want to increase
    // the critical time step of the model, while a value of 1 would do nothing.
    double dt_factor;

    // @sect4{Acceptability criteria}

    // CMS changes the model's inertia and therefore, its physical response.
    // Therefore, CMS is only acceptable if the inertia forces that appear due to
    // the artificially added mass are of second order, and you have not scaled
    // the majority of your mesh elements. 
    //
    // How to assess the impact of CMS and how to set proper acceptability criteria
    // for your own FE model is not a triviality. Firstly, you have to consider
    // the mass of your own model. If your model already has a mass of 6 lb, adding
    // 3 lb will impact your model's inertia drastically. But what about adding
    // 0.5 lb or even only 0.1 lb? It can be tricky to know when you have an acceptable
    // mass ratio. Another important criterion for acceptability of CMS is how
    // much of your mesh needs scaling. If a large portion of the mesh needs it,
    // then your CMS is probably not appropriate, since CMS is meant to be only
    // locally implemented to a few tiny cells. But then again, how do you set an
    // appropriate value?
    //
    // Both questions about the exact values for these acceptability criteria do not
    // have a straightforward answer. CMS is not a technique that you can apply blindly.
    // Instead, you need to understand your model. How does your model exactly react?
    // Maybe your model has too tiny cells and you cannot compute this reaction.
    // Therefore, to use CMS you need to have some understanding of the reaction of
    // your model even without performing a simulation. You then can conduct CMS, and
    // it may look great or it may look completely different to what you expected.
    // This is where convergence studies come into play. If you have little experience
    // with CMS, you may need various iterations of your simulation, so that in each
    // iteration you vary your acceptability criteria. Then, you look at the physical
    // responses of all of your iterations, and you conduct a convergence analysis.
    // You will see that setting too high values for the acceptability criteria will
    // lead to adding a lot of mass, which drastically affects the response. But the
    // smaller the values for your acceptability criteria are, the response will be
    // altered less and less. There will then be a point in which you can confidently
    // say that adding a certain amount of artificial mass to certain elements of your
    // mesh barely affects your physical response. And if this point is computationally
    // effective, you may use those acceptability criteria. You may also need to look
    // at energies, and forces, it will all depend on your exact application of the model,
    // and what artificial impact you are willing to accept.
    //
    // Nevertheless, not all models are suitable for CMS. After conducting several
    // convergence studies, you may find out that you cannot reach an acceptable trade-off
    // between an efficient computation time and realistic results. Then, you may
    // need to remesh again your model with a coarser mesh or accept high computation
    // time. Yet, CMS can be a powerful tool for explicit FEM simulations, especially for models
    // that need fine meshes or that are particularly stiff, since their critical time
    // step will be strictly limited to a very small number.
    //
    // This tutorial introduces the use of the CMS technique in a FEM simulation,
    // and does not focus on solving a FEM simulation. Therefore, we use only two
    // acceptability criteria, <code>max_added_mass_fraction</code> and
    // <code>max_scaled_cell_fraction</code>, which we both set to 0.10, but you could
    // set the values as you wish.
    double max_added_mass_fraction;
    double max_scaled_cell_fraction;

    // We do not want to increase the critical time step of our model
    // at any cost. If increasing the critical time step to our selected time step
    // means that we will violate our acceptability criteria, we will reduce
    // the target critical time step until the acceptability criteria are not longer violated.
    bool auto_reduce_dt_target;

    // Therefore, if the acceptability criteria are violated, we will reduce
    // <code>dt_target</code> by the factor <code>dt_shrink_factor</code> and retry.
    // Here, you could choose your own value for the factor.
    double dt_shrink_factor;

    // In this tutorial we have also set a maximum number of iterations in the loop that
    // reduces the target critical step. You may need to remesh your model if you want
    // to achieve a specific time step and after the maximum number of iteratively
    // reducing the target critical step, the CMS still violates the acceptability criteria.
    unsigned int max_dt_adjustment_steps;

    // @sect4{Variable for output field}

    // Output a VTK field with the resulting density scaling factor.
    bool write_density_scaling_output;
  };



  // @sect4{Default constructor}

  // All parameters are initialized to reasonable tutorial values that illustrate the CMS
  // algorithm without introducing additional material complexity.
  inline CmsParameters::CmsParameters()
    : lambda(1.0)
    , mu(1.0)
    , rho(1.0)
    , thickness(1.0)
    , dt_factor(1.1)
    , max_added_mass_fraction(0.10)
    , max_scaled_cell_fraction(0.10)
    , auto_reduce_dt_target(true)
    , dt_shrink_factor(0.95)
    , max_dt_adjustment_steps(25)
    , write_density_scaling_output(true)
  {}



  // @sect3{The <code>ElasticProblem</code> class template}

  // This class is named after the standard deal.II tutorial naming. However, note again:
  // we do not solve the PDE system here.
  //
  // The intent is to keep the familiar FEM scaffolding, so later it can be extended to
  // an explicit code, where you would assemble or matrix-free apply stiffness and mass.
  template <int dim>
  class ElasticProblem
  {
  public:
     
    // We construct the problem using CMS parameters.
    ElasticProblem(const CmsParameters &parameters);
    void run();

  private:

  // @sect4{A result struct to return all CMS diagnostics}

  // We group computed outputs in a struct to avoid long argument lists and to keep run()
  // readable. We first set all of the outputs to 0.0, since the values will be updated
  // when running the code.
    struct CmsResult
    {
      CmsResult();

      // <code>dt_crit_global_before</code> is the global critical time step of the entire
      // model computed on the original density.
      double dt_crit_global_before;

      // <code>dt_target</code> is the critical time step we aim for after scaling, but it
      // may be reduced automatically.
      double dt_target;

      // <code>dt_crit_global_after</code> is the global critical time step of the entire
      // model after performing CMS.
      double dt_crit_global_after;

      // We need some parameters for mass bookkeeping. Keep in mind that we have a 2D problem
      // and no actual thickness, which means that we also have no mass. This is why we
      // already introduced a thickness of 1.0 to calculate this non-existent mass, which is
      // later only used to assess the acceptability of our application of CMS.
      double total_mass;
      double added_mass;
      double added_mass_ratio;

      // We need some parameters to keep track of how many cells were scaled and what fraction
      // this represents.
      unsigned int n_scaled;
      double scaled_fraction;

      // We also need the cached characteristic length per cell.
      std::vector<double> h_char;

      // <code>density_scaling</code> represents the value of the new density divided by
      // the old density per cell.
      std::vector<double> density_scaling;
    };

    // @sect4{Mesh and geometry}

    // The geometry of the model in this tutorial is fairly simple, it is just a rectangle.
    // If we meshed the rectangle in a "normal" way, we would get a mesh with elements with
    // the exact same sizes. Even if we refined the mesh, the mesh would only be finer, but
    // all of the elements would still have the exact same size. However, we use CMS when
    // we have some very small cells in comparison with other larger cells, since those few
    // very tiny cells are what are limiting our critical time step in a very strict way.
    //
    // Therefore, for this tutorial, we need a new function that lets us mesh our rectangular
    // domain with elements with various sizes without hanging nodes and with the exact same
    // <code>FE_Q(1)</code> elements used in step-8.
    void make_conforming_graded_mesh();

    // We also need a function that computes the minimum edge of every mesh element, also known
    // as cell, so that we can define their characteristic lengths. Since we are using
    // <code>FE_Q(1)</code> elements in a 2D problem, the characteristic length <code>h_char</code>
    // of a mesh element corresponds to its minimum edge length.
    //
    // If we wanted to use other types of mesh elements, we would need other formula to compute
    // the characteristic length. In the introduction of this tutorial there are several
    // references to FEM books which contain these formulas for various elements.
    static double compute_min_edge_length_of_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell
    );
    void compute_characteristic_length(std::vector<double> &h_char) const;

    // @sect4{Computation of critical time step}

    // We first compute the critical time step of every cell. Afterwards, we set the global critical
    // time step of the complete model to the minimum critical time step per cell.
    double compute_dt_cell(
      const double h, const double rho_local
    ) const;
    double compute_dt_crit_global(
      const std::vector<double> &h_char, const std::vector<double> &density_scaling
    ) const;

    // @sect4{CMS}

    // For a chosen <code>dt_target</code>, we scale the densities only where needed.
    // We also reduce <code>dt_target</code> where needed to satisfy the acceptability
    // criteria.
    CmsResult apply_cms_for_target_dt(
      const std::vector<double> &h_char, const double dt_target
    ) const;
    CmsResult choose_dt_and_apply_cms(
      const std::vector<double> &h_char
    ) const;

    // @sect4{Sanity checks and debug output}

    // Since we create the mesh in an unconventional way in this tutorial due to the
    // geometry of the model, we need to assert that the mesh is indeed a valid mesh
    // and that the domain area matches the expectations. If we were to apply CMS to
    // a model, which we mesh in a more conventional way, we would not need to assert this.
    void verify_mesh_is_quad_and_non_degenerate() const;

    // We print per-cell information.
    void print_cell_debug_table(
      const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target
    ) const;

    // @sect4{Output}

    // We write the per-cell density scaling factor to a VTK file.
    void output_density_scaling(
      const std::vector<double> &density_scaling
    ) const;

    // @sect4{Input parameters and FEM core objects}

    const CmsParameters parameters;

    Triangulation<dim> triangulation;
    DoFHandler<dim> dof_handler;

    // Even though we do not solve elasticity here, this mirrors real usage of CMS.
    const FESystem<dim> fe;
    };



    // @sect4{Constructors}
    template <int dim>
    inline ElasticProblem<dim>::CmsResult::CmsResult()
      : dt_crit_global_before(0.0)
      , dt_target(0.0)
      , dt_crit_global_after(0.0)
      , total_mass(0.0)
      , added_mass(0.0)
      , added_mass_ratio(0.0)
      , n_scaled(0)
      , scaled_fraction(0.0)
    {}



    template <int dim>
    ElasticProblem<dim>::ElasticProblem(const CmsParameters &parameters)
      : parameters(parameters)
      , dof_handler(triangulation)
      , fe(FE_Q<dim>(1) ^ dim)
    {}



    // @sect3{Mesh generation without hanging nodes}

    // This function would not be needed in a real application of CMS. We only need it
    // because we have a model with a simple rectangular geometry, and we need a mesh
    // that would mirror a mesh that would benefit from the application of CMS.
    //
    // In this mesh we need a few small cells that constrain the critical time step, but a globally conforming
    // mesh without adaptive refinement complications. We therefore build a structured
    // quad mesh on [-1,1]^2 with non-uniform coordinate spacing in x and in y. This
    // produces graded cell sizes and guarantees conformity without hanging modes and
    // T-junctions because we never refine adaptively.
    template <int dim>
    void ElasticProblem<dim>::make_conforming_graded_mesh()
    {
      Assert(dim == 2, ExcMessage("This demo is implemented for dim=2 only."));

      // Coordinates with clustering in several regions. This creates tiny cells that will
      // dominate the global critical time step of the complete model.
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

      // Vertex indexing helper for a tensor-product grid.
      auto vid = [&](unsigned int i, unsigned int j){
        return j * x.size() + i;
      };

      // We create all vertices.
      for (unsigned int j = 0; j < y.size(); ++j)
        for (unsigned int i = 0; i < x.size(); ++i)
          vertices.emplace_back(x[i],y[j]);

      // We create quad cells, i.e. <code>FE_Q(1)</code> elements by connecting vertex indices.
      std::vector<CellData<2>> cells;
      cells.reserve((x.size()-1) * (y.size()-1));

      for (unsigned int j = 0; j < y.size()-1; ++j)
        for (unsigned int i = 0; i < x.size()-1; ++i)
          {
            CellData<2> c;
            c.vertices[0] = vid(i,j);
            c.vertices[1] = vid(i+1,j);
            c.vertices[2] = vid(i,j+1);
            c.vertices[3] = vid(i+1,j+1);
            c.material_id = 0;
            cells.push_back(c);
          }

      SubCellData subcell_data;
      triangulation.create_triangulation(vertices,cells,subcell_data);
    }



  // @sect3{Characteristic element length}

  // We compute the minimum edge length of our mesh elements, since for <code>FE_Q(1)</code>
  // elements this corresponds to their characteristic element length.
  template <int dim>
  double ElasticProblem<dim>::compute_min_edge_length_of_cell(
    const typename Triangulation<dim>::active_cell_iterator &cell
  )
  {
    double min_length = std::numeric_limits<double>::max();

    // <code>GeometryInfo<dim>::lines_per_cell</code> gives the number of edges for the
    // reference cell/element. This for-loop finds out the minimum edge length of a
    // mesh element.
    for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell; ++line)
    {
      const unsigned int v0_index = GeometryInfo<dim>::line_to_cell_vertices(line,0);
      const unsigned int v1_index = GeometryInfo<dim>::line_to_cell_vertices(line,1);
      
      const Point<dim> &v0 = cell->vertex(v0_index);
      const Point<dim> &v1 = cell->vertex(v1_index);

      min_length = std::min(min_length, v0.distance(v1));
    }

    // A valid element must have positive edge lengths.
    Assert(min_length > 0.0, ExcInternalError());

    return min_length;
  }



  // We then set the characteristic element length <code>h_char</code> of each element
  // to be their corresponding minimum element edge length:
  // @f[
  // h_{\text{char}} = l_{\text{min}}
  // @f]
  template <int dim>
  void ElasticProblem<dim>::compute_characteristic_length(
    std::vector<double> &h_char
  ) const
  {
    h_char.clear();
    h_char.reserve(triangulation.n_active_cells());

    // We iterate over active cells via the DoFHandler to reflect typical FEM structure.
    for (const auto &cell : dof_handler.active_cell_iterators())
      h_char.push_back(compute_min_edge_length_of_cell(cell));

    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
  }



  // @sect3{Compute global critical time step}
  
  // We first compute the critical step time per cell with the formula:
  // @f[
  // \Delta t \le h_{\text{char}} \, \sqrt{\frac{\rho}{\lambda + 2 \, \mu}}
  // ]
  template <int dim>
  double ElasticProblem<dim>::compute_dt_cell(
    const double h,
    const double rho_local
  ) const
  {
    const double denom = parameters.lambda + 2.0 * parameters.mu;

    Assert(denom > 0.0, ExcMessage("(lambda + 2 * mu) must be positive."));
    Assert(rho_local > 0.0, ExcMessage("rho_local must be positive."));
    Assert(h > 0.0, ExcMessage("h must be positive."));

    return h * std::sqrt(rho_local / denom);
  }



  // We then set the global critical time step of the overall model.
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

  // CMS increases the density of the cell where needed. This increases inertia
  // locally. In a real transient problem, that changes acceleration and thus the
  // dynamic response. The acceptability criteria aim to keep this change small.
  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::apply_cms_for_target_dt(
    const std::vector<double> &h_char, const double dt_target
  ) const
  {
    Assert(dt_target > 0.0, ExcMessage("dt_target must be positive."));
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());

    const double rho0 = parameters.rho;
    const double denom = parameters.lambda + 2.0 * parameters.mu;

    Assert(rho0 > 0.0, ExcMessage("Material density must be positive."));
    Assert(denom > 0.0, ExcMessage("(lambda + 2*mu) must be positive to define a wave speed."));

    // The result object collects all diagnostics we want to report.
    CmsResult result;

    // Store inputs for later output/debugging.
    result.h_char = h_char;
    result.dt_target = dt_target;

    // Start with no scaling per-cell.
    result.density_scaling.assign(h_char.size(), 1.0);

    // The global critical time step before scaling will be updated by <code>choose_dt_and_apply_cms()</code>.
    result.dt_crit_global_before = 0.0;

    // In a real 3D , the mass of a cell would be the density multiplied by the volume.
    // Here we work in 2D, therefore, we interpret "mass" as the multiplication of density,
    // area and an artificial thickness parameter set to the number 1. This makes the "added
    // mass ratio" dimensionally consistent and easy to interpret, even though the model is 2D.
    double total_mass = 0.0;
    double added_mass = 0.0;

    unsigned int cell_index = 0;

    // We loop over active FE cells via DoFHandler.
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      AssertIndexRange(cell_index,h_char.size());

      // We add the mass of each cell to the parameter <code>total_mass</code>.
      const double cell_area = cell->measure();
      const double cell_volume = cell_area * parameters.thickness;

      total_mass += rho0 * cell_volume;

      // We compute the critical time step for each cell with the original density without scaling.
      const double dt_i = compute_dt_cell(h_char[cell_index],rho0);

      // If the critical time step is smaller than the desired critical time step
      // <code>dt_target</code>, we then scale the density of this cell, i.e. we add
      // artificial density to said cell.
      if (dt_i < dt_target)
      {

        // We compute the density needed so that the critical time step of this cell
        // corresponds to the desired critical time step of the model.
        const double h = h_char[cell_index];
        const double rho_required = denom * (dt_target * dt_target) / (h * h);

        // CMS only adds density, it never decreases density.
        if (rho_required > rho0)
        {
          const double alpha = rho_required / rho0;
          result.density_scaling[cell_index] = alpha;

          added_mass += (rho_required - rho0) * cell_volume;
        }
      }
      ++cell_index;
    }

    // We store mass diagnostics in the result object.
    result.total_mass = total_mass;
    result.added_mass = added_mass;
    result.added_mass_ratio = (total_mass > 0.0 ? added_mass / total_mass : 0.0);

    // We count how many cells were actually scaled.
    for (const double a : result.density_scaling)
      if (a > 1.0 + 1e-12)
        ++result.n_scaled;

    result.scaled_fraction = (
      result.density_scaling.empty() ? 0.0 : static_cast<double>(result.n_scaled) / 
      static_cast<double>(result.density_scaling.size()));

    // We finally compute the global critical time step of the whole model after scaling,
    // i.e. after performing CMS.
    result.dt_crit_global_after = compute_dt_crit_global(h_char, result.density_scaling);

    return result;
  }



  // @sect3{Choose the desired critical time step and ensure mass increase stays acceptable}

  // The function <code>choose_dt_and_apply_cms</code> implements the "outer loop" around
  // CMS. We compute the baseline global critical time step, we choose an initial target critical
  // time step and we apply CMS for that target and obtain diagnostics. If the acceptability
  // criteria are violated, we reduce <code>dt_target</code> and retry. This is important because
  // a too ambitious <code>dt_target</code> can require scaling many cells or adding too much mass.
  //
  // The acceptability checks used here are deliberately simple. <code>max_added_mass_fraction</code>
  // limits the total added mass relative to the original model mass, and <code>max_scaled_cell_fraction</code>
  // limits how much of the mesh is affected by scaling, because CMS should ideally remain a local
  // modification.
  //
  // In real applications, one may use additional criteria such as response metrics,
  // energy balance, or frequency shifts, depending on the exact application and goal of the simulation.
  template <int dim>
  typename ElasticProblem<dim>::CmsResult
  ElasticProblem<dim>::choose_dt_and_apply_cms(const std::vector<double> &h_char) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());

    // We compute the global critical time step before any scaling. We do this by calling
    // <code>compute_dt_crit_global()</code> with a scaling vector that is identically one.
    std::vector<double> unscaled_density(h_char.size(), 1.0);
    const double dt_global_before = compute_dt_crit_global(h_char, unscaled_density);

    // We choose the initial target critical time step.
    double dt_target = parameters.dt_factor * dt_global_before;

    // First attempt: We apply CMS for the requested target.
    CmsResult best = apply_cms_for_target_dt(h_char, dt_target);
    best.dt_crit_global_before = dt_global_before;

    // If the user disabled automatic reduction <code>auto_reduce_dt_target</code>, the code
    // returns the first attempt even if it violates the acceptability criteria. This is useful
    // for diagnostic runs where we explicitly want to see how intrusive a given
    // <code>dt_target</code> would be.
    if (!parameters.auto_reduce_dt_target)
      return best;

    // If we are allowed to reduce <code>dt_target</code>, we perform a simple backtracking loop.
    // We keep shrinking <code>dt_target</code> until the acceptability criteria are met or
    // after we reach the maximum number of iterations.
    Assert(parameters.dt_shrink_factor > 0.0 && parameters.dt_shrink_factor < 1.0, 
      ExcMessage("Shrink factor must be in (0,1) when we enable the automatic reduction of the target critical step."));

    for (unsigned int k = 0; k < parameters.max_dt_adjustment_steps; ++k)
    {
      const bool too_much_mass = (best.added_mass_ratio > parameters.max_added_mass_fraction);

      const bool too_many_cells = (best.scaled_fraction > parameters.max_scaled_cell_fraction);

      // If neither acceptability criterion is violated, we accept the current result.
      if(!too_much_mass && !too_many_cells)
        return best;

      // Otherwise, we reduce the target critical time step and retry.
      dt_target *= parameters.dt_shrink_factor;
      best = apply_cms_for_target_dt(h_char, dt_target);
      best.dt_crit_global_before = dt_global_before;
    }
    // We return the best effort after the maximum number of iterations. The user should then
    // critically assess the results to decide whether CMS is acceptable or not.
    return best;
  }



  // @sect3{Mesh verification}

  // We need this function because we have created a specific mesh in an unconventional way
  // to demonstrate CMS. In a real application of CMS, you would probably not need this function.
  template <int dim>
  void ElasticProblem<dim>::verify_mesh_is_quad_and_non_degenerate() const
  {
    Assert(dim == 2, ExcNotImplemented());

    // We check that each cell/element has 4 vertices.
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      Assert(cell->n_vertices() == 4, ExcInternalError());
      Assert(cell->measure() > 0.0, ExcMessage("Found a cell with non-positive area."));
    }

    // We check that the sum of all element areas matches the total area of the domain that was meshed,
    // which here is 4.0, since our model looks like this: [-1,1]x[-1,1].
    double area_sum = 0.0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      area_sum += cell->measure();

    const double expected_area = 4.0;
    std::cout << "Mesh check: sum(cell areas) = " << area_sum << " (expected ~ " << expected_area << ")\n";

    AssertThrow(std::abs(area_sum - expected_area) < 1e-12, 
    ExcMessage("Total mesh area does not match expected domain area."));
  }



  // @sect3{Per-cell debug output}

  // CMS decisions are local per cell. Printing the edge lengths, the characteristic element length,
  // the critical time step per cell, and whether the cell has been scaled and with which scaling factor,
  // lets us verify that CMS has been run correctly.
  // 
  // We keep this function in this tutorial to show how CMS works. In a real application
  // of CMS you would not print this information for each cell, since it would create
  // an excessively long output.
  template <int dim>
  void ElasticProblem<dim>::print_cell_debug_table
  (
    const std::vector<double> &h_char, const std::vector<double> &density_scaling, const double dt_target
  ) const
  {
    Assert(h_char.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(density_scaling.size() == triangulation.n_active_cells(), ExcInternalError());
    Assert(dt_target > 0.0, ExcMessage("The target critical time step must be positive."));

    const double denom = parameters.lambda + 2.0 * parameters.mu;
    const double rho0  = parameters.rho;

    // Save/restore stream state so we don't affect other output.
    const auto cout_flags = std::cout.flags();
    const auto cout_prec  = std::cout.precision();

    // Choose column widths once and use them consistently for header and rows.
    // If you ever add/remove a column, update these in one place.
    const int w_idx    = 4;
    const int w_area   = 8;
    const int w_edge   = 8;
    const int w_hchar  = 10;
    const int w_dt     = 10;
    const int w_scaled = 7;
    const int w_alpha  = 8;

    // 2 spaces between columns
    const std::string sep = "  ";

    std::cout << "\n--- Per-cell debug for a target critical step = " << dt_target << " ---\n";

    std::cout << std::left
              << std::setw(w_idx)    << "idx" << sep
              << std::setw(w_area)   << "area" << sep
              << std::setw(w_edge)   << "e0" << sep
              << std::setw(w_edge)   << "e1" << sep
              << std::setw(w_edge)   << "e2" << sep
              << std::setw(w_edge)   << "e3" << sep
              << std::setw(w_hchar)  << "h_char" << sep
              << std::setw(w_dt)     << "dt_cell" << sep
              << std::setw(w_scaled) << "scaled" << sep
              << std::setw(w_alpha)  << "factor" << sep
              << '\n';

    const int total_width =
      w_idx + w_area + 4*w_edge + w_hchar + w_dt + w_scaled + w_alpha + static_cast<int>(sep.size()) * 9;
    std::cout << std::string(total_width, '-') << '\n';

    unsigned int idx = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      double e[4];
      for (unsigned int l = 0; l < 4; ++l)
      {
        const unsigned int v0 = GeometryInfo<dim>::line_to_cell_vertices(l, 0);
        const unsigned int v1 = GeometryInfo<dim>::line_to_cell_vertices(l, 1);
        e[l] = cell->vertex(v0).distance(cell->vertex(v1));
      }

      const double h_char_i = *std::min_element(std::begin(e), std::end(e));
      const double dt_i     = h_char_i * std::sqrt(rho0 / denom);

      const double alpha  = density_scaling[idx];
      const bool   scaled = (alpha > 1.0 + 1e-12);

      // Row: numbers right-aligned, fixed precision for stable width.
      std::cout << std::right << std::fixed
                << std::setprecision(0) << std::setw(w_idx)  << idx << sep
                << std::setprecision(6) << std::setw(w_area) << cell->measure() << sep
                << std::setw(w_edge) << e[0] << sep
                << std::setw(w_edge) << e[1] << sep
                << std::setw(w_edge) << e[2] << sep
                << std::setw(w_edge) << e[3] << sep
                << std::setw(w_hchar) << h_char_i
                << std::setw(w_dt)    << dt_i << sep
                << std::left  << std::setw(w_scaled) << (scaled ? "YES" : "NO") << sep
                << std::right << std::setprecision(4) << std::setw(w_alpha) << alpha
                << '\n';

      AssertThrow(std::abs(h_char[idx] - h_char_i) < 1e-14,
                  ExcMessage("h_char mismatch: stored h_char != recomputed min edge length"));

      ++idx;
    }

    std::cout << "--- End per-cell debug ---\n";

    // Restore stream state.
    std::cout.flags(cout_flags);
    std::cout.precision(cout_prec);
  }



  // @sect3{Output density scaling to VTK}

  // We want to output the density scaling to a VTK file to see where CMS acts.
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

  // Lastly, this is the driver function. This function mirrors the structure
  // of many deal.II tutorial steps:
  //
  // - create mesh
  // - distribute DoFs
  // - compute something
  // - output results
  //
  // In an explicit FEM code, after the DoF distribution, you would:
  //
  // 1. Assemble stiffness K.
  // 2. Assemble mass M (lumped mass is common for explicit schemes).
  // 3. Integrate in time, e.g. using central differences, using a stable critical time step.
  //
  // CMS would enter between the steps 2 and 3. There, you would modify the mass
  // through adding artificial density to some cells to permit a larger
  // critical time step. This demonstration focuses exactly on the CMS step,
  // which can then later be implemented into an explicit FEM simulation code.
  template <int dim>
  void ElasticProblem<dim>::run()
  {
    Assert(dim == 2, ExcMessage("This step-101 is written for dim=2."));

    // As already mentioned before, we need to create an unconventional mesh
    // to apply CMS in a rectangular [-1,1]x[-1,1] domain. In a real application
    // of CMS we would substitute the following 2 lines with the meshing technique
    // or techniques suitable for our own specific application.
    make_conforming_graded_mesh();
    verify_mesh_is_quad_and_non_degenerate();

    // Even though we do not assemble anything, we distribute DoFs so the code
    // resembles a normal FEM application.
    dof_handler.distribute_dofs(fe);

    std::cout << "Step-101 (CMS) (CMS timestep targeting)\n";
    std::cout << "Active cells: " << triangulation.n_active_cells() << '\n';
    std::cout << "DoFs: " << dof_handler.n_dofs() << '\n';

    // We compute the characteristic element length.
    std::vector<double> h_char;
    compute_characteristic_length(h_char);

    // We choose a target critical time step, and apply CMS including optional
    // auto reduction of the target value.
    const auto cms = choose_dt_and_apply_cms(h_char);

    // We print cell-by-cell debug information to better show the functionality of CMS.
    print_cell_debug_table(h_char, cms.density_scaling, cms.dt_target);

    // We summarize time step targeting.
    std::cout << "\nSummary:\n";
    std::cout << "Target improvement factor = " << parameters.dt_factor << " ("<< (parameters.dt_factor - 1) * 100 <<"%)\n";
    std::cout << "Global critical time step before CMS = " << cms.dt_crit_global_before << '\n';
    std::cout << "Target critical time step = " << cms.dt_target << '\n';
    std::cout << "Global critical time step after CMS = " << cms.dt_crit_global_after << "\n";
    std::cout << "Real improvement factor = " << (cms.dt_crit_global_after / cms.dt_crit_global_before) 
              << " ("<< ((cms.dt_crit_global_after / cms.dt_crit_global_before) - 1) * 100 <<"%)\n";

    // We summarize CMS impact.
    std::cout << "\nCMS impact:\n";
    std::cout << "Added mass ratio: " << cms.added_mass_ratio * 100 << "%\n";
    std::cout << "Scaled cells: " << cms.n_scaled << " / " << cms.density_scaling.size() << " (" << 100.0 * cms.scaled_fraction << "%)\n";

    // We check the acceptability criteria.
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
      std::cout << " - If your application tolerates it, relax the acceptability limits, but only after checking whether the dynamic response of the model is still acceptable.\n";
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
// following is the main function of the program. Here we simply set parameters,
// construct the problem, and run it.
//
// You can vary the values of the parameters and observe how the results change.
int main()
{
  try
    {
      Step101::CmsParameters prm;

      // Basic isotropic parameters with demo values.
      prm.lambda = 1.0;
      prm.mu = 1.0;
      prm.rho = 1.0;

      prm.thickness = 1.0;

      // We ask for a modest increase in global stable dt.
      prm.dt_factor = 1.05;

      // We set the acceptability criteria.
      prm.max_added_mass_fraction = 0.10;
      prm.max_scaled_cell_fraction = 0.10;

      // We enable the auto-reduction of the critical time step if the acceptability
      // constraints are violated.
      prm.auto_reduce_dt_target = true;
      prm.dt_shrink_factor = 0.95;
      prm.max_dt_adjustment_steps = 25;

      prm.write_density_scaling_output = true;

      Step101::ElasticProblem<2> problem(prm);
      problem.run();
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