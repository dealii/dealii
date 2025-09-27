/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * This test is quasi identical to step-68, with the following exceptions:
 * - There is no load balancing
 * - A simplex mesh is used for the background mesh
 * - The Euler Analytical integration using the analytically defined velocity
 * field is not used because it would not test anything relevant
 * - Step parameters are hardcoded instead of drawn from a parameter file
 *
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Peter Munch
 */

// @sect3{Include files}

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>

#include <cmath>
#include <iostream>

#include "../tests.h"



namespace Step68
{



  template <int dim>
  class Vortex : public Function<dim>
  {
  public:
    Vortex()
      : Function<dim>(dim)
    {}


    virtual void
    vector_value(const Point<dim> &point,
                 Vector<double>   &values) const override;
  };


  template <int dim>
  void
  Vortex<dim>::vector_value(const Point<dim> &point,
                            Vector<double>   &values) const
  {
    const double T = 4;
    const double t = this->get_time();

    const double px = numbers::PI * point[0];
    const double py = numbers::PI * point[1];
    const double pt = numbers::PI / T * t;

    values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py);
    values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px);
    if (dim == 3)
      {
        values[2] = 0;
      }
  }



  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking();
    void
    run();

  private:
    void
    generate_particles();

    void
    setup_background_dofs();

    void
    interpolate_function_to_field();

    void
    euler_step_interpolated(const double dt);
    void
    euler_step_analytical(const double dt);

    void
    output_particles(const unsigned int it);
    void
    output_background(const unsigned int it);

    void
    log_particles();

    MPI_Comm                                       mpi_communicator;
    parallel::fullydistributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim>                particle_handler;

    DoFHandler<dim>                            fluid_dh;
    FESystem<dim>                              fluid_fe;
    MappingFE<dim>                             mapping;
    LinearAlgebra::distributed::Vector<double> velocity_field;

    Vortex<dim> velocity;

    // Simulation parameters. In step-68 drawn from prm file, here hardcoded.
    std::string output_directory = "./";

    static constexpr unsigned int velocity_degree  = 1;
    static constexpr double       time_step        = 0.004;
    static constexpr double       final_time       = 4.0;
    static constexpr unsigned int output_frequency = 1000;
    static constexpr unsigned int fluid_refinement = 8;
  };

  template <int dim>
  ParticleTracking<dim>::ParticleTracking()
    : mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(mpi_communicator)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_SimplexP<dim>(velocity_degree), dim)
    , mapping(FE_SimplexP<dim>(velocity_degree))
  {}

  // @sect4{Particles generation}

  // This function generates the tracer particles and the background
  // triangulation on which these particles evolve.
  template <int dim>
  void
  ParticleTracking<dim>::generate_particles()
  {
    // We create a hyper cube triangulation which we globally refine. This
    // triangulation covers the full trajectory of the particles.
    // parallel::distributed::Triangulation<dim> tria_pdt(mpi_communicator);
    //
    // GridGenerator::hyper_cube(tria_pdt, 0, 1);
    // tria_pdt.refine_global(fluid_refinement);

    Triangulation<dim> temporary_triangulation;
    GridGenerator::subdivided_hyper_cube_with_simplices(temporary_triangulation,
                                                        fluid_refinement);


    // extract relevant information from distributed triangulation
    auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(temporary_triangulation,
                                            mpi_communicator);
    background_triangulation.create_triangulation(construction_data);


    // This initializes the background triangulation where the particles are
    // living and the number of properties of the particles.
    particle_handler.initialize(background_triangulation, mapping, 1 + dim);

    // We create a particle triangulation which is solely used to generate
    // the points which will be used to insert the particles. This
    // triangulation is a hyper shell which is offset from the
    // center of the simulation domain. This will be used to generate a
    // disk filled with particles which will allow an easy monitoring
    // of the motion due to the vortex.
    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      center[2] = 0.5;

    const double outer_radius = 0.15;
    const double inner_radius = 0.01;

    Triangulation<dim> temporary_quad_particle_triangulation;

    GridGenerator::hyper_shell(temporary_quad_particle_triangulation,
                               center,
                               inner_radius,
                               outer_radius,
                               6);


    Triangulation<dim> temporary_tri_particle_triangulation;
    GridGenerator::convert_hypercube_to_simplex_mesh(
      temporary_quad_particle_triangulation,
      temporary_tri_particle_triangulation);
    temporary_tri_particle_triangulation.set_manifold(0, FlatManifold<dim>());


    // extract relevant information from distributed triangulation
    auto particle_construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(
        temporary_tri_particle_triangulation, mpi_communicator);



    parallel::fullydistributed::Triangulation<dim> particle_triangulation(
      mpi_communicator);
    particle_triangulation.create_triangulation(particle_construction_data);
    particle_triangulation.set_manifold(0, FlatManifold<dim>());

    // We generate the necessary bounding boxes for the particles generator.
    // These bounding boxes are required to quickly identify in which
    // process's subdomain the inserted particle lies, and which cell owns it.

    std::vector<BoundingBox<dim>> all_boxes;
    all_boxes.reserve(background_triangulation.n_locally_owned_active_cells());
    for (const auto &cell : background_triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        all_boxes.emplace_back(cell->bounding_box());
    const auto tree        = pack_rtree(all_boxes);
    const auto local_boxes = extract_rtree_level(tree, 1);

    std::vector<std::vector<BoundingBox<dim>>> global_bounding_boxes;
    global_bounding_boxes =
      Utilities::MPI::all_gather(mpi_communicator, local_boxes);

    // We generate an empty vector of properties. We will attribute the
    // properties to the particles once they are generated.
    std::vector<std::vector<double>> properties(
      particle_triangulation.n_locally_owned_active_cells(),
      std::vector<double>(dim + 1, 0.));

    // We generate the particles at the position of a single
    // point quadrature. Consequently, one particle will be generated
    // at the centroid of each cell.
    QGaussSimplex<dim> quadrature_formula(1);

    Particles::Generators::quadrature_points(particle_triangulation,
                                             quadrature_formula,
                                             global_bounding_boxes,
                                             particle_handler,
                                             mapping,
                                             properties);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      deallog << "Number of particles inserted: "
              << particle_handler.n_global_particles() << std::endl;
  }



  // @sect4{Background DOFs and interpolation}

  // This function sets up the background degrees of freedom used for the
  // velocity interpolation and allocates the field vector where the entire
  // solution of the velocity field is stored.
  template <int dim>
  void
  ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(fluid_dh);

    velocity_field.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }



  // This function takes care of interpolating the
  // vortex velocity field to the field vector. This is achieved rather easily
  // by using the VectorTools::interpolate() function.
  template <int dim>
  void
  ParticleTracking<dim>::interpolate_function_to_field()
  {
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
    velocity_field.update_ghost_values();
  }



  // @sect4{Time integration of the trajectories}

  // In contrast to the previous function in this function we
  // integrate the particle trajectories by interpolating the value of
  // the velocity field at the degrees of freedom to the position of
  // the particles.
  template <int dim>
  void
  ParticleTracking<dim>::euler_step_interpolated(const double dt)
  {
    Vector<double> local_dof_values(fluid_fe.dofs_per_cell);

    // We loop over all the local particles. Although this could be achieved
    // directly by looping over all the cells, this would force us
    // to loop over numerous cells which do not contain particles.
    // Rather, we loop over all the particles, but, we get the reference
    // of the cell in which the particle lies and then loop over all particles
    // within that cell. This enables us to gather the values of the velocity
    // out of the `velocity_field` vector once and use them for all particles
    // that lie within the cell.
    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto cell = particle->get_surrounding_cell();
        const auto dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);

        dh_cell->get_dof_values(velocity_field, local_dof_values);

        // Next, compute the velocity at the particle locations by evaluating
        // the finite element solution at the position of the particles.
        // This is essentially an optimized version of the particle advection
        // functionality in step 19, but instead of creating quadrature
        // objects and FEValues objects for each cell, we do the
        // evaluation by hand, which is somewhat more efficient and only
        // matters for this tutorial, because the particle work is the
        // dominant cost of the whole program.
        const auto pic = particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (auto &p : pic)
          {
            const Point<dim> reference_location = p.get_reference_location();
            Tensor<1, dim>   particle_velocity;
            for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
              {
                const auto comp_j = fluid_fe.system_to_component_index(j);

                particle_velocity[comp_j.first] +=
                  fluid_fe.shape_value(j, reference_location) *
                  local_dof_values[j];
              }

            Point<dim> &particle_location = particle->get_location();
            for (int d = 0; d < dim; ++d)
              particle_location[d] += particle_velocity[d] * dt;

            // Again, we store the particle velocity and the processor id in the
            // particle properties for visualization purposes.
            ArrayView<double> properties = p.get_properties();
            for (int d = 0; d < dim; ++d)
              properties[d] = particle_velocity[d];

            properties[dim] =
              Utilities::MPI::this_mpi_process(mpi_communicator);

            ++particle;
          }
      }
  }



  // @sect4{Data output}

  // The next two functions take care of writing both the particles
  // and the background mesh to vtu with a pvtu record. This ensures
  // that the simulation results can be visualized when the simulation is
  // launched in parallel.
  template <int dim>
  void
  ParticleTracking<dim>::output_particles(const unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("process_id");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    particle_output.build_patches(particle_handler,
                                  solution_names,
                                  data_component_interpretation);
    const std::string output_folder(output_directory);
    const std::string file_name("interpolated-particles");

    particle_output.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void
  ParticleTracking<dim>::log_particles()
  {
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      deallog << "Particles location" << std::endl;
    for (unsigned int proc = 0;
         proc < Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++proc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == proc)
          {
            for (auto part : particle_handler)
              {
                deallog << part.get_location() << std::endl;
              }
          }
        MPI_Barrier(mpi_communicator);
      }
  }



  template <int dim>
  void
  ParticleTracking<dim>::output_background(const unsigned int it)
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    // Attach the solution data to data_out object
    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(velocity_field,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(background_triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = background_triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping);

    const std::string output_folder(output_directory);
    const std::string file_name("background");

    data_out.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }


  template <int dim>
  void
  ParticleTracking<dim>::run()
  {
    DiscreteTime discrete_time(0, final_time, time_step);

    generate_particles();

    // We set the initial property of the particles by doing an
    // explicit Euler iteration with a time-step of 0 both in the case
    // of the analytical and the interpolated approach.
    setup_background_dofs();
    interpolate_function_to_field();
    euler_step_interpolated(0.);

    output_particles(discrete_time.get_step_number());
    output_background(discrete_time.get_step_number());

    // The particles are advected by looping over time.
    while (!discrete_time.is_at_end())
      {
        discrete_time.advance_time();
        velocity.set_time(discrete_time.get_previous_time());

        interpolate_function_to_field();
        euler_step_interpolated(discrete_time.get_previous_step_size());

        // After the particles have been moved, it is necessary to identify
        // in which cell they now reside. This is achieved by calling
        // <code>sort_particles_into_subdomains_and_cells</code>
        particle_handler.sort_particles_into_subdomains_and_cells();

        if ((discrete_time.get_step_number() % output_frequency) == 0)
          {
            output_particles(discrete_time.get_step_number());
            output_background(discrete_time.get_step_number());
          }
      }
    log_particles();
  }

} // namespace Step68



// @sect3{The main() function}

// The remainder of the code, the `main()` function, is standard.
// We note that we run the particle tracking with the analytical velocity
// and the interpolated velocity and produce both results
int
main(int argc, char *argv[])
{
  using namespace Step68;
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  deallog.depth_console(1);


  try
    {
      Step68::ParticleTracking<2> particle_tracking;
      particle_tracking.run();
    }
  catch (const std::exception &exc)
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
