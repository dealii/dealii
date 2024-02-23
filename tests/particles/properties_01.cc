// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// When inserting multiple particles at once, we forgot to set their
// property pool pointer.

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"


template <unsigned int dim>
class SolutionFunction : public Function<dim>
{
public:
  inline SolutionFunction()
    : Function<dim>(dim + 1)
  {}

  inline virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    for (unsigned int i = 0; i < dim; ++i)
      {
        values[i] = p[i];
      }
    values[dim] = 1.0;
  }
};


template <int dim>
void
test()
{
  // number of "unknowns"
  const unsigned int ncomponents = dim + 1;

  // number of particles to randomly place
  const unsigned int nparticles = 10;

  // create the mesh
  const MappingQ1<dim> mapping;
  Triangulation<dim>   triangulation;
  DoFHandler<dim>      dofHandler(triangulation);
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  triangulation.refine_global(5);
  const dealii::FESystem<dim> fe(dealii::FE_Q<dim>(1), ncomponents);
  dofHandler.distribute_dofs(fe);

  // "solve" the pde
  SolutionFunction<dim> function;
  Vector<double>        solution(dofHandler.n_dofs());
  VectorTools::interpolate(dofHandler, function, solution);

  // create the particle handler
  Particles::ParticleHandler<dim> particleHandler(triangulation,
                                                  mapping,
                                                  ncomponents);

  // randomly place particles in the domain
  Functions::ConstantFunction<dim> pdf(1.0);
  Particles::Generators::probabilistic_locations(
    triangulation, pdf, true, nparticles, particleHandler, mapping);

  // map the conserved properties onto the particles
  dealii::Vector<double> interpolatedParticleQuantities(ncomponents *
                                                        nparticles);
  Particles::Utilities::interpolate_field_on_particles(
    dofHandler, particleHandler, solution, interpolatedParticleQuantities);

  // check to make sure we get what we expect
  unsigned int part = 0;
  for (auto iter = particleHandler.begin(); iter != particleHandler.end();
       ++iter, ++part)
    {
      deallog << "particle " << part << " quantities: ";
      for (unsigned int i = 0; i < ncomponents; ++i)
        {
          deallog << interpolatedParticleQuantities[part * ncomponents + i]
                  << ' ';
        }
      deallog << std::endl;
      deallog << "expected quantities: ";
      for (unsigned int i = 0; i < dim; ++i)
        {
          deallog << iter->get_location()[i] << ' ';
        }
      deallog << " 1" << std::endl;
      deallog << std::endl;
    }

  // Now try to set the interpolated quantities as particle
  // properties. This used to fail before because we had forgotten to
  // the set the property_pool pointer of the particles that were
  // bulk-inserted before.
  part = 0;
  for (auto iter = particleHandler.begin(); iter != particleHandler.end();
       ++iter, ++part)
    {
      std::vector<double> quantities(ncomponents);
      for (unsigned int i = 0; i < ncomponents; ++i)
        {
          quantities[i] =
            interpolatedParticleQuantities[part * ncomponents + i];
        }

      iter->set_properties(quantities);
    }

  // Finally output these properties
  for (const auto particle : particleHandler)
    {
      deallog << "Particle " << particle.get_id() << ": ";
      for (const auto p : particle.get_properties())
        deallog << p << ' ';
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  test<1>();
  test<2>();
  test<3>();
}
