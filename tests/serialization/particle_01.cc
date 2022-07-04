// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2022 by the deal.II authors
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


// check and illustrate the serialization process for individual particles
// outside a ParticleHandler class

#include <deal.II/particles/particle.h>
#include <deal.II/particles/property_pool.h>

#include "serialization.h"



template <int dim, int spacedim>
void
test()
{
  // save data to archive
  std::ostringstream oss;
  {
    Point<spacedim> location;
    Point<dim>      reference_location;

    location[spacedim - 1]                    = 0.3;
    reference_location[dim - 1]               = 0.5;
    const unsigned int                     id = 6;
    Particles::PropertyPool<dim, spacedim> property_pool(2);

    Particles::Particle<dim, spacedim> particle(location,
                                                reference_location,
                                                id);
    particle.set_property_pool(property_pool);
    particle.get_properties()[0] = 0.1;
    particle.get_properties()[1] = 0.7;

    deallog << "Before serialization particle id " << particle.get_id()
            << " has location " << particle.get_location()
            << ", has reference location " << particle.get_reference_location()
            << ", and has properties " << particle.get_properties()[0] << ' '
            << particle.get_properties()[1] << std::endl;

    boost::archive::text_oarchive oa(oss, boost::archive::no_header);
    oa << particle;

    // archive and stream closed when
    // destructors are called
  }
  deallog << "Serialized string: " << oss.str() << std::endl;

  // verify correctness of the serialization.
  {
    std::istringstream            iss(oss.str());
    boost::archive::text_iarchive ia(iss, boost::archive::no_header);

    Particles::PropertyPool<dim, spacedim> property_pool(2);

    Particles::Particle<dim, spacedim> particle;
    particle.set_property_pool(property_pool);

    ia >> particle;

    deallog << "After serialization particle id " << particle.get_id()
            << " has location " << particle.get_location()
            << ", has reference location " << particle.get_reference_location()
            << ", and has properties " << particle.get_properties()[0] << ' '
            << particle.get_properties()[1] << std::endl;
  }

  deallog << "OK" << std::endl << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  deallog.push("2d/2d");
  test<2, 2>();
  deallog.pop();
  deallog.push("2d/3d");
  test<2, 3>();
  deallog.pop();
  deallog.push("3d/3d");
  test<3, 3>();
  deallog.pop();
}
