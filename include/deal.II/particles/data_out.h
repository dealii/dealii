// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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
#ifndef dealii_particles_data_out_h
#define dealii_particles_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  class ParticleHandler;

  /**
   * This class generates graphical output for the particles stored by a
   * ParticleHandler object. From a particle handler, it generates patches which
   * can then be used to write traditional output files. This class currently
   * only supports witing the particle position and their ID and does not allow
   * to write the properties attached to the particles
   *
   * @ingroup Particle
   *
   * @author Bruno Blais, Luca Heltai 2019
   */
  template <int dim, int spacedim = dim>
  class DataOut : public dealii::DataOutInterface<0, spacedim>
  {
  public:
    /**
     * Default constructor for the Particles::DataOut class.
     */
    DataOut() = default;

    /**
     * Destructor for the Particles::DataOut class.
     */
    ~DataOut() = default;


    /**
     * Build the patches for a given particle handler.
     *
     * @param [in] particles A particle handler for which the patches will be built.
     * A dim=0 patch is built for each particle. The position of the particle is
     * used to build the node position and the ID of the particle is added as a
     * single data element.
     *
     * @author Bruno Blais, Luca Heltai 2019
     */
    void
    build_patches(const Particles::ParticleHandler<dim, spacedim> &particles);

  protected:
    /**
     * Returns the patches built by the data_out class which was previously
     * built using a particle handler
     */
    virtual const std::vector<DataOutBase::Patch<0, spacedim>> &
    get_patches() const override;

    /**
     * Virtual function through which the names of data sets are obtained from
     * this class
     */
    virtual std::vector<std::string>
    get_dataset_names() const override;

  private:
    /**
     * This is a list of patches that is created each time build_patches() is
     * called. These patches are used in the output routines of the base
     * classes.
     */
    std::vector<DataOutBase::Patch<0, spacedim>> patches;

    /**
     * A list of field names for all data components stored in patches.
     */
    std::vector<std::string> dataset_names;
  };

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
