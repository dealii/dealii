// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_particles_data_out_h
#define dealii_particles_data_out_h

#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>

#include <deal.II/numerics/data_component_interpretation.h>

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
   * can then be used to write traditional output files.
   *
   * The class does, in essence, for particles what ::DataOut does for
   * finite element fields, using a similar interface. It is used in
   * step-19, for example.
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
     * @param [in] data_component_names An optional vector of strings that
     * describe the properties of each particle. Particle properties will only
     * be written if this vector
     * is provided.
     * @param [in] data_component_interpretations An optional vector that
     * controls if the particle properties are interpreted as scalars, vectors,
     * or tensors. Has to be of the same length as @p data_component_names.
     */
    void
    build_patches(const Particles::ParticleHandler<dim, spacedim> &particles,
                  const std::vector<std::string> &data_component_names = {},
                  const std::vector<
                    DataComponentInterpretation::DataComponentInterpretation>
                    &data_component_interpretations = {});

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


    /**
     * Overload of the respective DataOutInterface::get_nonscalar_data_ranges()
     * function. See there for a more extensive documentation.
     * This function is a reimplementation of the function
     * DataOut_DoFData::get_nonscalar_data_ranges().
     */
    virtual std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
    get_nonscalar_data_ranges() const override;

  private:
    /**
     * This is a vector of patches that is created each time build_patches() is
     * called. These patches are used in the output routines of the base
     * classes.
     */
    std::vector<DataOutBase::Patch<0, spacedim>> patches;

    /**
     * A vector of field names for all data components stored in patches.
     */
    std::vector<std::string> dataset_names;

    /**
     * A vector that for each of the data components of the
     * current data set indicates whether they are scalar fields, parts of a
     * vector-field, or any of the other supported kinds of data.
     */
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretations;
  };

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
