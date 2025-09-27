// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_resample_h
#define dealii_data_out_resample_h



#include <deal.II/base/config.h>

#include <deal.II/base/mpi_remote_point_evaluation.h>
#include <deal.II/base/partitioner.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * A DataOut-like class which interpolates values defined on one Triangulation
 * onto a second potentially unrelated Triangulation. By using this class,
 * one can output the result obtained on an unstructured mesh onto a
 * structured one or one can create a slice in 3d.
 *
 * The following code snippet shows the steps how to use the class when the
 * solution is given for a three dimensional triangulation and the result
 * should be outputted on a (2d) slice:
 * @code
 * // setup and first usage
 * DataOutResample<3, 2, 3> data_out(patch_tria,patch_mapping);
 * data_out.add_data_vector(dof_handler, vector, "solution");
 * data_out.build_patches(mapping);
 *
 * // ... no changes in triangulation and mapping -> reuse internal data
 * // structures
 * data_out.build_patches();
 *
 * // ... changes in triangulation or mapping -> reinitialize internal data
 * // structures
 * data_out.build_patches(mapping);
 * @endcode
 *
 * @note While the dimension of the two triangulations might differ, their
 *   space dimension need to coincide.
 */
template <int dim, int patch_dim, int spacedim>
class DataOutResample
  : public DataOut_DoFData<dim, patch_dim, spacedim, spacedim>
{
public:
  /**
   * Constructor taking the triangulation and mapping for which the patches
   * should be generated.
   */
  DataOutResample(const Triangulation<patch_dim, spacedim> &patch_tria,
                  const Mapping<patch_dim, spacedim>       &patch_mapping);

  /**
   * Update the @p mapping of original triangulation. One needs to call this
   * function if the mapping has changed. Just like in the DataOut context,
   * @p n_subdivisions determines how many "patches" this function will build
   * out of every cell.
   *
   * This function involves an expensive setup: evaluation points are generated
   * and their owners are determined, which is used to set up the communication
   * pattern.
   *
   * @note If you use the version of build_patches() that does not take a
   *   mapping, this function has to be called before its first usage.
   */
  void
  update_mapping(const Mapping<dim, spacedim> &mapping,
                 const unsigned int            n_subdivisions = 0);

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. A
   * patch is, in essence, some intermediate representation of the data on
   * each cell of a triangulation and DoFHandler object that can then be used
   * to write files in some format that is readable by visualization programs.
   *
   * Since this function calls internally at the beginning update_mapping(),
   * this function also involves an expensive setup.
   *
   * Just like in the DataOut::build_patches() context, @p n_subdivisions
   * determines how many "patches" this function will build out of every cell.
   */
  void
  build_patches(
    const Mapping<dim, spacedim> &mapping,
    const unsigned int            n_subdivisions = 0,
    const typename DataOut<patch_dim, spacedim>::CurvedCellRegion
      curved_region =
        DataOut<patch_dim, spacedim>::CurvedCellRegion::curved_boundary);

  /**
   * Just like the above function, this function builds a list of
   * patches to be written by the low-level functions of the base class.
   * However it skips the update of the mapping and reuses the one registered
   * via update_mapping(). This allows to skip the expensive setup of the
   * internal communication routines.
   *
   * @note This function can be only used if a mapping has been registered via
   *   update_mapping() or the other build_patches() function. The same
   *   n_subdivisions previously passed to these functions is used.
   */
  void
  build_patches(
    const typename DataOut<patch_dim, spacedim>::CurvedCellRegion
      curved_region =
        DataOut<patch_dim, spacedim>::CurvedCellRegion::curved_boundary);

protected:
  virtual const std::vector<typename DataOutBase::Patch<patch_dim, spacedim>> &
  get_patches() const override;

private:
  /**
   * Intermediate DoFHandler
   */
  DoFHandler<patch_dim, spacedim> patch_dof_handler;

  /**
   * Mapping used in connection with patch_tria.
   */
  const ObserverPointer<const Mapping<patch_dim, spacedim>> patch_mapping;

  /**
   * DataOut object that does the actual building of the patches.
   */
  DataOut<patch_dim, spacedim> patch_data_out;

  /**
   * Object to evaluate
   */
  Utilities::MPI::RemotePointEvaluation<dim, spacedim> rpe;

  /**
   * Partitioner to create internally distributed vectors.
   */
  std::shared_ptr<Utilities::MPI::Partitioner> partitioner;

  /**
   * Process local indices to access efficiently internal distributed vectors.
   */
  std::vector<types::global_dof_index> point_to_local_vector_indices;

  /**
   * Mapping of the original triangulation provided in update_mapping().
   */
  ObserverPointer<const Mapping<dim, spacedim>> mapping;
};

DEAL_II_NAMESPACE_CLOSE

#endif
