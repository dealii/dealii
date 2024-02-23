// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_rotation_h
#define dealii_data_out_rotation_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutRotationImplementation
  {
    /**
     * A derived class for use in the DataOutFaces class. This is a class for
     * the AdditionalData kind of data structure discussed in the
     * documentation of the WorkStream class.
     */
    template <int dim, int spacedim>
    struct ParallelData
      : public internal::DataOutImplementation::ParallelDataBase<dim, spacedim>
    {
      ParallelData(const unsigned int               n_datasets,
                   const unsigned int               n_subdivisions,
                   const unsigned int               n_patches_per_circle,
                   const std::vector<unsigned int> &n_postprocessor_outputs,
                   const Mapping<dim, spacedim>    &mapping,
                   const std::vector<
                     std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                                    &finite_elements,
                   const UpdateFlags update_flags);

      const unsigned int n_patches_per_circle;

      std::vector<Point<spacedim>> patch_evaluation_points;
    };
  } // namespace DataOutRotationImplementation
} // namespace internal



/**
 * This class generates output in the full domain of computations that were
 * done using rotational symmetry of domain and solution. In particular, if a
 * computation of a three dimensional problem with rotational symmetry around
 * the @p z-axis (i.e. in the @p r-z-plane) was done, then this class can be
 * used to generate the output in the original @p x-y-z space. In order to do
 * so, it generates from each cell in the computational mesh a cell in the
 * space with dimension one greater than that of the DoFHandler object. The
 * resulting output will then consist of hexahedra forming an object that has
 * rotational symmetry around the z-axis. As most graphical programs can not
 * represent ring-like structures, the angular (rotation) variable is
 * discretized into a finite number of intervals as well; the number of these
 * intervals must be given to the @p build_patches function. It is noted,
 * however, that while this function generates nice pictures of the whole
 * domain, it often produces <em>very</em> large output files.
 *
 *
 * <h3>Interface</h3>
 *
 * The interface of this class is copied from the DataOut class. Furthermore,
 * they share the common parent class DataOut_DoFData(). See the reference of
 * these two classes for a discussion of the interface and how to extend it by
 * deriving further classes from this class.
 *
 *
 * <h3>Details for 1d computations</h3>
 *
 * The one coordinate in the triangulation used by the DoFHandler object
 * passed to this class is taken as the radial variable, and the output will
 * then be either a circle or a ring domain. It is in the user's
 * responsibility to assure that the radial coordinate only attains
 * non-negative values.
 *
 *
 * <h3>Details for 2d computations</h3>
 *
 * We consider the computation (represented by the DoFHandler object that is
 * attached to this class) to have happened in the @p r-z-plane, where @p r is
 * the radial variable and @p z denotes the axis of revolution around which
 * the solution is symmetric. The output is in @p x-y-z space, where the
 * radial dependence is transformed to the @p x-y plane. At present, it is not
 * possible to exchange the meaning of the first and second variable of the
 * plane in which the simulation was made, i.e. generate output from a
 * simulation where the first variable denoted the symmetry axis, and the
 * second denoted the radial variable. You have to take that into account when
 * first programming your application.
 *
 * It is in the responsibility of the user to make sure that the radial
 * variable attains only non-negative values.
 *
 * @ingroup output
 */
template <int dim, int spacedim = dim>
class DataOutRotation
  : public DataOut_DoFData<dim, dim + 1, spacedim, spacedim + 1>
{
  static_assert(dim == spacedim, "Not implemented for dim != spacedim.");

public:
  /**
   * Dimension parameters for the patches.
   */
  static constexpr int patch_dim      = dim + 1;
  static constexpr int patch_spacedim = spacedim + 1;

  /**
   * Typedef to the iterator type of the dof handler class under
   * consideration.
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
      cell_iterator;

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. A
   * patch is, in essence, some intermediate representation of the data on
   * each cell of a triangulation and DoFHandler object that can then be used
   * to write files in some format that is readable by visualization programs.
   *
   * You can find an overview of the use of this function in the general
   * documentation of this class. An example is also provided in the
   * documentation of this class's base class DataOut_DoFData.
   *
   * @param n_patches_per_circle Denotes into how many intervals the angular
   * (rotation) variable is to be subdivided.
   *
   * @param n_subdivisions See DataOut::build_patches() for an extensive
   * description of this parameter.
   */
  virtual void
  build_patches(const unsigned int n_patches_per_circle,
                const unsigned int n_subdivisions = 0);

  /**
   * Return the first cell which we want output for. The default
   * implementation returns the first
   * @ref GlossActive "active cell",
   * but you might want to return other cells in a derived class.
   */
  virtual cell_iterator
  first_cell();

  /**
   * Return the next cell after @p cell which we want output for. If there are
   * no more cells, <tt>dofs->end()</tt> shall be returned.
   *
   * The default implementation returns the next active cell, but you might
   * want to return other cells in a derived class. Note that the default
   * implementation assumes that the given @p cell is active, which is
   * guaranteed as long as @p first_cell is also used from the default
   * implementation. Overloading only one of the two functions might not be a
   * good idea.
   */
  virtual cell_iterator
  next_cell(const cell_iterator &cell);

  /**
   * Exception
   */
  DeclException1(ExcRadialVariableHasNegativeValues,
                 double,
                 << "You are attempting to use this class on a triangulation "
                    "in which some vertices have a negative radial coordinate "
                    "value of "
                 << arg1
                 << ". If you rotate such a triangulation around an "
                    "axis, you will get (dim+1)-dimensional meshes "
                    "that are not likely what you hoped to see.");

private:
  /**
   * Build all of the patches that correspond to the cell given in the first
   * argument. Use the second argument as scratch space for parallel
   * invocation in WorkStream, and put the results into the last argument.
   */
  void
  build_one_patch(
    const cell_iterator                                                  *cell,
    internal::DataOutRotationImplementation::ParallelData<dim, spacedim> &data,
    std::vector<DataOutBase::Patch<patch_dim, patch_spacedim>> &my_patches);
};


DEAL_II_NAMESPACE_CLOSE

#endif
