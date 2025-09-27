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

#ifndef dealii_data_out_faces_h
#define dealii_data_out_faces_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/data_out.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutFacesImplementation
  {
    /**
     * A derived class for use in the DataOutFaces class. This is a class for
     * the AdditionalData kind of data structure discussed in the
     * documentation of the WorkStream context.
     */
    template <int dim, int spacedim>
    struct ParallelData
      : public internal::DataOutImplementation::ParallelDataBase<dim, spacedim>
    {
      ParallelData(const unsigned int               n_datasets,
                   const unsigned int               n_subdivisions,
                   const std::vector<unsigned int> &n_postprocessor_outputs,
                   const Mapping<dim, spacedim>    &mapping,
                   const std::vector<
                     std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                                    &finite_elements,
                   const UpdateFlags update_flags);

      std::vector<Point<spacedim>> patch_evaluation_points;
    };
  } // namespace DataOutFacesImplementation
} // namespace internal


/**
 * This class generates output from faces of a triangulation. It might be used
 * to generate output only for the surface of the triangulation (this is the
 * default of this class), or for all faces of active cells, as specified in
 * the constructor. The output of this class is a set of patches (as defined
 * by the class DataOutBase::Patch()), one for each face for which output is
 * to be generated. These patches can then be written in several graphical
 * data formats by the functions of the underlying classes.
 *
 * <h3>Interface</h3>
 *
 * The interface of this class is copied from the DataOut class. Furthermore,
 * they share the common parent class DataOut_DoFData. See the reference of
 * these two classes for a discussion of the interface.
 *
 *
 * <h3>Extending this class</h3>
 *
 * The sequence of faces to generate patches from is generated in the same way
 * as in the DataOut class; see there for a description of the respective
 * interface. The functions generating the sequence of faces which shall be
 * used to generate output, are called first_face() and next_face().
 *
 * Since we need to initialize objects of type FEValues with the faces
 * generated from these functions, it is not sufficient that they only return
 * face iterators. Rather, we need a pair of cell and the number of the face,
 * as the values of finite element fields need not necessarily be unique on a
 * face (think of discontinuous finite elements, where the value of the finite
 * element field depend on the direction from which you approach a face, thus
 * it is necessary to use a pair of cell and face, rather than only a face
 * iterator). Therefore, this class defines an @p alias which creates a type
 * @p FaceDescriptor that is an abbreviation for a pair of cell iterator and
 * face number. The functions @p first_face and @p next_face operate on
 * objects of this type.
 *
 * Extending this class might, for example, be useful if you only want output
 * from certain portions of the boundary, e.g. as indicated by the boundary
 * indicator of the respective faces. However, it is also conceivable that one
 * generates patches not from boundary faces, but from interior faces that are
 * selected due to other criteria; one application might be to use only those
 * faces where one component of the solution attains a certain value, in order
 * to display the values of other solution components on these faces. Other
 * applications certainly exist, for which the author is not imaginative
 * enough.
 *
 * @todo Reimplement this whole class using actual FEFaceValues and
 * MeshWorker.
 *
 * @ingroup output
 */
template <int dim, int spacedim = dim>
class DataOutFaces : public DataOut_DoFData<dim, dim - 1, spacedim, spacedim>
{
  static_assert(dim == spacedim, "Not implemented for dim != spacedim.");

public:
  /**
   * Dimension parameters for the patches.
   */
  static constexpr int patch_dim      = dim - 1;
  static constexpr int patch_spacedim = spacedim;

  /**
   * Alias to the iterator type of the dof handler class under
   * consideration.
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
      cell_iterator;

  /**
   * Constructor.
   *
   * @param[in] surface_only If `true`, then this class only generates
   *   output on faces that lie on the boundary of the domain. This
   *   is typically what this class is used for: To output information
   *   about the solution, fluxes, and other quantities that live on
   *   the boundary of the domain. On the other hand, it is sometimes
   *   useful to also visualize data on internal faces. This is
   *   facilitated by setting this argument to `false`.
   */
  DataOutFaces(const bool surface_only = true);

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. A
   * patch is, in essence, some intermediate representation of the data on
   * each face of a triangulation and DoFHandler object that can then be used
   * to write files in some format that is readable by visualization programs.
   *
   * You can find an overview of the use of this function in the general
   * documentation of this class. An example is also provided in the
   * documentation of this class's base class DataOut_DoFData.
   *
   * @param n_subdivisions See DataOut::build_patches() for an extensive
   * description of this parameter.
   */
  virtual void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * Same as above, except that the additional first parameter defines a
   * mapping that is to be used in the generation of output. If
   * <tt>n_subdivisions>1</tt>, the points interior of subdivided patches
   * which originate from cells at the boundary of the domain can be computed
   * using the mapping, i.e. a higher order mapping leads to a representation
   * of a curved boundary by using more subdivisions.
   *
   * Even for non-curved cells the mapping argument can be used for the
   * Eulerian mappings (see class MappingQ1Eulerian) where a mapping is used
   * not only to determine the position of points interior to a cell, but also
   * of the vertices.  It offers an opportunity to watch the solution on a
   * deformed triangulation on which the computation was actually carried out,
   * even if the mesh is internally stored in its undeformed configuration and
   * the deformation is only tracked by an additional vector that holds the
   * deformation of each vertex.
   *
   * @todo The @p mapping argument should be replaced by a
   * hp::MappingCollection in case of a DoFHandler with hp-capabilities.
   */
  virtual void
  build_patches(const Mapping<dim, spacedim> &mapping,
                const unsigned int            n_subdivisions = 0);

  /**
   * Declare a way to describe a face which we would like to generate output
   * for. The usual way would, of course, be to use an object of type
   * <tt>DoFHandler<dim>::face_iterator</tt>, but since we have to describe
   * faces to objects of type FEValues, we can only represent faces by pairs
   * of a cell and the number of the face. This pair is here aliased to a name
   * that is better to type.
   */
  using FaceDescriptor = typename std::pair<cell_iterator, unsigned int>;


  /**
   * Return the first face which we want output for. The default
   * implementation returns the first face of a (locally owned) active cell
   * or, if the @p surface_only option was set in the destructor (as is the
   * default), the first such face that is located on the boundary.
   *
   * If you want to use a different logic to determine which faces should
   * contribute to the creation of graphical output, you can overload this
   * function in a derived class.
   */
  virtual FaceDescriptor
  first_face();

  /**
   * Return the next face after which we want output for. If there are no more
   * faces, <tt>dofs->end()</tt> is returned as the first component of the
   * return value.
   *
   * The default implementation returns the next face of a (locally owned)
   * active cell, or the next such on the boundary (depending on whether the
   * @p surface_only option was provided to the constructor).
   *
   * This function traverses the mesh active cell by active cell (restricted to
   * locally owned cells), and then through all faces of the cell. As a result,
   * interior faces are output twice, a feature that is useful for
   * discontinuous Galerkin methods or if a DataPostprocessor is used that
   * might produce results that are discontinuous between cells).
   *
   * This function can be overloaded in a derived class to select a
   * different set of faces. Note that the default implementation assumes that
   * the given @p face is active, which is guaranteed as long as first_face()
   * is also used from the default implementation. Overloading only one of the
   * two functions should be done with care.
   */
  virtual FaceDescriptor
  next_face(const FaceDescriptor &face);

private:
  /**
   * Parameter deciding between surface meshes and full wire basket.
   */
  const bool surface_only;

  /**
   * Build one patch. This function is called in a WorkStream context.
   */
  void
  build_one_patch(
    const FaceDescriptor *cell_and_face,
    internal::DataOutFacesImplementation::ParallelData<dim, spacedim> &data,
    DataOutBase::Patch<patch_dim, patch_spacedim>                     &patch);
};


DEAL_II_NAMESPACE_CLOSE

#endif
