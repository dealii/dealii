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

#ifndef dealii_non_matching_mesh_classifier
#define dealii_non_matching_mesh_classifier

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  namespace internal
  {
    namespace MeshClassifierImplementation
    {
      template <int dim>
      class LevelSetDescription;
    } // namespace MeshClassifierImplementation
  }   // namespace internal


  /**
   * Type describing how a cell or a face is located relative to the zero
   * contour of a level set function, $\psi$. The values of the type correspond
   * to:
   *
   * inside      if $\psi(x) < 0$,
   * outside     if $\psi(x) > 0$,
   * intersected if $\psi(x)$ varies in sign,
   *
   * over the cell/face. The value "unassigned" is used to describe that the
   * location of a cell/face has not yet been determined.
   */
  enum class LocationToLevelSet
  {
    inside,
    outside,
    intersected,
    unassigned
  };


  /**
   * Class responsible for determining how the active cells and faces of a
   * triangulation relate to the sign of a level set function. When calling the
   * reclassify() function each of the active cells and faces are categorized as
   * one of the values of LocationToLevelSet: inside, outside or intersected,
   * depending on the sign of the level set function over the cell/face. This
   * information is typically required in immersed/cut finite element methods,
   * both when distributing degrees of freedom over the triangulation and when
   * the system is assembled. The given class would then be used in the
   * following way:
   *
   * @code
   * Vector<double> &level_set = ...
   *
   * MeshClassifier<dim> classifier(dof_handler, level_set);
   * classifier.reclassify();
   *
   * LocationToLevelSet location = classifier.location_to_level_set(cell);
   * @endcode
   *
   * The level set function can either be described as a discrete function by a
   * (DoFHandler, Vector)-pair or as a general Function. In the case of a
   * discrete function, LocationToLevelSet for a given face is determined
   * by looking at the local degrees of freedom on the face. Since the Lagrange
   * basis functions are not positive definite, positive/negative definite
   * dof values do not imply that the interpolated function is
   * positive/negative definite. Thus, to classify a face this class internally
   * transforms the local dofs to a basis spanning the same polynomial space
   * over the face but where definite dof values imply a definite function.
   * Currently, only the case of FE_Q-elements is implemented, where we
   * internally change basis to FE_Bernstein. For cells, LocationToLevelSet is
   * determined from the faces of the cell. That is, if all faces of the cell
   * are inside/outside the LocationToLevelSet of the cell is set to
   * inside/outside. LocationToLevelSet of the cell is set to intersected if at
   * least one face is intersected or if its faces have different
   * LocationToLevelSet. Note that, this procedure will incorrectly classify the
   * cell as inside/outside, if the mesh refinement is so low that the whole
   * zero-contour is contained in a single cell (so that none of its faces are
   * intersected).
   *
   * When the level set function is described as a Function, the level set
   * function is locally interpolated to an FE_Q element and we proceed in the
   * same way as for the discrete level set function.
   */
  template <int dim>
  class MeshClassifier : public EnableObserverPointer
  {
  public:
    /**
     * Constructor. Takes a level set function described as a DoFHandler and a
     * Vector. The triangulation attached to DoFHandler is the one that will be
     * classified.
     */
    template <typename VectorType>
    MeshClassifier(const DoFHandler<dim> &level_set_dof_handler,
                   const VectorType      &level_set);

    /**
     * Constructor. Takes the triangulation that should be classified, a
     * level set function described as a Function, and a scalar element that we
     * interpolate the Function to in order to classify each cell/face.
     *
     * @note The Function and the FiniteElement must both have a single component.
     */
    MeshClassifier(const Triangulation<dim> &triangulation,
                   const Function<dim>      &level_set,
                   const FiniteElement<dim> &element);

    /**
     * Perform the classification of the non artificial cells and faces in the
     * triangulation.
     */
    void
    reclassify();

    /**
     * Return how the incoming cell is located relative to the level set
     * function.
     */
    LocationToLevelSet
    location_to_level_set(
      const typename Triangulation<dim>::cell_iterator &cell) const;

    /**
     * Return how a face of the incoming cell is located relative to the level
     * set function.
     */
    LocationToLevelSet
    location_to_level_set(
      const typename Triangulation<dim>::cell_iterator &cell,
      const unsigned int                                face_index) const;

  private:
    /**
     * For each element in the hp::FECollection returned by
     * level_set_description, sets up the local transformation matrices.
     */
    void
    initialize();

    /**
     * Computes how the face with the given index on the incoming cell is
     * located relative to the level set function.
     */
    LocationToLevelSet
    determine_face_location_to_levelset(
      const typename Triangulation<dim>::active_cell_iterator &cell,
      const unsigned int                                       face_index);

    /**
     * Pointer to the triangulation that should be classified.
     */
    const ObserverPointer<const Triangulation<dim>> triangulation;

    /**
     * Pointer to an object that describes what we need to know about the
     * level set function. The underlying object will be of different type
     * depending on whether the level set function is discrete
     * (DoFHandler, Vector) or described by a Function.
     */
    const std::unique_ptr<
      internal::MeshClassifierImplementation::LevelSetDescription<dim>>
      level_set_description;

    /**
     * A vector that stores how each active cell is located relative to the
     * level set function, based on the cells active index.
     */
    std::vector<LocationToLevelSet> cell_locations;

    /**
     * A vector that stores how each active face is located relative to the
     * level set function, based on the face's global index.
     */
    std::vector<LocationToLevelSet> face_locations;

    /**
     * For each element in the hp::FECollection returned by the
     * LevelSetDescription, and for each local face, this vector stores a
     * transformation matrix to a basis where positive/negative definite
     * face dofs implies that the underlying function is positive/negative
     * definite over the face.
     */
    std::vector<
      std::array<LAPACKFullMatrix<double>, GeometryInfo<dim>::faces_per_cell>>
      lagrange_to_bernstein_face;
  };


  namespace internal
  {
    namespace MeshClassifierImplementation
    {
      /**
       * Abstract class that describes what we need to know about the level set
       * function independently of whether it is a Function or a
       * (DoFHandler, Vector)-pair.
       */
      template <int dim>
      class LevelSetDescription
      {
      public:
        /**
         * Destructor, declared to mark it virtual.
         */
        virtual ~LevelSetDescription() = default;

        /**
         * Return a collection to all the elements that are used to locally
         * describe the level set function.
         */
        virtual const hp::FECollection<dim> &
        get_fe_collection() const = 0;

        /**
         * Return the index of the element in the FECollection that we associate
         * with the level set function on the incoming cell.
         */
        virtual unsigned int
        active_fe_index(const typename Triangulation<dim>::active_cell_iterator
                          &cell) const = 0;

        /**
         * Fill the DoF values of the associated level set representation on the
         * face of the incoming cell into the vector provided in the last
         * argument.
         *
         * @note Since this function extracts the dofs on the face of the cell,
         * it assumes that the underlying element has face support points.
         */
        virtual void
        get_local_level_set_values(
          const typename Triangulation<dim>::active_cell_iterator &cell,
          const unsigned int                                       face_index,
          Vector<double> &local_dofs) = 0;
      };

    } // namespace MeshClassifierImplementation
  }   // namespace internal
} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif
