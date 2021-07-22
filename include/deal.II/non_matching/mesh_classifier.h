// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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

#ifndef dealii_non_matching_mesh_classifier
#define dealii_non_matching_mesh_classifier

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <memory>

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
   * one of values of LocationToLevelSet: inside, outside or intersected,
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
   * LocationToLevelSet location =
   *   classifier.location_to_level_set(triangulation.begin_active());
   * @endcode
   *
   * The level set function can either be described as a discrete function by a
   * (DoFHandler, Vector)-pair or as a general Function. In the case of a
   * discrete function, LocationToLevelSet for a given cell/face is determined
   * by looking at the local degrees of freedom. Since the Lagrange basis
   * functions are not positive definite, positive/negative definite
   * dof-values does not imply that the interpolated function is
   * positive/negative definite. Thus, to classify a cell/face this class
   * internally transforms the local dofs to a basis spanning the same
   * polynomial space but where definite dof-values implies a definite function.
   * Currently only the case of FE_Q-elements is implemented, where we
   * internally change basis to FE_Bernstein. When the level set function is
   * described as a Function, the level set function is locally interpolated to
   * an FE_Q element and we proceed in the same way as for the discrete level
   * set function.
   */
  template <int dim>
  class MeshClassifier : public Subscriptor
  {
  public:
    /**
     * Constructor. Takes a level set function described as a DoFHandler and a
     * Vector. The triangulation attached to DoFHandler is the one that will be
     * classified.
     */
    template <class VECTOR>
    MeshClassifier(const DoFHandler<dim> &level_set_dof_handler,
                   const VECTOR &         level_set);


    /**
     * Constructor. Takes the triangulation that should be classified, a
     * level set function described as a Function, and a scalar element that we
     * interpolate the Function to in order to classify each cell/face.
     *
     * @note The Function and the FiniteElement must both have a single component.
     */
    MeshClassifier(const Triangulation<dim> &triangulation,
                   const Function<dim> &     level_set,
                   const FiniteElement<dim> &element);

    /**
     * Perform the classification of the active cells and faces in the
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
     * level_set_description, sets up the local transformation matrices
     */
    void
    initialize();

    /**
     * Computes how the incoming cell is located relative to the level set
     * function.
     */
    LocationToLevelSet
    determine_cell_location_to_levelset(
      const typename Triangulation<dim>::active_cell_iterator &cell);

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
    const SmartPointer<const Triangulation<dim>> triangulation;

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
     * A map that for each cell stores how each active cell is located relative
     * to the level set function.
     */
    std::map<typename Triangulation<dim>::cell_iterator, LocationToLevelSet>
      cell_locations;

    /**
     * A map that stores how each active face is located relative to the level
     * set function.
     */
    std::map<typename Triangulation<dim>::face_iterator, LocationToLevelSet>
      face_locations;

    /**
     * For each element in the hp::FECollection returned by the
     * LevelSetDescription, this vector stores a transformation matrix to a
     * basis where positive/negative definite cell dofs implies that the
     * underlying function is positive/negative definite over the cell.
     */
    std::vector<LAPACKFullMatrix<double>> lagrange_to_bernstein_cell;

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
       * function independently of weather it is a Function or a
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
         * Return the local level set dofs of the incoming cell.
         */
        virtual void
        get_local_level_set_values(
          const typename Triangulation<dim>::active_cell_iterator &cell,
          Vector<double> &local_dofs) = 0;

        /**
         * Return the local level set dofs of the face with the given index on
         * the incoming cell.
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
