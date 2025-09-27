// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_element_access.h>

#include <deal.II/non_matching/mesh_classifier.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  namespace internal
  {
    namespace MeshClassifierImplementation
    {
      DeclExceptionMsg(
        ExcReclassifyNotCalled,
        "The Triangulation has not been classified. You need to call the "
        "reclassify()-function before using this function.");

      DeclExceptionMsg(
        ExcTriangulationMismatch,
        "The incoming cell does not belong to the triangulation passed to "
        "the constructor.");

      /**
       * Return LocationToLevelSet::inside/outside if all values in incoming
       * vector are negative/positive, otherwise return
       * LocationToLevelSet::intersected.
       */
      template <typename VectorType>
      LocationToLevelSet
      location_from_dof_signs(const VectorType &local_levelset_values)
      {
        const auto min_max_element =
          std::minmax_element(local_levelset_values.begin(),
                              local_levelset_values.end());

        if (*min_max_element.second < 0)
          return LocationToLevelSet::inside;
        if (0 < *min_max_element.first)
          return LocationToLevelSet::outside;

        return LocationToLevelSet::intersected;
      }



      /**
       * The concrete LevelSetDescription used when the level set function is
       * described as a (DoFHandler, Vector)-pair.
       */
      template <int dim, typename VectorType>
      class DiscreteLevelSetDescription : public LevelSetDescription<dim>
      {
      public:
        /**
         * Constructor.
         */
        DiscreteLevelSetDescription(const DoFHandler<dim> &dof_handler,
                                    const VectorType      &level_set);

        /**
         * Return the FECollection of the DoFHandler passed to the constructor.
         */
        const hp::FECollection<dim> &
        get_fe_collection() const override;

        /**
         * Return the active FE index of the DoFCellAccessor associated with the
         * DoFHandler and the incoming cell in the triangulation.
         */
        unsigned int
        active_fe_index(const typename Triangulation<dim>::active_cell_iterator
                          &cell) const override;

        /**
         * Writes the local face dofs of the global level set vector to
         * @p local_levelset_values.
         */
        void
        get_local_level_set_values(
          const typename Triangulation<dim>::active_cell_iterator &cell,
          const unsigned int                                       face_index,
          Vector<double> &local_levelset_values) override;

      private:
        /**
         * Pointer to the DoFHandler associated with the level set function.
         */
        const ObserverPointer<const DoFHandler<dim>> dof_handler;

        /**
         * Pointer to the vector containing the level set function's global dof
         * values.
         */
        const ObserverPointer<const VectorType> level_set;
      };



      template <int dim, typename VectorType>
      DiscreteLevelSetDescription<dim, VectorType>::DiscreteLevelSetDescription(
        const DoFHandler<dim> &dof_handler,
        const VectorType      &level_set)
        : dof_handler(&dof_handler)
        , level_set(&level_set)
      {}



      template <int dim, typename VectorType>
      const hp::FECollection<dim> &
      DiscreteLevelSetDescription<dim, VectorType>::get_fe_collection() const
      {
        return dof_handler->get_fe_collection();
      }



      template <int dim, typename VectorType>
      void
      DiscreteLevelSetDescription<dim, VectorType>::get_local_level_set_values(
        const typename Triangulation<dim>::active_cell_iterator &cell,
        const unsigned int                                       face_index,
        Vector<double> &local_levelset_values)
      {
        const auto cell_with_dofs = cell->as_dof_handler_iterator(*dof_handler);

        const unsigned int n_dofs_per_face =
          dof_handler->get_fe().n_dofs_per_face();
        std::vector<types::global_dof_index> dof_indices(n_dofs_per_face);
        cell_with_dofs->face(face_index)->get_dof_indices(dof_indices);

        local_levelset_values.reinit(dof_indices.size());

        for (unsigned int i = 0; i < dof_indices.size(); i++)
          local_levelset_values[i] =
            dealii::internal::ElementAccess<VectorType>::get(*level_set,
                                                             dof_indices[i]);
      }



      template <int dim, typename VectorType>
      unsigned int
      DiscreteLevelSetDescription<dim, VectorType>::active_fe_index(
        const typename Triangulation<dim>::active_cell_iterator &cell) const
      {
        const auto cell_with_dofs = cell->as_dof_handler_iterator(*dof_handler);

        return cell_with_dofs->active_fe_index();
      }


      /**
       * The concrete LevelSetDescription used when the level set function is
       * described by a Function.
       */
      template <int dim>
      class AnalyticLevelSetDescription : public LevelSetDescription<dim>
      {
      public:
        /**
         * Constructor. Takes the Function that describes the geometry and the
         * element that this function should be interpolated to.
         */
        AnalyticLevelSetDescription(const Function<dim>      &level_set,
                                    const FiniteElement<dim> &element);

        /**
         * Returns the finite element passed to the constructor wrapped in a
         * collection.
         */
        const hp::FECollection<dim> &
        get_fe_collection() const override;

        /**
         * Returns 0, since there is always a single element in the
         * FECollection.
         */
        unsigned int
        active_fe_index(const typename Triangulation<dim>::active_cell_iterator
                          &cell) const override;

        /**
         * Return the level set function evaluated at the real space face
         * support points of the finite element passed to the constructor.
         */
        void
        get_local_level_set_values(
          const typename Triangulation<dim>::active_cell_iterator &cell,
          const unsigned int                                       face_index,
          Vector<double> &local_levelset_values) override;

      private:
        /**
         * Pointer to the level set function.
         */
        const ObserverPointer<const Function<dim>> level_set;

        /**
         * Collection containing the single element which we locally interpolate
         * the level set function to.
         */
        const hp::FECollection<dim> fe_collection;

        /**
         * FEFaceValues object used to transform the support points on a face to
         * real space.
         */
        FEFaceValues<dim> fe_face_values;
      };



      template <int dim>
      AnalyticLevelSetDescription<dim>::AnalyticLevelSetDescription(
        const Function<dim>      &level_set,
        const FiniteElement<dim> &element)
        : level_set(&level_set)
        , fe_collection(element)
        , fe_face_values(element,
                         Quadrature<dim - 1>(
                           element.get_unit_face_support_points()),
                         update_quadrature_points)
      {}



      template <int dim>
      void
      AnalyticLevelSetDescription<dim>::get_local_level_set_values(
        const typename Triangulation<dim>::active_cell_iterator &cell,
        const unsigned int                                       face_index,
        Vector<double> &local_levelset_values)
      {
        AssertDimension(local_levelset_values.size(),
                        fe_face_values.n_quadrature_points);

        fe_face_values.reinit(cell, face_index);
        const std::vector<Point<dim>> &points =
          fe_face_values.get_quadrature_points();

        for (unsigned int i = 0; i < points.size(); i++)
          local_levelset_values[i] = level_set->value(points[i]);
      }



      template <int dim>
      const hp::FECollection<dim> &
      AnalyticLevelSetDescription<dim>::get_fe_collection() const
      {
        return fe_collection;
      }



      template <int dim>
      unsigned int
      AnalyticLevelSetDescription<dim>::active_fe_index(
        const typename Triangulation<dim>::active_cell_iterator &) const
      {
        return 0;
      }
    } // namespace MeshClassifierImplementation
  }   // namespace internal



  template <int dim>
  template <typename VectorType>
  MeshClassifier<dim>::MeshClassifier(const DoFHandler<dim> &dof_handler,
                                      const VectorType      &level_set)
    : triangulation(&dof_handler.get_triangulation())
    , level_set_description(
        std::make_unique<internal::MeshClassifierImplementation::
                           DiscreteLevelSetDescription<dim, VectorType>>(
          dof_handler,
          level_set))
  {
#ifdef DEAL_II_WITH_LAPACK
    const hp::FECollection<dim> &fe_collection =
      dof_handler.get_fe_collection();
    for (unsigned int i = 0; i < fe_collection.size(); i++)
      {
        // The level set function must be scalar.
        AssertDimension(fe_collection[i].n_components(), 1);

        Assert(fe_collection[i].has_face_support_points(),
               ExcMessage(
                 "The elements in the FECollection of the incoming DoFHandler "
                 "must have face support points."));
      }
#else
    AssertThrow(false, ExcNeedsLAPACK());
#endif
  }



  template <int dim>
  MeshClassifier<dim>::MeshClassifier(const Triangulation<dim> &triangulation,
                                      const Function<dim>      &level_set,
                                      const FiniteElement<dim> &element)
    : triangulation(&triangulation)
    , level_set_description(
        std::make_unique<internal::MeshClassifierImplementation::
                           AnalyticLevelSetDescription<dim>>(level_set,
                                                             element))
  {
    // The level set function must be scalar.
    AssertDimension(level_set.n_components, 1);
    AssertDimension(element.n_components(), 1);
  }



  template <int dim>
  void
  MeshClassifier<dim>::reclassify()
  {
    initialize();
    cell_locations.assign(triangulation->n_active_cells(),
                          LocationToLevelSet::unassigned);
    face_locations.assign(triangulation->n_raw_faces(),
                          LocationToLevelSet::unassigned);

    // Loop over all cells and determine the location of all non artificial
    // cells and faces.
    for (const auto &cell : triangulation->active_cell_iterators())
      if (!cell->is_artificial())
        {
          const LocationToLevelSet face0_location =
            determine_face_location_to_levelset(cell, 0);

          face_locations[cell->face(0)->index()] = face0_location;
          LocationToLevelSet cell_location       = face0_location;

          for (unsigned int f = 1; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
              const LocationToLevelSet face_location =
                determine_face_location_to_levelset(cell, f);

              face_locations[cell->face(f)->index()] = face_location;

              if (face_location != face0_location)
                cell_location = LocationToLevelSet::intersected;
            }
          cell_locations[cell->active_cell_index()] = cell_location;
        }
  }



  template <int dim>
  LocationToLevelSet
  MeshClassifier<dim>::determine_face_location_to_levelset(
    const typename Triangulation<dim>::active_cell_iterator &cell,
    const unsigned int                                       face_index)
  {
    // The location of the face might already be computed on the neighboring
    // cell. If this is the case we just return the value.
    const LocationToLevelSet location =
      face_locations.at(cell->face(face_index)->index());
    if (location != LocationToLevelSet::unassigned)
      return location;

    // Determine the location by changing basis to FE_Bernstein and checking
    // the signs of the dofs.
    const unsigned int fe_index = level_set_description->active_fe_index(cell);
    const unsigned int n_local_dofs =
      lagrange_to_bernstein_face[fe_index][face_index].m();

    Vector<double> local_levelset_values(n_local_dofs);
    level_set_description->get_local_level_set_values(cell,
                                                      face_index,
                                                      local_levelset_values);

    const FiniteElement<dim> &fe =
      level_set_description->get_fe_collection()[fe_index];

    const FE_Q_iso_Q1<dim> *fe_q_iso_q1 =
      dynamic_cast<const FE_Q_iso_Q1<dim> *>(&fe);

    const FE_Poly<dim> *fe_poly = dynamic_cast<const FE_Poly<dim> *>(&fe);

    const bool is_linear = fe_q_iso_q1 != nullptr ||
                           (fe_poly != nullptr && fe_poly->get_degree() == 1);

    // shortcut for linear elements
    if (is_linear)
      {
        return internal::MeshClassifierImplementation::location_from_dof_signs(
          local_levelset_values);
      }

    lagrange_to_bernstein_face[fe_index][face_index].solve(
      local_levelset_values);

    return internal::MeshClassifierImplementation::location_from_dof_signs(
      local_levelset_values);
  }



  template <int dim>
  LocationToLevelSet
  MeshClassifier<dim>::location_to_level_set(
    const typename Triangulation<dim>::cell_iterator &cell) const
  {
    Assert(cell_locations.size() == triangulation->n_active_cells(),
           internal::MeshClassifierImplementation::ExcReclassifyNotCalled());
    Assert(&cell->get_triangulation() == triangulation,
           internal::MeshClassifierImplementation::ExcTriangulationMismatch());

    return cell_locations.at(cell->active_cell_index());
  }



  template <int dim>
  LocationToLevelSet
  MeshClassifier<dim>::location_to_level_set(
    const typename Triangulation<dim>::cell_iterator &cell,
    const unsigned int                                face_index) const
  {
    AssertIndexRange(face_index, GeometryInfo<dim>::faces_per_cell);
    Assert(face_locations.size() == triangulation->n_raw_faces(),
           internal::MeshClassifierImplementation::ExcReclassifyNotCalled());
    Assert(&cell->get_triangulation() == triangulation,
           internal::MeshClassifierImplementation::ExcTriangulationMismatch());

    return face_locations.at(cell->face(face_index)->index());
  }



  template <int dim>
  void
  MeshClassifier<dim>::initialize()
  {
    const hp::FECollection<dim> &fe_collection =
      level_set_description->get_fe_collection();

    // The level set function must be scalar.
    AssertDimension(fe_collection.n_components(), 1);

    lagrange_to_bernstein_face.resize(fe_collection.size());

    for (unsigned int i = 0; i < fe_collection.size(); i++)
      {
        const FiniteElement<dim> &element = fe_collection[i];
        const FE_Q_Base<dim>     *fe_q =
          dynamic_cast<const FE_Q_Base<dim> *>(&element);
        Assert(fe_q != nullptr, ExcNotImplemented());

        const FE_Bernstein<dim> fe_bernstein(fe_q->get_degree());

        const unsigned int dofs_per_face = fe_q->dofs_per_face;
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++)
          {
            FullMatrix<double> face_interpolation_matrix(dofs_per_face,
                                                         dofs_per_face);

            fe_bernstein.get_face_interpolation_matrix(
              *fe_q, face_interpolation_matrix, f);
            lagrange_to_bernstein_face[i][f].reinit(dofs_per_face);
            lagrange_to_bernstein_face[i][f] = face_interpolation_matrix;
            lagrange_to_bernstein_face[i][f].compute_lu_factorization();
          }
      }
  }

} // namespace NonMatching

#include "non_matching/mesh_classifier.inst"

DEAL_II_NAMESPACE_CLOSE
