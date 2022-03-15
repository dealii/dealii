// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <deal.II/non_matching/fe_values.h>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  RegionUpdateFlags::RegionUpdateFlags()
    : inside(update_default)
    , outside(update_default)
    , surface(update_default)
  {}


  namespace internal
  {
    namespace FEValuesImplementation
    {
      DeclExceptionMsg(
        ExcCellNotSet,
        "The set_active_cell function has to be called before calling this function.");


      /**
       * This class evaluates a function defined by a solution vector and a
       * DoFHandler transformed to reference space. To be precise, if we let
       * $\hat{x}$ be a point on the reference cell, this class implements the
       * function
       *
       * $\hat{f}(\hat{x}) = \sum_{j=0}^{n-1} f_j \hat{\phi}_j(\hat{x})$,
       *
       * where $f_j$ are the local solution values and $\hat{\phi}_j(\hat(x))$
       * are the local reference space shape functions. The gradient and Hessian
       * of this function are thus derivatives with respect to the reference
       * space coordinates, $\hat{x}_0, \hat{x}_1, \ldots$.
       *
       * Note that this class is similar to FEFieldFunction, but that
       * FEFieldFunction implements the following function on a given cell, $K$,
       *
       * $f(x) = \sum_{j=0}^{n-1} f_j \hat{\phi}_j(F_K^{-1}(x))$,
       *
       * which has the same coefficients but uses real space basis functions.
       * Here, $F_K$ is the mapping from the reference cell to the real cell.
       *
       * Before calling the value/gradient/hessian function, the set_active_cell
       * function must be called to specify which cell the function should be
       * evaluated on.
       */
      template <int dim, class VectorType = Vector<double>>
      class RefSpaceFEFieldFunction : public CellWiseFunction<dim>
      {
      public:
        /**
         * Constructor. Takes the solution vector and the associated DoFHandler.
         *
         * Pointers to the input arguments are stored internally, so they must
         * have a longer lifetime than the created RefSpaceFEFieldFunction
         * object.
         */
        RefSpaceFEFieldFunction(const DoFHandler<dim> &dof_handler,
                                const VectorType &     dof_values);

        /**
         * @copydoc CellWiseFunction::set_active_cell()
         */
        void
        set_active_cell(const typename Triangulation<dim>::active_cell_iterator
                          &cell) override;

        /**
         * @copydoc Function::value()
         *
         * @note The set_active_cell function must be called before this function.
         * The incoming point should be on the reference cell, but this is not
         * checked.
         */
        double
        value(const Point<dim> & point,
              const unsigned int component = 0) const override;

        /**
         * @copydoc Function::gradient()
         *
         * @note The set_active_cell function must be called before this function.
         * The incoming point should be on the reference cell, but this is not
         * checked.
         */
        Tensor<1, dim>
        gradient(const Point<dim> & point,
                 const unsigned int component = 0) const override;

        /**
         * @copydoc Function::hessian()
         *
         * @note The set_active_cell function must be called before this function.
         * The incoming point should be on the reference cell, but this is not
         * checked.
         */
        SymmetricTensor<2, dim>
        hessian(const Point<dim> & point,
                const unsigned int component = 0) const override;

      private:
        /**
         * Return whether the set_active_cell function has been called.
         */
        bool
        cell_is_set() const;

        /**
         * Pointer to the DoFHandler passed to the constructor.
         */
        const SmartPointer<const DoFHandler<dim>> dof_handler;

        /**
         * Pointer to the vector of solution coefficients passed to the
         * constructor.
         */
        const SmartPointer<const VectorType> global_dof_values;

        /**
         * Pointer to the element associated with the cell in the last call to
         * set_active_cell().
         */
        SmartPointer<const FiniteElement<dim>> element;

        /**
         * DOF-indices of the cell in the last call to set_active_cell()
         */
        std::vector<types::global_dof_index> local_dof_indices;

        /**
         * Local solution values of the cell in the last call to
         * set_active_cell()
         */
        std::vector<typename VectorType::value_type> local_dof_values;

        /**
         * Description of the 1D polynomial basis for tensor product elements
         * used for the fast path of this class using tensor product
         * evaluators.
         */
        std::vector<Polynomials::Polynomial<double>> poly;

        /**
         * Renumbering for the tensor-product evaluator in the fast path.
         */
        std::vector<unsigned int> renumber;

        /**
         * Check whether the shape functions are linear.
         */
        bool polynomials_are_hat_functions;
      };



      template <int dim, class VectorType>
      RefSpaceFEFieldFunction<dim, VectorType>::RefSpaceFEFieldFunction(
        const DoFHandler<dim> &dof_handler,
        const VectorType &     dof_values)
        : dof_handler(&dof_handler)
        , global_dof_values(&dof_values)
      {
        Assert(dof_handler.n_dofs() == dof_values.size(),
               ExcDimensionMismatch(dof_handler.n_dofs(), dof_values.size()));
      }



      template <int dim, class VectorType>
      void
      RefSpaceFEFieldFunction<dim, VectorType>::set_active_cell(
        const typename Triangulation<dim>::active_cell_iterator &cell)
      {
        Assert(
          &cell->get_triangulation() == &dof_handler->get_triangulation(),
          ExcMessage(
            "The incoming cell must belong to the triangulation associated with "
            "the DoFHandler passed to the constructor."));

        const typename DoFHandler<dim>::active_cell_iterator dof_handler_cell(
          &dof_handler->get_triangulation(),
          cell->level(),
          cell->index(),
          dof_handler);

        // Save the element and the local dof values, since this is what we need
        // to evaluate the function.

        // Check if we can use the fast path. In case we have a different
        // element from the one used before we need to set up the data
        // structures again.
        if (element != &dof_handler_cell->get_fe())
          {
            poly.clear();
            element = &dof_handler_cell->get_fe();

            if (element->n_base_elements() == 1 &&
                dealii::internal::FEPointEvaluation::is_fast_path_supported(
                  *element, 0))
              {
                dealii::internal::MatrixFreeFunctions::ShapeInfo<double>
                  shape_info;

                shape_info.reinit(QMidpoint<1>(), *element, 0);
                renumber = shape_info.lexicographic_numbering;
                poly =
                  dealii::internal::FEPointEvaluation::get_polynomial_space(
                    element->base_element(0));

                polynomials_are_hat_functions =
                  (poly.size() == 2 && poly[0].value(0.) == 1. &&
                   poly[0].value(1.) == 0. && poly[1].value(0.) == 0. &&
                   poly[1].value(1.) == 1.);
              }
          }
        else
          element = &dof_handler_cell->get_fe();

        local_dof_indices.resize(element->dofs_per_cell);
        dof_handler_cell->get_dof_indices(local_dof_indices);

        local_dof_values.resize(element->dofs_per_cell);

        for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
          local_dof_values[i] =
            dealii::internal::ElementAccess<VectorType>::get(
              *global_dof_values, local_dof_indices[i]);
      }



      template <int dim, class VectorType>
      bool
      RefSpaceFEFieldFunction<dim, VectorType>::cell_is_set() const
      {
        // If set cell hasn't been called the size of local_dof_values will be
        // zero.
        return local_dof_values.size() > 0;
      }



      template <int dim, class VectorType>
      double
      RefSpaceFEFieldFunction<dim, VectorType>::value(
        const Point<dim> & point,
        const unsigned int component) const
      {
        AssertIndexRange(component, this->n_components);
        Assert(cell_is_set(), ExcCellNotSet());

        if (!poly.empty() && component == 0)
          {
            // TODO: this could be extended to a component that is not zero
            return dealii::internal::evaluate_tensor_product_value_and_gradient(
                     poly,
                     local_dof_values,
                     point,
                     polynomials_are_hat_functions,
                     renumber)
              .first;
          }
        else
          {
            double value = 0;
            for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
              value += local_dof_values[i] *
                       element->shape_value_component(i, point, component);

            return value;
          }
      }



      template <int dim, class VectorType>
      Tensor<1, dim>
      RefSpaceFEFieldFunction<dim, VectorType>::gradient(
        const Point<dim> & point,
        const unsigned int component) const
      {
        AssertIndexRange(component, this->n_components);
        Assert(cell_is_set(), ExcCellNotSet());

        if (!poly.empty() && component == 0)
          {
            // TODO: this could be extended to a component that is not zero
            return dealii::internal::evaluate_tensor_product_value_and_gradient(
                     poly,
                     local_dof_values,
                     point,
                     polynomials_are_hat_functions,
                     renumber)
              .second;
          }
        else
          {
            Tensor<1, dim> gradient;
            for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
              gradient += local_dof_values[i] *
                          element->shape_grad_component(i, point, component);

            return gradient;
          }
      }



      template <int dim, class VectorType>
      SymmetricTensor<2, dim>
      RefSpaceFEFieldFunction<dim, VectorType>::hessian(
        const Point<dim> & point,
        const unsigned int component) const
      {
        AssertIndexRange(component, this->n_components);
        Assert(cell_is_set(), ExcCellNotSet());

        if (!poly.empty() && component == 0)
          {
            // TODO: this could be extended to a component that is not zero
            return dealii::internal::evaluate_tensor_product_hessian(
              poly, local_dof_values, point, renumber);
          }
        else
          {
            Tensor<2, dim> hessian;
            for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
              hessian +=
                local_dof_values[i] *
                element->shape_grad_grad_component(i, point, component);

            return symmetrize(hessian);
          }
      }
    } // namespace FEValuesImplementation
  }   // namespace internal



  template <int dim>
  template <class VectorType>
  FEValues<dim>::FEValues(const hp::FECollection<dim> &fe_collection,
                          const Quadrature<1> &        quadrature,
                          const RegionUpdateFlags      region_update_flags,
                          const MeshClassifier<dim> &  mesh_classifier,
                          const DoFHandler<dim> &      dof_handler,
                          const VectorType &           level_set,
                          const AdditionalData &       additional_data)
    : mapping_collection(&dealii::hp::StaticMappingQ1<dim>::mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(quadrature)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , quadrature_generator(q_collection_1D, additional_data)
    , reference_space_level_set(
        std::make_unique<internal::FEValuesImplementation::
                           RefSpaceFEFieldFunction<dim, VectorType>>(
          dof_handler,
          level_set))
  {
    // Tensor products of each quadrature in q_collection_1D. Used on the
    // non-intersected cells.
    hp::QCollection<dim> q_collection;
    for (unsigned int i = 0; i < q_collection_1D.size(); ++i)
      q_collection.push_back(Quadrature<dim>(q_collection_1D[i]));

    initialize(q_collection);
  }



  template <int dim>
  template <class VectorType>
  FEValues<dim>::FEValues(const hp::MappingCollection<dim> &mapping_collection,
                          const hp::FECollection<dim> &     fe_collection,
                          const hp::QCollection<dim> &      q_collection,
                          const hp::QCollection<1> &        q_collection_1D,
                          const RegionUpdateFlags           region_update_flags,
                          const MeshClassifier<dim> &       mesh_classifier,
                          const DoFHandler<dim> &           dof_handler,
                          const VectorType &                level_set,
                          const AdditionalData &            additional_data)
    : mapping_collection(&mapping_collection)
    , fe_collection(&fe_collection)
    , q_collection_1D(q_collection_1D)
    , region_update_flags(region_update_flags)
    , mesh_classifier(&mesh_classifier)
    , quadrature_generator(q_collection_1D, additional_data)
    , reference_space_level_set(
        std::make_unique<internal::FEValuesImplementation::
                           RefSpaceFEFieldFunction<dim, VectorType>>(
          dof_handler,
          level_set))
  {
    initialize(q_collection);
  }



  template <int dim>
  void
  FEValues<dim>::initialize(const hp::QCollection<dim> &q_collection)
  {
    current_cell_location = LocationToLevelSet::unassigned;
    active_fe_index       = numbers::invalid_unsigned_int;

    Assert(fe_collection->size() > 0,
           ExcMessage("Incoming hp::FECollection can not be empty."));
    Assert(mapping_collection->size() == fe_collection->size() ||
             mapping_collection->size() == 1,
           ExcMessage("Size of hp::MappingCollection must be "
                      "the same as hp::FECollection or 1."));
    Assert(q_collection.size() == fe_collection->size() ||
             q_collection.size() == 1,
           ExcMessage("Size of hp::QCollection<dim> must be the "
                      "same as hp::FECollection or 1."));
    Assert(q_collection_1D.size() == fe_collection->size() ||
             q_collection_1D.size() == 1,
           ExcMessage("Size of hp::QCollection<1> must be the "
                      "same as hp::FECollection or 1."));

    // For each element in fe_collection, create dealii::FEValues objects to use
    // on the non-intersected cells.
    fe_values_inside_full_quadrature.resize(fe_collection->size());
    fe_values_outside_full_quadrature.resize(fe_collection->size());
    for (unsigned int fe_index = 0; fe_index < fe_collection->size();
         ++fe_index)
      {
        const unsigned int mapping_index =
          mapping_collection->size() > 1 ? fe_index : 0;
        const unsigned int q_index = q_collection.size() > 1 ? fe_index : 0;

        fe_values_inside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.inside);
        fe_values_outside_full_quadrature[fe_index].emplace(
          (*mapping_collection)[mapping_index],
          (*fe_collection)[fe_index],
          q_collection[q_index],
          region_update_flags.outside);
      }
  }



  template <int dim>
  template <bool level_dof_access>
  void
  FEValues<dim>::reinit(
    const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell)
  {
    current_cell_location = mesh_classifier->location_to_level_set(cell);
    active_fe_index       = cell->active_fe_index();

    // These objects were created with a quadrature based on the previous cell
    // and are thus no longer valid.
    fe_values_inside.reset();
    fe_values_surface.reset();
    fe_values_outside.reset();

    switch (current_cell_location)
      {
        case LocationToLevelSet::inside:
          {
            fe_values_inside_full_quadrature.at(active_fe_index)->reinit(cell);
            break;
          }
        case LocationToLevelSet::outside:
          {
            fe_values_outside_full_quadrature.at(active_fe_index)->reinit(cell);
            break;
          }
        case LocationToLevelSet::intersected:
          {
            const unsigned int mapping_index =
              mapping_collection->size() > 1 ? active_fe_index : 0;

            const unsigned int q1D_index =
              q_collection_1D.size() > 1 ? active_fe_index : 0;
            quadrature_generator.set_1D_quadrature(q1D_index);

            reference_space_level_set->set_active_cell(cell);
            const BoundingBox<dim> unit_box = create_unit_bounding_box<dim>();
            quadrature_generator.generate(*reference_space_level_set, unit_box);

            const Quadrature<dim> &inside_quadrature =
              quadrature_generator.get_inside_quadrature();
            const Quadrature<dim> &outside_quadrature =
              quadrature_generator.get_outside_quadrature();
            const ImmersedSurfaceQuadrature<dim> &surface_quadrature =
              quadrature_generator.get_surface_quadrature();

            // Even if a cell is formally intersected the number of created
            // quadrature points can be 0. Avoid creating an FEValues object
            // if that is the case.
            if (inside_quadrature.size() > 0)
              {
                fe_values_inside.emplace((*mapping_collection)[mapping_index],
                                         (*fe_collection)[active_fe_index],
                                         inside_quadrature,
                                         region_update_flags.inside);

                fe_values_inside->reinit(cell);
              }

            if (outside_quadrature.size() > 0)
              {
                fe_values_outside.emplace((*mapping_collection)[mapping_index],
                                          (*fe_collection)[active_fe_index],
                                          outside_quadrature,
                                          region_update_flags.outside);

                fe_values_outside->reinit(cell);
              }

            if (surface_quadrature.size() > 0)
              {
                fe_values_surface.emplace((*mapping_collection)[mapping_index],
                                          (*fe_collection)[active_fe_index],
                                          surface_quadrature,
                                          region_update_flags.surface);
                fe_values_surface->reinit(cell);
              }

            break;
          }
        default:
          {
            Assert(false, ExcInternalError());
            break;
          }
      }
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEValues<dim>> &
  FEValues<dim>::get_inside_fe_values() const
  {
    if (current_cell_location == LocationToLevelSet::inside)
      return fe_values_inside_full_quadrature.at(active_fe_index);
    else
      return fe_values_inside;
  }



  template <int dim>
  const std_cxx17::optional<dealii::FEValues<dim>> &
  FEValues<dim>::get_outside_fe_values() const
  {
    if (current_cell_location == LocationToLevelSet::outside)
      return fe_values_outside_full_quadrature.at(active_fe_index);
    else
      return fe_values_outside;
  }



  template <int dim>
  const std_cxx17::optional<FEImmersedSurfaceValues<dim>> &
  FEValues<dim>::get_surface_fe_values() const
  {
    return fe_values_surface;
  }


#include "fe_values.inst"

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE
