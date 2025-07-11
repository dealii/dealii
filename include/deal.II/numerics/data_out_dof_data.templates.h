// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_data_out_dof_data_templates_h
#define dealii_data_out_dof_data_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/read_write_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutImplementation
  {
    template <int dim, int spacedim>
    ParallelDataBase<dim, spacedim>::ParallelDataBase(
      const unsigned int                    n_datasets,
      const unsigned int                    n_subdivisions,
      const std::vector<unsigned int>      &n_postprocessor_outputs,
      const dealii::Mapping<dim, spacedim> &mapping,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                       &finite_elements,
      const UpdateFlags update_flags,
      const bool        use_face_values)
      : ParallelDataBase<dim, spacedim>(
          n_datasets,
          n_subdivisions,
          n_postprocessor_outputs,
          dealii::hp::MappingCollection<dim, spacedim>(mapping),
          finite_elements,
          update_flags,
          use_face_values)
    {}



    /**
     * Generate evaluation points on a simplex with arbitrary number of
     * subdivisions.
     */
    template <int dim>
    inline std::vector<Point<dim>>
    generate_simplex_evaluation_points(const unsigned int /*n_subdivisions*/)
    {
      DEAL_II_NOT_IMPLEMENTED();

      return {};
    }



    /**
     * Helper function to create evaluation points recursively with
     * subdivisions=0,1,2 being the base case:
     *                                        +
     *                                        |\
     *                            +           +-+
     *                            |\          |\|\
     *                  +         +-+         +-+-+
     *                  |\        |\|\        |\|\|\
     *          +       +-+       +-+-+       +-+-+-+
     *          |\      |\|\      |\|\|\      |\|\|\|\
     *    +     +-+     +-+-+     +-+-+-+     +-+-+-+-+
     *
     *    0      1        2          3            4
     *    ^      ^                   |            |
     *    |      |                   |            |
     *    +--------------------------+            |
     *           |                                |
     *           +--------------------------------+
     */
    inline void
    generate_simplex_evaluation_points_recursively(
      const std::vector<Point<2>> &bounding_vertices,
      const unsigned int           n_subdivisions,
      std::vector<Point<2>>       &evaluation_points)
    {
      if (n_subdivisions == 0)
        {
          evaluation_points.push_back(bounding_vertices[0]);
          return;
        }

      for (const auto &p : bounding_vertices)
        evaluation_points.push_back(p);

      if (n_subdivisions == 1)
        return;

      // Helper functions to create intermediate points between p0 and p1 with
      // n_subdivisions. This function appends these new points to the vector
      // evaluation_points but also returns the points as a vector to be able
      // to easily access specific points on the line.
      const auto generate_inbetween_points = [&](const Point<2> &p0,
                                                 const Point<2> &p1) {
        std::vector<Point<2>> line_points;

        for (unsigned int i = 1; i < n_subdivisions; ++i)
          line_points.push_back(p0 + (p1 - p0) / n_subdivisions * i);

        evaluation_points.insert(evaluation_points.end(),
                                 line_points.begin(),
                                 line_points.end());

        return line_points;
      };

      const auto line_points_0 =
        generate_inbetween_points(bounding_vertices[0], bounding_vertices[1]);
      const auto line_points_1 =
        generate_inbetween_points(bounding_vertices[1], bounding_vertices[2]);
      const auto line_points_2 =
        generate_inbetween_points(bounding_vertices[2], bounding_vertices[0]);

      if (n_subdivisions == 2)
        return;

      generate_simplex_evaluation_points_recursively(
        // create new inner triangle (see ASCII art above)
        {{Point<2>(line_points_0[line_points_0.size() - 2][0],
                   line_points_1[0][1]),
          Point<2>(line_points_0[0][0],
                   line_points_1[line_points_1.size() - 2][1]),
          Point<2>(line_points_0[0][0], line_points_1[0][1])}},
        n_subdivisions - 3,
        evaluation_points);
    }



    /**
     * Specialization for triangles.
     */
    template <>
    inline std::vector<Point<2>>
    generate_simplex_evaluation_points(const unsigned int n_subdivisions)
    {
      std::vector<Point<2>> evaluation_points;

      generate_simplex_evaluation_points_recursively(
        {{Point<2>(0.0, 0.0), Point<2>(1.0, 0.0), Point<2>(0.0, 1.0)}},
        n_subdivisions,
        evaluation_points);

      return evaluation_points;
    }



    /**
     * Set up vectors of FEValues and FEFaceValues needed inside of
     * ParallelDataBase and return the maximum number of quadrature points
     * needed to allocate enough memory for the scratch data.
     */
    template <int dim, int spacedim>
    unsigned int
    setup_parallel_data_base_internal(
      const unsigned int                                  n_subdivisions,
      const dealii::hp::MappingCollection<dim, spacedim> &mapping_collection,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                       &finite_elements,
      const UpdateFlags update_flags,
      const bool        use_face_values,
      std::vector<std::shared_ptr<dealii::hp::FEValues<dim, spacedim>>>
        &x_fe_values,
      std::vector<std::shared_ptr<dealii::hp::FEFaceValues<dim, spacedim>>>
        &x_fe_face_values)
    {
      // First figure out which reference cell types are present in the
      // FECollections we got as arguments. To this end, use a lambda
      // function that for a FECollection object tests whether one of the
      // elements uses a specific reference cell. Then use this
      // lambda function in a std::any_of() call that accumulates over
      // all of the FECollection objects we were given.
      //
      // Note that in 1d, we will count a line segment as a hypercube
      // even though it is *also* a simplex. Furthermore, wedges and
      // pyramids can only appear in 3d set ups, so we do not need
      // to test there if dim<3.
      static const auto has_fe_with_reference_cell =
        [](const dealii::hp::FECollection<dim, spacedim> &fe_collection,
           const ReferenceCell                           &reference_cell) {
          for (unsigned int i = 0; i < fe_collection.size(); ++i)
            if (fe_collection[i].reference_cell() == reference_cell)
              return true;
          return false;
        };

      const bool needs_hypercube_setup = std::any_of(
        finite_elements.begin(),
        finite_elements.end(),
        [](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
             &fe_collection) {
          return has_fe_with_reference_cell(
            *fe_collection, ReferenceCells::get_hypercube<dim>());
        });
      const bool needs_simplex_setup =
        (dim > 1 &&
         std::any_of(
           finite_elements.begin(),
           finite_elements.end(),
           [](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
                &fe_collection) {
             return has_fe_with_reference_cell(
               *fe_collection, ReferenceCells::get_simplex<dim>());
           }));
      const bool needs_wedge_setup =
        (dim == 3 &&
         std::any_of(
           finite_elements.begin(),
           finite_elements.end(),
           [](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
                &fe_collection) {
             return has_fe_with_reference_cell(*fe_collection,
                                               ReferenceCells::Wedge);
           }));
      const bool needs_pyramid_setup =
        (dim == 3 &&
         std::any_of(
           finite_elements.begin(),
           finite_elements.end(),
           [](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
                &fe_collection) {
             return has_fe_with_reference_cell(*fe_collection,
                                               ReferenceCells::Pyramid);
           }));

      // Decide whether we want to work on cell or face FE(Face)Values objects:
      if (use_face_values == false)
        {
          std::unique_ptr<dealii::Quadrature<dim>> quadrature_simplex;
          std::unique_ptr<dealii::Quadrature<dim>> quadrature_hypercube;
          std::unique_ptr<dealii::Quadrature<dim>> quadrature_wedge;
          std::unique_ptr<dealii::Quadrature<dim>> quadrature_pyramid;

          if (needs_simplex_setup)
            {
              if (dim == 2)
                quadrature_simplex = std::make_unique<Quadrature<dim>>(
                  generate_simplex_evaluation_points<dim>(n_subdivisions));
              else
                quadrature_simplex = std::make_unique<Quadrature<dim>>(
                  FE_SimplexP<dim>(n_subdivisions).get_unit_support_points());
            }

          if (needs_hypercube_setup)
            {
              quadrature_hypercube =
                std::make_unique<QIterated<dim>>(QTrapezoid<1>(),
                                                 n_subdivisions);
            }

          if (needs_wedge_setup)
            {
              Assert(n_subdivisions == 1, ExcNotImplemented());

              quadrature_wedge = std::make_unique<Quadrature<dim>>(
                ReferenceCells::Wedge.get_nodal_type_quadrature<dim>());
            }

          if (needs_pyramid_setup)
            {
              Assert(n_subdivisions == 1, ExcNotImplemented());

              quadrature_pyramid = std::make_unique<Quadrature<dim>>(
                ReferenceCells::Pyramid.get_nodal_type_quadrature<dim>());
            }

          x_fe_values.resize(finite_elements.size());
          for (unsigned int i = 0; i < finite_elements.size(); ++i)
            {
              // Check if one of the previous finite elements is equal to the
              // present one. If so, re-use the FEValues object.
              for (unsigned int j = 0; j < i; ++j)
                if (finite_elements[i].get() == finite_elements[j].get())
                  {
                    x_fe_values[i] = x_fe_values[j];
                    break;
                  }

              // If none was found, create an FEValues object:
              if (x_fe_values[i].get() == nullptr)
                {
                  dealii::hp::QCollection<dim> quadrature;

                  for (unsigned int j = 0; j < finite_elements[i]->size(); ++j)
                    {
                      const ReferenceCell reference_cell =
                        (*finite_elements[i])[j].reference_cell();

                      if (reference_cell.is_hyper_cube())
                        quadrature.push_back(*quadrature_hypercube);
                      else if (reference_cell.is_simplex())
                        quadrature.push_back(*quadrature_simplex);
                      else if (reference_cell == ReferenceCells::Wedge)
                        quadrature.push_back(*quadrature_wedge);
                      else if (reference_cell == ReferenceCells::Pyramid)
                        quadrature.push_back(*quadrature_pyramid);
                      else
                        DEAL_II_NOT_IMPLEMENTED();
                    }

                  x_fe_values[i] =
                    std::make_shared<dealii::hp::FEValues<dim, spacedim>>(
                      mapping_collection,
                      *finite_elements[i],
                      quadrature,
                      // certain reference-cell kinds take by default the
                      // curved code path in DataOut::build_one_patch, for
                      // which the quadrature points need to be enabled
                      (needs_simplex_setup || needs_wedge_setup ||
                       needs_pyramid_setup) ?
                        (update_flags | update_quadrature_points) :
                        update_flags);
                }
            }

          // Return maximal number of evaluation points:
          return std::max(
            {needs_wedge_setup ? quadrature_wedge->size() : 0,
             needs_simplex_setup ? quadrature_simplex->size() : 0,
             needs_hypercube_setup ? quadrature_hypercube->size() : 0,
             needs_pyramid_setup ? quadrature_pyramid->size() : 0});
        }
      else // build FEFaceValues objects instead
        {
          // The code following is not quite right for wedges and pyramids.
          // Assert that we don't have these kinds of meshes.
          Assert(needs_pyramid_setup == false && needs_wedge_setup == false,
                 ExcNotImplemented());

          std::unique_ptr<dealii::Quadrature<dim - 1>> quadrature_simplex;
          std::unique_ptr<dealii::Quadrature<dim - 1>> quadrature_hypercube;

          // See whether we need simplex or hypercube quadrature formulas.
          // This is only an issue in 3d (in 2d every face integral is just
          // a line integral, so we can deal with that via the usual hypercube
          // quadrature rule).
          if ((dim == 3) &&
              (needs_simplex_setup || needs_pyramid_setup || needs_wedge_setup))
            {
              quadrature_simplex = std::make_unique<Quadrature<dim - 1>>(
                generate_simplex_evaluation_points<dim - 1>(n_subdivisions));
            }

          if ((dim < 3) || (needs_hypercube_setup || needs_pyramid_setup ||
                            needs_wedge_setup))
            {
              quadrature_hypercube =
                std::make_unique<QIterated<dim - 1>>(QTrapezoid<1>(),
                                                     n_subdivisions);
            }

          x_fe_face_values.resize(finite_elements.size());
          for (unsigned int i = 0; i < finite_elements.size(); ++i)
            {
              // Check if one of the previous finite elements is equal to the
              // present one. If so, re-use the FEValues object.
              for (unsigned int j = 0; j < i; ++j)
                if (finite_elements[i].get() == finite_elements[j].get())
                  {
                    x_fe_face_values[i] = x_fe_face_values[j];
                    break;
                  }

              // If none was found, create an FEFaceValues object:
              if (x_fe_face_values[i].get() == nullptr)
                {
                  dealii::hp::QCollection<dim - 1> quadrature;

                  for (unsigned int j = 0; j < finite_elements[i]->size(); ++j)
                    {
                      const ReferenceCell reference_cell =
                        (*finite_elements[i])[j].reference_cell();

                      // In 1d/2d and for hypercube/wedge/pyramid elements, we
                      // need hypercube quadratures.
                      if ((dim < 3) ||
                          (reference_cell.is_hyper_cube() ||
                           (reference_cell == ReferenceCells::Wedge) ||
                           (reference_cell == ReferenceCells::Pyramid)))
                        quadrature.push_back(*quadrature_hypercube);

                      // In 3d, if the element is for simplex/wedge/pyramid
                      // cells, then we also need simplex quadratures
                      if ((dim == 3) &&
                          (reference_cell.is_simplex() ||
                           (reference_cell == ReferenceCells::Wedge) ||
                           (reference_cell == ReferenceCells::Pyramid)))
                        quadrature.push_back(*quadrature_simplex);
                    }

                  x_fe_face_values[i] =
                    std::make_shared<dealii::hp::FEFaceValues<dim, spacedim>>(
                      mapping_collection,
                      *finite_elements[i],
                      quadrature,
                      update_flags);
                }
            }

          // Return maximal number of evaluation points:
          return std::max(
            {(dim == 3) && (needs_simplex_setup || needs_pyramid_setup ||
                            needs_wedge_setup) ?
               quadrature_simplex->size() :
               0,
             (dim < 3) || needs_hypercube_setup ? quadrature_hypercube->size() :
                                                  0});
        }
    }



    template <int dim, int spacedim>
    ParallelDataBase<dim, spacedim>::ParallelDataBase(
      const unsigned int               n_datasets,
      const unsigned int               n_subdivisions,
      const std::vector<unsigned int> &n_postprocessor_outputs,
      const dealii::hp::MappingCollection<dim, spacedim> &mapping,
      const std::vector<
        std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                       &finite_elements,
      const UpdateFlags update_flags,
      const bool        use_face_values)
      : n_datasets(n_datasets)
      , n_subdivisions(n_subdivisions)
      , postprocessed_values(n_postprocessor_outputs.size())
      , mapping_collection(mapping)
      , finite_elements(finite_elements)
      , update_flags(update_flags)
    {
      const unsigned int n_q_points =
        setup_parallel_data_base_internal(n_subdivisions,
                                          mapping,
                                          finite_elements,
                                          update_flags,
                                          use_face_values,
                                          x_fe_values,
                                          x_fe_face_values);

      patch_values_scalar.solution_values.resize(n_q_points);
      patch_values_scalar.solution_gradients.resize(n_q_points);
      patch_values_scalar.solution_hessians.resize(n_q_points);
      patch_values_system.solution_values.resize(n_q_points);
      patch_values_system.solution_gradients.resize(n_q_points);
      patch_values_system.solution_hessians.resize(n_q_points);

      for (unsigned int dataset = 0; dataset < n_postprocessor_outputs.size();
           ++dataset)
        if (n_postprocessor_outputs[dataset] != 0)
          postprocessed_values[dataset].resize(
            n_q_points,
            dealii::Vector<double>(n_postprocessor_outputs[dataset]));
    }



    // implement copy constructor to create a thread's own version of
    // x_fe_values
    template <int dim, int spacedim>
    ParallelDataBase<dim, spacedim>::ParallelDataBase(
      const ParallelDataBase<dim, spacedim> &data)
      : n_datasets(data.n_datasets)
      , n_subdivisions(data.n_subdivisions)
      , patch_values_scalar(data.patch_values_scalar)
      , patch_values_system(data.patch_values_system)
      , postprocessed_values(data.postprocessed_values)
      , mapping_collection(data.mapping_collection)
      , finite_elements(data.finite_elements)
      , update_flags(data.update_flags)
    {
      setup_parallel_data_base_internal(n_subdivisions,
                                        mapping_collection,
                                        finite_elements,
                                        update_flags,
                                        data.x_fe_values.empty(),
                                        x_fe_values,
                                        x_fe_face_values);
    }



    template <int dim, int spacedim>
    void
    ParallelDataBase<dim, spacedim>::reinit_all_fe_values(
      std::vector<std::shared_ptr<DataEntryBase<dim, spacedim>>> &dof_data,
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const unsigned int                                                  face)
    {
      for (unsigned int dataset = 0; dataset < dof_data.size(); ++dataset)
        {
          const bool is_duplicate = std::any_of(
            finite_elements.cbegin(),
            finite_elements.cbegin() + dataset,
            [&](const std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>
                  &fe) { return finite_elements[dataset].get() == fe.get(); });
          if (is_duplicate == false)
            {
              if (cell->is_active())
                {
                  const typename DoFHandler<dim, spacedim>::active_cell_iterator
                    dh_cell(&cell->get_triangulation(),
                            cell->level(),
                            cell->index(),
                            dof_data[dataset]->dof_handler);

                  // Check whether we need cell or face FEValues objects by
                  // testing which of the two arrays actually has any content.
                  if (x_fe_values.empty())
                    {
                      AssertIndexRange(face, GeometryInfo<dim>::faces_per_cell);
                      x_fe_face_values[dataset]->reinit(dh_cell, face);
                    }
                  else
                    x_fe_values[dataset]->reinit(dh_cell);
                }
              else
                x_fe_values[dataset]->reinit(cell);
            }
        }

      // If there is are no DoF-associated data (just cell-associated ones),
      // then the loop above will not execute any iterations. In that case,
      // do the initialization for the first FE(Face)Values object by
      // hand, using only the (triangulation) cell without a DoFHandler.
      if (dof_data.empty())
        {
          if (x_fe_values.empty())
            {
              AssertIndexRange(face, GeometryInfo<dim>::faces_per_cell);
              x_fe_face_values[0]->reinit(cell, face);
            }
          else
            x_fe_values[0]->reinit(cell);
        }
    }



    template <int dim, int spacedim>
    const FEValuesBase<dim, spacedim> &
    ParallelDataBase<dim, spacedim>::get_present_fe_values(
      const unsigned int dataset) const
    {
      AssertIndexRange(dataset, finite_elements.size());

      // Check whether we need cell or face FEValues objects by testing
      // which of the two arrays actually has any content.
      if (x_fe_values.empty())
        return x_fe_face_values[dataset]->get_present_fe_values();
      else
        return x_fe_values[dataset]->get_present_fe_values();
    }



    template <int dim, int spacedim>
    void
    ParallelDataBase<dim, spacedim>::resize_system_vectors(
      const unsigned int n_components)
    {
      Assert(patch_values_system.solution_values.size() > 0,
             ExcInternalError());
      AssertDimension(patch_values_system.solution_values.size(),
                      patch_values_system.solution_gradients.size());
      AssertDimension(patch_values_system.solution_values.size(),
                      patch_values_system.solution_hessians.size());
      for (unsigned int k = 0; k < patch_values_system.solution_values.size();
           ++k)
        {
          patch_values_system.solution_values[k].reinit(n_components);
          patch_values_system.solution_gradients[k].resize(n_components);
          patch_values_system.solution_hessians[k].resize(n_components);
        }
    }



    /**
     * In a WorkStream context, use this function to append the patch computed
     * by the parallel stage to the array of patches.
     */
    template <int dim, int spacedim>
    void
    append_patch_to_list(
      const DataOutBase::Patch<dim, spacedim>        &patch,
      std::vector<DataOutBase::Patch<dim, spacedim>> &patches)
    {
      patches.push_back(patch);
      patches.back().patch_index = patches.size() - 1;
    }
  } // namespace DataOutImplementation
} // namespace internal

namespace internal
{
  namespace DataOutImplementation
  {
    /**
     * Extract the specified component of a number. This template is used when
     * the given value is assumed to be a real scalar, so asking for the real
     * part is the only valid choice for the second argument.
     */
    template <typename NumberType>
    double
    get_component(const NumberType         value,
                  const ComponentExtractor extract_component)
    {
      (void)extract_component;
      static_assert(
        numbers::NumberTraits<NumberType>::is_complex == false,
        "This function must not be called for complex-valued data types.");
      Assert(extract_component == ComponentExtractor::real_part,
             ExcMessage("You cannot extract anything other than the real "
                        "part from a real number."));
      return value;
    }



    /**
     * Extract the specified component of a number. This template is used when
     * the given value is a complex number
     */
    template <typename NumberType>
    double
    get_component(const std::complex<NumberType> &value,
                  const ComponentExtractor        extract_component)
    {
      switch (extract_component)
        {
          case ComponentExtractor::real_part:
            return value.real();

          case ComponentExtractor::imaginary_part:
            return value.imag();

          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }

      return numbers::signaling_nan<double>();
    }



    template <int rank, int dim, typename NumberType>
    Tensor<rank, dim>
    get_component(const Tensor<rank, dim, NumberType> &value,
                  const ComponentExtractor             extract_component)
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcMessage("You cannot extract anything other than the real "
                        "part from a real number."));

      Tensor<rank, dim, double> t;
      for (unsigned int d = 0; d < dim; ++d)
        t[d] = get_component(value[d], extract_component);

      return t;
    }



    template <int dim, int spacedim>
    DataEntryBase<dim, spacedim>::DataEntryBase(
      const DoFHandler<dim, spacedim> *dofs,
      const std::vector<std::string>  &names_in,
      const std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
        &data_component_interpretation)
      : dof_handler(
          dofs,
          typeid(dealii::DataOut_DoFData<dim, dim, spacedim, spacedim>).name())
      , names(names_in)
      , data_component_interpretation(data_component_interpretation)
      , postprocessor(nullptr, typeid(*this).name())
      , n_output_variables(names.size())
    {
      Assert(names.size() == data_component_interpretation.size(),
             ExcDimensionMismatch(data_component_interpretation.size(),
                                  names.size()));

      // check that the names use only allowed characters
      for (const auto &name : names)
        {
          (void)name;
          Assert(name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                        "0123456789_<>()") == std::string::npos,
                 Exceptions::DataOutImplementation::ExcInvalidCharacter(
                   name,
                   name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                          "0123456789_<>()")));
        }
    }



    template <int dim, int spacedim>
    DataEntryBase<dim, spacedim>::DataEntryBase(
      const DoFHandler<dim, spacedim>   *dofs,
      const DataPostprocessor<spacedim> *data_postprocessor)
      : dof_handler(
          dofs,
          typeid(dealii::DataOut_DoFData<dim, dim, spacedim, spacedim>).name())
      , names(data_postprocessor->get_names())
      , data_component_interpretation(
          data_postprocessor->get_data_component_interpretation())
      , postprocessor(data_postprocessor, typeid(*this).name())
      , n_output_variables(names.size())
    {
      Assert(data_postprocessor->get_names().size() ==
               data_postprocessor->get_data_component_interpretation().size(),
             ExcDimensionMismatch(
               data_postprocessor->get_names().size(),
               data_postprocessor->get_data_component_interpretation().size()));

      // check that the names use only allowed characters
      for (const auto &name : names)
        {
          (void)name;
          Assert(name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                        "0123456789_<>()") == std::string::npos,
                 Exceptions::DataOutImplementation::ExcInvalidCharacter(
                   name,
                   name.find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                          "0123456789_<>()")));
        }
    }



    namespace CreateVectors
    {
      /**
       * Copy the data from an arbitrary non-block vector to a
       * LinearAlgebra::distributed::Vector.
       */
      template <typename VectorType, typename Number>
      void
      copy_locally_owned_data_from(
        const VectorType                           &src,
        LinearAlgebra::distributed::Vector<Number> &dst)
      {
        // If source and destination vector have the same underlying scalar,
        // we can directly import elements by using only one temporary vector:
        if constexpr (std::is_same_v<typename VectorType::value_type, Number>)
          {
            LinearAlgebra::ReadWriteVector<typename VectorType::value_type>
              temp;
            temp.reinit(src.locally_owned_elements());
            temp.import_elements(src, VectorOperation::insert);

            dst.import_elements(temp, VectorOperation::insert);
          }
        else
          // The source and destination vector have different scalar types. We
          // need to split the parallel import and local copy operations into
          // two phases
          {
            LinearAlgebra::ReadWriteVector<typename VectorType::value_type>
              temp;
            temp.reinit(src.locally_owned_elements());
            temp.import_elements(src, VectorOperation::insert);

            LinearAlgebra::ReadWriteVector<Number> temp2;
            temp2.reinit(temp, true);
            temp2 = temp;

            dst.import_elements(temp2, VectorOperation::insert);
          }
      }

#ifdef DEAL_II_WITH_TRILINOS
      template <typename Number>
      void
      copy_locally_owned_data_from(
        const TrilinosWrappers::MPI::Vector        &src,
        LinearAlgebra::distributed::Vector<Number> &dst)
      {
        // ReadWriteVector does not work for ghosted
        // TrilinosWrappers::MPI::Vector objects. Fall back to copy the
        // entries manually.
        for (const auto i : dst.locally_owned_elements())
          dst[i] = src[i];
      }
#endif

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
      template <typename Number, typename MemorySpace>
      void
      copy_locally_owned_data_from(
        const LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &src,
        LinearAlgebra::distributed::Vector<Number>                       &dst)
      {
        // ReadWriteVector does not work for ghosted
        // TrilinosWrappers::MPI::Vector objects. Fall back to copy the
        // entries manually.
        for (const auto i : dst.locally_owned_elements())
          dst[i] = src[i];
      }
#endif

      /**
       * Create a ghosted-copy of a block dof vector.
       */
      template <int dim,
                int spacedim,
                typename VectorType,
                typename Number,
                std::enable_if_t<IsBlockVector<VectorType>::value, VectorType>
                  * = nullptr>
      void
      create_dof_vector(
        const DoFHandler<dim, spacedim>                 &dof_handler,
        const VectorType                                &src,
        LinearAlgebra::distributed::BlockVector<Number> &dst,
        const unsigned int level = numbers::invalid_unsigned_int)
      {
        const IndexSet &locally_owned_dofs =
          (level == numbers::invalid_unsigned_int) ?
            dof_handler.locally_owned_dofs() :
            dof_handler.locally_owned_mg_dofs(level);

        const IndexSet locally_relevant_dofs =
          (level == numbers::invalid_unsigned_int) ?
            DoFTools::extract_locally_relevant_dofs(dof_handler) :
            DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

        std::vector<types::global_dof_index> n_indices_per_block(
          src.n_blocks());

        for (unsigned int b = 0; b < src.n_blocks(); ++b)
          n_indices_per_block[b] = src.get_block_indices().block_size(b);

        const auto locally_owned_dofs_b =
          locally_owned_dofs.split_by_block(n_indices_per_block);
        const auto locally_relevant_dofs_b =
          locally_relevant_dofs.split_by_block(n_indices_per_block);

        dst.reinit(src.n_blocks());

        for (unsigned int b = 0; b < src.n_blocks(); ++b)
          {
            Assert(src.block(b).locally_owned_elements().is_contiguous(),
                   ExcMessage(
                     "Using non-contiguous vector blocks is not currently "
                     "supported by DataOut and related classes. The typical "
                     "way you may end up with such vectors is if you order "
                     "degrees of freedom via DoFRenumber::component_wise() "
                     "but then group several vector components into one "
                     "vector block -- for example, all 'dim' components "
                     "of a velocity are grouped into the same vector block. "
                     "If you do this, you don't want to renumber degrees "
                     "of freedom based on vector component, but instead "
                     "based on the 'blocks' you will later use for grouping "
                     "things into vectors. Take a look at step-32 and step-55 "
                     "and how they use the second argument of "
                     "DoFRenumber::component_wise()."));

            dst.block(b).reinit(locally_owned_dofs_b[b],
                                locally_relevant_dofs_b[b],
                                dof_handler.get_mpi_communicator());
            copy_locally_owned_data_from(src.block(b), dst.block(b));
          }

        dst.collect_sizes();

        dst.update_ghost_values();
      }

      /**
       * Create a ghosted-copy of a non-block dof vector.
       */
      template <int dim,
                int spacedim,
                typename VectorType,
                typename Number,
                std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                  * = nullptr>
      void
      create_dof_vector(
        const DoFHandler<dim, spacedim>                 &dof_handler,
        const VectorType                                &src,
        LinearAlgebra::distributed::BlockVector<Number> &dst,
        const unsigned int level = numbers::invalid_unsigned_int)
      {
        const IndexSet &locally_owned_dofs =
          (level == numbers::invalid_unsigned_int) ?
            dof_handler.locally_owned_dofs() :
            dof_handler.locally_owned_mg_dofs(level);

        const IndexSet locally_relevant_dofs =
          (level == numbers::invalid_unsigned_int) ?
            DoFTools::extract_locally_relevant_dofs(dof_handler) :
            DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);


        Assert(locally_owned_dofs.is_contiguous(),
               ExcMessage(
                 "You are trying to add a non-block vector with non-contiguous "
                 "locally-owned index sets. This is not possible. Please "
                 "consider to use an adequate block vector!"));

        dst.reinit(1);

        dst.block(0).reinit(locally_owned_dofs,
                            locally_relevant_dofs,
                            dof_handler.get_mpi_communicator());
        copy_locally_owned_data_from(src, dst.block(0));

        dst.collect_sizes();

        dst.update_ghost_values();
      }

      /**
       * Create a ghosted-copy of a block cell vector.
       */
      template <typename VectorType,
                typename Number,
                std::enable_if_t<IsBlockVector<VectorType>::value, VectorType>
                  * = nullptr>
      void
      create_cell_vector(const VectorType                                &src,
                         LinearAlgebra::distributed::BlockVector<Number> &dst)
      {
        dst.reinit(src.n_blocks());

        for (unsigned int b = 0; b < src.n_blocks(); ++b)
          {
            dst.block(b).reinit(src.get_block_indices().block_size(b));
            copy_locally_owned_data_from(src.block(b), dst.block(b));
          }

        dst.collect_sizes();
      }


      /**
       * Create a ghosted-copy of a non-block cell vector.
       */
      template <typename VectorType,
                typename Number,
                std::enable_if_t<!IsBlockVector<VectorType>::value, VectorType>
                  * = nullptr>
      void
      create_cell_vector(const VectorType                                &src,
                         LinearAlgebra::distributed::BlockVector<Number> &dst)
      {
        dst.reinit(1);

        dst.block(0).reinit(src.size());
        copy_locally_owned_data_from(src, dst.block(0));

        dst.collect_sizes();

        dst.update_ghost_values();
      }
    } // namespace CreateVectors



    /**
     * Class that stores a pointer to a vector of type equal to the template
     * argument, and provides the functions to extract data from it.
     */
    template <int dim, int spacedim, typename ScalarType>
    class DataEntry : public DataEntryBase<dim, spacedim>
    {
    public:
      /**
       * Constructor. Give a list of names for the individual components of
       * the vector and their interpretation as scalar or vector data. This
       * constructor assumes that no postprocessor is going to be used.
       */
      template <typename DataVectorType, typename VectorType>
      DataEntry(const DoFHandler<dim, spacedim> *dofs,
                const VectorType                *data,
                const std::vector<std::string>  &names,
                const std::vector<
                  DataComponentInterpretation::DataComponentInterpretation>
                                    &data_component_interpretation,
                const DataVectorType actual_type);

      /**
       * Constructor when a data postprocessor is going to be used. In that
       * case, the names and vector declarations are going to be acquired from
       * the postprocessor.
       */
      template <typename VectorType>
      DataEntry(const DoFHandler<dim, spacedim>   *dofs,
                const VectorType                  *data,
                const DataPostprocessor<spacedim> *data_postprocessor);

      /**
       * Assuming that the stored vector is a cell vector, extract the given
       * element from it.
       */
      virtual double
      get_cell_data_value(
        const unsigned int       cell_number,
        const ComponentExtractor extract_component) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store.
       */
      virtual void
      get_function_values(const FEValuesBase<dim, spacedim> &fe_patch_values,
                          const ComponentExtractor           extract_component,
                          std::vector<double> &patch_values) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual void
      get_function_values(const FEValuesBase<dim, spacedim> &fe_patch_values,
                          const ComponentExtractor           extract_component,
                          std::vector<dealii::Vector<double>>
                            &patch_values_system) const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<dim, spacedim> &fe_patch_values,
        const ComponentExtractor           extract_component,
        std::vector<Tensor<1, spacedim>>  &patch_gradients) const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_gradients(const FEValuesBase<dim, spacedim> &fe_patch_values,
                             const ComponentExtractor extract_component,
                             std::vector<std::vector<Tensor<1, spacedim>>>
                               &patch_gradients_system) const override;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<dim, spacedim> &fe_patch_values,
        const ComponentExtractor           extract_component,
        std::vector<Tensor<2, spacedim>>  &patch_hessians) const override;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_hessians(const FEValuesBase<dim, spacedim> &fe_patch_values,
                            const ComponentExtractor extract_component,
                            std::vector<std::vector<Tensor<2, spacedim>>>
                              &patch_hessians_system) const override;

      /**
       * Return whether the data represented by (a derived class of) this object
       * represents a complex-valued (as opposed to real-valued) information.
       */
      virtual bool
      is_complex_valued() const override;

      /**
       * Clear all references to the vectors.
       */
      virtual void
      clear() override;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * Copy of the data of the vector passed to DataOutDoFData::add_vector().
       * This vector imports all elements necessary to create output from the
       * source vector and stores it until we no longer need it. No reference
       * to the original source vector is necessary nor stored.
       */
      LinearAlgebra::distributed::BlockVector<ScalarType> vector;
    };



    template <int dim, int spacedim, typename ScalarType>
    template <typename DataVectorType, typename VectorType>
    DataEntry<dim, spacedim, ScalarType>::DataEntry(
      const DoFHandler<dim, spacedim> *dofs,
      const VectorType                *data,
      const std::vector<std::string>  &names,
      const std::vector<
        DataComponentInterpretation::DataComponentInterpretation>
                          &data_component_interpretation,
      const DataVectorType actual_type)
      : DataEntryBase<dim, spacedim>(dofs, names, data_component_interpretation)
    {
      if (actual_type == DataVectorType::type_dof_data)
        CreateVectors::create_dof_vector(*dofs, *data, vector);
      else if (actual_type == DataVectorType::type_cell_data)
        CreateVectors::create_cell_vector(*data, vector);
      else
        DEAL_II_ASSERT_UNREACHABLE();
    }



    template <int dim, int spacedim, typename ScalarType>
    template <typename VectorType>
    DataEntry<dim, spacedim, ScalarType>::DataEntry(
      const DoFHandler<dim, spacedim>   *dofs,
      const VectorType                  *data,
      const DataPostprocessor<spacedim> *data_postprocessor)
      : DataEntryBase<dim, spacedim>(dofs, data_postprocessor)
    {
      CreateVectors::create_dof_vector(*dofs, *data, vector);
    }



    template <int dim, int spacedim, typename ScalarType>
    double
    DataEntry<dim, spacedim, ScalarType>::get_cell_data_value(
      const unsigned int       cell_number,
      const ComponentExtractor extract_component) const
    {
      return get_component(
        internal::ElementAccess<LinearAlgebra::distributed::BlockVector<
          ScalarType>>::get(vector, cell_number),
        extract_component);
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_values(
      const FEValuesBase<dim, spacedim>   &fe_patch_values,
      const ComponentExtractor             extract_component,
      std::vector<dealii::Vector<double>> &patch_values_system) const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));
          fe_patch_values.get_function_values(vector, patch_values_system);
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<dealii::Vector<ScalarType>> tmp(n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            tmp[i].reinit(n_components);

          fe_patch_values.get_function_values(vector, tmp);

          AssertDimension(patch_values_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            {
              AssertDimension(patch_values_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; ++j)
                patch_values_system[i](j) =
                  get_component(tmp[i](j), extract_component);
            }
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_values(
      const FEValuesBase<dim, spacedim> &fe_patch_values,
      const ComponentExtractor           extract_component,
      std::vector<double>               &patch_values) const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_values(vector, patch_values);
        }
      else
        {
          std::vector<ScalarType> tmp(patch_values.size());

          fe_patch_values.get_function_values(vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); ++i)
            patch_values[i] = get_component(tmp[i], extract_component);
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_gradients(
      const FEValuesBase<dim, spacedim>             &fe_patch_values,
      const ComponentExtractor                       extract_component,
      std::vector<std::vector<Tensor<1, spacedim>>> &patch_gradients_system)
      const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_gradients(vector,
                                                 patch_gradients_system);
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<std::vector<Tensor<1, spacedim, ScalarType>>> tmp(
            n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            tmp[i].resize(n_components);

          fe_patch_values.get_function_gradients(vector, tmp);

          AssertDimension(patch_gradients_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            {
              AssertDimension(patch_gradients_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; ++j)
                patch_gradients_system[i][j] =
                  get_component(tmp[i][j], extract_component);
            }
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_gradients(
      const FEValuesBase<dim, spacedim> &fe_patch_values,
      const ComponentExtractor           extract_component,
      std::vector<Tensor<1, spacedim>>  &patch_gradients) const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_gradients(vector, patch_gradients);
        }
      else
        {
          std::vector<Tensor<1, spacedim, ScalarType>> tmp;
          tmp.resize(patch_gradients.size());

          fe_patch_values.get_function_gradients(vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); ++i)
            patch_gradients[i] = get_component(tmp[i], extract_component);
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_hessians(
      const FEValuesBase<dim, spacedim>             &fe_patch_values,
      const ComponentExtractor                       extract_component,
      std::vector<std::vector<Tensor<2, spacedim>>> &patch_hessians_system)
      const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_hessians(vector, patch_hessians_system);
        }
      else
        {
          // OK, so we know that the data type stored by the user-provided
          // vector is not simply a double. In that case, we need to ask the
          // FEValuesBase object to extract function values for us from the
          // evaluation points in the provided data type, and then copy them by
          // hand into the output location.
          const unsigned int n_components =
            fe_patch_values.get_fe().n_components();
          const unsigned int n_eval_points =
            fe_patch_values.n_quadrature_points;

          std::vector<std::vector<Tensor<2, spacedim, ScalarType>>> tmp(
            n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            tmp[i].resize(n_components);

          fe_patch_values.get_function_hessians(vector, tmp);

          AssertDimension(patch_hessians_system.size(), n_eval_points);
          for (unsigned int i = 0; i < n_eval_points; ++i)
            {
              AssertDimension(patch_hessians_system[i].size(), n_components);

              for (unsigned int j = 0; j < n_components; ++j)
                patch_hessians_system[i][j] =
                  get_component(tmp[i][j], extract_component);
            }
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::get_function_hessians(
      const FEValuesBase<dim, spacedim> &fe_patch_values,
      const ComponentExtractor           extract_component,
      std::vector<Tensor<2, spacedim>>  &patch_hessians) const
    {
      if constexpr (std::is_same_v<ScalarType, double>)
        {
          (void)extract_component;
          Assert(extract_component == ComponentExtractor::real_part,
                 ExcMessage("You cannot extract anything other than the real "
                            "part from a real number."));

          fe_patch_values.get_function_hessians(vector, patch_hessians);
        }
      else
        {
          std::vector<Tensor<2, spacedim, ScalarType>> tmp(
            patch_hessians.size());

          fe_patch_values.get_function_hessians(vector, tmp);

          for (unsigned int i = 0; i < tmp.size(); ++i)
            patch_hessians[i] = get_component(tmp[i], extract_component);
        }
    }



    template <int dim, int spacedim, typename ScalarType>
    bool
    DataEntry<dim, spacedim, ScalarType>::is_complex_valued() const
    {
      return numbers::NumberTraits<ScalarType>::is_complex;
    }



    template <int dim, int spacedim, typename ScalarType>
    std::size_t
    DataEntry<dim, spacedim, ScalarType>::memory_consumption() const
    {
      return (vector.memory_consumption() +
              MemoryConsumption::memory_consumption(this->names));
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    DataEntry<dim, spacedim, ScalarType>::clear()
    {
      vector.reinit(0, 0);
      this->dof_handler = nullptr;
    }



    /**
     * Like DataEntry, but used to look up data from multigrid computations.
     * Data will use level-DoF indices to look up in a
     * MGLevelObject<VectorType> given on the specific level instead of
     * interpolating data to coarser cells.
     */
    template <int dim, int spacedim, typename ScalarType>
    class MGDataEntry : public DataEntryBase<dim, spacedim>
    {
    public:
      template <typename VectorType>
      MGDataEntry(const DoFHandler<dim, spacedim> *dofs,
                  const MGLevelObject<VectorType> *vectors,
                  const std::vector<std::string>  &names,
                  const std::vector<
                    DataComponentInterpretation::DataComponentInterpretation>
                    &data_component_interpretation)
        : DataEntryBase<dim, spacedim>(dofs,
                                       names,
                                       data_component_interpretation)
      {
        AssertDimension(vectors->min_level(), 0);

        this->vectors.resize(vectors->min_level(), vectors->max_level());

        for (unsigned int l = vectors->min_level(); l <= vectors->max_level();
             ++l)
          CreateVectors::create_dof_vector(*dofs,
                                           (*vectors)[l],
                                           this->vectors[l],
                                           l);
      }

      virtual double
      get_cell_data_value(
        const unsigned int       cell_number,
        const ComponentExtractor extract_component) const override;

      virtual void
      get_function_values(const FEValuesBase<dim, spacedim> &fe_patch_values,
                          const ComponentExtractor           extract_component,
                          std::vector<double> &patch_values) const override;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual void
      get_function_values(const FEValuesBase<dim, spacedim> &fe_patch_values,
                          const ComponentExtractor           extract_component,
                          std::vector<dealii::Vector<double>>
                            &patch_values_system) const override;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<dim, spacedim> & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<Tensor<1, spacedim>> & /*patch_gradients*/) const override
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_gradients(
        const FEValuesBase<dim, spacedim> & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<std::vector<Tensor<1, spacedim>>>
          & /*patch_gradients_system*/) const override
      {
        DEAL_II_NOT_IMPLEMENTED();
      }


      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<dim, spacedim> & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<Tensor<2, spacedim>> & /*patch_hessians*/) const override
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual void
      get_function_hessians(
        const FEValuesBase<dim, spacedim> & /*fe_patch_values*/,
        const ComponentExtractor /*extract_component*/,
        std::vector<std::vector<Tensor<2, spacedim>>>
          & /*patch_hessians_system*/) const override
      {
        DEAL_II_NOT_IMPLEMENTED();
      }

      /**
       * Return whether the data represented by (a derived class of) this object
       * represents a complex-valued (as opposed to real-valued) information.
       */
      virtual bool
      is_complex_valued() const override
      {
        Assert(numbers::NumberTraits<ScalarType>::is_complex == false,
               ExcNotImplemented());
        return numbers::NumberTraits<ScalarType>::is_complex;
      }

      /**
       * Clear all references to the vectors.
       */
      virtual void
      clear() override
      {
        vectors.clear();
      }

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t
      memory_consumption() const override
      {
        return sizeof(vectors);
      }

    private:
      MGLevelObject<LinearAlgebra::distributed::BlockVector<ScalarType>>
        vectors;

      /**
       * Extract the @p indices from @p vector and put them into @p values.
       */
      void
      extract(const LinearAlgebra::distributed::BlockVector<ScalarType> &vector,
              const std::vector<types::global_dof_index> &indices,
              const ComponentExtractor                    extract_component,
              std::vector<double>                        &values) const
      {
        for (unsigned int i = 0; i < values.size(); ++i)
          values[i] = get_component(vector[indices[i]], extract_component);
      }
    };



    template <int dim, int spacedim, typename ScalarType>
    double
    MGDataEntry<dim, spacedim, ScalarType>::get_cell_data_value(
      const unsigned int /*cell_number*/,
      const ComponentExtractor /*extract_component*/) const
    {
      DEAL_II_NOT_IMPLEMENTED();

      return 0.0;
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    MGDataEntry<dim, spacedim, ScalarType>::get_function_values(
      const FEValuesBase<dim, spacedim> &fe_patch_values,
      const ComponentExtractor           extract_component,
      std::vector<double>               &patch_values) const
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcNotImplemented());

      const typename DoFHandler<dim, spacedim>::level_cell_iterator dof_cell(
        &fe_patch_values.get_cell()->get_triangulation(),
        fe_patch_values.get_cell()->level(),
        fe_patch_values.get_cell()->index(),
        this->dof_handler);

      const auto *vector = &(vectors[dof_cell->level()]);

      const unsigned int dofs_per_cell =
        this->dof_handler->get_fe(0).n_dofs_per_cell();

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      dof_cell->get_mg_dof_indices(dof_indices);
      std::vector<double> values(dofs_per_cell);
      extract(*vector, dof_indices, extract_component, values);
      std::fill(patch_values.begin(), patch_values.end(), 0.0);

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int q_point = 0; q_point < patch_values.size(); ++q_point)
          patch_values[q_point] +=
            values[i] * fe_patch_values.shape_value(i, q_point);
    }



    template <int dim, int spacedim, typename ScalarType>
    void
    MGDataEntry<dim, spacedim, ScalarType>::get_function_values(
      const FEValuesBase<dim, spacedim>   &fe_patch_values,
      const ComponentExtractor             extract_component,
      std::vector<dealii::Vector<double>> &patch_values_system) const
    {
      Assert(extract_component == ComponentExtractor::real_part,
             ExcNotImplemented());

      typename DoFHandler<dim, spacedim>::level_cell_iterator dof_cell(
        &fe_patch_values.get_cell()->get_triangulation(),
        fe_patch_values.get_cell()->level(),
        fe_patch_values.get_cell()->index(),
        this->dof_handler);

      const auto *vector = &(vectors[dof_cell->level()]);

      const unsigned int dofs_per_cell =
        this->dof_handler->get_fe(0).n_dofs_per_cell();

      std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
      dof_cell->get_mg_dof_indices(dof_indices);
      std::vector<double> values(dofs_per_cell);
      extract(*vector, dof_indices, extract_component, values);

      const unsigned int n_components = fe_patch_values.get_fe().n_components();
      const unsigned int n_eval_points = fe_patch_values.n_quadrature_points;

      AssertDimension(patch_values_system.size(), n_eval_points);
      for (unsigned int q = 0; q < n_eval_points; ++q)
        {
          AssertDimension(patch_values_system[q].size(), n_components);
          patch_values_system[q] = 0.0;

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int c = 0; c < n_components; ++c)
              patch_values_system[q](c) +=
                values[i] * fe_patch_values.shape_value_component(i, q, c);
        }
    }

  } // namespace DataOutImplementation
} // namespace internal



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::DataOut_DoFData()
  : triangulation(nullptr, typeid(*this).name())
  , dofs(nullptr, typeid(*this).name())
{}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::~DataOut_DoFData()
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class for clarity be explicit on which function is
  // called
  DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::clear();
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::attach_dof_handler(
  const DoFHandler<dim, spacedim> &d)
{
  Assert(dof_data.empty(),
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());
  Assert(cell_data.empty(),
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());

  triangulation =
    ObserverPointer<const Triangulation<dim, spacedim>>(&d.get_triangulation(),
                                                        typeid(*this).name());
  dofs =
    ObserverPointer<const DoFHandler<dim, spacedim>>(&d, typeid(*this).name());
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::attach_triangulation(
  const Triangulation<dim, spacedim> &tria)
{
  Assert(dof_data.empty(),
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());
  Assert(cell_data.empty(),
         Exceptions::DataOutImplementation::ExcOldDataStillPresent());

  triangulation =
    ObserverPointer<const Triangulation<dim, spacedim>>(&tria,
                                                        typeid(*this).name());
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_data_vector(
  const DoFHandler<dim, spacedim>   &dof_handler,
  const VectorType                  &vec,
  const DataPostprocessor<spacedim> &data_postprocessor)
{
  // this is a specialized version of the other function where we have a
  // postprocessor. if we do, we know that we have type_dof_data, which makes
  // things a bit simpler, we also don't need to deal with some of the other
  // stuff and use a different constructor of DataEntry
  if (triangulation != nullptr)
    {
      Assert(&dof_handler.get_triangulation() == triangulation,
             ExcMessage("The triangulation attached to the DoFHandler does not "
                        "match with the one set previously"));
    }
  else
    {
      triangulation = ObserverPointer<const Triangulation<dim, spacedim>>(
        &dof_handler.get_triangulation(), typeid(*this).name());
    }

  Assert(vec.size() == dof_handler.n_dofs(),
         Exceptions::DataOutImplementation::ExcInvalidVectorSize(
           vec.size(),
           dof_handler.n_dofs(),
           dof_handler.get_triangulation().n_active_cells()));


  auto new_entry = std::make_unique<internal::DataOutImplementation::DataEntry<
    dim,
    spacedim,
    typename numbers::NumberTraits<
      typename VectorType::value_type>::double_type>>(&dof_handler,
                                                      &vec,
                                                      &data_postprocessor);
  dof_data.emplace_back(std::move(new_entry));
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
  add_data_vector_internal(
    const DoFHandler<dim, spacedim> *dof_handler,
    const VectorType                &data_vector,
    const std::vector<std::string>  &names,
    const DataVectorType             type,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
              &data_component_interpretation_,
    const bool deduce_output_names)
{
  // Check available mesh information:
  if (triangulation == nullptr)
    {
      Assert(dof_handler != nullptr, ExcInternalError());
      triangulation = ObserverPointer<const Triangulation<dim, spacedim>>(
        &dof_handler->get_triangulation(), typeid(*this).name());
    }

  if (dof_handler != nullptr)
    {
      Assert(&dof_handler->get_triangulation() == triangulation,
             ExcMessage("The triangulation attached to the DoFHandler does not "
                        "match with the one set previously"));
    }

  // Figure out the data type:
  DataVectorType actual_type = type;
  if (type == type_automatic)
    {
      Assert(
        (dof_handler == nullptr) ||
          (triangulation->n_active_cells() != dof_handler->n_dofs()),
        ExcMessage(
          "Unable to determine the type of vector automatically because the number of DoFs "
          "is equal to the number of cells. Please specify DataVectorType."));

      if (data_vector.size() == triangulation->n_active_cells())
        actual_type = type_cell_data;
      else
        actual_type = type_dof_data;
    }
  Assert(actual_type == type_cell_data || actual_type == type_dof_data,
         ExcInternalError());

  // If necessary, append '_1', '_2', etc. to component names:
  std::vector<std::string> deduced_names;
  if (deduce_output_names && actual_type == type_dof_data)
    {
      Assert(names.size() == 1, ExcInternalError());
      Assert(dof_handler != nullptr, ExcInternalError());
      Assert(dof_handler->n_dofs() > 0,
             ExcMessage("The DoF handler attached to the current output vector "
                        "does not have any degrees of freedom, so it is not "
                        "possible to output DoF data in this context."));
      const std::string &name         = names[0];
      const unsigned int n_components = dof_handler->get_fe(0).n_components();
      deduced_names.resize(n_components);
      if (n_components > 1)
        {
          for (unsigned int i = 0; i < n_components; ++i)
            {
              deduced_names[i] = name + '_' + std::to_string(i);
            }
        }
      else
        {
          deduced_names[0] = name;
        }
    }
  else
    {
      deduced_names = names;
    }

  // Check that things have the right sizes for the data type:
  switch (actual_type)
    {
      case type_cell_data:
        Assert(data_vector.size() == triangulation->n_active_cells(),
               ExcDimensionMismatch(data_vector.size(),
                                    triangulation->n_active_cells()));
        Assert(deduced_names.size() == 1,
               Exceptions::DataOutImplementation::ExcInvalidNumberOfNames(
                 deduced_names.size(), 1));
        break;
      case type_dof_data:
        Assert(dof_handler != nullptr,
               Exceptions::DataOutImplementation::ExcNoDoFHandlerSelected());
        Assert(data_vector.size() == dof_handler->n_dofs(),
               Exceptions::DataOutImplementation::ExcInvalidVectorSize(
                 data_vector.size(),
                 dof_handler->n_dofs(),
                 triangulation->n_active_cells()));
        Assert(deduced_names.size() == dof_handler->get_fe(0).n_components(),
               Exceptions::DataOutImplementation::ExcInvalidNumberOfNames(
                 deduced_names.size(), dof_handler->get_fe(0).n_components()));
        break;
      default:
        DEAL_II_ASSERT_UNREACHABLE();
    }

  const auto &data_component_interpretation =
    (data_component_interpretation_.size() != 0 ?
       data_component_interpretation_ :
       std::vector<DataComponentInterpretation::DataComponentInterpretation>(
         deduced_names.size(),
         DataComponentInterpretation::component_is_scalar));

  // finally, add the data vector:
  auto new_entry = std::make_unique<internal::DataOutImplementation::DataEntry<
    dim,
    spacedim,
    typename numbers::NumberTraits<
      typename VectorType::value_type>::double_type>>(
    dof_handler,
    &data_vector,
    deduced_names,
    data_component_interpretation,
    actual_type);

  if (actual_type == type_dof_data)
    dof_data.emplace_back(std::move(new_entry));
  else
    cell_data.emplace_back(std::move(new_entry));
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_mg_data_vector(
  const DoFHandler<dim, spacedim> &dof_handler,
  const MGLevelObject<VectorType> &data,
  const std::string               &name)
{
  // forward the call to the vector version:
  std::vector<std::string> names(1, name);
  add_mg_data_vector(dof_handler, data, names);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
template <typename VectorType>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::add_mg_data_vector(
  const DoFHandler<dim, spacedim> &dof_handler,
  const MGLevelObject<VectorType> &data,
  const std::vector<std::string>  &names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation_)
{
  if (triangulation == nullptr)
    triangulation = ObserverPointer<const Triangulation<dim, spacedim>>(
      &dof_handler.get_triangulation(), typeid(*this).name());

  Assert(&dof_handler.get_triangulation() == triangulation,
         ExcMessage("The triangulation attached to the DoFHandler does not "
                    "match with the one set previously"));

  const unsigned int       n_components  = dof_handler.get_fe(0).n_components();
  std::vector<std::string> deduced_names = names;

  if (names.size() == 1 && n_components > 1)
    {
      deduced_names.resize(n_components);
      for (unsigned int i = 0; i < n_components; ++i)
        {
          deduced_names[i] = names[0] + '_' + std::to_string(i);
        }
    }

  Assert(deduced_names.size() == n_components,
         ExcMessage("Invalid number of names given."));

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    &data_component_interpretation =
      (data_component_interpretation_.size() != 0 ?
         data_component_interpretation_ :
         std::vector<DataComponentInterpretation::DataComponentInterpretation>(
           n_components, DataComponentInterpretation::component_is_scalar));

  Assert(data_component_interpretation.size() == n_components,
         ExcMessage(
           "Invalid number of entries in data_component_interpretation."));

  auto new_entry = std::make_unique<
    internal::DataOutImplementation::
      MGDataEntry<dim, spacedim, typename VectorType::value_type>>(
    &dof_handler, &data, deduced_names, data_component_interpretation);
  dof_data.emplace_back(std::move(new_entry));
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::clear_data_vectors()
{
  dof_data.erase(dof_data.begin(), dof_data.end());
  cell_data.erase(cell_data.begin(), cell_data.end());

  // delete patches
  std::vector<Patch> dummy;
  patches.swap(dummy);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
  clear_input_data_references()
{
  for (unsigned int i = 0; i < dof_data.size(); ++i)
    dof_data[i]->clear();

  for (unsigned int i = 0; i < cell_data.size(); ++i)
    cell_data[i]->clear();

  if (dofs != nullptr)
    dofs = nullptr;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
void
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::clear()
{
  dof_data.erase(dof_data.begin(), dof_data.end());
  cell_data.erase(cell_data.begin(), cell_data.end());

  if (dofs != nullptr)
    dofs = nullptr;

  // delete patches
  std::vector<Patch> dummy;
  patches.swap(dummy);
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
std::vector<std::string>
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::get_dataset_names()
  const
{
  std::vector<std::string> names;

  // Loop over all DoF-data datasets and push the names. If the
  // vector underlying a data set is complex-valued, then
  // expand it into its real and imaginary part. Note, however,
  // that what comes back from a postprocessor is *always*
  // real-valued, regardless of what goes in, so we don't
  // have this to do this name expansion for data sets that
  // have a postprocessor.
  //
  // The situation is made complicated when considering vector- and
  // tensor-valued component sets. This is because if, for example, we have a
  // complex-valued vector, we don't want to output Re(u_x), then Im(u_x), then
  // Re(u_y), etc. That's because if we did this, then visualization programs
  // will not easily be able to patch together the 1st, 3rd, 5th components into
  // the vector representing the real part of a vector field, and similarly for
  // the 2nd, 4th, 6th component for the imaginary part of the vector field.
  // Rather, we need to put all real components of the same vector field into
  // consecutive components.
  //
  // This sort of logic is also explained in some detail in
  //   DataOut::build_one_patch().
  for (const auto &input_data : dof_data)
    if (input_data->is_complex_valued() == false ||
        (input_data->postprocessor != nullptr))
      {
        for (const auto &name : input_data->names)
          names.push_back(name);
      }
    else
      {
        // OK, so we have a complex-valued vector. We then need to go through
        // all components and order them appropriately
        for (unsigned int i = 0; i < input_data->names.size();
             /* increment of i happens below */)
          {
            switch (input_data->data_component_interpretation[i])
              {
                case DataComponentInterpretation::component_is_scalar:
                  {
                    // It's a scalar. Just output real and imaginary parts one
                    // after the other:
                    names.push_back(input_data->names[i] + "_re");
                    names.push_back(input_data->names[i] + "_im");

                    // Move forward by one component
                    ++i;

                    break;
                  }

                case DataComponentInterpretation::component_is_part_of_vector:
                  {
                    // It's a vector. First output all real parts, then all
                    // imaginary parts:
                    const unsigned int size = patch_spacedim;
                    for (unsigned int vec_comp = 0; vec_comp < size; ++vec_comp)
                      names.push_back(input_data->names[i + vec_comp] + "_re");
                    for (unsigned int vec_comp = 0; vec_comp < size; ++vec_comp)
                      names.push_back(input_data->names[i + vec_comp] + "_im");

                    // Move forward by dim components
                    i += size;

                    break;
                  }

                case DataComponentInterpretation::component_is_part_of_tensor:
                  {
                    // It's a tensor. First output all real parts, then all
                    // imaginary parts:
                    const unsigned int size = patch_spacedim * patch_spacedim;
                    for (unsigned int tensor_comp = 0; tensor_comp < size;
                         ++tensor_comp)
                      names.push_back(input_data->names[i + tensor_comp] +
                                      "_re");
                    for (unsigned int tensor_comp = 0; tensor_comp < size;
                         ++tensor_comp)
                      names.push_back(input_data->names[i + tensor_comp] +
                                      "_im");

                    // Move forward by dim*dim components
                    i += size;

                    break;
                  }

                default:
                  DEAL_II_ASSERT_UNREACHABLE();
              }
          }
      }

  // Do the same as above for cell-type data. This is simpler because it
  // is always scalar, and so we don't have to worry about whether some
  // components together form vectors or tensors.
  for (const auto &input_data : cell_data)
    {
      Assert(input_data->names.size() == 1, ExcInternalError());
      if ((input_data->is_complex_valued() == false) ||
          (input_data->postprocessor != nullptr))
        names.push_back(input_data->names[0]);
      else
        {
          names.push_back(input_data->names[0] + "_re");
          names.push_back(input_data->names[0] + "_im");
        }
    }

  return names;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
  get_nonscalar_data_ranges() const
{
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    ranges;

  // collect the ranges of dof and cell data
  unsigned int output_component = 0;
  for (const auto &input_data : dof_data)
    for (unsigned int i = 0; i < input_data->n_output_variables;
         /* i is updated below */)
      // see what kind of data we have here. note that for the purpose of the
      // current function all we care about is vector data
      switch (input_data->data_component_interpretation[i])
        {
          case DataComponentInterpretation::component_is_scalar:
            {
              // Just move component forward by one (or two if the
              // component happens to be complex-valued and we don't use a
              // postprocessor)
              // -- postprocessors always return real-valued things)
              ++i;
              output_component += (input_data->is_complex_valued() &&
                                       (input_data->postprocessor == nullptr) ?
                                     2 :
                                     1);

              break;
            }

          case DataComponentInterpretation::component_is_part_of_vector:
            {
              // ensure that there is a continuous number of next space_dim
              // components that all deal with vectors
              Assert(
                i + patch_spacedim <= input_data->n_output_variables,
                Exceptions::DataOutImplementation::ExcInvalidVectorDeclaration(
                  i, input_data->names[i]));
              for (unsigned int dd = 1; dd < patch_spacedim; ++dd)
                Assert(
                  input_data->data_component_interpretation[i + dd] ==
                    DataComponentInterpretation::component_is_part_of_vector,
                  Exceptions::DataOutImplementation::
                    ExcInvalidVectorDeclaration(i, input_data->names[i]));

              // all seems right, so figure out whether there is a common
              // name to these components. if not, leave the name empty and
              // let the output format writer decide what to do here
              std::string name = input_data->names[i];
              for (unsigned int dd = 1; dd < patch_spacedim; ++dd)
                if (name != input_data->names[i + dd])
                  {
                    name = "";
                    break;
                  }

              // Finally add a corresponding range. If this is a real-valued
              // vector, then we only need to do this once. But if it is a
              // complex-valued vector and it is not postprocessed, then we need
              // to do it twice -- once for the real parts and once for the
              // imaginary parts
              //
              // This sort of logic is also explained in some detail in
              //   DataOut::build_one_patch().
              if (input_data->is_complex_valued() == false ||
                  (input_data->postprocessor != nullptr))
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_spacedim - 1,
                    name,
                    DataComponentInterpretation::component_is_part_of_vector));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += patch_spacedim;
                  i += patch_spacedim;
                }
              else
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_spacedim - 1,
                    name + "_re",
                    DataComponentInterpretation::component_is_part_of_vector));
                  output_component += patch_spacedim;

                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + patch_spacedim - 1,
                    name + "_im",
                    DataComponentInterpretation::component_is_part_of_vector));
                  output_component += patch_spacedim;

                  i += patch_spacedim;
                }


              break;
            }

          case DataComponentInterpretation::component_is_part_of_tensor:
            {
              const unsigned int size = patch_spacedim * patch_spacedim;
              // ensure that there is a continuous number of next
              // space_dim*space_dim components that all deal with tensors
              Assert(
                i + size <= input_data->n_output_variables,
                Exceptions::DataOutImplementation::ExcInvalidTensorDeclaration(
                  i, input_data->names[i]));
              for (unsigned int dd = 1; dd < size; ++dd)
                Assert(
                  input_data->data_component_interpretation[i + dd] ==
                    DataComponentInterpretation::component_is_part_of_tensor,
                  Exceptions::DataOutImplementation::
                    ExcInvalidTensorDeclaration(i, input_data->names[i]));

              // all seems right, so figure out whether there is a common
              // name to these components. if not, leave the name empty and
              // let the output format writer decide what to do here
              std::string name = input_data->names[i];
              for (unsigned int dd = 1; dd < size; ++dd)
                if (name != input_data->names[i + dd])
                  {
                    name = "";
                    break;
                  }

              // Finally add a corresponding range. If this is a real-valued
              // tensor, then we only need to do this once. But if it is a
              // complex-valued tensor and it is not postprocessed, then we need
              // to do it twice -- once for the real parts and once for the
              // imaginary parts
              //
              // This sort of logic is also explained in some detail in
              //   DataOut::build_one_patch().
              if (input_data->is_complex_valued() == false ||
                  (input_data->postprocessor != nullptr))
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name,
                    DataComponentInterpretation::component_is_part_of_tensor));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += size;
                  i += size;
                }
              else
                {
                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name + "_re",
                    DataComponentInterpretation::component_is_part_of_tensor));
                  output_component += size;

                  ranges.emplace_back(std::forward_as_tuple(
                    output_component,
                    output_component + size - 1,
                    name + "_im",
                    DataComponentInterpretation::component_is_part_of_tensor));
                  output_component += size;

                  i += size;
                }
              break;
            }

          default:
            DEAL_II_NOT_IMPLEMENTED();
        }

  // note that we do not have to traverse the list of cell data here because
  // cell data is one value per (logical) cell and therefore cannot be a
  // vector

  return ranges;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
const std::vector<dealii::DataOutBase::Patch<patch_dim, patch_spacedim>> &
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::get_patches() const
{
  return patches;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
std::vector<std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::get_fes() const
{
  std::vector<std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
    finite_elements(this->dof_data.size());
  for (unsigned int i = 0; i < this->dof_data.size(); ++i)
    {
      Assert(dof_data[i]->dof_handler != nullptr,
             Exceptions::DataOutImplementation::ExcNoDoFHandlerSelected());

      // avoid creating too many finite elements and doing a lot of work on
      // initializing FEValues downstream: if two DoFHandlers are the same
      // (checked by pointer comparison), we can re-use the shared_ptr object
      // for the second one. We cannot check for finite element equalities
      // because we need different FEValues objects for different dof
      // handlers.
      bool duplicate = false;
      for (unsigned int j = 0; j < i; ++j)
        if (dof_data[i]->dof_handler == dof_data[j]->dof_handler)
          {
            finite_elements[i] = finite_elements[j];
            duplicate          = true;
          }
      if (duplicate == false)
        finite_elements[i] =
          std::make_shared<dealii::hp::FECollection<dim, spacedim>>(
            this->dof_data[i]->dof_handler->get_fe_collection());
    }
  if (this->dof_data.empty())
    {
      Assert(triangulation != nullptr, ExcNotImplemented());

      const auto &reference_cells = triangulation->get_reference_cells();
      finite_elements.reserve(reference_cells.size());

      // TODO: below we select linear, discontinuous elements for simplex,
      // wedge, and pyramid; we would like to select constant functions, but
      // that is not implemented yet
      for (const auto &reference_cell : reference_cells)
        {
          if (reference_cell.is_hyper_cube())
            finite_elements.emplace_back(
              std::make_shared<dealii::hp::FECollection<dim, spacedim>>(
                FE_DGQ<dim, spacedim>(0)));
          else if (reference_cell.is_simplex())
            finite_elements.emplace_back(
              std::make_shared<dealii::hp::FECollection<dim, spacedim>>(
                FE_SimplexDGP<dim, spacedim>(1)));
          else if (reference_cell == ReferenceCells::Wedge)
            finite_elements.emplace_back(
              std::make_shared<dealii::hp::FECollection<dim, spacedim>>(
                FE_WedgeDGP<dim, spacedim>(1)));
          else if (reference_cell == ReferenceCells::Pyramid)
            finite_elements.emplace_back(
              std::make_shared<dealii::hp::FECollection<dim, spacedim>>(
                FE_PyramidDGP<dim, spacedim>(1)));
          else
            DEAL_II_NOT_IMPLEMENTED();
        }
    }
  return finite_elements;
}



template <int dim, int patch_dim, int spacedim, int patch_spacedim>
std::size_t
DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::memory_consumption()
  const
{
  return (DataOutInterface<patch_dim, patch_spacedim>::memory_consumption() +
          MemoryConsumption::memory_consumption(dofs) +
          MemoryConsumption::memory_consumption(patches));
}

DEAL_II_NAMESPACE_CLOSE

#endif
