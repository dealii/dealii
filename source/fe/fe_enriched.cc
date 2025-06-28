// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_enriched.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <memory>
#include <set>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FE_Enriched
  {
    namespace
    {
      /**
       * Auxiliary function to create multiplicity vector from input enrichment
       * functions.
       */
      template <typename T>
      std::vector<unsigned int>
      build_multiplicities(const std::vector<std::vector<T>> &functions)
      {
        std::vector<unsigned int> multiplicities;
        multiplicities.push_back(1); // the first one is non-enriched FE
        for (unsigned int i = 0; i < functions.size(); ++i)
          multiplicities.push_back(functions[i].size());

        return multiplicities;
      }


      /**
       * Auxiliary function to build FiniteElement's vector
       */
      template <int dim, int spacedim>
      std::vector<const FiniteElement<dim, spacedim> *>
      build_fes(
        const FiniteElement<dim, spacedim>                      *fe_base,
        const std::vector<const FiniteElement<dim, spacedim> *> &fe_enriched)
      {
        std::vector<const FiniteElement<dim, spacedim> *> fes;
        fes.push_back(fe_base);
        for (unsigned int i = 0; i < fe_enriched.size(); ++i)
          fes.push_back(fe_enriched[i]);

        return fes;
      }


      /**
       * Auxiliary function which check consistency of the input parameters.
       * Returns true if everything is ok.
       */
      template <int dim, int spacedim>
      bool
      consistency_check(
        const std::vector<const FiniteElement<dim, spacedim> *> &fes,
        const std::vector<unsigned int>                         &multiplicities,
        const std::vector<std::vector<std::function<const Function<spacedim> *(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator
            &)>>>                                               &functions)
      {
        AssertThrow(fes.size() > 0, ExcMessage("FEs size should be >=1"));
        AssertThrow(fes.size() == multiplicities.size(),
                    ExcMessage(
                      "FEs and multiplicities should have the same size"));

        AssertThrow(functions.size() == fes.size() - 1,
                    ExcDimensionMismatch(functions.size(), fes.size() - 1));

        AssertThrow(multiplicities[0] == 1,
                    ExcMessage("First multiplicity should be 1"));

        const unsigned int n_comp_base = fes[0]->n_components();

        // start from fe=1 as 0th is always non-enriched FE.
        for (unsigned int fe = 1; fe < fes.size(); ++fe)
          {
            const FE_Nothing<dim> *fe_nothing =
              dynamic_cast<const FE_Nothing<dim> *>(fes[fe]);
            if (fe_nothing)
              AssertThrow(
                fe_nothing->is_dominating(),
                ExcMessage(
                  "Only dominating FE_Nothing can be used in FE_Enriched"));

            AssertThrow(
              fes[fe]->n_components() == n_comp_base,
              ExcMessage(
                "All elements must have the same number of components"));
          }
        return true;
      }


      /**
       * Auxiliary function which determines whether the FiniteElement will be
       * enriched.
       */
      template <int dim, int spacedim>
      bool
      check_if_enriched(
        const std::vector<const FiniteElement<dim, spacedim> *> &fes)
      {
        // start from fe=1 as 0th is always non-enriched FE.
        for (unsigned int fe = 1; fe < fes.size(); ++fe)
          if (dynamic_cast<const FE_Nothing<dim> *>(fes[fe]) == nullptr)
            // this is not FE_Nothing => there will be enrichment
            return true;

        return false;
      }
    } // namespace
  }   // namespace FE_Enriched
} // namespace internal


template <int dim, int spacedim>
FE_Enriched<dim, spacedim>::FE_Enriched(
  const FiniteElement<dim, spacedim> &fe_base)
  : FE_Enriched<dim, spacedim>(fe_base,
                               FE_Nothing<dim, spacedim>(fe_base.n_components(),
                                                         true),
                               nullptr)
{}


template <int dim, int spacedim>
FE_Enriched<dim, spacedim>::FE_Enriched(
  const FiniteElement<dim, spacedim> &fe_base,
  const FiniteElement<dim, spacedim> &fe_enriched,
  const Function<spacedim>           *enrichment_function)
  : FE_Enriched<dim, spacedim>(
      &fe_base,
      std::vector<const FiniteElement<dim, spacedim> *>(1, &fe_enriched),
      std::vector<std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>>(
        1,
        std::vector<std::function<const Function<spacedim> *(
          const typename Triangulation<dim, spacedim>::cell_iterator &)>>(
          1,
          [=](const typename Triangulation<dim, spacedim>::cell_iterator &)
            -> const Function<spacedim> * { return enrichment_function; })))
{}


template <int dim, int spacedim>
FE_Enriched<dim, spacedim>::FE_Enriched(
  const FiniteElement<dim, spacedim>                      *fe_base,
  const std::vector<const FiniteElement<dim, spacedim> *> &fe_enriched,
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>> &functions)
  : FE_Enriched<dim, spacedim>(
      internal::FE_Enriched::build_fes(fe_base, fe_enriched),
      internal::FE_Enriched::build_multiplicities(functions),
      functions)
{}


template <int dim, int spacedim>
FE_Enriched<dim, spacedim>::FE_Enriched(
  const std::vector<const FiniteElement<dim, spacedim> *> &fes,
  const std::vector<unsigned int>                         &multiplicities,
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>> &functions)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(fes, multiplicities, false),
      FETools::Compositing::compute_restriction_is_additive_flags(
        fes,
        multiplicities),
      FETools::Compositing::compute_nonzero_components(fes,
                                                       multiplicities,
                                                       false))
  , enrichments(functions)
  , is_enriched(internal::FE_Enriched::check_if_enriched(fes))
  , fe_system(std::make_unique<FESystem<dim, spacedim>>(fes, multiplicities))
{
  // descriptive error are thrown within the function.
  Assert(internal::FE_Enriched::consistency_check(fes,
                                                  multiplicities,
                                                  functions),
         ExcInternalError());

  initialize(fes, multiplicities);

  // resize to be consistent with all FEs used to construct the FE_Enriched,
  // even though we will never use the 0th element.
  base_no_mult_local_enriched_dofs.resize(fes.size());
  for (unsigned int fe = 1; fe < fes.size(); ++fe)
    base_no_mult_local_enriched_dofs[fe].resize(multiplicities[fe]);

  Assert(base_no_mult_local_enriched_dofs.size() == this->n_base_elements(),
         ExcDimensionMismatch(base_no_mult_local_enriched_dofs.size(),
                              this->n_base_elements()));

  // build the map: (base_no, base_m) -> vector of local element DoFs
  for (unsigned int system_index = 0; system_index < this->n_dofs_per_cell();
       ++system_index)
    {
      const unsigned int base_no =
        this->system_to_base_table[system_index].first.first;
      if (base_no == 0) // 0th is always non-enriched FE
        continue;

      const unsigned int base_m =
        this->system_to_base_table[system_index].first.second;

      Assert(base_m < base_no_mult_local_enriched_dofs[base_no].size(),
             ExcMessage(
               "Size mismatch for base_no_mult_local_enriched_dofs: "
               "base_index = " +
               std::to_string(this->system_to_base_table[system_index].second) +
               "; base_no = " + std::to_string(base_no) +
               "; base_m = " + std::to_string(base_m) +
               "; system_index = " + std::to_string(system_index)));

      Assert(base_m < base_no_mult_local_enriched_dofs[base_no].size(),
             ExcDimensionMismatch(
               base_m, base_no_mult_local_enriched_dofs[base_no].size()));

      base_no_mult_local_enriched_dofs[base_no][base_m].push_back(system_index);
    }

  // make sure that local_enriched_dofs.size() is correct, that is equals to
  // DoFs per cell of the corresponding FE.
  for (unsigned int base_no = 1;
       base_no < base_no_mult_local_enriched_dofs.size();
       base_no++)
    {
      for (unsigned int m = 0;
           m < base_no_mult_local_enriched_dofs[base_no].size();
           m++)
        Assert(base_no_mult_local_enriched_dofs[base_no][m].size() ==
                 fes[base_no]->n_dofs_per_cell(),
               ExcDimensionMismatch(
                 base_no_mult_local_enriched_dofs[base_no][m].size(),
                 fes[base_no]->n_dofs_per_cell()));
    }
}


template <int dim, int spacedim>
std::vector<std::vector<std::function<const Function<spacedim> *(
  const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
FE_Enriched<dim, spacedim>::get_enrichments() const
{
  return enrichments;
}


template <int dim, int spacedim>
double
FE_Enriched<dim, spacedim>::shape_value(const unsigned int i,
                                        const Point<dim>  &p) const
{
  Assert(
    !is_enriched,
    ExcMessage(
      "For enriched finite elements shape_value() can not be defined on the reference element."));
  return fe_system->shape_value(i, p);
}


template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FE_Enriched<dim, spacedim>::clone() const
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  std::vector<unsigned int>                         multiplicities;

  for (unsigned int i = 0; i < this->n_base_elements(); ++i)
    {
      fes.push_back(&base_element(i));
      multiplicities.push_back(this->element_multiplicity(i));
    }

  return std::unique_ptr<FE_Enriched<dim, spacedim>>(
    new FE_Enriched<dim, spacedim>(fes, multiplicities, get_enrichments()));
}


template <int dim, int spacedim>
UpdateFlags
FE_Enriched<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = fe_system->requires_update_flags(flags);

  if (is_enriched)
    {
      // if we ask for values or gradients, then we would need quadrature points
      if (flags & (update_values | update_gradients))
        out |= update_quadrature_points;

      // if need gradients, add update_values due to product rule
      if (out & update_gradients)
        out |= update_values;
    }

  Assert(!(flags & update_3rd_derivatives), ExcNotImplemented());

  return out;
}


template <int dim, int spacedim>
template <int dim_1>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_Enriched<dim, spacedim>::setup_data(
  std::unique_ptr<typename FESystem<dim, spacedim>::InternalData> fes_data,
  const UpdateFlags                                               flags,
  const Quadrature<dim_1> &quadrature) const
{
  // Pass ownership of the FiniteElement::InternalDataBase object
  // that fes_data points to, to the new InternalData object.
  auto update_each_flags = fes_data->update_each;
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
        data_ptr = std::make_unique<InternalData>(std::move(fes_data));
  auto &data     = dynamic_cast<InternalData &>(*data_ptr);

  // copy update_each from FESystem data:
  data.update_each = update_each_flags;

  // resize cache array according to requested flags
  data.enrichment.resize(this->n_base_elements());

  const unsigned int n_q_points = quadrature.size();

  for (unsigned int base = 0; base < this->n_base_elements(); ++base)
    {
      data.enrichment[base].resize(this->element_multiplicity(base));
      for (unsigned int m = 0; m < this->element_multiplicity(base); ++m)
        {
          if (flags & update_values)
            data.enrichment[base][m].values.resize(n_q_points);

          if (flags & update_gradients)
            data.enrichment[base][m].gradients.resize(n_q_points);

          if (flags & update_hessians)
            data.enrichment[base][m].hessians.resize(n_q_points);
        }
    }

  return data_ptr;
}


template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_Enriched<dim, spacedim>::get_face_data(
  const UpdateFlags               update_flags,
  const Mapping<dim, spacedim>   &mapping,
  const hp::QCollection<dim - 1> &quadrature,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  AssertDimension(quadrature.size(), 1);

  auto data =
    fe_system->get_face_data(update_flags, mapping, quadrature, output_data);
  return setup_data(Utilities::dynamic_unique_cast<
                      typename FESystem<dim, spacedim>::InternalData>(
                      std::move(data)),
                    update_flags,
                    quadrature[0]);
}


template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_Enriched<dim, spacedim>::get_subface_data(
  const UpdateFlags             update_flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1>    &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  auto data =
    fe_system->get_subface_data(update_flags, mapping, quadrature, output_data);
  return setup_data(Utilities::dynamic_unique_cast<
                      typename FESystem<dim, spacedim>::InternalData>(
                      std::move(data)),
                    update_flags,
                    quadrature);
}


template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FE_Enriched<dim, spacedim>::get_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim>        &quadrature,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  auto data = fe_system->get_data(flags, mapping, quadrature, output_data);
  return setup_data(Utilities::dynamic_unique_cast<
                      typename FESystem<dim, spacedim>::InternalData>(
                      std::move(data)),
                    flags,
                    quadrature);
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::initialize(
  const std::vector<const FiniteElement<dim, spacedim> *> &fes,
  const std::vector<unsigned int>                         &multiplicities)
{
  Assert(fes.size() == multiplicities.size(),
         ExcDimensionMismatch(fes.size(), multiplicities.size()));

  // Note that we need to skip every FE with multiplicity 0 in the following
  // block of code
  this->base_to_block_indices.reinit(0, 0);

  for (unsigned int i = 0; i < fes.size(); ++i)
    if (multiplicities[i] > 0)
      this->base_to_block_indices.push_back(multiplicities[i]);

  {
    // If the system is not primitive, these have not been initialized by
    // FiniteElement
    this->system_to_component_table.resize(this->n_dofs_per_cell());

    FETools::Compositing::build_cell_tables(this->system_to_base_table,
                                            this->system_to_component_table,
                                            this->component_to_base_table,
                                            *this,
                                            false);

    this->face_system_to_component_table.resize(this->n_unique_faces());

    for (unsigned int face_no = 0; face_no < this->n_unique_faces(); ++face_no)
      {
        this->face_system_to_component_table[0].resize(
          this->n_dofs_per_face(face_no));


        FETools::Compositing::build_face_tables(
          this->face_system_to_base_table[face_no],
          this->face_system_to_component_table[face_no],
          *this,
          false,
          face_no);
      }
  }

  // restriction and prolongation matrices are built on demand

  // now set up the interface constraints for h-refinement.
  // take them from fe_system:
  this->interface_constraints = fe_system->interface_constraints;

  // if we just wrap another FE (i.e. use FE_Nothing as a second FE)
  // then it makes sense to have support points.
  // However, functions like interpolate_boundary_values() need all FEs inside
  // FECollection to be able to provide support points irrespectively whether
  // this FE sits on the boundary or not. Thus for moment just copy support
  // points from FE system:
  {
    this->unit_support_points      = fe_system->unit_support_points;
    this->unit_face_support_points = fe_system->unit_face_support_points;
  }

  // take adjust_quad_dof_index_for_face_orientation_table from FESystem:
  {
    this->adjust_line_dof_index_for_line_orientation_table =
      fe_system->adjust_line_dof_index_for_line_orientation_table;
  }
}


template <int dim, int spacedim>
std::string
FE_Enriched<dim, spacedim>::get_name() const
{
  std::ostringstream namebuf;

  namebuf << "FE_Enriched<" << Utilities::dim_string(dim, spacedim) << ">[";
  for (unsigned int i = 0; i < this->n_base_elements(); ++i)
    {
      namebuf << base_element(i).get_name();
      if (this->element_multiplicity(i) != 1)
        namebuf << '^' << this->element_multiplicity(i);
      if (i != this->n_base_elements() - 1)
        namebuf << '-';
    }
  namebuf << ']';

  return namebuf.str();
}


template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FE_Enriched<dim, spacedim>::base_element(const unsigned int index) const
{
  return fe_system->base_element(index);
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim>                                      &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // call FESystem's method to fill everything without enrichment function
  fe_system->fill_fe_values(cell,
                            cell_similarity,
                            quadrature,
                            mapping,
                            mapping_internal,
                            mapping_data,
                            *fe_data.fesystem_data,
                            output_data);

  if (is_enriched)
    multiply_by_enrichment(
      quadrature, fe_data, mapping_data, cell, output_data);
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1>                             &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  AssertDimension(quadrature.size(), 1);

  // call FESystem's method to fill everything without enrichment function
  fe_system->fill_fe_face_values(cell,
                                 face_no,
                                 quadrature,
                                 mapping,
                                 mapping_internal,
                                 mapping_data,
                                 *fe_data.fesystem_data,
                                 output_data);

  if (is_enriched)
    multiply_by_enrichment(
      quadrature[0], fe_data, mapping_data, cell, output_data);
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const Quadrature<dim - 1>                                  &quadrature,
  const Mapping<dim, spacedim>                               &mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase    &mapping_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                                &mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // call FESystem's method to fill everything without enrichment function
  fe_system->fill_fe_subface_values(cell,
                                    face_no,
                                    sub_no,
                                    quadrature,
                                    mapping,
                                    mapping_internal,
                                    mapping_data,
                                    *fe_data.fesystem_data,
                                    output_data);

  if (is_enriched)
    multiply_by_enrichment(
      quadrature, fe_data, mapping_data, cell, output_data);
}


template <int dim, int spacedim>
template <int dim_1>
void
FE_Enriched<dim, spacedim>::multiply_by_enrichment(
  const Quadrature<dim_1> &quadrature,
  const InternalData      &fe_data,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
                                                             &mapping_data,
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  // mapping_data will contain quadrature points on the real element.
  // fe_internal is needed to get update flags
  // finally, output_data should store all the results we need.

  // Either dim_1==dim
  // (fill_fe_values) or dim_1==dim-1
  // (fill_fe_(sub)face_values)
  Assert(dim_1 == dim || dim_1 == dim - 1, ExcInternalError());
  const UpdateFlags flags = fe_data.update_each;

  const unsigned int n_q_points = quadrature.size();

  // First, populate output_data object (that shall hold everything requested
  // such as shape value, gradients, hessians, etc) from each base element. That
  // is almost identical to FESystem::compute_fill_one_base(), the difference
  // being that we do it irrespectively of cell_similarity and use
  // base_fe_data.update_flags

  // TODO: do we need it only for dim_1 == dim (i.e. fill_fe_values)?
  if (dim_1 == dim)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); ++base_no)
      {
        const FiniteElement<dim, spacedim> &base_fe = base_element(base_no);
        typename FiniteElement<dim, spacedim>::InternalDataBase &base_fe_data =
          fe_data.get_fe_data(base_no);
        internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                   spacedim>
          &base_data = fe_data.get_fe_output_object(base_no);

        const UpdateFlags base_flags = base_fe_data.update_each;

        for (unsigned int system_index = 0;
             system_index < this->n_dofs_per_cell();
             ++system_index)
          if (this->system_to_base_table[system_index].first.first == base_no)
            {
              const unsigned int base_index =
                this->system_to_base_table[system_index].second;
              Assert(base_index < base_fe.n_dofs_per_cell(),
                     ExcInternalError());

              // now copy. if the shape function is primitive, then there
              // is only one value to be copied, but for non-primitive
              // elements, there might be more values to be copied
              //
              // so, find out from which index to take this one value, and
              // to which index to put
              unsigned int out_index = 0;
              for (unsigned int i = 0; i < system_index; ++i)
                out_index += this->n_nonzero_components(i);
              unsigned int in_index = 0;
              for (unsigned int i = 0; i < base_index; ++i)
                in_index += base_fe.n_nonzero_components(i);

              // then loop over the number of components to be copied
              Assert(this->n_nonzero_components(system_index) ==
                       base_fe.n_nonzero_components(base_index),
                     ExcInternalError());
              for (unsigned int s = 0;
                   s < this->n_nonzero_components(system_index);
                   ++s)
                {
                  if (base_flags & update_values)
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      output_data.shape_values[out_index + s][q] =
                        base_data.shape_values(in_index + s, q);

                  if (base_flags & update_gradients)
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      output_data.shape_gradients[out_index + s][q] =
                        base_data.shape_gradients[in_index + s][q];

                  if (base_flags & update_hessians)
                    for (unsigned int q = 0; q < n_q_points; ++q)
                      output_data.shape_hessians[out_index + s][q] =
                        base_data.shape_hessians[in_index + s][q];
                }
            }
      }

  Assert(base_no_mult_local_enriched_dofs.size() == fe_data.enrichment.size(),
         ExcDimensionMismatch(base_no_mult_local_enriched_dofs.size(),
                              fe_data.enrichment.size()));
  // calculate hessians, gradients and values for each function
  for (unsigned int base_no = 1; base_no < this->n_base_elements(); ++base_no)
    {
      Assert(
        base_no_mult_local_enriched_dofs[base_no].size() ==
          fe_data.enrichment[base_no].size(),
        ExcDimensionMismatch(base_no_mult_local_enriched_dofs[base_no].size(),
                             fe_data.enrichment[base_no].size()));
      for (unsigned int m = 0;
           m < base_no_mult_local_enriched_dofs[base_no].size();
           m++)
        {
          // Avoid evaluating quadrature points if no dofs are assigned. This
          // happens when FE_Nothing is used together with other FE (i.e. FE_Q)
          // as enrichments.
          if (base_no_mult_local_enriched_dofs[base_no][m].empty())
            continue;

          Assert(enrichments[base_no - 1][m](cell) != nullptr,
                 ExcMessage(
                   "The pointer to the enrichment function is not set"));

          Assert(enrichments[base_no - 1][m](cell)->n_components == 1,
                 ExcMessage(
                   "Only scalar-valued enrichment functions are allowed"));

          if (flags & update_hessians)
            {
              Assert(fe_data.enrichment[base_no][m].hessians.size() ==
                       n_q_points,
                     ExcDimensionMismatch(
                       fe_data.enrichment[base_no][m].hessians.size(),
                       n_q_points));
              for (unsigned int q = 0; q < n_q_points; ++q)
                fe_data.enrichment[base_no][m].hessians[q] =
                  enrichments[base_no - 1][m](cell)->hessian(
                    mapping_data.quadrature_points[q]);
            }

          if (flags & update_gradients)
            {
              Assert(fe_data.enrichment[base_no][m].gradients.size() ==
                       n_q_points,
                     ExcDimensionMismatch(
                       fe_data.enrichment[base_no][m].gradients.size(),
                       n_q_points));
              for (unsigned int q = 0; q < n_q_points; ++q)
                fe_data.enrichment[base_no][m].gradients[q] =
                  enrichments[base_no - 1][m](cell)->gradient(
                    mapping_data.quadrature_points[q]);
            }

          if (flags & update_values)
            {
              Assert(fe_data.enrichment[base_no][m].values.size() == n_q_points,
                     ExcDimensionMismatch(
                       fe_data.enrichment[base_no][m].values.size(),
                       n_q_points));
              for (unsigned int q = 0; q < n_q_points; ++q)
                fe_data.enrichment[base_no][m].values[q] =
                  enrichments[base_no - 1][m](cell)->value(
                    mapping_data.quadrature_points[q]);
            }
        }
    }

  // Finally, update the standard data stored in output_data
  // by expanding the product rule for enrichment function.
  // note that the order if important, namely
  // output_data.shape_XYZ contains values of standard FEM and we overwrite
  // it with the updated one in the following order: hessians -> gradients ->
  // values
  if (flags & update_hessians)
    {
      for (unsigned int base_no = 1; base_no < this->n_base_elements();
           base_no++)
        {
          for (unsigned int m = 0;
               m < base_no_mult_local_enriched_dofs[base_no].size();
               m++)
            for (unsigned int i = 0;
                 i < base_no_mult_local_enriched_dofs[base_no][m].size();
                 i++)
              {
                const unsigned int enriched_dof =
                  base_no_mult_local_enriched_dofs[base_no][m][i];
                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    const Tensor<2, spacedim> grad_grad = outer_product(
                      output_data.shape_gradients[enriched_dof][q],
                      fe_data.enrichment[base_no][m].gradients[q]);
                    const Tensor<2, spacedim, double> sym_grad_grad =
                      symmetrize(grad_grad) * 2.0; // symmetrize does [s+s^T]/2

                    output_data.shape_hessians[enriched_dof][q] *=
                      fe_data.enrichment[base_no][m].values[q];
                    output_data.shape_hessians[enriched_dof][q] +=
                      sym_grad_grad +
                      output_data.shape_values[enriched_dof][q] *
                        fe_data.enrichment[base_no][m].hessians[q];
                  }
              }
        }
    }

  if (flags & update_gradients)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); ++base_no)
      {
        for (unsigned int m = 0;
             m < base_no_mult_local_enriched_dofs[base_no].size();
             m++)
          for (unsigned int i = 0;
               i < base_no_mult_local_enriched_dofs[base_no][m].size();
               i++)
            {
              const unsigned int enriched_dof =
                base_no_mult_local_enriched_dofs[base_no][m][i];
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  output_data.shape_gradients[enriched_dof][q] *=
                    fe_data.enrichment[base_no][m].values[q];
                  output_data.shape_gradients[enriched_dof][q] +=
                    output_data.shape_values[enriched_dof][q] *
                    fe_data.enrichment[base_no][m].gradients[q];
                }
            }
      }

  if (flags & update_values)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); ++base_no)
      {
        for (unsigned int m = 0;
             m < base_no_mult_local_enriched_dofs[base_no].size();
             m++)
          for (unsigned int i = 0;
               i < base_no_mult_local_enriched_dofs[base_no][m].size();
               i++)
            {
              const unsigned int enriched_dof =
                base_no_mult_local_enriched_dofs[base_no][m][i];
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  output_data.shape_values[enriched_dof][q] *=
                    fe_data.enrichment[base_no][m].values[q];
                }
            }
      }
}


template <int dim, int spacedim>
const FESystem<dim, spacedim> &
FE_Enriched<dim, spacedim>::get_fe_system() const
{
  return *fe_system;
}


template <int dim, int spacedim>
bool
FE_Enriched<dim, spacedim>::hp_constraints_are_implemented() const
{
  return true;
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source,
  FullMatrix<double>                 &matrix,
  const unsigned int                  face_no) const
{
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&source))
    {
      fe_system->get_face_interpolation_matrix(fe_enr_other->get_fe_system(),
                                               matrix,
                                               face_no);
    }
  else
    {
      AssertThrow(
        false,
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));
    }
}


template <int dim, int spacedim>
void
FE_Enriched<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &source,
  const unsigned int                  subface,
  FullMatrix<double>                 &matrix,
  const unsigned int                  face_no) const
{
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&source))
    {
      fe_system->get_subface_interpolation_matrix(fe_enr_other->get_fe_system(),
                                                  subface,
                                                  matrix,
                                                  face_no);
    }
  else
    {
      AssertThrow(
        false,
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Enriched<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&fe_other))
    {
      return fe_system->hp_vertex_dof_identities(fe_enr_other->get_fe_system());
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Enriched<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&fe_other))
    {
      return fe_system->hp_line_dof_identities(fe_enr_other->get_fe_system());
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_Enriched<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&fe_other))
    {
      return fe_system->hp_quad_dof_identities(fe_enr_other->get_fe_system(),
                                               face_no);
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Enriched<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face/cell domination
  // --------------------------------
  // need to decide which element constrain another.
  // for example Q(2) dominate Q(4) and thus some DoFs of Q(4) will be
  // constrained. If we have Q(2) and Q(4)+POU, then it's clear that Q(2)
  // dominates, namely our DoFs will be constrained to make field continuous.
  // However, we need to check for situations like Q(4) vs Q(2)+POU.
  // In that case the domination for the underlying FEs should be the other way,
  // but this implies that we can't constrain POU dofs to make the field
  // continuous. In that case, throw an error

  // if it's also enriched, do domination based on each one's FESystem
  if (const FE_Enriched<dim, spacedim> *fe_enr_other =
        dynamic_cast<const FE_Enriched<dim, spacedim> *>(&fe_other))
    {
      return fe_system->compare_for_domination(fe_enr_other->get_fe_system(),
                                               codim);
    }
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return FiniteElementDomination::neither_element_dominates;
    }
}


template <int dim, int spacedim>
const FullMatrix<double> &
FE_Enriched<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  return fe_system->get_prolongation_matrix(child, refinement_case);
}


template <int dim, int spacedim>
const FullMatrix<double> &
FE_Enriched<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  return fe_system->get_restriction_matrix(child, refinement_case);
}


/* ----------------------- FESystem::InternalData ------------------- */


template <int dim, int spacedim>
FE_Enriched<dim, spacedim>::InternalData::InternalData(
  std::unique_ptr<typename FESystem<dim, spacedim>::InternalData> fesystem_data)
  : fesystem_data(std::move(fesystem_data))
{}


template <int dim, int spacedim>
typename FiniteElement<dim, spacedim>::InternalDataBase &
FE_Enriched<dim, spacedim>::InternalData::get_fe_data(
  const unsigned int base_no) const
{
  return fesystem_data->get_fe_data(base_no);
}


template <int dim, int spacedim>
internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
FE_Enriched<dim, spacedim>::InternalData::get_fe_output_object(
  const unsigned int base_no) const
{
  return fesystem_data->get_fe_output_object(base_no);
}


namespace ColorEnriched
{
  namespace internal
  {
    template <int dim, int spacedim>
    bool
    find_connection_between_subdomains(
      const DoFHandler<dim, spacedim>         &dof_handler,
      const predicate_function<dim, spacedim> &predicate_1,
      const predicate_function<dim, spacedim> &predicate_2)
    {
      // Use a vector to mark vertices
      std::vector<bool> vertices_subdomain_1(
        dof_handler.get_triangulation().n_vertices(), false);

      // Mark vertices that belong to cells in subdomain 1
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (predicate_1(cell)) // True ==> part of subdomain 1
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            vertices_subdomain_1[cell->vertex_index(v)] = true;

      // Find if cells in subdomain 2 and subdomain 1 share vertices.
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (predicate_2(cell)) // True ==> part of subdomain 2
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            if (vertices_subdomain_1[cell->vertex_index(v)] == true)
              {
                return true;
              }
      return false;
    }



    template <int dim, int spacedim>
    unsigned int
    color_predicates(
      const DoFHandler<dim, spacedim>                      &mesh,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      std::vector<unsigned int>                            &predicate_colors)
    {
      const unsigned int num_indices = predicates.size();

      // Use sparsity pattern to represent connections between subdomains.
      // Each predicate (i.e a subdomain) is a node in the graph.
      DynamicSparsityPattern dsp;
      dsp.reinit(num_indices, num_indices);

      /*
       * Find connections between subdomains taken two at a time.
       * If the connection exists, add it to a graph object represented
       * by dynamic sparsity pattern.
       */
      for (unsigned int i = 0; i < num_indices; ++i)
        for (unsigned int j = i + 1; j < num_indices; ++j)
          if (internal::find_connection_between_subdomains(mesh,
                                                           predicates[i],
                                                           predicates[j]))
            dsp.add(i, j);

      dsp.symmetrize();

      // Copy dynamic sparsity pattern to sparsity pattern needed by
      // coloring function
      SparsityPattern sp_graph;
      sp_graph.copy_from(dsp);
      predicate_colors.resize(num_indices);

      // Assign each predicate with a color and return number of colors
      return SparsityTools::color_sparsity_pattern(sp_graph, predicate_colors);
    }



    template <int dim, int spacedim>
    void
    set_cellwise_color_set_and_fe_index(
      DoFHandler<dim, spacedim>                            &dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      const std::vector<unsigned int>                      &predicate_colors,
      std::map<unsigned int, std::map<unsigned int, unsigned int>>
                                          &cellwise_color_predicate_map,
      std::vector<std::set<unsigned int>> &fe_sets)
    {
      // clear output variables first
      fe_sets.clear();
      cellwise_color_predicate_map.clear();

      /*
       * Add first element of fe_sets which is empty by default. This means that
       * the default, FE index = 0 is associated with an empty set, since no
       * predicate is active in these regions.
       */
      fe_sets.resize(1);

      /*
       * Loop through cells and find set of predicate colors associated
       * with the cell. As an example, a cell with an FE index associated
       * with colors {a,b} means that predicates active in the cell have
       * colors a or b.
       *
       * Create new active FE index in case of the color
       * set is not already listed in fe_sets. If the set already exists,
       * find index of the set in fe_sets. In either case, use the id in
       * fe_sets to modify cell->active_fe_index.
       *
       * Associate each cell_id with a set of pairs. The pair represents
       * predicate color and the active predicate with that color.
       * Each color can only correspond to a single predicate since
       * predicates with the same color correspond to disjoint domains.
       * This is what the graph coloring in color_predicates
       * function ensures. The number of pairs is equal to the number
       * of predicates active in the given cell.
       */
      unsigned int map_index = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          // set default FE index ==> no enrichment and no active predicates
          cell->set_active_fe_index(0);

          // Give each cell a unique id, which the cellwise_color_predicate_map
          // can later use to associate a colors of active predicates with
          // the actual predicate id.
          //
          // When the grid is refined, material id is inherited to
          // children cells. So, the cellwise_color_predicate_map stays
          // relevant.
          cell->set_material_id(map_index);
          std::set<unsigned int> color_list;

          // loop through active predicates for the cell and insert map.
          // Eg: if the cell with material id 100 has active
          // predicates 4 (color = 1) and 5 (color = 2), the map will insert
          // pairs (1, 4) and (2, 5) at key 100 (i.e unique id of cell is
          // mapped with a map which associates color with predicate id)
          // Note that color list for the cell would be {1,2}.
          std::map<unsigned int, unsigned int> &cell_map =
            cellwise_color_predicate_map[map_index];
          for (unsigned int i = 0; i < predicates.size(); ++i)
            {
              if (predicates[i](cell))
                {
                  /*
                   * create a pair predicate color and predicate id and add it
                   * to a map associated with each enriched cell
                   */
                  auto ret = cell_map.insert(
                    std::pair<unsigned int, unsigned int>(predicate_colors[i],
                                                          i));

                  AssertThrow(ret.second == 1,
                              ExcMessage(
                                "Only one enrichment function per color"));

                  color_list.insert(predicate_colors[i]);
                }
            }


          /*
           * check if color combination is already added.
           * If already added, set the active FE index based on
           * its index in the fe_sets. If the combination doesn't
           * exist, add the set to fe_sets and once again set the
           * active FE index as last index in fe_sets.
           *
           * Eg: if cell has color list {1,2} associated and
           * fe_sets = { {}, {1}, {2} } for now, we need to add a new set {1,2}
           * to fe_sets and a new active FE index 3 because 0 to 2 FE indices
           * are already taken by existing sets in fe_sets.
           */
          if (!color_list.empty())
            {
              const auto it =
                std::find(fe_sets.begin(), fe_sets.end(), color_list);
              // when entry is not found
              if (it == fe_sets.end())
                {
                  fe_sets.push_back(color_list);
                  cell->set_active_fe_index(fe_sets.size() - 1);
                }
              // when entry is found
              else
                {
                  cell->set_active_fe_index(std::distance(fe_sets.begin(), it));
                }
            }
          /*
           * map_index is used to identify cells and should be unique. The
           * index is stored in the material_id of the cell and hence
           * stays relevant even when the cells are refined (which is
           * why cell_id is not used).
           * Two distant cells can have the same set of colors but different
           * enrichment functions can be associated with any given
           * color. So, in order to figure which enrichment function
           * belongs to a color, we use a map that uses this index.
           */
          ++map_index;
        }


      /*
       * Treat interface between enriched cells specially,
       * until #1496 (https://github.com/dealii/dealii/issues/1496) is resolved.
       * Each time we build constraints at the interface between two different
       * FE_Enriched, we look for the least dominating FE of their common
       * subspace via hp::FECollection::find_dominating_fe_extended().
       * If we don't take further actions, we may find a dominating FE that is
       * too restrictive, i.e. enriched FE consisting of only FE_Nothing. New
       * elements needs to be added to FECollection object to help find the
       * correct enriched FE underlying the spaces in the adjacent cells. This
       * is done by creating an appropriate set in fe_sets and a call to the
       * function make_fe_collection_from_colored_enrichments at a later stage.
       *
       * Consider a domain with three predicates and hence with three different
       * enrichment functions. Let the enriched finite element of a cell with
       * enrichment functions 1 and 2 be represented by [1 1 0], with the last
       * entry as zero since the 3rd enrichment function is not relevant for
       * the cell. If the interface has enriched FE [1 0 1] and [0 1 1]
       * on adjacent cells, an enriched FE [0 0 1] should exist and is
       * found as the least dominating finite element for the two cells by
       * DoFTools::make_hanging_node_constraints, using the above mentioned
       * hp::FECollection functions. Denoting the FE set in adjacent cells as
       * {1,3} and {2,3}, this implies that an FE set {3} needs to be added!
       * Based on the predicate configuration, this may not be automatically
       * done without the following special treatment.
       */

      // loop through faces
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          const unsigned int           fe_index = cell->active_fe_index();
          const std::set<unsigned int> fe_set   = fe_sets.at(fe_index);
          for (const unsigned int face : GeometryInfo<dim>::face_indices())
            {
              // cell shouldn't be at the boundary and
              // neighboring cell is not already visited (to avoid visiting
              // same face twice). Note that the cells' material ids are
              // labeled according to their order in dof_handler previously.
              if (!cell->at_boundary(face) &&
                  cell->material_id() < cell->neighbor(face)->material_id())
                {
                  const auto nbr_fe_index =
                    cell->neighbor(face)->active_fe_index();

                  // find corresponding FE set
                  const auto nbr_fe_set = fe_sets.at(nbr_fe_index);

                  // find intersection of the FE sets: fe_set and nbr_fe_set
                  std::set<unsigned int> intersection_set;
                  std::set_intersection(
                    fe_set.begin(),
                    fe_set.end(),
                    nbr_fe_set.begin(),
                    nbr_fe_set.end(),
                    std::inserter(intersection_set, intersection_set.begin()));

                  // add only the new sets
                  if (!intersection_set.empty())
                    {
                      const auto it = std::find(fe_sets.begin(),
                                                fe_sets.end(),
                                                intersection_set);
                      // add the set if it is not found
                      if (it == fe_sets.end())
                        {
                          fe_sets.push_back(intersection_set);
                        }
                    }
                }
            }
        }
    }



    template <int dim, int spacedim>
    void
    make_colorwise_enrichment_functions(
      const unsigned int                                      n_colors,
      const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments,
      const std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &cellwise_color_predicate_map,
      std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments)
    {
      color_enrichments.clear();

      // Each color should be associated with a single enrichment function
      // called color enrichment function which calls the correct enrichment
      // function for a given cell.
      //
      // Assume that a cell has a active predicates with ids 4 (color = 1) and
      // 5 (color = 2). cellwise_color_predicate_map has this information,
      // provided we know the material id.
      //
      // The constructed color_enrichments is such that
      // color_enrichments[1](cell) will return return a pointer to
      // function with id=4, i.e. enrichments[4].
      // In other words, using the previously collected information in
      // this function we translate a vector of user provided enrichment
      // functions into a vector of functions suitable for FE_Enriched class.
      color_enrichments.resize(n_colors);
      for (unsigned int i = 0; i < n_colors; ++i)
        {
          color_enrichments[i] =
            [&, i](const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) {
              const unsigned int id = cell->material_id();

              /*
               * i'th color_enrichment function corresponds to (i+1)'th color.
               * Since FE_Enriched takes function pointers, we return a
               * function pointer.
               */
              return enrichments[cellwise_color_predicate_map.at(id).at(i + 1)]
                .get();
            };
        }
    }



    template <int dim, int spacedim>
    void
    make_fe_collection_from_colored_enrichments(
      const unsigned int                         n_colors,
      const std::vector<std::set<unsigned int>> &fe_sets,
      const std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
                                         &color_enrichments,
      const FiniteElement<dim, spacedim> &fe_base,
      const FiniteElement<dim, spacedim> &fe_enriched,
      const FE_Nothing<dim, spacedim>    &fe_nothing,
      hp::FECollection<dim, spacedim>    &fe_collection)
    {
      // define dummy function which is associated with FE_Nothing
      const std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>
        dummy_function =
          [=](const typename Triangulation<dim, spacedim>::cell_iterator &)
        -> const Function<spacedim> * {
        AssertThrow(false,
                    ExcMessage("Called enrichment function for FE_Nothing"));
        return nullptr;
      };


      // loop through color sets and create FE_enriched element for each
      // of them provided before calling this function, we have color
      // enrichment function associated with each color.
      for (const auto &fe_set : fe_sets)
        {
          std::vector<const FiniteElement<dim, spacedim> *> vec_fe_enriched(
            n_colors, &fe_nothing);
          std::vector<std::vector<std::function<const Function<spacedim> *(
            const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
            functions(n_colors, {dummy_function});

          for (const unsigned int color_id : fe_set)
            {
              // Given a color id, corresponding color enrichment
              // function is at index id-1 because color_enrichments are
              // indexed from zero and colors are indexed from 1.
              const unsigned int ind = color_id - 1;

              AssertIndexRange(ind, vec_fe_enriched.size());
              AssertIndexRange(ind, functions.size());
              AssertIndexRange(ind, color_enrichments.size());

              // Assume an active predicate colors {1,2} for a cell.
              // We then need to create a vector of FE enriched elements
              // with vec_fe_enriched[0] = vec_fe_enriched[1] = &fe_enriched
              // which can later be associated with enrichment functions.
              vec_fe_enriched[ind] = &fe_enriched;

              // color_set_id'th color function is (color_set_id-1)
              // element of color wise enrichments
              functions[ind][0] = color_enrichments[ind];
            }

          AssertDimension(vec_fe_enriched.size(), functions.size());

          FE_Enriched<dim, spacedim> fe_component(&fe_base,
                                                  vec_fe_enriched,
                                                  functions);
          fe_collection.push_back(fe_component);
        }
    }
  } // namespace internal



  template <int dim, int spacedim>
  Helper<dim, spacedim>::Helper(
    const FiniteElement<dim, spacedim>                     &fe_base,
    const FiniteElement<dim, spacedim>                     &fe_enriched,
    const std::vector<predicate_function<dim, spacedim>>   &predicates,
    const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments)
    : fe_base(fe_base)
    , fe_enriched(fe_enriched)
    , fe_nothing(fe_base.n_components(), true)
    , predicates(predicates)
    , enrichments(enrichments)
    , n_colors(numbers::invalid_unsigned_int)
  {
    AssertDimension(predicates.size(), enrichments.size());
    AssertDimension(fe_base.n_components(), fe_enriched.n_components());
    AssertThrow(predicates.size() > 0,
                ExcMessage("Number of predicates should be positive"));
  }



  template <int dim, int spacedim>
  const hp::FECollection<dim, spacedim> &
  Helper<dim, spacedim>::build_fe_collection(
    DoFHandler<dim, spacedim> &dof_handler)
  {
    // color the predicates based on connections between corresponding
    // subdomains
    n_colors =
      internal::color_predicates(dof_handler, predicates, predicate_colors);

    // create color maps and color list for each cell
    internal::set_cellwise_color_set_and_fe_index(dof_handler,
                                                  predicates,
                                                  predicate_colors,
                                                  cellwise_color_predicate_map,
                                                  fe_sets);
    // setup color wise enrichment functions
    // i'th function corresponds to (i+1) color!
    internal::make_colorwise_enrichment_functions<dim, spacedim>(
      n_colors, enrichments, cellwise_color_predicate_map, color_enrichments);

    // make FE_Collection
    internal::make_fe_collection_from_colored_enrichments(n_colors,
                                                          fe_sets,
                                                          color_enrichments,
                                                          fe_base,
                                                          fe_enriched,
                                                          fe_nothing,
                                                          fe_collection);

    return fe_collection;
  }
} // namespace ColorEnriched


// explicit instantiations
#include "fe/fe_enriched.inst"

DEAL_II_NAMESPACE_CLOSE
