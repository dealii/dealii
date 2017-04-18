// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/fe/fe_enriched.h>

#include <deal.II/fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  /**
   * Auxiliary function to create multiplicity vector from input enrichment functions.
   */
  template <typename T>
  std::vector<unsigned int>
  build_multiplicities(const std::vector<std::vector<T > > &functions )
  {
    std::vector<unsigned int> multiplicities;
    multiplicities.push_back(1); // the first one is non-enriched FE
    for (unsigned int i = 0; i < functions.size(); i++)
      multiplicities.push_back(functions[i].size());

    return multiplicities;
  }


  /**
   * Auxiliary function to build FiniteElement's vector
   */
  template <int dim, int spacedim>
  std::vector< const FiniteElement< dim, spacedim > * >
  build_fes(const FiniteElement<dim,spacedim> *fe_base,
            const std::vector<const FiniteElement<dim,spacedim> * > &fe_enriched)
  {
    std::vector< const FiniteElement< dim, spacedim > * > fes;
    fes.push_back(fe_base);
    for (unsigned int i = 0; i < fe_enriched.size(); i++)
      fes.push_back(fe_enriched[i]);

    return fes;
  }


  /**
   * Auxiliary function which check consistency of the input parameters.
   * Returns true if everything is ok.
   */
  template <int dim, int spacedim>
  bool
  consistency_check (const std::vector< const FiniteElement< dim, spacedim > * > &fes,
                     const std::vector< unsigned int > &multiplicities,
                     const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > &functions)
  {
    AssertThrow(fes.size() > 0,
                ExcMessage("FEs size should be >=1"));
    AssertThrow(fes.size() == multiplicities.size(),
                ExcMessage("FEs and multiplicities should have the same size"));

    AssertThrow (functions.size() == fes.size() - 1,
                 ExcDimensionMismatch(functions.size(), fes.size()-1));

    AssertThrow(multiplicities[0] == 1,
                ExcMessage("First multiplicity should be 1"));

    const unsigned int n_comp_base = fes[0]->n_components();

    // start from fe=1 as 0th is always non-enriched FE.
    for (unsigned int fe=1; fe < fes.size(); fe++)
      {
        const FE_Nothing<dim> *fe_nothing = dynamic_cast<const FE_Nothing<dim>*>(fes[fe]);
        if (fe_nothing)
          AssertThrow (fe_nothing->is_dominating(),
                       ExcMessage("Only dominating FE_Nothing can be used in FE_Enriched"));

        AssertThrow (fes[fe]->n_components() == n_comp_base,
                     ExcMessage("All elements must have the same number of components"));
      }
    return true;
  }


  /**
   * Auxiliary function which determines whether the FiniteElement will be enriched.
   */
  template <int dim, int spacedim>
  bool
  check_if_enriched (const std::vector< const FiniteElement< dim, spacedim > * > &fes)
  {
    // start from fe=1 as 0th is always non-enriched FE.
    for (unsigned int fe=1; fe < fes.size(); fe++)
      if (dynamic_cast<const FE_Nothing<dim>*>(fes[fe]) == nullptr)
        // this is not FE_Nothing => there will be enrichment
        return true;

    return false;
  }
}


template <int dim, int spacedim>
FE_Enriched<dim,spacedim>::FE_Enriched (const FiniteElement<dim,spacedim> &fe_base)
  :
  FE_Enriched<dim,spacedim>(fe_base,
                            FE_Nothing<dim,spacedim>(fe_base.n_components(),true),
                            nullptr)
{
}


template <int dim, int spacedim>
FE_Enriched<dim,spacedim>::FE_Enriched (const FiniteElement<dim,spacedim> &fe_base,
                                        const FiniteElement<dim,spacedim> &fe_enriched,
                                        const Function<spacedim>      *enrichment_function)
  :
  FE_Enriched<dim,spacedim>
  (&fe_base,
   std::vector<const FiniteElement<dim,spacedim>*>(1, &fe_enriched),
   std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > >
   (1,
    std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > >
    (1,
     [=] (const typename Triangulation<dim, spacedim>::cell_iterator &) -> const Function<spacedim> *
{
  return enrichment_function;
})
 )
)
{}


template <int dim, int spacedim>
FE_Enriched<dim,spacedim>::FE_Enriched (const FiniteElement<dim,spacedim> *fe_base,
                                        const std::vector<const FiniteElement<dim,spacedim> * > &fe_enriched,
                                        const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > &functions)
  :
  FE_Enriched<dim,spacedim> (build_fes(fe_base,fe_enriched),
                             build_multiplicities(functions),
                             functions)
{}


template <int dim, int spacedim>
FE_Enriched<dim,spacedim>::FE_Enriched (const std::vector< const FiniteElement< dim, spacedim > * > &fes,
                                        const std::vector< unsigned int > &multiplicities,
                                        const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > > &functions)
  :
  FiniteElement<dim,spacedim> (FETools::Compositing::multiply_dof_numbers(fes,multiplicities,false),
                               FETools::Compositing::compute_restriction_is_additive_flags(fes,multiplicities),
                               FETools::Compositing::compute_nonzero_components(fes,multiplicities,false)),
  enrichments(functions),
  is_enriched(check_if_enriched(fes)),
  fe_system(fes,multiplicities)
{
  // descriptive error are thrown within the function.
  Assert(consistency_check(fes,multiplicities,functions),
         ExcInternalError());

  initialize(fes, multiplicities);

  // resize to be consistent with all FEs used to construct the FE_Enriched,
  // even though we will never use the 0th element.
  base_no_mult_local_enriched_dofs.resize(fes.size());
  for (unsigned int fe=1; fe < fes.size(); fe++)
    base_no_mult_local_enriched_dofs[fe].resize(multiplicities[fe]);

  Assert (base_no_mult_local_enriched_dofs.size() == this->n_base_elements(),
          ExcDimensionMismatch(base_no_mult_local_enriched_dofs.size(),
                               this->n_base_elements()));

  // build the map: (base_no, base_m) -> vector of local element DoFs
  for (unsigned int system_index=0; system_index<this->dofs_per_cell;
       ++system_index)
    {
      const unsigned int base_no = this->system_to_base_table[system_index].first.first;
      if (base_no == 0) // 0th is always non-enriched FE
        continue;

      const unsigned int base_m  = this->system_to_base_table[system_index].first.second;

      Assert (base_m < base_no_mult_local_enriched_dofs[base_no].size(),
              ExcMessage("Size mismatch for base_no_mult_local_enriched_dofs: "
                         "base_index = " + std::to_string(this->system_to_base_table[system_index].second) +
                         "; base_no = " + std::to_string(base_no) +
                         "; base_m = " + std::to_string(base_m) +
                         "; system_index = " + std::to_string(system_index)));

      Assert (base_m < base_no_mult_local_enriched_dofs[base_no].size(),
              ExcDimensionMismatch(base_m,
                                   base_no_mult_local_enriched_dofs[base_no].size()));

      base_no_mult_local_enriched_dofs[base_no][base_m].push_back(system_index);
    }

  // make sure that local_enriched_dofs.size() is correct, that is equals to DoFs
  // per cell of the corresponding FE.
  for (unsigned int base_no = 1; base_no < base_no_mult_local_enriched_dofs.size(); base_no++)
    {
      for (unsigned int m=0; m < base_no_mult_local_enriched_dofs[base_no].size(); m++)
        Assert ( base_no_mult_local_enriched_dofs[base_no][m].size() == fes[base_no]->dofs_per_cell,
                 ExcDimensionMismatch(base_no_mult_local_enriched_dofs[base_no][m].size(),
                                      fes[base_no]->dofs_per_cell));
    }
}


template <int dim, int spacedim>
const std::vector<std::vector<std::function<const Function<spacedim> *(const typename Triangulation<dim, spacedim>::cell_iterator &) > > >
FE_Enriched<dim,spacedim>::get_enrichments() const
{
  return enrichments;
}


template <int dim, int spacedim>
double
FE_Enriched<dim,spacedim>::shape_value(const unsigned int   i,
                                       const Point< dim > &p) const
{
  Assert(!is_enriched, ExcMessage("For enriched finite elements shape_value() can not be defined on the reference element."));
  return fe_system.shape_value(i,p);
}


template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Enriched<dim,spacedim>::clone() const
{
  std::vector< const FiniteElement< dim, spacedim > * > fes;
  std::vector< unsigned int > multiplicities;

  for (unsigned int i=0; i<this->n_base_elements(); i++)
    {
      fes.push_back( & base_element(i) );
      multiplicities.push_back(this->element_multiplicity(i) );
    }

  return new FE_Enriched<dim,spacedim>(fes, multiplicities, get_enrichments());
}


template <int dim, int spacedim>
UpdateFlags
FE_Enriched<dim,spacedim>::requires_update_flags (const UpdateFlags flags) const
{
  UpdateFlags out = fe_system.requires_update_flags(flags);

  if (is_enriched)
    {
      // if we ask for values or gradients, then we would need quadrature points
      if (flags & (update_values | update_gradients))
        out |= update_quadrature_points;

      // if need gradients, add update_values due to product rule
      if (out & update_gradients)
        out |= update_values;
    }

  Assert (!(flags & update_3rd_derivatives),
          ExcNotImplemented());

  return out;
}


template <int dim, int spacedim>
template <int dim_1>
typename FiniteElement<dim,spacedim>::InternalDataBase *
FE_Enriched<dim,spacedim>::setup_data (std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase> fes_data,
                                       const UpdateFlags      flags,
                                       const Quadrature<dim_1> &quadrature) const
{
  Assert ((dynamic_cast<typename FESystem<dim,spacedim>::InternalData *> (fes_data.get()) != nullptr),
          ExcInternalError());
  typename FESystem<dim,spacedim>::InternalData *data_fesystem =
    static_cast<typename FESystem<dim,spacedim>::InternalData *> (fes_data.get());

  // FESystem::InternalData will be aggregated (owned) by
  // our InternalData.
  fes_data.release();
  InternalData *data = new InternalData(std::unique_ptr<typename FESystem<dim,spacedim>::InternalData>(data_fesystem));

  // copy update_each from FESystem data:
  data->update_each = data_fesystem->update_each;

  // resize cache array according to requested flags
  data->enrichment.resize(this->n_base_elements());

  const unsigned int n_q_points = quadrature.size();

  for (unsigned int base=0; base < this->n_base_elements(); ++base)
    {
      data->enrichment[base].resize(this->element_multiplicity(base));
      for (unsigned int m = 0; m < this->element_multiplicity(base); ++m)
        {
          if (flags & update_values)
            data->enrichment[base][m].values.resize(n_q_points);

          if (flags & update_gradients)
            data->enrichment[base][m].gradients.resize(n_q_points);

          if (flags & update_hessians)
            data->enrichment[base][m].hessians.resize(n_q_points);
        }
    }

  return data;
}


template <int dim, int spacedim>
typename FiniteElement<dim,spacedim>::InternalDataBase *
FE_Enriched<dim,spacedim>::get_face_data (const UpdateFlags      update_flags,
                                          const Mapping<dim,spacedim>    &mapping,
                                          const Quadrature<dim-1> &quadrature,
                                          internal::FEValues::FiniteElementRelatedData< dim, spacedim >        &output_data) const
{
  return setup_data(std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>(fe_system.get_face_data(update_flags,mapping,quadrature,output_data)),
                    update_flags,
                    quadrature);
}


template <int dim, int spacedim>
typename FiniteElement<dim,spacedim>::InternalDataBase *
FE_Enriched<dim,spacedim>::get_subface_data (const UpdateFlags      update_flags,
                                             const Mapping<dim,spacedim>    &mapping,
                                             const Quadrature<dim-1> &quadrature,
                                             dealii::internal::FEValues::FiniteElementRelatedData<dim, spacedim> &output_data) const
{
  return setup_data(std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>(fe_system.get_subface_data(update_flags,mapping,quadrature,output_data)),
                    update_flags,
                    quadrature);
}


template <int dim, int spacedim>
typename FiniteElement<dim,spacedim>::InternalDataBase *
FE_Enriched<dim,spacedim>::get_data (const UpdateFlags      flags,
                                     const Mapping<dim,spacedim>    &mapping,
                                     const Quadrature<dim> &quadrature,
                                     internal::FEValues::FiniteElementRelatedData< dim, spacedim >   &output_data) const
{
  return setup_data(std::unique_ptr<typename FiniteElement<dim,spacedim>::InternalDataBase>(fe_system.get_data(flags,mapping,quadrature,output_data)),
                    flags,
                    quadrature);
}


template <int dim, int spacedim>
void FE_Enriched<dim,spacedim>::initialize (const std::vector<const FiniteElement<dim,spacedim>*> &fes,
                                            const std::vector<unsigned int> &multiplicities)
{
  Assert (fes.size() == multiplicities.size(),
          ExcDimensionMismatch (fes.size(), multiplicities.size()) );

  // Note that we need to skip every fe with multiplicity 0 in the following block of code
  this->base_to_block_indices.reinit(0, 0);

  for (unsigned int i=0; i<fes.size(); i++)
    if (multiplicities[i]>0)
      this->base_to_block_indices.push_back( multiplicities[i] );

  {
    // If the system is not primitive, these have not been initialized by
    // FiniteElement
    this->system_to_component_table.resize(this->dofs_per_cell);
    this->face_system_to_component_table.resize(this->dofs_per_face);

    FETools::Compositing::build_cell_tables(this->system_to_base_table,
                                            this->system_to_component_table,
                                            this->component_to_base_table,
                                            *this,
                                            false);

    FETools::Compositing::build_face_tables(this->face_system_to_base_table,
                                            this->face_system_to_component_table,
                                            *this,
                                            false);
  }

  // restriction and prolongation matrices are built on demand

  // now set up the interface constraints for h-refinement.
  // take them from fe_system:
  this->interface_constraints = fe_system.interface_constraints;

  // if we just wrap another FE (i.e. use FE_Nothing as a second FE)
  // then it makes sense to have support points.
  // However, functions like interpolate_boundary_values() need all FEs inside
  // FECollection to be able to provide support points irrespectively whether
  // this FE sits on the boundary or not. Thus for moment just copy support
  // points from fe system:
  {
    this->unit_support_points      = fe_system.unit_support_points;
    this->unit_face_support_points = fe_system.unit_face_support_points;
  }

  // take adjust_quad_dof_index_for_face_orientation_table from FESystem:
  {
    this->adjust_line_dof_index_for_line_orientation_table = fe_system.adjust_line_dof_index_for_line_orientation_table;
  }
}


template <int dim, int spacedim>
std::string
FE_Enriched<dim,spacedim>::get_name () const
{
  std::ostringstream namebuf;

  namebuf << "FE_Enriched<"
          << Utilities::dim_string(dim,spacedim)
          << ">[";
  for (unsigned int i=0; i< this->n_base_elements(); ++i)
    {
      namebuf << base_element(i).get_name();
      if (this->element_multiplicity(i) != 1)
        namebuf << '^' << this->element_multiplicity(i);
      if (i != this->n_base_elements()-1)
        namebuf << '-';
    }
  namebuf << ']';

  return namebuf.str();
}


template <int dim, int spacedim>
const FiniteElement<dim,spacedim> &
FE_Enriched<dim,spacedim>::base_element (const unsigned int index) const
{
  return fe_system.base_element(index);
}


template <int dim, int spacedim>
void
FE_Enriched<dim,spacedim>::fill_fe_values (const typename Triangulation< dim, spacedim >::cell_iterator &cell,
                                           const CellSimilarity::Similarity cell_similarity,
                                           const Quadrature< dim > &quadrature,
                                           const Mapping< dim, spacedim > &mapping,
                                           const typename Mapping< dim, spacedim >::InternalDataBase &mapping_internal,
                                           const dealii::internal::FEValues::MappingRelatedData< dim, spacedim > &mapping_data,
                                           const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
                                           internal::FEValues::FiniteElementRelatedData< dim, spacedim > &output_data
                                          ) const
{
  Assert (dynamic_cast<const InternalData *> (&fe_internal) != nullptr,
          ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  // call FESystem's method to fill everything without enrichment function
  fe_system.fill_fe_values(cell,
                           cell_similarity,
                           quadrature,
                           mapping,
                           mapping_internal,
                           mapping_data,
                           *fe_data.fesystem_data,
                           output_data);

  if (is_enriched)
    multiply_by_enrichment(quadrature,
                           fe_data,
                           mapping_data,
                           cell,
                           output_data);
}


template <int dim, int spacedim>
void
FE_Enriched<dim,spacedim>::fill_fe_face_values
(const typename Triangulation< dim, spacedim >::cell_iterator &cell,
 const unsigned int face_no,
 const Quadrature< dim-1 > &quadrature,
 const Mapping< dim, spacedim > &mapping,
 const typename Mapping< dim, spacedim >::InternalDataBase &mapping_internal,
 const dealii::internal::FEValues::MappingRelatedData< dim, spacedim > &mapping_data,
 const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
 internal::FEValues::FiniteElementRelatedData< dim, spacedim > &output_data
) const
{
  Assert (dynamic_cast<const InternalData *> (&fe_internal) != nullptr,
          ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  // call FESystem's method to fill everything without enrichment function
  fe_system.fill_fe_face_values(cell,
                                face_no,
                                quadrature,
                                mapping,
                                mapping_internal,
                                mapping_data,
                                *fe_data.fesystem_data,
                                output_data);

  if (is_enriched)
    multiply_by_enrichment(quadrature,
                           fe_data,
                           mapping_data,
                           cell,
                           output_data);
}


template <int dim, int spacedim>
void
FE_Enriched<dim,spacedim>::fill_fe_subface_values
(const typename Triangulation< dim, spacedim >::cell_iterator &cell,
 const unsigned int face_no,
 const unsigned int sub_no,
 const Quadrature< dim-1 > &quadrature,
 const Mapping< dim, spacedim > &mapping,
 const typename Mapping< dim, spacedim >::InternalDataBase &mapping_internal,
 const dealii::internal::FEValues::MappingRelatedData< dim, spacedim > &mapping_data,
 const typename FiniteElement<dim,spacedim>::InternalDataBase &fe_internal,
 internal::FEValues::FiniteElementRelatedData< dim, spacedim > &output_data
) const
{
  Assert (dynamic_cast<const InternalData *> (&fe_internal) != nullptr,
          ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  // call FESystem's method to fill everything without enrichment function
  fe_system.fill_fe_subface_values(cell,
                                   face_no,
                                   sub_no,
                                   quadrature,
                                   mapping,
                                   mapping_internal,
                                   mapping_data,
                                   *fe_data.fesystem_data,
                                   output_data);

  if (is_enriched)
    multiply_by_enrichment(quadrature,
                           fe_data,
                           mapping_data,
                           cell,
                           output_data);

}


template <int dim, int spacedim>
template <int dim_1>
void
FE_Enriched<dim,spacedim>::multiply_by_enrichment
(const Quadrature<dim_1> &quadrature,
 const InternalData &fe_data,
 const internal::FEValues::MappingRelatedData<dim,spacedim> &mapping_data,
 const typename Triangulation< dim, spacedim >::cell_iterator &cell,
 internal::FEValues::FiniteElementRelatedData<dim,spacedim> &output_data) const
{
  // mapping_data will contain quadrature points on the real element.
  // fe_internal is needed to get update flags
  // finally, output_data should store all the results we need.

  // Either dim_1==dim
  // (fill_fe_values) or dim_1==dim-1
  // (fill_fe_(sub)face_values)
  Assert(dim_1==dim || dim_1==dim-1, ExcInternalError());
  const UpdateFlags flags = fe_data.update_each;

  const unsigned int n_q_points  = quadrature.size();

  // First, populate output_data object (that shall hold everything requested such
  // as shape value, gradients, hessians, etc) from each base element.
  // That is almost identical to FESystem::compute_fill_one_base(),
  // the difference being that we do it irrespectively of cell_similarity
  // and use base_fe_data.update_flags

  // TODO: do we need it only for dim_1 == dim (i.e. fill_fe_values)?
  if (dim_1 == dim)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); base_no++)
      {
        const FiniteElement<dim,spacedim> &
        base_fe      = base_element(base_no);
        typename FiniteElement<dim,spacedim>::InternalDataBase &
        base_fe_data = fe_data.get_fe_data(base_no);
        internal::FEValues::FiniteElementRelatedData<dim,spacedim> &
        base_data    = fe_data.get_fe_output_object(base_no);

        const UpdateFlags base_flags = base_fe_data.update_each;

        for (unsigned int system_index=0; system_index<this->dofs_per_cell;
             ++system_index)
          if (this->system_to_base_table[system_index].first.first == base_no)
            {
              const unsigned int
              base_index = this->system_to_base_table[system_index].second;
              Assert (base_index<base_fe.dofs_per_cell, ExcInternalError());

              // now copy. if the shape function is primitive, then there
              // is only one value to be copied, but for non-primitive
              // elements, there might be more values to be copied
              //
              // so, find out from which index to take this one value, and
              // to which index to put
              unsigned int out_index = 0;
              for (unsigned int i=0; i<system_index; ++i)
                out_index += this->n_nonzero_components(i);
              unsigned int in_index = 0;
              for (unsigned int i=0; i<base_index; ++i)
                in_index += base_fe.n_nonzero_components(i);

              // then loop over the number of components to be copied
              Assert (this->n_nonzero_components(system_index) ==
                      base_fe.n_nonzero_components(base_index),
                      ExcInternalError());
              for (unsigned int s=0; s<this->n_nonzero_components(system_index); ++s)
                {
                  if (base_flags & update_values)
                    for (unsigned int q=0; q<n_q_points; ++q)
                      output_data.shape_values[out_index+s][q] =
                        base_data.shape_values(in_index+s,q);

                  if (base_flags & update_gradients)
                    for (unsigned int q=0; q<n_q_points; ++q)
                      output_data.shape_gradients[out_index+s][q] =
                        base_data.shape_gradients[in_index+s][q];

                  if (base_flags & update_hessians)
                    for (unsigned int q=0; q<n_q_points; ++q)
                      output_data.shape_hessians[out_index+s][q] =
                        base_data.shape_hessians[in_index+s][q];
                }
            }
      }

  Assert (base_no_mult_local_enriched_dofs.size() == fe_data.enrichment.size(),
          ExcDimensionMismatch(base_no_mult_local_enriched_dofs.size(),
                               fe_data.enrichment.size()));
  // calculate hessians, gradients and values for each function
  for (unsigned int base_no = 1; base_no < this->n_base_elements(); base_no++)
    {
      Assert (base_no_mult_local_enriched_dofs[base_no].size() == fe_data.enrichment[base_no].size(),
              ExcDimensionMismatch(base_no_mult_local_enriched_dofs[base_no].size(),
                                   fe_data.enrichment[base_no].size()));
      for (unsigned int m=0; m < base_no_mult_local_enriched_dofs[base_no].size(); m++)
        {
          Assert (enrichments[base_no-1][m](cell) != nullptr,
                  ExcMessage("The pointer to the enrichment function is NULL"));

          Assert (enrichments[base_no-1][m](cell)->n_components == 1,
                  ExcMessage("Only scalar-valued enrichment functions are allowed"));

          if (flags & update_hessians)
            {
              Assert (fe_data.enrichment[base_no][m].hessians.size() == n_q_points,
                      ExcDimensionMismatch(fe_data.enrichment[base_no][m].hessians.size(),
                                           n_q_points));
              for (unsigned int q=0; q<n_q_points; q++)
                fe_data.enrichment[base_no][m].hessians[q]  = enrichments[base_no-1][m](cell)->hessian (mapping_data.quadrature_points[q]);
            }

          if (flags & update_gradients)
            {
              Assert (fe_data.enrichment[base_no][m].gradients.size() == n_q_points,
                      ExcDimensionMismatch(fe_data.enrichment[base_no][m].gradients.size(),
                                           n_q_points));
              for (unsigned int q=0; q<n_q_points; q++)
                fe_data.enrichment[base_no][m].gradients[q] = enrichments[base_no-1][m](cell)->gradient(mapping_data.quadrature_points[q]);
            }

          if (flags & update_values)
            {
              Assert (fe_data.enrichment[base_no][m].values.size() == n_q_points,
                      ExcDimensionMismatch(fe_data.enrichment[base_no][m].values.size(),
                                           n_q_points));
              for (unsigned int q=0; q<n_q_points; q++)
                fe_data.enrichment[base_no][m].values[q]    = enrichments[base_no-1][m](cell)->value   (mapping_data.quadrature_points[q]);
            }
        }
    }

  // Finally, update the standard data stored in output_data
  // by expanding the product rule for enrichment function.
  // note that the order if important, namely
  // output_data.shape_XYZ contains values of standard FEM and we overwrite
  // it with the updated one in the following order: hessians -> gradients -> values
  if (flags & update_hessians)
    {
      for (unsigned int base_no = 1; base_no < this->n_base_elements(); base_no++)
        {
          for (unsigned int m=0; m < base_no_mult_local_enriched_dofs[base_no].size(); m++)
            for (unsigned int i=0; i < base_no_mult_local_enriched_dofs[base_no][m].size(); i++)
              {
                const unsigned int enriched_dof = base_no_mult_local_enriched_dofs[base_no][m][i];
                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    const Tensor<2, spacedim> grad_grad =  outer_product(output_data.shape_gradients[enriched_dof][q],
                                                                         fe_data.enrichment[base_no][m].gradients[q]);
                    const Tensor<2,spacedim,double> sym_grad_grad = symmetrize (grad_grad) * 2.0; // symmetrize does [s+s^T]/2

                    output_data.shape_hessians   [enriched_dof][q]*= fe_data.enrichment[base_no][m].values[q];
                    output_data.shape_hessians   [enriched_dof][q]+=
                      sym_grad_grad +
                      output_data.shape_values   [enriched_dof][q] * fe_data.enrichment[base_no][m].hessians[q];
                  }
              }
        }
    }

  if (flags & update_gradients)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); base_no++)
      {
        for (unsigned int m=0; m < base_no_mult_local_enriched_dofs[base_no].size(); m++)
          for (unsigned int i=0; i < base_no_mult_local_enriched_dofs[base_no][m].size(); i++)
            {
              const unsigned int enriched_dof = base_no_mult_local_enriched_dofs[base_no][m][i];
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  output_data.shape_gradients[enriched_dof][q]*=  fe_data.enrichment[base_no][m].values[q];
                  output_data.shape_gradients[enriched_dof][q]+=
                    output_data.shape_values [enriched_dof][q] *  fe_data.enrichment[base_no][m].gradients[q];
                }
            }
      }

  if (flags & update_values)
    for (unsigned int base_no = 1; base_no < this->n_base_elements(); base_no++)
      {
        for (unsigned int m=0; m < base_no_mult_local_enriched_dofs[base_no].size(); m++)
          for (unsigned int i=0; i < base_no_mult_local_enriched_dofs[base_no][m].size(); i++)
            {
              const unsigned int enriched_dof = base_no_mult_local_enriched_dofs[base_no][m][i];
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  output_data.shape_values[enriched_dof][q] *= fe_data.enrichment[base_no][m].values[q];
                }
            }
      }
}


template <int dim, int spacedim>
const FESystem<dim,spacedim> &
FE_Enriched<dim,spacedim>::
get_fe_system() const
{
  return fe_system;
}


template <int dim, int spacedim>
bool
FE_Enriched<dim,spacedim>::
hp_constraints_are_implemented () const
{
  return true;
}


template <int dim, int spacedim>
void
FE_Enriched<dim,spacedim>::
get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                               FullMatrix<double>                &matrix) const
{
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&source))
    {
      fe_system.get_face_interpolation_matrix(fe_enr_other->get_fe_system(),
                                              matrix);
    }
  else
    {
      typedef FiniteElement<dim,spacedim> FEL;
      AssertThrow(false,typename FEL::ExcInterpolationNotImplemented());
    }
}


template <int dim, int spacedim>
void
FE_Enriched<dim,spacedim>::
get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
                                  const unsigned int        subface,
                                  FullMatrix<double>       &matrix) const
{
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&source))
    {
      fe_system.get_subface_interpolation_matrix(fe_enr_other->get_fe_system(),
                                                 subface,
                                                 matrix);
    }
  else
    {
      typedef FiniteElement<dim,spacedim> FEL;
      AssertThrow(false,typename FEL::ExcInterpolationNotImplemented());
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Enriched<dim,spacedim>::
hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&fe_other))
    {
      return fe_system.hp_vertex_dof_identities(fe_enr_other->get_fe_system());
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> >();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Enriched<dim,spacedim>::
hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&fe_other))
    {
      return fe_system.hp_line_dof_identities(fe_enr_other->get_fe_system());
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> >();
    }
}


template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Enriched<dim,spacedim>::
hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const
{
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&fe_other))
    {
      return fe_system.hp_quad_dof_identities(fe_enr_other->get_fe_system());
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int> >();
    }
}


template <int dim, int spacedim>
FiniteElementDomination::Domination
FE_Enriched<dim,spacedim>::
compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const
{
  // need to decide which element constrain another.
  // for example Q(2) dominate Q(4) and thus some DoFs of Q(4) will be constrained.
  // If we have Q(2) and Q(4)+POU, then it's clear that Q(2) dominates,
  // namely our DoFs will be constrained to make field continuous.
  // However, we need to check for situations like Q(4) vs Q(2)+POU.
  // In that case the domination for the underlying FEs should be the other way,
  // but this implies that we can't constrain POU dofs to make the field continuous.
  // In that case, through an error

  // if it's also enriched, do domination based on each one's FESystem
  if (const FE_Enriched<dim,spacedim> *fe_enr_other
      = dynamic_cast<const FE_Enriched<dim,spacedim>*>(&fe_other))
    {
      return fe_system.compare_for_face_domination(fe_enr_other->get_fe_system());
    }
  else
    {
      Assert (false, ExcNotImplemented());
      return FiniteElementDomination::neither_element_dominates;
    }
}

template <int dim, int spacedim>
const FullMatrix<double> &
FE_Enriched<dim,spacedim>::get_prolongation_matrix (const unsigned int child,
                                                    const RefinementCase<dim> &refinement_case) const
{
  return fe_system.get_prolongation_matrix(child, refinement_case);
}


template <int dim, int spacedim>
const FullMatrix<double> &
FE_Enriched<dim,spacedim>::get_restriction_matrix (const unsigned int child,
                                                   const RefinementCase<dim> &refinement_case) const
{
  return fe_system.get_restriction_matrix(child, refinement_case);
}


/* ----------------------- FESystem::InternalData ------------------- */


template <int dim, int spacedim>
FE_Enriched<dim,spacedim>::InternalData::InternalData(std::unique_ptr<typename FESystem<dim,spacedim>::InternalData> fesystem_data)
  :
  fesystem_data(std::move(fesystem_data))
{}


template <int dim, int spacedim>
typename FiniteElement<dim,spacedim>::InternalDataBase &
FE_Enriched<dim,spacedim>::
InternalData::get_fe_data (const unsigned int base_no) const
{
  return fesystem_data->get_fe_data(base_no);
}


template <int dim, int spacedim>
internal::FEValues::FiniteElementRelatedData<dim,spacedim> &
FE_Enriched<dim,spacedim>::
InternalData::get_fe_output_object (const unsigned int base_no) const
{
  return fesystem_data->get_fe_output_object(base_no);
}


// explicit instantiations
#include "fe_enriched.inst"

DEAL_II_NAMESPACE_CLOSE
