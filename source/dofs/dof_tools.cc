// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>

#include <algorithm>
#include <numeric>

DEAL_II_NAMESPACE_OPEN



namespace DoFTools
{
  namespace internal
  {
    // return an array that for each dof on the reference cell
    // lists the corresponding vector component.
    //
    // if an element is non-primitive then we assign to each degree of freedom
    // the following component:
    // - if the nonzero components that belong to a shape function are not
    //   selected in the component_mask, then the shape function is assigned
    //   to the first nonzero vector component that corresponds to this
    //   shape function
    // - otherwise, the shape function is assigned the first component selected
    //   in the component_mask that corresponds to this shape function
    template <int dim, int spacedim>
    std::vector<unsigned char>
    get_local_component_association (const FiniteElement<dim,spacedim>  &fe,
                                     const ComponentMask        &component_mask)
    {
      std::vector<unsigned char> local_component_association (fe.dofs_per_cell,
                                                              (unsigned char)(-1));

      // compute the component each local dof belongs to.
      // if the shape function is primitive, then this
      // is simple and we can just associate it with
      // what system_to_component_index gives us
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        if (fe.is_primitive(i))
          local_component_association[i] =
            fe.system_to_component_index(i).first;
        else
          // if the shape function is not primitive, then either use the
          // component of the first nonzero component corresponding
          // to this shape function (if the component is not specified
          // in the component_mask), or the first component of this block
          // that is listed in the component_mask (if the block this
          // component corresponds to is indeed specified in the component
          // mask)
          {
            const unsigned int first_comp =
              fe.get_nonzero_components(i).first_selected_component();

            if ((fe.get_nonzero_components(i)
                 &
                 component_mask).n_selected_components(fe.n_components()) == 0)
              local_component_association[i] = first_comp;
            else
              // pick the component selected. we know from the previous 'if'
              // that within the components that are nonzero for this
              // shape function there must be at least one for which the
              // mask is true, so we will for sure run into the break()
              // at one point
              for (unsigned int c=first_comp; c<fe.n_components(); ++c)
                if (component_mask[c] == true)
                  {
                    local_component_association[i] = c;
                    break;
                  }
          }

      Assert (std::find (local_component_association.begin(),
                         local_component_association.end(),
                         (unsigned char)(-1))
              ==
              local_component_association.end(),
              ExcInternalError());

      return local_component_association;
    }


    // this internal function assigns to each dof the respective component
    // of the vector system.
    //
    // the output array dofs_by_component lists for each dof the
    // corresponding vector component. if the DoFHandler is based on a
    // parallel distributed triangulation then the output array is index by
    // dof.locally_owned_dofs().index_within_set(indices[i])
    //
    // if an element is non-primitive then we assign to each degree of
    // freedom the following component:
    // - if the nonzero components that belong to a shape function are not
    //   selected in the component_mask, then the shape function is assigned
    //   to the first nonzero vector component that corresponds to this
    //   shape function
    // - otherwise, the shape function is assigned the first component selected
    //   in the component_mask that corresponds to this shape function
    template <typename DoFHandlerType>
    void
    get_component_association (const DoFHandlerType       &dof,
                               const ComponentMask        &component_mask,
                               std::vector<unsigned char> &dofs_by_component)
    {
      const dealii::hp::FECollection<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
      fe_collection (dof.get_fe());
      Assert (fe_collection.n_components() < 256, ExcNotImplemented());
      Assert (dofs_by_component.size() == dof.n_locally_owned_dofs(),
              ExcDimensionMismatch(dofs_by_component.size(),
                                   dof.n_locally_owned_dofs()));

      // next set up a table for the degrees of freedom on each of the
      // cells (regardless of the fact whether it is listed in the
      // component_select argument or not)
      //
      // for each element 'f' of the FECollection,
      // local_component_association[f][d] then returns the vector
      // component that degree of freedom 'd' belongs to
      std::vector<std::vector<unsigned char> >
      local_component_association (fe_collection.size());
      for (unsigned int f=0; f<fe_collection.size(); ++f)
        {
          const FiniteElement<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &fe =
            fe_collection[f];
          local_component_association[f]
            = get_local_component_association (fe, component_mask);
        }

      // then loop over all cells and do the work
      std::vector<types::global_dof_index> indices;
      for (typename DoFHandlerType::active_cell_iterator c=dof.begin_active();
           c!=dof.end(); ++ c)
        if (c->is_locally_owned())
          {
            const unsigned int fe_index = c->active_fe_index();
            const unsigned int dofs_per_cell = c->get_fe().dofs_per_cell;
            indices.resize(dofs_per_cell);
            c->get_dof_indices(indices);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              if (dof.locally_owned_dofs().is_element(indices[i]))
                dofs_by_component[dof.locally_owned_dofs().index_within_set(indices[i])]
                  = local_component_association[fe_index][i];
          }
    }


    // this is the function corresponding to the one above but working on
    // blocks instead of components.
    //
    // the output array dofs_by_block lists for each dof the corresponding
    // vector block. if the DoFHandler is based on a parallel distributed
    // triangulation then the output array is index by
    // dof.locally_owned_dofs().index_within_set(indices[i])
    template <typename DoFHandlerType>
    inline
    void
    get_block_association (const DoFHandlerType       &dof,
                           std::vector<unsigned char> &dofs_by_block)
    {
      const dealii::hp::FECollection<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
      fe_collection (dof.get_fe());
      Assert (fe_collection.n_components() < 256, ExcNotImplemented());
      Assert (dofs_by_block.size() == dof.n_locally_owned_dofs(),
              ExcDimensionMismatch(dofs_by_block.size(),
                                   dof.n_locally_owned_dofs()));

      // next set up a table for the degrees of freedom on each of the
      // cells (regardless of the fact whether it is listed in the
      // component_select argument or not)
      //
      // for each element 'f' of the FECollection,
      // local_block_association[f][d] then returns the vector block that
      // degree of freedom 'd' belongs to
      std::vector<std::vector<unsigned char> > local_block_association
      (fe_collection.size());
      for (unsigned int f=0; f<fe_collection.size(); ++f)
        {
          const FiniteElement<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &fe =
            fe_collection[f];
          local_block_association[f].resize(fe.dofs_per_cell,
                                            (unsigned char)(-1));
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            local_block_association[f][i] = fe.system_to_block_index(i).first;

          Assert (std::find (local_block_association[f].begin(),
                             local_block_association[f].end(),
                             (unsigned char)(-1))
                  ==
                  local_block_association[f].end(),
                  ExcInternalError());
        }

      // then loop over all cells and do the work
      std::vector<types::global_dof_index> indices;
      for (typename DoFHandlerType::active_cell_iterator c=dof.begin_active();
           c!=dof.end(); ++ c)
        if (c->is_locally_owned())
          {
            const unsigned int fe_index = c->active_fe_index();
            const unsigned int dofs_per_cell = c->get_fe().dofs_per_cell;
            indices.resize(dofs_per_cell);
            c->get_dof_indices(indices);
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              if (dof.locally_owned_dofs().is_element(indices[i]))
                dofs_by_block[dof.locally_owned_dofs().index_within_set(indices[i])]
                  = local_block_association[fe_index][i];
          }
    }
  }



  template <typename DoFHandlerType, typename Number>
  void distribute_cell_to_dof_vector (const DoFHandlerType &dof_handler,
                                      const Vector<Number> &cell_data,
                                      Vector<double>       &dof_data,
                                      const unsigned int    component)
  {
    const unsigned int dim = DoFHandlerType::dimension;
    const unsigned int spacedim = DoFHandlerType::space_dimension;
    const Triangulation<dim,spacedim> &tria = dof_handler.get_triangulation();
    (void)tria;

    AssertDimension (cell_data.size(), tria.n_active_cells());
    AssertDimension (dof_data.size(), dof_handler.n_dofs());
    AssertIndexRange (component, n_components(dof_handler));
    Assert (fe_is_primitive(dof_handler) == true,
            typename FiniteElement<dim>::ExcFENotPrimitive());

    // store a flag whether we should care about different components. this
    // is just a simplification, we could ask for this at every single
    // place equally well
    const bool consider_components = (n_components(dof_handler) != 1);

    // zero out the components that we will touch
    if (consider_components == false)
      dof_data = 0;
    else
      {
        std::vector<unsigned char> component_dofs (dof_handler.n_locally_owned_dofs());
        internal::get_component_association (dof_handler,
                                             dof_handler.get_fe().component_mask
                                             (FEValuesExtractors::Scalar(component)),
                                             component_dofs);

        for (unsigned int i=0; i<dof_data.size(); ++i)
          if (component_dofs[i] == static_cast<unsigned char>(component))
            dof_data(i) = 0;
      }

    // count how often we have added a value in the sum for each dof
    std::vector<unsigned char> touch_count (dof_handler.n_dofs(), 0);

    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    std::vector<types::global_dof_index> dof_indices;
    dof_indices.reserve (max_dofs_per_cell(dof_handler));

    for (unsigned int present_cell = 0; cell!=endc; ++cell, ++present_cell)
      {
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        dof_indices.resize (dofs_per_cell);
        cell->get_dof_indices (dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          // consider this dof only if it is the right component. if there
          // is only one component, short cut the test
          if (!consider_components ||
              (cell->get_fe().system_to_component_index(i).first == component))
            {
              // sum up contribution of the present_cell to this dof
              dof_data(dof_indices[i]) += cell_data(present_cell);

              // note that we added another summand
              ++touch_count[dof_indices[i]];
            }
      }

    // compute the mean value on all the dofs by dividing with the number
    // of summands.
    for (types::global_dof_index i=0; i<dof_handler.n_dofs(); ++i)
      {
        // assert that each dof was used at least once. this needs not be
        // the case if the vector has more than one component
        Assert (consider_components || (touch_count[i]!=0),
                ExcInternalError());
        if (touch_count[i] != 0)
          dof_data(i) /=  touch_count[i];
      }
  }



  template <int dim, int spacedim>
  void
  extract_dofs (const DoFHandler<dim,spacedim> &dof,
                const ComponentMask            &component_mask,
                std::vector<bool>              &selected_dofs)
  {
    const FiniteElement<dim,spacedim> &fe = dof.get_fe();
    (void)fe;

    Assert(component_mask.represents_n_components(fe.n_components()),
           ExcMessage ("The given component mask is not sized correctly to represent the "
                       "components of the given finite element."));
    Assert(selected_dofs.size() == dof.n_locally_owned_dofs(),
           ExcDimensionMismatch(selected_dofs.size(), dof.n_locally_owned_dofs()));

    // two special cases: no component is selected, and all components are
    // selected; both rather stupid, but easy to catch
    if (component_mask.n_selected_components(n_components(dof)) == 0)
      {
        std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), false);
        return;
      }
    else if (component_mask.n_selected_components(n_components(dof)) == n_components(dof))
      {
        std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), true);
        return;
      }


    // preset all values by false
    std::fill_n (selected_dofs.begin(), dof.n_locally_owned_dofs(), false);

    // get the component association of each DoF and then select the ones
    // that match the given set of blocks
    std::vector<unsigned char> dofs_by_component (dof.n_locally_owned_dofs());
    internal::get_component_association (dof, component_mask,
                                         dofs_by_component);

    for (types::global_dof_index i=0; i<dof.n_locally_owned_dofs(); ++i)
      if (component_mask[dofs_by_component[i]] == true)
        selected_dofs[i] = true;
  }


  // TODO: Unify the following two functions with the non-hp case

  template <int dim, int spacedim>
  void
  extract_dofs (const hp::DoFHandler<dim,spacedim> &dof,
                const ComponentMask                &component_mask,
                std::vector<bool>                  &selected_dofs)
  {
    const FiniteElement<dim,spacedim> &fe = dof.begin_active()->get_fe();
    (void)fe;

    Assert(component_mask.represents_n_components(fe.n_components()),
           ExcMessage ("The given component mask is not sized correctly to represent the "
                       "components of the given finite element."));
    Assert(selected_dofs.size() == dof.n_dofs(),
           ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs()));

    // two special cases: no component is selected, and all components are
    // selected; both rather stupid, but easy to catch
    if (component_mask.n_selected_components(n_components(dof)) == 0)
      {
        std::fill_n (selected_dofs.begin(), dof.n_dofs(), false);
        return;
      }
    else if (component_mask.n_selected_components(n_components(dof)) == n_components(dof))
      {
        std::fill_n (selected_dofs.begin(), dof.n_dofs(), true);
        return;
      }


    // preset all values by false
    std::fill_n (selected_dofs.begin(), dof.n_dofs(), false);

    // get the component association of each DoF and then select the ones
    // that match the given set of components
    std::vector<unsigned char> dofs_by_component (dof.n_dofs());
    internal::get_component_association (dof, component_mask,
                                         dofs_by_component);

    for (types::global_dof_index i=0; i<dof.n_dofs(); ++i)
      if (component_mask[dofs_by_component[i]] == true)
        selected_dofs[i] = true;
  }



  template <int dim, int spacedim>
  void
  extract_dofs (const DoFHandler<dim,spacedim>   &dof,
                const BlockMask     &block_mask,
                std::vector<bool>       &selected_dofs)
  {
    // simply forward to the function that works based on a component mask
    extract_dofs (dof, dof.get_fe().component_mask (block_mask),
                  selected_dofs);
  }



  template <int dim, int spacedim>
  void
  extract_dofs (const hp::DoFHandler<dim,spacedim>   &dof,
                const BlockMask     &block_mask,
                std::vector<bool>       &selected_dofs)
  {
    // simply forward to the function that works based on a component mask
    extract_dofs (dof, dof.get_fe().component_mask (block_mask),
                  selected_dofs);
  }



  template<typename DoFHandlerType>
  void
  extract_level_dofs (const unsigned int    level,
                      const DoFHandlerType &dof,
                      const ComponentMask  &component_mask,
                      std::vector<bool>    &selected_dofs)
  {
    const FiniteElement<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &fe = dof.get_fe();

    Assert(component_mask.represents_n_components(n_components(dof)),
           ExcMessage ("The given component mask is not sized correctly to represent the "
                       "components of the given finite element."));
    Assert(selected_dofs.size() == dof.n_dofs(level),
           ExcDimensionMismatch(selected_dofs.size(), dof.n_dofs(level)));

    // two special cases: no component is selected, and all components are
    // selected, both rather stupid, but easy to catch
    if (component_mask.n_selected_components(n_components(dof)) == 0)
      {
        std::fill_n (selected_dofs.begin(), dof.n_dofs(level), false);
        return;
      }
    else if (component_mask.n_selected_components(n_components(dof)) == n_components(dof))
      {
        std::fill_n (selected_dofs.begin(), dof.n_dofs(level), true);
        return;
      }

    // preset all values by false
    std::fill_n (selected_dofs.begin(), dof.n_dofs(level), false);

    // next set up a table for the degrees of freedom on each of the cells
    // whether it is something interesting or not
    std::vector<unsigned char> local_component_asssociation
      = internal::get_local_component_association (fe, component_mask);
    std::vector<bool> local_selected_dofs (fe.dofs_per_cell);
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      local_selected_dofs[i] = component_mask[local_component_asssociation[i]];

    // then loop over all cells and do work
    std::vector<types::global_dof_index> indices(fe.dofs_per_cell);
    typename DoFHandlerType::level_cell_iterator c;
    for (c = dof.begin(level) ; c != dof.end(level) ; ++ c)
      {
        c->get_mg_dof_indices(indices);
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          selected_dofs[indices[i]] = local_selected_dofs[i];
      }
  }



  template<typename DoFHandlerType>
  void
  extract_level_dofs (const unsigned int    level,
                      const DoFHandlerType &dof,
                      const BlockMask      &block_mask,
                      std::vector<bool>    &selected_dofs)
  {
    // simply defer to the other extract_level_dofs() function
    extract_level_dofs (level, dof, dof.get_fe().component_mask(block_mask),
                        selected_dofs);
  }



  template <typename DoFHandlerType>
  void
  extract_boundary_dofs (const DoFHandlerType               &dof_handler,
                         const ComponentMask                &component_mask,
                         std::vector<bool>                  &selected_dofs,
                         const std::set<types::boundary_id> &boundary_ids)
  {
    Assert ((dynamic_cast<const parallel::distributed::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension>*>
             (&dof_handler.get_triangulation())
             == 0),
            ExcMessage ("This function can not be used with distributed triangulations."
                        "See the documentation for more information."));

    IndexSet indices;
    extract_boundary_dofs (dof_handler, component_mask,
                           indices, boundary_ids);

    // clear and reset array by default values
    selected_dofs.clear ();
    selected_dofs.resize (dof_handler.n_dofs(), false);

    // then convert the values computed above to the binary vector
    indices.fill_binary_vector(selected_dofs);
  }


  template <typename DoFHandlerType>
  void
  extract_boundary_dofs (const DoFHandlerType               &dof_handler,
                         const ComponentMask                &component_mask,
                         IndexSet                           &selected_dofs,
                         const std::set<types::boundary_id> &boundary_ids)
  {
    Assert (component_mask.represents_n_components(n_components(dof_handler)),
            ExcMessage ("Component mask has invalid size."));
    Assert (boundary_ids.find (numbers::internal_face_boundary_id) == boundary_ids.end(),
            ExcInvalidBoundaryIndicator());
    const unsigned int dim=DoFHandlerType::dimension;

    // first reset output argument
    selected_dofs.clear ();
    selected_dofs.set_size(dof_handler.n_dofs());

    // let's see whether we have to check for certain boundary indicators
    // or whether we can accept all
    const bool check_boundary_id = (boundary_ids.size() != 0);

    // also see whether we have to check whether a certain vector component
    // is selected, or all
    const bool check_vector_component
      = ((component_mask.represents_the_all_selected_mask() == false)
         ||
         (component_mask.n_selected_components(n_components(dof_handler)) !=
          n_components(dof_handler)));

    std::vector<types::global_dof_index> face_dof_indices;
    face_dof_indices.reserve (max_dofs_per_face(dof_handler));

    // now loop over all cells and check whether their faces are at the
    // boundary. note that we need not take special care of single lines
    // being at the boundary (using @p{cell->has_boundary_lines}), since we
    // do not support boundaries of dimension dim-2, and so every isolated
    // boundary line is also part of a boundary face which we will be
    // visiting sooner or later
    for (typename DoFHandlerType::active_cell_iterator cell=dof_handler.begin_active();
         cell!=dof_handler.end(); ++cell)

      // only work on cells that are either locally owned or at least ghost
      // cells
      if (cell->is_artificial() == false)
        for (unsigned int face=0;
             face<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++face)
          if (cell->at_boundary(face))
            if (! check_boundary_id ||
                (boundary_ids.find (cell->face(face)->boundary_id())
                 != boundary_ids.end()))
              {
                const FiniteElement<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
                &fe = cell->get_fe();

                const unsigned int dofs_per_face = fe.dofs_per_face;
                face_dof_indices.resize (dofs_per_face);
                cell->face(face)->get_dof_indices (face_dof_indices,
                                                   cell->active_fe_index());

                for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                  if (!check_vector_component)
                    selected_dofs.add_index (face_dof_indices[i]);
                  else
                    // check for component is required. somewhat tricky as
                    // usual for the case that the shape function is
                    // non-primitive, but use usual convention (see docs)
                    {
                      // first get at the cell-global number of a face dof,
                      // to ask the fe certain questions
                      const unsigned int cell_index
                        = (dim == 1 ?
                           i
                           :
                           (dim == 2 ?
                            (i<2*fe.dofs_per_vertex ? i : i+2*fe.dofs_per_vertex)
                            :
                            (dim == 3 ?
                             (i<4*fe.dofs_per_vertex ?
                              i
                              :
                              (i<4*fe.dofs_per_vertex+4*fe.dofs_per_line ?
                               i+4*fe.dofs_per_vertex
                               :
                               i+4*fe.dofs_per_vertex+8*fe.dofs_per_line))
                             :
                             numbers::invalid_unsigned_int)));
                      if (fe.is_primitive (cell_index))
                        {
                          if (component_mask[fe.face_system_to_component_index(i).first]
                              == true)
                            selected_dofs.add_index (face_dof_indices[i]);
                        }
                      else // not primitive
                        {
                          const unsigned int first_nonzero_comp
                            = fe.get_nonzero_components(cell_index).first_selected_component();
                          Assert (first_nonzero_comp < fe.n_components(),
                                  ExcInternalError());

                          if (component_mask[first_nonzero_comp] == true)
                            selected_dofs.add_index (face_dof_indices[i]);
                        }
                    }
              }
  }



  template <typename DoFHandlerType>
  void
  extract_dofs_with_support_on_boundary (const DoFHandlerType               &dof_handler,
                                         const ComponentMask                &component_mask,
                                         std::vector<bool>                  &selected_dofs,
                                         const std::set<types::boundary_id> &boundary_ids)
  {
    Assert (component_mask.represents_n_components (n_components(dof_handler)),
            ExcMessage ("This component mask has the wrong size."));
    Assert (boundary_ids.find (numbers::internal_face_boundary_id) == boundary_ids.end(),
            ExcInvalidBoundaryIndicator());

    // let's see whether we have to check for certain boundary indicators
    // or whether we can accept all
    const bool check_boundary_id = (boundary_ids.size() != 0);

    // also see whether we have to check whether a certain vector component
    // is selected, or all
    const bool check_vector_component
      = (component_mask.represents_the_all_selected_mask() == false);

    // clear and reset array by default values
    selected_dofs.clear ();
    selected_dofs.resize (dof_handler.n_dofs(), false);
    std::vector<types::global_dof_index> cell_dof_indices;
    cell_dof_indices.reserve (max_dofs_per_cell(dof_handler));

    // now loop over all cells and check whether their faces are at the
    // boundary. note that we need not take special care of single lines
    // being at the boundary (using @p{cell->has_boundary_lines}), since we
    // do not support boundaries of dimension dim-2, and so every isolated
    // boundary line is also part of a boundary face which we will be
    // visiting sooner or later
    for (typename DoFHandlerType::active_cell_iterator cell=dof_handler.begin_active();
         cell!=dof_handler.end(); ++cell)
      for (unsigned int face=0;
           face<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++face)
        if (cell->at_boundary(face))
          if (! check_boundary_id ||
              (boundary_ids.find (cell->face(face)->boundary_id())
               != boundary_ids.end()))
            {
              const FiniteElement<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe
                = cell->get_fe();

              const unsigned int dofs_per_cell = fe.dofs_per_cell;
              cell_dof_indices.resize (dofs_per_cell);
              cell->get_dof_indices (cell_dof_indices);

              for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
                if (fe.has_support_on_face(i,face))
                  {
                    if (!check_vector_component)
                      selected_dofs[cell_dof_indices[i]] = true;
                    else
                      // check for component is required. somewhat tricky
                      // as usual for the case that the shape function is
                      // non-primitive, but use usual convention (see docs)
                      {
                        if (fe.is_primitive (i))
                          selected_dofs[cell_dof_indices[i]]
                            = (component_mask[fe.system_to_component_index(i).first]
                               == true);
                        else // not primitive
                          {
                            const unsigned int first_nonzero_comp
                              = fe.get_nonzero_components(i).first_selected_component();
                            Assert (first_nonzero_comp < fe.n_components(),
                                    ExcInternalError());

                            selected_dofs[cell_dof_indices[i]]
                              = (component_mask[first_nonzero_comp]
                                 == true);
                          }
                      }
                  }
            }
  }




  namespace internal
  {
    namespace
    {
      template <int spacedim>
      void extract_hanging_node_dofs (const dealii::DoFHandler<1,spacedim> &dof_handler,
                                      std::vector<bool>           &selected_dofs)
      {
        Assert(selected_dofs.size() == dof_handler.n_dofs(),
               ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
        // preset all values by false
        std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

        // there are no hanging nodes in 1d
      }


      template <int spacedim>
      void extract_hanging_node_dofs (const dealii::DoFHandler<2,spacedim> &dof_handler,
                                      std::vector<bool>           &selected_dofs)
      {
        const unsigned int dim = 2;

        Assert(selected_dofs.size() == dof_handler.n_dofs(),
               ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
        // preset all values by false
        std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

        const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();

        // this function is similar to the make_sparsity_pattern function,
        // see there for more information
        typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->has_children())
              {
                const typename dealii::DoFHandler<dim,spacedim>::line_iterator
                line = cell->face(face);

                for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
                  selected_dofs[line->child(0)->vertex_dof_index(1,dof)] = true;

                for (unsigned int child=0; child<2; ++child)
                  for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
                    selected_dofs[line->child(child)->dof_index(dof)] = true;
              }
      }


      template <int spacedim>
      void extract_hanging_node_dofs (const dealii::DoFHandler<3,spacedim> &dof_handler,
                                      std::vector<bool>           &selected_dofs)
      {
        const unsigned int dim = 3;

        Assert(selected_dofs.size() == dof_handler.n_dofs(),
               ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));
        // preset all values by false
        std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

        const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();

        // this function is similar to the make_sparsity_pattern function,
        // see there for more information

        typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->face(f)->has_children())
              {
                const typename dealii::DoFHandler<dim,spacedim>::face_iterator
                face = cell->face(f);

                for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
                  selected_dofs[face->child(0)->vertex_dof_index(2,dof)] = true;

                // dof numbers on the centers of the lines bounding this
                // face
                for (unsigned int line=0; line<4; ++line)
                  for (unsigned int dof=0; dof!=fe.dofs_per_vertex; ++dof)
                    selected_dofs[face->line(line)->child(0)->vertex_dof_index(1,dof)] = true;

                // next the dofs on the lines interior to the face; the
                // order of these lines is laid down in the FiniteElement
                // class documentation
                for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
                  selected_dofs[face->child(0)->line(1)->dof_index(dof)] = true;
                for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
                  selected_dofs[face->child(1)->line(2)->dof_index(dof)] = true;
                for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
                  selected_dofs[face->child(2)->line(3)->dof_index(dof)] = true;
                for (unsigned int dof=0; dof<fe.dofs_per_line; ++dof)
                  selected_dofs[face->child(3)->line(0)->dof_index(dof)] = true;

                // dofs on the bordering lines
                for (unsigned int line=0; line<4; ++line)
                  for (unsigned int child=0; child<2; ++child)
                    for (unsigned int dof=0; dof!=fe.dofs_per_line; ++dof)
                      selected_dofs[face->line(line)->child(child)->dof_index(dof)] = true;

                // finally, for the dofs interior to the four child faces
                for (unsigned int child=0; child<4; ++child)
                  for (unsigned int dof=0; dof!=fe.dofs_per_quad; ++dof)
                    selected_dofs[face->child(child)->dof_index(dof)] = true;
              }
      }
    }
  }



  template <int dim, int spacedim>
  void

  extract_hanging_node_dofs (const DoFHandler<dim,spacedim> &dof_handler,
                             std::vector<bool>              &selected_dofs)
  {
    internal::extract_hanging_node_dofs (dof_handler,
                                         selected_dofs);
  }



  template <typename DoFHandlerType>
  void
  extract_subdomain_dofs (const DoFHandlerType      &dof_handler,
                          const types::subdomain_id  subdomain_id,
                          std::vector<bool>         &selected_dofs)
  {
    Assert(selected_dofs.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(selected_dofs.size(), dof_handler.n_dofs()));

    // preset all values by false
    std::fill_n (selected_dofs.begin(), dof_handler.n_dofs(), false);

    std::vector<types::global_dof_index> local_dof_indices;
    local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

    // this function is similar to the make_sparsity_pattern function, see
    // there for more information
    typename DoFHandlerType::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->subdomain_id() == subdomain_id)
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          local_dof_indices.resize (dofs_per_cell);
          cell->get_dof_indices (local_dof_indices);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            selected_dofs[local_dof_indices[i]] = true;
        };
  }



  template <typename DoFHandlerType>
  void
  extract_locally_owned_dofs (const DoFHandlerType &dof_handler,
                              IndexSet             &dof_set)
  {
    // collect all the locally owned dofs
    dof_set = dof_handler.locally_owned_dofs();
    dof_set.compress ();
  }



  template <typename DoFHandlerType>
  void
  extract_locally_active_dofs (const DoFHandlerType &dof_handler,
                               IndexSet             &dof_set)
  {
    // collect all the locally owned dofs
    dof_set = dof_handler.locally_owned_dofs();

    // add the DoF on the adjacent ghost cells to the IndexSet, cache them
    // in a set. need to check each dof manually because we can't be sure
    // that the dof range of locally_owned_dofs is really contiguous.
    std::vector<types::global_dof_index> dof_indices;
    std::set<types::global_dof_index> global_dof_indices;

    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          dof_indices.resize(cell->get_fe().dofs_per_cell);
          cell->get_dof_indices(dof_indices);

          for (std::vector<types::global_dof_index>::iterator it=dof_indices.begin();
               it!=dof_indices.end();
               ++it)
            if (!dof_set.is_element(*it))
              global_dof_indices.insert(*it);
        }

    dof_set.add_indices(global_dof_indices.begin(), global_dof_indices.end());

    dof_set.compress();
  }



  template <typename DoFHandlerType>
  void
  extract_locally_relevant_dofs (const DoFHandlerType &dof_handler,
                                 IndexSet             &dof_set)
  {
    // collect all the locally owned dofs
    dof_set = dof_handler.locally_owned_dofs();

    // now add the DoF on the adjacent ghost cells to the IndexSet

    // Note: For certain meshes (in particular in 3D and with many
    // processors), it is really necessary to cache intermediate data. After
    // trying several objects such as std::set, a vector that is always kept
    // sorted, and a vector that is initially unsorted and sorted once at the
    // end, the latter has been identified to provide the best performance.
    // Martin Kronbichler
    std::vector<types::global_dof_index> dof_indices;
    std::vector<types::global_dof_index> dofs_on_ghosts;

    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_ghost())
        {
          dof_indices.resize(cell->get_fe().dofs_per_cell);
          cell->get_dof_indices(dof_indices);
          for (unsigned int i=0; i<dof_indices.size(); ++i)
            if (!dof_set.is_element(dof_indices[i]))
              dofs_on_ghosts.push_back(dof_indices[i]);
        }

    // sort, compress out duplicates, fill into index set
    std::sort(dofs_on_ghosts.begin(), dofs_on_ghosts.end());
    dof_set.add_indices(dofs_on_ghosts.begin(), std::unique(dofs_on_ghosts.begin(),
                                                            dofs_on_ghosts.end()));
    dof_set.compress();
  }


  template <typename DoFHandlerType>
  void
  extract_locally_relevant_level_dofs (const DoFHandlerType &dof_handler,
                                       const unsigned int    level,
                                       IndexSet             &dof_set)
  {
    // collect all the locally owned dofs
    dof_set = dof_handler.locally_owned_mg_dofs(level);

    // add the DoF on the adjacent ghost cells to the IndexSet

    // Note: For certain meshes (in particular in 3D and with many
    // processors), it is really necessary to cache intermediate data. After
    // trying several objects such as std::set, a vector that is always kept
    // sorted, and a vector that is initially unsorted and sorted once at the
    // end, the latter has been identified to provide the best performance.
    // Martin Kronbichler
    std::vector<types::global_dof_index> dof_indices;
    std::vector<types::global_dof_index> dofs_on_ghosts;

    typename DoFHandlerType::cell_iterator cell = dof_handler.begin(level),
                                           endc = dof_handler.end(level);
    for (; cell!=endc; ++cell)
      {
        const types::subdomain_id id = cell->level_subdomain_id();

        // skip artificial and own cells (only look at ghost cells)
        if (id == dof_handler.get_triangulation().locally_owned_subdomain()
            || id == numbers::artificial_subdomain_id)
          continue;

        dof_indices.resize(cell->get_fe().dofs_per_cell);
        cell->get_mg_dof_indices(dof_indices);
        for (unsigned int i=0; i<dof_indices.size(); ++i)
          if (!dof_set.is_element(dof_indices[i]))
            dofs_on_ghosts.push_back(dof_indices[i]);
      }

    // sort, compress out duplicates, fill into index set
    std::sort(dofs_on_ghosts.begin(), dofs_on_ghosts.end());
    dof_set.add_indices(dofs_on_ghosts.begin(), std::unique(dofs_on_ghosts.begin(),
                                                            dofs_on_ghosts.end()));

    dof_set.compress();
  }

  template <typename DoFHandlerType>
  void
  extract_constant_modes (const DoFHandlerType            &dof_handler,
                          const ComponentMask             &component_mask,
                          std::vector<std::vector<bool> > &constant_modes)
  {
    const unsigned int n_components = dof_handler.get_fe().n_components();
    Assert (component_mask.represents_n_components(n_components),
            ExcDimensionMismatch(n_components,
                                 component_mask.size()));

    std::vector<unsigned char> dofs_by_component (dof_handler.n_locally_owned_dofs());
    internal::get_component_association (dof_handler, component_mask,
                                         dofs_by_component);
    unsigned int n_selected_dofs = 0;
    for (unsigned int i=0; i<n_components; ++i)
      if (component_mask[i] == true)
        n_selected_dofs += std::count (dofs_by_component.begin(),
                                       dofs_by_component.end(), i);

    // Find local numbering within the selected components
    const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
    std::vector<unsigned int> component_numbering(locally_owned_dofs.n_elements(),
                                                  numbers::invalid_unsigned_int);
    for (unsigned int i=0, count=0; i<locally_owned_dofs.n_elements(); ++i)
      if (component_mask[dofs_by_component[i]])
        component_numbering[i] = count++;

    // get the element constant modes and find a translation table between
    // index in the constant modes and the components.
    //
    // TODO: We might be able to extend this also for elements which do not
    // have the same constant modes, but that is messy...
    const dealii::hp::FECollection<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
    fe_collection (dof_handler.get_fe());
    std::vector<Table<2,bool> > element_constant_modes;
    std::vector<std::vector<std::pair<unsigned int, unsigned int> > >
    constant_mode_to_component_translation(n_components);
    unsigned int n_constant_modes = 0;
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
        std::pair<Table<2,bool>, std::vector<unsigned int> > data
          = fe_collection[f].get_constant_modes();
        element_constant_modes.push_back(data.first);
        if (f==0)
          for (unsigned int i=0; i<data.second.size(); ++i)
            if (component_mask[data.second[i]])
              constant_mode_to_component_translation[data.second[i]].
              push_back(std::make_pair(n_constant_modes++,i));
        AssertDimension(element_constant_modes.back().n_rows(),
                        element_constant_modes[0].n_rows());
      }

    // First count the number of dofs in the current component.
    constant_modes.clear ();
    constant_modes.resize (n_constant_modes, std::vector<bool>(n_selected_dofs,
                                                               false));

    // Loop over all owned cells and ask the element for the constant modes

    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    std::vector<types::global_dof_index> dof_indices;
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          dof_indices.resize(cell->get_fe().dofs_per_cell);
          cell->get_dof_indices(dof_indices);

          for (unsigned int i=0; i<dof_indices.size(); ++i)
            if (locally_owned_dofs.is_element(dof_indices[i]))
              {
                const unsigned int loc_index =
                  locally_owned_dofs.index_within_set(dof_indices[i]);
                const unsigned int comp = dofs_by_component[loc_index];
                if (component_mask[comp])
                  for (unsigned int j=0; j<constant_mode_to_component_translation[comp].size(); ++j)
                    constant_modes[constant_mode_to_component_translation[comp][j].first]
                    [component_numbering[loc_index]] =
                      element_constant_modes[cell->active_fe_index()]
                      (constant_mode_to_component_translation[comp][j].second,i);
              }
        }
  }



  template <typename DoFHandlerType>
  void
  get_active_fe_indices (const DoFHandlerType      &dof_handler,
                         std::vector<unsigned int> &active_fe_indices)
  {
    AssertDimension (active_fe_indices.size(), dof_handler.get_triangulation().n_active_cells());

    typename DoFHandlerType::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      active_fe_indices[cell->active_cell_index()] = cell->active_fe_index();
  }

  template <typename DoFHandlerType>
  std::vector<IndexSet>
  locally_owned_dofs_per_subdomain (const DoFHandlerType  &dof_handler)
  {
    // If the Triangulation is distributed, the only thing we can usefully
    // ask is for its locally owned subdomain
    Assert ((dynamic_cast<const parallel::distributed::
             Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *>
             (&dof_handler.get_triangulation()) == 0),
            ExcMessage ("For parallel::distributed::Triangulation objects and "
                        "associated DoF handler objects, asking for any information "
                        "related to a subdomain other than the locally owned one does "
                        "not make sense."));

    //the following is a random process (flip of a coin), thus should be called once only.
    std::vector< dealii::types::subdomain_id > subdomain_association (dof_handler.n_dofs ());
    dealii::DoFTools::get_subdomain_association (dof_handler, subdomain_association);

    const unsigned int n_subdomains
      = (dynamic_cast<const parallel::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *>
         (&dof_handler.get_triangulation()) == 0
         ?
         1
         :
         Utilities::MPI::n_mpi_processes
         (dynamic_cast<const parallel::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *>
          (&dof_handler.get_triangulation())->get_communicator()));
    Assert (n_subdomains >
            *std::max_element (subdomain_association.begin (),
                               subdomain_association.end ()),
            ExcInternalError());

    std::vector<dealii::IndexSet> index_sets (n_subdomains,
                                              dealii::IndexSet(dof_handler.n_dofs()));

    // loop over subdomain_association and populate IndexSet when a
    // change in subdomain ID is found
    dealii::types::global_dof_index i_min          = 0;
    dealii::types::global_dof_index this_subdomain = subdomain_association[0];

    for (dealii::types::global_dof_index index = 1;
         index < subdomain_association.size (); ++index)
      {
        //found index different from the current one
        if (subdomain_association[index] != this_subdomain)
          {
            index_sets[this_subdomain].add_range (i_min, index);
            i_min = index;
            this_subdomain = subdomain_association[index];
          }
      }

    // the very last element is of different index
    if (i_min == subdomain_association.size () - 1)
      {
        index_sets[this_subdomain].add_index (i_min);
      }

    // otherwise there are at least two different indices
    else
      {
        index_sets[this_subdomain].add_range (
          i_min, subdomain_association.size ());
      }

    for (unsigned int i = 0; i < n_subdomains; i++)
      index_sets[i].compress ();

    return index_sets;
  }

  template <typename DoFHandlerType>
  std::vector<IndexSet>
  locally_relevant_dofs_per_subdomain (const DoFHandlerType  &dof_handler)
  {
    // If the Triangulation is distributed, the only thing we can usefully
    // ask is for its locally owned subdomain
    Assert ((dynamic_cast<const parallel::distributed::
             Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *>
             (&dof_handler.get_triangulation()) == 0),
            ExcMessage ("For parallel::distributed::Triangulation objects and "
                        "associated DoF handler objects, asking for any information "
                        "related to a subdomain other than the locally owned one does "
                        "not make sense."));

    // Collect all the locally owned DoFs
    // Note: Even though the distribution of DoFs by the locally_owned_dofs_per_subdomain
    // function is pseudo-random, we will collect all the DoFs on the subdomain
    // and its layer cell. Therefore, the random nature of this function does
    // not play a role in the extraction of the locally relevant DoFs
    std::vector<IndexSet> dof_set = locally_owned_dofs_per_subdomain(dof_handler);
    const dealii::types::subdomain_id n_subdomains = dof_set.size();

    // Add the DoFs on the adjacent (equivalent ghost) cells to the IndexSet,
    // cache them in a set. Need to check each DoF manually because we can't
    // be sure that the DoF range of locally_owned_dofs is really contiguous.
    for (dealii::types::subdomain_id subdomain_id = 0;
         subdomain_id < n_subdomains; ++subdomain_id)
      {
        // Extract the layer of cells around this subdomain
        std_cxx11::function<bool (const typename DoFHandlerType::active_cell_iterator &)> predicate
          = IteratorFilters::SubdomainEqualTo(subdomain_id);
        const std::vector<typename DoFHandlerType::active_cell_iterator>
        active_halo_layer = GridTools::compute_active_cell_halo_layer (dof_handler,
                            predicate);

        // Extract DoFs associated with halo layer
        std::vector<types::global_dof_index> local_dof_indices;
        std::set<types::global_dof_index> subdomain_halo_global_dof_indices;
        for (typename std::vector<typename DoFHandlerType::active_cell_iterator>::const_iterator
             it_cell = active_halo_layer.begin(); it_cell!=active_halo_layer.end(); ++it_cell)
          {
            const typename DoFHandlerType::active_cell_iterator &cell = *it_cell;
            Assert(cell->subdomain_id() != subdomain_id,
                   ExcMessage("The subdomain ID of the halo cell should not match that of the vector entry."));

            local_dof_indices.resize(cell->get_fe().dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);

            for (std::vector<types::global_dof_index>::iterator it=local_dof_indices.begin();
                 it!=local_dof_indices.end();
                 ++it)
              subdomain_halo_global_dof_indices.insert(*it);
          }

        dof_set[subdomain_id].add_indices(subdomain_halo_global_dof_indices.begin(),
                                          subdomain_halo_global_dof_indices.end());

        dof_set[subdomain_id].compress();
      }

    return dof_set;
  }

  template <typename DoFHandlerType>
  void
  get_subdomain_association (const DoFHandlerType &dof_handler,
                             std::vector<types::subdomain_id> &subdomain_association)
  {
    // if the Triangulation is distributed, the only thing we can usefully
    // ask is for its locally owned subdomain
    Assert ((dynamic_cast<const parallel::distributed::
             Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *>
             (&dof_handler.get_triangulation()) == 0),
            ExcMessage ("For parallel::distributed::Triangulation objects and "
                        "associated DoF handler objects, asking for any subdomain other "
                        "than the locally owned one does not make sense."));

    Assert(subdomain_association.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(subdomain_association.size(),
                                dof_handler.n_dofs()));

    Assert(dof_handler.n_dofs() > 0,
           ExcMessage("Number of DoF is not positive. "
                      "This could happen when the function is called before NumberCache is written."));

    // In case this function is executed with parallel::shared::Triangulation
    // with possibly artifical cells, we need to take "true" subdomain IDs (i.e. without
    // artificial cells). Otherwise we are good to use subdomain_id as stored
    // in cell->subdomain_id().
    std::vector<types::subdomain_id> cell_owners (dof_handler.get_triangulation().n_active_cells());
    if (const parallel::shared::Triangulation<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
        *tr = (dynamic_cast<const parallel::shared::Triangulation<DoFHandlerType::dimension,
               DoFHandlerType::space_dimension>*> (&dof_handler.get_triangulation())))
      {
        cell_owners = tr->get_true_subdomain_ids_of_cells();
        Assert (tr->get_true_subdomain_ids_of_cells().size() == tr->n_active_cells(),
                ExcInternalError());
      }
    else
      {
        for (typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active();
             cell!= dof_handler.end(); cell++)
          if (cell->is_locally_owned())
            cell_owners[cell->active_cell_index()] = cell->subdomain_id();
      }

    // preset all values by an invalid value
    std::fill_n (subdomain_association.begin(), dof_handler.n_dofs(),
                 numbers::invalid_subdomain_id);

    std::vector<types::global_dof_index> local_dof_indices;
    local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

    // pseudo-randomly assign variables which lie on the interface between
    // subdomains to each of the two or more
    bool coin_flip = true;

    // loop over all cells and record which subdomain a DoF belongs to.
    // toss a coin in case it is on an interface
    typename DoFHandlerType::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        const types::subdomain_id subdomain_id = cell_owners[cell->active_cell_index()];
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        local_dof_indices.resize (dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);

        // set subdomain ids. if dofs already have their values set then
        // they must be on partition interfaces. in that case randomly
        // assign them to either the previous association or the current
        // one, where we take "random" to be "once this way once that way"
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (subdomain_association[local_dof_indices[i]] ==
              numbers::invalid_unsigned_int)
            subdomain_association[local_dof_indices[i]] = subdomain_id;
          else
            {
              if (coin_flip == true)
                subdomain_association[local_dof_indices[i]] = subdomain_id;
              coin_flip = !coin_flip;
            }
      }

    Assert (std::find (subdomain_association.begin(),
                       subdomain_association.end(),
                       numbers::invalid_subdomain_id)
            == subdomain_association.end(),
            ExcInternalError());
  }



  template <typename DoFHandlerType>
  unsigned int
  count_dofs_with_subdomain_association (const DoFHandlerType      &dof_handler,
                                         const types::subdomain_id  subdomain)
  {
    std::vector<types::subdomain_id> subdomain_association (dof_handler.n_dofs());
    get_subdomain_association (dof_handler, subdomain_association);

    return std::count (subdomain_association.begin(),
                       subdomain_association.end(),
                       subdomain);
  }



  template <typename DoFHandlerType>
  IndexSet
  dof_indices_with_subdomain_association (const DoFHandlerType      &dof_handler,
                                          const types::subdomain_id  subdomain)
  {

    // If we have a distributed::Triangulation only allow locally_owned
    // subdomain.
    Assert (
      (dof_handler.get_triangulation().locally_owned_subdomain() == numbers::invalid_subdomain_id)
      ||
      (subdomain == dof_handler.get_triangulation().locally_owned_subdomain()),
      ExcMessage ("For parallel::distributed::Triangulation objects and "
                  "associated DoF handler objects, asking for any subdomain other "
                  "than the locally owned one does not make sense."));

    IndexSet index_set (dof_handler.n_dofs());

    std::vector<types::global_dof_index> local_dof_indices;
    local_dof_indices.reserve (max_dofs_per_cell(dof_handler));

    // first generate an unsorted list of all indices which we fill from
    // the back. could also insert them directly into the IndexSet, but
    // that inserts indices in the middle, which is an O(n^2) algorithm and
    // hence too expensive. Could also use std::set, but that is in general
    // more expensive than a vector
    std::vector<types::global_dof_index> subdomain_indices;

    typename DoFHandlerType::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if ((cell->is_artificial() == false)
          &&
          (cell->subdomain_id() == subdomain))
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          local_dof_indices.resize (dofs_per_cell);
          cell->get_dof_indices (local_dof_indices);
          subdomain_indices.insert(subdomain_indices.end(),
                                   local_dof_indices.begin(),
                                   local_dof_indices.end());
        }
    // sort indices and remove duplicates
    std::sort (subdomain_indices.begin(), subdomain_indices.end());
    subdomain_indices.erase (std::unique(subdomain_indices.begin(),
                                         subdomain_indices.end()),
                             subdomain_indices.end());

    // insert into IndexSet
    index_set.add_indices (subdomain_indices.begin(), subdomain_indices.end());
    index_set.compress ();

    return index_set;
  }



  template <typename DoFHandlerType>
  void
  count_dofs_with_subdomain_association (const DoFHandlerType      &dof_handler,
                                         const types::subdomain_id  subdomain,
                                         std::vector<unsigned int> &n_dofs_on_subdomain)
  {
    Assert (n_dofs_on_subdomain.size() == dof_handler.get_fe().n_components(),
            ExcDimensionMismatch (n_dofs_on_subdomain.size(),
                                  dof_handler.get_fe().n_components()));
    std::fill (n_dofs_on_subdomain.begin(), n_dofs_on_subdomain.end(), 0);

    // in debug mode, make sure that there are some cells at least with
    // this subdomain id
#ifdef DEBUG
    {
      bool found = false;
      for (typename Triangulation<DoFHandlerType::dimension,
           DoFHandlerType::space_dimension>::active_cell_iterator
           cell=dof_handler.get_triangulation().begin_active();
           cell!=dof_handler.get_triangulation().end(); ++cell)
        if (cell->subdomain_id() == subdomain)
          {
            found = true;
            break;
          }
      Assert (found == true,
              ExcMessage ("There are no cells for the given subdomain!"));
    }
#endif

    std::vector<types::subdomain_id> subdomain_association (dof_handler.n_dofs());
    get_subdomain_association (dof_handler, subdomain_association);

    std::vector<unsigned char> component_association (dof_handler.n_dofs());
    internal::get_component_association (dof_handler, std::vector<bool>(),
                                         component_association);

    for (unsigned int c=0; c<dof_handler.get_fe().n_components(); ++c)
      {
        for (types::global_dof_index i=0; i<dof_handler.n_dofs(); ++i)
          if ((subdomain_association[i] == subdomain) &&
              (component_association[i] == static_cast<unsigned char>(c)))
            ++n_dofs_on_subdomain[c];
      }
  }



  namespace internal
  {
    // TODO: why is this function so complicated? It would be nice to have
    // comments that explain why we can't just loop over all components and
    // count the entries in dofs_by_component that have this component's
    // index
    template <int dim, int spacedim>
    void
    resolve_components (const FiniteElement<dim,spacedim> &fe,
                        const std::vector<unsigned char> &dofs_by_component,
                        const std::vector<unsigned int>  &target_component,
                        const bool                        only_once,
                        std::vector<types::global_dof_index> &dofs_per_component,
                        unsigned int                     &component)
    {
      for (unsigned int b=0; b<fe.n_base_elements(); ++b)
        {
          const FiniteElement<dim,spacedim> &base = fe.base_element(b);
          // Dimension of base element
          unsigned int d = base.n_components();

          for (unsigned int m=0; m<fe.element_multiplicity(b); ++m)
            {
              if (base.n_base_elements() > 1)
                resolve_components(base, dofs_by_component, target_component,
                                   only_once, dofs_per_component, component);
              else
                {
                  for (unsigned int dd=0; dd<d; ++dd,++component)
                    dofs_per_component[target_component[component]]
                    += std::count(dofs_by_component.begin(),
                                  dofs_by_component.end(),
                                  component);

                  // if we have non-primitive FEs and want all components
                  // to show the number of dofs, need to copy the result to
                  // those components
                  if (!base.is_primitive() && !only_once)
                    for (unsigned int dd=1; dd<d; ++dd)
                      dofs_per_component[target_component[component-d+dd]] =
                        dofs_per_component[target_component[component-d]];
                }
            }
        }
    }


    template <int dim, int spacedim>
    void
    resolve_components (const hp::FECollection<dim,spacedim> &fe_collection,
                        const std::vector<unsigned char> &dofs_by_component,
                        const std::vector<unsigned int>  &target_component,
                        const bool                        only_once,
                        std::vector<types::global_dof_index> &dofs_per_component,
                        unsigned int                     &component)
    {
      // assert that all elements in the collection have the same structure
      // (base elements and multiplicity, components per base element) and
      // then simply call the function above
      for (unsigned int fe=1; fe<fe_collection.size(); ++fe)
        {
          Assert (fe_collection[fe].n_components() == fe_collection[0].n_components(),
                  ExcNotImplemented());
          Assert (fe_collection[fe].n_base_elements() == fe_collection[0].n_base_elements(),
                  ExcNotImplemented());
          for (unsigned int b=0; b<fe_collection[0].n_base_elements(); ++b)
            {
              Assert (fe_collection[fe].base_element(b).n_components() == fe_collection[0].base_element(b).n_components(),
                      ExcNotImplemented());
              Assert (fe_collection[fe].base_element(b).n_base_elements() == fe_collection[0].base_element(b).n_base_elements(),
                      ExcNotImplemented());
            }
        }

      resolve_components (fe_collection[0], dofs_by_component,
                          target_component, only_once, dofs_per_component,
                          component);
    }
  }



  namespace internal
  {
    namespace
    {
      /**
       * Return true if the given element is primitive.
       */
      template <int dim, int spacedim>
      bool all_elements_are_primitive (const FiniteElement<dim,spacedim> &fe)
      {
        return fe.is_primitive();
      }


      /**
       * Return true if each element of the given element collection is primitive.
       */
      template <int dim, int spacedim>
      bool all_elements_are_primitive (const dealii::hp::FECollection<dim,spacedim> &fe_collection)
      {
        for (unsigned int i=0; i<fe_collection.size(); ++i)
          if (fe_collection[i].is_primitive() == false)
            return false;

        return true;
      }
    }
  }

  template <typename DoFHandlerType>
  void
  count_dofs_per_component (const DoFHandlerType                 &dof_handler,
                            std::vector<types::global_dof_index> &dofs_per_component,
                            bool                                  only_once,
                            std::vector<unsigned int>             target_component)
  {
    const unsigned int n_components = dof_handler.get_fe().n_components();

    std::fill (dofs_per_component.begin(), dofs_per_component.end(),
               types::global_dof_index(0));

    // If the empty vector was given as default argument, set up this
    // vector as identity.
    if (target_component.size()==0)
      {
        target_component.resize(n_components);
        for (unsigned int i=0; i<n_components; ++i)
          target_component[i] = i;
      }
    else
      Assert (target_component.size()==n_components,
              ExcDimensionMismatch(target_component.size(),
                                   n_components));


    const unsigned int max_component
      = *std::max_element (target_component.begin(),
                           target_component.end());
    const unsigned int n_target_components = max_component + 1;
    (void)n_target_components; // silence possible warning about unused variable

    AssertDimension (dofs_per_component.size(), n_target_components);

    // special case for only one component. treat this first since it does
    // not require any computations
    if (n_components == 1)
      {
        dofs_per_component[0] = dof_handler.n_locally_owned_dofs();
        return;
      }


    // otherwise determine the number of dofs in each component separately.
    // do so in parallel
    std::vector<unsigned char> dofs_by_component (dof_handler.n_locally_owned_dofs());
    internal::get_component_association (dof_handler, ComponentMask(),
                                         dofs_by_component);

    // next count what we got
    unsigned int component = 0;
    internal::resolve_components(dof_handler.get_fe(),
                                 dofs_by_component, target_component,
                                 only_once, dofs_per_component, component);
    Assert (n_components == component, ExcInternalError());

    // finally sanity check. this is only valid if the finite element is
    // actually primitive, so exclude other elements from this
    Assert ((internal::all_elements_are_primitive(dof_handler.get_fe()) == false)
            ||
            (std::accumulate (dofs_per_component.begin(),
                              dofs_per_component.end(),
                              types::global_dof_index(0))
             == dof_handler.n_locally_owned_dofs()),
            ExcInternalError());

    // reduce information from all CPUs
#ifdef DEAL_II_WITH_MPI
    const unsigned int dim = DoFHandlerType::dimension;
    const unsigned int spacedim = DoFHandlerType::space_dimension;

    if (const parallel::Triangulation<dim,spacedim> *tria
        = (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
           (&dof_handler.get_triangulation())))
      {
        std::vector<types::global_dof_index> local_dof_count = dofs_per_component;

        MPI_Allreduce ( &local_dof_count[0], &dofs_per_component[0], n_target_components,
                        DEAL_II_DOF_INDEX_MPI_TYPE,
                        MPI_SUM, tria->get_communicator());
      }
#endif
  }



  template <typename DoFHandlerType>
  void
  count_dofs_per_block (const DoFHandlerType                 &dof_handler,
                        std::vector<types::global_dof_index> &dofs_per_block,
                        const std::vector<unsigned int>      &target_block_)
  {
    std::vector<unsigned int>  target_block = target_block_;

    const dealii::hp::FECollection<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
    fe_collection (dof_handler.get_fe());
    Assert (fe_collection.size() < 256, ExcNotImplemented());

    for (unsigned int this_fe=0; this_fe<fe_collection.size(); ++this_fe)
      {
        const FiniteElement<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &fe = fe_collection[this_fe];
        std::fill (dofs_per_block.begin(), dofs_per_block.end(),
                   types::global_dof_index(0));

        // If the empty vector was given as default argument, set up this
        // vector as identity.
        if (target_block.size()==0)
          {
            target_block.resize(fe.n_blocks());
            for (unsigned int i=0; i<fe.n_blocks(); ++i)
              target_block[i] = i;
          }
        else
          Assert (target_block.size()==fe.n_blocks(),
                  ExcDimensionMismatch(target_block.size(),
                                       fe.n_blocks()));



        const unsigned int max_block
          = *std::max_element (target_block.begin(),
                               target_block.end());
        const unsigned int n_target_blocks = max_block + 1;
        (void)n_target_blocks; // silence possible warning about unused variable

        const unsigned int n_blocks = fe.n_blocks();

        AssertDimension (dofs_per_block.size(), n_target_blocks);

        // special case for only one block. treat this first since it does
        // not require any computations
        if (n_blocks == 1)
          {
            dofs_per_block[0] = dof_handler.n_dofs();
            return;
          }
        // otherwise determine the number of dofs in each block separately.
        std::vector<unsigned char> dofs_by_block (dof_handler.n_locally_owned_dofs());
        internal::get_block_association (dof_handler, dofs_by_block);

        // next count what we got
        for (unsigned int block=0; block<fe.n_blocks(); ++block)
          dofs_per_block[target_block[block]]
          += std::count(dofs_by_block.begin(), dofs_by_block.end(),
                        block);

#ifdef DEAL_II_WITH_MPI
        // if we are working on a parallel mesh, we now need to collect
        // this information from all processors
        if (const parallel::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> *tria
            = (dynamic_cast<const parallel::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension>*>
               (&dof_handler.get_triangulation())))
          {
            std::vector<types::global_dof_index> local_dof_count = dofs_per_block;
            MPI_Allreduce ( &local_dof_count[0], &dofs_per_block[0],
                            n_target_blocks,
                            DEAL_II_DOF_INDEX_MPI_TYPE,
                            MPI_SUM, tria->get_communicator());
          }
#endif
      }
  }



  template <typename DoFHandlerType>
  void
  map_dof_to_boundary_indices (const DoFHandlerType &dof_handler,
                               std::vector<types::global_dof_index> &mapping)
  {
    mapping.clear ();
    mapping.insert (mapping.end(), dof_handler.n_dofs(),
                    DoFHandlerType::invalid_dof_index);

    std::vector<types::global_dof_index> dofs_on_face;
    dofs_on_face.reserve (max_dofs_per_face(dof_handler));
    types::global_dof_index next_boundary_index = 0;

    // now loop over all cells and check whether their faces are at the
    // boundary. note that we need not take special care of single lines
    // being at the boundary (using @p{cell->has_boundary_lines}), since we
    // do not support boundaries of dimension dim-2, and so every isolated
    // boundary line is also part of a boundary face which we will be
    // visiting sooner or later
    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);
            cell->face(f)->get_dof_indices (dofs_on_face,
                                            cell->active_fe_index());
            for (unsigned int i=0; i<dofs_per_face; ++i)
              if (mapping[dofs_on_face[i]] == DoFHandlerType::invalid_dof_index)
                mapping[dofs_on_face[i]] = next_boundary_index++;
          }

    AssertDimension (next_boundary_index, dof_handler.n_boundary_dofs());
  }



  template <typename DoFHandlerType>
  void map_dof_to_boundary_indices
  (const DoFHandlerType                 &dof_handler,
   const std::set<types::boundary_id>   &boundary_ids,
   std::vector<types::global_dof_index> &mapping)
  {
    Assert (boundary_ids.find (numbers::internal_face_boundary_id) == boundary_ids.end(),
            ExcInvalidBoundaryIndicator());

    mapping.clear ();
    mapping.insert (mapping.end(), dof_handler.n_dofs(),
                    DoFHandlerType::invalid_dof_index);

    // return if there is nothing to do
    if (boundary_ids.size() == 0)
      return;

    std::vector<types::global_dof_index> dofs_on_face;
    dofs_on_face.reserve (max_dofs_per_face(dof_handler));
    types::global_dof_index next_boundary_index = 0;

    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
        if (boundary_ids.find (cell->face(f)->boundary_id()) !=
            boundary_ids.end())
          {
            const unsigned int dofs_per_face = cell->get_fe().dofs_per_face;
            dofs_on_face.resize (dofs_per_face);
            cell->face(f)->get_dof_indices (dofs_on_face, cell->active_fe_index());
            for (unsigned int i=0; i<dofs_per_face; ++i)
              if (mapping[dofs_on_face[i]] == DoFHandlerType::invalid_dof_index)
                mapping[dofs_on_face[i]] = next_boundary_index++;
          }

    AssertDimension (next_boundary_index,
                     dof_handler.n_boundary_dofs (boundary_ids));
  }

  namespace internal
  {
    namespace
    {
      template <typename DoFHandlerType>
      void
      map_dofs_to_support_points
      (const hp::MappingCollection<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &mapping,
       const DoFHandlerType                                                      &dof_handler,
       std::map<types::global_dof_index,Point<DoFHandlerType::space_dimension> > &support_points)
      {
        const unsigned int dim = DoFHandlerType::dimension;
        const unsigned int spacedim = DoFHandlerType::space_dimension;

        hp::FECollection<dim, spacedim> fe_collection(dof_handler.get_fe());
        hp::QCollection<dim> q_coll_dummy;

        for (unsigned int fe_index = 0; fe_index < fe_collection.size(); ++fe_index)
          {
            // check whether every fe in the collection has support points
            Assert(fe_collection[fe_index].has_support_points(),
                   typename FiniteElement<dim>::ExcFEHasNoSupportPoints());
            q_coll_dummy.push_back(
              Quadrature<dim> (
                fe_collection[fe_index].get_unit_support_points()));
          }

        // Now loop over all cells and enquire the support points on each
        // of these. we use dummy quadrature formulas where the quadrature
        // points are located at the unit support points to enquire the
        // location of the support points in real space.
        //
        // The weights of the quadrature rule have been set to invalid
        // values by the used constructor.
        hp::FEValues<dim, spacedim> hp_fe_values(mapping, fe_collection,
                                                 q_coll_dummy, update_quadrature_points);
        typename DoFHandlerType::active_cell_iterator cell =
          dof_handler.begin_active(), endc = dof_handler.end();

        std::vector<types::global_dof_index> local_dof_indices;
        for (; cell != endc; ++cell)
          // only work on locally relevant cells
          if (cell->is_artificial() == false)
            {
              hp_fe_values.reinit(cell);
              const FEValues<dim, spacedim> &fe_values = hp_fe_values.get_present_fe_values();

              local_dof_indices.resize(cell->get_fe().dofs_per_cell);
              cell->get_dof_indices(local_dof_indices);

              const std::vector<Point<spacedim> > &points =
                fe_values.get_quadrature_points();
              for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell; ++i)
                // insert the values into the map
                support_points[local_dof_indices[i]] = points[i];
            }
      }


      template <typename DoFHandlerType>
      void
      map_dofs_to_support_points
      (const hp::MappingCollection<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &mapping,
       const DoFHandlerType                                 &dof_handler,
       std::vector<Point<DoFHandlerType::space_dimension> > &support_points)
      {
        // get the data in the form of the map as above
        std::map<types::global_dof_index,Point<DoFHandlerType::space_dimension> >  x_support_points;
        map_dofs_to_support_points(mapping, dof_handler, x_support_points);

        // now convert from the map to the linear vector. make sure every
        // entry really appeared in the map
        for (types::global_dof_index i=0; i<dof_handler.n_dofs(); ++i)
          {
            Assert (x_support_points.find(i) != x_support_points.end(),
                    ExcInternalError());
            support_points[i] = x_support_points[i];
          }
      }
    }
  }

  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const Mapping<dim,spacedim>    &mapping,
                              const DoFHandler<dim,spacedim> &dof_handler,
                              std::vector<Point<spacedim> >  &support_points)
  {
    AssertDimension(support_points.size(), dof_handler.n_dofs());
    Assert ((dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
             (&dof_handler.get_triangulation())
             ==
             0),
            ExcMessage ("This function can not be used with distributed triangulations."
                        "See the documentation for more information."));

    // Let the internal function do all the work, just make sure that it
    // gets a MappingCollection
    const hp::MappingCollection<dim, spacedim> mapping_collection(mapping);

    internal::map_dofs_to_support_points (mapping_collection,
                                          dof_handler,
                                          support_points);
  }


  template<int dim, int spacedim>
  void
  map_dofs_to_support_points(const hp::MappingCollection<dim, spacedim> &mapping,
                             const hp::DoFHandler<dim, spacedim>        &dof_handler,
                             std::vector<Point<spacedim> >              &support_points)
  {
    AssertDimension(support_points.size(), dof_handler.n_dofs());
    Assert ((dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
             (&dof_handler.get_triangulation())
             ==
             0),
            ExcMessage ("This function can not be used with distributed triangulations."
                        "See the documentation for more information."));

    // Let the internal function do all the work, just make sure that it
    // gets a MappingCollection
    internal::map_dofs_to_support_points (mapping,
                                          dof_handler,
                                          support_points);
  }


  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const Mapping<dim,spacedim>    &mapping,
                              const DoFHandler<dim,spacedim> &dof_handler,
                              std::map<types::global_dof_index, Point<spacedim> > &support_points)
  {
    support_points.clear();

    // Let the internal function do all the work, just make sure that it
    // gets a MappingCollection
    const hp::MappingCollection<dim, spacedim> mapping_collection(mapping);

    internal::map_dofs_to_support_points (mapping_collection,
                                          dof_handler,
                                          support_points);
  }


  template<int dim, int spacedim>
  void
  map_dofs_to_support_points(const hp::MappingCollection<dim, spacedim> &mapping,
                             const hp::DoFHandler<dim, spacedim>        &dof_handler,
                             std::map<types::global_dof_index, Point<spacedim> > &support_points)
  {
    support_points.clear();

    // Let the internal function do all the work, just make sure that it
    // gets a MappingCollection
    internal::map_dofs_to_support_points (mapping,
                                          dof_handler,
                                          support_points);
  }


  template<int dim, int spacedim>
  void
  convert_couplings_to_blocks
  (const DoFHandler<dim,spacedim>  &dof_handler,
   const Table<2, Coupling>        &table,
   std::vector<Table<2,Coupling> > &tables_by_block)
  {
    const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();
    const unsigned int nb = fe.n_blocks();

    tables_by_block.resize(1);
    tables_by_block[0].reinit(nb, nb);
    tables_by_block[0].fill(none);

    for (unsigned int i=0; i<fe.n_components(); ++i)
      {
        const unsigned int ib = fe.component_to_block_index(i);
        for (unsigned int j=0; j<fe.n_components(); ++j)
          {
            const unsigned int jb = fe.component_to_block_index(j);
            tables_by_block[0](ib,jb) |= table(i,j);
          }
      }
  }


  template<int dim, int spacedim>
  void
  convert_couplings_to_blocks (const hp::DoFHandler<dim,spacedim> &dof_handler,
                               const Table<2, Coupling>           &table,
                               std::vector<Table<2,Coupling> >    &tables_by_block)
  {
    const hp::FECollection<dim> &fe_collection = dof_handler.get_fe();
    tables_by_block.resize(fe_collection.size());

    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
        const FiniteElement<dim,spacedim> &fe = fe_collection[f];

        const unsigned int nb = fe.n_blocks();
        tables_by_block[f].reinit(nb, nb);
        tables_by_block[f].fill(none);
        for (unsigned int i=0; i<fe.n_components(); ++i)
          {
            const unsigned int ib = fe.component_to_block_index(i);
            for (unsigned int j=0; j<fe.n_components(); ++j)
              {
                const unsigned int jb = fe.component_to_block_index(j);
                tables_by_block[f](ib,jb) |= table(i,j);
              }
          }
      }
  }



  template <typename DoFHandlerType, class Sparsity>
  void make_cell_patches(Sparsity                &block_list,
                         const DoFHandlerType    &dof_handler,
                         const unsigned int       level,
                         const std::vector<bool> &selected_dofs,
                         types::global_dof_index  offset)
  {
    typename DoFHandlerType::level_cell_iterator cell;
    typename DoFHandlerType::level_cell_iterator endc = dof_handler.end(level);
    std::vector<types::global_dof_index> indices;

    unsigned int i=0;
    for (cell=dof_handler.begin(level); cell != endc; ++i, ++cell)
      {
        indices.resize(cell->get_fe().dofs_per_cell);
        cell->get_mg_dof_indices(indices);

        if (selected_dofs.size()!=0)
          AssertDimension(indices.size(), selected_dofs.size());

        for (types::global_dof_index j=0; j<indices.size(); ++j)
          {
            if (selected_dofs.size() == 0)
              block_list.add(i,indices[j]-offset);
            else
              {
                if (selected_dofs[j])
                  block_list.add(i,indices[j]-offset);
              }
          }
      }
  }


  template <typename DoFHandlerType>
  void make_single_patch(SparsityPattern      &block_list,
                         const DoFHandlerType &dof_handler,
                         const unsigned int    level,
                         const bool            interior_only)
  {
    const FiniteElement<DoFHandlerType::dimension> &fe = dof_handler.get_fe();
    block_list.reinit(1, dof_handler.n_dofs(level), dof_handler.n_dofs(level));
    typename DoFHandlerType::level_cell_iterator cell;
    typename DoFHandlerType::level_cell_iterator endc = dof_handler.end(level);

    std::vector<types::global_dof_index> indices;
    std::vector<bool> exclude;

    for (cell=dof_handler.begin(level); cell != endc; ++cell)
      {
        indices.resize(cell->get_fe().dofs_per_cell);
        cell->get_mg_dof_indices(indices);

        if (interior_only)
          {
            // Exclude degrees of freedom on faces opposite to the vertex
            exclude.resize(fe.dofs_per_cell);
            std::fill(exclude.begin(), exclude.end(), false);
            const unsigned int dpf = fe.dofs_per_face;

            for (unsigned int face=0; face<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++face)
              if (cell->at_boundary(face) || cell->neighbor(face)->level() != cell->level())
                for (unsigned int i=0; i<dpf; ++i)
                  exclude[fe.face_to_cell_index(i,face)] = true;
            for (types::global_dof_index j=0; j<indices.size(); ++j)
              if (!exclude[j])
                block_list.add(0, indices[j]);
          }
        else
          {
            for (types::global_dof_index j=0; j<indices.size(); ++j)
              block_list.add(0, indices[j]);
          }
      }
  }


  template <typename DoFHandlerType>
  void make_child_patches (SparsityPattern      &block_list,
                           const DoFHandlerType &dof_handler,
                           const unsigned int    level,
                           const bool            interior_dofs_only,
                           const bool            boundary_dofs)
  {
    Assert(level > 0 && level < dof_handler.get_triangulation().n_levels(),
           ExcIndexRange(level, 1, dof_handler.get_triangulation().n_levels()));

    typename DoFHandlerType::level_cell_iterator pcell = dof_handler.begin(level-1);
    typename DoFHandlerType::level_cell_iterator endc = dof_handler.end(level-1);

    std::vector<types::global_dof_index> indices;
    std::vector<bool> exclude;

    for (unsigned int block = 0; pcell != endc; ++pcell)
      {
        if (!pcell->has_children())
          continue;

        for (unsigned int child=0; child<pcell->n_children(); ++child)
          {
            const typename DoFHandlerType::level_cell_iterator cell = pcell->child(child);

            // For hp, only this line here would have to be replaced.
            const FiniteElement<DoFHandlerType::dimension> &fe = dof_handler.get_fe();
            const unsigned int n_dofs = fe.dofs_per_cell;
            indices.resize(n_dofs);
            exclude.resize(n_dofs);
            std::fill(exclude.begin(), exclude.end(), false);
            cell->get_mg_dof_indices(indices);

            if (interior_dofs_only)
              {
                // Eliminate dofs on faces of the child which are on faces
                // of the parent
                const unsigned int dpf = fe.dofs_per_face;

                for (unsigned int d=0; d<DoFHandlerType::dimension; ++d)
                  {
                    const unsigned int face = GeometryInfo<DoFHandlerType::dimension>::vertex_to_face[child][d];
                    for (unsigned int i=0; i<dpf; ++i)
                      exclude[fe.face_to_cell_index(i,face)] = true;
                  }

                // Now remove all degrees of freedom on the domain boundary
                // from the exclusion list
                if (boundary_dofs)
                  for (unsigned int face=0; face< GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++face)
                    if (cell->at_boundary(face))
                      for (unsigned int i=0; i<dpf; ++i)
                        exclude[fe.face_to_cell_index(i,face)] = false;
              }

            for (unsigned int i=0; i<n_dofs; ++i)
              if (!exclude[i])
                block_list.add(block, indices[i]);
          }
        ++block;
      }
  }


  template <typename DoFHandlerType>
  std::vector<unsigned int>
  make_vertex_patches (SparsityPattern      &block_list,
                       const DoFHandlerType &dof_handler,
                       const unsigned int    level,
                       const bool            interior_only,
                       const bool            boundary_patches,
                       const bool            level_boundary_patches,
                       const bool            single_cell_patches,
                       const bool            invert_vertex_mapping)
  {
    typename DoFHandlerType::level_cell_iterator cell;
    typename DoFHandlerType::level_cell_iterator endc = dof_handler.end(level);

    // Vector mapping from vertex index in the triangulation to consecutive
    // block indices on this level The number of cells at a vertex
    std::vector<unsigned int> vertex_cell_count(dof_handler.get_triangulation().n_vertices(), 0);

    // Is a vertex at the boundary?
    std::vector<bool> vertex_boundary(dof_handler.get_triangulation().n_vertices(), false);

    std::vector<unsigned int> vertex_mapping(dof_handler.get_triangulation().n_vertices(),
                                             numbers::invalid_unsigned_int);

    // Estimate for the number of dofs at this point
    std::vector<unsigned int> vertex_dof_count(dof_handler.get_triangulation().n_vertices(), 0);

    // Identify all vertices active on this level and remember some data
    // about them
    for (cell=dof_handler.begin(level); cell != endc; ++cell)
      for (unsigned int v=0; v<GeometryInfo<DoFHandlerType::dimension>::vertices_per_cell; ++v)
        {
          const unsigned int vg = cell->vertex_index(v);
          vertex_dof_count[vg] += cell->get_fe().dofs_per_cell;
          ++vertex_cell_count[vg];
          for (unsigned int d=0; d<DoFHandlerType::dimension; ++d)
            {
              const unsigned int face = GeometryInfo<DoFHandlerType::dimension>::vertex_to_face[v][d];
              if (cell->at_boundary(face))
                vertex_boundary[vg] = true;
              else if ((!level_boundary_patches)
                       && (cell->neighbor(face)->level() != (int) level))
                vertex_boundary[vg] = true;
            }
        }
    // From now on, only vertices with positive dof count are "in".

    // Remove vertices at boundaries or in corners
    for (unsigned int vg=0; vg<vertex_dof_count.size(); ++vg)
      if ((!single_cell_patches && vertex_cell_count[vg] < 2)
          ||
          (!boundary_patches && vertex_boundary[vg]))
        vertex_dof_count[vg] = 0;

    // Create a mapping from all vertices to the ones used here
    unsigned int n_vertex_count=0;
    for (unsigned int vg=0; vg<vertex_mapping.size(); ++vg)
      if (vertex_dof_count[vg] != 0)
        vertex_mapping[vg] = n_vertex_count++;

    // Compactify dof count
    for (unsigned int vg=0; vg<vertex_mapping.size(); ++vg)
      if (vertex_dof_count[vg] != 0)
        vertex_dof_count[vertex_mapping[vg]] = vertex_dof_count[vg];

    // Now that we have all the data, we reduce it to the part we actually
    // want
    vertex_dof_count.resize(n_vertex_count);

    // At this point, the list of patches is ready. Now we enter the dofs
    // into the sparsity pattern.
    block_list.reinit(vertex_dof_count.size(), dof_handler.n_dofs(level), vertex_dof_count);

    std::vector<types::global_dof_index> indices;
    std::vector<bool> exclude;

    for (cell=dof_handler.begin(level); cell != endc; ++cell)
      {
        const FiniteElement<DoFHandlerType::dimension> &fe = cell->get_fe();
        indices.resize(fe.dofs_per_cell);
        cell->get_mg_dof_indices(indices);

        for (unsigned int v=0; v<GeometryInfo<DoFHandlerType::dimension>::vertices_per_cell; ++v)
          {
            const unsigned int vg = cell->vertex_index(v);
            const unsigned int block = vertex_mapping[vg];
            if (block == numbers::invalid_unsigned_int)
              continue;

            if (interior_only)
              {
                // Exclude degrees of freedom on faces opposite to the
                // vertex
                exclude.resize(fe.dofs_per_cell);
                std::fill(exclude.begin(), exclude.end(), false);
                const unsigned int dpf = fe.dofs_per_face;

                for (unsigned int d=0; d<DoFHandlerType::dimension; ++d)
                  {
                    const unsigned int a_face = GeometryInfo<DoFHandlerType::dimension>::vertex_to_face[v][d];
                    const unsigned int face = GeometryInfo<DoFHandlerType::dimension>::opposite_face[a_face];
                    for (unsigned int i=0; i<dpf; ++i)
                      exclude[fe.face_to_cell_index(i,face)] = true;
                  }
                for (unsigned int j=0; j<indices.size(); ++j)
                  if (!exclude[j])
                    block_list.add(block, indices[j]);
              }
            else
              {
                for (unsigned int j=0; j<indices.size(); ++j)
                  block_list.add(block, indices[j]);
              }
          }
      }

    if (invert_vertex_mapping)
      {
        // Compress vertex mapping
        unsigned int n_vertex_count = 0;
        for (unsigned int vg = 0; vg < vertex_mapping.size(); ++vg)
          if (vertex_mapping[vg] != numbers::invalid_unsigned_int)
            vertex_mapping[n_vertex_count++] = vg;

        // Now we reduce it to the part we actually want
        vertex_mapping.resize(n_vertex_count);
      }

    return vertex_mapping;
  }


  template <typename DoFHandlerType>
  unsigned int
  count_dofs_on_patch (const std::vector<typename DoFHandlerType::active_cell_iterator> &patch)
  {
    std::set<types::global_dof_index> dofs_on_patch;
    std::vector<types::global_dof_index> local_dof_indices;

    // loop over the cells in the patch and get the DoFs on each.
    // add all of them to a std::set which automatically makes sure
    // all duplicates are ignored
    for (unsigned int i=0; i<patch.size(); ++i)
      {
        const typename DoFHandlerType::active_cell_iterator cell = patch[i];
        Assert (cell->is_artificial() == false,
                ExcMessage("This function can not be called with cells that are "
                           "not either locally owned or ghost cells."));
        local_dof_indices.resize (cell->get_fe().dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);
        dofs_on_patch.insert (local_dof_indices.begin(),
                              local_dof_indices.end());
      }

    // now return the number of DoFs (duplicates were ignored)
    return dofs_on_patch.size();
  }



  template <typename DoFHandlerType>
  std::vector<types::global_dof_index>
  get_dofs_on_patch (const std::vector<typename DoFHandlerType::active_cell_iterator> &patch)
  {
    std::set<types::global_dof_index> dofs_on_patch;
    std::vector<types::global_dof_index> local_dof_indices;

    // loop over the cells in the patch and get the DoFs on each.
    // add all of them to a std::set which automatically makes sure
    // all duplicates are ignored
    for (unsigned int i=0; i<patch.size(); ++i)
      {
        const typename DoFHandlerType::active_cell_iterator cell = patch[i];
        Assert (cell->is_artificial() == false,
                ExcMessage("This function can not be called with cells that are "
                           "not either locally owned or ghost cells."));
        local_dof_indices.resize (cell->get_fe().dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);
        dofs_on_patch.insert (local_dof_indices.begin(),
                              local_dof_indices.end());
      }

    Assert (dofs_on_patch.size() == count_dofs_on_patch<DoFHandlerType>(patch),
            ExcInternalError());

    // return a vector with the content of the set above. copying
    // also ensures that we retain sortedness as promised in the
    // documentation and as necessary to retain the block structure
    // also on the local system
    return std::vector<types::global_dof_index> (dofs_on_patch.begin(),
                                                 dofs_on_patch.end());
  }


} // end of namespace DoFTools



// explicit instantiations

#include "dof_tools.inst"



DEAL_II_NAMESPACE_CLOSE
