// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_fe_tools_templates_H
#define dealii_fe_tools_templates_H


#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_rannacher_turek.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_rt_bubbles.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>

#include <cctype>
#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace FETools
{
  namespace Compositing
  {
    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities,
      const bool do_tensor_product)
    {
      AssertDimension(fes.size(), multiplicities.size());

      unsigned int multiplied_dofs_per_vertex = 0;
      unsigned int multiplied_dofs_per_line   = 0;
      unsigned int multiplied_dofs_per_quad   = 0;
      unsigned int multiplied_dofs_per_hex    = 0;

      unsigned int multiplied_n_components = 0;

      unsigned int degree = 0; // degree is the maximal degree of the components

      unsigned int n_components = 0;
      // Get the number of components from the first given finite element.
      for (unsigned int i = 0; i < fes.size(); i++)
        if (multiplicities[i] > 0)
          {
            n_components = fes[i]->n_components();
            break;
          }

      for (unsigned int i = 0; i < fes.size(); i++)
        if (multiplicities[i] > 0)
          {
            // TODO: the implementation makes the assumption that all faces have
            // the same number of dofs -> don't construct DPO but
            // PrecomputedFiniteElementData
            AssertDimension(fes[i]->n_unique_quads(), 1);

            multiplied_dofs_per_vertex +=
              fes[i]->n_dofs_per_vertex() * multiplicities[i];
            multiplied_dofs_per_line +=
              fes[i]->n_dofs_per_line() * multiplicities[i];
            multiplied_dofs_per_quad +=
              fes[i]->n_dofs_per_quad(0) * multiplicities[i];
            multiplied_dofs_per_hex +=
              fes[i]->n_dofs_per_hex() * multiplicities[i];

            multiplied_n_components +=
              fes[i]->n_components() * multiplicities[i];

            Assert(do_tensor_product ||
                     (n_components == fes[i]->n_components()),
                   ExcDimensionMismatch(n_components, fes[i]->n_components()));

            degree = std::max(degree, fes[i]->tensor_degree());
          }

      // assume conformity of the first finite element and then take away
      // bits as indicated by the base elements. if all multiplicities
      // happen to be zero, then it doesn't matter what we set it to.
      typename FiniteElementData<dim>::Conformity total_conformity =
        typename FiniteElementData<dim>::Conformity();
      {
        unsigned int index = 0;
        for (index = 0; index < fes.size(); ++index)
          if (multiplicities[index] > 0)
            {
              total_conformity = fes[index]->conforming_space;
              break;
            }

        for (; index < fes.size(); ++index)
          if (multiplicities[index] > 0)
            total_conformity = typename FiniteElementData<dim>::Conformity(
              total_conformity & fes[index]->conforming_space);
      }

      std::vector<unsigned int> dpo;
      dpo.push_back(multiplied_dofs_per_vertex);
      dpo.push_back(multiplied_dofs_per_line);
      if (dim > 1)
        dpo.push_back(multiplied_dofs_per_quad);
      if (dim > 2)
        dpo.push_back(multiplied_dofs_per_hex);

      BlockIndices block_indices(0, 0);

      for (unsigned int base = 0; base < fes.size(); ++base)
        for (unsigned int m = 0; m < multiplicities[base]; ++m)
          block_indices.push_back(fes[base]->n_dofs_per_cell());

      return FiniteElementData<dim>(dpo,
                                    fes.front()->reference_cell(),
                                    (do_tensor_product ?
                                       multiplied_n_components :
                                       n_components),
                                    degree,
                                    total_conformity,
                                    block_indices);
    }



    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(const FiniteElement<dim, spacedim> *fe1,
                         const unsigned int                  N1,
                         const FiniteElement<dim, spacedim> *fe2,
                         const unsigned int                  N2,
                         const FiniteElement<dim, spacedim> *fe3,
                         const unsigned int                  N3,
                         const FiniteElement<dim, spacedim> *fe4,
                         const unsigned int                  N4,
                         const FiniteElement<dim, spacedim> *fe5,
                         const unsigned int                  N5)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fes;
      fes.push_back(fe1);
      fes.push_back(fe2);
      fes.push_back(fe3);
      fes.push_back(fe4);
      fes.push_back(fe5);

      std::vector<unsigned int> mult;
      mult.push_back(N1);
      mult.push_back(N2);
      mult.push_back(N3);
      mult.push_back(N4);
      mult.push_back(N5);
      return multiply_dof_numbers(fes, mult);
    }



    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities)
    {
      AssertDimension(fes.size(), multiplicities.size());

      // first count the number of dofs and components that will emerge from the
      // given FEs
      unsigned int n_shape_functions = 0;
      for (unsigned int i = 0; i < fes.size(); ++i)
        if (multiplicities[i] > 0) // check needed as FE might be nullptr
          n_shape_functions += fes[i]->n_dofs_per_cell() * multiplicities[i];

      // generate the array that will hold the output
      std::vector<bool> retval(n_shape_functions, false);

      // finally go through all the shape functions of the base elements, and
      // copy their flags. this somehow copies the code in build_cell_table,
      // which is not nice as it uses too much implicit knowledge about the
      // layout of the individual bases in the composed FE, but there seems no
      // way around...
      //
      // for each shape function, copy the flags from the base element to this
      // one, taking into account multiplicities, and other complications
      unsigned int total_index = 0;
      for (const unsigned int vertex_number :
           fes.front()->reference_cell().vertex_indices())
        {
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base]; ++m)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_vertex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fes[base]->n_dofs_per_vertex() * vertex_number +
                     local_index);

                  Assert(index_in_base < fes[base]->n_dofs_per_cell(),
                         ExcInternalError());
                  retval[total_index] =
                    fes[base]->restriction_is_additive(index_in_base);
                }
        }

      // 2. Lines
      for (const unsigned int line_number :
           fes.front()->reference_cell().line_indices())
        {
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base]; ++m)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_line();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fes[base]->n_dofs_per_line() * line_number + local_index +
                     fes[base]->get_first_line_index());

                  Assert(index_in_base < fes[base]->n_dofs_per_cell(),
                         ExcInternalError());
                  retval[total_index] =
                    fes[base]->restriction_is_additive(index_in_base);
                }
        }

      // 3. Quads
      for (unsigned int quad_number = 0;
           quad_number <
           (dim == 2 ?
              1 :
              (dim == 3 ? fes.front()->reference_cell().n_faces() : 0));
           ++quad_number)
        {
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base]; ++m)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_quad(quad_number);
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index + fes[base]->get_first_quad_index(quad_number);

                  Assert(index_in_base < fes[base]->n_dofs_per_cell(),
                         ExcInternalError());
                  retval[total_index] =
                    fes[base]->restriction_is_additive(index_in_base);
                }
        }

      // 4. Hexes
      if (dim == 3)
        {
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base]; ++m)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_hex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index + fes[base]->get_first_hex_index();

                  Assert(index_in_base < fes[base]->n_dofs_per_cell(),
                         ExcInternalError());
                  retval[total_index] =
                    fes[base]->restriction_is_additive(index_in_base);
                }
        }

      Assert(total_index == n_shape_functions, ExcInternalError());

      return retval;
    }



    /**
     * Take a @p FiniteElement object
     * and return an boolean vector including the @p
     * restriction_is_additive_flags of the mixed element consisting of @p N
     * elements of the sub-element @p fe.
     */
    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const FiniteElement<dim, spacedim> *fe1,
      const unsigned int                  N1,
      const FiniteElement<dim, spacedim> *fe2,
      const unsigned int                  N2,
      const FiniteElement<dim, spacedim> *fe3,
      const unsigned int                  N3,
      const FiniteElement<dim, spacedim> *fe4,
      const unsigned int                  N4,
      const FiniteElement<dim, spacedim> *fe5,
      const unsigned int                  N5)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fe_list;
      std::vector<unsigned int>                         multiplicities;

      fe_list.push_back(fe1);
      multiplicities.push_back(N1);

      fe_list.push_back(fe2);
      multiplicities.push_back(N2);

      fe_list.push_back(fe3);
      multiplicities.push_back(N3);

      fe_list.push_back(fe4);
      multiplicities.push_back(N4);

      fe_list.push_back(fe5);
      multiplicities.push_back(N5);
      return compute_restriction_is_additive_flags(fe_list, multiplicities);
    }



    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities,
      const bool do_tensor_product)
    {
      AssertDimension(fes.size(), multiplicities.size());
      Assert(
        fes.size() > 0,
        ExcMessage(
          "This function only makes sense if at least one FiniteElement is provided."));
      // first count the number of dofs and components that will emerge from the
      // given FEs
      unsigned int n_shape_functions = 0;
      for (unsigned int i = 0; i < fes.size(); ++i)
        if (multiplicities[i] > 0) // needed because FE might be nullptr
          n_shape_functions += fes[i]->n_dofs_per_cell() * multiplicities[i];

      unsigned int n_components = 0;
      if (do_tensor_product)
        {
          for (unsigned int i = 0; i < fes.size(); ++i)
            if (multiplicities[i] > 0) // needed because FE might be nullptr
              n_components += fes[i]->n_components() * multiplicities[i];
        }
      else
        {
          for (unsigned int i = 0; i < fes.size(); ++i)
            if (multiplicities[i] > 0) // needed because FE might be nullptr
              {
                n_components = fes[i]->n_components();
                break;
              }
          // Now check that all FEs have the same number of components:
          for (unsigned int i = 0; i < fes.size(); ++i)
            if (multiplicities[i] > 0) // needed because FE might be nullptr
              Assert(n_components == fes[i]->n_components(),
                     ExcDimensionMismatch(n_components,
                                          fes[i]->n_components()));
        }

      // generate the array that will hold the output
      std::vector<std::vector<bool>> retval(n_shape_functions,
                                            std::vector<bool>(n_components,
                                                              false));

      // finally go through all the shape functions of the base elements, and
      // copy their flags. this somehow copies the code in build_cell_table,
      // which is not nice as it uses too much implicit knowledge about the
      // layout of the individual bases in the composed FE, but there seems no
      // way around...
      //
      // for each shape function, copy the non-zero flags from the base element
      // to this one, taking into account multiplicities, multiple components in
      // base elements, and other complications
      unsigned int total_index = 0;
      for (const unsigned int vertex_number :
           fes.front()->reference_cell().vertex_indices())
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base];
                 ++m,
                              comp_start +=
                              fes[base]->n_components() * do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_vertex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fes[base]->n_dofs_per_vertex() * vertex_number +
                     local_index);

                  Assert(comp_start + fes[base]->n_components() <=
                           retval[total_index].size(),
                         ExcInternalError());
                  for (unsigned int c = 0; c < fes[base]->n_components(); ++c)
                    {
                      Assert(c < fes[base]
                                   ->get_nonzero_components(index_in_base)
                                   .size(),
                             ExcInternalError());
                      retval[total_index][comp_start + c] =
                        fes[base]->get_nonzero_components(index_in_base)[c];
                    }
                }
        }

      // 2. Lines
      for (const unsigned int line_number :
           fes.front()->reference_cell().line_indices())
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base];
                 ++m,
                              comp_start +=
                              fes[base]->n_components() * do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_line();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fes[base]->n_dofs_per_line() * line_number + local_index +
                     fes[base]->get_first_line_index());

                  Assert(comp_start + fes[base]->n_components() <=
                           retval[total_index].size(),
                         ExcInternalError());
                  for (unsigned int c = 0; c < fes[base]->n_components(); ++c)
                    {
                      Assert(c < fes[base]
                                   ->get_nonzero_components(index_in_base)
                                   .size(),
                             ExcInternalError());
                      retval[total_index][comp_start + c] =
                        fes[base]->get_nonzero_components(index_in_base)[c];
                    }
                }
        }

      // 3. Quads
      for (unsigned int quad_number = 0;
           quad_number <
           (dim == 2 ?
              1 :
              (dim == 3 ? fes.front()->reference_cell().n_faces() : 0));
           ++quad_number)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base];
                 ++m,
                              comp_start +=
                              fes[base]->n_components() * do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_quad(quad_number);
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index + fes[base]->get_first_quad_index(quad_number);

                  Assert(comp_start + fes[base]->n_components() <=
                           retval[total_index].size(),
                         ExcInternalError());
                  for (unsigned int c = 0; c < fes[base]->n_components(); ++c)
                    {
                      Assert(c < fes[base]
                                   ->get_nonzero_components(index_in_base)
                                   .size(),
                             ExcInternalError());
                      retval[total_index][comp_start + c] =
                        fes[base]->get_nonzero_components(index_in_base)[c];
                    }
                }
        }

      // 4. Hexes
      if (dim == 3)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fes.size(); ++base)
            for (unsigned int m = 0; m < multiplicities[base];
                 ++m,
                              comp_start +=
                              fes[base]->n_components() * do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fes[base]->n_dofs_per_hex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index + fes[base]->get_first_hex_index();

                  Assert(comp_start + fes[base]->n_components() <=
                           retval[total_index].size(),
                         ExcInternalError());
                  for (unsigned int c = 0; c < fes[base]->n_components(); ++c)
                    {
                      Assert(c < fes[base]
                                   ->get_nonzero_components(index_in_base)
                                   .size(),
                             ExcInternalError());
                      retval[total_index][comp_start + c] =
                        fes[base]->get_nonzero_components(index_in_base)[c];
                    }
                }
        }

      Assert(total_index == n_shape_functions, ExcInternalError());

      // now copy the vector<vector<bool> > into a vector<ComponentMask>.
      // this appears complicated but we do it this way since it's just
      // awkward to generate ComponentMasks directly and so we need the
      // recourse of the inner vector<bool> anyway.
      std::vector<ComponentMask> xretval(retval.size());
      for (unsigned int i = 0; i < retval.size(); ++i)
        xretval[i] = ComponentMask(retval[i]);
      return xretval;
    }



    /**
     * Compute the non-zero vector components of a composed finite element.
     */
    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(const FiniteElement<dim, spacedim> *fe1,
                               const unsigned int                  N1,
                               const FiniteElement<dim, spacedim> *fe2,
                               const unsigned int                  N2,
                               const FiniteElement<dim, spacedim> *fe3,
                               const unsigned int                  N3,
                               const FiniteElement<dim, spacedim> *fe4,
                               const unsigned int                  N4,
                               const FiniteElement<dim, spacedim> *fe5,
                               const unsigned int                  N5,
                               const bool do_tensor_product)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fe_list;
      std::vector<unsigned int>                         multiplicities;

      fe_list.push_back(fe1);
      multiplicities.push_back(N1);

      fe_list.push_back(fe2);
      multiplicities.push_back(N2);

      fe_list.push_back(fe3);
      multiplicities.push_back(N3);

      fe_list.push_back(fe4);
      multiplicities.push_back(N4);

      fe_list.push_back(fe5);
      multiplicities.push_back(N5);

      return compute_nonzero_components(fe_list,
                                        multiplicities,
                                        do_tensor_product);
    }



    template <int dim, int spacedim>
    void
    build_cell_tables(
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &system_to_base_table,
      std::vector<std::pair<unsigned int, unsigned int>>
        &                                   system_to_component_table,
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &component_to_base_table,
      const FiniteElement<dim, spacedim> &  fe,
      const bool                            do_tensor_product)
    {
      unsigned int total_index = 0;

      if (do_tensor_product)
        {
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base); ++m)
              {
                for (unsigned int k = 0;
                     k < fe.base_element(base).n_components();
                     ++k)
                  component_to_base_table[total_index++] =
                    std::make_pair(std::make_pair(base, k), m);
              }
          Assert(total_index == component_to_base_table.size(),
                 ExcInternalError());
        }
      else
        {
          // The base element establishing a component does not make sense in
          // this case. Set up to something meaningless:
          std::fill(
            component_to_base_table.begin(),
            component_to_base_table.end(),
            std::make_pair(std::make_pair(numbers::invalid_unsigned_int,
                                          numbers::invalid_unsigned_int),
                           numbers::invalid_unsigned_int));
        }


      // Initialize index tables.  Multi-component base elements have to be
      // thought of. For non-primitive shape functions, have a special invalid
      // index.
      const std::pair<unsigned int, unsigned int> non_primitive_index(
        numbers::invalid_unsigned_int, numbers::invalid_unsigned_int);

      // First enumerate vertex indices, where we first enumerate all indices on
      // the first vertex in the order of the base elements, then of the second
      // vertex, etc
      total_index = 0;
      for (const unsigned int vertex_number :
           fe.reference_cell().vertex_indices())
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_vertex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fe.base_element(base).n_dofs_per_vertex() * vertex_number +
                     local_index);

                  system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int index_in_comp =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .second;
                      system_to_component_table[total_index] =
                        std::make_pair(comp, index_in_comp);
                    }
                  else
                    system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }

      // 2. Lines
      for (const unsigned int line_number : fe.reference_cell().line_indices())
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_line();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    (fe.base_element(base).n_dofs_per_line() * line_number +
                     local_index +
                     fe.base_element(base).get_first_line_index());

                  system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int index_in_comp =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .second;
                      system_to_component_table[total_index] =
                        std::make_pair(comp, index_in_comp);
                    }
                  else
                    system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }

      // 3. Quads
      for (unsigned int quad_number = 0;
           quad_number <
           (dim == 2 ? 1 : (dim == 3 ? fe.reference_cell().n_faces() : 0));
           ++quad_number)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index <
                   fe.base_element(base).n_dofs_per_quad(quad_number);
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index +
                    fe.base_element(base).get_first_quad_index(quad_number);

                  system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int index_in_comp =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .second;
                      system_to_component_table[total_index] =
                        std::make_pair(comp, index_in_comp);
                    }
                  else
                    system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }

      // 4. Hexes
      if (dim == 3)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_hex();
                   ++local_index, ++total_index)
                {
                  const unsigned int index_in_base =
                    local_index + fe.base_element(base).get_first_hex_index();

                  system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int index_in_comp =
                        fe.base_element(base)
                          .system_to_component_index(index_in_base)
                          .second;
                      system_to_component_table[total_index] =
                        std::make_pair(comp, index_in_comp);
                    }
                  else
                    system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }
    }



    template <int dim, int spacedim>
    void
    build_face_tables(
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &face_system_to_base_table,
      std::vector<std::pair<unsigned int, unsigned int>>
        &                                 face_system_to_component_table,
      const FiniteElement<dim, spacedim> &fe,
      const bool                          do_tensor_product,
      const unsigned int                  face_no)
    {
      // Initialize index tables. do this in the same way as done for the cell
      // tables, except that we now loop over the objects of faces

      // For non-primitive shape functions, have a special invalid index
      const std::pair<unsigned int, unsigned int> non_primitive_index(
        numbers::invalid_unsigned_int, numbers::invalid_unsigned_int);

      // 1. Vertices
      unsigned int total_index = 0;
      for (unsigned int vertex_number = 0;
           vertex_number <
           fe.reference_cell().face_reference_cell(face_no).n_vertices();
           ++vertex_number)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_vertex();
                   ++local_index, ++total_index)
                {
                  // get (cell) index of this shape function inside the base
                  // element to see whether the shape function is primitive
                  // (assume that all shape functions on vertices share the same
                  // primitivity property; assume likewise for all shape
                  // functions located on lines, quads, etc. this way, we can
                  // ask for primitivity of only _one_ shape function, which is
                  // taken as representative for all others located on the same
                  // type of object):
                  const unsigned int index_in_base =
                    (fe.base_element(base).n_dofs_per_vertex() * vertex_number +
                     local_index);

                  const unsigned int face_index_in_base =
                    (fe.base_element(base).n_dofs_per_vertex() * vertex_number +
                     local_index);

                  face_system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), face_index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int face_index_in_comp =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .second;
                      face_system_to_component_table[total_index] =
                        std::make_pair(comp, face_index_in_comp);
                    }
                  else
                    face_system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }

      // 2. Lines
      for (unsigned int line_number = 0;
           line_number <
           fe.reference_cell().face_reference_cell(face_no).n_lines();
           ++line_number)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_line();
                   ++local_index, ++total_index)
                {
                  // do everything alike for this type of object
                  const unsigned int index_in_base =
                    (fe.base_element(base).n_dofs_per_line() * line_number +
                     local_index +
                     fe.base_element(base).get_first_line_index());

                  const unsigned int face_index_in_base =
                    (fe.base_element(base).get_first_face_line_index(face_no) +
                     fe.base_element(base).n_dofs_per_line() * line_number +
                     local_index);

                  face_system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), face_index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int face_index_in_comp =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .second;
                      face_system_to_component_table[total_index] =
                        std::make_pair(comp, face_index_in_comp);
                    }
                  else
                    face_system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }

      // 3. Quads
      if (dim == 3)
        {
          unsigned int comp_start = 0;
          for (unsigned int base = 0; base < fe.n_base_elements(); ++base)
            for (unsigned int m = 0; m < fe.element_multiplicity(base);
                 ++m,
                              comp_start +=
                              fe.base_element(base).n_components() *
                              do_tensor_product)
              for (unsigned int local_index = 0;
                   local_index < fe.base_element(base).n_dofs_per_quad(face_no);
                   ++local_index, ++total_index)
                {
                  // do everything alike for this type of object
                  const unsigned int index_in_base =
                    (local_index +
                     fe.base_element(base).get_first_quad_index(face_no));

                  const unsigned int face_index_in_base =
                    (fe.base_element(base).get_first_face_quad_index(face_no) +
                     local_index);

                  face_system_to_base_table[total_index] =
                    std::make_pair(std::make_pair(base, m), face_index_in_base);

                  if (fe.base_element(base).is_primitive(index_in_base))
                    {
                      const unsigned int comp_in_base =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .first;
                      const unsigned int comp = comp_start + comp_in_base;
                      const unsigned int face_index_in_comp =
                        fe.base_element(base)
                          .face_system_to_component_index(face_index_in_base,
                                                          face_no)
                          .second;
                      face_system_to_component_table[total_index] =
                        std::make_pair(comp, face_index_in_comp);
                    }
                  else
                    face_system_to_component_table[total_index] =
                      non_primitive_index;
                }
        }
      Assert(total_index == fe.n_dofs_per_face(face_no), ExcInternalError());
      Assert(total_index == face_system_to_component_table.size(),
             ExcInternalError());
      Assert(total_index == face_system_to_base_table.size(),
             ExcInternalError());
    }
  } // namespace Compositing



  // Not implemented in the general case.
  template <class FE>
  std::unique_ptr<FiniteElement<FE::dimension, FE::space_dimension>>
  FEFactory<FE>::get(const Quadrature<1> &) const
  {
    Assert(false, ExcNotImplemented());
    return nullptr;
  }

  // Specializations for FE_Q.
  template <>
  std::unique_ptr<FiniteElement<1, 1>>
  FEFactory<FE_Q<1, 1>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q<1>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<2, 2>>
  FEFactory<FE_Q<2, 2>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q<2>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<3, 3>>
  FEFactory<FE_Q<3, 3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q<3>>(quad);
  }

  // Specializations for FE_Q_DG0.
  template <>
  std::unique_ptr<FiniteElement<1, 1>>
  FEFactory<FE_Q_DG0<1, 1>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_DG0<1>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<2, 2>>
  FEFactory<FE_Q_DG0<2, 2>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_DG0<2>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<3, 3>>
  FEFactory<FE_Q_DG0<3, 3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_DG0<3>>(quad);
  }

  // Specializations for FE_Q_Bubbles.
  template <>
  std::unique_ptr<FiniteElement<1, 1>>
  FEFactory<FE_Q_Bubbles<1, 1>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_Bubbles<1>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<2, 2>>
  FEFactory<FE_Q_Bubbles<2, 2>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_Bubbles<2>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<3, 3>>
  FEFactory<FE_Q_Bubbles<3, 3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_Q_Bubbles<3>>(quad);
  }

  // Specializations for FE_DGQArbitraryNodes.
  template <>
  std::unique_ptr<FiniteElement<1, 1>>
  FEFactory<FE_DGQ<1>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<1>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<1, 2>>
  FEFactory<FE_DGQ<1, 2>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<1, 2>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<1, 3>>
  FEFactory<FE_DGQ<1, 3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<1, 3>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<2, 2>>
  FEFactory<FE_DGQ<2>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<2>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<2, 3>>
  FEFactory<FE_DGQ<2, 3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<2, 3>>(quad);
  }

  template <>
  std::unique_ptr<FiniteElement<3, 3>>
  FEFactory<FE_DGQ<3>>::get(const Quadrature<1> &quad) const
  {
    return std::make_unique<FE_DGQArbitraryNodes<3>>(quad);
  }



  namespace internal
  {
    // The following three functions serve to fill the maps from element
    // names to elements fe_name_map below. The first one exists because
    // we have finite elements which are not implemented for nonzero
    // codimension. These should be transferred to the second function
    // eventually.
    namespace FEToolsAddFENameHelper
    {
      template <int dim>
      void
      fill_no_codim_fe_names(
        std::map<std::string, std::unique_ptr<const Subscriptor>> &result)
      {
        result["FE_Q_Hierarchical"] =
          std::make_unique<FETools::FEFactory<FE_Q_Hierarchical<dim>>>();
        result["FE_ABF"] = std::make_unique<FETools::FEFactory<FE_ABF<dim>>>();
        result["FE_Bernstein"] =
          std::make_unique<FETools::FEFactory<FE_Bernstein<dim>>>();
        result["FE_BDM"] = std::make_unique<FETools::FEFactory<FE_BDM<dim>>>();
        result["FE_DGBDM"] =
          std::make_unique<FETools::FEFactory<FE_DGBDM<dim>>>();
        result["FE_DGNedelec"] =
          std::make_unique<FETools::FEFactory<FE_DGNedelec<dim>>>();
        result["FE_DGRaviartThomas"] =
          std::make_unique<FETools::FEFactory<FE_DGRaviartThomas<dim>>>();
        result["FE_RaviartThomas"] =
          std::make_unique<FETools::FEFactory<FE_RaviartThomas<dim>>>();
        result["FE_RaviartThomasNodal"] =
          std::make_unique<FETools::FEFactory<FE_RaviartThomasNodal<dim>>>();
        result["FE_RT_Bubbles"] =
          std::make_unique<FETools::FEFactory<FE_RT_Bubbles<dim>>>();
        result["FE_Nedelec"] =
          std::make_unique<FETools::FEFactory<FE_Nedelec<dim>>>();
        result["FE_NedelecSZ"] =
          std::make_unique<FETools::FEFactory<FE_NedelecSZ<dim>>>();
        result["FE_DGPNonparametric"] =
          std::make_unique<FETools::FEFactory<FE_DGPNonparametric<dim>>>();
        result["FE_DGP"] = std::make_unique<FETools::FEFactory<FE_DGP<dim>>>();
        result["FE_DGPMonomial"] =
          std::make_unique<FETools::FEFactory<FE_DGPMonomial<dim>>>();
        result["FE_DGQ"] = std::make_unique<FETools::FEFactory<FE_DGQ<dim>>>();
        result["FE_DGQArbitraryNodes"] =
          std::make_unique<FETools::FEFactory<FE_DGQ<dim>>>();
        result["FE_DGQLegendre"] =
          std::make_unique<FETools::FEFactory<FE_DGQLegendre<dim>>>();
        result["FE_DGQHermite"] =
          std::make_unique<FETools::FEFactory<FE_DGQHermite<dim>>>();
        result["FE_FaceQ"] =
          std::make_unique<FETools::FEFactory<FE_FaceQ<dim>>>();
        result["FE_FaceP"] =
          std::make_unique<FETools::FEFactory<FE_FaceP<dim>>>();
        result["FE_Q"] = std::make_unique<FETools::FEFactory<FE_Q<dim>>>();
        result["FE_Q_DG0"] =
          std::make_unique<FETools::FEFactory<FE_Q_DG0<dim>>>();
        result["FE_Q_Bubbles"] =
          std::make_unique<FETools::FEFactory<FE_Q_Bubbles<dim>>>();
        result["FE_Q_iso_Q1"] =
          std::make_unique<FETools::FEFactory<FE_Q_iso_Q1<dim>>>();
        result["FE_Nothing"] =
          std::make_unique<FETools::FEFactory<FE_Nothing<dim>>>();
        result["FE_RannacherTurek"] =
          std::make_unique<FETools::FEFactory<FE_RannacherTurek<dim>>>();
      }



      // This function fills a map from names to finite elements for any
      // dimension and codimension for those elements which support
      // nonzero codimension.
      template <int dim, int spacedim>
      void
      fill_codim_fe_names(
        std::map<std::string, std::unique_ptr<const Subscriptor>> &result)
      {
        result["FE_Bernstein"] =
          std::make_unique<FETools::FEFactory<FE_Bernstein<dim, spacedim>>>();
        result["FE_DGP"] =
          std::make_unique<FETools::FEFactory<FE_DGP<dim, spacedim>>>();
        result["FE_DGQ"] =
          std::make_unique<FETools::FEFactory<FE_DGQ<dim, spacedim>>>();
        result["FE_Nothing"] =
          std::make_unique<FETools::FEFactory<FE_Nothing<dim, spacedim>>>();
        result["FE_DGQArbitraryNodes"] =
          std::make_unique<FETools::FEFactory<FE_DGQ<dim, spacedim>>>();
        result["FE_DGQLegendre"] =
          std::make_unique<FETools::FEFactory<FE_DGQLegendre<dim, spacedim>>>();
        result["FE_DGQHermite"] =
          std::make_unique<FETools::FEFactory<FE_DGQHermite<dim, spacedim>>>();
        result["FE_Q_Bubbles"] =
          std::make_unique<FETools::FEFactory<FE_Q_Bubbles<dim, spacedim>>>();
        result["FE_Q_DG0"] =
          std::make_unique<FETools::FEFactory<FE_Q_DG0<dim, spacedim>>>();
        result["FE_Q_iso_Q1"] =
          std::make_unique<FETools::FEFactory<FE_Q_iso_Q1<dim, spacedim>>>();
        result["FE_Q"] =
          std::make_unique<FETools::FEFactory<FE_Q<dim, spacedim>>>();
        result["FE_Bernstein"] =
          std::make_unique<FETools::FEFactory<FE_Bernstein<dim, spacedim>>>();
      }

      // The function filling the vector fe_name_map below. It iterates
      // through all legal dimension/spacedimension pairs and fills
      // fe_name_map[dimension][spacedimension] with the maps generated
      // by the functions above.
      std::array<
        std::array<std::map<std::string, std::unique_ptr<const Subscriptor>>,
                   4>,
        4>
      fill_default_map()
      {
        std::array<
          std::array<std::map<std::string, std::unique_ptr<const Subscriptor>>,
                     4>,
          4>
          result;

        fill_no_codim_fe_names<1>(result[1][1]);
        fill_no_codim_fe_names<2>(result[2][2]);
        fill_no_codim_fe_names<3>(result[3][3]);

        fill_codim_fe_names<1, 2>(result[1][2]);
        fill_codim_fe_names<1, 3>(result[1][3]);
        fill_codim_fe_names<2, 3>(result[2][3]);

        return result;
      }


      // have a lock that guarantees that at most one thread is changing
      // and accessing the fe_name_map variable. make this lock local to
      // this file.
      //
      // This variable is declared static (even though
      // it belongs to an internal namespace) in order to make icc happy
      // (which otherwise reports a multiply defined symbol when linking
      // libraries for more than one space dimension together)
      static std::mutex fe_name_map_lock;

      // This is the map used by FETools::get_fe_by_name and
      // FETools::add_fe_name. It is only accessed by functions in this
      // file, so it is safe to make it a static variable here. It must be
      // static so that we can link several dimensions together.

      // The organization of this storage is such that
      // fe_name_map[dim][spacedim][name] points to an
      // FEFactoryBase<dim,spacedim> with the name given. Since
      // all entries of this vector are of different type, we store
      // pointers to generic objects and cast them when needed.

      // We use a unique pointer to factory objects, to ensure that they
      // get deleted at the end of the program run and don't end up as
      // apparent memory leaks to programs like valgrind.

      // This vector is initialized at program start time using the
      // function above. because at this time there are no threads
      // running, there are no thread-safety issues here. since this is
      // compiled for all dimensions at once, need to create objects for
      // each dimension and then separate between them further down
      inline std::array<
        std::array<std::map<std::string, std::unique_ptr<const Subscriptor>>,
                   4>,
        4> &
      get_fe_name_map()
      {
        static std::array<
          std::array<std::map<std::string, std::unique_ptr<const Subscriptor>>,
                     4>,
          4>
          fe_name_map = fill_default_map();
        return fe_name_map;
      }
    } // namespace FEToolsAddFENameHelper

    namespace FEToolsGetInterpolationMatrixHelper
    {
      // forwarder function for
      // FE::get_interpolation_matrix. we
      // will want to call that function
      // for arbitrary FullMatrix<T>
      // types, but it only accepts
      // double arguments. since it is a
      // virtual function, this can also
      // not be changed. so have a
      // forwarder function that calls
      // that function directly if
      // T==double, and otherwise uses a
      // temporary
      template <int dim, int spacedim>
      inline void
      gim_forwarder(const FiniteElement<dim, spacedim> &fe1,
                    const FiniteElement<dim, spacedim> &fe2,
                    FullMatrix<double> &                interpolation_matrix)
      {
        fe2.get_interpolation_matrix(fe1, interpolation_matrix);
      }



      template <int dim, typename number, int spacedim>
      inline void
      gim_forwarder(const FiniteElement<dim, spacedim> &fe1,
                    const FiniteElement<dim, spacedim> &fe2,
                    FullMatrix<number> &                interpolation_matrix)
      {
        FullMatrix<double> tmp(interpolation_matrix.m(),
                               interpolation_matrix.n());
        fe2.get_interpolation_matrix(fe1, tmp);
        interpolation_matrix = tmp;
      }
    } // namespace FEToolsGetInterpolationMatrixHelper
  }   // namespace internal



  template <int dim, int spacedim>
  void
  compute_component_wise(const FiniteElement<dim, spacedim> &    element,
                         std::vector<unsigned int> &             renumbering,
                         std::vector<std::vector<unsigned int>> &comp_start)
  {
    Assert(renumbering.size() == element.n_dofs_per_cell(),
           ExcDimensionMismatch(renumbering.size(), element.n_dofs_per_cell()));

    comp_start.resize(element.n_base_elements());

    unsigned int index = 0;
    for (unsigned int i = 0; i < comp_start.size(); ++i)
      {
        comp_start[i].resize(element.element_multiplicity(i));
        const unsigned int increment =
          element.base_element(i).n_dofs_per_cell();

        for (unsigned int &first_index_of_component : comp_start[i])
          {
            first_index_of_component = index;
            index += increment;
          }
      }

    // For each index i of the
    // unstructured cellwise
    // numbering, renumbering
    // contains the index of the
    // cell-block numbering
    for (unsigned int i = 0; i < element.n_dofs_per_cell(); ++i)
      {
        std::pair<std::pair<unsigned int, unsigned int>, unsigned int> indices =
          element.system_to_base_index(i);
        renumbering[i] = comp_start[indices.first.first][indices.first.second] +
                         indices.second;
      }
  }



  template <int dim, int spacedim>
  void
  compute_block_renumbering(const FiniteElement<dim, spacedim> &  element,
                            std::vector<types::global_dof_index> &renumbering,
                            std::vector<types::global_dof_index> &block_data,
                            bool return_start_indices)
  {
    Assert(renumbering.size() == element.n_dofs_per_cell(),
           ExcDimensionMismatch(renumbering.size(), element.n_dofs_per_cell()));
    Assert(block_data.size() == element.n_blocks(),
           ExcDimensionMismatch(block_data.size(), element.n_blocks()));

    types::global_dof_index k     = 0;
    unsigned int            count = 0;
    for (unsigned int b = 0; b < element.n_base_elements(); ++b)
      for (unsigned int m = 0; m < element.element_multiplicity(b); ++m)
        {
          block_data[count++] = (return_start_indices) ?
                                  k :
                                  (element.base_element(b).n_dofs_per_cell());
          k += element.base_element(b).n_dofs_per_cell();
        }
    Assert(count == element.n_blocks(), ExcInternalError());

    std::vector<types::global_dof_index> start_indices(block_data.size());
    k = 0;
    for (unsigned int i = 0; i < block_data.size(); ++i)
      if (return_start_indices)
        start_indices[i] = block_data[i];
      else
        {
          start_indices[i] = k;
          k += block_data[i];
        }

    for (unsigned int i = 0; i < element.n_dofs_per_cell(); ++i)
      {
        std::pair<unsigned int, types::global_dof_index> indices =
          element.system_to_block_index(i);
        renumbering[i] = start_indices[indices.first] + indices.second;
      }
  }



  template <int dim, typename number, int spacedim>
  void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &fe1,
                           const FiniteElement<dim, spacedim> &fe2,
                           FullMatrix<number> &interpolation_matrix)
  {
    Assert(fe1.n_components() == fe2.n_components(),
           ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
    Assert(interpolation_matrix.m() == fe2.n_dofs_per_cell() &&
             interpolation_matrix.n() == fe1.n_dofs_per_cell(),
           ExcMatrixDimensionMismatch(interpolation_matrix.m(),
                                      interpolation_matrix.n(),
                                      fe2.n_dofs_per_cell(),
                                      fe1.n_dofs_per_cell()));

    // first try the easy way: maybe the FE wants to implement things itself:
    try
      {
        internal::FEToolsGetInterpolationMatrixHelper::gim_forwarder(
          fe1, fe2, interpolation_matrix);
        return;
      }
    catch (
      typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented &)
      {
        // too bad....
      }

    // uh, so this was not the
    // case. hm. then do it the hard
    // way. note that this will only
    // work if the element is
    // primitive, so check this first
    Assert(fe1.is_primitive() == true, ExcFENotPrimitive());
    Assert(fe2.is_primitive() == true, ExcFENotPrimitive());

    // Initialize FEValues for fe1 at
    // the unit support points of the
    // fe2 element.
    const std::vector<Point<dim>> &fe2_support_points =
      fe2.get_unit_support_points();

    Assert(fe2_support_points.size() == fe2.n_dofs_per_cell(),
           (typename FiniteElement<dim, spacedim>::ExcFEHasNoSupportPoints()));

    for (unsigned int i = 0; i < fe2.n_dofs_per_cell(); ++i)
      {
        const unsigned int i1 = fe2.system_to_component_index(i).first;
        for (unsigned int j = 0; j < fe1.n_dofs_per_cell(); ++j)
          {
            const unsigned int j1 = fe1.system_to_component_index(j).first;
            if (i1 == j1)
              interpolation_matrix(i, j) =
                fe1.shape_value(j, fe2_support_points[i]);
            else
              interpolation_matrix(i, j) = 0.;
          }
      }
  }



  template <int dim, typename number, int spacedim>
  void
  get_back_interpolation_matrix(const FiniteElement<dim, spacedim> &fe1,
                                const FiniteElement<dim, spacedim> &fe2,
                                FullMatrix<number> &interpolation_matrix)
  {
    Assert(fe1.n_components() == fe2.n_components(),
           ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
    Assert(interpolation_matrix.m() == fe1.n_dofs_per_cell() &&
             interpolation_matrix.n() == fe1.n_dofs_per_cell(),
           ExcMatrixDimensionMismatch(interpolation_matrix.m(),
                                      interpolation_matrix.n(),
                                      fe1.n_dofs_per_cell(),
                                      fe1.n_dofs_per_cell()));

    FullMatrix<number> first_matrix(fe2.n_dofs_per_cell(),
                                    fe1.n_dofs_per_cell());
    FullMatrix<number> second_matrix(fe1.n_dofs_per_cell(),
                                     fe2.n_dofs_per_cell());

    get_interpolation_matrix(fe1, fe2, first_matrix);
    get_interpolation_matrix(fe2, fe1, second_matrix);

    // int_matrix=second_matrix*first_matrix
    second_matrix.mmult(interpolation_matrix, first_matrix);
  }



  template <int dim, typename number, int spacedim>
  void
  get_interpolation_difference_matrix(const FiniteElement<dim, spacedim> &fe1,
                                      const FiniteElement<dim, spacedim> &fe2,
                                      FullMatrix<number> &difference_matrix)
  {
    Assert(fe1.n_components() == fe2.n_components(),
           ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
    Assert(difference_matrix.m() == fe1.n_dofs_per_cell() &&
             difference_matrix.n() == fe1.n_dofs_per_cell(),
           ExcMatrixDimensionMismatch(difference_matrix.m(),
                                      difference_matrix.n(),
                                      fe1.n_dofs_per_cell(),
                                      fe1.n_dofs_per_cell()));

    FullMatrix<number> interpolation_matrix(fe1.n_dofs_per_cell());
    get_back_interpolation_matrix(fe1, fe2, interpolation_matrix);

    // compute difference
    difference_matrix = IdentityMatrix(fe1.n_dofs_per_cell());
    difference_matrix.add(-1, interpolation_matrix);
  }



  template <int dim, typename number, int spacedim>
  void
  get_projection_matrix(const FiniteElement<dim, spacedim> &fe1,
                        const FiniteElement<dim, spacedim> &fe2,
                        FullMatrix<number> &                matrix)
  {
    Assert(fe1.n_components() == 1, ExcNotImplemented());
    Assert(fe1.n_components() == fe2.n_components(),
           ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
    Assert(matrix.m() == fe2.n_dofs_per_cell() &&
             matrix.n() == fe1.n_dofs_per_cell(),
           ExcMatrixDimensionMismatch(matrix.m(),
                                      matrix.n(),
                                      fe2.n_dofs_per_cell(),
                                      fe1.n_dofs_per_cell()));
    matrix = 0;

    const unsigned int n1 = fe1.n_dofs_per_cell();
    const unsigned int n2 = fe2.n_dofs_per_cell();

    const ReferenceCell reference_cell = fe1.reference_cell();

    Assert(fe1.reference_cell() == fe2.reference_cell(), ExcNotImplemented());

    // First, create a local mass matrix for the unit cell
    Triangulation<dim, spacedim> tr;
    GridGenerator::reference_cell(tr, reference_cell);

    const auto &mapping =
      reference_cell.template get_default_linear_mapping<dim, spacedim>();

    // Choose a Gauss quadrature rule that is exact up to degree 2n-1
    const unsigned int degree =
      std::max(fe1.tensor_degree(), fe2.tensor_degree());
    Assert(degree != numbers::invalid_unsigned_int, ExcNotImplemented());
    const auto quadrature =
      reference_cell.get_gauss_type_quadrature<dim>(degree + 1);

    // Set up FEValues.
    const UpdateFlags flags =
      update_values | update_quadrature_points | update_JxW_values;
    FEValues<dim> val1(mapping, fe1, quadrature, update_values);
    val1.reinit(tr.begin_active());

    FEValues<dim> val2(mapping, fe2, quadrature, flags);
    val2.reinit(tr.begin_active());

    // Integrate and invert mass matrix. This happens in the target space
    FullMatrix<double> mass(n2, n2);

    for (unsigned int k = 0; k < quadrature.size(); ++k)
      {
        const double dx = val2.JxW(k);
        for (unsigned int i = 0; i < n2; ++i)
          {
            const double v = val2.shape_value(i, k);
            for (unsigned int j = 0; j < n2; ++j)
              mass(i, j) += v * val2.shape_value(j, k) * dx;
          }
      }
    // Invert the matrix. Gauss-Jordan should be sufficient since we expect
    // the mass matrix to be well-conditioned
    mass.gauss_jordan();

    // Now, test every function of fe1 with test functions of fe2 and
    // compute the projection of each unit vector.
    Vector<double> b(n2);
    Vector<double> x(n2);

    for (unsigned int j = 0; j < n1; ++j)
      {
        b = 0.;
        for (unsigned int i = 0; i < n2; ++i)
          for (unsigned int k = 0; k < quadrature.size(); ++k)
            {
              const double dx = val2.JxW(k);
              const double u  = val1.shape_value(j, k);
              const double v  = val2.shape_value(i, k);
              b(i) += u * v * dx;
            }

        // Multiply by the inverse
        mass.vmult(x, b);
        for (unsigned int i = 0; i < n2; ++i)
          matrix(i, j) = x(i);
      }
  }



  template <int dim, int spacedim>
  FullMatrix<double>
  compute_node_matrix(const FiniteElement<dim, spacedim> &fe)
  {
    const unsigned int n_dofs = fe.n_dofs_per_cell();

    FullMatrix<double> N(n_dofs, n_dofs);

    Assert(fe.has_generalized_support_points(), ExcNotInitialized());
    Assert(fe.n_components() == dim, ExcNotImplemented());

    const std::vector<Point<dim>> &points = fe.get_generalized_support_points();

    // We need the values of the polynomials in all generalized support
    // points. This function specifically works for the case where shape
    // functions have 'dim' vector components, so allocate that much space
    std::vector<Vector<double>> support_point_values(points.size(),
                                                     Vector<double>(dim));

    // In this vector, we store the
    // result of the interpolation
    std::vector<double> nodal_values(n_dofs);

    // Get the values of each shape function in turn. Remember that these
    // are the 'raw' shape functions (i.e., where the element has not yet
    // computed the expansion coefficients with regard to the basis
    // provided by the polynomial space).
    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        // get the values of the current set of shape functions
        // at the generalized support points
        for (unsigned int k = 0; k < points.size(); ++k)
          for (unsigned int d = 0; d < dim; ++d)
            {
              support_point_values[k][d] =
                fe.shape_value_component(i, points[k], d);
              Assert(numbers::is_finite(support_point_values[k][d]),
                     ExcInternalError());
            }

        fe.convert_generalized_support_point_values_to_dof_values(
          support_point_values, nodal_values);

        // Enter the interpolated dofs into the matrix
        for (unsigned int j = 0; j < n_dofs; ++j)
          {
            N(j, i) = nodal_values[j];
            Assert(numbers::is_finite(nodal_values[j]), ExcInternalError());
          }
      }

    return N;
  }



  namespace internal
  {
    namespace FEToolsComputeEmbeddingMatricesHelper
    {
      template <int dim, typename number, int spacedim>
      void
      compute_embedding_for_shape_function(
        const unsigned int                  i,
        const FiniteElement<dim, spacedim> &fe,
        const FEValues<dim, spacedim> &     coarse,
        const Householder<double> &         H,
        FullMatrix<number> &                this_matrix,
        const double                        threshold)
      {
        const unsigned int n  = fe.n_dofs_per_cell();
        const unsigned int nd = fe.n_components();
        const unsigned int nq = coarse.n_quadrature_points;

        Vector<number> v_coarse(nq * nd);
        Vector<number> v_fine(n);

        // The right hand side of
        // the least squares
        // problem consists of the
        // function values of the
        // coarse grid function in
        // each quadrature point.
        if (fe.is_primitive())
          {
            const unsigned int d     = fe.system_to_component_index(i).first;
            const double *     phi_i = &coarse.shape_value(i, 0);

            for (unsigned int k = 0; k < nq; ++k)
              v_coarse(k * nd + d) = phi_i[k];
          }

        else
          for (unsigned int d = 0; d < nd; ++d)
            for (unsigned int k = 0; k < nq; ++k)
              v_coarse(k * nd + d) = coarse.shape_value_component(i, k, d);

        // solve the least squares
        // problem.
        const double result = H.least_squares(v_fine, v_coarse);
        Assert(result <= threshold, FETools::ExcLeastSquaresError(result));
        // Avoid warnings in release mode
        (void)result;
        (void)threshold;

        // Copy into the result
        // matrix. Since the matrix
        // maps a coarse grid
        // function to a fine grid
        // function, the columns
        // are fine grid.
        for (unsigned int j = 0; j < n; ++j)
          this_matrix(j, i) = v_fine(j);
      }



      template <int dim, typename number, int spacedim>
      void
      compute_embedding_matrices_for_refinement_case(
        const FiniteElement<dim, spacedim> &fe,
        std::vector<FullMatrix<number>> &   matrices,
        const unsigned int                  ref_case,
        const double                        threshold)
      {
        const unsigned int n = fe.n_dofs_per_cell();
        const unsigned int nc =
          GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

        AssertDimension(matrices.size(), nc);

        for (unsigned int i = 0; i < nc; ++i)
          {
            Assert(matrices[i].n() == n,
                   ExcDimensionMismatch(matrices[i].n(), n));
            Assert(matrices[i].m() == n,
                   ExcDimensionMismatch(matrices[i].m(), n));
          }

        const ReferenceCell reference_cell = fe.reference_cell();

        // Set up meshes, one with a single
        // reference cell and refine it once
        Triangulation<dim, spacedim> tria;
        GridGenerator::reference_cell(tria, reference_cell);
        tria.begin_active()->set_refine_flag(RefinementCase<dim>(ref_case));
        tria.execute_coarsening_and_refinement();

        const unsigned int degree = fe.degree;

        const auto &mapping =
          reference_cell.template get_default_linear_mapping<dim, spacedim>();
        const auto &q_fine =
          reference_cell.get_gauss_type_quadrature<dim>(degree + 1);
        const unsigned int nq = q_fine.size();

        FEValues<dim, spacedim> fine(mapping,
                                     fe,
                                     q_fine,
                                     update_quadrature_points |
                                       update_JxW_values | update_values);

        // We search for the polynomial on
        // the small cell, being equal to
        // the coarse polynomial in all
        // quadrature points.

        // First build the matrix for this
        // least squares problem. This
        // contains the values of the fine
        // cell polynomials in the fine
        // cell grid points.

        // This matrix is the same for all
        // children.
        fine.reinit(tria.begin_active());
        const unsigned int nd = fe.n_components();
        FullMatrix<number> A(nq * nd, n);

        for (unsigned int j = 0; j < n; ++j)
          for (unsigned int d = 0; d < nd; ++d)
            for (unsigned int k = 0; k < nq; ++k)
              A(k * nd + d, j) = fine.shape_value_component(j, k, d);

        Householder<double> H(A);

        Threads::TaskGroup<void> task_group;

        for (const auto &fine_cell : tria.active_cell_iterators())
          {
            fine.reinit(fine_cell);

            // evaluate on the coarse cell (which
            // is the first -- inactive -- cell on
            // the lowest level of the
            // triangulation we have created)
            const std::vector<Point<spacedim>> &q_points_fine =
              fine.get_quadrature_points();
            std::vector<Point<dim>> q_points_coarse(q_points_fine.size());
            for (unsigned int i = 0; i < q_points_fine.size(); ++i)
              for (unsigned int j = 0; j < dim; ++j)
                q_points_coarse[i](j) = q_points_fine[i](j);
            const Quadrature<dim>   q_coarse(q_points_coarse,
                                           fine.get_JxW_values());
            FEValues<dim, spacedim> coarse(mapping,
                                           fe,
                                           q_coarse,
                                           update_values);

            coarse.reinit(tria.begin(0));

            FullMatrix<double> &this_matrix =
              matrices[fine_cell->active_cell_index()];

            // Compute this once for each
            // coarse grid basis function. can
            // spawn subtasks if n is
            // sufficiently large so that there
            // are more than about 5000
            // operations in the inner loop
            // (which is basically const * n^2
            // operations).
            if (n > 30)
              {
                for (unsigned int i = 0; i < n; ++i)
                  {
                    task_group += Threads::new_task(
                      &compute_embedding_for_shape_function<dim,
                                                            number,
                                                            spacedim>,
                      i,
                      fe,
                      coarse,
                      H,
                      this_matrix,
                      threshold);
                  }
                task_group.join_all();
              }
            else
              {
                for (unsigned int i = 0; i < n; ++i)
                  {
                    compute_embedding_for_shape_function<dim, number, spacedim>(
                      i, fe, coarse, H, this_matrix, threshold);
                  }
              }

            // Remove small entries from
            // the matrix
            for (unsigned int i = 0; i < this_matrix.m(); ++i)
              for (unsigned int j = 0; j < this_matrix.n(); ++j)
                if (std::fabs(this_matrix(i, j)) < 1e-12)
                  this_matrix(i, j) = 0.;
          }
      }
    } // namespace FEToolsComputeEmbeddingMatricesHelper
  }   // namespace internal



  template <int dim, typename number, int spacedim>
  void
  compute_embedding_matrices(
    const FiniteElement<dim, spacedim> &          fe,
    std::vector<std::vector<FullMatrix<number>>> &matrices,
    const bool                                    isotropic_only,
    const double                                  threshold)
  {
    Threads::TaskGroup<void> task_group;

    // loop over all possible refinement cases
    unsigned int ref_case = (isotropic_only) ?
                              RefinementCase<dim>::isotropic_refinement :
                              RefinementCase<dim>::cut_x;

    for (; ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
      task_group += Threads::new_task(
        &internal::FEToolsComputeEmbeddingMatricesHelper::
          compute_embedding_matrices_for_refinement_case<dim, number, spacedim>,
        fe,
        matrices[ref_case - 1],
        ref_case,
        threshold);

    task_group.join_all();
  }


  template <int dim, typename number, int spacedim>
  void
  compute_face_embedding_matrices(
    const FiniteElement<dim, spacedim> &fe,
    FullMatrix<number> (&matrices)[GeometryInfo<dim>::max_children_per_face],
    const unsigned int face_coarse,
    const unsigned int face_fine,
    const double       threshold)
  {
    const unsigned int face_no = face_coarse;

    Assert(face_coarse == 0, ExcNotImplemented());
    Assert(face_fine == 0, ExcNotImplemented());

    const unsigned int nc     = GeometryInfo<dim>::max_children_per_face;
    const unsigned int n      = fe.n_dofs_per_face(face_no);
    const unsigned int nd     = fe.n_components();
    const unsigned int degree = fe.degree;

    const bool normal     = fe.conforms(FiniteElementData<dim>::Hdiv);
    const bool tangential = fe.conforms(FiniteElementData<dim>::Hcurl);

    for (unsigned int i = 0; i < nc; ++i)
      {
        Assert(matrices[i].n() == n, ExcDimensionMismatch(matrices[i].n(), n));
        Assert(matrices[i].m() == n, ExcDimensionMismatch(matrices[i].m(), n));
      }

    // In order to make the loops below
    // simpler, we introduce vectors
    // containing for indices 0-n the
    // number of the corresponding
    // shape value on the cell.
    std::vector<unsigned int> face_c_dofs(n);
    std::vector<unsigned int> face_f_dofs(n);
    {
      unsigned int face_dof = 0;
      for (unsigned int i = 0;
           i < fe.reference_cell().face_reference_cell(face_no).n_vertices();
           ++i)
        {
          const unsigned int offset_c =
            GeometryInfo<dim>::face_to_cell_vertices(face_coarse, i) *
            fe.n_dofs_per_vertex();
          const unsigned int offset_f =
            GeometryInfo<dim>::face_to_cell_vertices(face_fine, i) *
            fe.n_dofs_per_vertex();
          for (unsigned int j = 0; j < fe.n_dofs_per_vertex(); ++j)
            {
              face_c_dofs[face_dof] = offset_c + j;
              face_f_dofs[face_dof] = offset_f + j;
              ++face_dof;
            }
        }

      for (unsigned int i = 1;
           i <= fe.reference_cell().face_reference_cell(face_no).n_lines();
           ++i)
        {
          const unsigned int offset_c =
            fe.get_first_line_index() +
            GeometryInfo<dim>::face_to_cell_lines(face_coarse, i - 1) *
              fe.n_dofs_per_line();
          const unsigned int offset_f =
            fe.get_first_line_index() +
            GeometryInfo<dim>::face_to_cell_lines(face_fine, i - 1) *
              fe.n_dofs_per_line();
          for (unsigned int j = 0; j < fe.n_dofs_per_line(); ++j)
            {
              face_c_dofs[face_dof] = offset_c + j;
              face_f_dofs[face_dof] = offset_f + j;
              ++face_dof;
            }
        }

      if (dim == 3)
        {
          const unsigned int offset_c = fe.get_first_quad_index(face_coarse);
          const unsigned int offset_f = fe.get_first_quad_index(face_fine);
          for (unsigned int j = 0; j < fe.n_dofs_per_quad(face_no); ++j)
            {
              face_c_dofs[face_dof] = offset_c + j;
              face_f_dofs[face_dof] = offset_f + j;
              ++face_dof;
            }
        }
      Assert(face_dof == fe.n_dofs_per_face(face_no), ExcInternalError());
    }

    // Set up meshes, one with a single
    // reference cell and refine it once
    Triangulation<dim, spacedim> tria;
    GridGenerator::hyper_cube(tria, 0, 1);
    tria.refine_global(1);
    MappingCartesian<dim> mapping;

    // Setup quadrature and FEValues
    // for a face. We cannot use
    // FEFaceValues and
    // FESubfaceValues because of
    // some nifty handling of
    // refinement cases. Guido stops
    // disliking and instead starts
    // hating the anisotropic implementation
    QGauss<dim - 1>       q_gauss(degree + 1);
    const Quadrature<dim> q_fine =
      QProjector<dim>::project_to_face(fe.reference_cell(), q_gauss, face_fine);
    const unsigned int nq = q_fine.size();

    FEValues<dim> fine(mapping,
                       fe,
                       q_fine,
                       update_quadrature_points | update_JxW_values |
                         update_values);

    // We search for the polynomial on
    // the small cell, being equal to
    // the coarse polynomial in all
    // quadrature points.

    // First build the matrix for this
    // least squares problem. This
    // contains the values of the fine
    // cell polynomials in the fine
    // cell grid points.

    // This matrix is the same for all
    // children.
    fine.reinit(tria.begin_active());
    FullMatrix<number> A(nq * nd, n);
    for (unsigned int j = 0; j < n; ++j)
      for (unsigned int k = 0; k < nq; ++k)
        if (nd != dim)
          for (unsigned int d = 0; d < nd; ++d)
            A(k * nd + d, j) = fine.shape_value_component(face_f_dofs[j], k, d);
        else
          {
            if (normal)
              A(k * nd, j) = fine.shape_value_component(face_f_dofs[j], k, 0);
            if (tangential)
              for (unsigned int d = 1; d < dim; ++d)
                A(k * nd + d, j) =
                  fine.shape_value_component(face_f_dofs[j], k, d);
          }

    Householder<double> H(A);

    Vector<number> v_coarse(nq * nd);
    Vector<number> v_fine(n);


    for (unsigned int cell_number = 0;
         cell_number < GeometryInfo<dim>::max_children_per_face;
         ++cell_number)
      {
        const Quadrature<dim> q_coarse = QProjector<dim>::project_to_subface(
          fe.reference_cell(), q_gauss, face_coarse, cell_number);
        FEValues<dim> coarse(mapping, fe, q_coarse, update_values);

        typename Triangulation<dim, spacedim>::active_cell_iterator fine_cell =
          tria.begin(0)->child(GeometryInfo<dim>::child_cell_on_face(
            tria.begin(0)->refinement_case(), face_coarse, cell_number));
        fine.reinit(fine_cell);
        coarse.reinit(tria.begin(0));

        FullMatrix<double> &this_matrix = matrices[cell_number];

        // Compute this once for each
        // coarse grid basis function
        for (unsigned int i = 0; i < n; ++i)
          {
            // The right hand side of
            // the least squares
            // problem consists of the
            // function values of the
            // coarse grid function in
            // each quadrature point.
            for (unsigned int k = 0; k < nq; ++k)
              if (nd != dim)
                for (unsigned int d = 0; d < nd; ++d)
                  v_coarse(k * nd + d) =
                    coarse.shape_value_component(face_c_dofs[i], k, d);
              else
                {
                  if (normal)
                    v_coarse(k * nd) =
                      coarse.shape_value_component(face_c_dofs[i], k, 0);
                  if (tangential)
                    for (unsigned int d = 1; d < dim; ++d)
                      v_coarse(k * nd + d) =
                        coarse.shape_value_component(face_c_dofs[i], k, d);
                }
            // solve the least squares
            // problem.
            const double result = H.least_squares(v_fine, v_coarse);
            Assert(result <= threshold, FETools::ExcLeastSquaresError(result));
            // Avoid compiler warnings in Release mode
            (void)result;
            (void)threshold;

            // Copy into the result
            // matrix. Since the matrix
            // maps a coarse grid
            // function to a fine grid
            // function, the columns
            // are fine grid.
            for (unsigned int j = 0; j < n; ++j)
              this_matrix(j, i) = v_fine(j);
          }
        // Remove small entries from
        // the matrix
        for (unsigned int i = 0; i < this_matrix.m(); ++i)
          for (unsigned int j = 0; j < this_matrix.n(); ++j)
            if (std::fabs(this_matrix(i, j)) < 1e-12)
              this_matrix(i, j) = 0.;
      }
  }


  template <int dim, typename number, int spacedim>
  void
  compute_projection_matrices(
    const FiniteElement<dim, spacedim> &          fe,
    std::vector<std::vector<FullMatrix<number>>> &matrices,
    const bool                                    isotropic_only)
  {
    const unsigned int n      = fe.n_dofs_per_cell();
    const unsigned int nd     = fe.n_components();
    const unsigned int degree = fe.degree;

    const ReferenceCell reference_cell = fe.reference_cell();

    const auto &mapping =
      reference_cell.template get_default_linear_mapping<dim, spacedim>();
    const auto &q_fine =
      reference_cell.get_gauss_type_quadrature<dim>(degree + 1);

    // prepare FEValues, quadrature etc on
    // coarse cell
    const unsigned int nq = q_fine.size();

    // create mass matrix on coarse cell.
    FullMatrix<number> mass(n, n);
    {
      // set up a triangulation for coarse cell
      Triangulation<dim, spacedim> tr;
      GridGenerator::reference_cell(tr, reference_cell);

      FEValues<dim, spacedim> coarse(mapping,
                                     fe,
                                     q_fine,
                                     update_JxW_values | update_values);

      typename Triangulation<dim, spacedim>::cell_iterator coarse_cell =
        tr.begin(0);
      coarse.reinit(coarse_cell);

      const std::vector<double> &JxW = coarse.get_JxW_values();
      for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = 0; j < n; ++j)
          if (fe.is_primitive())
            {
              const double *coarse_i = &coarse.shape_value(i, 0);
              const double *coarse_j = &coarse.shape_value(j, 0);
              double        mass_ij  = 0;
              for (unsigned int k = 0; k < nq; ++k)
                mass_ij += JxW[k] * coarse_i[k] * coarse_j[k];
              mass(i, j) = mass_ij;
            }
          else
            {
              double mass_ij = 0;
              for (unsigned int d = 0; d < nd; ++d)
                for (unsigned int k = 0; k < nq; ++k)
                  mass_ij += JxW[k] * coarse.shape_value_component(i, k, d) *
                             coarse.shape_value_component(j, k, d);
              mass(i, j) = mass_ij;
            }

      // invert mass matrix
      mass.gauss_jordan();
    }


    const auto compute_one_case =
      [&reference_cell, &mapping, &fe, &q_fine, n, nd, nq](
        const unsigned int               ref_case,
        const FullMatrix<double> &       inverse_mass_matrix,
        std::vector<FullMatrix<double>> &matrices) {
        const unsigned int nc =
          GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

        for (unsigned int i = 0; i < nc; ++i)
          {
            Assert(matrices[i].n() == n,
                   ExcDimensionMismatch(matrices[i].n(), n));
            Assert(matrices[i].m() == n,
                   ExcDimensionMismatch(matrices[i].m(), n));
          }

        // create a respective refinement on the triangulation
        Triangulation<dim, spacedim> tr;
        GridGenerator::reference_cell(tr, reference_cell);
        tr.begin_active()->set_refine_flag(RefinementCase<dim>(ref_case));
        tr.execute_coarsening_and_refinement();

        FEValues<dim, spacedim> fine(mapping,
                                     fe,
                                     q_fine,
                                     update_quadrature_points |
                                       update_JxW_values | update_values);

        typename Triangulation<dim, spacedim>::cell_iterator coarse_cell =
          tr.begin(0);

        Vector<number> v_coarse(n);
        Vector<number> v_fine(n);

        for (unsigned int cell_number = 0; cell_number < nc; ++cell_number)
          {
            FullMatrix<double> &this_matrix = matrices[cell_number];

            // Compute right hand side, which is a fine level basis
            // function tested with the coarse level functions.
            fine.reinit(coarse_cell->child(cell_number));
            const std::vector<Point<spacedim>> &q_points_fine =
              fine.get_quadrature_points();
            std::vector<Point<dim>> q_points_coarse(q_points_fine.size());
            for (unsigned int q = 0; q < q_points_fine.size(); ++q)
              for (unsigned int j = 0; j < dim; ++j)
                q_points_coarse[q](j) = q_points_fine[q](j);
            Quadrature<dim> q_coarse(q_points_coarse, fine.get_JxW_values());
            FEValues<dim, spacedim> coarse(
              coarse_cell->reference_cell()
                .template get_default_linear_mapping<dim, spacedim>(),
              fe,
              q_coarse,
              update_values);
            coarse.reinit(coarse_cell);

            // Build RHS

            const std::vector<double> &JxW = fine.get_JxW_values();

            // Outer loop over all fine grid shape functions phi_j
            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
              {
                for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
                  {
                    if (fe.is_primitive())
                      {
                        const double *coarse_i = &coarse.shape_value(i, 0);
                        const double *fine_j   = &fine.shape_value(j, 0);

                        double update = 0;
                        for (unsigned int k = 0; k < nq; ++k)
                          update += JxW[k] * coarse_i[k] * fine_j[k];
                        v_fine(i) = update;
                      }
                    else
                      {
                        double update = 0;
                        for (unsigned int d = 0; d < nd; ++d)
                          for (unsigned int k = 0; k < nq; ++k)
                            update += JxW[k] *
                                      coarse.shape_value_component(i, k, d) *
                                      fine.shape_value_component(j, k, d);
                        v_fine(i) = update;
                      }
                  }

                // RHS ready. Solve system and enter row into matrix
                inverse_mass_matrix.vmult(v_coarse, v_fine);
                for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
                  this_matrix(i, j) = v_coarse(i);
              }

            // Remove small entries from the matrix
            for (unsigned int i = 0; i < this_matrix.m(); ++i)
              for (unsigned int j = 0; j < this_matrix.n(); ++j)
                if (std::fabs(this_matrix(i, j)) < 1e-12)
                  this_matrix(i, j) = 0.;
          }
      };


    // finally loop over all possible refinement cases
    Threads::TaskGroup<> tasks;
    unsigned int         ref_case = (isotropic_only) ?
                              RefinementCase<dim>::isotropic_refinement :
                              RefinementCase<dim>::cut_x;
    for (; ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
      tasks += Threads::new_task([&, ref_case]() {
        compute_one_case(ref_case, mass, matrices[ref_case - 1]);
      });

    tasks.join_all();
  }


  template <int dim, int spacedim>
  void
  add_fe_name(const std::string &                 parameter_name,
              const FEFactoryBase<dim, spacedim> *factory)
  {
    // Erase everything after the
    // actual class name
    std::string  name     = parameter_name;
    unsigned int name_end = name.find_first_not_of(std::string(
      "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
    if (name_end < name.size())
      name.erase(name_end);
    // first make sure that no other
    // thread intercepts the
    // operation of this function;
    // for this, acquire the lock
    // until we quit this function
    std::lock_guard<std::mutex> lock(
      internal::FEToolsAddFENameHelper::fe_name_map_lock);

    Assert(
      internal::FEToolsAddFENameHelper::get_fe_name_map()[dim][spacedim].find(
        name) ==
        internal::FEToolsAddFENameHelper::get_fe_name_map()[dim][spacedim]
          .end(),
      ExcMessage("Cannot change existing element in finite element name list"));

    // Insert the normalized name into
    // the map
    internal::FEToolsAddFENameHelper::get_fe_name_map()[dim][spacedim][name] =
      std::unique_ptr<const Subscriptor>(factory);
  }



  namespace internal
  {
    namespace FEToolsGetFEHelper
    {
      // TODO: this encapsulates the call to the
      // dimension-dependent fe_name_map so that we
      // have a unique interface. could be done
      // smarter?
      template <int dim, int spacedim>
      std::unique_ptr<FiniteElement<dim, spacedim>>
      get_fe_by_name_ext(
        std::string &name,
        const std::map<std::string, std::unique_ptr<const Subscriptor>>
          &fe_name_map)
      {
        // Extract the name of the
        // finite element class, which only
        // contains characters, numbers and
        // underscores.
        unsigned int      name_end = name.find_first_not_of(std::string(
          "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
        const std::string name_part(name, 0, name_end);
        name.erase(0, name_part.size());

        // now things get a little more
        // complicated: FESystem. it's
        // more complicated, since we
        // have to figure out what the
        // base elements are. this can
        // only be done recursively
        if (name_part == "FESystem")
          {
            // next we have to get at the
            // base elements. start with
            // the first. wrap the whole
            // block into try-catch to
            // make sure we destroy the
            // pointers we got from
            // recursive calls if one of
            // these calls should throw
            // an exception
            std::vector<std::unique_ptr<const FiniteElement<dim, spacedim>>>
                                      base_fes;
            std::vector<unsigned int> base_multiplicities;

            // Now, just the [...]
            // part should be left.
            if (name.size() == 0 || name[0] != '[')
              throw std::string("Invalid first character in ") + name;
            do
              {
                // Erase the
                // leading '[' or '-'
                name.erase(0, 1);
                // Now, the name of the
                // first base element is
                // first... Let's get it
                base_fes.push_back(
                  get_fe_by_name_ext<dim, spacedim>(name, fe_name_map));
                // next check whether
                // FESystem placed a
                // multiplicity after
                // the element name
                if (name[0] == '^')
                  {
                    // yes. Delete the '^'
                    // and read this
                    // multiplicity
                    name.erase(0, 1);

                    const std::pair<int, unsigned int> tmp =
                      Utilities::get_integer_at_position(name, 0);
                    name.erase(0, tmp.second);
                    // add to length,
                    // including the '^'
                    base_multiplicities.push_back(tmp.first);
                  }
                else
                  // no, so
                  // multiplicity is
                  // 1
                  base_multiplicities.push_back(1);

                // so that's it for
                // this base
                // element. base
                // elements are
                // separated by '-',
                // and the list is
                // terminated by ']',
                // so loop while the
                // next character is
                // '-'
              }
            while (name[0] == '-');

            // so we got to the end
            // of the '-' separated
            // list. make sure that
            // we actually had a ']'
            // there
            if (name.size() == 0 || name[0] != ']')
              throw std::string("Invalid first character in ") + name;
            name.erase(0, 1);
            // just one more sanity check
            Assert((base_fes.size() == base_multiplicities.size()) &&
                     (base_fes.size() > 0),
                   ExcInternalError());

            // this is a workaround since the constructor for FESystem
            // only accepts raw pointers.
            std::vector<const FiniteElement<dim, spacedim> *> raw_base_fes;
            raw_base_fes.reserve(base_fes.size());
            for (const auto &base_fe : base_fes)
              raw_base_fes.push_back(base_fe.get());

            // ok, apparently everything went ok. so generate the composed
            // element. and return it.
            return std::make_unique<FESystem<dim, spacedim>>(
              raw_base_fes, base_multiplicities);
          }
        else if (name_part == "FE_Nothing")
          {
            // remove the () from FE_Nothing()
            name.erase(0, 2);

            // this is a bit of a hack, as
            // FE_Nothing does not take a
            // degree, but it does take an
            // argument, which defaults to 1,
            // so this properly returns
            // FE_Nothing()
            const Subscriptor *ptr = fe_name_map.find(name_part)->second.get();
            const FETools::FEFactoryBase<dim, spacedim> *fef =
              dynamic_cast<const FETools::FEFactoryBase<dim, spacedim> *>(ptr);
            return fef->get(1);
          }
        else
          {
            // Make sure no other thread
            // is just adding an element
            std::lock_guard<std::mutex> lock(
              internal::FEToolsAddFENameHelper::fe_name_map_lock);
            AssertThrow(fe_name_map.find(name_part) != fe_name_map.end(),
                        FETools::ExcInvalidFEName(name));

            // Now, just the (degree)
            // or (Quadrature<1>(degree+1))
            // part should be left.
            if (name.size() == 0 || name[0] != '(')
              throw std::string("Invalid first character in ") + name;
            name.erase(0, 1);
            if (name[0] != 'Q')
              {
                const std::pair<int, unsigned int> tmp =
                  Utilities::get_integer_at_position(name, 0);
                name.erase(0, tmp.second + 1);
                const Subscriptor *ptr =
                  fe_name_map.find(name_part)->second.get();
                const FETools::FEFactoryBase<dim, spacedim> *fef =
                  dynamic_cast<const FETools::FEFactoryBase<dim, spacedim> *>(
                    ptr);
                return fef->get(tmp.first);
              }
            else
              {
                unsigned int      position = name.find('(');
                const std::string quadrature_name(name, 0, position);
                name.erase(0, position + 1);
                if (quadrature_name.compare("QGaussLobatto") == 0)
                  {
                    const std::pair<int, unsigned int> tmp =
                      Utilities::get_integer_at_position(name, 0);
                    // delete "))"
                    name.erase(0, tmp.second + 2);
                    const Subscriptor *ptr =
                      fe_name_map.find(name_part)->second.get();
                    const FETools::FEFactoryBase<dim, spacedim> *fef =
                      dynamic_cast<
                        const FETools::FEFactoryBase<dim, spacedim> *>(ptr);
                    return fef->get(QGaussLobatto<1>(tmp.first));
                  }
                else if (quadrature_name.compare("QGauss") == 0)
                  {
                    const std::pair<int, unsigned int> tmp =
                      Utilities::get_integer_at_position(name, 0);
                    // delete "))"
                    name.erase(0, tmp.second + 2);
                    const Subscriptor *ptr =
                      fe_name_map.find(name_part)->second.get();
                    const FETools::FEFactoryBase<dim, spacedim> *fef =
                      dynamic_cast<
                        const FETools::FEFactoryBase<dim, spacedim> *>(ptr);
                    return fef->get(QGauss<1>(tmp.first));
                  }
                else if (quadrature_name.compare("QIterated") == 0)
                  {
                    // find sub-quadrature
                    position = name.find('(');
                    const std::string subquadrature_name(name, 0, position);
                    AssertThrow(subquadrature_name == "QTrapez" ||
                                  subquadrature_name == "QTrapezoid",
                                ExcNotImplemented(
                                  "Could not detect quadrature of name " +
                                  subquadrature_name));
                    // delete "QTrapezoid(),"
                    name.erase(0, position + 3);
                    const std::pair<int, unsigned int> tmp =
                      Utilities::get_integer_at_position(name, 0);
                    // delete "))"
                    name.erase(0, tmp.second + 2);
                    const Subscriptor *ptr =
                      fe_name_map.find(name_part)->second.get();
                    const FETools::FEFactoryBase<dim, spacedim> *fef =
                      dynamic_cast<
                        const FETools::FEFactoryBase<dim, spacedim> *>(ptr);
                    return fef->get(QIterated<1>(QTrapezoid<1>(), tmp.first));
                  }
                else
                  {
                    AssertThrow(false, ExcNotImplemented());
                  }
              }
          }


        // hm, if we have come thus far, we
        // didn't know what to do with the
        // string we got. so do as the docs
        // say: raise an exception
        AssertThrow(false, FETools::ExcInvalidFEName(name));

        // make some compilers happy that
        // do not realize that we can't get
        // here after throwing
        return nullptr;
      }



      template <int dim, int spacedim>
      std::unique_ptr<FiniteElement<dim, spacedim>>
      get_fe_by_name(std::string &name)
      {
        return get_fe_by_name_ext<dim, spacedim>(
          name, FEToolsAddFENameHelper::get_fe_name_map()[dim][spacedim]);
      }
    } // namespace FEToolsGetFEHelper
  }   // namespace internal



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  get_fe_by_name(const std::string &parameter_name)
  {
    std::string name  = Utilities::trim(parameter_name);
    std::size_t index = 1;
    // remove spaces that are not between two word (things that match the
    // regular expression [A-Za-z0-9_]) characters.
    while (2 < name.size() && index < name.size() - 1)
      {
        if (name[index] == ' ' &&
            (!(std::isalnum(name[index - 1]) || name[index - 1] == '_') ||
             !(std::isalnum(name[index + 1]) || name[index + 1] == '_')))
          {
            name.erase(index, 1);
          }
        else
          {
            ++index;
          }
      }

    // Create a version of the name
    // string where all template
    // parameters are eliminated.
    for (unsigned int pos1 = name.find('<'); pos1 < name.size();
         pos1              = name.find('<'))
      {
        const unsigned int pos2 = name.find('>');
        // If there is only a single
        // character between those two,
        // it should be 'd' or the number
        // representing the dimension.
        if (pos2 - pos1 == 2)
          {
            const char dimchar = '0' + dim;
            (void)dimchar;
            if (name.at(pos1 + 1) != 'd')
              Assert(name.at(pos1 + 1) == dimchar,
                     ExcInvalidFEDimension(name.at(pos1 + 1), dim));
          }
        else
          Assert(pos2 - pos1 == 4, ExcInvalidFEName(name));

        // If pos1==pos2, then we are
        // probably at the end of the
        // string
        if (pos2 != pos1)
          name.erase(pos1, pos2 - pos1 + 1);
      }
    // Replace all occurrences of "^dim"
    // by "^d" to be handled by the
    // next loop
    for (unsigned int pos = name.find("^dim"); pos < name.size();
         pos              = name.find("^dim"))
      name.erase(pos + 2, 2);

    // Replace all occurrences of "^d"
    // by using the actual dimension
    for (unsigned int pos = name.find("^d"); pos < name.size();
         pos              = name.find("^d"))
      name.at(pos + 1) = '0' + dim;

    try
      {
        auto fe =
          internal::FEToolsGetFEHelper::get_fe_by_name<dim, spacedim>(name);

        // Make sure the auxiliary function
        // ate up all characters of the name.
        AssertThrow(name.size() == 0,
                    ExcInvalidFEName(parameter_name +
                                     std::string(" extra characters after "
                                                 "end of name")));
        return fe;
      }
    catch (const std::string &errline)
      {
        AssertThrow(false,
                    ExcInvalidFEName(parameter_name + std::string(" at ") +
                                     errline));
        return nullptr;
      }
  }



  template <int dim, int spacedim>
  void
  compute_projection_from_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             lhs_quadrature,
    const Quadrature<dim> &             rhs_quadrature,
    FullMatrix<double> &                X)
  {
    Assert(fe.n_components() == 1, ExcNotImplemented());

    // first build the matrices M and Q
    // described in the documentation
    FullMatrix<double> M(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
    FullMatrix<double> Q(fe.n_dofs_per_cell(), rhs_quadrature.size());

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
        for (unsigned int q = 0; q < lhs_quadrature.size(); ++q)
          M(i, j) += fe.shape_value(i, lhs_quadrature.point(q)) *
                     fe.shape_value(j, lhs_quadrature.point(q)) *
                     lhs_quadrature.weight(q);

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      for (unsigned int q = 0; q < rhs_quadrature.size(); ++q)
        Q(i, q) +=
          fe.shape_value(i, rhs_quadrature.point(q)) * rhs_quadrature.weight(q);

    // then invert M
    FullMatrix<double> M_inverse(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
    M_inverse.invert(M);

    // finally compute the result
    X.reinit(fe.n_dofs_per_cell(), rhs_quadrature.size());
    M_inverse.mmult(X, Q);

    Assert(X.m() == fe.n_dofs_per_cell(), ExcInternalError());
    Assert(X.n() == rhs_quadrature.size(), ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  compute_interpolation_to_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    FullMatrix<double> &                I_q)
  {
    Assert(fe.n_components() == 1, ExcNotImplemented());
    Assert(I_q.m() == quadrature.size(), ExcMessage("Wrong matrix size"));
    Assert(I_q.n() == fe.n_dofs_per_cell(), ExcMessage("Wrong matrix size"));

    for (unsigned int q = 0; q < quadrature.size(); ++q)
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        I_q(q, i) = fe.shape_value(i, quadrature.point(q));
  }



  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double> &         projection_matrix,
    const std::vector<Tensor<1, dim>> &vector_of_tensors_at_qp,
    std::vector<Tensor<1, dim>> &      vector_of_tensors_at_nodes)
  {
    // check that the number columns of the projection_matrix
    // matches the size of the vector_of_tensors_at_qp
    Assert(projection_matrix.n_cols() == vector_of_tensors_at_qp.size(),
           ExcDimensionMismatch(projection_matrix.n_cols(),
                                vector_of_tensors_at_qp.size()));

    // check that the number rows of the projection_matrix
    // matches the size of the vector_of_tensors_at_nodes
    Assert(projection_matrix.n_rows() == vector_of_tensors_at_nodes.size(),
           ExcDimensionMismatch(projection_matrix.n_rows(),
                                vector_of_tensors_at_nodes.size()));

    // number of support points (nodes) to project to
    const unsigned int n_support_points = projection_matrix.n_rows();
    // number of quadrature points to project from
    const unsigned int n_quad_points = projection_matrix.n_cols();

    // component projected to the nodes
    Vector<double> component_at_node(n_support_points);
    // component at the quadrature point
    Vector<double> component_at_qp(n_quad_points);

    for (unsigned int ii = 0; ii < dim; ++ii)
      {
        component_at_qp = 0;

        // populate the vector of components at the qps
        // from vector_of_tensors_at_qp
        // vector_of_tensors_at_qp data is in form:
        //      columns:        0, 1, ...,  dim
        //      rows:           0,1,....,  n_quad_points
        // so extract the ii'th column of vector_of_tensors_at_qp
        for (unsigned int q = 0; q < n_quad_points; ++q)
          {
            component_at_qp(q) = vector_of_tensors_at_qp[q][ii];
          }

        // project from the qps -> nodes
        // component_at_node = projection_matrix_u * component_at_qp
        projection_matrix.vmult(component_at_node, component_at_qp);

        // rewrite the projection of the components
        // back into the vector of tensors
        for (unsigned int nn = 0; nn < n_support_points; ++nn)
          {
            vector_of_tensors_at_nodes[nn][ii] = component_at_node(nn);
          }
      }
  }



  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double> &                  projection_matrix,
    const std::vector<SymmetricTensor<2, dim>> &vector_of_tensors_at_qp,
    std::vector<SymmetricTensor<2, dim>> &      vector_of_tensors_at_nodes)
  {
    // check that the number columns of the projection_matrix
    // matches the size of the vector_of_tensors_at_qp
    Assert(projection_matrix.n_cols() == vector_of_tensors_at_qp.size(),
           ExcDimensionMismatch(projection_matrix.n_cols(),
                                vector_of_tensors_at_qp.size()));

    // check that the number rows of the projection_matrix
    // matches the size of the vector_of_tensors_at_nodes
    Assert(projection_matrix.n_rows() == vector_of_tensors_at_nodes.size(),
           ExcDimensionMismatch(projection_matrix.n_rows(),
                                vector_of_tensors_at_nodes.size()));

    // number of support points (nodes)
    const unsigned int n_support_points = projection_matrix.n_rows();
    // number of quadrature points to project from
    const unsigned int n_quad_points = projection_matrix.n_cols();

    // number of unique entries in a symmetric second-order tensor
    const unsigned int n_independent_components =
      SymmetricTensor<2, dim>::n_independent_components;

    // component projected to the nodes
    Vector<double> component_at_node(n_support_points);
    // component at the quadrature point
    Vector<double> component_at_qp(n_quad_points);

    // loop over the number of unique dimensions of the tensor
    for (unsigned int ii = 0; ii < n_independent_components; ++ii)
      {
        component_at_qp = 0;

        // row-column entry of tensor corresponding the unrolled index
        TableIndices<2> row_column_index =
          SymmetricTensor<2, dim>::unrolled_to_component_indices(ii);
        const unsigned int row    = row_column_index[0];
        const unsigned int column = row_column_index[1];

        //  populate the vector of components at the qps
        //  from vector_of_tensors_at_qp
        //  vector_of_tensors_at_qp is in form:
        //      columns:       0, 1, ..., n_independent_components
        //      rows:           0,1,....,  n_quad_points
        //  so extract the ii'th column of vector_of_tensors_at_qp
        for (unsigned int q = 0; q < n_quad_points; ++q)
          {
            component_at_qp(q) = (vector_of_tensors_at_qp[q])[row][column];
          }

        // project from the qps -> nodes
        // component_at_node = projection_matrix_u * component_at_qp
        projection_matrix.vmult(component_at_node, component_at_qp);

        // rewrite the projection of the components back into the vector of
        // tensors
        for (unsigned int nn = 0; nn < n_support_points; ++nn)
          {
            (vector_of_tensors_at_nodes[nn])[row][column] =
              component_at_node(nn);
          }
      }
  }



  template <int dim, int spacedim>
  void
  compute_projection_from_face_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim - 1> &         lhs_quadrature,
    const Quadrature<dim - 1> &         rhs_quadrature,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face,
    FullMatrix<double> &                                            X)
  {
    Assert(fe.n_components() == 1, ExcNotImplemented());
    Assert(lhs_quadrature.size() > fe.degree,
           ExcNotGreaterThan(lhs_quadrature.size(), fe.degree));



    // build the matrices M and Q
    // described in the documentation
    FullMatrix<double> M(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
    FullMatrix<double> Q(fe.n_dofs_per_cell(), rhs_quadrature.size());

    {
      // need an FEFaceValues object to evaluate shape function
      // values on the specified face.
      FEFaceValues<dim> fe_face_values(fe, lhs_quadrature, update_values);
      fe_face_values.reinit(cell, face); // setup shape_value on this face.

      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
          for (unsigned int q = 0; q < lhs_quadrature.size(); ++q)
            M(i, j) += fe_face_values.shape_value(i, q) *
                       fe_face_values.shape_value(j, q) *
                       lhs_quadrature.weight(q);
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        {
          M(i, i) = (M(i, i) == 0 ? 1 : M(i, i));
        }
    }

    {
      FEFaceValues<dim> fe_face_values(fe, rhs_quadrature, update_values);
      fe_face_values.reinit(cell, face); // setup shape_value on this face.

      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        for (unsigned int q = 0; q < rhs_quadrature.size(); ++q)
          Q(i, q) +=
            fe_face_values.shape_value(i, q) * rhs_quadrature.weight(q);
    }
    // then invert M
    FullMatrix<double> M_inverse(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
    M_inverse.invert(M);

    // finally compute the result
    X.reinit(fe.n_dofs_per_cell(), rhs_quadrature.size());
    M_inverse.mmult(X, Q);

    Assert(X.m() == fe.n_dofs_per_cell(), ExcInternalError());
    Assert(X.n() == rhs_quadrature.size(), ExcInternalError());
  }



  namespace internal
  {
    namespace FEToolsConvertHelper
    {
      // Helper functions for
      // FETools::convert_generalized_support_point_values_to_dof_values

      template <int dim, int spacedim, typename number>
      static void
      convert_helper(const FiniteElement<dim, spacedim> &finite_element,
                     const std::vector<Vector<number>> & support_point_values,
                     std::vector<number> &               dof_values)
      {
        static Threads::ThreadLocalStorage<std::vector<Vector<double>>>
          double_support_point_values;
        static Threads::ThreadLocalStorage<std::vector<double>>
          double_dof_values;

        double_support_point_values.get().resize(support_point_values.size());
        double_dof_values.get().resize(dof_values.size());

        for (unsigned int i = 0; i < support_point_values.size(); ++i)
          {
            double_support_point_values.get()[i].reinit(
              finite_element.n_components(), false);
            std::copy(std::begin(support_point_values[i]),
                      std::end(support_point_values[i]),
                      std::begin(double_support_point_values.get()[i]));
          }

        finite_element.convert_generalized_support_point_values_to_dof_values(
          double_support_point_values.get(), double_dof_values.get());

        std::copy(std::begin(double_dof_values.get()),
                  std::end(double_dof_values.get()),
                  std::begin(dof_values));
      }


      template <int dim, int spacedim, typename number>
      static void
      convert_helper(
        const FiniteElement<dim, spacedim> &             finite_element,
        const std::vector<Vector<std::complex<number>>> &support_point_values,
        std::vector<std::complex<number>> &              dof_values)
      {
        static Threads::ThreadLocalStorage<std::vector<Vector<double>>>
          double_support_point_values_real;
        static Threads::ThreadLocalStorage<std::vector<double>>
          double_dof_values_real;
        static Threads::ThreadLocalStorage<std::vector<Vector<double>>>
          double_support_point_values_imag;
        static Threads::ThreadLocalStorage<std::vector<double>>
          double_dof_values_imag;

        double_support_point_values_real.get().resize(
          support_point_values.size());
        double_dof_values_real.get().resize(dof_values.size());
        double_support_point_values_imag.get().resize(
          support_point_values.size());
        double_dof_values_imag.get().resize(dof_values.size());

        for (unsigned int i = 0; i < support_point_values.size(); ++i)
          {
            double_support_point_values_real.get()[i].reinit(
              finite_element.n_components(), false);
            double_support_point_values_imag.get()[i].reinit(
              finite_element.n_components(), false);

            std::transform(
              std::begin(support_point_values[i]),
              std::end(support_point_values[i]),
              std::begin(double_support_point_values_real.get()[i]),
              [](std::complex<number> c) -> double { return c.real(); });

            std::transform(
              std::begin(support_point_values[i]),
              std::end(support_point_values[i]),
              std::begin(double_support_point_values_imag.get()[i]),
              [](std::complex<number> c) -> double { return c.imag(); });
          }

        finite_element.convert_generalized_support_point_values_to_dof_values(
          double_support_point_values_real.get(), double_dof_values_real.get());
        finite_element.convert_generalized_support_point_values_to_dof_values(
          double_support_point_values_imag.get(), double_dof_values_imag.get());

        std::transform(std::begin(double_dof_values_real.get()),
                       std::end(double_dof_values_real.get()),
                       std::begin(double_dof_values_imag.get()),
                       std::begin(dof_values),
                       [](number real, number imag) -> std::complex<number> {
                         return {real, imag};
                       });
      }


      template <int dim, int spacedim>
      static void
      convert_helper(const FiniteElement<dim, spacedim> &finite_element,
                     const std::vector<Vector<double>> & support_point_values,
                     std::vector<double> &               dof_values)
      {
        finite_element.convert_generalized_support_point_values_to_dof_values(
          support_point_values, dof_values);
      }

    } // namespace FEToolsConvertHelper
  }   // namespace internal



  template <int dim, int spacedim, typename number>
  void
  convert_generalized_support_point_values_to_dof_values(
    const FiniteElement<dim, spacedim> &finite_element,
    const std::vector<Vector<number>> & support_point_values,
    std::vector<number> &               dof_values)
  {
    AssertDimension(support_point_values.size(),
                    finite_element.get_generalized_support_points().size());
    AssertDimension(dof_values.size(), finite_element.n_dofs_per_cell());

    internal::FEToolsConvertHelper::convert_helper<dim, spacedim>(
      finite_element, support_point_values, dof_values);
  }



  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering(const unsigned int degree)
  {
    // number of support points in each direction
    const unsigned int n = degree + 1;

    const unsigned int dofs_per_cell = Utilities::fixed_power<dim>(n);

    std::vector<unsigned int> h2l(dofs_per_cell);

    // polynomial degree
    const unsigned int dofs_per_line = degree - 1;

    // the following lines of code are somewhat odd, due to the way the
    // hierarchic numbering is organized. if someone would really want to
    // understand these lines, you better draw some pictures where you
    // indicate the indices and orders of vertices, lines, etc, along with the
    // numbers of the degrees of freedom in hierarchical and lexicographical
    // order
    switch (dim)
      {
        case 0:
          {
            h2l[0] = 0;
            break;
          }
        case 1:
          {
            h2l[0] = 0;
            h2l[1] = dofs_per_cell - 1;
            for (unsigned int i = 2; i < dofs_per_cell; ++i)
              h2l[i] = i - 1;

            break;
          }

        case 2:
          {
            unsigned int next_index = 0;
            // first the four vertices
            h2l[next_index++] = 0;
            h2l[next_index++] = n - 1;
            h2l[next_index++] = n * (n - 1);
            h2l[next_index++] = n * n - 1;

            // left   line
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (1 + i) * n;

            // right  line
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (2 + i) * n - 1;

            // bottom line
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = 1 + i;

            // top    line
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n * (n - 1) + i + 1;

            // inside quad
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = n * (i + 1) + j + 1;

            Assert(next_index == dofs_per_cell, ExcInternalError());

            break;
          }

        case 3:
          {
            unsigned int next_index = 0;
            // first the eight vertices
            h2l[next_index++] = 0;                        // 0
            h2l[next_index++] = (1) * degree;             // 1
            h2l[next_index++] = (n)*degree;               // 2
            h2l[next_index++] = (n + 1) * degree;         // 3
            h2l[next_index++] = (n * n) * degree;         // 4
            h2l[next_index++] = (n * n + 1) * degree;     // 5
            h2l[next_index++] = (n * n + n) * degree;     // 6
            h2l[next_index++] = (n * n + n + 1) * degree; // 7

            // line 0
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (i + 1) * n;
            // line 1
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n - 1 + (i + 1) * n;
            // line 2
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = 1 + i;
            // line 3
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = 1 + i + n * (n - 1);

            // line 4
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (n - 1) * n * n + (i + 1) * n;
            // line 5
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (n - 1) * (n * n + 1) + (i + 1) * n;
            // line 6
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n * n * (n - 1) + i + 1;
            // line 7
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n * n * (n - 1) + i + 1 + n * (n - 1);

            // line 8
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (i + 1) * n * n;
            // line 9
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n - 1 + (i + 1) * n * n;
            // line 10
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = (i + 1) * n * n + n * (n - 1);
            // line 11
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              h2l[next_index++] = n - 1 + (i + 1) * n * n + n * (n - 1);


            // inside quads
            // face 0
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = (i + 1) * n * n + n * (j + 1);
            // face 1
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = (i + 1) * n * n + n - 1 + n * (j + 1);
            // face 2, note the orientation!
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = (j + 1) * n * n + i + 1;
            // face 3, note the orientation!
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = (j + 1) * n * n + n * (n - 1) + i + 1;
            // face 4
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = n * (i + 1) + j + 1;
            // face 5
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                h2l[next_index++] = (n - 1) * n * n + n * (i + 1) + j + 1;

            // inside hex
            for (unsigned int i = 0; i < dofs_per_line; ++i)
              for (unsigned int j = 0; j < dofs_per_line; ++j)
                for (unsigned int k = 0; k < dofs_per_line; ++k)
                  h2l[next_index++] = n * n * (i + 1) + n * (j + 1) + k + 1;

            Assert(next_index == dofs_per_cell, ExcInternalError());

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }

    return h2l;
  }



  template <int dim>
  void
  hierarchic_to_lexicographic_numbering(const unsigned int         degree,
                                        std::vector<unsigned int> &h2l)
  {
    AssertDimension(h2l.size(), Utilities::fixed_power<dim>(degree + 1));
    h2l = hierarchic_to_lexicographic_numbering<dim>(degree);
  }



  template <int dim>
  void
  hierarchic_to_lexicographic_numbering(const FiniteElementData<dim> &fe,
                                        std::vector<unsigned int> &   h2l)
  {
    Assert(h2l.size() == fe.n_dofs_per_cell(),
           ExcDimensionMismatch(h2l.size(), fe.n_dofs_per_cell()));
    hierarchic_to_lexicographic_numbering<dim>(fe.n_dofs_per_line() + 1, h2l);
  }



  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering(const FiniteElementData<dim> &fe)
  {
    Assert(fe.n_components() == 1, ExcInvalidFE());
    return hierarchic_to_lexicographic_numbering<dim>(fe.n_dofs_per_line() + 1);
  }



  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering(const unsigned int degree)
  {
    return Utilities::invert_permutation(
      hierarchic_to_lexicographic_numbering<dim>(degree));
  }



  template <int dim>
  void
  lexicographic_to_hierarchic_numbering(const FiniteElementData<dim> &fe,
                                        std::vector<unsigned int> &   l2h)
  {
    l2h = lexicographic_to_hierarchic_numbering(fe);
  }



  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering(const FiniteElementData<dim> &fe)
  {
    return Utilities::invert_permutation(
      hierarchic_to_lexicographic_numbering(fe));
  }

} // namespace FETools


DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_fe_tools_templates_H */
