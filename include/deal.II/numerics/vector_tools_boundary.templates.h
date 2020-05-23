// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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


#ifndef dealii_vector_tools_boundary_templates_h
#define dealii_vector_tools_boundary_templates_h

#include <deal.II/base/qprojector.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  // ----------- interpolate_boundary_values for std::map --------------------

  namespace internal
  {
    template <int dim,
              int spacedim,
              typename number,
              template <int, int> class DoFHandlerType,
              template <int, int> class M_or_MC>
    static inline void
    do_interpolate_boundary_values(
      const M_or_MC<dim, spacedim> &       mapping,
      const DoFHandlerType<dim, spacedim> &dof,
      const std::map<types::boundary_id, const Function<spacedim, number> *>
        &                                        function_map,
      std::map<types::global_dof_index, number> &boundary_values,
      const ComponentMask &                      component_mask)
    {
      Assert(
        component_mask.represents_n_components(dof.get_fe(0).n_components()),
        ExcMessage("The number of components in the mask has to be either "
                   "zero or equal to the number of components in the finite "
                   "element."));


      // if for whatever reason we were passed an empty map, return
      // immediately
      if (function_map.size() == 0)
        return;

      Assert(function_map.find(numbers::internal_face_boundary_id) ==
               function_map.end(),
             ExcMessage("You cannot specify the special boundary indicator "
                        "for interior faces in your function map."));

      const unsigned int n_components = DoFTools::n_components(dof);
      for (typename std::map<types::boundary_id,
                             const Function<spacedim, number> *>::const_iterator
             i = function_map.begin();
           i != function_map.end();
           ++i)
        Assert(n_components == i->second->n_components,
               ExcDimensionMismatch(n_components, i->second->n_components));


      // interpolate boundary values in 1d. in higher dimensions, we
      // use FEValues to figure out what to do on faces, but in 1d
      // faces are points and it is far easier to simply work on
      // individual vertices
      if (dim == 1)
        {
          for (const auto &cell : dof.active_cell_iterators())
            for (const unsigned int direction :
                 GeometryInfo<dim>::face_indices())
              if (cell->at_boundary(direction) &&
                  (function_map.find(cell->face(direction)->boundary_id()) !=
                   function_map.end()))
                {
                  const Function<spacedim, number> &boundary_function =
                    *function_map.find(cell->face(direction)->boundary_id())
                       ->second;

                  // get the FE corresponding to this cell
                  const FiniteElement<dim, spacedim> &fe = cell->get_fe();
                  Assert(fe.n_components() == boundary_function.n_components,
                         ExcDimensionMismatch(fe.n_components(),
                                              boundary_function.n_components));

                  Assert(component_mask.n_selected_components(
                           fe.n_components()) > 0,
                         ComponentMask::ExcNoComponentSelected());

                  // now set the value of the vertex degree of
                  // freedom. setting also creates the entry in the
                  // map if it did not exist beforehand
                  //
                  // save some time by requesting values only once for
                  // each point, irrespective of the number of
                  // components of the function
                  Vector<number> function_values(fe.n_components());
                  if (fe.n_components() == 1)
                    function_values(0) =
                      boundary_function.value(cell->vertex(direction));
                  else
                    boundary_function.vector_value(cell->vertex(direction),
                                                   function_values);

                  for (unsigned int i = 0; i < fe.dofs_per_vertex; ++i)
                    if (component_mask[fe.face_system_to_component_index(i)
                                         .first])
                      boundary_values[cell->vertex_dof_index(
                        direction, i, cell->active_fe_index())] =
                        function_values(
                          fe.face_system_to_component_index(i).first);
                }
        }
      else // dim > 1
        {
          const bool fe_is_system = (n_components != 1);

          // field to store the indices
          std::vector<types::global_dof_index> face_dofs;
          face_dofs.reserve(DoFTools::max_dofs_per_face(dof));

          // array to store the values of the boundary function at the boundary
          // points. have two arrays for scalar and vector functions to use the
          // more efficient one respectively
          std::vector<number>         dof_values_scalar;
          std::vector<Vector<number>> dof_values_system;
          dof_values_scalar.reserve(DoFTools::max_dofs_per_face(dof));
          dof_values_system.reserve(DoFTools::max_dofs_per_face(dof));

          // before we start with the loop over all cells create an hp::FEValues
          // object that holds the interpolation points of all finite elements
          // that may ever be in use
          const dealii::hp::FECollection<dim, spacedim> &finite_elements =
            dof.get_fe_collection();
          dealii::hp::QCollection<dim - 1> q_collection;
          for (unsigned int f = 0; f < finite_elements.size(); ++f)
            {
              const FiniteElement<dim, spacedim> &fe = finite_elements[f];

              // generate a quadrature rule on the face from the unit support
              // points. this will be used to obtain the quadrature points on
              // the real cell's face
              //
              // to do this, we check whether the FE has support points on the
              // face at all:
              if (fe.has_face_support_points())
                q_collection.push_back(
                  Quadrature<dim - 1>(fe.get_unit_face_support_points()));
              else
                {
                  // if not, then we should try a more clever way. the idea is
                  // that a finite element may not offer support points for all
                  // its shape functions, but maybe only some. if it offers
                  // support points for the components we are interested in in
                  // this function, then that's fine. if not, the function we
                  // call in the finite element will raise an exception. the
                  // support points for the other shape functions are left
                  // uninitialized (well, initialized by the default
                  // constructor), since we don't need them anyway.
                  //
                  // As a detour, we must make sure we only query
                  // face_system_to_component_index if the index corresponds to
                  // a primitive shape function. since we know that all the
                  // components we are interested in are primitive (by the above
                  // check), we can safely put such a check in front
                  std::vector<Point<dim - 1>> unit_support_points(
                    fe.dofs_per_face);

                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (fe.is_primitive(fe.face_to_cell_index(i, 0)))
                      if (component_mask[fe.face_system_to_component_index(i)
                                           .first] == true)
                        unit_support_points[i] = fe.unit_face_support_point(i);

                  q_collection.push_back(
                    Quadrature<dim - 1>(unit_support_points));
                }
            }
          // now that we have a q_collection object with all the right
          // quadrature points, create an hp::FEFaceValues object that we can
          // use to evaluate the boundary values at
          const auto mapping_collection =
            dealii::hp::MappingCollection<dim, spacedim>(mapping);
          dealii::hp::FEFaceValues<dim, spacedim> x_fe_values(
            mapping_collection,
            finite_elements,
            q_collection,
            update_quadrature_points);

          typename DoFHandlerType<dim, spacedim>::active_cell_iterator
            cell = dof.begin_active(),
            endc = dof.end();
          for (; cell != endc; ++cell)
            if (!cell->is_artificial())
              for (const unsigned int face_no :
                   GeometryInfo<dim>::face_indices())
                {
                  const FiniteElement<dim, spacedim> &fe = cell->get_fe();

                  // we can presently deal only with primitive elements for
                  // boundary values. this does not preclude us using
                  // non-primitive elements in components that we aren't
                  // interested in, however. make sure that all shape functions
                  // that are non-zero for the components we are interested in,
                  // are in fact primitive
                  for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell;
                       ++i)
                    {
                      const ComponentMask &nonzero_component_array =
                        cell->get_fe().get_nonzero_components(i);
                      for (unsigned int c = 0; c < n_components; ++c)
                        if ((nonzero_component_array[c] == true) &&
                            (component_mask[c] == true))
                          Assert(
                            cell->get_fe().is_primitive(i),
                            ExcMessage(
                              "This function can only deal with requested boundary "
                              "values that correspond to primitive (scalar) base "
                              "elements. You may want to look up in the deal.II "
                              "glossary what the term 'primitive' means."
                              "\n\n"
                              "There are alternative boundary value interpolation "
                              "functions in namespace 'VectorTools' that you can "
                              "use for non-primitive finite elements."));
                    }

                  const typename DoFHandlerType<dim, spacedim>::face_iterator
                                           face = cell->face(face_no);
                  const types::boundary_id boundary_component =
                    face->boundary_id();

                  // see if this face is part of the boundaries for which we are
                  // supposed to do something, and also see if the finite
                  // element in use here has DoFs on the face at all
                  if ((function_map.find(boundary_component) !=
                       function_map.end()) &&
                      (cell->get_fe().dofs_per_face > 0))
                    {
                      // face is of the right component
                      x_fe_values.reinit(cell, face_no);
                      const dealii::FEFaceValues<dim, spacedim> &fe_values =
                        x_fe_values.get_present_fe_values();

                      // get indices, physical location and boundary values of
                      // dofs on this face
                      face_dofs.resize(fe.dofs_per_face);
                      face->get_dof_indices(face_dofs, cell->active_fe_index());
                      const std::vector<Point<spacedim>> &dof_locations =
                        fe_values.get_quadrature_points();

                      if (fe_is_system)
                        {
                          // resize array. avoid construction of a memory
                          // allocating temporary if possible
                          if (dof_values_system.size() < fe.dofs_per_face)
                            dof_values_system.resize(fe.dofs_per_face,
                                                     Vector<number>(
                                                       fe.n_components()));
                          else
                            dof_values_system.resize(fe.dofs_per_face);

                          function_map.find(boundary_component)
                            ->second->vector_value_list(dof_locations,
                                                        dof_values_system);

                          // enter those dofs into the list that match the
                          // component signature. avoid the usual complication
                          // that we can't just use *_system_to_component_index
                          // for non-primitive FEs
                          for (unsigned int i = 0; i < face_dofs.size(); ++i)
                            {
                              unsigned int component;
                              if (fe.is_primitive())
                                component =
                                  fe.face_system_to_component_index(i).first;
                              else
                                {
                                  // non-primitive case. make sure that this
                                  // particular shape function _is_ primitive,
                                  // and get at it's component. use usual trick
                                  // to transfer face dof index to cell dof
                                  // index
                                  const unsigned int cell_i =
                                    (dim == 1 ?
                                       i :
                                       (dim == 2 ?
                                          (i < 2 * fe.dofs_per_vertex ?
                                             i :
                                             i + 2 * fe.dofs_per_vertex) :
                                          (dim == 3 ?
                                             (i < 4 * fe.dofs_per_vertex ?
                                                i :
                                                (i < 4 * fe.dofs_per_vertex +
                                                       4 * fe.dofs_per_line ?
                                                   i + 4 * fe.dofs_per_vertex :
                                                   i + 4 * fe.dofs_per_vertex +
                                                     8 * fe.dofs_per_line)) :
                                             numbers::invalid_unsigned_int)));
                                  Assert(cell_i < fe.dofs_per_cell,
                                         ExcInternalError());

                                  // make sure that if this is not a primitive
                                  // shape function, then all the corresponding
                                  // components in the mask are not set
                                  if (!fe.is_primitive(cell_i))
                                    for (unsigned int c = 0; c < n_components;
                                         ++c)
                                      if (fe.get_nonzero_components(cell_i)[c])
                                        Assert(component_mask[c] == false,
                                               FETools::ExcFENotPrimitive());

                                  // let's pick the first of possibly more than
                                  // one non-zero components. if shape function
                                  // is non-primitive, then we will ignore the
                                  // result in the following anyway, otherwise
                                  // there's only one non-zero component which
                                  // we will use
                                  component = fe.get_nonzero_components(cell_i)
                                                .first_selected_component();
                                }

                              if (component_mask[component] == true)
                                boundary_values[face_dofs[i]] =
                                  dof_values_system[i](component);
                            }
                        }
                      else
                        // fe has only one component, so save some computations
                        {
                          // get only the one component that this function has
                          dof_values_scalar.resize(fe.dofs_per_face);
                          function_map.find(boundary_component)
                            ->second->value_list(dof_locations,
                                                 dof_values_scalar,
                                                 0);

                          // enter into list

                          for (unsigned int i = 0; i < face_dofs.size(); ++i)
                            boundary_values[face_dofs[i]] =
                              dof_values_scalar[i];
                        }
                    }
                }
        }
    } // end of interpolate_boundary_values
  }   // namespace internal



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &                      component_mask_)
  {
    internal::do_interpolate_boundary_values(
      mapping, dof, function_map, boundary_values, component_mask_);
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &             mapping,
    const DoFHandlerType<dim, spacedim> &      dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &                      component_mask)
  {
    std::map<types::boundary_id, const Function<spacedim, number> *>
      function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values(
      mapping, dof, function_map, boundary_values, component_mask);
  }


  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &                      component_mask_)
  {
    internal::do_interpolate_boundary_values(
      mapping, dof, function_map, boundary_values, component_mask_);
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &      dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &                      component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                                dof,
                                boundary_component,
                                boundary_function,
                                boundary_values,
                                component_mask);
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &                      component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                                dof,
                                function_map,
                                boundary_values,
                                component_mask);
  }



  // ----------- interpolate_boundary_values for AffineConstraints
  // --------------



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask_)
  {
    std::map<types::global_dof_index, number> boundary_values;
    interpolate_boundary_values(
      mapping, dof, function_map, boundary_values, component_mask_);
    typename std::map<types::global_dof_index, number>::const_iterator
      boundary_value = boundary_values.begin();
    for (; boundary_value != boundary_values.end(); ++boundary_value)
      {
        if (constraints.can_store_line(boundary_value->first) &&
            !constraints.is_constrained(boundary_value->first))
          {
            constraints.add_line(boundary_value->first);
            constraints.set_inhomogeneity(boundary_value->first,
                                          boundary_value->second);
          }
      }
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const types::boundary_id             boundary_component,
    const Function<spacedim, number> &   boundary_function,
    AffineConstraints<number> &          constraints,
    const ComponentMask &                component_mask)
  {
    std::map<types::boundary_id, const Function<spacedim, number> *>
      function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values(
      mapping, dof, function_map, constraints, component_mask);
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const types::boundary_id             boundary_component,
    const Function<spacedim, number> &   boundary_function,
    AffineConstraints<number> &          constraints,
    const ComponentMask &                component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                                dof,
                                boundary_component,
                                boundary_function,
                                constraints,
                                component_mask);
  }



  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                                dof,
                                function_map,
                                constraints,
                                component_mask);
  }



  // -------- implementation for project_boundary_values with std::map --------


  namespace internal
  {
    // keep the first argument non-reference since we use it
    // with 1e-8 * number
    template <typename number1, typename number2>
    bool
    real_part_bigger_than(const number1 a, const number2 &b)
    {
      return a > b;
    }

    template <typename number1, typename number2>
    bool
    real_part_bigger_than(const number1 a, const std::complex<number2> b)
    {
      Assert(std::abs(b.imag()) <= 1e-15 * std::abs(b), ExcInternalError());
      return a > b.real();
    }

    template <typename number1, typename number2>
    bool
    real_part_bigger_than(const std::complex<number1> a, const number2 b)
    {
      Assert(std::abs(a.imag()) <= 1e-15 * std::abs(a), ExcInternalError());
      return a.real() > b;
    }

    template <typename number1, typename number2>
    bool
    real_part_bigger_than(const std::complex<number1> a,
                          const std::complex<number2> b)
    {
      Assert(std::abs(a.imag()) <= 1e-15 * std::abs(a), ExcInternalError());
      Assert(std::abs(b.imag()) <= 1e-15 * std::abs(b), ExcInternalError());
      return a.real() > b.real();
    }

    // this function is needed to get an idea where
    // rhs.norm_sqr()  is too small for a given type.
    template <typename number>
    number
    min_number(const number & /*dummy*/)
    {
      return std::numeric_limits<number>::min();
    }

    // Sine rhs.norm_sqr() is non-negative real, in complex case we
    // take the numeric limits of the underlying type used in std::complex<>.
    template <typename number>
    number
    min_number(const std::complex<number> & /*dummy*/)
    {
      return std::numeric_limits<number>::min();
    }

    template <int dim,
              int spacedim,
              template <int, int> class DoFHandlerType,
              template <int, int> class M_or_MC,
              template <int> class Q_or_QC,
              typename number>
    void
    do_project_boundary_values(
      const M_or_MC<dim, spacedim> &       mapping,
      const DoFHandlerType<dim, spacedim> &dof,
      const std::map<types::boundary_id, const Function<spacedim, number> *>
        &                                        boundary_functions,
      const Q_or_QC<dim - 1> &                   q,
      std::map<types::global_dof_index, number> &boundary_values,
      std::vector<unsigned int>                  component_mapping)
    {
      // in 1d, projection onto the 0d end points == interpolation
      if (dim == 1)
        {
          Assert(component_mapping.size() == 0, ExcNotImplemented());
          interpolate_boundary_values(
            mapping, dof, boundary_functions, boundary_values, ComponentMask());
          return;
        }

      // TODO:[?] In project_boundary_values, no condensation of sparsity
      //    structures, matrices and right hand sides or distribution of
      //    solution vectors is performed. This is ok for dim<3 because then
      //    there are no constrained nodes on the boundary, but is not
      //    acceptable for higher dimensions. Fix this.

      if (component_mapping.size() == 0)
        {
          AssertDimension(dof.get_fe(0).n_components(),
                          boundary_functions.begin()->second->n_components);
          // I still do not see why i
          // should create another copy
          // here
          component_mapping.resize(dof.get_fe(0).n_components());
          for (unsigned int i = 0; i < component_mapping.size(); ++i)
            component_mapping[i] = i;
        }
      else
        AssertDimension(dof.get_fe(0).n_components(), component_mapping.size());

      std::vector<types::global_dof_index> dof_to_boundary_mapping;
      std::set<types::boundary_id>         selected_boundary_components;
      for (typename std::map<types::boundary_id,
                             const Function<spacedim, number> *>::const_iterator
             i = boundary_functions.begin();
           i != boundary_functions.end();
           ++i)
        selected_boundary_components.insert(i->first);

      DoFTools::map_dof_to_boundary_indices(dof,
                                            selected_boundary_components,
                                            dof_to_boundary_mapping);

      // Done if no degrees of freedom on the boundary
      if (dof.n_boundary_dofs(boundary_functions) == 0)
        return;

      // set up sparsity structure
      DynamicSparsityPattern dsp(dof.n_boundary_dofs(boundary_functions),
                                 dof.n_boundary_dofs(boundary_functions));
      DoFTools::make_boundary_sparsity_pattern(dof,
                                               boundary_functions,
                                               dof_to_boundary_mapping,
                                               dsp);
      SparsityPattern sparsity;
      sparsity.copy_from(dsp);



      // note: for three or more dimensions, there
      // may be constrained nodes on the boundary
      // in this case the boundary mass matrix has
      // to be condensed and the solution is to
      // be distributed afterwards, which is not
      // yet implemented. The reason for this is
      // that we cannot simply use the condense
      // family of functions, since the matrices
      // and vectors do not use the global
      // numbering but rather the boundary
      // numbering, i.e. the condense function
      // needs to use another indirection. There
      // should be not many technical problems,
      // but it needs to be implemented
      if (dim >= 3)
        {
#ifdef DEBUG
          // Assert that there are no hanging nodes at the boundary
          int level = -1;
          for (const auto &cell : dof.active_cell_iterators())
            for (auto f : GeometryInfo<dim>::face_indices())
              {
                if (cell->at_boundary(f))
                  {
                    if (level == -1)
                      level = cell->level();
                    else
                      {
                        Assert(
                          level == cell->level(),
                          ExcMessage(
                            "The mesh you use in projecting boundary values "
                            "has hanging nodes at the boundary. This would require "
                            "dealing with hanging node constraints when solving "
                            "the linear system on the boundary, but this is not "
                            "currently implemented."));
                      }
                  }
              }
#endif
        }
      sparsity.compress();


      // make mass matrix and right hand side
      SparseMatrix<number> mass_matrix(sparsity);
      Vector<number>       rhs(sparsity.n_rows());


      MatrixCreator::create_boundary_mass_matrix(
        mapping,
        dof,
        q,
        mass_matrix,
        boundary_functions,
        rhs,
        dof_to_boundary_mapping,
        static_cast<const Function<spacedim, number> *>(nullptr),
        component_mapping);

      Vector<number> boundary_projection(rhs.size());

      // cannot reduce residual in a useful way if we are close to the square
      // root of the minimal double value
      if (rhs.norm_sqr() < 1e28 * min_number(number()))
        boundary_projection = 0;
      else
        {
          // Allow for a maximum of 5*n steps to reduce the residual by 10^-12.
          // n steps may not be sufficient, since roundoff errors may accumulate
          // for badly conditioned matrices
          ReductionControl control(5 * rhs.size(), 0., 1e-12, false, false);
          GrowingVectorMemory<Vector<number>> memory;
          SolverCG<Vector<number>>            cg(control, memory);

          PreconditionSSOR<SparseMatrix<number>> prec;
          prec.initialize(mass_matrix, 1.2);

          cg.solve(mass_matrix, boundary_projection, rhs, prec);
        }
      // fill in boundary values
      for (unsigned int i = 0; i < dof_to_boundary_mapping.size(); ++i)
        if (dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
          {
            AssertIsFinite(boundary_projection(dof_to_boundary_mapping[i]));

            // this dof is on one of the
            // interesting boundary parts
            //
            // remember: i is the global dof
            // number, dof_to_boundary_mapping[i]
            // is the number on the boundary and
            // thus in the solution vector
            boundary_values[i] =
              boundary_projection(dof_to_boundary_mapping[i]);
          }
    }
  } // namespace internal

  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    internal::do_project_boundary_values(
      mapping, dof, boundary_functions, q, boundary_values, component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                            dof,
                            boundary_functions,
                            q,
                            boundary_values,
                            component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const hp::QCollection<dim - 1> &           q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    internal::do_project_boundary_values(
      mapping, dof, boundary_functions, q, boundary_values, component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_function,
    const hp::QCollection<dim - 1> &           q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    project_boundary_values(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof,
      boundary_function,
      q,
      boundary_values,
      component_mapping);
  }


  // ---- implementation for project_boundary_values with AffineConstraints ----



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_functions,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping)
  {
    std::map<types::global_dof_index, number> boundary_values;
    project_boundary_values(
      mapping, dof, boundary_functions, q, boundary_values, component_mapping);
    typename std::map<types::global_dof_index, number>::const_iterator
      boundary_value = boundary_values.begin();
    for (; boundary_value != boundary_values.end(); ++boundary_value)
      {
        if (!constraints.is_constrained(boundary_value->first))
          {
            constraints.add_line(boundary_value->first);
            constraints.set_inhomogeneity(boundary_value->first,
                                          boundary_value->second);
          }
      }
  }



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_functions,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim, spacedim>::mapping,
                            dof,
                            boundary_functions,
                            q,
                            constraints,
                            component_mapping);
  }


  namespace internals
  {
    // This function computes the
    // projection of the boundary
    // function on edges for 3D.
    template <typename cell_iterator>
    void
    compute_edge_projection(const cell_iterator &cell,
                            const unsigned int   face,
                            const unsigned int   line,
                            hp::FEValues<3> &    hp_fe_values,
                            const Function<3> &  boundary_function,
                            const unsigned int   first_vector_component,
                            std::vector<double> &dof_values,
                            std::vector<bool> &  dofs_processed)
    {
      const double tol =
        0.5 * cell->face(face)->line(line)->diameter() / cell->get_fe().degree;
      const unsigned int dim      = 3;
      const unsigned int spacedim = 3;

      hp_fe_values.reinit(
        cell,
        (cell->active_fe_index() * GeometryInfo<dim>::faces_per_cell + face) *
            GeometryInfo<dim>::lines_per_face +
          line);

      // Initialize the required
      // objects.
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
      const FiniteElement<dim> &                           fe = cell->get_fe();
      const std::vector<DerivativeForm<1, dim, spacedim>> &jacobians =
        fe_values.get_jacobians();
      const std::vector<Point<dim>> &quadrature_points =
        fe_values.get_quadrature_points();

      std::vector<Tensor<1, dim>> tangentials(fe_values.n_quadrature_points);
      std::vector<Vector<double>> values(fe_values.n_quadrature_points,
                                         Vector<double>(fe.n_components()));

      // Get boundary function values
      // at quadrature points.
      boundary_function.vector_value_list(quadrature_points, values);

      const std::vector<Point<dim>> &reference_quadrature_points =
        fe_values.get_quadrature().get_points();
      std::pair<unsigned int, unsigned int> base_indices(0, 0);

      if (dynamic_cast<const FESystem<dim> *>(&cell->get_fe()) != nullptr)
        {
          unsigned int fe_index     = 0;
          unsigned int fe_index_old = 0;
          unsigned int i            = 0;

          for (; i < fe.n_base_elements(); ++i)
            {
              fe_index_old = fe_index;
              fe_index +=
                fe.element_multiplicity(i) * fe.base_element(i).n_components();

              if (fe_index > first_vector_component)
                break;
            }

          base_indices.first  = i;
          base_indices.second = (first_vector_component - fe_index_old) /
                                fe.base_element(i).n_components();
        }

      // coordinate directions of
      // the edges of the face.
      const unsigned int
        edge_coordinate_direction[GeometryInfo<dim>::faces_per_cell]
                                 [GeometryInfo<dim>::lines_per_face] = {
                                   {2, 2, 1, 1},
                                   {2, 2, 1, 1},
                                   {0, 0, 2, 2},
                                   {0, 0, 2, 2},
                                   {1, 1, 0, 0},
                                   {1, 1, 0, 0}};
      const FEValuesExtractors::Vector vec(first_vector_component);

      // The interpolation for the
      // lowest order edge shape
      // functions is just the mean
      // value of the tangential
      // components of the boundary
      // function on the edge.
      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
           ++q_point)
        {
          // Therefore compute the
          // tangential of the edge at
          // the quadrature point.
          Point<dim> shifted_reference_point_1 =
            reference_quadrature_points[q_point];
          Point<dim> shifted_reference_point_2 =
            reference_quadrature_points[q_point];

          shifted_reference_point_1(edge_coordinate_direction[face][line]) +=
            tol;
          shifted_reference_point_2(edge_coordinate_direction[face][line]) -=
            tol;
          tangentials[q_point] =
            (0.5 *
             (fe_values.get_mapping().transform_unit_to_real_cell(
                cell, shifted_reference_point_1) -
              fe_values.get_mapping().transform_unit_to_real_cell(
                cell, shifted_reference_point_2)) /
             tol);
          tangentials[q_point] /= tangentials[q_point].norm();

          // Compute the degrees of
          // freedom.
          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            if (((dynamic_cast<const FESystem<dim> *>(&fe) != nullptr) &&
                 (fe.system_to_base_index(fe.face_to_cell_index(i, face))
                    .first == base_indices) &&
                 (fe.base_element(base_indices.first)
                    .face_to_cell_index(line * fe.degree, face) <=
                  fe.system_to_base_index(fe.face_to_cell_index(i, face))
                    .second) &&
                 (fe.system_to_base_index(fe.face_to_cell_index(i, face))
                    .second <=
                  fe.base_element(base_indices.first)
                    .face_to_cell_index((line + 1) * fe.degree - 1, face))) ||
                ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr) &&
                 (line * fe.degree <= i) && (i < (line + 1) * fe.degree)))
              {
                const double tangential_solution_component =
                  (values[q_point](first_vector_component) *
                     tangentials[q_point][0] +
                   values[q_point](first_vector_component + 1) *
                     tangentials[q_point][1] +
                   values[q_point](first_vector_component + 2) *
                     tangentials[q_point][2]);
                dof_values[i] +=
                  (fe_values.JxW(q_point) * tangential_solution_component *
                   (fe_values[vec].value(fe.face_to_cell_index(i, face),
                                         q_point) *
                    tangentials[q_point]) /
                   std::sqrt(
                     jacobians[q_point][0]
                              [edge_coordinate_direction[face][line]] *
                       jacobians[q_point][0]
                                [edge_coordinate_direction[face][line]] +
                     jacobians[q_point][1]
                              [edge_coordinate_direction[face][line]] *
                       jacobians[q_point][1]
                                [edge_coordinate_direction[face][line]] +
                     jacobians[q_point][2]
                              [edge_coordinate_direction[face][line]] *
                       jacobians[q_point][2]
                                [edge_coordinate_direction[face][line]]));

                if (q_point == 0)
                  dofs_processed[i] = true;
              }
        }
    }

    // dummy implementation of above
    // function for all other
    // dimensions
    template <int dim, typename cell_iterator>
    void
    compute_edge_projection(const cell_iterator &,
                            const unsigned int,
                            const unsigned int,
                            hp::FEValues<dim> &,
                            const Function<dim> &,
                            const unsigned int,
                            std::vector<double> &,
                            std::vector<bool> &)
    {
      Assert(false, ExcInternalError());
    }

    // This function computes the
    // projection of the boundary
    // function on the interior of
    // faces.
    template <int dim, typename cell_iterator, typename number>
    void
    compute_face_projection_curl_conforming(
      const cell_iterator &        cell,
      const unsigned int           face,
      hp::FEValues<dim> &          hp_fe_values,
      const Function<dim, number> &boundary_function,
      const unsigned int           first_vector_component,
      std::vector<double> &        dof_values,
      std::vector<bool> &          dofs_processed)
    {
      const unsigned int spacedim = dim;
      hp_fe_values.reinit(cell,
                          cell->active_fe_index() *
                              GeometryInfo<dim>::faces_per_cell +
                            face);
      // Initialize the required
      // objects.
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
      const FiniteElement<dim> &                           fe = cell->get_fe();
      const std::vector<DerivativeForm<1, dim, spacedim>> &jacobians =
        fe_values.get_jacobians();
      const std::vector<Point<dim>> &quadrature_points =
        fe_values.get_quadrature_points();
      const unsigned int                    degree = fe.degree - 1;
      std::pair<unsigned int, unsigned int> base_indices(0, 0);

      if (dynamic_cast<const FESystem<dim> *>(&cell->get_fe()) != nullptr)
        {
          unsigned int fe_index     = 0;
          unsigned int fe_index_old = 0;
          unsigned int i            = 0;

          for (; i < fe.n_base_elements(); ++i)
            {
              fe_index_old = fe_index;
              fe_index +=
                fe.element_multiplicity(i) * fe.base_element(i).n_components();

              if (fe_index > first_vector_component)
                break;
            }

          base_indices.first  = i;
          base_indices.second = (first_vector_component - fe_index_old) /
                                fe.base_element(i).n_components();
        }

      std::vector<Vector<double>> values(fe_values.n_quadrature_points,
                                         Vector<double>(fe.n_components()));

      // Get boundary function
      // values at quadrature
      // points.
      boundary_function.vector_value_list(quadrature_points, values);

      switch (dim)
        {
          case 2:
            {
              const double tol =
                0.5 * cell->face(face)->diameter() / cell->get_fe().degree;
              std::vector<Tensor<1, dim>> tangentials(
                fe_values.n_quadrature_points);

              const std::vector<Point<dim>> &reference_quadrature_points =
                fe_values.get_quadrature().get_points();

              // coordinate directions
              // of the face.
              const unsigned int
                face_coordinate_direction[GeometryInfo<dim>::faces_per_cell] = {
                  1, 1, 0, 0};
              const FEValuesExtractors::Vector vec(first_vector_component);

              // The interpolation for
              // the lowest order face
              // shape functions is just
              // the mean value of the
              // tangential  components
              // of the boundary function
              // on the edge.
              for (unsigned int q_point = 0;
                   q_point < fe_values.n_quadrature_points;
                   ++q_point)
                {
                  // Therefore compute the
                  // tangential of the
                  // face at the quadrature
                  // point.
                  Point<dim> shifted_reference_point_1 =
                    reference_quadrature_points[q_point];
                  Point<dim> shifted_reference_point_2 =
                    reference_quadrature_points[q_point];

                  shifted_reference_point_1(face_coordinate_direction[face]) +=
                    tol;
                  shifted_reference_point_2(face_coordinate_direction[face]) -=
                    tol;
                  tangentials[q_point] =
                    (fe_values.get_mapping().transform_unit_to_real_cell(
                       cell, shifted_reference_point_1) -
                     fe_values.get_mapping().transform_unit_to_real_cell(
                       cell, shifted_reference_point_2)) /
                    tol;
                  tangentials[q_point] /= tangentials[q_point].norm();

                  // Compute the degrees
                  // of freedom.
                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (((dynamic_cast<const FESystem<dim> *>(&fe) !=
                          nullptr) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .first == base_indices)) ||
                        (dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr))
                      {
                        dof_values[i] +=
                          fe_values.JxW(q_point) *
                          (values[q_point](first_vector_component) *
                             tangentials[q_point][0] +
                           values[q_point](first_vector_component + 1) *
                             tangentials[q_point][1]) *
                          (fe_values[vec].value(fe.face_to_cell_index(i, face),
                                                q_point) *
                           tangentials[q_point]);

                        if (q_point == 0)
                          dofs_processed[i] = true;
                      }
                }

              break;
            }

          case 3:
            {
              const FEValuesExtractors::Vector vec(first_vector_component);
              FullMatrix<double>               assembling_matrix(
                degree * fe.degree, dim * fe_values.n_quadrature_points);
              Vector<double>     assembling_vector(assembling_matrix.n());
              Vector<double>     cell_rhs(assembling_matrix.m());
              FullMatrix<double> cell_matrix(assembling_matrix.m(),
                                             assembling_matrix.m());
              FullMatrix<double> cell_matrix_inv(assembling_matrix.m(),
                                                 assembling_matrix.m());
              Vector<double>     solution(cell_matrix.m());

              // Get coordinate directions
              // of the face.
              const unsigned int global_face_coordinate_directions
                [GeometryInfo<3>::faces_per_cell][2] = {
                  {1, 2}, {1, 2}, {2, 0}, {2, 0}, {0, 1}, {0, 1}};

              // The projection is divided into two steps.  In the first step we
              // project the boundary function on the horizontal shape
              // functions. Then the boundary function is projected on the
              // vertical shape functions.  We begin with the horizontal shape
              // functions and set up a linear system of equations to get the
              // values for degrees of freedom associated with the interior of
              // the face.
              for (unsigned int q_point = 0;
                   q_point < fe_values.n_quadrature_points;
                   ++q_point)
                {
                  // The right hand
                  // side of the
                  // corresponding problem
                  // is the residual
                  // of the boundary
                  // function and
                  // the already
                  // interpolated part
                  // on the edges.
                  Tensor<1, dim> tmp;

                  for (unsigned int d = 0; d < dim; ++d)
                    tmp[d] = values[q_point](first_vector_component + d);

                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (((dynamic_cast<const FESystem<dim> *>(&fe) !=
                          nullptr) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .first == base_indices) &&
                         (fe.base_element(base_indices.first)
                            .face_to_cell_index(2 * fe.degree, face) <=
                          fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second <=
                          fe.base_element(base_indices.first)
                            .face_to_cell_index(4 * fe.degree - 1, face))) ||
                        ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                          nullptr) &&
                         (2 * fe.degree <= i) && (i < 4 * fe.degree)))
                      tmp -=
                        dof_values[i] *
                        fe_values[vec].value(fe.face_to_cell_index(i, face),
                                             q_point);

                  const double JxW = std::sqrt(
                    fe_values.JxW(q_point) /
                    ((jacobians[q_point][0]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][0]
                                 [global_face_coordinate_directions[face][0]] +
                      jacobians[q_point][1]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][1]
                                 [global_face_coordinate_directions[face][0]] +
                      jacobians[q_point][2]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][2]
                                 [global_face_coordinate_directions[face][0]]) *
                     (jacobians[q_point][0]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][0]
                                 [global_face_coordinate_directions[face][1]] +
                      jacobians[q_point][1]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][1]
                                 [global_face_coordinate_directions[face][1]] +
                      jacobians[q_point][2]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][2]
                                 [global_face_coordinate_directions[face]
                                                                   [1]])));

                  // In the weak form
                  // the right hand
                  // side function
                  // is multiplicated
                  // by the horizontal
                  // shape functions
                  // defined in the
                  // interior of
                  // the face.
                  for (unsigned int d = 0; d < dim; ++d)
                    assembling_vector(dim * q_point + d) = JxW * tmp[d];

                  unsigned int index = 0;

                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (((dynamic_cast<const FESystem<dim> *>(&fe) !=
                          nullptr) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .first == base_indices) &&
                         (fe.base_element(base_indices.first)
                            .face_to_cell_index(
                              GeometryInfo<dim>::lines_per_face * fe.degree,
                              face) <=
                          fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second <
                          fe.base_element(base_indices.first)
                            .face_to_cell_index(
                              (degree + GeometryInfo<dim>::lines_per_face) *
                                fe.degree,
                              face))) ||
                        ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                          nullptr) &&
                         (GeometryInfo<dim>::lines_per_face * fe.degree <= i) &&
                         (i < (degree + GeometryInfo<dim>::lines_per_face) *
                                fe.degree)))
                      {
                        const Tensor<1, dim> shape_value =
                          (JxW *
                           fe_values[vec].value(fe.face_to_cell_index(i, face),
                                                q_point));

                        for (unsigned int d = 0; d < dim; ++d)
                          assembling_matrix(index, dim * q_point + d) =
                            shape_value[d];

                        ++index;
                      }
                }

              // Create the system matrix by multiplying the assembling matrix
              // with its transposed and the right hand side vector by
              // multiplying the assembling matrix with the assembling vector.
              // Invert the system matrix.
              assembling_matrix.mTmult(cell_matrix, assembling_matrix);
              cell_matrix_inv.invert(cell_matrix);
              assembling_matrix.vmult(cell_rhs, assembling_vector);
              cell_matrix_inv.vmult(solution, cell_rhs);

              // Store the computed
              // values.
              {
                unsigned int index = 0;

                for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                  if (((dynamic_cast<const FESystem<dim> *>(&fe) != nullptr) &&
                       (fe.system_to_base_index(fe.face_to_cell_index(i, face))
                          .first == base_indices) &&
                       (fe.base_element(base_indices.first)
                          .face_to_cell_index(
                            GeometryInfo<dim>::lines_per_face * fe.degree,
                            face) <=
                        fe.system_to_base_index(fe.face_to_cell_index(i, face))
                          .second) &&
                       (fe.system_to_base_index(fe.face_to_cell_index(i, face))
                          .second <
                        fe.base_element(base_indices.first)
                          .face_to_cell_index(
                            (degree + GeometryInfo<dim>::lines_per_face) *
                              fe.degree,
                            face))) ||
                      ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                        nullptr) &&
                       (GeometryInfo<dim>::lines_per_face * fe.degree <= i) &&
                       (i < (degree + GeometryInfo<dim>::lines_per_face) *
                              fe.degree)))
                    {
                      dof_values[i]     = solution(index);
                      dofs_processed[i] = true;
                      ++index;
                    }
              }

              // Now we do the same as above with the vertical shape functions
              // instead of the horizontal ones.
              for (unsigned int q_point = 0;
                   q_point < fe_values.n_quadrature_points;
                   ++q_point)
                {
                  Tensor<1, dim> tmp;

                  for (unsigned int d = 0; d < dim; ++d)
                    tmp[d] = values[q_point](first_vector_component + d);

                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (((dynamic_cast<const FESystem<dim> *>(&fe) !=
                          nullptr) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .first == base_indices) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second <=
                          fe.base_element(base_indices.first)
                            .face_to_cell_index(2 * fe.degree - 1, face)) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second >= fe.base_element(base_indices.first)
                                         .face_to_cell_index(0, face))) ||
                        ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                          nullptr) &&
                         (i < 2 * fe.degree)))
                      tmp -=
                        dof_values[i] *
                        fe_values[vec].value(fe.face_to_cell_index(i, face),
                                             q_point);

                  const double JxW = std::sqrt(
                    fe_values.JxW(q_point) /
                    ((jacobians[q_point][0]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][0]
                                 [global_face_coordinate_directions[face][0]] +
                      jacobians[q_point][1]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][1]
                                 [global_face_coordinate_directions[face][0]] +
                      jacobians[q_point][2]
                               [global_face_coordinate_directions[face][0]] *
                        jacobians[q_point][2]
                                 [global_face_coordinate_directions[face][0]]) *
                     (jacobians[q_point][0]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][0]
                                 [global_face_coordinate_directions[face][1]] +
                      jacobians[q_point][1]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][1]
                                 [global_face_coordinate_directions[face][1]] +
                      jacobians[q_point][2]
                               [global_face_coordinate_directions[face][1]] *
                        jacobians[q_point][2]
                                 [global_face_coordinate_directions[face]
                                                                   [1]])));

                  for (unsigned int d = 0; d < dim; ++d)
                    assembling_vector(dim * q_point + d) = JxW * tmp[d];

                  unsigned int index = 0;

                  for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                    if (((dynamic_cast<const FESystem<dim> *>(&fe) !=
                          nullptr) &&
                         (fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .first == base_indices) &&
                         (fe.base_element(base_indices.first)
                            .face_to_cell_index(
                              (degree + GeometryInfo<dim>::lines_per_face) *
                                fe.degree,
                              face) <=
                          fe.system_to_base_index(
                              fe.face_to_cell_index(i, face))
                            .second)) ||
                        ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                          nullptr) &&
                         ((degree + GeometryInfo<dim>::lines_per_face) *
                            fe.degree <=
                          i)))
                      {
                        const Tensor<1, dim> shape_value =
                          JxW *
                          fe_values[vec].value(fe.face_to_cell_index(i, face),
                                               q_point);

                        for (unsigned int d = 0; d < dim; ++d)
                          assembling_matrix(index, dim * q_point + d) =
                            shape_value[d];

                        ++index;
                      }
                }

              assembling_matrix.mTmult(cell_matrix, assembling_matrix);
              cell_matrix_inv.invert(cell_matrix);
              assembling_matrix.vmult(cell_rhs, assembling_vector);
              cell_matrix_inv.vmult(solution, cell_rhs);

              unsigned int index = 0;

              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim> *>(&fe) != nullptr) &&
                     (fe.system_to_base_index(fe.face_to_cell_index(i, face))
                        .first == base_indices) &&
                     (fe.base_element(base_indices.first)
                        .face_to_cell_index(
                          (degree + GeometryInfo<dim>::lines_per_face) *
                            fe.degree,
                          face) <=
                      fe.system_to_base_index(fe.face_to_cell_index(i, face))
                        .second)) ||
                    ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr) &&
                     ((degree + GeometryInfo<dim>::lines_per_face) *
                        fe.degree <=
                      i)))
                  {
                    dof_values[i]     = solution(index);
                    dofs_processed[i] = true;
                    ++index;
                  }

              break;
            }

          default:
            Assert(false, ExcNotImplemented());
        }
    }
  } // namespace internals



  template <int dim>
  void
  project_boundary_values_curl_conforming(
    const DoFHandler<dim> &    dof_handler,
    const unsigned int         first_vector_component,
    const Function<dim> &      boundary_function,
    const types::boundary_id   boundary_component,
    AffineConstraints<double> &constraints,
    const Mapping<dim> &       mapping)
  {
    // Projection-based interpolation is performed in two (in 2D) respectively
    // three (in 3D) steps. First the tangential component of the function is
    // interpolated on each edge.  This gives the values for the degrees of
    // freedom corresponding to the edge shape functions. Now we are done for
    // 2D, but in 3D we possibly have also degrees of freedom, which are
    // located in the interior of the faces. Therefore we compute the residual
    // of the function describing the boundary values and the interpolated
    // part, which we have computed in the last step. On the faces there are
    // two kinds of shape functions, the horizontal and the vertical
    // ones. Thus we have to solve two linear systems of equations of size
    // <tt>degree * (degree + 1)<tt> to obtain the values for the
    // corresponding degrees of freedom.
    const unsigned int    superdegree = dof_handler.get_fe().degree;
    const QGauss<dim - 1> reference_face_quadrature(2 * superdegree);
    const unsigned int    dofs_per_face = dof_handler.get_fe().dofs_per_face;
    const hp::FECollection<dim> &fe_collection(dof_handler.get_fe_collection());
    const hp::MappingCollection<dim> mapping_collection(mapping);
    hp::QCollection<dim>             face_quadrature_collection;

    for (unsigned int face : GeometryInfo<dim>::face_indices())
      face_quadrature_collection.push_back(
        QProjector<dim>::project_to_face(reference_face_quadrature, face));

    hp::FEValues<dim> fe_face_values(mapping_collection,
                                     fe_collection,
                                     face_quadrature_collection,
                                     update_jacobians | update_JxW_values |
                                       update_quadrature_points |
                                       update_values);

    std::vector<bool>                    dofs_processed(dofs_per_face);
    std::vector<double>                  dof_values(dofs_per_face);
    std::vector<types::global_dof_index> face_dof_indices(dofs_per_face);
    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active();

    switch (dim)
      {
        case 2:
          {
            for (; cell != dof_handler.end(); ++cell)
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // if the FE is a
                      // FE_Nothing object
                      // there is no work to
                      // do
                      if (dynamic_cast<const FE_Nothing<dim> *>(
                            &cell->get_fe()) != nullptr)
                        return;

                      // This is only
                      // implemented, if the
                      // FE is a Nedelec
                      // element. If the FE
                      // is a FESystem, we
                      // cannot check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(
                            dynamic_cast<const FE_Nedelec<dim> *>(
                              &cell->get_fe()) != nullptr,
                            (typename FiniteElement<
                              dim>::ExcInterpolationNotImplemented()));
                        }

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        {
                          dof_values[dof]     = 0.0;
                          dofs_processed[dof] = false;
                        }

                      // Compute the
                      // projection of the
                      // boundary function on
                      // the edge.
                      internals::compute_face_projection_curl_conforming(
                        cell,
                        face,
                        fe_face_values,
                        boundary_function,
                        first_vector_component,
                        dof_values,
                        dofs_processed);
                      cell->face(face)->get_dof_indices(
                        face_dof_indices, cell->active_fe_index());

                      // Add the computed constraints to the constraints
                      // object, if the degree of freedom is not already
                      // constrained.
                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        if (dofs_processed[dof] &&
                            constraints.can_store_line(face_dof_indices[dof]) &&
                            !(constraints.is_constrained(
                              face_dof_indices[dof])))
                          {
                            constraints.add_line(face_dof_indices[dof]);

                            if (std::abs(dof_values[dof]) > 1e-13)
                              constraints.set_inhomogeneity(
                                face_dof_indices[dof], dof_values[dof]);
                          }
                    }

            break;
          }

        case 3:
          {
            const QGauss<dim - 2> reference_edge_quadrature(2 * superdegree);
            const unsigned int    degree = superdegree - 1;
            hp::QCollection<dim>  edge_quadrature_collection;

            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              for (unsigned int line = 0;
                   line < GeometryInfo<dim>::lines_per_face;
                   ++line)
                edge_quadrature_collection.push_back(
                  QProjector<dim>::project_to_face(
                    QProjector<dim - 1>::project_to_face(
                      reference_edge_quadrature, line),
                    face));

            hp::FEValues<dim> fe_edge_values(mapping_collection,
                                             fe_collection,
                                             edge_quadrature_collection,
                                             update_jacobians |
                                               update_JxW_values |
                                               update_quadrature_points |
                                               update_values);

            for (; cell != dof_handler.end(); ++cell)
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // if the FE is a
                      // FE_Nothing object
                      // there is no work to
                      // do
                      if (dynamic_cast<const FE_Nothing<dim> *>(
                            &cell->get_fe()) != nullptr)
                        return;

                      // This is only
                      // implemented, if the
                      // FE is a Nedelec
                      // element. If the FE is
                      // a FESystem we cannot
                      // check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(dynamic_cast<const FE_Nedelec<dim> *>(
                                        &cell->get_fe()) != nullptr,
                                      typename FiniteElement<
                                        dim>::ExcInterpolationNotImplemented());
                        }

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        {
                          dof_values[dof]     = 0.0;
                          dofs_processed[dof] = false;
                        }

                      // First we compute the
                      // projection on the
                      // edges.
                      for (unsigned int line = 0;
                           line < GeometryInfo<3>::lines_per_face;
                           ++line)
                        internals::compute_edge_projection(
                          cell,
                          face,
                          line,
                          fe_edge_values,
                          boundary_function,
                          first_vector_component,
                          dof_values,
                          dofs_processed);

                      // If there are higher
                      // order shape
                      // functions, there is
                      // still some work
                      // left.
                      if (degree > 0)
                        internals::compute_face_projection_curl_conforming(
                          cell,
                          face,
                          fe_face_values,
                          boundary_function,
                          first_vector_component,
                          dof_values,
                          dofs_processed);

                      // Store the computed
                      // values in the global
                      // vector.
                      cell->face(face)->get_dof_indices(
                        face_dof_indices, cell->active_fe_index());

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        if (dofs_processed[dof] &&
                            constraints.can_store_line(face_dof_indices[dof]) &&
                            !(constraints.is_constrained(
                              face_dof_indices[dof])))
                          {
                            constraints.add_line(face_dof_indices[dof]);

                            if (std::abs(dof_values[dof]) > 1e-13)
                              constraints.set_inhomogeneity(
                                face_dof_indices[dof], dof_values[dof]);
                          }
                    }

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  }



  template <int dim>
  void

  project_boundary_values_curl_conforming(
    const hp::DoFHandler<dim> &       dof_handler,
    const unsigned int                first_vector_component,
    const Function<dim> &             boundary_function,
    const types::boundary_id          boundary_component,
    AffineConstraints<double> &       constraints,
    const hp::MappingCollection<dim> &mapping_collection)
  {
    const hp::FECollection<dim> &fe_collection(dof_handler.get_fe_collection());
    hp::QCollection<dim>         face_quadrature_collection;

    for (unsigned int i = 0; i < fe_collection.size(); ++i)
      {
        const QGauss<dim - 1> reference_face_quadrature(
          2 * fe_collection[i].degree);

        for (unsigned int face : GeometryInfo<dim>::face_indices())
          face_quadrature_collection.push_back(
            QProjector<dim>::project_to_face(reference_face_quadrature, face));
      }

    hp::FEValues<dim>                    fe_face_values(mapping_collection,
                                     fe_collection,
                                     face_quadrature_collection,
                                     update_jacobians | update_JxW_values |
                                       update_quadrature_points |
                                       update_values);
    std::vector<bool>                    dofs_processed;
    std::vector<double>                  dof_values;
    std::vector<types::global_dof_index> face_dof_indices;
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active();

    switch (dim)
      {
        case 2:
          {
            for (; cell != dof_handler.end(); ++cell)
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // if the FE is a FE_Nothing object there is no work to do
                      if (dynamic_cast<const FE_Nothing<dim> *>(
                            &cell->get_fe()) != nullptr)
                        return;

                      // This is only implemented, if the FE is a Nedelec
                      // element. If the FE is a FESystem we cannot check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(dynamic_cast<const FE_Nedelec<dim> *>(
                                        &cell->get_fe()) != nullptr,
                                      typename FiniteElement<
                                        dim>::ExcInterpolationNotImplemented());
                        }

                      const unsigned int dofs_per_face =
                        cell->get_fe().dofs_per_face;

                      dofs_processed.resize(dofs_per_face);
                      dof_values.resize(dofs_per_face);

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        {
                          dof_values[dof]     = 0.0;
                          dofs_processed[dof] = false;
                        }

                      internals::compute_face_projection_curl_conforming(
                        cell,
                        face,
                        fe_face_values,
                        boundary_function,
                        first_vector_component,
                        dof_values,
                        dofs_processed);
                      face_dof_indices.resize(dofs_per_face);
                      cell->face(face)->get_dof_indices(
                        face_dof_indices, cell->active_fe_index());

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        if (dofs_processed[dof] &&
                            constraints.can_store_line(face_dof_indices[dof]) &&
                            !(constraints.is_constrained(
                              face_dof_indices[dof])))
                          {
                            constraints.add_line(face_dof_indices[dof]);

                            if (std::abs(dof_values[dof]) > 1e-13)
                              constraints.set_inhomogeneity(
                                face_dof_indices[dof], dof_values[dof]);
                          }
                    }

            break;
          }

        case 3:
          {
            hp::QCollection<dim> edge_quadrature_collection;

            for (unsigned int i = 0; i < fe_collection.size(); ++i)
              {
                const QGauss<dim - 2> reference_edge_quadrature(
                  2 * fe_collection[i].degree);

                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  for (unsigned int line = 0;
                       line < GeometryInfo<dim>::lines_per_face;
                       ++line)
                    edge_quadrature_collection.push_back(
                      QProjector<dim>::project_to_face(
                        QProjector<dim - 1>::project_to_face(
                          reference_edge_quadrature, line),
                        face));
              }

            hp::FEValues<dim> fe_edge_values(mapping_collection,
                                             fe_collection,
                                             edge_quadrature_collection,
                                             update_jacobians |
                                               update_JxW_values |
                                               update_quadrature_points |
                                               update_values);

            for (; cell != dof_handler.end(); ++cell)
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // if the FE is a FE_Nothing object there is no work to do
                      if (dynamic_cast<const FE_Nothing<dim> *>(
                            &cell->get_fe()) != nullptr)
                        return;

                      // This is only implemented, if the FE is a Nedelec
                      // element. If the FE is a FESystem we cannot check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(dynamic_cast<const FE_Nedelec<dim> *>(
                                        &cell->get_fe()) != nullptr,
                                      typename FiniteElement<
                                        dim>::ExcInterpolationNotImplemented());
                        }

                      const unsigned int superdegree = cell->get_fe().degree;
                      const unsigned int degree      = superdegree - 1;
                      const unsigned int dofs_per_face =
                        cell->get_fe().dofs_per_face;

                      dofs_processed.resize(dofs_per_face);
                      dof_values.resize(dofs_per_face);

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        {
                          dof_values[dof]     = 0.0;
                          dofs_processed[dof] = false;
                        }

                      for (unsigned int line = 0;
                           line < GeometryInfo<dim>::lines_per_face;
                           ++line)
                        internals::compute_edge_projection(
                          cell,
                          face,
                          line,
                          fe_edge_values,
                          boundary_function,
                          first_vector_component,
                          dof_values,
                          dofs_processed);

                      // If there are higher order shape functions, there is
                      // still some work left.
                      if (degree > 0)
                        internals::compute_face_projection_curl_conforming(
                          cell,
                          face,
                          fe_face_values,
                          boundary_function,
                          first_vector_component,
                          dof_values,
                          dofs_processed);


                      face_dof_indices.resize(dofs_per_face);
                      cell->face(face)->get_dof_indices(
                        face_dof_indices, cell->active_fe_index());

                      for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                        if (dofs_processed[dof] &&
                            constraints.can_store_line(face_dof_indices[dof]) &&
                            !(constraints.is_constrained(
                              face_dof_indices[dof])))
                          {
                            constraints.add_line(face_dof_indices[dof]);

                            if (std::abs(dof_values[dof]) > 1e-13)
                              constraints.set_inhomogeneity(
                                face_dof_indices[dof], dof_values[dof]);
                          }
                    }

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  }


  namespace internals
  {
    template <int dim, typename cell_iterator, typename number>
    typename std::enable_if<dim == 3>::type
    compute_edge_projection_l2(const cell_iterator &        cell,
                               const unsigned int           face,
                               const unsigned int           line,
                               hp::FEValues<dim> &          hp_fe_values,
                               const Function<dim, number> &boundary_function,
                               const unsigned int   first_vector_component,
                               std::vector<number> &dof_values,
                               std::vector<bool> &  dofs_processed)
    {
      // This function computes the L2-projection of the given
      // boundary function on 3D edges and returns the constraints
      // associated with the edge functions for the given cell.
      //
      // In the context of this function, by associated DoFs we mean:
      // the DoFs corresponding to the group of components making up the vector
      // with first component first_vector_component (length dim).
      const FiniteElement<dim> &fe = cell->get_fe();

      // reinit for this cell, face and line.
      hp_fe_values.reinit(
        cell,
        (cell->active_fe_index() * GeometryInfo<dim>::faces_per_cell + face) *
            GeometryInfo<dim>::lines_per_face +
          line);

      // Initialize the required objects.
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

      const std::vector<Point<dim>> &quadrature_points =
        fe_values.get_quadrature_points();
      std::vector<Vector<number>> values(fe_values.n_quadrature_points,
                                         Vector<number>(fe.n_components()));

      // Get boundary function values
      // at quadrature points.
      boundary_function.vector_value_list(quadrature_points, values);

      // Find the group of vector components we want to project onto
      // (dim of them, starting at first_vector_component) within the
      // overall finite element (which may be an FESystem).
      std::pair<unsigned int, unsigned int> base_indices(0, 0);
      if (dynamic_cast<const FESystem<dim> *>(&cell->get_fe()) != nullptr)
        {
          unsigned int fe_index     = 0;
          unsigned int fe_index_old = 0;
          unsigned int i            = 0;

          // Find base element:
          // base_indices.first
          //
          // Then select which copy of that base element
          // [ each copy is of length
          // fe.base_element(base_indices.first).n_components() ] corresponds to
          // first_vector_component: base_index.second
          for (; i < fe.n_base_elements(); ++i)
            {
              fe_index_old = fe_index;
              fe_index +=
                fe.element_multiplicity(i) * fe.base_element(i).n_components();

              if (fe_index > first_vector_component)
                break;
            }

          base_indices.first  = i;
          base_indices.second = (first_vector_component - fe_index_old) /
                                fe.base_element(i).n_components();
        }
      else
        // The only other element we know how to deal with (so far) is
        // FE_Nedelec, which has one base element and one copy of it
        // (with 3 components). In that case, the values of
        // 'base_indices' as initialized above are correct.
        Assert((dynamic_cast<const FE_Nedelec<dim> *>(&cell->get_fe()) !=
                nullptr) ||
                 (dynamic_cast<const FE_NedelecSZ<dim> *>(&cell->get_fe()) !=
                  nullptr),
               ExcNotImplemented());


      // Store the 'degree' of the Nedelec element as fe.degree-1. For
      // Nedelec elements, FE_Nedelec<dim>(0) returns fe.degree = 1
      // because fe.degree stores the *polynomial* degree, not the
      // degree of the element (which is typically defined based on
      // the largest polynomial space that is *complete* within the
      // finite element).
      const unsigned int degree =
        fe.base_element(base_indices.first).degree - 1;

      // Find DoFs we want to constrain: There are
      // fe.base_element(base_indices.first).dofs_per_line DoFs
      // associated with the given line on the given face on the given
      // cell.
      //
      // We need to know which of these DoFs (there are degree+1 of interest)
      // are associated with the components given by first_vector_component.
      // Then we can make a map from the associated line DoFs to the face DoFs.
      //
      // For a single FE_Nedelec<3> element this is simple:
      //    We know the ordering of local DoFs goes
      //    lines -> faces -> cells
      //
      // For a set of FESystem<3> elements we need to pick out the matching base
      // element and the index within this ordering.
      //
      // We call the map associated_edge_dof_to_face_dof
      std::vector<unsigned int> associated_edge_dof_to_face_dof(
        degree + 1, numbers::invalid_unsigned_int);

      // Lowest DoF in the base element allowed for this edge:
      const unsigned int lower_bound =
        fe.base_element(base_indices.first)
          .face_to_cell_index(line * (degree + 1), face);
      // Highest DoF in the base element allowed for this edge:
      const unsigned int upper_bound =
        fe.base_element(base_indices.first)
          .face_to_cell_index((line + 1) * (degree + 1) - 1, face);

      unsigned int associated_edge_dof_index = 0;
      for (unsigned int line_dof_idx = 0; line_dof_idx < fe.dofs_per_line;
           ++line_dof_idx)
        {
          // For each DoF associated with the (interior of) the line, we need
          // to figure out which base element it belongs to and then if
          // that's the correct base element. This is complicated by the
          // fact that the FiniteElement class has functions that translate
          // from face to cell, but not from edge to cell index systems. So
          // we have to do that step by step.
          //
          // DoFs on a face in 3d are numbered in order by vertices then lines
          // then faces.
          // i.e. line 0 has degree+1 dofs numbered 0,..,degree
          //      line 1 has degree+1 dofs numbered (degree+1),..,2*(degree+1)
          //      and so on.

          const unsigned int face_dof_idx =
            GeometryInfo<dim>::vertices_per_face * fe.dofs_per_vertex +
            line * fe.dofs_per_line + line_dof_idx;

          // Note, assuming that the edge orientations are "standard"
          //       i.e. cell->line_orientation(line) = true.
          Assert(cell->line_orientation(line),
                 ExcMessage("Edge orientation does not meet expectation."));
          // Next, translate from face to cell. Note, this might be assuming
          // that the edge orientations are "standard" (not sure any more at
          // this time), i.e.
          //       cell->line_orientation(line) = true.
          const unsigned int cell_dof_idx =
            fe.face_to_cell_index(face_dof_idx, face);

          // Check that this cell_idx belongs to the correct base_element,
          // component and line. We do this for each of the supported elements
          // separately
          bool dof_is_of_interest = false;
          if (dynamic_cast<const FESystem<dim> *>(&fe) != nullptr)
            {
              dof_is_of_interest =
                (fe.system_to_base_index(cell_dof_idx).first == base_indices) &&
                (lower_bound <= fe.system_to_base_index(cell_dof_idx).second) &&
                (fe.system_to_base_index(cell_dof_idx).second <= upper_bound);
            }
          else if ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr) ||
                   (dynamic_cast<const FE_NedelecSZ<dim> *>(&fe) != nullptr))
            {
              Assert((line * (degree + 1) <= face_dof_idx) &&
                       (face_dof_idx < (line + 1) * (degree + 1)),
                     ExcInternalError());
              dof_is_of_interest = true;
            }
          else
            Assert(false, ExcNotImplemented());

          if (dof_is_of_interest)
            {
              associated_edge_dof_to_face_dof[associated_edge_dof_index] =
                face_dof_idx;
              ++associated_edge_dof_index;
            }
        }
      // Sanity check:
      const unsigned int n_associated_edge_dofs = associated_edge_dof_index;
      Assert(n_associated_edge_dofs == degree + 1,
             ExcMessage("Error: Unexpected number of 3D edge DoFs"));

      // Matrix and RHS vectors to store linear system:
      // We have (degree+1) basis functions for an edge
      FullMatrix<number> edge_matrix(degree + 1, degree + 1);
      FullMatrix<number> edge_matrix_inv(degree + 1, degree + 1);
      Vector<number>     edge_rhs(degree + 1);
      Vector<number>     edge_solution(degree + 1);

      const FEValuesExtractors::Vector vec(first_vector_component);

      // coordinate directions of
      // the edges of the face.
      const unsigned int
        edge_coordinate_direction[GeometryInfo<dim>::faces_per_cell]
                                 [GeometryInfo<dim>::lines_per_face] = {
                                   {2, 2, 1, 1},
                                   {2, 2, 1, 1},
                                   {0, 0, 2, 2},
                                   {0, 0, 2, 2},
                                   {1, 1, 0, 0},
                                   {1, 1, 0, 0}};

      const double tol =
        0.5 * cell->face(face)->line(line)->diameter() / fe.degree;
      const std::vector<Point<dim>> &reference_quadrature_points =
        fe_values.get_quadrature().get_points();

      // Project the boundary function onto the shape functions for this edge
      // and set up a linear system of equations to get the values for the DoFs
      // associated with this edge.
      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
           ++q_point)
        {
          // Compute the tangential
          // of the edge at
          // the quadrature point.
          Point<dim> shifted_reference_point_1 =
            reference_quadrature_points[q_point];
          Point<dim> shifted_reference_point_2 =
            reference_quadrature_points[q_point];

          shifted_reference_point_1(edge_coordinate_direction[face][line]) +=
            tol;
          shifted_reference_point_2(edge_coordinate_direction[face][line]) -=
            tol;
          Tensor<1, dim> tangential =
            (0.5 *
             (fe_values.get_mapping().transform_unit_to_real_cell(
                cell, shifted_reference_point_1) -
              fe_values.get_mapping().transform_unit_to_real_cell(
                cell, shifted_reference_point_2)) /
             tol);
          tangential /= tangential.norm();

          // Compute the entries of the linear system
          // Note the system is symmetric so we could only compute the
          // lower/upper triangle.
          //
          // The matrix entries are
          // \int_{edge}
          // (tangential*edge_shape_function_i)*(tangential*edge_shape_function_j)
          // dS
          //
          // The RHS entries are:
          // \int_{edge}
          // (tangential*boundary_value)*(tangential*edge_shape_function_i) dS.
          for (unsigned int j = 0; j < n_associated_edge_dofs; ++j)
            {
              const unsigned int j_face_idx =
                associated_edge_dof_to_face_dof[j];
              const unsigned int j_cell_idx =
                fe.face_to_cell_index(j_face_idx, face);
              for (unsigned int i = 0; i < n_associated_edge_dofs; ++i)
                {
                  const unsigned int i_face_idx =
                    associated_edge_dof_to_face_dof[i];
                  const unsigned int i_cell_idx =
                    fe.face_to_cell_index(i_face_idx, face);

                  edge_matrix(i, j) +=
                    fe_values.JxW(q_point) *
                    (fe_values[vec].value(i_cell_idx, q_point) * tangential) *
                    (fe_values[vec].value(j_cell_idx, q_point) * tangential);
                }
              // Compute the RHS entries:
              edge_rhs(j) +=
                fe_values.JxW(q_point) *
                (values[q_point](first_vector_component) * tangential[0] +
                 values[q_point](first_vector_component + 1) * tangential[1] +
                 values[q_point](first_vector_component + 2) * tangential[2]) *
                (fe_values[vec].value(j_cell_idx, q_point) * tangential);
            }
        }

      // Invert linear system
      edge_matrix_inv.invert(edge_matrix);
      edge_matrix_inv.vmult(edge_solution, edge_rhs);

      // Store computed DoFs
      for (unsigned int i = 0; i < n_associated_edge_dofs; ++i)
        {
          dof_values[associated_edge_dof_to_face_dof[i]]     = edge_solution(i);
          dofs_processed[associated_edge_dof_to_face_dof[i]] = true;
        }
    }


    template <int dim, typename cell_iterator, typename number>
    typename std::enable_if<dim != 3>::type
    compute_edge_projection_l2(const cell_iterator &,
                               const unsigned int,
                               const unsigned int,
                               hp::FEValues<dim> &,
                               const Function<dim, number> &,
                               const unsigned int,
                               std::vector<number> &,
                               std::vector<bool> &)
    {
      // dummy implementation of above function
      // for all other dimensions
      Assert(false, ExcInternalError());
    }


    template <int dim, typename cell_iterator, typename number>
    void
    compute_face_projection_curl_conforming_l2(
      const cell_iterator &        cell,
      const unsigned int           face,
      hp::FEFaceValues<dim> &      hp_fe_face_values,
      const Function<dim, number> &boundary_function,
      const unsigned int           first_vector_component,
      std::vector<number> &        dof_values,
      std::vector<bool> &          dofs_processed)
    {
      // This function computes the L2-projection of the boundary
      // function on the interior of faces only. In 3D, this should only be
      // called after first calling compute_edge_projection_l2, as it relies on
      // edge constraints which are found.

      // In the context of this function, by associated DoFs we mean:
      // the DoFs corresponding to the group of components making up the vector
      // with first component first_vector_component (with total components
      // dim).

      // Copy to the standard FEFaceValues object:
      hp_fe_face_values.reinit(cell, face);
      const FEFaceValues<dim> &fe_face_values =
        hp_fe_face_values.get_present_fe_values();

      // Initialize the required objects.
      const FiniteElement<dim> &     fe = cell->get_fe();
      const std::vector<Point<dim>> &quadrature_points =
        fe_face_values.get_quadrature_points();

      std::vector<Vector<number>> values(fe_face_values.n_quadrature_points,
                                         Vector<number>(fe.n_components()));

      // Get boundary function values at quadrature points.
      boundary_function.vector_value_list(quadrature_points, values);

      // Find where the group of vector components (dim of them,
      // starting at first_vector_component) are within an FESystem.
      //
      // If not using FESystem then must be using FE_Nedelec,
      // which has one base element and one copy of it (with 3 components).
      std::pair<unsigned int, unsigned int> base_indices(0, 0);
      if (dynamic_cast<const FESystem<dim> *>(&cell->get_fe()) != nullptr)
        {
          unsigned int fe_index     = 0;
          unsigned int fe_index_old = 0;
          unsigned int i            = 0;

          // Find base element:
          // base_indices.first
          //
          // Then select which copy of that base element
          // [ each copy is of length
          // fe.base_element(base_indices.first).n_components() ] corresponds to
          // first_vector_component: base_index.second
          for (; i < fe.n_base_elements(); ++i)
            {
              fe_index_old = fe_index;
              fe_index +=
                fe.element_multiplicity(i) * fe.base_element(i).n_components();

              if (fe_index > first_vector_component)
                break;
            }
          base_indices.first  = i;
          base_indices.second = (first_vector_component - fe_index_old) /
                                fe.base_element(i).n_components();
        }
      else
        {
          // Assert that the FE is in fact an FE_Nedelec, so that the default
          // base_indices == (0,0) is correct.
          Assert((dynamic_cast<const FE_Nedelec<dim> *>(&cell->get_fe()) !=
                  nullptr) ||
                   (dynamic_cast<const FE_NedelecSZ<dim> *>(&cell->get_fe()) !=
                    nullptr),
                 ExcNotImplemented());
        }
      const unsigned int degree =
        fe.base_element(base_indices.first).degree - 1;

      switch (dim)
        {
          case 2:
            // NOTE: This is very similar to compute_edge_projection as used in
            // 3D,
            //       and contains a lot of overlap with that function.
            {
              // Find the DoFs we want to constrain. There are degree+1 in
              // total. Create a map from these to the face index Note:
              //    - for a single FE_Nedelec<2> element this is
              //      simply 0 to fe.dofs_per_face
              //    - for FESystem<2> this just requires matching the
              //      base element, fe.system_to_base_index.first.first
              //      and the copy of the base element we're interested
              //      in, fe.system_to_base_index.first.second
              std::vector<unsigned int> associated_edge_dof_to_face_dof(degree +
                                                                        1);

              unsigned int associated_edge_dof_index = 0;
              for (unsigned int face_idx = 0; face_idx < fe.dofs_per_face;
                   ++face_idx)
                {
                  const unsigned int cell_idx =
                    fe.face_to_cell_index(face_idx, face);
                  if (((dynamic_cast<const FESystem<dim> *>(&fe) != nullptr) &&
                       (fe.system_to_base_index(cell_idx).first ==
                        base_indices)) ||
                      (dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr) ||
                      (dynamic_cast<const FE_NedelecSZ<dim> *>(&fe) != nullptr))
                    {
                      associated_edge_dof_to_face_dof
                        [associated_edge_dof_index] = face_idx;
                      ++associated_edge_dof_index;
                    }
                }
              // Sanity check:
              const unsigned int associated_edge_dofs =
                associated_edge_dof_index;
              Assert(associated_edge_dofs == degree + 1,
                     ExcMessage("Error: Unexpected number of 2D edge DoFs"));

              // Matrix and RHS vectors to store:
              // We have (degree+1) edge basis functions
              FullMatrix<number> edge_matrix(degree + 1, degree + 1);
              FullMatrix<number> edge_matrix_inv(degree + 1, degree + 1);
              Vector<number>     edge_rhs(degree + 1);
              Vector<number>     edge_solution(degree + 1);

              const FEValuesExtractors::Vector vec(first_vector_component);

              // Project the boundary function onto the shape functions for this
              // edge and set up a linear system of equations to get the values
              // for the DoFs associated with this edge.
              for (unsigned int q_point = 0;
                   q_point < fe_face_values.n_quadrature_points;
                   ++q_point)
                {
                  // Compute the entries of the linear system
                  // Note the system is symmetric so we could only compute the
                  // lower/upper triangle.
                  //
                  // The matrix entries are
                  // \int_{edge} (tangential * edge_shape_function_i) *
                  // (tangential * edge_shape_function_j) dS
                  //
                  // The RHS entries are:
                  // \int_{edge} (tangential* boundary_value) * (tangential *
                  // edge_shape_function_i) dS.
                  //
                  // In 2D, tangential*vector is equivalent to
                  // cross_product_3d(normal, vector), so we use this instead.
                  // This avoids possible issues with the computation of the
                  // tangent.

                  // Store the normal at this quad point:
                  Tensor<1, dim> normal_at_q_point =
                    fe_face_values.normal_vector(q_point);
                  for (unsigned int j = 0; j < associated_edge_dofs; ++j)
                    {
                      const unsigned int j_face_idx =
                        associated_edge_dof_to_face_dof[j];
                      const unsigned int j_cell_idx =
                        fe.face_to_cell_index(j_face_idx, face);

                      Tensor<1, dim> phi_j =
                        fe_face_values[vec].value(j_cell_idx, q_point);
                      for (unsigned int i = 0; i < associated_edge_dofs; ++i)
                        {
                          const unsigned int i_face_idx =
                            associated_edge_dof_to_face_dof[i];
                          const unsigned int i_cell_idx =
                            fe.face_to_cell_index(i_face_idx, face);

                          Tensor<1, dim> phi_i =
                            fe_face_values[vec].value(i_cell_idx, q_point);

                          // Using n cross phi
                          edge_matrix(i, j) +=
                            fe_face_values.JxW(q_point) *
                            ((phi_i[1] * normal_at_q_point[0] -
                              phi_i[0] * normal_at_q_point[1]) *
                             (phi_j[1] * normal_at_q_point[0] -
                              phi_j[0] * normal_at_q_point[1]));
                        }
                      // Using n cross phi
                      edge_rhs(j) +=
                        fe_face_values.JxW(q_point) *
                        ((values[q_point](first_vector_component + 1) *
                            normal_at_q_point[0] -
                          values[q_point](first_vector_component) *
                            normal_at_q_point[1]) *
                         (phi_j[1] * normal_at_q_point[0] -
                          phi_j[0] * normal_at_q_point[1]));
                    }
                }

              // Invert linear system
              edge_matrix_inv.invert(edge_matrix);
              edge_matrix_inv.vmult(edge_solution, edge_rhs);

              // Store computed DoFs
              for (unsigned int associated_edge_dof_index = 0;
                   associated_edge_dof_index < associated_edge_dofs;
                   ++associated_edge_dof_index)
                {
                  dof_values[associated_edge_dof_to_face_dof
                               [associated_edge_dof_index]] =
                    edge_solution(associated_edge_dof_index);
                  dofs_processed[associated_edge_dof_to_face_dof
                                   [associated_edge_dof_index]] = true;
                }
              break;
            }

          case 3:
            {
              const FEValuesExtractors::Vector vec(first_vector_component);

              // First group DoFs associated with edges which we already know.
              // Sort these into groups of dofs (0 -> degree+1 of them) by each
              // edge. This will help when computing the residual for the face
              // projections.
              //
              // This matches with the search done in compute_edge_projection.
              const unsigned int lines_per_face =
                GeometryInfo<dim>::lines_per_face;
              std::vector<std::vector<unsigned int>>
                                        associated_edge_dof_to_face_dof(lines_per_face,
                                                                        std::vector<unsigned int>(degree +
                                                                        1));
              std::vector<unsigned int> associated_edge_dofs(lines_per_face);

              for (unsigned int line = 0; line < lines_per_face; ++line)
                {
                  // Lowest DoF in the base element allowed for this edge:
                  const unsigned int lower_bound =
                    fe.base_element(base_indices.first)
                      .face_to_cell_index(line * (degree + 1), face);
                  // Highest DoF in the base element allowed for this edge:
                  const unsigned int upper_bound =
                    fe.base_element(base_indices.first)
                      .face_to_cell_index((line + 1) * (degree + 1) - 1, face);
                  unsigned int associated_edge_dof_index = 0;

                  for (unsigned int line_dof_idx = 0;
                       line_dof_idx < fe.dofs_per_line;
                       ++line_dof_idx)
                    {
                      // For each DoF associated with the (interior of) the
                      // line, we need to figure out which base element it
                      // belongs to and then if that's the correct base element.
                      // This is complicated by the fact that the FiniteElement
                      // class has functions that translate from face to cell,
                      // but not from edge to cell index systems. So we have to
                      // do that step by step.
                      //
                      // DoFs on a face in 3d are numbered in order by vertices
                      // then lines then faces. i.e. line 0 has degree+1 dofs
                      // numbered 0,..,degree
                      //      line 1 has degree+1 dofs numbered
                      //      (degree+1),..,2*(degree+1) and so on.
                      const unsigned int face_dof_idx =
                        GeometryInfo<dim>::vertices_per_face *
                          fe.dofs_per_vertex +
                        line * fe.dofs_per_line + line_dof_idx;

                      // Next, translate from face to cell. Note, this might be
                      // assuming that the edge orientations are "standard" (not
                      // sure any more at this time), i.e.
                      //       cell->line_orientation(line) = true.
                      const unsigned int cell_dof_idx =
                        fe.face_to_cell_index(face_dof_idx, face);

                      // Check that this cell_idx belongs to the correct
                      // base_element, component and line. We do this for each
                      // of the supported elements separately
                      bool dof_is_of_interest = false;
                      if (dynamic_cast<const FESystem<dim> *>(&fe) != nullptr)
                        {
                          dof_is_of_interest =
                            (fe.system_to_base_index(cell_dof_idx).first ==
                             base_indices) &&
                            (lower_bound <=
                             fe.system_to_base_index(cell_dof_idx).second) &&
                            (fe.system_to_base_index(cell_dof_idx).second <=
                             upper_bound);
                        }
                      else if ((dynamic_cast<const FE_Nedelec<dim> *>(&fe) !=
                                nullptr) ||
                               (dynamic_cast<const FE_NedelecSZ<dim> *>(&fe) !=
                                nullptr))
                        {
                          Assert((line * (degree + 1) <= face_dof_idx) &&
                                   (face_dof_idx < (line + 1) * (degree + 1)),
                                 ExcInternalError());
                          dof_is_of_interest = true;
                        }
                      else
                        Assert(false, ExcNotImplemented());

                      if (dof_is_of_interest)
                        {
                          associated_edge_dof_to_face_dof
                            [line][associated_edge_dof_index] = face_dof_idx;
                          ++associated_edge_dof_index;
                        }
                    }
                  // Sanity check:
                  associated_edge_dofs[line] = associated_edge_dof_index;
                  Assert(associated_edge_dofs[line] == degree + 1,
                         ExcInternalError());
                }

              // Next find the face DoFs associated with the vector components
              // we're interested in. There are 2*degree*(degree+1) DoFs
              // associated with the interior of each face (not including
              // edges!).
              //
              // Create a map mapping from the consecutively numbered
              // associated_dofs to the face DoF (which can be transferred to a
              // local cell index).
              //
              // For FE_Nedelec<3> we just need to have a face numbering greater
              // than the number of edge DoFs (=lines_per_face*(degree+1).
              //
              // For FESystem<3> we need to base the base_indices (base element
              // and copy within that base element) and ensure we're above the
              // number of edge DoFs within that base element.
              std::vector<unsigned int> associated_face_dof_to_face_dof(
                2 * degree * (degree + 1));

              // Loop over these quad-interior dofs.
              unsigned int associated_face_dof_index = 0;
              for (unsigned int quad_dof_idx = 0;
                   quad_dof_idx < fe.dofs_per_quad;
                   ++quad_dof_idx)
                {
                  const unsigned int face_idx =
                    GeometryInfo<dim>::vertices_per_face * fe.dofs_per_vertex +
                    lines_per_face * fe.dofs_per_line + quad_dof_idx;
                  const unsigned int cell_idx =
                    fe.face_to_cell_index(face_idx, face);
                  if (((dynamic_cast<const FESystem<dim> *>(&fe) != nullptr) &&
                       (fe.system_to_base_index(cell_idx).first ==
                        base_indices)) ||
                      (dynamic_cast<const FE_Nedelec<dim> *>(&fe) != nullptr) ||
                      (dynamic_cast<const FE_NedelecSZ<dim> *>(&fe) != nullptr))
                    {
                      AssertIndexRange(associated_face_dof_index,
                                       associated_face_dof_to_face_dof.size());
                      associated_face_dof_to_face_dof
                        [associated_face_dof_index] = face_idx;
                      ++associated_face_dof_index;
                    }
                }
              // Sanity check:
              const unsigned int associated_face_dofs =
                associated_face_dof_index;
              Assert(associated_face_dofs == 2 * degree * (degree + 1),
                     ExcMessage("Error: Unexpected number of 3D face DoFs"));

              // Storage for the linear system.
              // There are 2*degree*(degree+1) DoFs associated with a face in
              // 3D. Note this doesn't include the DoFs associated with edges on
              // that face.
              FullMatrix<number> face_matrix(2 * degree * (degree + 1));
              FullMatrix<number> face_matrix_inv(2 * degree * (degree + 1));
              Vector<number>     face_rhs(2 * degree * (degree + 1));
              Vector<number>     face_solution(2 * degree * (degree + 1));

              // Project the boundary function onto the shape functions for this
              // face and set up a linear system of equations to get the values
              // for the DoFs associated with this face. We also must include
              // the residuals from the shape functions associated with edges.
              Tensor<1, dim, number> tmp;
              Tensor<1, dim>         cross_product_i;
              Tensor<1, dim>         cross_product_j;
              Tensor<1, dim, number> cross_product_rhs;

              // Loop to construct face linear system.
              for (unsigned int q_point = 0;
                   q_point < fe_face_values.n_quadrature_points;
                   ++q_point)
                {
                  // First calculate the residual from the edge functions
                  // store the result in tmp.
                  //
                  // Edge_residual =
                  //        boundary_value - (
                  //            \sum_(edges on face)
                  //                 \sum_(DoFs on edge)
                  //                 edge_dof_value*edge_shape_function
                  //                   )
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      tmp[d] = 0.0;
                    }
                  for (unsigned int line = 0; line < lines_per_face; ++line)
                    {
                      for (unsigned int associated_edge_dof = 0;
                           associated_edge_dof < associated_edge_dofs[line];
                           ++associated_edge_dof)
                        {
                          const unsigned int face_idx =
                            associated_edge_dof_to_face_dof
                              [line][associated_edge_dof];
                          const unsigned int cell_idx =
                            fe.face_to_cell_index(face_idx, face);
                          tmp -= dof_values[face_idx] *
                                 fe_face_values[vec].value(cell_idx, q_point);
                        }
                    }

                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      tmp[d] += values[q_point](first_vector_component + d);
                    }

                  // Tensor of normal vector on the face at q_point;
                  const Tensor<1, dim> normal_vector =
                    fe_face_values.normal_vector(q_point);

                  // Now compute the linear system:
                  // On a face:
                  // The matrix entries are:
                  // \int_{face} (n x face_shape_function_i) \cdot ( n x
                  // face_shape_function_j) dS
                  //
                  // The RHS entries are:
                  // \int_{face} (n x (Edge_residual) \cdot (n x
                  // face_shape_function_i) dS

                  for (unsigned int j = 0; j < associated_face_dofs; ++j)
                    {
                      const unsigned int j_face_idx =
                        associated_face_dof_to_face_dof[j];
                      const unsigned int cell_j =
                        fe.face_to_cell_index(j_face_idx, face);

                      cross_product_j =
                        cross_product_3d(normal_vector,
                                         fe_face_values[vec].value(cell_j,
                                                                   q_point));

                      for (unsigned int i = 0; i < associated_face_dofs; ++i)
                        {
                          const unsigned int i_face_idx =
                            associated_face_dof_to_face_dof[i];
                          const unsigned int cell_i =
                            fe.face_to_cell_index(i_face_idx, face);
                          cross_product_i = cross_product_3d(
                            normal_vector,
                            fe_face_values[vec].value(cell_i, q_point));

                          face_matrix(i, j) += fe_face_values.JxW(q_point) *
                                               cross_product_i *
                                               cross_product_j;
                        }
                      // compute rhs
                      cross_product_rhs = cross_product_3d(normal_vector, tmp);
                      face_rhs(j) += fe_face_values.JxW(q_point) *
                                     cross_product_rhs * cross_product_j;
                    }
                }

              // Solve linear system:
              if (associated_face_dofs > 0)
                {
                  face_matrix_inv.invert(face_matrix);
                  face_matrix_inv.vmult(face_solution, face_rhs);
                }

              // Store computed DoFs:
              for (unsigned int associated_face_dof = 0;
                   associated_face_dof < associated_face_dofs;
                   ++associated_face_dof)
                {
                  dof_values
                    [associated_face_dof_to_face_dof[associated_face_dof]] =
                      face_solution(associated_face_dof);
                  dofs_processed
                    [associated_face_dof_to_face_dof[associated_face_dof]] =
                      true;
                }
              break;
            }
          default:
            Assert(false, ExcNotImplemented());
        }
    }


    template <int dim, typename DoFHandlerType, typename number>
    void
    compute_project_boundary_values_curl_conforming_l2(
      const DoFHandlerType &                 dof_handler,
      const unsigned int                     first_vector_component,
      const Function<dim, number> &          boundary_function,
      const types::boundary_id               boundary_component,
      AffineConstraints<number> &            constraints,
      const hp::MappingCollection<dim, dim> &mapping_collection)
    {
      // L2-projection based interpolation formed in one (in 2D) or two (in 3D)
      // steps.
      //
      // In 2D we only need to constrain edge DoFs.
      //
      // In 3D we need to constrain both edge and face DoFs. This is done in two
      // parts.
      //
      // For edges, since the face shape functions are zero here ("bubble
      // functions"), we project the tangential component of the boundary
      // function and compute the L2-projection. This returns the values for the
      // DoFs associated with each edge shape function. In 3D, this is computed
      // by internals::compute_edge_projection_l2, in 2D, it is handled by
      // compute_face_projection_curl_conforming_l2.
      //
      // For faces we compute the residual of the boundary function which is
      // satisfied by the edge shape functions alone. Which can then be used to
      // calculate the remaining face DoF values via a projection which leads to
      // a linear system to solve. This is handled by
      // compute_face_projection_curl_conforming_l2
      //
      // For details see (for example) section 4.2:
      // Electromagnetic scattering simulation using an H (curl) conforming hp
      // finite element method in three dimensions, PD Ledger, K Morgan, O
      // Hassan, Int. J.  Num. Meth. Fluids, Volume 53, Issue 8, pages
      // 1267-1296, 20 March 2007:
      // http://onlinelibrary.wiley.com/doi/10.1002/fld.1223/abstract

      // Create hp FEcollection, dof_handler can be either hp or standard type.
      // From here on we can treat it like a hp-namespace object.
      const hp::FECollection<dim> &fe_collection(
        dof_handler.get_fe_collection());

      // Create face quadrature collection
      hp::QCollection<dim - 1> face_quadrature_collection;
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        {
          const QGauss<dim - 1> reference_face_quadrature(
            2 * fe_collection[i].degree + 1);
          face_quadrature_collection.push_back(reference_face_quadrature);
        }

      hp::FEFaceValues<dim> fe_face_values(mapping_collection,
                                           fe_collection,
                                           face_quadrature_collection,
                                           update_values |
                                             update_quadrature_points |
                                             update_normal_vectors |
                                             update_JxW_values);

      // Storage for dof values found and whether they have been processed:
      std::vector<bool>                             dofs_processed;
      std::vector<number>                           dof_values;
      std::vector<types::global_dof_index>          face_dof_indices;
      typename DoFHandlerType::active_cell_iterator cell =
        dof_handler.begin_active();

      switch (dim)
        {
          case 2:
            {
              for (; cell != dof_handler.end(); ++cell)
                {
                  if (cell->at_boundary() && cell->is_locally_owned())
                    {
                      for (const unsigned int face :
                           GeometryInfo<dim>::face_indices())
                        {
                          if (cell->face(face)->boundary_id() ==
                              boundary_component)
                            {
                              // If the FE is an FE_Nothing object there is no
                              // work to do
                              if (dynamic_cast<const FE_Nothing<dim> *>(
                                    &cell->get_fe()) != nullptr)
                                {
                                  return;
                                }

                              // This is only implemented for FE_Nedelec
                              // elements. If the FE is a FESystem we cannot
                              // check this.
                              if (dynamic_cast<const FESystem<dim> *>(
                                    &cell->get_fe()) == nullptr)
                                {
                                  AssertThrow(
                                    (dynamic_cast<const FE_Nedelec<dim> *>(
                                       &cell->get_fe()) != nullptr) ||
                                      (dynamic_cast<const FE_NedelecSZ<dim> *>(
                                         &cell->get_fe()) != nullptr),
                                    typename FiniteElement<
                                      dim>::ExcInterpolationNotImplemented());
                                }

                              const unsigned int dofs_per_face =
                                cell->get_fe().dofs_per_face;

                              dofs_processed.resize(dofs_per_face);
                              dof_values.resize(dofs_per_face);

                              for (unsigned int dof = 0; dof < dofs_per_face;
                                   ++dof)
                                {
                                  dof_values[dof]     = 0.0;
                                  dofs_processed[dof] = false;
                                }

                              // Compute the projection of the boundary function
                              // on the edge. In 2D this is all that's required.
                              compute_face_projection_curl_conforming_l2(
                                cell,
                                face,
                                fe_face_values,
                                boundary_function,
                                first_vector_component,
                                dof_values,
                                dofs_processed);

                              // store the local->global map:
                              face_dof_indices.resize(dofs_per_face);
                              cell->face(face)->get_dof_indices(
                                face_dof_indices, cell->active_fe_index());

                              // Add the computed constraints to the
                              // AffineConstraints object, assuming the degree
                              // of freedom is not already constrained.
                              for (unsigned int dof = 0; dof < dofs_per_face;
                                   ++dof)
                                {
                                  if (dofs_processed[dof] &&
                                      constraints.can_store_line(
                                        face_dof_indices[dof]) &&
                                      !(constraints.is_constrained(
                                        face_dof_indices[dof])))
                                    {
                                      constraints.add_line(
                                        face_dof_indices[dof]);
                                      if (std::abs(dof_values[dof]) > 1e-13)
                                        {
                                          constraints.set_inhomogeneity(
                                            face_dof_indices[dof],
                                            dof_values[dof]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
              break;
            }

          case 3:
            {
              hp::QCollection<dim> edge_quadrature_collection;

              // Create equivalent of FEEdgeValues:
              for (unsigned int i = 0; i < fe_collection.size(); ++i)
                {
                  const QGauss<dim - 2> reference_edge_quadrature(
                    2 * fe_collection[i].degree + 1);
                  for (const unsigned int face :
                       GeometryInfo<dim>::face_indices())
                    {
                      for (unsigned int line = 0;
                           line < GeometryInfo<dim>::lines_per_face;
                           ++line)
                        {
                          edge_quadrature_collection.push_back(
                            QProjector<dim>::project_to_face(
                              QProjector<dim - 1>::project_to_face(
                                reference_edge_quadrature, line),
                              face));
                        }
                    }
                }

              hp::FEValues<dim> fe_edge_values(mapping_collection,
                                               fe_collection,
                                               edge_quadrature_collection,
                                               update_jacobians |
                                                 update_JxW_values |
                                                 update_quadrature_points |
                                                 update_values);

              for (; cell != dof_handler.end(); ++cell)
                {
                  if (cell->at_boundary() && cell->is_locally_owned())
                    {
                      for (const unsigned int face :
                           GeometryInfo<dim>::face_indices())
                        {
                          if (cell->face(face)->boundary_id() ==
                              boundary_component)
                            {
                              // If the FE is an FE_Nothing object there is no
                              // work to do
                              if (dynamic_cast<const FE_Nothing<dim> *>(
                                    &cell->get_fe()) != nullptr)
                                {
                                  return;
                                }

                              // This is only implemented for FE_Nedelec
                              // elements. If the FE is a FESystem we cannot
                              // check this.
                              if (dynamic_cast<const FESystem<dim> *>(
                                    &cell->get_fe()) == nullptr)
                                {
                                  AssertThrow(
                                    (dynamic_cast<const FE_Nedelec<dim> *>(
                                       &cell->get_fe()) != nullptr) ||
                                      (dynamic_cast<const FE_NedelecSZ<dim> *>(
                                         &cell->get_fe()) != nullptr),
                                    typename FiniteElement<
                                      dim>::ExcInterpolationNotImplemented());
                                }

                              const unsigned int superdegree =
                                cell->get_fe().degree;
                              const unsigned int degree = superdegree - 1;
                              const unsigned int dofs_per_face =
                                cell->get_fe().dofs_per_face;

                              dofs_processed.resize(dofs_per_face);
                              dof_values.resize(dofs_per_face);
                              for (unsigned int dof = 0; dof < dofs_per_face;
                                   ++dof)
                                {
                                  dof_values[dof]     = 0.0;
                                  dofs_processed[dof] = false;
                                }

                              // First compute the projection on the edges.
                              for (unsigned int line = 0;
                                   line < GeometryInfo<3>::lines_per_face;
                                   ++line)
                                {
                                  compute_edge_projection_l2(
                                    cell,
                                    face,
                                    line,
                                    fe_edge_values,
                                    boundary_function,
                                    first_vector_component,
                                    dof_values,
                                    dofs_processed);
                                }

                              // If there are higher order shape functions, then
                              // we still need to compute the face projection
                              if (degree > 0)
                                {
                                  compute_face_projection_curl_conforming_l2(
                                    cell,
                                    face,
                                    fe_face_values,
                                    boundary_function,
                                    first_vector_component,
                                    dof_values,
                                    dofs_processed);
                                }

                              // Store the computed values in the global vector.
                              face_dof_indices.resize(dofs_per_face);
                              cell->face(face)->get_dof_indices(
                                face_dof_indices, cell->active_fe_index());

                              for (unsigned int dof = 0; dof < dofs_per_face;
                                   ++dof)
                                {
                                  if (dofs_processed[dof] &&
                                      constraints.can_store_line(
                                        face_dof_indices[dof]) &&
                                      !(constraints.is_constrained(
                                        face_dof_indices[dof])))
                                    {
                                      constraints.add_line(
                                        face_dof_indices[dof]);

                                      if (std::abs(dof_values[dof]) > 1e-13)
                                        {
                                          constraints.set_inhomogeneity(
                                            face_dof_indices[dof],
                                            dof_values[dof]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
              break;
            }
          default:
            Assert(false, ExcNotImplemented());
        }
    }

  } // namespace internals


  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const DoFHandler<dim> &      dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, number> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<number> &  constraints,
    const Mapping<dim> &         mapping)
  {
    // non-hp version - calls the internal
    // compute_project_boundary_values_curl_conforming_l2() function
    // above after recasting the mapping.

    const hp::MappingCollection<dim> mapping_collection(mapping);
    internals::compute_project_boundary_values_curl_conforming_l2(
      dof_handler,
      first_vector_component,
      boundary_function,
      boundary_component,
      constraints,
      mapping_collection);
  }

  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const hp::DoFHandler<dim> &            dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, number> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<number> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection)
  {
    // hp version - calls the internal
    // compute_project_boundary_values_curl_conforming_l2() function above.
    internals::compute_project_boundary_values_curl_conforming_l2(
      dof_handler,
      first_vector_component,
      boundary_function,
      boundary_component,
      constraints,
      mapping_collection);
  }



  namespace internals
  {
    // This function computes the projection of the boundary function on the
    // boundary in 2d.
    template <typename cell_iterator>
    void
    compute_face_projection_div_conforming(
      const cell_iterator &                       cell,
      const unsigned int                          face,
      const FEFaceValues<2> &                     fe_values,
      const unsigned int                          first_vector_component,
      const Function<2> &                         boundary_function,
      const std::vector<DerivativeForm<1, 2, 2>> &jacobians,
      AffineConstraints<double> &                 constraints)
    {
      // Compute the integral over the product of the normal components of
      // the boundary function times the normal components of the shape
      // functions supported on the boundary.
      const FEValuesExtractors::Vector vec(first_vector_component);
      const FiniteElement<2> &         fe      = cell->get_fe();
      const std::vector<Tensor<1, 2>> &normals = fe_values.get_normal_vectors();
      const unsigned int
                                  face_coordinate_direction[GeometryInfo<2>::faces_per_cell] = {1,
                                                                      1,
                                                                      0,
                                                                      0};
      std::vector<Vector<double>> values(fe_values.n_quadrature_points,
                                         Vector<double>(2));
      Vector<double>              dof_values(fe.dofs_per_face);

      // Get the values of the boundary function at the quadrature points.
      {
        const std::vector<Point<2>> &quadrature_points =
          fe_values.get_quadrature_points();

        boundary_function.vector_value_list(quadrature_points, values);
      }

      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
           ++q_point)
        {
          double tmp = 0.0;

          for (unsigned int d = 0; d < 2; ++d)
            tmp += normals[q_point][d] * values[q_point](d);

          tmp *=
            fe_values.JxW(q_point) *
            std::sqrt(jacobians[q_point][0][face_coordinate_direction[face]] *
                        jacobians[q_point][0][face_coordinate_direction[face]] +
                      jacobians[q_point][1][face_coordinate_direction[face]] *
                        jacobians[q_point][1][face_coordinate_direction[face]]);

          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            dof_values(i) +=
              tmp * (normals[q_point] *
                     fe_values[vec].value(
                       fe.face_to_cell_index(i,
                                             face,
                                             cell->face_orientation(face),
                                             cell->face_flip(face),
                                             cell->face_rotation(face)),
                       q_point));
        }

      std::vector<types::global_dof_index> face_dof_indices(fe.dofs_per_face);

      cell->face(face)->get_dof_indices(face_dof_indices,
                                        cell->active_fe_index());

      // Copy the computed values in the AffineConstraints only, if the degree
      // of freedom is not already constrained.
      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        if (!(constraints.is_constrained(face_dof_indices[i])) &&
            fe.get_nonzero_components(fe.face_to_cell_index(
              i,
              face,
              cell->face_orientation(face),
              cell->face_flip(face),
              cell->face_rotation(face)))[first_vector_component])
          {
            constraints.add_line(face_dof_indices[i]);

            if (std::abs(dof_values(i)) > 1e-14)
              constraints.set_inhomogeneity(face_dof_indices[i], dof_values(i));
          }
    }

    // dummy implementation of above function for all other dimensions
    template <int dim, typename cell_iterator>
    void
    compute_face_projection_div_conforming(
      const cell_iterator &,
      const unsigned int,
      const FEFaceValues<dim> &,
      const unsigned int,
      const Function<dim> &,
      const std::vector<DerivativeForm<1, dim, dim>> &,
      AffineConstraints<double> &)
    {
      Assert(false, ExcNotImplemented());
    }

    // This function computes the projection of the boundary function on the
    // boundary in 3d.
    template <typename cell_iterator>
    void
    compute_face_projection_div_conforming(
      const cell_iterator &                       cell,
      const unsigned int                          face,
      const FEFaceValues<3> &                     fe_values,
      const unsigned int                          first_vector_component,
      const Function<3> &                         boundary_function,
      const std::vector<DerivativeForm<1, 3, 3>> &jacobians,
      std::vector<double> &                       dof_values,
      std::vector<types::global_dof_index> &      projected_dofs)
    {
      // Compute the intergral over the product of the normal components of
      // the boundary function times the normal components of the shape
      // functions supported on the boundary.
      const FEValuesExtractors::Vector vec(first_vector_component);
      const FiniteElement<3> &         fe      = cell->get_fe();
      const std::vector<Tensor<1, 3>> &normals = fe_values.get_normal_vectors();
      const unsigned int
        face_coordinate_directions[GeometryInfo<3>::faces_per_cell][2] = {
          {1, 2}, {1, 2}, {2, 0}, {2, 0}, {0, 1}, {0, 1}};
      std::vector<Vector<double>> values(fe_values.n_quadrature_points,
                                         Vector<double>(3));
      Vector<double>              dof_values_local(fe.dofs_per_face);

      {
        const std::vector<Point<3>> &quadrature_points =
          fe_values.get_quadrature_points();

        boundary_function.vector_value_list(quadrature_points, values);
      }

      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points;
           ++q_point)
        {
          double tmp = 0.0;

          for (unsigned int d = 0; d < 3; ++d)
            tmp += normals[q_point][d] * values[q_point](d);

          tmp *=
            fe_values.JxW(q_point) *
            std::sqrt(
              (jacobians[q_point][0][face_coordinate_directions[face][0]] *
                 jacobians[q_point][0][face_coordinate_directions[face][0]] +
               jacobians[q_point][1][face_coordinate_directions[face][0]] *
                 jacobians[q_point][1][face_coordinate_directions[face][0]] +
               jacobians[q_point][2][face_coordinate_directions[face][0]] *
                 jacobians[q_point][2][face_coordinate_directions[face][0]]) *
              (jacobians[q_point][0][face_coordinate_directions[face][1]] *
                 jacobians[q_point][0][face_coordinate_directions[face][1]] +
               jacobians[q_point][1][face_coordinate_directions[face][1]] *
                 jacobians[q_point][1][face_coordinate_directions[face][1]] +
               jacobians[q_point][2][face_coordinate_directions[face][1]] *
                 jacobians[q_point][2][face_coordinate_directions[face][1]]));

          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            dof_values_local(i) +=
              tmp * (normals[q_point] *
                     fe_values[vec].value(
                       fe.face_to_cell_index(i,
                                             face,
                                             cell->face_orientation(face),
                                             cell->face_flip(face),
                                             cell->face_rotation(face)),
                       q_point));
        }

      std::vector<types::global_dof_index> face_dof_indices(fe.dofs_per_face);

      cell->face(face)->get_dof_indices(face_dof_indices,
                                        cell->active_fe_index());

      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        if (projected_dofs[face_dof_indices[i]] < fe.degree &&
            fe.get_nonzero_components(fe.face_to_cell_index(
              i,
              face,
              cell->face_orientation(face),
              cell->face_flip(face),
              cell->face_rotation(face)))[first_vector_component])
          {
            dof_values[face_dof_indices[i]]     = dof_values_local(i);
            projected_dofs[face_dof_indices[i]] = fe.degree;
          }
    }

    // dummy implementation of above
    // function for all other
    // dimensions
    template <int dim, typename cell_iterator>
    void
    compute_face_projection_div_conforming(
      const cell_iterator &,
      const unsigned int,
      const FEFaceValues<dim> &,
      const unsigned int,
      const Function<dim> &,
      const std::vector<DerivativeForm<1, dim, dim>> &,
      std::vector<double> &,
      std::vector<types::global_dof_index> &)
    {
      Assert(false, ExcNotImplemented());
    }
  } // namespace internals


  template <int dim>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim> &    dof_handler,
    const unsigned int         first_vector_component,
    const Function<dim> &      boundary_function,
    const types::boundary_id   boundary_component,
    AffineConstraints<double> &constraints,
    const Mapping<dim> &       mapping)
  {
    const unsigned int spacedim = dim;
    // Interpolate the normal components
    // of the boundary functions. Since
    // the Raviart-Thomas elements are
    // constructed from a Lagrangian
    // basis, it suffices to compute
    // the integral over the product
    // of the normal components of the
    // boundary function times the
    // normal components of the shape
    // functions supported on the
    // boundary.
    const FiniteElement<dim> &       fe = dof_handler.get_fe();
    QGauss<dim - 1>                  face_quadrature(fe.degree + 1);
    FEFaceValues<dim>                fe_face_values(mapping,
                                     fe,
                                     face_quadrature,
                                     update_JxW_values | update_normal_vectors |
                                       update_quadrature_points |
                                       update_values);
    hp::FECollection<dim>            fe_collection(fe);
    const hp::MappingCollection<dim> mapping_collection(mapping);
    hp::QCollection<dim>             quadrature_collection;

    for (unsigned int face : GeometryInfo<dim>::face_indices())
      quadrature_collection.push_back(
        QProjector<dim>::project_to_face(face_quadrature, face));

    hp::FEValues<dim> fe_values(mapping_collection,
                                fe_collection,
                                quadrature_collection,
                                update_jacobians);

    switch (dim)
      {
        case 2:
          {
            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // if the FE is a
                      // FE_Nothing object
                      // there is no work to
                      // do
                      if (dynamic_cast<const FE_Nothing<dim> *>(
                            &cell->get_fe()) != nullptr)
                        return;

                      // This is only
                      // implemented, if the
                      // FE is a Raviart-Thomas
                      // element. If the FE is
                      // a FESystem we cannot
                      // check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(
                            dynamic_cast<const FE_RaviartThomas<dim> *>(
                              &cell->get_fe()) != nullptr,
                            typename FiniteElement<
                              dim>::ExcInterpolationNotImplemented());
                        }

                      fe_values.reinit(cell,
                                       face +
                                         cell->active_fe_index() *
                                           GeometryInfo<dim>::faces_per_cell);

                      const std::vector<DerivativeForm<1, dim, spacedim>>
                        &jacobians =
                          fe_values.get_present_fe_values().get_jacobians();

                      fe_face_values.reinit(cell, face);
                      internals::compute_face_projection_div_conforming(
                        cell,
                        face,
                        fe_face_values,
                        first_vector_component,
                        boundary_function,
                        jacobians,
                        constraints);
                    }

            break;
          }

        case 3:
          {
            // In three dimensions the edges between two faces are treated
            // twice. Therefore we store the computed values in a vector
            // and copy them over in the AffineConstraints after all values
            // have been computed. If we have two values for one edge, we
            // choose the one, which was computed with the higher order
            // element. If both elements are of the same order, we just
            // keep the first value and do not compute a second one.
            const unsigned int                   n_dofs = dof_handler.n_dofs();
            std::vector<double>                  dof_values(n_dofs);
            std::vector<types::global_dof_index> projected_dofs(n_dofs);

            for (unsigned int dof = 0; dof < n_dofs; ++dof)
              projected_dofs[dof] = 0;

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // This is only implemented, if the FE is a
                      // Raviart-Thomas element. If the FE is a FESystem we
                      // cannot check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(
                            dynamic_cast<const FE_RaviartThomas<dim> *>(
                              &cell->get_fe()) != nullptr,
                            typename FiniteElement<
                              dim>::ExcInterpolationNotImplemented());
                        }

                      fe_values.reinit(cell,
                                       face +
                                         cell->active_fe_index() *
                                           GeometryInfo<dim>::faces_per_cell);

                      const std::vector<DerivativeForm<1, dim, spacedim>>
                        &jacobians =
                          fe_values.get_present_fe_values().get_jacobians();

                      fe_face_values.reinit(cell, face);
                      internals::compute_face_projection_div_conforming(
                        cell,
                        face,
                        fe_face_values,
                        first_vector_component,
                        boundary_function,
                        jacobians,
                        dof_values,
                        projected_dofs);
                    }

            for (unsigned int dof = 0; dof < n_dofs; ++dof)
              if ((projected_dofs[dof] != 0) &&
                  !(constraints.is_constrained(dof)))
                {
                  constraints.add_line(dof);

                  if (std::abs(dof_values[dof]) > 1e-14)
                    constraints.set_inhomogeneity(dof, dof_values[dof]);
                }

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  }


  template <int dim>
  void
  project_boundary_values_div_conforming(
    const hp::DoFHandler<dim> &            dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim> &                  boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<double> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection)
  {
    const unsigned int           spacedim = dim;
    const hp::FECollection<dim> &fe_collection =
      dof_handler.get_fe_collection();
    hp::QCollection<dim - 1> face_quadrature_collection;
    hp::QCollection<dim>     quadrature_collection;

    for (unsigned int i = 0; i < fe_collection.size(); ++i)
      {
        const QGauss<dim - 1> quadrature(fe_collection[i].degree + 1);

        face_quadrature_collection.push_back(quadrature);

        for (unsigned int face : GeometryInfo<dim>::face_indices())
          quadrature_collection.push_back(
            QProjector<dim>::project_to_face(quadrature, face));
      }

    hp::FEFaceValues<dim> fe_face_values(mapping_collection,
                                         fe_collection,
                                         face_quadrature_collection,
                                         update_JxW_values |
                                           update_normal_vectors |
                                           update_quadrature_points |
                                           update_values);
    hp::FEValues<dim>     fe_values(mapping_collection,
                                fe_collection,
                                quadrature_collection,
                                update_jacobians);

    switch (dim)
      {
        case 2:
          {
            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // This is only
                      // implemented, if the
                      // FE is a Raviart-Thomas
                      // element. If the FE is
                      // a FESystem we cannot
                      // check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(
                            dynamic_cast<const FE_RaviartThomas<dim> *>(
                              &cell->get_fe()) != nullptr,
                            typename FiniteElement<
                              dim>::ExcInterpolationNotImplemented());
                        }

                      fe_values.reinit(cell,
                                       face +
                                         cell->active_fe_index() *
                                           GeometryInfo<dim>::faces_per_cell);

                      const std::vector<DerivativeForm<1, dim, spacedim>>
                        &jacobians =
                          fe_values.get_present_fe_values().get_jacobians();

                      fe_face_values.reinit(cell, face);
                      internals::compute_face_projection_div_conforming(
                        cell,
                        face,
                        fe_face_values.get_present_fe_values(),
                        first_vector_component,
                        boundary_function,
                        jacobians,
                        constraints);
                    }

            break;
          }

        case 3:
          {
            const unsigned int                   n_dofs = dof_handler.n_dofs();
            std::vector<double>                  dof_values(n_dofs);
            std::vector<types::global_dof_index> projected_dofs(n_dofs);

            for (unsigned int dof = 0; dof < n_dofs; ++dof)
              projected_dofs[dof] = 0;

            for (const auto &cell : dof_handler.active_cell_iterators())
              if (cell->at_boundary() && cell->is_locally_owned())
                for (const unsigned int face :
                     GeometryInfo<dim>::face_indices())
                  if (cell->face(face)->boundary_id() == boundary_component)
                    {
                      // This is only
                      // implemented, if the
                      // FE is a Raviart-Thomas
                      // element. If the FE is
                      // a FESystem we cannot
                      // check this.
                      if (dynamic_cast<const FESystem<dim> *>(
                            &cell->get_fe()) == nullptr)
                        {
                          AssertThrow(
                            dynamic_cast<const FE_RaviartThomas<dim> *>(
                              &cell->get_fe()) != nullptr,
                            typename FiniteElement<
                              dim>::ExcInterpolationNotImplemented());
                        }

                      fe_values.reinit(cell,
                                       face +
                                         cell->active_fe_index() *
                                           GeometryInfo<dim>::faces_per_cell);

                      const std::vector<DerivativeForm<1, dim, spacedim>>
                        &jacobians =
                          fe_values.get_present_fe_values().get_jacobians();

                      fe_face_values.reinit(cell, face);
                      internals::compute_face_projection_div_conforming(
                        cell,
                        face,
                        fe_face_values.get_present_fe_values(),
                        first_vector_component,
                        boundary_function,
                        jacobians,
                        dof_values,
                        projected_dofs);
                    }

            for (unsigned int dof = 0; dof < n_dofs; ++dof)
              if ((projected_dofs[dof] != 0) &&
                  !(constraints.is_constrained(dof)))
                {
                  constraints.add_line(dof);

                  if (std::abs(dof_values[dof]) > 1e-14)
                    constraints.set_inhomogeneity(dof, dof_values[dof]);
                }

            break;
          }

        default:
          Assert(false, ExcNotImplemented());
      }
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_templates_h
