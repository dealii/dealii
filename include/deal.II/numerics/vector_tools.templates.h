// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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


#ifndef _deal2__vector_tools_templates_h
#define _deal2__vector_tools_templates_h

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/base/std_cxx11/array.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <set>
#include <list>

DEAL_II_NAMESPACE_OPEN


namespace VectorTools
{

  template <class VECTOR, int dim, int spacedim, template <int,int> class DH>
  void interpolate (const Mapping<dim,spacedim>    &mapping,
                    const DH<dim,spacedim>         &dof,
                    const Function<spacedim>       &function,
                    VECTOR                         &vec)
  {
    Assert (dof.get_fe().n_components() == function.n_components,
            ExcDimensionMismatch(dof.get_fe().n_components(),
                                 function.n_components));

    const hp::FECollection<dim,spacedim> fe (dof.get_fe());
    const unsigned int          n_components = fe.n_components();
    const bool                  fe_is_system = (n_components != 1);

    typename DH<dim,spacedim>::active_cell_iterator cell = dof.begin_active(),
                                                    endc = dof.end();

    // For FESystems many of the
    // unit_support_points will appear
    // multiple times, as a point may be
    // unit_support_point for several of the
    // components of the system.  The following
    // is rather complicated, but at least
    // attempts to avoid evaluating the
    // vectorfunction multiple times at the
    // same point on a cell.
    //
    // note that we have to set up all of the
    // following arrays for each of the
    // elements in the FECollection (which
    // means only once if this is for a regular
    // DoFHandler)
    std::vector<std::vector<Point<dim> > > unit_support_points (fe.size());
    for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
      {
        unit_support_points[fe_index] = fe[fe_index].get_unit_support_points();
        Assert (unit_support_points[fe_index].size() != 0,
                ExcNonInterpolatingFE());
      }


    // Find the support points on a cell that
    // are mentioned multiple times in
    // unit_support_points.  Mark the first
    // representative of each support point
    // mentioned multiple times by appending
    // its dof index to dofs_of_rep_points.
    // Each multiple point gets to know the dof
    // index of its representative point by the
    // dof_to_rep_dof_table.

    // the following vector collects all dofs i,
    // 0<=i<fe.dofs_per_cell, for that
    // unit_support_points[i]
    // is a representative one. i.e.
    // the following vector collects all rep dofs.
    // the position of a rep dof within this vector
    // is called rep index.
    std::vector<std::vector<types::global_dof_index> > dofs_of_rep_points(fe.size());
    // the following table converts a dof i
    // to the rep index.
    std::vector<std::vector<types::global_dof_index> > dof_to_rep_index_table(fe.size());

    std::vector<unsigned int> n_rep_points (fe.size(), 0);

    for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
      {
        for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
          {
            bool representative=true;
            // the following loop is looped
            // the other way round to get
            // the minimal effort of
            // O(fe.dofs_per_cell) for multiple
            // support points that are placed
            // one after the other.
            for (unsigned int j=dofs_of_rep_points[fe_index].size(); j>0; --j)
              if (unit_support_points[fe_index][i]
                  == unit_support_points[fe_index][dofs_of_rep_points[fe_index][j-1]])
                {
                  dof_to_rep_index_table[fe_index].push_back(j-1);
                  representative=false;
                  break;
                }

            if (representative)
              {
                // rep_index=dofs_of_rep_points.size()
                dof_to_rep_index_table[fe_index].push_back(dofs_of_rep_points[fe_index].size());
                // dofs_of_rep_points[rep_index]=i
                dofs_of_rep_points[fe_index].push_back(i);
                ++n_rep_points[fe_index];
              }
          }

        Assert(dofs_of_rep_points[fe_index].size()==n_rep_points[fe_index],
               ExcInternalError());
        Assert(dof_to_rep_index_table[fe_index].size()==fe[fe_index].dofs_per_cell,
               ExcInternalError());
      }

    const unsigned int max_rep_points = *std::max_element (n_rep_points.begin(),
                                                           n_rep_points.end());
    std::vector<types::global_dof_index> dofs_on_cell (fe.max_dofs_per_cell());
    std::vector<Point<spacedim> >  rep_points (max_rep_points);

    // get space for the values of the
    // function at the rep support points.
    //
    // have two versions, one for system fe
    // and one for scalar ones, to take the
    // more efficient one respectively
    std::vector<std::vector<double> >         function_values_scalar(fe.size());
    std::vector<std::vector<Vector<double> > > function_values_system(fe.size());

    // Make a quadrature rule from support points
    // to feed it into FEValues
    hp::QCollection<dim> support_quadrature;
    for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
      support_quadrature.push_back (Quadrature<dim>(unit_support_points[fe_index]));

    // Transformed support points are computed by
    // FEValues
    hp::MappingCollection<dim,spacedim> mapping_collection (mapping);

    hp::FEValues<dim,spacedim> fe_values (mapping_collection,
                                          fe, support_quadrature, update_quadrature_points);

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          const unsigned int fe_index = cell->active_fe_index();

          // for each cell:
          // get location of finite element
          // support_points
          fe_values.reinit(cell);
          const std::vector<Point<spacedim> > &support_points =
            fe_values.get_present_fe_values().get_quadrature_points();

          // pick out the representative
          // support points
          rep_points.resize (dofs_of_rep_points[fe_index].size());
          for (unsigned int j=0; j<dofs_of_rep_points[fe_index].size(); ++j)
            rep_points[j] = support_points[dofs_of_rep_points[fe_index][j]];

          // get indices of the dofs on this cell
          dofs_on_cell.resize (fe[fe_index].dofs_per_cell);
          cell->get_dof_indices (dofs_on_cell);


          if (fe_is_system)
            {
              // get function values at
              // these points. Here: get
              // all components
              function_values_system[fe_index]
              .resize (n_rep_points[fe_index],
                       Vector<double> (fe[fe_index].n_components()));
              function.vector_value_list (rep_points,
                                          function_values_system[fe_index]);
              // distribute the function
              // values to the global
              // vector
              for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
                {
                  const unsigned int component
                    = fe[fe_index].system_to_component_index(i).first;
                  const unsigned int rep_dof=dof_to_rep_index_table[fe_index][i];
                  vec(dofs_on_cell[i])
                    = function_values_system[fe_index][rep_dof](component);
                }
            }
          else
            {
              // get first component only,
              // which is the only component
              // in the function anyway
              function_values_scalar[fe_index].resize (n_rep_points[fe_index]);
              function.value_list (rep_points,
                                   function_values_scalar[fe_index],
                                   0);
              // distribute the function
              // values to the global
              // vector
              for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
                vec(dofs_on_cell[i])
                  = function_values_scalar[fe_index][dof_to_rep_index_table[fe_index][i]];
            }
        }
    vec.compress(VectorOperation::insert);
  }


  template <class VECTOR, class DH>
  void interpolate (const DH              &dof,
                    const Function<DH::space_dimension>   &function,
                    VECTOR                &vec)
  {
    interpolate(StaticMappingQ1<DH::dimension, DH::space_dimension>::mapping,
                dof, function, vec);
  }




  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolate (const DoFHandler<dim,spacedim>           &dof_1,
               const DoFHandler<dim,spacedim>           &dof_2,
               const FullMatrix<double>        &transfer,
               const InVector                  &data_1,
               OutVector                       &data_2)
  {
    Vector<double> cell_data_1(dof_1.get_fe().dofs_per_cell);
    Vector<double> cell_data_2(dof_2.get_fe().dofs_per_cell);

    std::vector<short unsigned int> touch_count (dof_2.n_dofs(), 0); //TODO: check on datatype... kinda strange (UK)
    std::vector<types::global_dof_index>       local_dof_indices (dof_2.get_fe().dofs_per_cell);

    typename DoFHandler<dim,spacedim>::active_cell_iterator h = dof_1.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator l = dof_2.begin_active();
    const typename DoFHandler<dim,spacedim>::cell_iterator endh = dof_1.end();

    for (; h != endh; ++h, ++l)
      {
        h->get_dof_values(data_1, cell_data_1);
        transfer.vmult(cell_data_2, cell_data_1);

        l->get_dof_indices (local_dof_indices);

        // distribute cell vector
        for (unsigned int j=0; j<dof_2.get_fe().dofs_per_cell; ++j)
          {
            data_2(local_dof_indices[j]) += cell_data_2(j);

            // count, how often we have
            // added to this dof
            Assert (touch_count[local_dof_indices[j]] < numbers::internal_face_boundary_id,
                    ExcInternalError());
            ++touch_count[local_dof_indices[j]];
          }
      }

    // compute the mean value of the
    // sum which we have placed in each
    // entry of the output vector
    for (unsigned int i=0; i<dof_2.n_dofs(); ++i)
      {
        Assert (touch_count[i] != 0,
                ExcInternalError());

        data_2(i) /= touch_count[i];
      }
  }


  template<typename VECTOR, typename DH>
  void
  interpolate_based_on_material_id(const Mapping<DH::dimension, DH::space_dimension>                          &mapping,
                                   const DH                                                                   &dof,
                                   const std::map< types::material_id, const Function<DH::space_dimension>* > &function_map,
                                   VECTOR                                                                     &dst,
                                   const ComponentMask                                                        &component_mask)
  {
    const unsigned int dim = DH::dimension;

    Assert( component_mask.represents_n_components(dof.get_fe().n_components()),
            ExcMessage("The number of components in the mask has to be either "
                       "zero or equal to the number of components in the finite "
                       "element.") );

    if ( function_map.size() == 0 )
      return;

    Assert( function_map.find(numbers::invalid_material_id) == function_map.end(),
            ExcInvalidMaterialIndicator() );

    for ( typename std::map< types::material_id, const Function<DH::space_dimension>* >::const_iterator
          iter  = function_map.begin();
          iter != function_map.end();
          ++iter )
      {
        Assert( dof.get_fe().n_components() == iter->second->n_components,
                ExcDimensionMismatch(dof.get_fe().n_components(), iter->second->n_components) );
      }

    const hp::FECollection<DH::dimension, DH::space_dimension> fe(dof.get_fe());
    const unsigned int n_components =  fe.n_components();
    const bool         fe_is_system = (n_components != 1);

    typename DH::active_cell_iterator cell = dof.begin_active(),
                                      endc = dof.end();

    std::vector< std::vector< Point<dim> > > unit_support_points(fe.size());
    for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
      {
        unit_support_points[fe_index] = fe[fe_index].get_unit_support_points();
        Assert( unit_support_points[fe_index].size() != 0,
                ExcNonInterpolatingFE() );
      }

    std::vector< std::vector<unsigned int> > dofs_of_rep_points(fe.size());
    std::vector< std::vector<unsigned int> > dof_to_rep_index_table(fe.size());
    std::vector<unsigned int>                n_rep_points(fe.size(), 0);

    for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
      {
        for (unsigned int i = 0; i < fe[fe_index].dofs_per_cell; ++i)
          {
            bool representative = true;

            for (unsigned int j = dofs_of_rep_points[fe_index].size(); j > 0; --j)
              if ( unit_support_points[fe_index][i] == unit_support_points[fe_index][dofs_of_rep_points[fe_index][j-1]] )
                {
                  dof_to_rep_index_table[fe_index].push_back(j-1);
                  representative = false;
                  break;
                }

            if (representative)
              {
                dof_to_rep_index_table[fe_index].push_back(dofs_of_rep_points[fe_index].size());
                dofs_of_rep_points[fe_index].push_back(i);
                ++n_rep_points[fe_index];
              }
          }

        Assert( dofs_of_rep_points[fe_index].size() == n_rep_points[fe_index],
                ExcInternalError() );
        Assert( dof_to_rep_index_table[fe_index].size() == fe[fe_index].dofs_per_cell,
                ExcInternalError() );
      }

    const unsigned int max_rep_points = *std::max_element(n_rep_points.begin(),
                                                          n_rep_points.end());
    std::vector< types::global_dof_index>     dofs_on_cell(fe.max_dofs_per_cell());
    std::vector< Point<DH::space_dimension> > rep_points(max_rep_points);

    std::vector< std::vector<double> >           function_values_scalar(fe.size());
    std::vector< std::vector< Vector<double> > > function_values_system(fe.size());

    hp::QCollection<dim> support_quadrature;
    for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
      support_quadrature.push_back( Quadrature<dim>(unit_support_points[fe_index]) );

    hp::MappingCollection<dim, DH::space_dimension> mapping_collection(mapping);
    hp::FEValues<dim, DH::space_dimension> fe_values(mapping_collection,
                                                     fe,
                                                     support_quadrature,
                                                     update_quadrature_points);

    for ( ; cell != endc; ++cell)
      if ( cell->is_locally_owned() )
        if ( function_map.find(cell->material_id()) != function_map.end() )
          {
            const unsigned int fe_index = cell->active_fe_index();

            fe_values.reinit(cell);

            const std::vector< Point<DH::space_dimension> > &support_points = fe_values.get_present_fe_values().get_quadrature_points();

            rep_points.resize( dofs_of_rep_points[fe_index].size() );
            for (unsigned int i = 0; i < dofs_of_rep_points[fe_index].size(); ++i)
              rep_points[i] = support_points[dofs_of_rep_points[fe_index][i]];

            dofs_on_cell.resize( fe[fe_index].dofs_per_cell );
            cell->get_dof_indices(dofs_on_cell);

            if (fe_is_system)
              {
                function_values_system[fe_index].resize( n_rep_points[fe_index],
                                                         Vector<double>(fe[fe_index].n_components()) );

                function_map.find(cell->material_id())->second->vector_value_list(rep_points,
                    function_values_system[fe_index]);

                for (unsigned int i = 0; i < fe[fe_index].dofs_per_cell; ++i)
                  {
                    const unsigned int component = fe[fe_index].system_to_component_index(i).first;

                    if ( component_mask[component] )
                      {
                        const unsigned int rep_dof = dof_to_rep_index_table[fe_index][i];
                        dst(dofs_on_cell[i])       = function_values_system[fe_index][rep_dof](component);
                      }
                  }
              }
            else
              {
                function_values_scalar[fe_index].resize(n_rep_points[fe_index]);

                function_map.find(cell->material_id())->second->value_list(rep_points,
                                                                           function_values_scalar[fe_index],
                                                                           0);

                for (unsigned int i = 0; i < fe[fe_index].dofs_per_cell; ++i)
                  dst(dofs_on_cell[i]) = function_values_scalar[fe_index][dof_to_rep_index_table[fe_index][i]];
              }
          }

    dst.compress (VectorOperation::insert);
  }


  namespace internal
  {
    /**
     * Interpolate zero boundary values. We don't need to worry about a mapping
     * here because the function we evaluate for the DoFs is zero in the mapped
     * locations as well as in the original, unmapped locations
     */
    template <class DH>
    void
    interpolate_zero_boundary_values (const DH                                 &dof_handler,
                                      std::map<types::global_dof_index,double> &boundary_values)
    {
      const unsigned int dim      = DH::dimension;

      // loop over all boundary faces
      // to get all dof indices of
      // dofs on the boundary. note
      // that in 3d there are cases
      // where a face is not at the
      // boundary, yet one of its
      // lines is, and we should
      // consider the degrees of
      // freedom on it as boundary
      // nodes. likewise, in 2d and
      // 3d there are cases where a
      // cell is only at the boundary
      // by one vertex. nevertheless,
      // since we do not support
      // boundaries with dimension
      // less or equal to dim-2, each
      // such boundary dof is also
      // found from some other face
      // that is actually wholly on
      // the boundary, not only by
      // one line or one vertex
      typename DH::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      std::vector<types::global_dof_index> face_dof_indices;
      for (; cell!=endc; ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f))
            {
              face_dof_indices.resize (cell->get_fe().dofs_per_face);
              cell->face(f)->get_dof_indices (face_dof_indices,
                                              cell->active_fe_index());
              for (unsigned int i=0; i<cell->get_fe().dofs_per_face; ++i)
                // enter zero boundary values
                // for all boundary nodes
                //
                // we need not care about
                // vector valued elements here,
                // since we set all components
                boundary_values[face_dof_indices[i]] = 0.;
            }
    }
  }



  template <int dim, int spacedim,
            template <int,int> class DH,
            class Vector>
  void
  interpolate_to_different_mesh (const DH<dim, spacedim> &dof1,
                                 const Vector            &u1,
                                 const DH<dim, spacedim> &dof2,
                                 Vector                  &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage ("The two containers must represent triangulations that "
                       "have the same coarse meshes"));

    InterGridMap<DH<dim, spacedim> > intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    ConstraintMatrix dummy;
    dummy.close();

    interpolate_to_different_mesh(intergridmap, u1, dummy, u2);
  }



  template <int dim, int spacedim,
            template <int,int> class DH,
            class Vector>
  void
  interpolate_to_different_mesh (const DH<dim, spacedim> &dof1,
                                 const Vector            &u1,
                                 const DH<dim, spacedim> &dof2,
                                 const ConstraintMatrix  &constraints,
                                 Vector                  &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage ("The two containers must represent triangulations that "
                       "have the same coarse meshes"));

    InterGridMap<DH<dim, spacedim> > intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    interpolate_to_different_mesh(intergridmap, u1, constraints, u2);
  }



  template <int dim, int spacedim,
            template <int,int> class DH,
            class Vector>
  void
  interpolate_to_different_mesh (const InterGridMap<DH<dim, spacedim> > &intergridmap,
                                 const Vector                           &u1,
                                 const ConstraintMatrix                 &constraints,
                                 Vector                                 &u2)
  {
    const DH<dim, spacedim> &dof1 = intergridmap.get_source_grid();
    const DH<dim, spacedim> &dof2 = intergridmap.get_destination_grid();

    Assert(u1.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(),dof2.n_dofs()));

    Vector cache;

    // Looping over the finest common
    // mesh, this means that source and
    // destination cells have to be on the
    // same level and at least one has to
    // be active.
    //
    // Therefor, loop over all cells
    // (active and inactive) of the source
    // grid ..
    typename DH<dim,spacedim>::cell_iterator       cell1 = dof1.begin();
    const typename DH<dim,spacedim>::cell_iterator endc1 = dof1.end();

    for (; cell1 != endc1; ++cell1)
      {
        const typename DH<dim,spacedim>::cell_iterator cell2 = intergridmap[cell1];

        // .. and skip if source and destination
        // cells are not on the same level ..
        if (cell1->level() != cell2->level())
          continue;
        // .. or none of them is active.
        if (!cell1->active() && !cell2->active())
          continue;

        Assert(cell1->get_fe().get_name() ==
               cell2->get_fe().get_name(),
               ExcMessage ("Source and destination cells need to use the same finite element"));

        cache.reinit(cell1->get_fe().dofs_per_cell);

        // Get and set the corresponding
        // dof_values by interpolation.
        cell1->get_interpolated_dof_values(u1, cache);
        cell2->set_dof_values_by_interpolation(cache, u2);
      }

    // Apply hanging node constraints.
    constraints.distribute(u2);
  }


  namespace
  {
    /**
     * Compute the boundary values to be used in the project() functions.
     */
    template <int dim, int spacedim,
              template <int,int> class DH,
              template <int,int> class M_or_MC,
              template <int> class Q_or_QC>
    void project_compute_b_v (const M_or_MC<dim, spacedim>   &mapping,
                              const DH<dim,spacedim> &dof,
                              const Function<spacedim> &function,
                              const bool                enforce_zero_boundary,
                              const Q_or_QC<dim-1>  &q_boundary,
                              const bool                project_to_boundary_first,
                              std::map<types::global_dof_index,double> &boundary_values)
    {
      if (enforce_zero_boundary == true)
        // no need to project boundary
        // values, but enforce
        // homogeneous boundary values
        // anyway
        internal::
        interpolate_zero_boundary_values (dof, boundary_values);

      else
        // no homogeneous boundary values
        if (project_to_boundary_first == true)
          // boundary projection required
          {
            // set up a list of boundary
            // functions for the
            // different boundary
            // parts. We want the
            // function to hold on
            // all parts of the boundary
            const std::vector<types::boundary_id>
            used_boundary_indicators = dof.get_tria().get_boundary_indicators();

            typename FunctionMap<spacedim>::type boundary_functions;
            for (unsigned int i=0; i<used_boundary_indicators.size(); ++i)
              boundary_functions[used_boundary_indicators[i]] = &function;
            project_boundary_values (mapping, dof, boundary_functions, q_boundary,
                                     boundary_values);
          }
    }


    /**
     * Return whether the boundary values try to constrain a degree of
     * freedom that is already constrained to something else
     */
    bool constraints_and_b_v_are_compatible (const ConstraintMatrix   &constraints,
                                             std::map<types::global_dof_index,double> &boundary_values)
    {
      for (std::map<types::global_dof_index,double>::iterator it=boundary_values.begin();
           it != boundary_values.end(); ++it)
        if (constraints.is_constrained(it->first))
//TODO: This looks wrong -- shouldn't it be ==0 in the first condition and && ?
          if (!(constraints.get_constraint_entries(it->first)->size() > 0
                ||
                (constraints.get_inhomogeneity(it->first) == it->second)))
            return false;

      return true;
    }


    /**
     * Generic implementation of the project() function
     */
    template <int dim, int spacedim,
              class Vector,
              template <int,int> class DH,
              template <int,int> class M_or_MC,
              template <int> class Q_or_QC>
    void do_project (const M_or_MC<dim, spacedim>   &mapping,
                     const DH<dim,spacedim> &dof,
                     const ConstraintMatrix   &constraints,
                     const Q_or_QC<dim>       &quadrature,
                     const Function<spacedim> &function,
                     Vector                   &vec_result,
                     const bool                enforce_zero_boundary,
                     const Q_or_QC<dim-1>  &q_boundary,
                     const bool                project_to_boundary_first)
    {
      Assert (dof.get_fe().n_components() == function.n_components,
              ExcDimensionMismatch(dof.get_fe().n_components(),
                                   function.n_components));
      Assert (vec_result.size() == dof.n_dofs(),
              ExcDimensionMismatch (vec_result.size(), dof.n_dofs()));

      // make up boundary values
      std::map<types::global_dof_index,double> boundary_values;
      project_compute_b_v(mapping, dof, function, enforce_zero_boundary,
                          q_boundary, project_to_boundary_first, boundary_values);

      // check if constraints are compatible (see below)
      const bool constraints_are_compatible =
        constraints_and_b_v_are_compatible(constraints, boundary_values);

      // set up mass matrix and right hand side
      dealii::Vector<double> vec (dof.n_dofs());
      SparsityPattern sparsity;
      {
        CompressedSimpleSparsityPattern csp (dof.n_dofs(), dof.n_dofs());
        DoFTools::make_sparsity_pattern (dof, csp, constraints,
                                         !constraints_are_compatible);

        sparsity.copy_from (csp);
      }
      SparseMatrix<double> mass_matrix (sparsity);
      dealii::Vector<double> tmp (mass_matrix.n());

      // If the constraint matrix does not conflict with the given boundary
      // values (i.e., it either does not contain boundary values or it contains
      // the same as boundary_values), we can let it call
      // distribute_local_to_global straight away, otherwise we need to first
      // interpolate the boundary values and then condense the matrix and vector
      if (constraints_are_compatible)
        {
          const Function<spacedim> *dummy = 0;
          MatrixCreator::create_mass_matrix (mapping, dof, quadrature,
                                             mass_matrix, function, tmp,
                                             dummy, constraints);
          if (boundary_values.size() > 0)
            MatrixTools::apply_boundary_values (boundary_values,
                                                mass_matrix, vec, tmp,
                                                true);
        }
      else
        {
          // create mass matrix and rhs at once, which is faster.
          MatrixCreator::create_mass_matrix (mapping, dof, quadrature,
                                             mass_matrix, function, tmp);
          MatrixTools::apply_boundary_values (boundary_values,
                                              mass_matrix, vec, tmp,
                                              true);
          constraints.condense(mass_matrix, tmp);
        }

      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices
      ReductionControl      control(5*tmp.size(), 0., 1e-12, false, false);
      GrowingVectorMemory<> memory;
      SolverCG<>            cg(control,memory);

      PreconditionSSOR<> prec;
      prec.initialize(mass_matrix, 1.2);

      cg.solve (mass_matrix, vec, tmp, prec);

      constraints.distribute (vec);

      // copy vec into vec_result. we can't use vec_result itself above, since
      // it may be of another type than Vector<double> and that wouldn't
      // necessarily go together with the matrix and other functions
      for (unsigned int i=0; i<vec.size(); ++i)
        vec_result(i) = vec(i);
    }
  }


  template <int dim, class VECTOR, int spacedim>
  void project (const Mapping<dim, spacedim>   &mapping,
                const DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix   &constraints,
                const Quadrature<dim>    &quadrature,
                const Function<spacedim> &function,
                VECTOR                   &vec_result,
                const bool                enforce_zero_boundary,
                const Quadrature<dim-1>  &q_boundary,
                const bool                project_to_boundary_first)
  {
    do_project (mapping, dof, constraints, quadrature,
                function, vec_result,
                enforce_zero_boundary, q_boundary,
                project_to_boundary_first);
  }


  template <int dim, class VECTOR, int spacedim>
  void project (const DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix   &constraints,
                const Quadrature<dim>    &quadrature,
                const Function<spacedim> &function,
                VECTOR                   &vec,
                const bool                enforce_zero_boundary,
                const Quadrature<dim-1>  &q_boundary,
                const bool                project_to_boundary_first)
  {
    project(StaticMappingQ1<dim,spacedim>::mapping, dof, constraints, quadrature, function, vec,
            enforce_zero_boundary, q_boundary, project_to_boundary_first);
  }



  template <int dim, class VECTOR, int spacedim>
  void project (const hp::MappingCollection<dim, spacedim>   &mapping,
                const hp::DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix   &constraints,
                const hp::QCollection<dim>    &quadrature,
                const Function<spacedim> &function,
                VECTOR                   &vec_result,
                const bool                enforce_zero_boundary,
                const hp::QCollection<dim-1>  &q_boundary,
                const bool                project_to_boundary_first)
  {
    do_project (mapping, dof, constraints, quadrature,
                function, vec_result,
                enforce_zero_boundary, q_boundary,
                project_to_boundary_first);
  }


  template <int dim, class VECTOR, int spacedim>
  void project (const hp::DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix   &constraints,
                const hp::QCollection<dim>    &quadrature,
                const Function<spacedim> &function,
                VECTOR                   &vec,
                const bool                enforce_zero_boundary,
                const hp::QCollection<dim-1>  &q_boundary,
                const bool                project_to_boundary_first)
  {
    project(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
            dof, constraints, quadrature, function, vec,
            enforce_zero_boundary, q_boundary, project_to_boundary_first);
  }


  template <int dim, int spacedim>
  void create_right_hand_side (const Mapping<dim, spacedim>    &mapping,
                               const DoFHandler<dim,spacedim> &dof_handler,
                               const Quadrature<dim> &quadrature,
                               const Function<spacedim>   &rhs_function,
                               Vector<double>        &rhs_vector)
  {
    const FiniteElement<dim,spacedim> &fe  = dof_handler.get_fe();
    Assert (fe.n_components() == rhs_function.n_components,
            ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    rhs_vector = 0;

    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    FEValues<dim,spacedim> fe_values (mapping, fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points,
                       n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs (dofs_per_cell);
    Vector<double> cell_vector (dofs_per_cell);

    typename DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    if (n_components==1)
      {
        std::vector<double> rhs_values(n_q_points);

        for (; cell!=endc; ++cell)
          {
            fe_values.reinit(cell);

            const std::vector<double> &weights   = fe_values.get_JxW_values ();
            rhs_function.value_list (fe_values.get_quadrature_points(),
                                     rhs_values);

            cell_vector = 0;
            for (unsigned int point=0; point<n_q_points; ++point)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                cell_vector(i) += rhs_values[point] *
                                  fe_values.shape_value(i,point) *
                                  weights[point];

            cell->get_dof_indices (dofs);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              rhs_vector(dofs[i]) += cell_vector(i);
          }

      }
    else
      {
        std::vector<Vector<double> > rhs_values(n_q_points,
                                                Vector<double>(n_components));

        for (; cell!=endc; ++cell)
          {
            fe_values.reinit(cell);

            const std::vector<double> &weights   = fe_values.get_JxW_values ();
            rhs_function.vector_value_list (fe_values.get_quadrature_points(),
                                            rhs_values);

            cell_vector = 0;
            // Use the faster code if the
            // FiniteElement is primitive
            if (fe.is_primitive ())
              {
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                      const unsigned int component
                        = fe.system_to_component_index(i).first;

                      cell_vector(i) += rhs_values[point](component) *
                                        fe_values.shape_value(i,point) *
                                        weights[point];
                    }
              }
            else
              {
                // Otherwise do it the way
                // proposed for vector valued
                // elements
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                      if (fe.get_nonzero_components(i)[comp_i])
                        {
                          cell_vector(i) += rhs_values[point](comp_i) *
                                            fe_values.shape_value_component(i,point,comp_i) *
                                            weights[point];
                        }
              }

            cell->get_dof_indices (dofs);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              rhs_vector(dofs[i]) += cell_vector(i);
          }
      }
  }



  template <int dim, int spacedim>
  void create_right_hand_side (const DoFHandler<dim,spacedim>    &dof_handler,
                               const Quadrature<dim>    &quadrature,
                               const Function<spacedim>      &rhs_function,
                               Vector<double>           &rhs_vector)
  {
    create_right_hand_side(StaticMappingQ1<dim,spacedim>::mapping, dof_handler, quadrature,
                           rhs_function, rhs_vector);
  }




  template <int dim, int spacedim>
  void create_right_hand_side (const hp::MappingCollection<dim,spacedim>    &mapping,
                               const hp::DoFHandler<dim,spacedim> &dof_handler,
                               const hp::QCollection<dim> &quadrature,
                               const Function<spacedim>   &rhs_function,
                               Vector<double>        &rhs_vector)
  {
    const hp::FECollection<dim,spacedim> &fe  = dof_handler.get_fe();
    Assert (fe.n_components() == rhs_function.n_components,
            ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    rhs_vector = 0;

    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    hp::FEValues<dim,spacedim> x_fe_values (mapping, fe, quadrature, update_flags);

    const unsigned int n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs (fe.max_dofs_per_cell());
    Vector<double> cell_vector (fe.max_dofs_per_cell());

    typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    if (n_components==1)
      {
        std::vector<double> rhs_values;

        for (; cell!=endc; ++cell)
          {
            x_fe_values.reinit(cell);

            const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values();

            const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                               n_q_points    = fe_values.n_quadrature_points;
            rhs_values.resize (n_q_points);
            dofs.resize (dofs_per_cell);
            cell_vector.reinit (dofs_per_cell);

            const std::vector<double> &weights   = fe_values.get_JxW_values ();
            rhs_function.value_list (fe_values.get_quadrature_points(),
                                     rhs_values);

            cell_vector = 0;
            for (unsigned int point=0; point<n_q_points; ++point)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                cell_vector(i) += rhs_values[point] *
                                  fe_values.shape_value(i,point) *
                                  weights[point];

            cell->get_dof_indices (dofs);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              rhs_vector(dofs[i]) += cell_vector(i);
          }

      }
    else
      {
        std::vector<Vector<double> > rhs_values;

        for (; cell!=endc; ++cell)
          {
            x_fe_values.reinit(cell);

            const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values();

            const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                               n_q_points    = fe_values.n_quadrature_points;
            rhs_values.resize (n_q_points,
                               Vector<double>(n_components));
            dofs.resize (dofs_per_cell);
            cell_vector.reinit (dofs_per_cell);

            const std::vector<double> &weights   = fe_values.get_JxW_values ();
            rhs_function.vector_value_list (fe_values.get_quadrature_points(),
                                            rhs_values);

            cell_vector = 0;

            // Use the faster code if the
            // FiniteElement is primitive
            if (cell->get_fe().is_primitive ())
              {
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                      const unsigned int component
                        = cell->get_fe().system_to_component_index(i).first;

                      cell_vector(i) += rhs_values[point](component) *
                                        fe_values.shape_value(i,point) *
                                        weights[point];
                    }
              }
            else
              {
                // Otherwise do it the way proposed
                // for vector valued elements
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                      if (cell->get_fe().get_nonzero_components(i)[comp_i])
                        {
                          cell_vector(i) += rhs_values[point](comp_i) *
                                            fe_values.shape_value_component(i,point,comp_i) *
                                            weights[point];
                        }
              }

            cell->get_dof_indices (dofs);

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              rhs_vector(dofs[i]) += cell_vector(i);
          }
      }
  }



  template <int dim, int spacedim>
  void create_right_hand_side (const hp::DoFHandler<dim,spacedim>    &dof_handler,
                               const hp::QCollection<dim>    &quadrature,
                               const Function<spacedim>      &rhs_function,
                               Vector<double>           &rhs_vector)
  {
    create_right_hand_side(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                           dof_handler, quadrature,
                           rhs_function, rhs_vector);
  }




  template <int dim, int spacedim>
  void create_point_source_vector (const Mapping<dim, spacedim>       &mapping,
                                   const DoFHandler<dim,spacedim>    &dof_handler,
                                   const Point<spacedim>         &p,
                                   Vector<double>           &rhs_vector)
  {
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert (dof_handler.get_fe().n_components() == 1,
            ExcMessage ("This function only works for scalar finite elements"));

    rhs_vector = 0;

    std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

    Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim,spacedim> fe_values(mapping, dof_handler.get_fe(),
                                     q, UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    cell_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  fe_values.shape_value(i,0);
  }



  template <int dim, int spacedim>
  void create_point_source_vector (const DoFHandler<dim,spacedim>    &dof_handler,
                                   const Point<spacedim>         &p,
                                   Vector<double>           &rhs_vector)
  {
    create_point_source_vector(StaticMappingQ1<dim,spacedim>::mapping, dof_handler,
                               p, rhs_vector);
  }


  template <int dim, int spacedim>
  void create_point_source_vector (const hp::MappingCollection<dim,spacedim>       &mapping,
                                   const hp::DoFHandler<dim,spacedim>    &dof_handler,
                                   const Point<spacedim>         &p,
                                   Vector<double>           &rhs_vector)
  {
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert (dof_handler.get_fe().n_components() == 1,
            ExcMessage ("This function only works for scalar finite elements"));

    rhs_vector = 0;

    std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

    Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
                            cell_point.first->get_fe(), q, UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    cell_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  fe_values.shape_value(i,0);
  }



  template <int dim, int spacedim>
  void create_point_source_vector (const hp::DoFHandler<dim,spacedim>    &dof_handler,
                                   const Point<spacedim>         &p,
                                   Vector<double>           &rhs_vector)
  {
    create_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
                               dof_handler,
                               p, rhs_vector);
  }




  template <int dim, int spacedim>
  void create_point_source_vector (const Mapping<dim, spacedim>       &mapping,
                                   const DoFHandler<dim,spacedim>     &dof_handler,
                                   const Point<spacedim>              &p,
                                   const Point<dim>                   &orientation,
                                   Vector<double>                     &rhs_vector)
  {
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert (dof_handler.get_fe().n_components() == dim,
            ExcMessage ("This function only works for vector-valued finite elements."));

    rhs_vector = 0;

    const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

    const Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    const FEValuesExtractors::Vector vec (0);
    FEValues<dim,spacedim> fe_values(mapping, dof_handler.get_fe(),
                                     q, UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    cell_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  orientation * fe_values[vec].value(i,0);
  }



  template <int dim, int spacedim>
  void create_point_source_vector (const DoFHandler<dim,spacedim>    &dof_handler,
                                   const Point<spacedim>             &p,
                                   const Point<dim>                  &orientation,
                                   Vector<double>                    &rhs_vector)
  {
    create_point_source_vector(StaticMappingQ1<dim,spacedim>::mapping, dof_handler,
                               p, orientation, rhs_vector);
  }


  template <int dim, int spacedim>
  void create_point_source_vector (const hp::MappingCollection<dim,spacedim> &mapping,
                                   const hp::DoFHandler<dim,spacedim>        &dof_handler,
                                   const Point<spacedim>                     &p,
                                   const Point<dim>                          &orientation,
                                   Vector<double>                            &rhs_vector)
  {
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert (dof_handler.get_fe().n_components() == dim,
            ExcMessage ("This function only works for vector-valued finite elements."));

    rhs_vector = 0;

    std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

    Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    const FEValuesExtractors::Vector vec (0);
    FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
                            cell_point.first->get_fe(), q, UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    cell_point.first->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  orientation * fe_values[vec].value(i,0);
  }



  template <int dim, int spacedim>
  void create_point_source_vector (const hp::DoFHandler<dim,spacedim>   &dof_handler,
                                   const Point<spacedim>                &p,
                                   const Point<dim>                     &orientation,
                                   Vector<double>                       &rhs_vector)
  {
    create_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
                               dof_handler,
                               p, orientation, rhs_vector);
  }



  template <int dim, int spacedim>
  void
  create_boundary_right_hand_side (const Mapping<dim, spacedim>      &mapping,
                                   const DoFHandler<dim,spacedim>   &dof_handler,
                                   const Quadrature<dim-1> &quadrature,
                                   const Function<spacedim>     &rhs_function,
                                   Vector<double>          &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_indicators)
  {
    const FiniteElement<dim> &fe  = dof_handler.get_fe();
    Assert (fe.n_components() == rhs_function.n_components,
            ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));

    rhs_vector = 0;

    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    FEFaceValues<dim> fe_values (mapping, fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points,
                       n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs (dofs_per_cell);
    Vector<double> cell_vector (dofs_per_cell);

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                            endc = dof_handler.end();

    if (n_components==1)
      {
        std::vector<double> rhs_values(n_q_points);

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_indicators.empty() ||
                 (boundary_indicators.find (cell->face(face)->boundary_indicator())
                  !=
                  boundary_indicators.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_vector(i) += rhs_values[point] *
                                      fe_values.shape_value(i,point) *
                                      weights[point];

                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
    else
      {
        std::vector<Vector<double> > rhs_values(n_q_points, Vector<double>(n_components));

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_indicators.empty() ||
                 (boundary_indicators.find (cell->face(face)->boundary_indicator())
                  !=
                  boundary_indicators.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;

                // Use the faster code if the
                // FiniteElement is primitive
                if (fe.is_primitive ())
                  {
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        {
                          const unsigned int component
                            = fe.system_to_component_index(i).first;

                          cell_vector(i) += rhs_values[point](component) *
                                            fe_values.shape_value(i,point) *
                                            weights[point];
                        }
                  }
                else
                  {
                    // And the full featured
                    // code, if vector valued
                    // FEs are used
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                          if (fe.get_nonzero_components(i)[comp_i])
                            {
                              cell_vector(i)
                              += rhs_values[point](comp_i) *
                                 fe_values.shape_value_component(i,point,comp_i) *
                                 weights[point];
                            }
                  }

                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
  }



  template <int dim, int spacedim>
  void
  create_boundary_right_hand_side (const DoFHandler<dim,spacedim>   &dof_handler,
                                   const Quadrature<dim-1> &quadrature,
                                   const Function<spacedim>     &rhs_function,
                                   Vector<double>          &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_indicators)
  {
    create_boundary_right_hand_side(StaticMappingQ1<dim>::mapping, dof_handler,
                                    quadrature,
                                    rhs_function, rhs_vector,
                                    boundary_indicators);
  }



  template <int dim, int spacedim>
  void
  create_boundary_right_hand_side (const hp::MappingCollection<dim,spacedim> &mapping,
                                   const hp::DoFHandler<dim,spacedim> &dof_handler,
                                   const hp::QCollection<dim-1>  &quadrature,
                                   const Function<spacedim>      &rhs_function,
                                   Vector<double>                &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_indicators)
  {
    const hp::FECollection<dim> &fe  = dof_handler.get_fe();
    Assert (fe.n_components() == rhs_function.n_components,
            ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));

    rhs_vector = 0;

    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    hp::FEFaceValues<dim> x_fe_values (mapping, fe, quadrature, update_flags);

    const unsigned int n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs (fe.max_dofs_per_cell());
    Vector<double> cell_vector (fe.max_dofs_per_cell());

    typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    if (n_components==1)
      {
        std::vector<double> rhs_values;

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_indicators.empty() ||
                 (boundary_indicators.find (cell->face(face)->boundary_indicator())
                  !=
                  boundary_indicators.end())))
              {
                x_fe_values.reinit(cell, face);

                const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

                const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                   n_q_points    = fe_values.n_quadrature_points;
                rhs_values.resize (n_q_points);

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_vector(i) += rhs_values[point] *
                                      fe_values.shape_value(i,point) *
                                      weights[point];

                dofs.resize(dofs_per_cell);
                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
    else
      {
        std::vector<Vector<double> > rhs_values;

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_indicators.empty() ||
                 (boundary_indicators.find (cell->face(face)->boundary_indicator())
                  !=
                  boundary_indicators.end())))
              {
                x_fe_values.reinit(cell, face);

                const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

                const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                   n_q_points    = fe_values.n_quadrature_points;
                rhs_values.resize (n_q_points, Vector<double>(n_components));

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;

                // Use the faster code if the
                // FiniteElement is primitive
                if (cell->get_fe().is_primitive ())
                  {
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        {
                          const unsigned int component
                            = cell->get_fe().system_to_component_index(i).first;

                          cell_vector(i) += rhs_values[point](component) *
                                            fe_values.shape_value(i,point) *
                                            weights[point];
                        }
                  }
                else
                  {
                    // And the full featured
                    // code, if vector valued
                    // FEs are used
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                          if (cell->get_fe().get_nonzero_components(i)[comp_i])
                            {
                              cell_vector(i)
                              += rhs_values[point](comp_i) *
                                 fe_values.shape_value_component(i,point,comp_i) *
                                 weights[point];
                            }
                  }
                dofs.resize(dofs_per_cell);
                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
  }



  template <int dim, int spacedim>
  void
  create_boundary_right_hand_side (const hp::DoFHandler<dim,spacedim> &dof_handler,
                                   const hp::QCollection<dim-1>  &quadrature,
                                   const Function<spacedim>      &rhs_function,
                                   Vector<double>                &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_indicators)
  {
    create_boundary_right_hand_side(hp::StaticMappingQ1<dim>::mapping_collection,
                                    dof_handler, quadrature,
                                    rhs_function, rhs_vector,
                                    boundary_indicators);
  }



// ----------- interpolate_boundary_values for std::map --------------------

  namespace
  {
    // interpolate boundary values in
    // 1d. in higher dimensions, we
    // use FEValues to figure out
    // what to do on faces, but in 1d
    // faces are points and it is far
    // easier to simply work on
    // individual vertices
    template <class DH,
              template <int,int> class M_or_MC>
    static inline
    void do_interpolate_boundary_values (const M_or_MC<DH::dimension, DH::space_dimension> &,
                                         const DH                 &dof,
                                         const typename FunctionMap<DH::space_dimension>::type &function_map,
                                         std::map<types::global_dof_index,double> &boundary_values,
                                         const ComponentMask       &component_mask,
                                         const dealii::internal::int2type<1>)
    {
      const unsigned int dim = DH::dimension;
      const unsigned int spacedim=DH::space_dimension;

      Assert (component_mask.represents_n_components(dof.get_fe().n_components()),
              ExcMessage ("The number of components in the mask has to be either "
                          "zero or equal to the number of components in the finite "
                          "element."));

      // if for whatever reason we were
      // passed an empty map, return
      // immediately
      if (function_map.size() == 0)
        return;

      for (typename DH::active_cell_iterator cell = dof.begin_active();
           cell != dof.end(); ++cell)
        for (unsigned int direction=0;
             direction<GeometryInfo<dim>::faces_per_cell; ++direction)
          if (cell->at_boundary(direction)
              &&
              (function_map.find(cell->face(direction)->boundary_indicator()) != function_map.end()))
            {
              const Function<DH::space_dimension> &boundary_function
                = *function_map.find(cell->face(direction)->boundary_indicator())->second;

              // get the FE corresponding to this
              // cell
              const FiniteElement<dim,spacedim> &fe = cell->get_fe();
              Assert (fe.n_components() == boundary_function.n_components,
                      ExcDimensionMismatch(fe.n_components(),
                                           boundary_function.n_components));

              Assert (component_mask.n_selected_components(fe.n_components()) > 0,
                      ComponentMask::ExcNoComponentSelected());

              // now set the value of
              // the vertex degree of
              // freedom. setting
              // also creates the
              // entry in the map if
              // it did not exist
              // beforehand
              //
              // save some time by
              // requesting values
              // only once for each
              // point, irrespective
              // of the number of
              // components of the
              // function
              dealii::Vector<double> function_values (fe.n_components());
              if (fe.n_components() == 1)
                function_values(0)
                  = boundary_function.value (cell->vertex(direction));
              else
                boundary_function.vector_value (cell->vertex(direction),
                                                function_values);

              for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
                if (component_mask[fe.face_system_to_component_index(i).first])
                  boundary_values[cell->
                                  vertex_dof_index(direction,i,
                                                   cell->active_fe_index())]
                    = function_values(fe.face_system_to_component_index(i).first);
            }
    }



    // template for the case dim!=1. Since the function has a template argument
    // dim_, it is clearly less specialized than the 1d function above and
    // whenever possible (i.e., if dim==1), the function template above
    // will be used
    template <class DH,
              template <int,int> class M_or_MC,
              int dim_>
    static inline
    void
    do_interpolate_boundary_values (const M_or_MC<DH::dimension, DH::space_dimension> &mapping,
                                    const DH                 &dof,
                                    const typename FunctionMap<DH::space_dimension>::type &function_map,
                                    std::map<types::global_dof_index,double> &boundary_values,
                                    const ComponentMask       &component_mask,
                                    const dealii::internal::int2type<dim_>)
    {
      const unsigned int dim = DH::dimension;
      const unsigned int spacedim=DH::space_dimension;

      Assert (component_mask.represents_n_components(dof.get_fe().n_components()),
              ExcMessage ("The number of components in the mask has to be either "
                          "zero or equal to the number of components in the finite "
                          "element."));


      // if for whatever reason we were passed an empty map, return
      // immediately
      if (function_map.size() == 0)
        return;

      Assert (function_map.find(numbers::internal_face_boundary_id) == function_map.end(),
              ExcInvalidBoundaryIndicator());

      const unsigned int        n_components = DoFTools::n_components(dof);
      const bool                fe_is_system = (n_components != 1);

      for (typename FunctionMap<spacedim>::type::const_iterator i=function_map.begin();
           i!=function_map.end(); ++i)
        Assert (n_components == i->second->n_components,
                ExcDimensionMismatch(n_components, i->second->n_components));

      // field to store the indices
      std::vector<types::global_dof_index> face_dofs;
      face_dofs.reserve (DoFTools::max_dofs_per_face(dof));

      std::vector<Point<spacedim> >  dof_locations;
      dof_locations.reserve (DoFTools::max_dofs_per_face(dof));

      // array to store the values of the boundary function at the boundary
      // points. have two arrays for scalar and vector functions to use the
      // more efficient one respectively
      std::vector<double>          dof_values_scalar;
      std::vector<dealii::Vector<double> > dof_values_system;
      dof_values_scalar.reserve (DoFTools::max_dofs_per_face (dof));
      dof_values_system.reserve (DoFTools::max_dofs_per_face (dof));

      // before we start with the loop over all cells create an hp::FEValues
      // object that holds the interpolation points of all finite elements
      // that may ever be in use
      dealii::hp::FECollection<dim,spacedim> finite_elements (dof.get_fe());
      dealii::hp::QCollection<dim-1>  q_collection;
      for (unsigned int f=0; f<finite_elements.size(); ++f)
        {
          const FiniteElement<dim,spacedim> &fe = finite_elements[f];

          // generate a quadrature rule on the face from the unit support
          // points. this will be used to obtain the quadrature points on the
          // real cell's face
          //
          // to do this, we check whether the FE has support points on the
          // face at all:
          if (fe.has_face_support_points())
            q_collection.push_back (Quadrature<dim-1>(fe.get_unit_face_support_points()));
          else
            {
              // if not, then we should try a more clever way. the idea is
              // that a finite element may not offer support points for all
              // its shape functions, but maybe only some. if it offers
              // support points for the components we are interested in in
              // this function, then that's fine. if not, the function we call
              // in the finite element will raise an exception. the support
              // points for the other shape functions are left uninitialized
              // (well, initialized by the default constructor), since we
              // don't need them anyway.
              //
              // As a detour, we must make sure we only query
              // face_system_to_component_index if the index corresponds to a
              // primitive shape function. since we know that all the
              // components we are interested in are primitive (by the above
              // check), we can safely put such a check in front
              std::vector<Point<dim-1> > unit_support_points (fe.dofs_per_face);

              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                if (fe.is_primitive (fe.face_to_cell_index(i,0)))
                  if (component_mask[fe.face_system_to_component_index(i).first]
                      == true)
                    unit_support_points[i] = fe.unit_face_support_point(i);

              q_collection.push_back (Quadrature<dim-1>(unit_support_points));
            }
        }
      // now that we have a q_collection object with all the right quadrature
      // points, create an hp::FEFaceValues object that we can use to evaluate
      // the boundary values at
      dealii::hp::MappingCollection<dim,spacedim> mapping_collection (mapping);
      dealii::hp::FEFaceValues<dim,spacedim> x_fe_values (mapping_collection, finite_elements, q_collection,
                                                          update_quadrature_points);

      typename DH::active_cell_iterator cell = dof.begin_active(),
                                        endc = dof.end();
      for (; cell!=endc; ++cell)
        if (!cell->is_artificial())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
               ++face_no)
            {
              const FiniteElement<dim,spacedim> &fe = cell->get_fe();

              // we can presently deal only with primitive elements for
              // boundary values. this does not preclude us using
              // non-primitive elements in components that we aren't
              // interested in, however. make sure that all shape functions
              // that are non-zero for the components we are interested in,
              // are in fact primitive
              for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
                {
                  const ComponentMask &nonzero_component_array
                    = cell->get_fe().get_nonzero_components (i);
                  for (unsigned int c=0; c<n_components; ++c)
                    if ((nonzero_component_array[c] == true)
                        &&
                        (component_mask[c] == true))
                      Assert (cell->get_fe().is_primitive (i),
                              ExcMessage ("This function can only deal with requested boundary "
                                          "values that correspond to primitive (scalar) base "
                                          "elements"));
                }

              const typename DH::face_iterator face = cell->face(face_no);
              const types::boundary_id boundary_component = face->boundary_indicator();

              // see if this face is part of the boundaries for which we are
              // supposed to do something, and also see if the finite element
              // in use here has DoFs on the face at all
              if ((function_map.find(boundary_component) != function_map.end())
                  &&
                  (cell->get_fe().dofs_per_face > 0))
                {
                  // face is of the right component
                  x_fe_values.reinit(cell, face_no);
                  const dealii::FEFaceValues<dim,spacedim> &fe_values =
                    x_fe_values.get_present_fe_values();

                  // get indices, physical location and boundary values of
                  // dofs on this face
                  face_dofs.resize (fe.dofs_per_face);
                  face->get_dof_indices (face_dofs, cell->active_fe_index());
                  const std::vector<Point<spacedim> > &dof_locations
                    = fe_values.get_quadrature_points ();

                  if (fe_is_system)
                    {
                      // resize array. avoid construction of a memory
                      // allocating temporary if possible
                      if (dof_values_system.size() < fe.dofs_per_face)
                        dof_values_system.resize (fe.dofs_per_face,
                                                  dealii::Vector<double>(fe.n_components()));
                      else
                        dof_values_system.resize (fe.dofs_per_face);

                      function_map.find(boundary_component)->second
                      ->vector_value_list (dof_locations, dof_values_system);

                      // enter those dofs into the list that match the
                      // component signature. avoid the usual complication
                      // that we can't just use *_system_to_component_index
                      // for non-primitive FEs
                      for (unsigned int i=0; i<face_dofs.size(); ++i)
                        {
                          unsigned int component;
                          if (fe.is_primitive())
                            component = fe.face_system_to_component_index(i).first;
                          else
                            {
                              // non-primitive case. make sure that this
                              // particular shape function _is_ primitive, and
                              // get at it's component. use usual trick to
                              // transfer face dof index to cell dof index
                              const unsigned int cell_i
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
                              Assert (cell_i < fe.dofs_per_cell, ExcInternalError());

                              // make sure that if this is not a primitive
                              // shape function, then all the corresponding
                              // components in the mask are not set
                              if (!fe.is_primitive(cell_i))
                                for (unsigned int c=0; c<n_components; ++c)
                                  if (fe.get_nonzero_components(cell_i)[c])
                                    Assert (component_mask[c] == false,
                                            FETools::ExcFENotPrimitive());

                              // let's pick the first of possibly more than
                              // one non-zero components. if shape function is
                              // non-primitive, then we will ignore the result
                              // in the following anyway, otherwise there's
                              // only one non-zero component which we will use
                              component = fe.get_nonzero_components(cell_i).first_selected_component();
                            }

                          if (component_mask[component] == true)
                            boundary_values[face_dofs[i]] = dof_values_system[i](component);
                        }
                    }
                  else
                    // fe has only one component, so save some computations
                    {
                      // get only the one component that this function has
                      dof_values_scalar.resize (fe.dofs_per_face);
                      function_map.find(boundary_component)->second
                      ->value_list (dof_locations, dof_values_scalar, 0);

                      // enter into list

                      for (unsigned int i=0; i<face_dofs.size(); ++i)
                        boundary_values[face_dofs[i]] = dof_values_scalar[i];
                    }
                }
            }
    } // end of interpolate_boundary_values
  } // end of namespace internal



  template <class DH>
  void

  interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension>            &mapping,
                               const DH                 &dof,
                               const typename FunctionMap<DH::space_dimension>::type &function_map,
                               std::map<types::global_dof_index,double> &boundary_values,
                               const ComponentMask       &component_mask_)
  {
    do_interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                    component_mask_,
                                    dealii::internal::int2type<DH::dimension>());
  }



  template <class DH>
  void
  interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension>            &mapping,
                               const DH                 &dof,
                               const types::boundary_id            boundary_component,
                               const Function<DH::space_dimension>           &boundary_function,
                               std::map<types::global_dof_index,double> &boundary_values,
                               const ComponentMask       &component_mask)
  {
    typename FunctionMap<DH::space_dimension>::type function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                 component_mask);
  }


  template <int dim, int spacedim>
  void
  interpolate_boundary_values (const hp::MappingCollection<dim,spacedim>            &mapping,
                               const hp::DoFHandler<dim,spacedim>                 &dof,
                               const typename FunctionMap<spacedim>::type &function_map,
                               std::map<types::global_dof_index,double> &boundary_values,
                               const ComponentMask       &component_mask_)
  {
    do_interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                    component_mask_,
                                    dealii::internal::int2type<dim>());
  }



  template <class DH>
  void
  interpolate_boundary_values (const DH                 &dof,
                               const types::boundary_id            boundary_component,
                               const Function<DH::space_dimension>           &boundary_function,
                               std::map<types::global_dof_index,double> &boundary_values,
                               const ComponentMask       &component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
                                dof, boundary_component,
                                boundary_function, boundary_values, component_mask);
  }



  template <class DH>
  void
  interpolate_boundary_values (const DH                 &dof,
                               const typename FunctionMap<DH::space_dimension>::type &function_map,
                               std::map<types::global_dof_index,double> &boundary_values,
                               const ComponentMask       &component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
                                dof, function_map,
                                boundary_values, component_mask);
  }




// ----------- interpolate_boundary_values for ConstraintMatrix --------------



  template <class DH>
  void
  interpolate_boundary_values
  (const Mapping<DH::dimension, DH::space_dimension>     &mapping,
   const DH                                              &dof,
   const typename FunctionMap<DH::space_dimension>::type &function_map,
   ConstraintMatrix                                      &constraints,
   const ComponentMask                               &component_mask_)
  {
    std::map<types::global_dof_index,double> boundary_values;
    interpolate_boundary_values (mapping, dof, function_map,
                                 boundary_values, component_mask_);
    std::map<types::global_dof_index,double>::const_iterator boundary_value =
      boundary_values.begin();
    for ( ; boundary_value !=boundary_values.end(); ++boundary_value)
      {
        if (constraints.can_store_line (boundary_value->first)
            &&
            !constraints.is_constrained(boundary_value->first))
          {
            constraints.add_line (boundary_value->first);
            constraints.set_inhomogeneity (boundary_value->first,
                                           boundary_value->second);
          }
      }
  }



  template <class DH>
  void
  interpolate_boundary_values
  (const Mapping<DH::dimension, DH::space_dimension> &mapping,
   const DH                                          &dof,
   const types::boundary_id                                boundary_component,
   const Function<DH::space_dimension>               &boundary_function,
   ConstraintMatrix                                  &constraints,
   const ComponentMask                           &component_mask)
  {
    typename FunctionMap<DH::space_dimension>::type function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values (mapping, dof, function_map, constraints,
                                 component_mask);
  }



  template <class DH>
  void
  interpolate_boundary_values
  (const DH                            &dof,
   const types::boundary_id                  boundary_component,
   const Function<DH::space_dimension> &boundary_function,
   ConstraintMatrix                    &constraints,
   const ComponentMask             &component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
                                dof, boundary_component,
                                boundary_function, constraints, component_mask);
  }



  template <class DH>
  void
  interpolate_boundary_values
  (const DH                                              &dof,
   const typename FunctionMap<DH::space_dimension>::type &function_map,
   ConstraintMatrix                                      &constraints,
   const ComponentMask                               &component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping,
                                dof, function_map,
                                constraints, component_mask);
  }




// -------- implementation for project_boundary_values with std::map --------


  namespace
  {
    template <int dim, int spacedim,
              template <int,int> class DH,
              template <int,int> class M_or_MC,
              template <int> class Q_or_QC>
    void
    do_project_boundary_values (const M_or_MC<dim, spacedim>   &mapping,
                                const DH<dim, spacedim> &dof,
                                const typename FunctionMap<spacedim>::type &boundary_functions,
                                const Q_or_QC<dim-1>        &q,
                                std::map<types::global_dof_index,double>  &boundary_values,
                                std::vector<unsigned int>       component_mapping)
    {
      // in 1d, projection onto the 0d end points == interpolation
      if (dim == 1)
        {
          Assert (component_mapping.size() == 0, ExcNotImplemented());
          interpolate_boundary_values (mapping, dof, boundary_functions,
                                       boundary_values, ComponentMask());
          return;
        }

      //TODO:[?] In project_boundary_values, no condensation of sparsity
      //    structures, matrices and right hand sides or distribution of
      //    solution vectors is performed. This is ok for dim<3 because then
      //    there are no constrained nodes on the boundary, but is not
      //    acceptable for higher dimensions. Fix this.

      if (component_mapping.size() == 0)
        {
          AssertDimension (dof.get_fe().n_components(), boundary_functions.begin()->second->n_components);
          // I still do not see why i
          // should create another copy
          // here
          component_mapping.resize(dof.get_fe().n_components());
          for (unsigned int i= 0 ; i < component_mapping.size() ; ++i)
            component_mapping[i] = i;
        }
      else
        AssertDimension (dof.get_fe().n_components(), component_mapping.size());

      std::vector<types::global_dof_index> dof_to_boundary_mapping;
      std::set<types::boundary_id> selected_boundary_components;
      for (typename FunctionMap<spacedim>::type::const_iterator i=boundary_functions.begin();
           i!=boundary_functions.end(); ++i)
        selected_boundary_components.insert (i->first);

      DoFTools::map_dof_to_boundary_indices (dof, selected_boundary_components,
                                             dof_to_boundary_mapping);

      // Done if no degrees of freedom on the boundary
      if (dof.n_boundary_dofs (boundary_functions) == 0)
        return;

      // set up sparsity structure
      CompressedSparsityPattern c_sparsity(dof.n_boundary_dofs (boundary_functions),
                                           dof.n_boundary_dofs (boundary_functions));
      DoFTools::make_boundary_sparsity_pattern (dof,
                                                boundary_functions,
                                                dof_to_boundary_mapping,
                                                c_sparsity);
      SparsityPattern sparsity;
      sparsity.copy_from(c_sparsity);



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
      if (dim>=3)
        {
#ifdef DEBUG
          // Assert that there are no hanging nodes at the boundary
          int level = -1;
          for (typename DH<dim,spacedim>::active_cell_iterator cell = dof.begin_active();
               cell != dof.end(); ++cell)
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (cell->at_boundary(f))
                  {
                    if (level == -1)
                      level = cell->level();
                    else
                      {
                        Assert (level == cell->level(), ExcNotImplemented());
                      }
                  }
              }
#endif
        }
      sparsity.compress();


      // make mass matrix and right hand side
      SparseMatrix<double> mass_matrix(sparsity);
      Vector<double>       rhs(sparsity.n_rows());


      MatrixCreator::create_boundary_mass_matrix (mapping, dof, q,
                                                  mass_matrix, boundary_functions,
                                                  rhs, dof_to_boundary_mapping, (const Function<spacedim> *) 0,
                                                  component_mapping);

      // For certain weird elements,
      // there might be degrees of
      // freedom on the boundary, but
      // their shape functions do not
      // have support there. Let's
      // eliminate them here.

      // The Bogner-Fox-Schmidt element
      // is an example for those.

//TODO: Maybe we should figure out if the element really needs this

      FilteredMatrix<Vector<double> > filtered_mass_matrix(mass_matrix, true);
      FilteredMatrix<Vector<double> > filtered_precondition;
      std::vector<bool> excluded_dofs(mass_matrix.m(), false);

      double max_element = 0.;
      for (unsigned int i=0; i<mass_matrix.m(); ++i)
        if (mass_matrix.diag_element(i) > max_element)
          max_element = mass_matrix.diag_element(i);

      for (unsigned int i=0; i<mass_matrix.m(); ++i)
        if (mass_matrix.diag_element(i) < 1.e-8 * max_element)
          {
            filtered_mass_matrix.add_constraint(i, 0.);
            filtered_precondition.add_constraint(i, 0.);
            mass_matrix.diag_element(i) = 1.;
            excluded_dofs[i] = true;
          }

      Vector<double> boundary_projection (rhs.size());

      // cannot reduce residual in a useful way if we are close to the square
      // root of the minimal double value
      if (rhs.norm_sqr() < 1e28 * std::numeric_limits<double>::min())
        boundary_projection = 0;
      else
        {
          // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
          // steps may not be sufficient, since roundoff errors may accumulate for
          // badly conditioned matrices
          ReductionControl        control(5*rhs.size(), 0., 1.e-12, false, false);
          GrowingVectorMemory<> memory;
          SolverCG<>              cg(control,memory);

          PreconditionSSOR<> prec;
          prec.initialize(mass_matrix, 1.2);
          filtered_precondition.initialize(prec, true);
          // solve
          cg.solve (filtered_mass_matrix, boundary_projection, rhs, filtered_precondition);
          filtered_precondition.apply_constraints(boundary_projection, true);
          filtered_precondition.clear();
        }
      // fill in boundary values
      for (unsigned int i=0; i<dof_to_boundary_mapping.size(); ++i)
        if (dof_to_boundary_mapping[i] != DoFHandler<dim,spacedim>::invalid_dof_index
            && ! excluded_dofs[dof_to_boundary_mapping[i]])
          {
            Assert(numbers::is_finite(boundary_projection(dof_to_boundary_mapping[i])), ExcNumberNotFinite());

            // this dof is on one of the
            // interesting boundary parts
            //
            // remember: i is the global dof
            // number, dof_to_boundary_mapping[i]
            // is the number on the boundary and
            // thus in the solution vector
            boundary_values[i] = boundary_projection(dof_to_boundary_mapping[i]);
          }
    }
  }

  template <int dim, int spacedim>
  void
  project_boundary_values (const Mapping<dim, spacedim>   &mapping,
                           const DoFHandler<dim, spacedim> &dof,
                           const typename FunctionMap<spacedim>::type &boundary_functions,
                           const Quadrature<dim-1>        &q,
                           std::map<types::global_dof_index,double>  &boundary_values,
                           std::vector<unsigned int>       component_mapping)
  {
    do_project_boundary_values(mapping, dof, boundary_functions, q,
                               boundary_values, component_mapping);
  }



  template <int dim, int spacedim>
  void
  project_boundary_values (const DoFHandler<dim,spacedim>    &dof,
                           const typename FunctionMap<spacedim>::type &boundary_functions,
                           const Quadrature<dim-1>  &q,
                           std::map<types::global_dof_index,double> &boundary_values,
                           std::vector<unsigned int> component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
                            boundary_values, component_mapping);
  }



  template <int dim, int spacedim>
  void project_boundary_values (const hp::MappingCollection<dim, spacedim>       &mapping,
                                const hp::DoFHandler<dim,spacedim>    &dof,
                                const typename FunctionMap<spacedim>::type &boundary_functions,
                                const hp::QCollection<dim-1>  &q,
                                std::map<types::global_dof_index,double> &boundary_values,
                                std::vector<unsigned int> component_mapping)
  {
    do_project_boundary_values (mapping, dof,
                                boundary_functions,
                                q, boundary_values,
                                component_mapping);
  }



  template <int dim, int spacedim>
  void project_boundary_values (const hp::DoFHandler<dim,spacedim>    &dof,
                                const typename FunctionMap<spacedim>::type &boundary_function,
                                const hp::QCollection<dim-1>  &q,
                                std::map<types::global_dof_index,double> &boundary_values,
                                std::vector<unsigned int> component_mapping)
  {
    project_boundary_values (hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof,
                             boundary_function,
                             q, boundary_values,
                             component_mapping);
  }


  // ----- implementation for project_boundary_values with ConstraintMatrix -----



  template <int dim, int spacedim>
  void
  project_boundary_values (const Mapping<dim, spacedim>       &mapping,
                           const DoFHandler<dim,spacedim>    &dof,
                           const typename FunctionMap<spacedim>::type &boundary_functions,
                           const Quadrature<dim-1>  &q,
                           ConstraintMatrix &constraints,
                           std::vector<unsigned int> component_mapping)
  {
    std::map<types::global_dof_index,double> boundary_values;
    project_boundary_values (mapping, dof, boundary_functions, q,
                             boundary_values, component_mapping);
    std::map<types::global_dof_index,double>::const_iterator boundary_value =
      boundary_values.begin();
    for ( ; boundary_value !=boundary_values.end(); ++boundary_value)
      {
        if (!constraints.is_constrained(boundary_value->first))
          {
            constraints.add_line (boundary_value->first);
            constraints.set_inhomogeneity (boundary_value->first,
                                           boundary_value->second);
          }
      }
  }



  template <int dim, int spacedim>
  void
  project_boundary_values (const DoFHandler<dim,spacedim>    &dof,
                           const typename FunctionMap<spacedim>::type &boundary_functions,
                           const Quadrature<dim-1>  &q,
                           ConstraintMatrix &constraints,
                           std::vector<unsigned int> component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
                            constraints, component_mapping);
  }




  namespace internal
  {
    /**
     * A structure that stores the dim DoF
     * indices that correspond to a
     * vector-valued quantity at a single
     * support point.
     */
    template <int dim>
    struct VectorDoFTuple
    {
      types::global_dof_index dof_indices[dim];

      VectorDoFTuple ()
      {
        for (unsigned int i=0; i<dim; ++i)
          dof_indices[i] = numbers::invalid_dof_index;
      }


      bool operator < (const VectorDoFTuple<dim> &other) const
      {
        for (unsigned int i=0; i<dim; ++i)
          if (dof_indices[i] < other.dof_indices[i])
            return true;
          else if (dof_indices[i] > other.dof_indices[i])
            return false;
        return false;
      }

      bool operator == (const VectorDoFTuple<dim> &other) const
      {
        for (unsigned int i=0; i<dim; ++i)
          if (dof_indices[i] != other.dof_indices[i])
            return false;

        return true;
      }

      bool operator != (const VectorDoFTuple<dim> &other) const
      {
        return ! (*this == other);
      }
    };


    template <int dim>
    std::ostream &operator << (std::ostream &out,
                               const VectorDoFTuple<dim> &vdt)
    {
      for (unsigned int d=0; d<dim; ++d)
        out << vdt.dof_indices[d] << (d < dim-1 ? " " : "");
      return out;
    }



    /**
     * Add the constraint
     * $\vec n \cdot \vec u = inhom$
     * to the list of constraints.
     *
     * Here, $\vec u$ is represented
     * by the set of given DoF
     * indices, and $\vec n$ by the
     * vector specified as the second
     * argument.
     *
     * The function does not add constraints
     * if a degree of freedom is already
     * constrained in the constraints object.
     */
    template <int dim>
    void
    add_constraint (const VectorDoFTuple<dim> &dof_indices,
                    const Tensor<1,dim>       &constraining_vector,
                    ConstraintMatrix          &constraints,
                    const double              inhomogeneity=0)
    {

      // choose the DoF that has the
      // largest component in the
      // constraining_vector as the
      // one to be constrained as
      // this makes the process
      // stable in cases where the
      // constraining_vector has the
      // form n=(1,0) or n=(0,1)
      //
      // we get constraints of the form
      //   x0 = a_1*x1 + a2*x2 + ...
      // if one of the weights is
      // essentially zero then skip
      // this part. the ConstraintMatrix
      // can also deal with cases like
      //   x0 = 0
      // if necessary
      //
      // there is a problem if we have a
      // normal vector of the form
      // (a,a,small) or (a,a,a). Depending on
      // round-off we may choose the first or
      // second component (or third, in the
      // latter case) as the largest one, and
      // depending on our choice one or
      // another degree of freedom will be
      // constrained. On a single processor
      // this is not much of a problem, but
      // it's a nightmare when we run in
      // parallel and two processors disagree
      // on which DoF should be
      // constrained. This led to an
      // incredibly difficult to find bug in
      // step-32 when running in parallel
      // with 9 or more processors.
      //
      // in practice, such normal vectors of
      // the form (a,a,small) or (a,a,a)
      // happen not infrequently since they
      // lie on the diagonals where vertices
      // frequently happen to land upon mesh
      // refinement if one starts from a
      // symmetric and regular body. we work
      // around this problem in the following
      // way: if we have a normal vector of
      // the form (a,b) (similarly algorithm
      // in 3d), we choose 'a' as the largest
      // coefficient not if a>b but if
      // a>b+1e-10. this shifts the problem
      // away from the frequently visited
      // diagonal to a line that is off the
      // diagonal. there will of course be
      // problems where the exact values of a
      // and b differ by exactly 1e-10 and we
      // get into the same instability, but
      // from a practical viewpoint such
      // problems should be much rarer. in
      // particular, meshes have to be very
      // very fine for a vertex to land on
      // this line if the original body had a
      // vertex on the diagonal as well
      switch (dim)
        {
        case 2:
        {
          if (std::fabs(constraining_vector[0]) > std::fabs(constraining_vector[1]) + 1e-10)
            {
              if (!constraints.is_constrained(dof_indices.dof_indices[0])
                  &&
                  constraints.can_store_line(dof_indices.dof_indices[0]))
                {
                  constraints.add_line (dof_indices.dof_indices[0]);

                  if (std::fabs (constraining_vector[1]/constraining_vector[0])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[0],
                                           dof_indices.dof_indices[1],
                                           -constraining_vector[1]/constraining_vector[0]);

                  if (std::fabs (inhomogeneity/constraining_vector[0])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(dof_indices.dof_indices[0],
                                                  inhomogeneity/constraining_vector[0]);
                }
            }
          else
            {
              if (!constraints.is_constrained(dof_indices.dof_indices[1])
                  &&
                  constraints.can_store_line(dof_indices.dof_indices[1]))
                {
                  constraints.add_line (dof_indices.dof_indices[1]);

                  if (std::fabs (constraining_vector[0]/constraining_vector[1])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[1],
                                           dof_indices.dof_indices[0],
                                           -constraining_vector[0]/constraining_vector[1]);

                  if (std::fabs (inhomogeneity/constraining_vector[1])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(dof_indices.dof_indices[1],
                                                  inhomogeneity/constraining_vector[1]);
                }
            }
          break;
        }

        case 3:
        {
          if ((std::fabs(constraining_vector[0]) >= std::fabs(constraining_vector[1])+1e-10)
              &&
              (std::fabs(constraining_vector[0]) >= std::fabs(constraining_vector[2])+2e-10))
            {
              if (!constraints.is_constrained(dof_indices.dof_indices[0])
                  &&
                  constraints.can_store_line(dof_indices.dof_indices[0]))
                {
                  constraints.add_line (dof_indices.dof_indices[0]);

                  if (std::fabs (constraining_vector[1]/constraining_vector[0])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[0],
                                           dof_indices.dof_indices[1],
                                           -constraining_vector[1]/constraining_vector[0]);

                  if (std::fabs (constraining_vector[2]/constraining_vector[0])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[0],
                                           dof_indices.dof_indices[2],
                                           -constraining_vector[2]/constraining_vector[0]);

                  if (std::fabs (inhomogeneity/constraining_vector[0])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(dof_indices.dof_indices[0],
                                                  inhomogeneity/constraining_vector[0]);
                }
            }
          else if ((std::fabs(constraining_vector[1])+1e-10 >= std::fabs(constraining_vector[0]))
                   &&
                   (std::fabs(constraining_vector[1]) >= std::fabs(constraining_vector[2])+1e-10))
            {
              if (!constraints.is_constrained(dof_indices.dof_indices[1])
                  &&
                  constraints.can_store_line(dof_indices.dof_indices[1]))
                {
                  constraints.add_line (dof_indices.dof_indices[1]);

                  if (std::fabs (constraining_vector[0]/constraining_vector[1])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[1],
                                           dof_indices.dof_indices[0],
                                           -constraining_vector[0]/constraining_vector[1]);

                  if (std::fabs (constraining_vector[2]/constraining_vector[1])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[1],
                                           dof_indices.dof_indices[2],
                                           -constraining_vector[2]/constraining_vector[1]);

                  if (std::fabs (inhomogeneity/constraining_vector[1])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(dof_indices.dof_indices[1],
                                                  inhomogeneity/constraining_vector[1]);
                }
            }
          else
            {
              if (!constraints.is_constrained(dof_indices.dof_indices[2])
                  &&
                  constraints.can_store_line(dof_indices.dof_indices[2]))
                {
                  constraints.add_line (dof_indices.dof_indices[2]);

                  if (std::fabs (constraining_vector[0]/constraining_vector[2])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[2],
                                           dof_indices.dof_indices[0],
                                           -constraining_vector[0]/constraining_vector[2]);

                  if (std::fabs (constraining_vector[1]/constraining_vector[2])
                      > std::numeric_limits<double>::epsilon())
                    constraints.add_entry (dof_indices.dof_indices[2],
                                           dof_indices.dof_indices[1],
                                           -constraining_vector[1]/constraining_vector[2]);

                  if (std::fabs (inhomogeneity/constraining_vector[2])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(dof_indices.dof_indices[2],
                                                  inhomogeneity/constraining_vector[2]);
                }
            }

          break;
        }

        default:
          Assert (false, ExcNotImplemented());
        }
    }


    /**
     * Add the constraint $(\vec u-\vec u_\Gamma) \|
     * \vec t$ to the list of
     * constraints. In 2d, this is a
     * single constraint, in 3d these
     * are two constraints
     *
     * Here, $\vec u$ is represented
     * by the set of given DoF
     * indices, and $\vec t$ by the
     * vector specified as the second
     * argument.
     *
     * The function does not add constraints
     * if a degree of freedom is already
     * constrained in the constraints object.
     */
    template <int dim>
    void
    add_tangentiality_constraints
    (const VectorDoFTuple<dim> &dof_indices,
     const Tensor<1,dim>       &tangent_vector,
     ConstraintMatrix          &constraints,
     const Vector<double>      &b_values = Vector<double>(dim))
    {

      // choose the DoF that has the
      // largest component in the
      // tangent_vector as the
      // independent component, and
      // then constrain the others to
      // it. specifically, if, say,
      // component 0 of the tangent
      // vector t is largest by
      // magnitude, then
      // x1=(b[1]*t[0]-b[0]*t[1])/t[0]+t[1]/t[0]*x_0, etc.
      unsigned int largest_component = 0;
      for (unsigned int d=1; d<dim; ++d)
        if (std::fabs(tangent_vector[d]) > std::fabs(tangent_vector[largest_component]) + 1e-10)
          largest_component = d;

      // then constrain all of the
      // other degrees of freedom in
      // terms of the one just found
      for (unsigned int d=0; d<dim; ++d)
        if (d != largest_component)
          if (!constraints.is_constrained(dof_indices.dof_indices[d])
              &&
              constraints.can_store_line(dof_indices.dof_indices[d]))
            {
              constraints.add_line (dof_indices.dof_indices[d]);

              if (std::fabs (tangent_vector[d]/tangent_vector[largest_component])
                  > std::numeric_limits<double>::epsilon())
                constraints.add_entry (dof_indices.dof_indices[d],
                                       dof_indices.dof_indices[largest_component],
                                       tangent_vector[d]/tangent_vector[largest_component]);

              const double inhomogeneity
                = (b_values(d)*tangent_vector[largest_component]
                   -b_values(largest_component)*tangent_vector[d])
                  /tangent_vector[largest_component];

              if (std::fabs(inhomogeneity)
                  > std::numeric_limits<double>::epsilon())
                constraints.set_inhomogeneity(dof_indices.dof_indices[d],
                                              inhomogeneity);
            }
    }



    /**
     * Given a vector, compute a set
     * of dim-1 vectors that are
     * orthogonal to the first one
     * and mutually orthonormal as
     * well.
     */
    template <int dim>
    void
    compute_orthonormal_vectors (const Tensor<1,dim> &vector,
                                 Tensor<1,dim> (&orthonormals)[dim-1])
    {
      switch (dim)
        {
        case 3:
        {
          // to do this in 3d, take
          // one vector that is
          // guaranteed to be not
          // aligned with the
          // average tangent and
          // form the cross
          // product. this yields
          // one vector that is
          // certainly
          // perpendicular to the
          // tangent; then take the
          // cross product between
          // this vector and the
          // tangent and get one
          // vector that is
          // perpendicular to both

          // construct a
          // temporary vector
          // by swapping the
          // larger two
          // components and
          // flipping one
          // sign; this can
          // not be collinear
          // with the average
          // tangent
          Tensor<1,dim> tmp = vector;
          if ((std::fabs(tmp[0]) > std::fabs(tmp[1]))
              &&
              (std::fabs(tmp[0]) > std::fabs(tmp[2])))
            {
              // entry zero
              // is the
              // largest
              if ((std::fabs(tmp[1]) > std::fabs(tmp[2])))
                std::swap (tmp[0], tmp[1]);
              else
                std::swap (tmp[0], tmp[2]);

              tmp[0] *= -1;
            }
          else if ((std::fabs(tmp[1]) > std::fabs(tmp[0]))
                   &&
                   (std::fabs(tmp[1]) > std::fabs(tmp[2])))
            {
              // entry one
              // is the
              // largest
              if ((std::fabs(tmp[0]) > std::fabs(tmp[2])))
                std::swap (tmp[1], tmp[0]);
              else
                std::swap (tmp[1], tmp[2]);

              tmp[1] *= -1;
            }
          else
            {
              // entry two
              // is the
              // largest
              if ((std::fabs(tmp[0]) > std::fabs(tmp[1])))
                std::swap (tmp[2], tmp[0]);
              else
                std::swap (tmp[2], tmp[1]);

              tmp[2] *= -1;
            }

          // make sure the two vectors
          // are indeed not collinear
          Assert (std::fabs(vector * tmp / vector.norm() / tmp.norm())
                  <
                  (1-1e-12),
                  ExcInternalError());

          // now compute the
          // two normals
          cross_product (orthonormals[0], vector, tmp);
          cross_product (orthonormals[1], vector, orthonormals[0]);

          break;
        }

        default:
          Assert (false, ExcNotImplemented());
        }
    }
  }


  namespace internals
  {
    // This function computes the
    // projection of the boundary
    // function on edges for 3D.
    template<typename cell_iterator>
    void
    compute_edge_projection (const cell_iterator &cell,
                             const unsigned int face,
                             const unsigned int line,
                             hp::FEValues<3> &hp_fe_values,
                             const Function<3> &boundary_function,
                             const unsigned int first_vector_component,
                             std::vector<double> &dof_values,
                             std::vector<bool> &dofs_processed)
    {
      const double tol = 0.5 * cell->face (face)->line (line)->diameter () / cell->get_fe ().degree;
      const unsigned int dim = 3;
      const unsigned int spacedim = 3;

      hp_fe_values.reinit
      (cell,
       (cell->active_fe_index () * GeometryInfo<dim>::faces_per_cell + face)
       * GeometryInfo<dim>::lines_per_face + line);

      // Initialize the required
      // objects.
      const FEValues<dim> &
      fe_values = hp_fe_values.get_present_fe_values ();
      const FiniteElement<dim> &fe = cell->get_fe ();
      const std::vector< DerivativeForm<1,dim,spacedim> > &
      jacobians = fe_values.get_jacobians ();
      const std::vector<Point<dim> > &
      quadrature_points = fe_values.get_quadrature_points ();

      std::vector<Point<dim> > tangentials (fe_values.n_quadrature_points);
      std::vector<Vector<double> > values (fe_values.n_quadrature_points,
                                           Vector<double> (fe.n_components ()));

      // Get boundary function values
      // at quadrature points.
      boundary_function.vector_value_list (quadrature_points, values);

      const std::vector<Point<dim> > &
      reference_quadrature_points = fe_values.get_quadrature ().get_points ();
      std::pair<unsigned int, unsigned int> base_indices (0, 0);

      if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) != 0)
        {
          unsigned int fe_index = 0;
          unsigned int fe_index_old = 0;
          unsigned int i = 0;

          for (; i < fe.n_base_elements (); ++i)
            {
              fe_index_old = fe_index;
              fe_index += fe.element_multiplicity (i) * fe.base_element (i).n_components ();

              if (fe_index > first_vector_component)
                break;
            }

          base_indices.first = i;
          base_indices.second = (first_vector_component - fe_index_old) / fe.base_element (i).n_components ();
        }

      // coordinate directions of
      // the edges of the face.
      const unsigned int
      edge_coordinate_direction
      [GeometryInfo<dim>::faces_per_cell]
      [GeometryInfo<dim>::lines_per_face]
      = { { 2, 2, 1, 1 },
        { 2, 2, 1, 1 },
        { 0, 0, 2, 2 },
        { 0, 0, 2, 2 },
        { 1, 1, 0, 0 },
        { 1, 1, 0, 0 }
      };
      const FEValuesExtractors::Vector vec (first_vector_component);

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
          Point<dim> shifted_reference_point_1 = reference_quadrature_points[q_point];
          Point<dim> shifted_reference_point_2 = reference_quadrature_points[q_point];

          shifted_reference_point_1 (edge_coordinate_direction[face][line]) += tol;
          shifted_reference_point_2 (edge_coordinate_direction[face][line]) -= tol;
          tangentials[q_point]
            = (0.5 *
               (fe_values.get_mapping ()
                .transform_unit_to_real_cell (cell,
                                              shifted_reference_point_1)
                -
                fe_values.get_mapping ()
                .transform_unit_to_real_cell (cell,
                                              shifted_reference_point_2))
               / tol);
          tangentials[q_point]
          /= std::sqrt (tangentials[q_point].square ());

          // Compute the degrees of
          // freedom.
          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                 && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                 && (fe.base_element (base_indices.first).face_to_cell_index (line * fe.degree, face)
                     <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second)
                 && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                     <= fe.base_element (base_indices.first).face_to_cell_index
                     ((line + 1) * fe.degree - 1, face)))
                || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0) && (line * fe.degree <= i)
                    && (i < (line + 1) * fe.degree)))
              {
                dof_values[i]
                += (fe_values.JxW (q_point)
                    * (values[q_point] (first_vector_component) * tangentials[q_point] (0)
                       + values[q_point] (first_vector_component + 1) * tangentials[q_point] (1)
                       + values[q_point] (first_vector_component + 2) * tangentials[q_point] (2))
                    * (fe_values[vec].value (fe.face_to_cell_index (i, face), q_point) * tangentials[q_point])
                    / std::sqrt (jacobians[q_point][0][edge_coordinate_direction[face][line]]
                                 * jacobians[q_point][0][edge_coordinate_direction[face][line]]
                                 + jacobians[q_point][1][edge_coordinate_direction[face][line]]
                                 * jacobians[q_point][1][edge_coordinate_direction[face][line]]
                                 + jacobians[q_point][2][edge_coordinate_direction[face][line]]
                                 * jacobians[q_point][2][edge_coordinate_direction[face][line]]));

                if (q_point == 0)
                  dofs_processed[i] = true;
              }
        }
    }

    // dummy implementation of above
    // function for all other
    // dimensions
    template<int dim, typename cell_iterator>
    void
    compute_edge_projection (const cell_iterator &,
                             const unsigned int,
                             const unsigned int,
                             hp::FEValues<dim> &,
                             const Function<dim> &,
                             const unsigned int,
                             std::vector<double> &,
                             std::vector<bool> &)
    {
      Assert (false, ExcInternalError ());
    }

    // This function computes the
    // projection of the boundary
    // function on the interior of
    // faces.
    template<int dim, typename cell_iterator>
    void
    compute_face_projection_curl_conforming (const cell_iterator &cell,
                                             const unsigned int face,
                                             hp::FEValues<dim> &hp_fe_values,
                                             const Function<dim> &boundary_function,
                                             const unsigned int first_vector_component,
                                             std::vector<double> &dof_values,
                                             std::vector<bool> &dofs_processed)
    {
      const unsigned int spacedim = dim;
      hp_fe_values.reinit (cell, cell->active_fe_index ()
                           * GeometryInfo<dim>::faces_per_cell + face);
      // Initialize the required
      // objects.
      const FEValues<dim> &
      fe_values = hp_fe_values.get_present_fe_values ();
      const FiniteElement<dim> &fe = cell->get_fe ();
      const std::vector< DerivativeForm<1,dim,spacedim> > &
      jacobians = fe_values.get_jacobians ();
      const std::vector<Point<dim> > &
      quadrature_points = fe_values.get_quadrature_points ();
      const unsigned int degree = fe.degree - 1;
      std::pair<unsigned int, unsigned int> base_indices (0, 0);

      if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) != 0)
        {
          unsigned int fe_index = 0;
          unsigned int fe_index_old = 0;
          unsigned int i = 0;

          for (; i < fe.n_base_elements (); ++i)
            {
              fe_index_old = fe_index;
              fe_index += fe.element_multiplicity (i) * fe.base_element (i).n_components ();

              if (fe_index > first_vector_component)
                break;
            }

          base_indices.first = i;
          base_indices.second = (first_vector_component - fe_index_old) / fe.base_element (i).n_components ();
        }

      std::vector<Vector<double> >
      values (fe_values.n_quadrature_points, Vector<double> (fe.n_components ()));

      // Get boundary function
      // values at quadrature
      // points.
      boundary_function.vector_value_list (quadrature_points, values);

      switch (dim)
        {
        case 2:
        {
          const double tol = 0.5 * cell->face (face)->diameter () / cell->get_fe ().degree;
          std::vector<Point<dim> >
          tangentials (fe_values.n_quadrature_points);

          const std::vector<Point<dim> > &
          reference_quadrature_points = fe_values.get_quadrature ().get_points ();

          // coordinate directions
          // of the face.
          const unsigned int
          face_coordinate_direction[GeometryInfo<dim>::faces_per_cell]
            = { 1, 1, 0, 0 };
          const FEValuesExtractors::Vector vec (first_vector_component);

          // The interpolation for
          // the lowest order face
          // shape functions is just
          // the mean value of the
          // tangential  components
          // of the boundary function
          // on the edge.
          for (unsigned int q_point = 0;
               q_point < fe_values.n_quadrature_points; ++q_point)
            {
              // Therefore compute the
              // tangential of the
              // face at the quadrature
              // point.
              Point<dim> shifted_reference_point_1
                = reference_quadrature_points[q_point];
              Point<dim> shifted_reference_point_2
                = reference_quadrature_points[q_point];

              shifted_reference_point_1 (face_coordinate_direction[face])
              += tol;
              shifted_reference_point_2 (face_coordinate_direction[face])
              -= tol;
              tangentials[q_point]
                = (fe_values.get_mapping ()
                   .transform_unit_to_real_cell (cell,
                                                 shifted_reference_point_1)
                   -
                   fe_values.get_mapping ()
                   .transform_unit_to_real_cell (cell,
                                                 shifted_reference_point_2))
                  / tol;
              tangentials[q_point]
              /= std::sqrt (tangentials[q_point].square ());
              // Compute the degrees
              // of freedom.
              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices))
                    || (dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0))
                  {
                    dof_values[i]
                    += fe_values.JxW (q_point)
                       * (values[q_point] (first_vector_component)
                          * tangentials[q_point] (0)
                          + values[q_point] (first_vector_component + 1)
                          * tangentials[q_point] (1))
                       * (fe_values[vec].value (fe.face_to_cell_index (i, face), q_point)
                          * tangentials[q_point]);

                    if (q_point == 0)
                      dofs_processed[i] = true;
                  }
            }

          break;
        }

        case 3:
        {
          const FEValuesExtractors::Vector vec (first_vector_component);
          FullMatrix<double>
          assembling_matrix (degree * fe.degree,
                             dim * fe_values.n_quadrature_points);
          Vector<double> assembling_vector (assembling_matrix.n ());
          Vector<double> cell_rhs (assembling_matrix.m ());
          FullMatrix<double> cell_matrix (assembling_matrix.m (),
                                          assembling_matrix.m ());
          FullMatrix<double> cell_matrix_inv (assembling_matrix.m (),
                                              assembling_matrix.m ());
          Vector<double> solution (cell_matrix.m ());

          // Get coordinate directions
          // of the face.
          const unsigned int
          global_face_coordinate_directions[GeometryInfo<3>::faces_per_cell][2]
          = { { 1, 2 },
            { 1, 2 },
            { 2, 0 },
            { 2, 0 },
            { 0, 1 },
            { 0, 1 }
          };

          // The projection is divided into two steps.  In the first step we
          // project the boundary function on the horizontal shape functions.
          // Then the boundary function is projected on the vertical shape
          // functions.  We begin with the horizontal shape functions and set
          // up a linear system of equations to get the values for degrees of
          // freedom associated with the interior of the face.
          for (unsigned int q_point = 0;
               q_point < fe_values.n_quadrature_points; ++q_point)
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
                tmp[d] = values[q_point] (first_vector_component + d);

              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                     && (fe.base_element (base_indices.first).face_to_cell_index (2 * fe.degree, face)
                         <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                         <= fe.base_element (base_indices.first).face_to_cell_index (4 * fe.degree - 1, face)))
                    || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0) && (2 * fe.degree <= i)
                        && (i < 4 * fe.degree)))
                  tmp -= dof_values[i] * fe_values[vec].value (fe.face_to_cell_index (i, face), q_point);

              const double JxW
                = std::sqrt (fe_values.JxW (q_point)
                             / ((jacobians[q_point][0][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][0][global_face_coordinate_directions[face][0]]
                                 +
                                 jacobians[q_point][1][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][1][global_face_coordinate_directions[face][0]]
                                 +
                                 jacobians[q_point][2][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][2][global_face_coordinate_directions[face][0]])
                                *
                                (jacobians[q_point][0][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][0][global_face_coordinate_directions[face][1]]
                                 +
                                 jacobians[q_point][1][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][1][global_face_coordinate_directions[face][1]]
                                 +
                                 jacobians[q_point][2][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][2][global_face_coordinate_directions[face][1]])));

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
                assembling_vector (dim * q_point + d) = JxW * tmp[d];

              unsigned int index = 0;

              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                     && (fe.base_element (base_indices.first).face_to_cell_index
                         (GeometryInfo<dim>::lines_per_face * fe.degree, face)
                         <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                         < fe.base_element (base_indices.first).face_to_cell_index
                         ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree, face)))
                    || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                        && (GeometryInfo<dim>::lines_per_face * fe.degree <= i)
                        && (i < (degree + GeometryInfo<dim>::lines_per_face) * fe.degree)))
                  {
                    const Tensor<1, dim> shape_value
                      = (JxW * fe_values[vec].value (fe.face_to_cell_index (i, face),
                                                     q_point));

                    for (unsigned int d = 0; d < dim; ++d)
                      assembling_matrix (index, dim * q_point + d) = shape_value[d];

                    ++index;
                  }
            }

          // Create the system matrix by multiplying the assembling matrix
          // with its transposed and the right hand side vector by mutliplying
          // the assembling matrix with the assembling vector.  Invert the
          // system matrix.
          assembling_matrix.mTmult (cell_matrix, assembling_matrix);
          cell_matrix_inv.invert (cell_matrix);
          assembling_matrix.vmult (cell_rhs, assembling_vector);
          cell_matrix_inv.vmult (solution, cell_rhs);

          // Store the computed
          // values.
          {
            unsigned int index = 0;

            for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
              if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                   && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                   && (fe.base_element (base_indices.first).face_to_cell_index
                       (GeometryInfo<dim>::lines_per_face * fe.degree, face)
                       <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second)
                   && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                       < fe.base_element (base_indices.first).face_to_cell_index
                       ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree, face)))
                  || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                      && (GeometryInfo<dim>::lines_per_face * fe.degree <= i)
                      && (i < (degree + GeometryInfo<dim>::lines_per_face) * fe.degree)))
                {
                  dof_values[i] = solution (index);
                  dofs_processed[i] = true;
                  ++index;
                }
          }

          // Now we do the same as above with the vertical shape functions
          // instead of the horizontal ones.
          for (unsigned int q_point = 0;
               q_point < fe_values.n_quadrature_points; ++q_point)
            {
              Tensor<1, dim> tmp;

              for (unsigned int d = 0; d < dim; ++d)
                tmp[d] = values[q_point] (first_vector_component + d);

              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                         <= fe.base_element (base_indices.first).face_to_cell_index (2 * fe.degree - 1, face))
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).second
                         >= fe.base_element (base_indices.first).face_to_cell_index (0, face)))
                    || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0) && (i < 2 * fe.degree)))
                  tmp -= dof_values[i] * fe_values[vec].value (fe.face_to_cell_index (i, face), q_point);

              const double JxW
                = std::sqrt (fe_values.JxW (q_point)
                             / ((jacobians[q_point][0][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][0][global_face_coordinate_directions[face][0]]
                                 +
                                 jacobians[q_point][1][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][1][global_face_coordinate_directions[face][0]]
                                 +
                                 jacobians[q_point][2][global_face_coordinate_directions[face][0]]
                                 * jacobians[q_point][2][global_face_coordinate_directions[face][0]])
                                *
                                (jacobians[q_point][0][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][0][global_face_coordinate_directions[face][1]]
                                 +
                                 jacobians[q_point][1][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][1][global_face_coordinate_directions[face][1]]
                                 +
                                 jacobians[q_point][2][global_face_coordinate_directions[face][1]]
                                 * jacobians[q_point][2][global_face_coordinate_directions[face][1]])));

              for (unsigned int d = 0; d < dim; ++d)
                assembling_vector (dim * q_point + d) = JxW * tmp[d];

              unsigned int index = 0;

              for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
                if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                     && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                     && (fe.base_element (base_indices.first).face_to_cell_index
                         ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree, face)
                         <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second))
                    || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                        && ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree <= i)))
                  {
                    const Tensor<1, dim> shape_value
                      = JxW * fe_values[vec].value (fe.face_to_cell_index (i, face),
                                                    q_point);

                    for (unsigned int d = 0; d < dim; ++d)
                      assembling_matrix (index, dim * q_point + d) = shape_value[d];

                    ++index;
                  }
            }

          assembling_matrix.mTmult (cell_matrix, assembling_matrix);
          cell_matrix_inv.invert (cell_matrix);
          assembling_matrix.vmult (cell_rhs, assembling_vector);
          cell_matrix_inv.vmult (solution, cell_rhs);

          unsigned int index = 0;

          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                 && (fe.system_to_base_index (fe.face_to_cell_index (i, face)).first == base_indices)
                 && (fe.base_element (base_indices.first).face_to_cell_index
                     ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree, face)
                     <= fe.system_to_base_index (fe.face_to_cell_index (i, face)).second))
                || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                    && ((degree + GeometryInfo<dim>::lines_per_face) * fe.degree <= i)))
              {
                dof_values[i] = solution (index);
                dofs_processed[i] = true;
                ++index;
              }

          break;
        }

        default:
          Assert (false, ExcNotImplemented ());
        }
    }
  }




  template <int dim>
  void

  project_boundary_values_curl_conforming (const DoFHandler<dim> &dof_handler,
                                           const unsigned int first_vector_component,
                                           const Function<dim> &boundary_function,
                                           const types::boundary_id boundary_component,
                                           ConstraintMatrix &constraints,
                                           const Mapping<dim> &mapping)
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
    const unsigned int superdegree = dof_handler.get_fe ().degree;
    const QGauss<dim - 1> reference_face_quadrature (2 * superdegree);
    const unsigned int dofs_per_face = dof_handler.get_fe ().dofs_per_face;
    hp::FECollection<dim> fe_collection (dof_handler.get_fe ());
    hp::MappingCollection<dim> mapping_collection (mapping);
    hp::QCollection<dim> face_quadrature_collection;

    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      face_quadrature_collection.push_back
      (QProjector<dim>::project_to_face (reference_face_quadrature, face));

    hp::FEValues<dim> fe_face_values (mapping_collection, fe_collection,
                                      face_quadrature_collection,
                                      update_jacobians |
                                      update_JxW_values |
                                      update_quadrature_points |
                                      update_values);

    std::vector<bool> dofs_processed (dofs_per_face);
    std::vector<double> dof_values (dofs_per_face);
    std::vector<types::global_dof_index> face_dof_indices (dofs_per_face);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();

    switch (dim)
      {
      case 2:
      {
        for (; cell != dof_handler.end (); ++cell)
          if (cell->at_boundary () && cell->is_locally_owned ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // if the FE is a
                  // FE_Nothing object
                  // there is no work to
                  // do
                  if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                    return;

                  // This is only
                  // implemented, if the
                  // FE is a Nedelec
                  // element. If the FE
                  // is a FESystem, we
                  // cannot check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;
                      AssertThrow (dynamic_cast<const FE_Nedelec<dim>*> (&cell->get_fe ()) != 0,

                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    {
                      dof_values[dof] = 0.0;
                      dofs_processed[dof] = false;
                    }

                  // Compute the
                  // projection of the
                  // boundary function on
                  // the edge.
                  internals
                  ::compute_face_projection_curl_conforming (cell, face, fe_face_values,
                                                             boundary_function,
                                                             first_vector_component,
                                                             dof_values, dofs_processed);
                  cell->face (face)->get_dof_indices (face_dof_indices,
                                                      cell->active_fe_index ());

                  // Add the computed
                  // constraints to the
                  // constraint matrix,
                  // if the degree of
                  // freedom is not
                  // already constrained.
                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                        && !(constraints.is_constrained (face_dof_indices[dof])))
                      {
                        constraints.add_line (face_dof_indices[dof]);

                        if (std::abs (dof_values[dof]) > 1e-13)
                          constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
                      }
                }

        break;
      }

      case 3:
      {
        const QGauss<dim - 2> reference_edge_quadrature (2 * superdegree);
        const unsigned int degree = superdegree - 1;
        hp::QCollection<dim> edge_quadrature_collection;

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_face; ++line)
            edge_quadrature_collection.push_back
            (QProjector<dim>::project_to_face
             (QProjector<dim - 1>::project_to_face
              (reference_edge_quadrature, line), face));

        hp::FEValues<dim> fe_edge_values (mapping_collection, fe_collection,
                                          edge_quadrature_collection,
                                          update_jacobians |
                                          update_JxW_values |
                                          update_quadrature_points |
                                          update_values);

        for (; cell != dof_handler.end (); ++cell)
          if (cell->at_boundary () && cell->is_locally_owned ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // if the FE is a
                  // FE_Nothing object
                  // there is no work to
                  // do
                  if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                    return;

                  // This is only
                  // implemented, if the
                  // FE is a Nedelec
                  // element. If the FE is
                  // a FESystem we cannot
                  // check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_Nedelec<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    {
                      dof_values[dof] = 0.0;
                      dofs_processed[dof] = false;
                    }

                  // First we compute the
                  // projection on the
                  // edges.
                  for (unsigned int line = 0;
                       line < GeometryInfo<3>::lines_per_face; ++line)
                    internals
                    ::compute_edge_projection (cell, face, line, fe_edge_values,
                                               boundary_function,
                                               first_vector_component,
                                               dof_values, dofs_processed);

                  // If there are higher
                  // order shape
                  // functions, there is
                  // still some work
                  // left.
                  if (degree > 0)
                    internals
                    ::compute_face_projection_curl_conforming (cell, face, fe_face_values,
                                                               boundary_function,
                                                               first_vector_component,
                                                               dof_values,
                                                               dofs_processed);

                  // Store the computed
                  // values in the global
                  // vector.
                  cell->face (face)->get_dof_indices (face_dof_indices,
                                                      cell->active_fe_index ());

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                        && !(constraints.is_constrained (face_dof_indices[dof])))
                      {
                        constraints.add_line (face_dof_indices[dof]);

                        if (std::abs (dof_values[dof]) > 1e-13)
                          constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
                      }
                }

        break;
      }

      default:
        Assert (false, ExcNotImplemented ());
      }
  }



  template <int dim>
  void

  project_boundary_values_curl_conforming (const hp::DoFHandler<dim> &dof_handler,
                                           const unsigned int first_vector_component,
                                           const Function<dim> &boundary_function,
                                           const types::boundary_id boundary_component,
                                           ConstraintMatrix &constraints,
                                           const hp::MappingCollection<dim> &mapping_collection)
  {
    hp::FECollection<dim> fe_collection (dof_handler.get_fe ());
    hp::QCollection<dim> face_quadrature_collection;

    for (unsigned int i = 0; i < fe_collection.size (); ++i)
      {
        const QGauss<dim - 1>
        reference_face_quadrature (2 * fe_collection[i].degree);

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          face_quadrature_collection.push_back
          (QProjector<dim>::project_to_face (reference_face_quadrature, face));
      }

    hp::FEValues<dim> fe_face_values (mapping_collection, fe_collection,
                                      face_quadrature_collection,
                                      update_jacobians |
                                      update_JxW_values |
                                      update_quadrature_points |
                                      update_values);
    std::vector<bool> dofs_processed;
    std::vector<double> dof_values;
    std::vector<types::global_dof_index> face_dof_indices;
    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();

    switch (dim)
      {
      case 2:
      {
        for (; cell != dof_handler.end (); ++cell)
          if (cell->at_boundary () && cell->is_locally_owned ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // if the FE is a FE_Nothing object there is no work to do
                  if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                    return;

                  // This is only implemented, if the FE is a Nedelec
                  // element. If the FE is a FESystem we cannot check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_Nedelec<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  const unsigned int dofs_per_face = cell->get_fe ().dofs_per_face;

                  dofs_processed.resize (dofs_per_face);
                  dof_values.resize (dofs_per_face);

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    {
                      dof_values[dof] = 0.0;
                      dofs_processed[dof] = false;
                    }

                  internals
                  ::compute_face_projection_curl_conforming (cell, face, fe_face_values,
                                                             boundary_function,
                                                             first_vector_component,
                                                             dof_values, dofs_processed);
                  face_dof_indices.resize (dofs_per_face);
                  cell->face (face)->get_dof_indices (face_dof_indices,
                                                      cell->active_fe_index ());

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                        && !(constraints.is_constrained (face_dof_indices[dof])))
                      {
                        constraints.add_line (face_dof_indices[dof]);

                        if (std::abs (dof_values[dof]) > 1e-13)
                          constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
                      }
                }

        break;
      }

      case 3:
      {
        hp::QCollection<dim> edge_quadrature_collection;

        for (unsigned int i = 0; i < fe_collection.size (); ++i)
          {
            const QGauss<dim - 2>
            reference_edge_quadrature (2 * fe_collection[i].degree);

            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_face; ++line)
                edge_quadrature_collection.push_back
                (QProjector<dim>::project_to_face
                 (QProjector<dim - 1>::project_to_face (reference_edge_quadrature, line),
                  face));
          }

        hp::FEValues<dim> fe_edge_values (mapping_collection, fe_collection,
                                          edge_quadrature_collection,
                                          update_jacobians |
                                          update_JxW_values |
                                          update_quadrature_points |
                                          update_values);

        for (; cell != dof_handler.end (); ++cell)
          if (cell->at_boundary () && cell->is_locally_owned ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // if the FE is a FE_Nothing object there is no work to do
                  if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                    return;

                  // This is only implemented, if the FE is a Nedelec
                  // element. If the FE is a FESystem we cannot check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_Nedelec<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  const unsigned int superdegree = cell->get_fe ().degree;
                  const unsigned int degree = superdegree - 1;
                  const unsigned int dofs_per_face = cell->get_fe ().dofs_per_face;

                  dofs_processed.resize (dofs_per_face);
                  dof_values.resize (dofs_per_face);

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    {
                      dof_values[dof] = 0.0;
                      dofs_processed[dof] = false;
                    }

                  for (unsigned int line = 0;
                       line < GeometryInfo<dim>::lines_per_face; ++line)
                    internals
                    ::compute_edge_projection (cell, face, line, fe_edge_values,
                                               boundary_function,
                                               first_vector_component,
                                               dof_values, dofs_processed);

                  // If there are higher order shape functions, there is still
                  // some work left.
                  if (degree > 0)
                    internals
                    ::compute_face_projection_curl_conforming (cell, face, fe_face_values,
                                                               boundary_function,
                                                               first_vector_component,
                                                               dof_values, dofs_processed);


                  face_dof_indices.resize (dofs_per_face);
                  cell->face (face)->get_dof_indices (face_dof_indices,
                                                      cell->active_fe_index ());

                  for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                    if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                        && !(constraints.is_constrained (face_dof_indices[dof])))
                      {
                        constraints.add_line (face_dof_indices[dof]);

                        if (std::abs (dof_values[dof]) > 1e-13)
                          constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
                      }
                }

        break;
      }

      default:
        Assert (false, ExcNotImplemented ());
      }
  }


  namespace internals
  {
    // This function computes the projection of the boundary function on the
    // boundary in 2d.
    template <typename cell_iterator>
    void
    compute_face_projection_div_conforming (const cell_iterator &cell,
                                            const unsigned int face,
                                            const FEFaceValues<2> &fe_values,
                                            const unsigned int first_vector_component,
                                            const Function<2> &boundary_function,
                                            const std::vector<DerivativeForm<1,2,2> > &jacobians,
                                            ConstraintMatrix &constraints)
    {
      // Compute the intergral over the product of the normal components of
      // the boundary function times the normal components of the shape
      // functions supported on the boundary.
      const FEValuesExtractors::Vector vec (first_vector_component);
      const FiniteElement<2> &fe = cell->get_fe ();
      const std::vector<Point<2> > &normals = fe_values.get_normal_vectors ();
      const unsigned int
      face_coordinate_direction[GeometryInfo<2>::faces_per_cell] = {1, 1, 0, 0};
      std::vector<Vector<double> >
      values (fe_values.n_quadrature_points, Vector<double> (2));
      Vector<double> dof_values (fe.dofs_per_face);

      // Get the values of the boundary function at the quadrature points.
      {
        const std::vector<Point<2> > &
        quadrature_points = fe_values.get_quadrature_points ();

        boundary_function.vector_value_list (quadrature_points, values);
      }

      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; ++q_point)
        {
          double tmp = 0.0;

          for (unsigned int d = 0; d < 2; ++d)
            tmp += normals[q_point][d] * values[q_point] (d);

          tmp *= fe_values.JxW (q_point)
                 * std::sqrt (jacobians[q_point][0][face_coordinate_direction[face]]
                              * jacobians[q_point][0][face_coordinate_direction[face]]
                              + jacobians[q_point][1][face_coordinate_direction[face]]
                              * jacobians[q_point][1][face_coordinate_direction[face]]);

          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            dof_values (i) += tmp * (normals[q_point]
                                     * fe_values[vec].value (fe.face_to_cell_index (i, face), q_point));
        }

      std::vector<types::global_dof_index> face_dof_indices (fe.dofs_per_face);

      cell->face (face)->get_dof_indices (face_dof_indices, cell->active_fe_index ());

      // Copy the computed values in the ConstraintMatrix only, if the degree
      // of freedom is not already constrained.
      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        if (!(constraints.is_constrained (face_dof_indices[i])))
          {
            constraints.add_line (face_dof_indices[i]);

            if (std::abs (dof_values (i)) > 1e-14)
              constraints.set_inhomogeneity (face_dof_indices[i], dof_values (i));
          }
    }

    // dummy implementation of above function for all other dimensions
    template<int dim, typename cell_iterator>
    void
    compute_face_projection_div_conforming (const cell_iterator &,
                                            const unsigned int,
                                            const FEFaceValues<dim> &,
                                            const unsigned int,
                                            const Function<dim> &,
                                            const std::vector<DerivativeForm<1,dim,dim> > &,
                                            ConstraintMatrix &)
    {
      Assert (false, ExcNotImplemented ());
    }

    // This function computes the projection of the boundary function on the
    // boundary in 3d.
    template<typename cell_iterator>
    void
    compute_face_projection_div_conforming (const cell_iterator &cell,
                                            const unsigned int face,
                                            const FEFaceValues<3> &fe_values,
                                            const unsigned int first_vector_component,
                                            const Function<3> &boundary_function,
                                            const std::vector<DerivativeForm<1,3,3> > &jacobians,
                                            std::vector<double> &dof_values,
                                            std::vector<types::global_dof_index> &projected_dofs)
    {
      // Compute the intergral over the product of the normal components of
      // the boundary function times the normal components of the shape
      // functions supported on the boundary.
      const FEValuesExtractors::Vector vec (first_vector_component);
      const FiniteElement<3> &fe = cell->get_fe ();
      const std::vector<Point<3> > &normals = fe_values.get_normal_vectors ();
      const unsigned int
      face_coordinate_directions[GeometryInfo<3>::faces_per_cell][2] = {{1, 2},
        {1, 2},
        {2, 0},
        {2, 0},
        {0, 1},
        {0, 1}
      };
      std::vector<Vector<double> >
      values (fe_values.n_quadrature_points, Vector<double> (3));
      Vector<double> dof_values_local (fe.dofs_per_face);

      {
        const std::vector<Point<3> > &
        quadrature_points = fe_values.get_quadrature_points ();

        boundary_function.vector_value_list (quadrature_points, values);
      }

      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; ++q_point)
        {
          double tmp = 0.0;

          for (unsigned int d = 0; d < 3; ++d)
            tmp += normals[q_point][d] * values[q_point] (d);

          tmp *= fe_values.JxW (q_point)
                 * std::sqrt ((jacobians[q_point][0][face_coordinate_directions[face][0]]
                               * jacobians[q_point][0][face_coordinate_directions[face][0]]
                               + jacobians[q_point][1][face_coordinate_directions[face][0]]
                               * jacobians[q_point][1][face_coordinate_directions[face][0]]
                               + jacobians[q_point][2][face_coordinate_directions[face][0]]
                               * jacobians[q_point][2][face_coordinate_directions[face][0]])
                              * (jacobians[q_point][0][face_coordinate_directions[face][1]]
                                 * jacobians[q_point][0][face_coordinate_directions[face][1]]
                                 + jacobians[q_point][1][face_coordinate_directions[face][1]]
                                 * jacobians[q_point][1][face_coordinate_directions[face][1]]
                                 + jacobians[q_point][2][face_coordinate_directions[face][1]]
                                 * jacobians[q_point][2][face_coordinate_directions[face][1]]));

          for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
            dof_values_local (i) += tmp * (normals[q_point]
                                           * fe_values[vec].value (fe.face_to_cell_index (i, face), q_point));
        }

      std::vector<types::global_dof_index> face_dof_indices (fe.dofs_per_face);

      cell->face (face)->get_dof_indices (face_dof_indices, cell->active_fe_index ());

      for (unsigned int i = 0; i < fe.dofs_per_face; ++i)
        if (projected_dofs[face_dof_indices[i]] < fe.degree)
          {
            dof_values[face_dof_indices[i]] = dof_values_local (i);
            projected_dofs[face_dof_indices[i]] = fe.degree;
          }
    }

    // dummy implementation of above
    // function for all other
    // dimensions
    template<int dim, typename cell_iterator>
    void
    compute_face_projection_div_conforming (const cell_iterator &,
                                            const unsigned int,
                                            const FEFaceValues<dim> &,
                                            const unsigned int,
                                            const Function<dim> &,
                                            const std::vector<DerivativeForm<1,dim,dim> > &,
                                            std::vector<double> &,
                                            std::vector<types::global_dof_index> &)
    {
      Assert (false, ExcNotImplemented ());
    }
  }


  template <int dim>
  void
  project_boundary_values_div_conforming (const DoFHandler<dim> &dof_handler,
                                          const unsigned int first_vector_component,
                                          const Function<dim> &boundary_function,
                                          const types::boundary_id boundary_component,
                                          ConstraintMatrix &constraints,
                                          const Mapping<dim> &mapping)
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
    const FiniteElement<dim> &fe = dof_handler.get_fe ();
    QGauss<dim - 1> face_quadrature (fe.degree + 1);
    FEFaceValues<dim> fe_face_values (mapping, fe, face_quadrature, update_JxW_values |
                                      update_normal_vectors |
                                      update_quadrature_points |
                                      update_values);
    hp::FECollection<dim> fe_collection (fe);
    hp::MappingCollection<dim> mapping_collection (mapping);
    hp::QCollection<dim> quadrature_collection;

    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      quadrature_collection.push_back (QProjector<dim>::project_to_face (face_quadrature,
                                       face));

    hp::FEValues<dim> fe_values (mapping_collection, fe_collection, quadrature_collection,
                                 update_jacobians);

    switch (dim)
      {
      case 2:
      {
        for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();
             cell != dof_handler.end (); ++cell)
          if (cell->at_boundary ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // if the FE is a
                  // FE_Nothing object
                  // there is no work to
                  // do
                  if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                    return;

                  // This is only
                  // implemented, if the
                  // FE is a Raviart-Thomas
                  // element. If the FE is
                  // a FESystem we cannot
                  // check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_RaviartThomas<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  fe_values.reinit (cell, face + cell->active_fe_index ()
                                    * GeometryInfo<dim>::faces_per_cell);

                  const std::vector<DerivativeForm<1,dim,spacedim> > &
                  jacobians = fe_values.get_present_fe_values ().get_jacobians ();

                  fe_face_values.reinit (cell, face);
                  internals::compute_face_projection_div_conforming (cell, face,
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
        // In three dimensions the
        // edges between two faces
        // are treated twice.
        // Therefore we store the
        // computed values in a
        // vector and copy them over
        // in the ConstraintMatrix
        // after all values have been
        // computed.
        // If we have two values for
        // one edge, we choose the one,
        // which was computed with the
        // higher order element.
        // If both elements are of the
        // same order, we just keep the
        // first value and do not
        // compute a second one.
        const unsigned int &n_dofs = dof_handler.n_dofs ();
        std::vector<double> dof_values (n_dofs);
        std::vector<types::global_dof_index> projected_dofs (n_dofs);

        for (unsigned int dof = 0; dof < n_dofs; ++dof)
          projected_dofs[dof] = 0;

        for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();
             cell != dof_handler.end (); ++cell)
          if (cell->at_boundary ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // This is only
                  // implemented, if the
                  // FE is a Raviart-Thomas
                  // element. If the FE is
                  // a FESystem we cannot
                  // check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_RaviartThomas<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  fe_values.reinit (cell, face + cell->active_fe_index ()
                                    * GeometryInfo<dim>::faces_per_cell);

                  const std::vector<DerivativeForm<1,dim ,spacedim> > &
                  jacobians = fe_values.get_present_fe_values ().get_jacobians ();

                  fe_face_values.reinit (cell, face);
                  internals::compute_face_projection_div_conforming (cell, face,
                                                                     fe_face_values,
                                                                     first_vector_component,
                                                                     boundary_function,
                                                                     jacobians, dof_values,
                                                                     projected_dofs);
                }

        for (unsigned int dof = 0; dof < n_dofs; ++dof)
          if ((projected_dofs[dof] != 0) && !(constraints.is_constrained (dof)))
            {
              constraints.add_line (dof);

              if (std::abs (dof_values[dof]) > 1e-14)
                constraints.set_inhomogeneity (dof, dof_values[dof]);
            }

        break;
      }

      default:
        Assert (false, ExcNotImplemented ());
      }
  }


  template <int dim>
  void
  project_boundary_values_div_conforming (const hp::DoFHandler<dim> &dof_handler,
                                          const unsigned int first_vector_component,
                                          const Function<dim> &boundary_function,
                                          const types::boundary_id boundary_component,
                                          ConstraintMatrix &constraints,
                                          const hp::MappingCollection<dim, dim> &mapping_collection)
  {
    const unsigned int spacedim = dim;
    const hp::FECollection<dim> &fe_collection = dof_handler.get_fe ();
    hp::QCollection<dim - 1> face_quadrature_collection;
    hp::QCollection<dim> quadrature_collection;

    for (unsigned int i = 0; i < fe_collection.size (); ++i)
      {
        const QGauss<dim - 1> quadrature (fe_collection[i].degree + 1);

        face_quadrature_collection.push_back (quadrature);

        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          quadrature_collection.push_back (QProjector<dim>::project_to_face (quadrature,
                                           face));
      }

    hp::FEFaceValues<dim> fe_face_values (mapping_collection, fe_collection,
                                          face_quadrature_collection, update_JxW_values |
                                          update_normal_vectors |
                                          update_quadrature_points |
                                          update_values);
    hp::FEValues<dim> fe_values (mapping_collection, fe_collection, quadrature_collection,
                                 update_jacobians);

    switch (dim)
      {
      case 2:
      {
        for (typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();
             cell != dof_handler.end (); ++cell)
          if (cell->at_boundary ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // This is only
                  // implemented, if the
                  // FE is a Raviart-Thomas
                  // element. If the FE is
                  // a FESystem we cannot
                  // check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_RaviartThomas<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  fe_values.reinit (cell, face + cell->active_fe_index ()
                                    * GeometryInfo<dim>::faces_per_cell);

                  const std::vector<DerivativeForm<1,dim,spacedim> > &
                  jacobians = fe_values.get_present_fe_values ().get_jacobians ();

                  fe_face_values.reinit (cell, face);
                  internals::compute_face_projection_div_conforming (cell, face,
                                                                     fe_face_values.get_present_fe_values (),
                                                                     first_vector_component,
                                                                     boundary_function,
                                                                     jacobians,
                                                                     constraints);
                }

        break;
      }

      case 3:
      {
        const unsigned int &n_dofs = dof_handler.n_dofs ();
        std::vector<double> dof_values (n_dofs);
        std::vector<types::global_dof_index> projected_dofs (n_dofs);

        for (unsigned int dof = 0; dof < n_dofs; ++dof)
          projected_dofs[dof] = 0;

        for (typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active ();
             cell != dof_handler.end (); ++cell)
          if (cell->at_boundary ())
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face (face)->boundary_indicator () == boundary_component)
                {
                  // This is only
                  // implemented, if the
                  // FE is a Raviart-Thomas
                  // element. If the FE is
                  // a FESystem we cannot
                  // check this.
                  if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) == 0)
                    {
                      typedef FiniteElement<dim> FEL;

                      AssertThrow (dynamic_cast<const FE_RaviartThomas<dim>*> (&cell->get_fe ()) != 0,
                                   typename FEL::ExcInterpolationNotImplemented ());
                    }

                  fe_values.reinit (cell, face + cell->active_fe_index ()
                                    * GeometryInfo<dim>::faces_per_cell);

                  const std::vector<DerivativeForm<1,dim,spacedim> > &
                  jacobians = fe_values.get_present_fe_values ().get_jacobians ();

                  fe_face_values.reinit (cell, face);
                  internals::compute_face_projection_div_conforming (cell, face,
                                                                     fe_face_values.get_present_fe_values (),
                                                                     first_vector_component,
                                                                     boundary_function,
                                                                     jacobians, dof_values,
                                                                     projected_dofs);
                }

        for (unsigned int dof = 0; dof < n_dofs; ++dof)
          if ((projected_dofs[dof] != 0) && !(constraints.is_constrained (dof)))
            {
              constraints.add_line (dof);

              if (std::abs (dof_values[dof]) > 1e-14)
                constraints.set_inhomogeneity (dof, dof_values[dof]);
            }

        break;
      }

      default:
        Assert (false, ExcNotImplemented ());
      }
  }



  template <int dim, template <int, int> class DH, int spacedim>
  void
  compute_no_normal_flux_constraints (const DH<dim,spacedim>         &dof_handler,
                                      const unsigned int     first_vector_component,
                                      const std::set<types::boundary_id> &boundary_ids,
                                      ConstraintMatrix      &constraints,
                                      const Mapping<dim, spacedim>    &mapping)
  {
    ZeroFunction<dim>zero_function(dim);
    typename FunctionMap<spacedim>::type function_map;
    std::set<types::boundary_id>::const_iterator it
      = boundary_ids.begin();
    for (; it != boundary_ids.end(); ++it)
      function_map[*it] = &zero_function;
    compute_nonzero_normal_flux_constraints(dof_handler,
                                            first_vector_component,
                                            boundary_ids,
                                            function_map,
                                            constraints,
                                            mapping);
  }

  template <int dim, template <int, int> class DH, int spacedim>
  void
  compute_nonzero_normal_flux_constraints (const DH<dim,spacedim>         &dof_handler,
                                           const unsigned int     first_vector_component,
                                           const std::set<types::boundary_id> &boundary_ids,
                                           typename FunctionMap<spacedim>::type &function_map,
                                           ConstraintMatrix      &constraints,
                                           const Mapping<dim, spacedim>    &mapping)
  {
    Assert (dim > 1,
            ExcMessage ("This function is not useful in 1d because it amounts "
                        "to imposing Dirichlet values on the vector-valued "
                        "quantity."));

    std::vector<types::global_dof_index> face_dofs;

    // create FE and mapping collections for all elements in use by this
    // DoFHandler
    hp::FECollection<dim,spacedim>      fe_collection (dof_handler.get_fe());
    hp::MappingCollection<dim,spacedim> mapping_collection;
    for (unsigned int i=0; i<fe_collection.size(); ++i)
      mapping_collection.push_back (mapping);

    // now also create a quadrature collection for the faces of a cell. fill
    // it with a quadrature formula with the support points on faces for each
    // FE
    hp::QCollection<dim-1> face_quadrature_collection;
    for (unsigned int i=0; i<fe_collection.size(); ++i)
      {
        const std::vector<Point<dim-1> > &
        unit_support_points = fe_collection[i].get_unit_face_support_points();

        Assert (unit_support_points.size() == fe_collection[i].dofs_per_face,
                ExcInternalError());

        face_quadrature_collection.push_back (Quadrature<dim-1> (unit_support_points));
      }

    // now create the object with which we will generate the normal vectors
    hp::FEFaceValues<dim,spacedim> x_fe_face_values (mapping_collection,
                                                     fe_collection,
                                                     face_quadrature_collection,
                                                     update_q_points |
                                                     update_normal_vectors);

    // have a map that stores normal vectors for each vector-dof tuple we want
    // to constrain. since we can get at the same vector dof tuple more than
    // once (for example if it is located at a vertex that we visit from all
    // adjacent cells), we will want to average later on the normal vectors
    // computed on different cells as described in the documentation of this
    // function. however, we can only average if the contributions came from
    // different cells, whereas we want to constrain twice or more in case the
    // contributions came from different faces of the same cell
    // (i.e. constrain not just the *average normal direction* but *all normal
    // directions* we find). consequently, we also have to store which cell a
    // normal vector was computed on
    typedef
    std::multimap<internal::VectorDoFTuple<dim>,
        std::pair<Tensor<1,dim>, typename DH<dim,spacedim>::active_cell_iterator> >
        DoFToNormalsMap;
    std::map<internal::VectorDoFTuple<dim>, Vector<double> >
    dof_vector_to_b_values;

    DoFToNormalsMap dof_to_normals_map;

    // now loop over all cells and all faces
    typename DH<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    std::set<types::boundary_id>::iterator b_id;
    for (; cell!=endc; ++cell)
      if (!cell->is_artificial())
        for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if ((b_id=boundary_ids.find(cell->face(face_no)->boundary_indicator()))
              != boundary_ids.end())
            {
              const FiniteElement<dim> &fe = cell->get_fe ();
              typename DH<dim,spacedim>::face_iterator face = cell->face(face_no);

              // get the indices of the dofs on this cell...
              face_dofs.resize (fe.dofs_per_face);
              face->get_dof_indices (face_dofs, cell->active_fe_index());

              x_fe_face_values.reinit (cell, face_no);
              const FEFaceValues<dim> &fe_values = x_fe_face_values.get_present_fe_values();

              // then identify which of them correspond to the selected set of
              // vector components
              for (unsigned int i=0; i<face_dofs.size(); ++i)
                if (fe.face_system_to_component_index(i).first ==
                    first_vector_component)
                  {
                    // find corresponding other components of vector
                    internal::VectorDoFTuple<dim> vector_dofs;
                    vector_dofs.dof_indices[0] = face_dofs[i];

                    Assert(first_vector_component+dim<=fe.n_components(),
                           ExcMessage("Error: the finite element does not have enough components "
                                      "to define a normal direction."));

                    for (unsigned int k=0; k<fe.dofs_per_face; ++k)
                      if ((k != i)
                          &&
                          (face_quadrature_collection[cell->active_fe_index()].point(k) ==
                           face_quadrature_collection[cell->active_fe_index()].point(i))
                          &&
                          (fe.face_system_to_component_index(k).first >=
                           first_vector_component)
                          &&
                          (fe.face_system_to_component_index(k).first <
                           first_vector_component + dim))
                        vector_dofs.dof_indices[fe.face_system_to_component_index(k).first -
                                                first_vector_component]
                          = face_dofs[k];

                    for (unsigned int d=0; d<dim; ++d)
                      Assert (vector_dofs.dof_indices[d] < dof_handler.n_dofs(),
                              ExcInternalError());

                    // we need the normal vector on this face. we know that it
                    // is a vector of length 1 but at least with higher order
                    // mappings it isn't always possible to guarantee that
                    // each component is exact up to zero tolerance. in
                    // particular, as shown in the deal.II/no_flux_06 test, if
                    // we just take the normal vector as given by the
                    // fe_values object, we can get entries in the normal
                    // vectors of the unit cube that have entries up to
                    // several times 1e-14.
                    //
                    // the problem with this is that this later yields
                    // constraints that are circular (e.g., in the testcase,
                    // we get constraints of the form
                    //
                    // x22 =  2.93099e-14*x21 + 2.93099e-14*x23
                    // x21 = -2.93099e-14*x22 + 2.93099e-14*x21
                    //
                    // in both of these constraints, the small numbers should
                    // be zero and the constraints should simply be
                    // x22 = x21 = 0
                    //
                    // to achieve this, we utilize that we know that the
                    // normal vector has (or should have) length 1 and that we
                    // can simply set small elements to zero (without having
                    // to check that they are small *relative to something
                    // else*). we do this and then normalize the length of the
                    // vector back to one, just to be on the safe side
                    //
                    // one more point: we would like to use the "real" normal
                    // vector here, as provided by the boundary description
                    // and as opposed to what we get from the FEValues object.
                    // we do this in the immediately next line, but as is
                    // obvious, the boundary only has a vague idea which side
                    // of a cell it is on -- indicated by the face number. in
                    // other words, it may provide the inner or outer normal.
                    // by and large, there is no harm from this, since the
                    // tangential vector we compute is still the same. however,
                    // we do average over normal vectors from adjacent cells
                    // and if they have recorded normal vectors from the inside
                    // once and from the outside the other time, then this
                    // averaging is going to run into trouble. as a consequence
                    // we ask the mapping after all for its normal vector,
                    // but we only ask it so that we can possibly correct the
                    // sign of the normal vector provided by the boundary
                    // if they should point in different directions. this is the
                    // case in tests/deal.II/no_flux_11.
                    Point<dim> normal_vector
                      = (cell->face(face_no)->get_boundary().normal_vector
                         (cell->face(face_no), fe_values.quadrature_point(i)));
                    if (normal_vector * fe_values.normal_vector(i) < 0)
                      normal_vector *= -1;
                    Assert (std::fabs(normal_vector.norm() - 1) < 1e-14,
                            ExcInternalError());
                    for (unsigned int d=0; d<dim; ++d)
                      if (std::fabs(normal_vector[d]) < 1e-13)
                        normal_vector[d] = 0;
                    normal_vector /= normal_vector.norm();

                    const Point<dim> point
                      = fe_values.quadrature_point(i);
                    Vector<double> b_values(dim);
                    function_map[*b_id]->vector_value(point, b_values);

                    // now enter the (dofs,(normal_vector,cell)) entry into
                    // the map
                    dof_to_normals_map.insert
                    (std::make_pair (vector_dofs,
                                     std::make_pair (normal_vector,cell)));
                    dof_vector_to_b_values.insert
                    (std::make_pair(vector_dofs, b_values));

#ifdef DEBUG_NO_NORMAL_FLUX
                    std::cout << "Adding normal vector:" << std::endl
                              << "   dofs=" << vector_dofs << std::endl
                              << "   cell=" << cell << " at " << cell->center() << std::endl
                              << "   normal=" << normal_vector << std::endl;
#endif
                  }
            }

    // Now do something with the collected information. To this end, loop
    // through all sets of pairs (dofs,normal_vector) and identify which
    // entries belong to the same set of dofs and then do as described in the
    // documentation, i.e. either average the normal vector or don't for this
    // particular set of dofs
    typename DoFToNormalsMap::const_iterator
    p = dof_to_normals_map.begin();

    while (p != dof_to_normals_map.end())
      {
        // first find the range of entries in the multimap that corresponds to
        // the same vector-dof tuple. as usual, we define the range
        // half-open. the first entry of course is 'p'
        typename DoFToNormalsMap::const_iterator same_dof_range[2] = { p };
        for (++p; p != dof_to_normals_map.end(); ++p)
          if (p->first != same_dof_range[0]->first)
            {
              same_dof_range[1] = p;
              break;
            }
        if (p == dof_to_normals_map.end())
          same_dof_range[1] = dof_to_normals_map.end();

#ifdef DEBUG_NO_NORMAL_FLUX
        std::cout << "For dof indices <" << p->first << ">, found the following normals"
                  << std::endl;
        for (typename DoFToNormalsMap::const_iterator
             q = same_dof_range[0];
             q != same_dof_range[1]; ++q)
          std::cout << "   " << q->second.first
                    << " from cell " << q->second.second
                    << std::endl;
#endif


        // now compute the reverse mapping: for each of the cells that
        // contributed to the current set of vector dofs, add up the normal
        // vectors. the values of the map are pairs of normal vectors and
        // number of cells that have contributed
        typedef std::map<typename DH<dim,spacedim>::active_cell_iterator,
                std::pair<Tensor<1,dim>, unsigned int> >
                CellToNormalsMap;

        CellToNormalsMap cell_to_normals_map;
        for (typename DoFToNormalsMap::const_iterator
             q = same_dof_range[0];
             q != same_dof_range[1]; ++q)
          if (cell_to_normals_map.find (q->second.second)
              == cell_to_normals_map.end())
            cell_to_normals_map[q->second.second]
              = std::make_pair (q->second.first, 1U);
          else
            {
              const Tensor<1,dim> old_normal
                = cell_to_normals_map[q->second.second].first;
              const unsigned int old_count
                = cell_to_normals_map[q->second.second].second;

              Assert (old_count > 0, ExcInternalError());

              // in the same entry, store again the now averaged normal vector
              // and the new count
              cell_to_normals_map[q->second.second]
                = std::make_pair ((old_normal * old_count + q->second.first) / (old_count + 1),
                                  old_count + 1);
            }
        Assert (cell_to_normals_map.size() >= 1, ExcInternalError());

#ifdef DEBUG_NO_NORMAL_FLUX
        std::cout << "   cell_to_normals_map:" << std::endl;
        for (typename CellToNormalsMap::const_iterator
             x = cell_to_normals_map.begin();
             x != cell_to_normals_map.end(); ++x)
          std::cout << "      " << x->first << " -> ("
                    << x->second.first << ',' << x->second.second << ')'
                    << std::endl;
#endif

        // count the maximum number of contributions from each cell
        unsigned int max_n_contributions_per_cell = 1;
        for (typename CellToNormalsMap::const_iterator
             x = cell_to_normals_map.begin();
             x != cell_to_normals_map.end(); ++x)
          max_n_contributions_per_cell
            = std::max (max_n_contributions_per_cell,
                        x->second.second);

        // verify that each cell can have only contributed at most dim times,
        // since that is the maximum number of faces that come together at a
        // single place
        Assert (max_n_contributions_per_cell <= dim, ExcInternalError());

        switch (max_n_contributions_per_cell)
          {
          // first deal with the case that a number of cells all have
          // registered that they have a normal vector defined at the
          // location of a given vector dof, and that each of them have
          // encountered this vector dof exactly once while looping over all
          // their faces. as stated in the documentation, this is the case
          // where we want to simply average over all normal vectors
          //
          // the typical case is in 2d where multiple cells meet at one
          // vertex sitting on the boundary. same in 3d for a vertex that
          // is associated with only one of the boundary indicators passed
          // to this function
          case 1:
          {
            // compute the average normal vector from all the ones that have
            // the same set of dofs. we could add them up and divide them by
            // the number of additions, or simply normalize them right away
            // since we want them to have unit length anyway
            Tensor<1,dim> normal;
            for (typename CellToNormalsMap::const_iterator
                 x = cell_to_normals_map.begin();
                 x != cell_to_normals_map.end(); ++x)
              normal += x->second.first;
            normal /= normal.norm();

            // normalize again
            for (unsigned int d=0; d<dim; ++d)
              if (std::fabs(normal[d]) < 1e-13)
                normal[d] = 0;
            normal /= normal.norm();

            // then construct constraints from this:
            const internal::VectorDoFTuple<dim> &
            dof_indices = same_dof_range[0]->first;
            double normal_value = 0.;
            const Vector<double> b_values = dof_vector_to_b_values[dof_indices];
            for (unsigned int i=0; i<dim; ++i)
              normal_value += b_values[i]*normal[i];
            internal::add_constraint (dof_indices, normal,
                                      constraints, normal_value);

            break;
          }

          // this is the slightly more complicated case that a single cell has
          // contributed with exactly DIM normal vectors to the same set of
          // vector dofs. this is what happens in a corner in 2d and 3d (but
          // not on an edge in 3d, where we have only 2, i.e. <DIM,
          // contributions. Here we do not want to average the normal
          // vectors. Since we have DIM contributions, let's assume (and
          // verify) that they are in fact all linearly independent; in that
          // case, all vector components are constrained and we need to set all
          // of them to the corresponding boundary values
          case dim:
          {
            // assert that indeed only a single cell has contributed
            Assert (cell_to_normals_map.size() == 1,
                    ExcInternalError());

            // check linear independence by computing the determinant of the
            // matrix created from all the normal vectors. if they are
            // linearly independent, then the determinant is nonzero. if they
            // are orthogonal, then the matrix is in fact equal to 1 (since
            // they are all unit vectors); make sure the determinant is larger
            // than 1e-3 to avoid cases where cells are degenerate
            {
              Tensor<2,dim> t;

              typename DoFToNormalsMap::const_iterator x = same_dof_range[0];
              for (unsigned int i=0; i<dim; ++i, ++x)
                for (unsigned int j=0; j<dim; ++j)
                  t[i][j] = x->second.first[j];

              Assert (std::fabs(determinant (t)) > 1e-3,
                      ExcMessage("Found a set of normal vectors that are nearly collinear."));
            }

            // so all components of this vector dof are constrained. enter
            // this into the constraint matrix
            //
            // ignore dofs already constrained
            const internal::VectorDoFTuple<dim> &
            dof_indices = same_dof_range[0]->first;
            const Vector<double> b_values = dof_vector_to_b_values[dof_indices];
            for (unsigned int i=0; i<dim; ++i)
              if (!constraints.is_constrained(same_dof_range[0]->first.dof_indices[i])
                  &&
                  constraints.can_store_line(same_dof_range[0]->first.dof_indices[i]))
                {
                  const types::global_dof_index line
                    = dof_indices.dof_indices[i];
                  constraints.add_line (line);
                  if (std::fabs(b_values[i])
                      > std::numeric_limits<double>::epsilon())
                    constraints.set_inhomogeneity(line, b_values[i]);
                  // no add_entries here
                }

            break;
          }

          // this is the case of an edge contribution in 3d, i.e. the vector
          // is constrained in two directions but not the third.
          default:
          {
            Assert (dim >= 3, ExcNotImplemented());
            Assert (max_n_contributions_per_cell == 2, ExcInternalError());

            // as described in the documentation, let us first collect what
            // each of the cells contributed at the current point. we use a
            // std::list instead of a std::set (which would be more natural)
            // because std::set requires that the stored elements are
            // comparable with operator<
            typedef std::map<typename DH<dim,spacedim>::active_cell_iterator,
                    std::list<Tensor<1,dim> > >
                    CellContributions;
            CellContributions cell_contributions;

            for (typename DoFToNormalsMap::const_iterator
                 q = same_dof_range[0];
                 q != same_dof_range[1]; ++q)
              cell_contributions[q->second.second].push_back (q->second.first);
            Assert (cell_contributions.size() >= 1, ExcInternalError());

            // now for each cell that has contributed determine the number of
            // normal vectors it has contributed. we currently only implement
            // if this is dim-1 for all cells (if a single cell has
            // contributed dim, or if all adjacent cells have contributed 1
            // normal vector, this is already handled above).
            //
            // we only implement the case that all cells contribute
            // dim-1 because we assume that we are following an edge
            // of the domain (think: we are looking at a vertex
            // located on one of the edges of a refined cube where the
            // boundary indicators of the two adjacent faces of the
            // cube are both listed in the set of boundary indicators
            // passed to this function). in that case, all cells along
            // that edge of the domain are assumed to have contributed
            // dim-1 normal vectors. however, there are cases where
            // this assumption is not justified (see the lengthy
            // explanation in test no_flux_12.cc) and in those cases
            // we simply ignore the cell that contributes only
            // once. this is also discussed at length in the
            // documentation of this function.
            //
            // for each contributing cell compute the tangential vector that
            // remains unconstrained
            std::list<Tensor<1,dim> > tangential_vectors;
            for (typename CellContributions::const_iterator
                 contribution = cell_contributions.begin();
                 contribution != cell_contributions.end();
                 ++contribution)
              {
#ifdef DEBUG_NO_NORMAL_FLUX
                std::cout << "   Treating edge case with dim-1 contributions." << std::endl
                          << "   Looking at cell " << contribution->first
                          << " which has contributed these normal vectors:"
                          << std::endl;
                for (typename std::list<Tensor<1,dim> >::const_iterator
                     t = contribution->second.begin();
                     t != contribution->second.end();
                     ++t)
                  std::cout << "      " << *t << std::endl;
#endif

                // as mentioned above, simply ignore cells that only
                // contribute once
                if (contribution->second.size() < dim-1)
                  continue;

                Tensor<1,dim> normals[dim-1];
                {
                  unsigned int index=0;
                  for (typename std::list<Tensor<1,dim> >::const_iterator
                       t = contribution->second.begin();
                       t != contribution->second.end();
                       ++t, ++index)
                    normals[index] = *t;
                  Assert (index == dim-1, ExcInternalError());
                }

                // calculate the tangent as the outer product of the normal
                // vectors. since these vectors do not need to be orthogonal
                // (think, for example, the case of the deal.II/no_flux_07
                // test: a sheared cube in 3d, with Q2 elements, where we have
                // constraints from the two normal vectors of two faces of the
                // sheared cube that are not perpendicular to each other), we
                // have to normalize the outer product
                Tensor<1,dim> tangent;
                switch (dim)
                  {
                  case 3:
                    // take cross product between normals[0] and
                    // normals[1]. write it in the current form (with [dim-2])
                    // to make sure that compilers don't warn about
                    // out-of-bounds accesses -- the warnings are bogus since
                    // we get here only for dim==3, but at least one isn't
                    // quite smart enough to notice this and warns when
                    // compiling the function in 2d
                    cross_product (tangent, normals[0], normals[dim-2]);
                    break;
                  default:
                    Assert (false, ExcNotImplemented());
                  }

                Assert (std::fabs (tangent.norm()) > 1e-12,
                        ExcMessage("Two normal vectors from adjacent faces are almost "
                                   "parallel."));
                tangent /= tangent.norm();

                tangential_vectors.push_back (tangent);
              }

            // go through the list of tangents and make sure that they all
            // roughly point in the same direction as the first one (i.e. have
            // an angle less than 90 degrees); if they don't then flip their
            // sign
            {
              const Tensor<1,dim> first_tangent = tangential_vectors.front();
              typename std::list<Tensor<1,dim> >::iterator
              t = tangential_vectors.begin();
              ++t;
              for (; t != tangential_vectors.end(); ++t)
                if (*t * first_tangent < 0)
                  *t *= -1;
            }

            // now compute the average tangent and normalize it
            Tensor<1,dim> average_tangent;
            for (typename std::list<Tensor<1,dim> >::const_iterator
                 t = tangential_vectors.begin();
                 t != tangential_vectors.end();
                 ++t)
              average_tangent += *t;
            average_tangent /= average_tangent.norm();

            // now all that is left is that we add the constraints that the
            // vector is parallel to the tangent
            const internal::VectorDoFTuple<dim> &
            dof_indices = same_dof_range[0]->first;
            const Vector<double> b_values = dof_vector_to_b_values[dof_indices];
            internal::add_tangentiality_constraints (dof_indices,
                                                     average_tangent,
                                                     constraints,
                                                     b_values);
          }
          }
      }
  }



  namespace
  {
    template <int dim>
    struct PointComparator
    {
      bool operator ()(const std_cxx11::array<types::global_dof_index,dim> &p1,
                       const std_cxx11::array<types::global_dof_index,dim> &p2)
      {
        for (unsigned int d=0; d<dim; ++d)
          if (p1[d] < p2[d])
            return true;
        return false;
      }
    };
  }



  template <int dim, template <int, int> class DH, int spacedim>
  void
  compute_normal_flux_constraints (const DH<dim,spacedim> &dof_handler,
                                   const unsigned int     first_vector_component,
                                   const std::set<types::boundary_id> &boundary_ids,
                                   ConstraintMatrix      &constraints,
                                   const Mapping<dim, spacedim> &mapping)
  {
    ZeroFunction<dim>zero_function(dim);
    typename FunctionMap<spacedim>::type function_map;
    std::set<types::boundary_id>::const_iterator it
      = boundary_ids.begin();
    for (; it != boundary_ids.end(); ++it)
      function_map[*it] = &zero_function;
    compute_nonzero_tangential_flux_constraints(dof_handler,
                                                first_vector_component,
                                                boundary_ids,
                                                function_map,
                                                constraints,
                                                mapping);
  }

  template <int dim, template <int, int> class DH, int spacedim>
  void
  compute_nonzero_tangential_flux_constraints (const DH<dim,spacedim> &dof_handler,
                                               const unsigned int     first_vector_component,
                                               const std::set<types::boundary_id> &boundary_ids,
                                               typename FunctionMap<spacedim>::type &function_map,
                                               ConstraintMatrix      &constraints,
                                               const Mapping<dim, spacedim> &mapping)
  {
    ConstraintMatrix no_normal_flux_constraints(constraints.get_local_lines());
    compute_nonzero_normal_flux_constraints (dof_handler,
                                             first_vector_component,
                                             boundary_ids,
                                             function_map,
                                             no_normal_flux_constraints,
                                             mapping);

    hp::FECollection<dim,spacedim>      fe_collection (dof_handler.get_fe());
    hp::MappingCollection<dim,spacedim> mapping_collection;
    for (unsigned int i=0; i<fe_collection.size(); ++i)
      mapping_collection.push_back (mapping);

    // now also create a quadrature collection for the faces of a cell. fill
    // it with a quadrature formula with the support points on faces for each
    // FE
    hp::QCollection<dim-1> face_quadrature_collection;
    for (unsigned int i=0; i<fe_collection.size(); ++i)
      {
        const std::vector<Point<dim-1> > &
        unit_support_points = fe_collection[i].get_unit_face_support_points();

        Assert (unit_support_points.size() == fe_collection[i].dofs_per_face,
                ExcInternalError());

        face_quadrature_collection.push_back (Quadrature<dim-1> (unit_support_points));
      }

    // now create the object with which we will generate the normal vectors
    hp::FEFaceValues<dim,spacedim> x_fe_face_values (mapping_collection,
                                                     fe_collection,
                                                     face_quadrature_collection,
                                                     update_q_points |
                                                     update_normal_vectors);

    // Extract a list that collects all vector components that belong to the
    // same node (scalar basis function). When creating that list, we use an
    // array of dim components that stores the global degree of freedom.
    std::set<std_cxx11::array<types::global_dof_index,dim>, PointComparator<dim> > vector_dofs;
    std::vector<types::global_dof_index> face_dofs;

    std::map<std_cxx11::array<types::global_dof_index,dim>, Vector<double> >
    dof_vector_to_b_values;

    std::set<types::boundary_id>::iterator b_id;
    std::vector<std_cxx11::array<types::global_dof_index,dim> > cell_vector_dofs;
    for (typename DH<dim,spacedim>::active_cell_iterator cell =
           dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      if (!cell->is_artificial())
        for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if ((b_id=boundary_ids.find(cell->face(face_no)->boundary_indicator()))
              != boundary_ids.end())
            {
              const FiniteElement<dim> &fe = cell->get_fe();
              typename DH<dim,spacedim>::face_iterator face=cell->face(face_no);

              // get the indices of the dofs on this cell...
              face_dofs.resize (fe.dofs_per_face);
              face->get_dof_indices (face_dofs, cell->active_fe_index());

              x_fe_face_values.reinit (cell, face_no);
              const FEFaceValues<dim> &fe_values = x_fe_face_values.get_present_fe_values();

              std::map<types::global_dof_index, double> dof_to_b_value;

              unsigned int n_scalar_indices = 0;
              cell_vector_dofs.resize(fe.dofs_per_face);
              for (unsigned int i=0; i<fe.dofs_per_face; ++i)
                {
                  if (fe.face_system_to_component_index(i).first >= first_vector_component &&
                      fe.face_system_to_component_index(i).first < first_vector_component + dim)
                    {
                      const unsigned int component
                        = fe.face_system_to_component_index(i).first
                          -first_vector_component;
                      n_scalar_indices =
                        std::max(n_scalar_indices,
                                 fe.face_system_to_component_index(i).second+1);
                      cell_vector_dofs[fe.face_system_to_component_index(i).second]
                      [component]
                        = face_dofs[i];

                      const Point<dim> point
                        = fe_values.quadrature_point(i);
                      const double b_value
                        = function_map[*b_id]->value(point, component);
                      dof_to_b_value.insert
                      (std::make_pair(face_dofs[i], b_value));
                    }
                }

              // now we identified the vector indices on the cell, so next
              // insert them into the set (it would be expensive to directly
              // insert incomplete points into the set)
              for (unsigned int i=0; i<n_scalar_indices; ++i)
                {
                  vector_dofs.insert(cell_vector_dofs[i]);
                  Vector<double> b_values(dim);
                  for (unsigned int j=0; j<dim; ++j)
                    b_values[j]=dof_to_b_value[cell_vector_dofs[i][j]];
                  dof_vector_to_b_values.insert
                  (std::make_pair(cell_vector_dofs[i], b_values));
                }

            }

    // iterate over the list of all vector components we found and see if we
    // can find constrained ones
    unsigned int n_total_constraints_found = 0;
    for (typename std::set<std_cxx11::array<types::global_dof_index,dim>,
         PointComparator<dim> >::const_iterator it=vector_dofs.begin();
         it!=vector_dofs.end(); ++it)
      {
        unsigned int n_constraints = 0;
        bool is_constrained[dim];
        for (unsigned int d=0; d<dim; ++d)
          if (no_normal_flux_constraints.is_constrained((*it)[d]))
            {
              is_constrained[d] = true;
              ++n_constraints;
              ++n_total_constraints_found;
            }
          else
            is_constrained[d] = false;
        if (n_constraints > 0)
          {
            // if more than one no-flux constraint is present, we need to
            // constrain all vector degrees of freedom (we are in a corner
            // where several faces meet and to get a continuous FE solution we
            // need to set all conditions corresponding to the boundary function.).
            if (n_constraints > 1)
              {
                const Vector<double> b_value = dof_vector_to_b_values[*it];
                for (unsigned int d=0; d<dim; ++d)
                  {
                    constraints.add_line((*it)[d]);
                    constraints.set_inhomogeneity((*it)[d], b_value(d));
                  }
                continue;
              }

            // ok, this is a no-flux constraint, so get the index of the dof
            // that is currently constrained and make it unconstrained. The
            // constraint indices will get the normal that contain the other
            // indices.
            Tensor<1,dim> normal;
            unsigned constrained_index = -1;
            for (unsigned int d=0; d<dim; ++d)
              if (is_constrained[d])
                {
                  constrained_index = d;
                  normal[d] = 1.;
                }
            AssertIndexRange(constrained_index, dim);
            const std::vector<std::pair<types::global_dof_index, double> > *constrained
              = no_normal_flux_constraints.get_constraint_entries((*it)[constrained_index]);
            // find components to which this index is constrained to
            Assert(constrained != 0, ExcInternalError());
            Assert(constrained->size() < dim, ExcInternalError());
            for (unsigned int c=0; c<constrained->size(); ++c)
              {
                int index = -1;
                for (unsigned int d=0; d<dim; ++d)
                  if ((*constrained)[c].first == (*it)[d])
                    index = d;
                Assert (index != -1, ExcInternalError());
                normal[index] = (*constrained)[c].second;
              }
            Vector<double> boundary_value = dof_vector_to_b_values[*it];
            for (unsigned int d=0; d<dim; ++d)
              {
                if (is_constrained[d])
                  continue;
                const unsigned int new_index = (*it)[d];
                if (!constraints.is_constrained(new_index))
                  {
                    constraints.add_line(new_index);
                    if (std::abs(normal[d]) > 1e-13)
                      constraints.add_entry(new_index, (*it)[constrained_index],
                                            -normal[d]);
                    constraints.set_inhomogeneity(new_index, boundary_value[d]);
                  }
              }
          }
      }
    AssertDimension(n_total_constraints_found,
                    no_normal_flux_constraints.n_constraints());
  }



  namespace internal
  {
    template <int dim, int spacedim>
    struct IDScratchData
    {
      IDScratchData (const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                     const dealii::hp::FECollection<dim,spacedim> &fe,
                     const dealii::hp::QCollection<dim> &q,
                     const UpdateFlags update_flags);

      IDScratchData (const IDScratchData &data);

      void resize_vectors (const unsigned int n_q_points,
                           const unsigned int n_components);

      std::vector<dealii::Vector<double> > function_values;
      std::vector<std::vector<Tensor<1,spacedim> > > function_grads;
      std::vector<double> weight_values;
      std::vector<dealii::Vector<double> > weight_vectors;

      std::vector<dealii::Vector<double> > psi_values;
      std::vector<std::vector<Tensor<1,spacedim> > > psi_grads;
      std::vector<double> psi_scalar;

      std::vector<double>         tmp_values;
      std::vector<Tensor<1,spacedim> > tmp_gradients;

      dealii::hp::FEValues<dim,spacedim> x_fe_values;
    };


    template <int dim, int spacedim>
    IDScratchData<dim,spacedim>
    ::IDScratchData(const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                    const dealii::hp::FECollection<dim,spacedim> &fe,
                    const dealii::hp::QCollection<dim> &q,
                    const UpdateFlags update_flags)
      :
      x_fe_values(mapping, fe, q, update_flags)
    {}

    template <int dim, int spacedim>
    IDScratchData<dim,spacedim>::IDScratchData (const IDScratchData &data)
      :
      x_fe_values(data.x_fe_values.get_mapping_collection(),
                  data.x_fe_values.get_fe_collection(),
                  data.x_fe_values.get_quadrature_collection(),
                  data.x_fe_values.get_update_flags())
    {}

    template <int dim, int spacedim>
    void
    IDScratchData<dim,spacedim>::resize_vectors (const unsigned int n_q_points,
                                                 const unsigned int n_components)
    {
      function_values.resize (n_q_points,
                              dealii::Vector<double>(n_components));
      function_grads.resize (n_q_points,
                             std::vector<Tensor<1,spacedim> >(n_components));

      weight_values.resize (n_q_points);
      weight_vectors.resize (n_q_points,
                             dealii::Vector<double>(n_components));

      psi_values.resize (n_q_points,
                         dealii::Vector<double>(n_components));
      psi_grads.resize (n_q_points,
                        std::vector<Tensor<1,spacedim> >(n_components));
      psi_scalar.resize (n_q_points);

      tmp_values.resize (n_q_points);
      tmp_gradients.resize (n_q_points);
    }


    // avoid compiling inner function for many vector types when we always
    // really do the same thing by putting the main work into this helper
    // function
    template <int dim, int spacedim>
    double
    integrate_difference_inner (const Function<spacedim>   &exact_solution,
                                const NormType              &norm,
                                const Function<spacedim>    *weight,
                                const UpdateFlags            update_flags,
                                const double                 exponent,
                                const unsigned int           n_components,
                                IDScratchData<dim,spacedim> &data)
    {
      const bool fe_is_system = (n_components != 1);
      const dealii::FEValues<dim, spacedim> &fe_values  = data.x_fe_values.get_present_fe_values ();
      const unsigned int n_q_points = fe_values.n_quadrature_points;

      if (weight!=0)
        {
          if (weight->n_components>1)
            weight->vector_value_list (fe_values.get_quadrature_points(),
                                       data.weight_vectors);
          else
            {
              weight->value_list (fe_values.get_quadrature_points(),
                                  data.weight_values);
              for (unsigned int k=0; k<n_q_points; ++k)
                data.weight_vectors[k] = data.weight_values[k];
            }
        }
      else
        {
          for (unsigned int k=0; k<n_q_points; ++k)
            data.weight_vectors[k] = 1.;
        }


      if (update_flags & update_values)
        {
          // first compute the exact solution (vectors) at the quadrature
          // points try to do this as efficient as possible by avoiding a
          // second virtual function call in case the function really has only
          // one component
          if (fe_is_system)
            exact_solution.vector_value_list (fe_values.get_quadrature_points(),
                                              data.psi_values);
          else
            {
              exact_solution.value_list (fe_values.get_quadrature_points(),
                                         data.tmp_values);
              for (unsigned int i=0; i<n_q_points; ++i)
                data.psi_values[i](0) = data.tmp_values[i];
            }

          // then subtract finite element fe_function
          for (unsigned int q=0; q<n_q_points; ++q)
            data.psi_values[q] -= data.function_values[q];
        }

      // Do the same for gradients, if required
      if (update_flags & update_gradients)
        {
          // try to be a little clever to avoid recursive virtual function
          // calls when calling gradient_list for functions that are really
          // scalar functions
          if (fe_is_system)
            exact_solution.vector_gradient_list (fe_values.get_quadrature_points(),
                                                 data.psi_grads);
          else
            {
              exact_solution.gradient_list (fe_values.get_quadrature_points(),
                                            data.tmp_gradients);
              for (unsigned int i=0; i<n_q_points; ++i)
                data.psi_grads[i][0] = data.tmp_gradients[i];
            }

          // then subtract finite element function_grads. We need to be
          // careful in the codimension one case, since there we only have
          // tangential gradients in the finite element function, not the full
          // gradient. This is taken care of, by subtracting the normal
          // component of the gradient from the exact function.
          if (update_flags & update_normal_vectors)
            for (unsigned int k=0; k<n_components; ++k)
              for (unsigned int q=0; q<n_q_points; ++q)
                data.psi_grads[q][k] -= (data.function_grads[q][k] +
                                         (data.psi_grads[q][k]* // (f.n) n
                                          fe_values.normal_vector(q))*
                                         fe_values.normal_vector(q));
          else
            for (unsigned int k=0; k<n_components; ++k)
              for (unsigned int q=0; q<n_q_points; ++q)
                data.psi_grads[q][k] -= data.function_grads[q][k];
        }

      double diff = 0;
      switch (norm)
        {
        case mean:
          // Compute values in quadrature points and integrate
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += data.psi_values[q](k) * data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }
          break;

        case Lp_norm:
        case L1_norm:
        case W1p_norm:
          // Compute values in quadrature points and integrate
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += std::pow(data.psi_values[q](k)*data.psi_values[q](k),
                                exponent/2.) * data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }

          // Compute the root only, if no derivative values are added later
          if (!(update_flags & update_gradients))
            diff = std::pow(diff, 1./exponent);
          break;

        case L2_norm:
        case H1_norm:
          // Compute values in quadrature points and integrate
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += data.psi_values[q](k) * data.psi_values[q](k) *
                       data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }
          // Compute the root only, if no derivative values are added later
          if (norm == L2_norm)
            diff=std::sqrt(diff);
          break;

        case Linfty_norm:
        case W1infty_norm:
          for (unsigned int q=0; q<n_q_points; ++q)
            for (unsigned int k=0; k<n_components; ++k)
              diff = std::max (diff, std::abs(data.psi_values[q](k)*
                                              data.weight_vectors[q](k)));
          break;

        case H1_seminorm:
        case W1p_seminorm:
        case W1infty_seminorm:
          break;

        default:
          Assert (false, ExcNotImplemented());
          break;
        }

      switch (norm)
        {
        case W1p_seminorm:
        case W1p_norm:
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += std::pow(data.psi_grads[q][k]*data.psi_grads[q][k],
                                exponent/2.) * data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }
          diff = std::pow(diff, 1./exponent);
          break;

        case H1_seminorm:
        case H1_norm:
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += (data.psi_grads[q][k] * data.psi_grads[q][k]) *
                       data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }
          diff = std::sqrt(diff);
          break;

        case W1infty_seminorm:
        case W1infty_norm:
        {
          double t = 0;
          for (unsigned int q=0; q<n_q_points; ++q)
            for (unsigned int k=0; k<n_components; ++k)
              for (unsigned int d=0; d<dim; ++d)
                t = std::max(t, std::abs(data.psi_grads[q][k][d]) *
                             data.weight_vectors[q](k));

          // then add seminorm to norm if that had previously been computed
          diff += t;
        }
        break;
        default:
          break;
        }

      // append result of this cell to the end of the vector
      Assert (numbers::is_finite(diff), ExcNumberNotFinite());
      return diff;
    }



    template <int dim, class InVector, class OutVector, class DH, int spacedim>
    static
    void
    do_integrate_difference (const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                             const DH              &dof,
                             const InVector        &fe_function,
                             const Function<spacedim>   &exact_solution,
                             OutVector             &difference,
                             const dealii::hp::QCollection<dim> &q,
                             const NormType &norm,
                             const Function<spacedim>   *weight,
                             const double           exponent_1)
    {
      // we mark the "exponent" parameter to this function "const" since it is
      // strictly incoming, but we need to set it to something different later
      // on, if necessary, so have a read-write version of it:
      double exponent = exponent_1;

      const unsigned int        n_components = dof.get_fe().n_components();

      if (weight!=0)
        {
          Assert ((weight->n_components==1) || (weight->n_components==n_components),
                  ExcDimensionMismatch(weight->n_components, n_components));
        }

      difference.reinit (dof.get_tria().n_active_cells());

      switch (norm)
        {
        case L2_norm:
        case H1_seminorm:
        case H1_norm:
          exponent = 2.;
          break;

        case L1_norm:
          exponent = 1.;
          break;

        default:
          break;
        }

      UpdateFlags update_flags = UpdateFlags (update_quadrature_points  |
                                              update_JxW_values);
      switch (norm)
        {
        case H1_seminorm:
        case W1p_seminorm:
        case W1infty_seminorm:
          update_flags |= UpdateFlags (update_gradients);
          if (spacedim == dim+1)
            update_flags |= UpdateFlags (update_normal_vectors);

          break;

        case H1_norm:
        case W1p_norm:
        case W1infty_norm:
          update_flags |= UpdateFlags (update_gradients);
          if (spacedim == dim+1)
            update_flags |= UpdateFlags (update_normal_vectors);
        // no break!

        default:
          update_flags |= UpdateFlags (update_values);
          break;
        }

      dealii::hp::FECollection<dim,spacedim> fe_collection (dof.get_fe());
      IDScratchData<dim,spacedim> data(mapping, fe_collection, q, update_flags);

      // loop over all cells
      typename DH::active_cell_iterator cell = dof.begin_active(),
                                        endc = dof.end();
      for (unsigned int index=0; cell != endc; ++cell, ++index)
        if (cell->is_locally_owned())
          {
            // initialize for this cell
            data.x_fe_values.reinit (cell);

            const dealii::FEValues<dim, spacedim> &fe_values  = data.x_fe_values.get_present_fe_values ();
            const unsigned int   n_q_points = fe_values.n_quadrature_points;
            data.resize_vectors (n_q_points, n_components);

            if (update_flags & update_values)
              fe_values.get_function_values (fe_function, data.function_values);
            if (update_flags & update_gradients)
              fe_values.get_function_grads (fe_function, data.function_grads);

            difference(index) =
              integrate_difference_inner (exact_solution, norm, weight,
                                          update_flags, exponent,
                                          n_components, data);
          }
        else
          // the cell is a ghost cell or is artificial. write a zero into the
          // corresponding value of the returned vector
          difference(index) = 0;
    }

  } // namespace internal




  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference (const Mapping<dim, spacedim>    &mapping,
                        const DoFHandler<dim,spacedim> &dof,
                        const InVector        &fe_function,
                        const Function<spacedim>   &exact_solution,
                        OutVector             &difference,
                        const Quadrature<dim> &q,
                        const NormType        &norm,
                        const Function<spacedim>   *weight,
                        const double           exponent)
  {
    internal
    ::do_integrate_difference (hp::MappingCollection<dim,spacedim>(mapping),
                               dof, fe_function, exact_solution,
                               difference, hp::QCollection<dim>(q),
                               norm, weight, exponent);
  }


  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference (const DoFHandler<dim,spacedim>    &dof,
                        const InVector           &fe_function,
                        const Function<spacedim>      &exact_solution,
                        OutVector                &difference,
                        const Quadrature<dim>    &q,
                        const NormType           &norm,
                        const Function<spacedim>      *weight,
                        const double              exponent)
  {
    internal
    ::do_integrate_difference(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                              dof, fe_function, exact_solution,
                              difference, hp::QCollection<dim>(q),
                              norm, weight, exponent);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference (const dealii::hp::MappingCollection<dim,spacedim>    &mapping,
                        const dealii::hp::DoFHandler<dim,spacedim> &dof,
                        const InVector        &fe_function,
                        const Function<spacedim>   &exact_solution,
                        OutVector             &difference,
                        const dealii::hp::QCollection<dim> &q,
                        const NormType        &norm,
                        const Function<spacedim>   *weight,
                        const double           exponent)
  {
    internal
    ::do_integrate_difference (hp::MappingCollection<dim,spacedim>(mapping),
                               dof, fe_function, exact_solution,
                               difference, q,
                               norm, weight, exponent);
  }


  template <int dim, class InVector, class OutVector, int spacedim>
  void
  integrate_difference (const dealii::hp::DoFHandler<dim,spacedim>    &dof,
                        const InVector           &fe_function,
                        const Function<spacedim>      &exact_solution,
                        OutVector                &difference,
                        const dealii::hp::QCollection<dim>    &q,
                        const NormType           &norm,
                        const Function<spacedim>      *weight,
                        const double              exponent)
  {
    internal
    ::do_integrate_difference(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                              dof, fe_function, exact_solution,
                              difference, q,
                              norm, weight, exponent);
  }



  template <int dim, class InVector, int spacedim>
  void
  point_difference (const DoFHandler<dim,spacedim> &dof,
                    const InVector        &fe_function,
                    const Function<spacedim>   &exact_function,
                    Vector<double>        &difference,
                    const Point<spacedim>      &point)
  {
    point_difference(StaticMappingQ1<dim>::mapping,
                     dof,
                     fe_function,
                     exact_function,
                     difference,
                     point);
  }


  template <int dim, class InVector, int spacedim>
  void
  point_difference (const Mapping<dim, spacedim>    &mapping,
                    const DoFHandler<dim,spacedim> &dof,
                    const InVector        &fe_function,
                    const Function<spacedim>   &exact_function,
                    Vector<double>        &difference,
                    const Point<spacedim>      &point)
  {
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(difference.size() == fe.n_components(),
           ExcDimensionMismatch(difference.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point = GridTools::find_active_cell_around_point (mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell_point.first);

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    if (fe.n_components() == 1)
      difference(0) = exact_function.value(point);
    else
      exact_function.vector_value(point, difference);

    for (unsigned int i=0; i<difference.size(); ++i)
      difference(i) -= u_value[0](i);
  }


  template <int dim, class InVector, int spacedim>
  void
  point_value (const DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point,
               Vector<double>        &value)
  {

    point_value (StaticMappingQ1<dim,spacedim>::mapping,
                 dof,
                 fe_function,
                 point,
                 value);
  }


  template <int dim, class InVector, int spacedim>
  void
  point_value (const hp::DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point,
               Vector<double>        &value)
  {
    point_value(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                dof,
                fe_function,
                point,
                value);
  }


  template <int dim, class InVector, int spacedim>
  double
  point_value (const DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point)
  {
    return point_value (StaticMappingQ1<dim,spacedim>::mapping,
                        dof,
                        fe_function,
                        point);
  }


  template <int dim, class InVector, int spacedim>
  double
  point_value (const hp::DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point)
  {
    return point_value(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                       dof,
                       fe_function,
                       point);
  }


  template <int dim, class InVector, int spacedim>
  void
  point_value (const Mapping<dim, spacedim>    &mapping,
               const DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point,
               Vector<double>        &value)
  {
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(value.size() == fe.n_components(),
           ExcDimensionMismatch(value.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point
      = GridTools::find_active_cell_around_point (mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell_point.first);

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, class InVector, int spacedim>
  void
  point_value (const hp::MappingCollection<dim, spacedim>    &mapping,
               const hp::DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point,
               Vector<double>        &value)
  {
    const hp::FECollection<dim, spacedim> &fe = dof.get_fe();

    Assert(value.size() == fe.n_components(),
           ExcDimensionMismatch(value.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =  GridTools::find_active_cell_around_point (mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    hp::FEValues<dim, spacedim> hp_fe_values(mapping, fe, hp::QCollection<dim>(quadrature), update_values);
    hp_fe_values.reinit(cell_point.first);
    const FEValues<dim, spacedim> &fe_values = hp_fe_values.get_present_fe_values();

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, class InVector, int spacedim>
  double
  point_value (const Mapping<dim, spacedim>    &mapping,
               const DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point)
  {
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(fe.n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point = GridTools::find_active_cell_around_point (mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell_point.first);

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<double> u_value(1);
    fe_values.get_function_values(fe_function, u_value);

    return u_value[0];
  }


  template <int dim, class InVector, int spacedim>
  double
  point_value (const hp::MappingCollection<dim, spacedim>    &mapping,
               const hp::DoFHandler<dim,spacedim> &dof,
               const InVector        &fe_function,
               const Point<spacedim>      &point)
  {
    const hp::FECollection<dim, spacedim> &fe = dof.get_fe();

    Assert(fe.n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point =  GridTools::find_active_cell_around_point (mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    hp::FEValues<dim, spacedim> hp_fe_values(mapping, fe, hp::QCollection<dim>(quadrature), update_values);
    hp_fe_values.reinit(cell_point.first);
    const FEValues<dim, spacedim> &fe_values = hp_fe_values.get_present_fe_values();

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<double> u_value(1);
    fe_values.get_function_values(fe_function, u_value);

    return u_value[0];
  }



  template <class VECTOR>
  void
  subtract_mean_value(VECTOR                  &v,
                      const std::vector<bool> &p_select)
  {
    if (p_select.size() == 0)
      {
        // In case of an empty boolean mask operate on the whole vector:
        v.add( - v.mean_value() );
      }
    else
      {
        // This function is not implemented for distributed vectors, so
        // if v is not a boring Vector or BlockVector:
        Assert(   dynamic_cast<Vector<double> *>(& v)
                  || dynamic_cast<Vector<float> *>(& v)
                  || dynamic_cast<Vector<long double> *>(& v)
                  || dynamic_cast<BlockVector<double> *>(& v)
                  || dynamic_cast<BlockVector<float> *>(& v)
                  || dynamic_cast<BlockVector<long double> *>(& v),
                  ExcNotImplemented());

        const unsigned int n = v.size();

        Assert(p_select.size() == n,
               ExcDimensionMismatch(p_select.size(), n));

        typename VECTOR::value_type s = 0.;
        unsigned int counter = 0;
        for (unsigned int i=0; i<n; ++i)
          if (p_select[i])
            {
              s += v(i);
              ++counter;
            }
        // Error out if we have not constrained anything. Note that in this
        // case the vector v is always nonempty.
        Assert (n == 0 || counter > 0, ComponentMask::ExcNoComponentSelected());

        s /= counter;

        for (unsigned int i=0; i<n; ++i)
          if (p_select[i])
            v(i) -= s;
      }
  }



  template <int dim, class InVector, int spacedim>
  double
  compute_mean_value (const Mapping<dim, spacedim>    &mapping,
                      const DoFHandler<dim,spacedim> &dof,
                      const Quadrature<dim> &quadrature,
                      const InVector        &v,
                      const unsigned int     component)
  {
    Assert (v.size() == dof.n_dofs(),
            ExcDimensionMismatch (v.size(), dof.n_dofs()));
    Assert (component < dof.get_fe().n_components(),
            ExcIndexRange(component, 0, dof.get_fe().n_components()));

    FEValues<dim,spacedim> fe(mapping, dof.get_fe(), quadrature,
                              UpdateFlags(update_JxW_values
                                          | update_values));

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell;
    std::vector<Vector<double> > values(quadrature.size(),
                                        Vector<double> (dof.get_fe().n_components()));

    double mean = 0.;
    double area = 0.;
    // Compute mean value
    for (cell = dof.begin_active(); cell != dof.end(); ++cell)
      if (cell->is_locally_owned())
        {
          fe.reinit (cell);
          fe.get_function_values(v, values);
          for (unsigned int k=0; k< quadrature.size(); ++k)
            {
              mean += fe.JxW(k) * values[k](component);
              area += fe.JxW(k);
            }
        }

#ifdef DEAL_II_WITH_P4EST
    // if this was a distributed DoFHandler, we need to do the reduction
    // over the entire domain
    if (const parallel::distributed::Triangulation<dim,spacedim> *
        p_d_triangulation
        = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim> *>(&dof.get_tria()))
      {
        double my_values[2] = { mean, area };
        double global_values[2];

        MPI_Allreduce (&my_values, &global_values, 2, MPI_DOUBLE,
                       MPI_SUM,
                       p_d_triangulation->get_communicator());

        mean = global_values[0];
        area = global_values[1];
      }
#endif

    return (mean/area);
  }


  template <int dim, class InVector, int spacedim>
  double
  compute_mean_value (const DoFHandler<dim,spacedim> &dof,
                      const Quadrature<dim> &quadrature,
                      const InVector        &v,
                      const unsigned int     component)
  {
    return compute_mean_value(StaticMappingQ1<dim,spacedim>::mapping, dof, quadrature, v, component);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
