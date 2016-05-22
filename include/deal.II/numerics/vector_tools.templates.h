// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2016 by the deal.II authors
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


#ifndef dealii__vector_tools_templates_h
#define dealii__vector_tools_templates_h

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
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
#include <deal.II/fe/mapping_q.h>
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

  template <typename VectorType, int dim, int spacedim,
            template <int, int> class DoFHandlerType>
  void interpolate (const Mapping<dim,spacedim>        &mapping,
                    const DoFHandlerType<dim,spacedim> &dof,
                    const Function<spacedim, typename VectorType::value_type>           &function,
                    VectorType                         &vec)
  {
    typedef typename VectorType::value_type number;
    Assert (vec.size() == dof.n_dofs(),
            ExcDimensionMismatch (vec.size(), dof.n_dofs()));
    Assert (dof.get_fe().n_components() == function.n_components,
            ExcDimensionMismatch(dof.get_fe().n_components(),
                                 function.n_components));

    const hp::FECollection<dim,spacedim> fe (dof.get_fe());
    const unsigned int          n_components = fe.n_components();
    const bool                  fe_is_system = (n_components != 1);

    typename DoFHandlerType<dim,spacedim>::active_cell_iterator
    cell = dof.begin_active(),
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
        Assert ((unit_support_points[fe_index].size() != 0)||(fe[fe_index].dofs_per_cell == 0),
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
                dof_to_rep_index_table[fe_index].push_back(dofs_of_rep_points[fe_index].size());
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
    std::vector<std::vector<number> >         function_values_scalar(fe.size());
    std::vector<std::vector<Vector<number> > > function_values_system(fe.size());

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
          if (fe[fe_index].dofs_per_cell != 0)
            {
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
                           Vector<number> (fe[fe_index].n_components()));
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
        }
    vec.compress(VectorOperation::insert);
  }


  template <typename VectorType, typename DoFHandlerType>
  void interpolate (const DoFHandlerType                            &dof,
                    const Function<DoFHandlerType::space_dimension,typename VectorType::value_type> &function,
                    VectorType                                      &vec)
  {
    interpolate(StaticMappingQ1<DoFHandlerType::dimension,
                DoFHandlerType::space_dimension>::mapping,
                dof,
                function,
                vec);
  }




  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolate (const DoFHandler<dim,spacedim>           &dof_1,
               const DoFHandler<dim,spacedim>           &dof_2,
               const FullMatrix<double>        &transfer,
               const InVector                  &data_1,
               OutVector                       &data_2)
  {
    typedef typename OutVector::value_type number;
    Vector<number> cell_data_1(dof_1.get_fe().dofs_per_cell);
    Vector<number> cell_data_2(dof_2.get_fe().dofs_per_cell);

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


  template<typename VectorType, typename DoFHandlerType>
  void
  interpolate_based_on_material_id
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &mapping,
   const DoFHandlerType                                                      &dof,
   const std::map<types::material_id, const Function<DoFHandlerType::space_dimension, typename VectorType::value_type> *> &function_map,
   VectorType                                                                &dst,
   const ComponentMask                                                       &component_mask)
  {
    typedef typename VectorType::value_type number;
    const unsigned int dim = DoFHandlerType::dimension;

    Assert( component_mask.represents_n_components(dof.get_fe().n_components()),
            ExcMessage("The number of components in the mask has to be either "
                       "zero or equal to the number of components in the finite "
                       "element.") );

    if ( function_map.size() == 0 )
      return;

    Assert( function_map.find(numbers::invalid_material_id) == function_map.end(),
            ExcMessage("You cannot specify the invalid material indicator "
                       "in your function map."));

    for (typename std::map<types::material_id, const Function<DoFHandlerType::space_dimension, number>* >
         ::const_iterator
         iter  = function_map.begin();
         iter != function_map.end();
         ++iter )
      {
        Assert( dof.get_fe().n_components() == iter->second->n_components,
                ExcDimensionMismatch(dof.get_fe().n_components(), iter->second->n_components) );
      }

    const hp::FECollection<DoFHandlerType::dimension, DoFHandlerType::space_dimension>
    fe(dof.get_fe());
    const unsigned int n_components =  fe.n_components();
    const bool         fe_is_system = (n_components != 1);

    typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
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
    std::vector< types::global_dof_index> dofs_on_cell(fe.max_dofs_per_cell());
    std::vector< Point<DoFHandlerType::space_dimension> > rep_points(max_rep_points);

    std::vector< std::vector<number> >           function_values_scalar(fe.size());
    std::vector< std::vector< Vector<number> > > function_values_system(fe.size());

    hp::QCollection<dim> support_quadrature;
    for (unsigned int fe_index = 0; fe_index < fe.size(); ++fe_index)
      support_quadrature.push_back( Quadrature<dim>(unit_support_points[fe_index]) );

    hp::MappingCollection<dim, DoFHandlerType::space_dimension> mapping_collection(mapping);
    hp::FEValues<dim, DoFHandlerType::space_dimension> fe_values(mapping_collection,
        fe,
        support_quadrature,
        update_quadrature_points);

    for ( ; cell != endc; ++cell)
      if ( cell->is_locally_owned() )
        if ( function_map.find(cell->material_id()) != function_map.end() )
          {
            const unsigned int fe_index = cell->active_fe_index();

            fe_values.reinit(cell);

            const std::vector< Point<DoFHandlerType::space_dimension> > &support_points
              = fe_values.get_present_fe_values().get_quadrature_points();

            rep_points.resize( dofs_of_rep_points[fe_index].size() );
            for (unsigned int i = 0; i < dofs_of_rep_points[fe_index].size(); ++i)
              rep_points[i] = support_points[dofs_of_rep_points[fe_index][i]];

            dofs_on_cell.resize( fe[fe_index].dofs_per_cell );
            cell->get_dof_indices(dofs_on_cell);

            if (fe_is_system)
              {
                function_values_system[fe_index].resize( n_rep_points[fe_index],
                                                         Vector<number>(fe[fe_index].n_components()) );

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
     * Interpolate zero boundary values. We don't need to worry about a
     * mapping here because the function we evaluate for the DoFs is zero in
     * the mapped locations as well as in the original, unmapped locations
     */
    template <typename DoFHandlerType, typename number>
    void
    interpolate_zero_boundary_values (const DoFHandlerType                     &dof_handler,
                                      std::map<types::global_dof_index,number> &boundary_values)
    {
      const unsigned int dim = DoFHandlerType::dimension;

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
      typename DoFHandlerType::active_cell_iterator
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



  template <int dim, int spacedim, template <int, int> class DoFHandlerType, typename VectorType>
  void
  interpolate_to_different_mesh (const DoFHandlerType<dim, spacedim> &dof1,
                                 const VectorType                    &u1,
                                 const DoFHandlerType<dim, spacedim> &dof2,
                                 VectorType                          &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage ("The two DoF handlers must represent triangulations that "
                       "have the same coarse meshes"));

    InterGridMap<DoFHandlerType<dim, spacedim> > intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    ConstraintMatrix dummy;
    dummy.close();

    interpolate_to_different_mesh(intergridmap, u1, dummy, u2);
  }



  template <int dim, int spacedim, template <int, int> class DoFHandlerType, typename VectorType>
  void
  interpolate_to_different_mesh (const DoFHandlerType<dim, spacedim> &dof1,
                                 const VectorType                    &u1,
                                 const DoFHandlerType<dim, spacedim> &dof2,
                                 const ConstraintMatrix              &constraints,
                                 VectorType                          &u2)
  {
    Assert(GridTools::have_same_coarse_mesh(dof1, dof2),
           ExcMessage ("The two DoF handlers must represent triangulations that "
                       "have the same coarse meshes"));

    InterGridMap<DoFHandlerType<dim, spacedim> > intergridmap;
    intergridmap.make_mapping(dof1, dof2);

    interpolate_to_different_mesh(intergridmap, u1, constraints, u2);
  }

  namespace internal
  {
    /**
     * Returns whether the cell and all of its descendants are locally owned.
     */
    template <typename cell_iterator>
    bool is_locally_owned(const cell_iterator &cell)
    {
      if (cell->active())
        return cell->is_locally_owned();

      for (unsigned int c=0; c<cell->n_children(); ++c)
        if (!is_locally_owned(cell->child(c)))
          return false;

      return true;
    }
  }


  template <int dim, int spacedim, template <int, int> class DoFHandlerType, typename VectorType>
  void
  interpolate_to_different_mesh
  (const InterGridMap<DoFHandlerType<dim, spacedim> > &intergridmap,
   const VectorType       &u1,
   const ConstraintMatrix &constraints,
   VectorType             &u2)
  {
    const DoFHandlerType<dim, spacedim> &dof1 = intergridmap.get_source_grid();
    const DoFHandlerType<dim, spacedim> &dof2 = intergridmap.get_destination_grid();
    (void)dof2;

    Assert(u1.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(),dof2.n_dofs()));

    Vector<typename VectorType::value_type> cache;

    // Looping over the finest common
    // mesh, this means that source and
    // destination cells have to be on the
    // same level and at least one has to
    // be active.
    //
    // Therefor, loop over all cells
    // (active and inactive) of the source
    // grid ..
    typename DoFHandlerType<dim,spacedim>::cell_iterator       cell1 = dof1.begin();
    const typename DoFHandlerType<dim,spacedim>::cell_iterator endc1 = dof1.end();

    for (; cell1 != endc1; ++cell1)
      {
        const typename DoFHandlerType<dim,spacedim>::cell_iterator cell2 = intergridmap[cell1];

        // .. and skip if source and destination
        // cells are not on the same level ..
        if (cell1->level() != cell2->level())
          continue;
        // .. or none of them is active.
        if (!cell1->active() && !cell2->active())
          continue;

        Assert(internal::is_locally_owned(cell1) == internal::is_locally_owned(cell2),
               ExcMessage("The two Triangulations are required to have the same parallel partitioning."));

        // Skip foreign cells.
        if (cell1->active() && !cell1->is_locally_owned())
          continue;
        if (cell2->active() && !cell2->is_locally_owned())
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

    // finish the work on parallel vectors
    u2.compress (VectorOperation::insert);
    // Apply hanging node constraints.
    constraints.distribute(u2);
  }


  namespace
  {
    /**
     * Compute the boundary values to be used in the project() functions.
     */
    template <int dim, int spacedim, template <int, int> class DoFHandlerType,
              template <int,int> class M_or_MC, template <int> class Q_or_QC,
              typename number>
    void project_compute_b_v
    (const M_or_MC<dim, spacedim>             &mapping,
     const DoFHandlerType<dim,spacedim>       &dof,
     const Function<spacedim,number>          &function,
     const bool                                enforce_zero_boundary,
     const Q_or_QC<dim-1>                     &q_boundary,
     const bool                                project_to_boundary_first,
     std::map<types::global_dof_index,number> &boundary_values)
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
            used_boundary_ids = dof.get_triangulation().get_boundary_ids();

            std::map<types::boundary_id, const Function<spacedim,number>*> boundary_functions;
            for (unsigned int i=0; i<used_boundary_ids.size(); ++i)
              boundary_functions[used_boundary_ids[i]] = &function;
            project_boundary_values (mapping, dof, boundary_functions, q_boundary,
                                     boundary_values);
          }
    }


    /**
     * Return whether the boundary values try to constrain a degree of freedom
     * that is already constrained to something else
     */
    template <typename number>
    bool constraints_and_b_v_are_compatible (const ConstraintMatrix   &constraints,
                                             std::map<types::global_dof_index,number> &boundary_values)
    {
      for (typename std::map<types::global_dof_index,number>::iterator it=boundary_values.begin();
           it != boundary_values.end(); ++it)
        if (constraints.is_constrained(it->first))
//TODO: This looks wrong -- shouldn't it be ==0 in the first condition and && ?
          if (!(constraints.get_constraint_entries(it->first)->size() > 0
                ||
                (constraints.get_inhomogeneity(it->first) == it->second)))
            return false;

      return true;
    }

    template <typename number>
    void invert_mass_matrix(const SparseMatrix<number> &mass_matrix,
                            const Vector<number> &rhs,
                            Vector<number> &solution)
    {
      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices
      ReductionControl      control(5*rhs.size(), 0., 1e-12, false, false);
      GrowingVectorMemory<Vector<number> > memory;
      SolverCG<Vector<number> >            cg(control,memory);

      PreconditionSSOR<SparseMatrix<number> > prec;
      prec.initialize(mass_matrix, 1.2);

      cg.solve (mass_matrix, solution, rhs, prec);
    }

    template <typename number>
    void invert_mass_matrix(const SparseMatrix<std::complex<number> > &/*mass_matrix*/,
                            const Vector<std::complex<number> > &/*rhs*/,
                            Vector<std::complex<number> > &/*solution*/)
    {
      Assert(false, ExcNotImplemented());
    }


    /**
     * Generic implementation of the project() function
     */
    template <int dim, int spacedim, typename VectorType,
              template <int, int> class DoFHandlerType,
              template <int,int> class M_or_MC, template <int> class Q_or_QC>
    void do_project (const M_or_MC<dim, spacedim>       &mapping,
                     const DoFHandlerType<dim,spacedim> &dof,
                     const ConstraintMatrix             &constraints,
                     const Q_or_QC<dim>                 &quadrature,
                     const Function<spacedim, typename VectorType::value_type>           &function,
                     VectorType                         &vec_result,
                     const bool                          enforce_zero_boundary,
                     const Q_or_QC<dim-1>               &q_boundary,
                     const bool                          project_to_boundary_first)
    {
      typedef typename VectorType::value_type number;
      Assert (dof.get_fe().n_components() == function.n_components,
              ExcDimensionMismatch(dof.get_fe().n_components(),
                                   function.n_components));
      Assert (vec_result.size() == dof.n_dofs(),
              ExcDimensionMismatch (vec_result.size(), dof.n_dofs()));

      // make up boundary values
      std::map<types::global_dof_index,number> boundary_values;
      project_compute_b_v(mapping, dof, function, enforce_zero_boundary,
                          q_boundary, project_to_boundary_first, boundary_values);

      // check if constraints are compatible (see below)
      const bool constraints_are_compatible =
        constraints_and_b_v_are_compatible<number>(constraints, boundary_values);

      // set up mass matrix and right hand side
      Vector<number> vec (dof.n_dofs());
      SparsityPattern sparsity;
      {
        DynamicSparsityPattern dsp (dof.n_dofs(), dof.n_dofs());
        DoFTools::make_sparsity_pattern (dof, dsp, constraints,
                                         !constraints_are_compatible);

        sparsity.copy_from (dsp);
      }
      SparseMatrix<number> mass_matrix (sparsity);
      Vector<number> tmp (mass_matrix.n());

      // If the constraint matrix does not conflict with the given boundary
      // values (i.e., it either does not contain boundary values or it contains
      // the same as boundary_values), we can let it call
      // distribute_local_to_global straight away, otherwise we need to first
      // interpolate the boundary values and then condense the matrix and vector
      if (constraints_are_compatible)
        {
          const Function<spacedim,number> *dummy = 0;
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

      invert_mass_matrix(mass_matrix,tmp,vec);
      constraints.distribute (vec);

      // copy vec into vec_result. we can't use vec_result itself above, since
      // it may be of another type than Vector<double> and that wouldn't
      // necessarily go together with the matrix and other functions
      for (unsigned int i=0; i<vec.size(); ++i)
        vec_result(i) = vec(i);
    }
  }


  template <int dim, typename VectorType, int spacedim>
  void project (const Mapping<dim, spacedim>   &mapping,
                const DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix         &constraints,
                const Quadrature<dim>          &quadrature,
                const Function<spacedim, typename VectorType::value_type>       &function,
                VectorType                     &vec_result,
                const bool                     enforce_zero_boundary,
                const Quadrature<dim-1>        &q_boundary,
                const bool                     project_to_boundary_first)
  {
    do_project (mapping, dof, constraints, quadrature,
                function, vec_result,
                enforce_zero_boundary, q_boundary,
                project_to_boundary_first);
  }


  template <int dim, typename VectorType, int spacedim>
  void project (const DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix         &constraints,
                const Quadrature<dim>          &quadrature,
                const Function<spacedim, typename VectorType::value_type>       &function,
                VectorType                     &vec,
                const bool                     enforce_zero_boundary,
                const Quadrature<dim-1>        &q_boundary,
                const bool                     project_to_boundary_first)
  {
    project(StaticMappingQ1<dim,spacedim>::mapping, dof, constraints, quadrature, function, vec,
            enforce_zero_boundary, q_boundary, project_to_boundary_first);
  }



  template <int dim, typename VectorType, int spacedim>
  void project (const hp::MappingCollection<dim, spacedim> &mapping,
                const hp::DoFHandler<dim,spacedim>         &dof,
                const ConstraintMatrix                     &constraints,
                const hp::QCollection<dim>                 &quadrature,
                const Function<spacedim, typename VectorType::value_type>                   &function,
                VectorType                                 &vec_result,
                const bool                                 enforce_zero_boundary,
                const hp::QCollection<dim-1>               &q_boundary,
                const bool                                 project_to_boundary_first)
  {
    do_project (mapping, dof, constraints, quadrature,
                function, vec_result,
                enforce_zero_boundary, q_boundary,
                project_to_boundary_first);
  }


  template <int dim, typename VectorType, int spacedim>
  void project (const hp::DoFHandler<dim,spacedim> &dof,
                const ConstraintMatrix             &constraints,
                const hp::QCollection<dim>         &quadrature,
                const Function<spacedim, typename VectorType::value_type>           &function,
                VectorType                         &vec,
                const bool                         enforce_zero_boundary,
                const hp::QCollection<dim-1>       &q_boundary,
                const bool                         project_to_boundary_first)
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
                                   const std::set<types::boundary_id> &boundary_ids)
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
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
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
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
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
                                   const std::set<types::boundary_id> &boundary_ids)
  {
    create_boundary_right_hand_side(StaticMappingQ1<dim>::mapping, dof_handler,
                                    quadrature,
                                    rhs_function, rhs_vector,
                                    boundary_ids);
  }



  template <int dim, int spacedim>
  void
  create_boundary_right_hand_side (const hp::MappingCollection<dim,spacedim> &mapping,
                                   const hp::DoFHandler<dim,spacedim> &dof_handler,
                                   const hp::QCollection<dim-1>  &quadrature,
                                   const Function<spacedim>      &rhs_function,
                                   Vector<double>                &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_ids)
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
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
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
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
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
                                   const std::set<types::boundary_id> &boundary_ids)
  {
    create_boundary_right_hand_side(hp::StaticMappingQ1<dim>::mapping_collection,
                                    dof_handler, quadrature,
                                    rhs_function, rhs_vector,
                                    boundary_ids);
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
    template <typename number, typename DoFHandlerType, template <int,int> class M_or_MC>
    static inline
    void do_interpolate_boundary_values
    (const M_or_MC<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &,
     const DoFHandlerType                                                                        &dof,
     const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
     std::map<types::global_dof_index,number>                                                    &boundary_values,
     const ComponentMask                                                                         &component_mask,
     const dealii::internal::int2type<1>)
    {
      const unsigned int dim = DoFHandlerType::dimension;
      const unsigned int spacedim = DoFHandlerType::space_dimension;

      Assert (component_mask.represents_n_components(dof.get_fe().n_components()),
              ExcMessage ("The number of components in the mask has to be either "
                          "zero or equal to the number of components in the finite "
                          "element."));

      // if for whatever reason we were
      // passed an empty map, return
      // immediately
      if (function_map.size() == 0)
        return;

      for (typename DoFHandlerType::active_cell_iterator cell = dof.begin_active();
           cell != dof.end(); ++cell)
        for (unsigned int direction=0;
             direction<GeometryInfo<dim>::faces_per_cell; ++direction)
          if (cell->at_boundary(direction)
              &&
              (function_map.find(cell->face(direction)->boundary_id()) != function_map.end()))
            {
              const Function<DoFHandlerType::space_dimension,number> &boundary_function
                = *function_map.find(cell->face(direction)->boundary_id())->second;

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
              Vector<number> function_values (fe.n_components());
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
    template <typename number, typename DoFHandlerType, template <int,int> class M_or_MC, int dim_>
    static inline
    void
    do_interpolate_boundary_values
    (const M_or_MC<DoFHandlerType::dimension, DoFHandlerType::space_dimension>                   &mapping,
     const DoFHandlerType                                                                        &dof,
     const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
     std::map<types::global_dof_index,number>                                                    &boundary_values,
     const ComponentMask                                                                         &component_mask,
     const dealii::internal::int2type<dim_>)
    {
      const unsigned int dim = DoFHandlerType::dimension;
      const unsigned int spacedim=DoFHandlerType::space_dimension;

      Assert (component_mask.represents_n_components(dof.get_fe().n_components()),
              ExcMessage ("The number of components in the mask has to be either "
                          "zero or equal to the number of components in the finite "
                          "element."));


      // if for whatever reason we were passed an empty map, return
      // immediately
      if (function_map.size() == 0)
        return;

      Assert (function_map.find(numbers::internal_face_boundary_id) == function_map.end(),
              ExcMessage("You cannot specify the special boundary indicator "
                         "for interior faces in your function map."));

      const unsigned int        n_components = DoFTools::n_components(dof);
      const bool                fe_is_system = (n_components != 1);

      for (typename std::map<types::boundary_id, const Function<spacedim,number>*>::const_iterator i=function_map.begin();
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
      std::vector<number>          dof_values_scalar;
      std::vector<Vector<number> > dof_values_system;
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

      typename DoFHandlerType::active_cell_iterator cell = dof.begin_active(),
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

              const typename DoFHandlerType::face_iterator face = cell->face(face_no);
              const types::boundary_id boundary_component = face->boundary_id();

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
                                                  Vector<number>(fe.n_components()));
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



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>                   &mapping,
   const DoFHandlerType                                                                        &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
   std::map<types::global_dof_index,number>                                                    &boundary_values,
   const ComponentMask                                                                         &component_mask_)
  {
    do_interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                    component_mask_,
                                    dealii::internal::int2type<DoFHandlerType::dimension>());
  }



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &mapping,
   const DoFHandlerType                                                      &dof,
   const types::boundary_id                                                   boundary_component,
   const Function<DoFHandlerType::space_dimension,number>                    &boundary_function,
   std::map<types::global_dof_index,number>                                  &boundary_values,
   const ComponentMask                                                       &component_mask)
  {
    std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                 component_mask);
  }


  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values
  (const hp::MappingCollection<dim,spacedim>                            &mapping,
   const hp::DoFHandler<dim,spacedim>                                   &dof,
   const std::map<types::boundary_id, const Function<spacedim,number>*> &function_map,
   std::map<types::global_dof_index,number>                             &boundary_values,
   const ComponentMask                                                  &component_mask_)
  {
    do_interpolate_boundary_values (mapping, dof, function_map, boundary_values,
                                    component_mask_,
                                    dealii::internal::int2type<dim>());
  }



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const DoFHandlerType                                   &dof,
   const types::boundary_id                                boundary_component,
   const Function<DoFHandlerType::space_dimension,number> &boundary_function,
   std::map<types::global_dof_index,number>               &boundary_values,
   const ComponentMask                                    &component_mask)
  {
    interpolate_boundary_values(StaticMappingQ1<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::mapping,
                                dof, boundary_component,
                                boundary_function, boundary_values, component_mask);
  }



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const DoFHandlerType                                                                        &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
   std::map<types::global_dof_index,number>                                                    &boundary_values,
   const ComponentMask                                                                         &component_mask)
  {
    interpolate_boundary_values
    (StaticMappingQ1<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::mapping,
     dof, function_map,
     boundary_values, component_mask);
  }




// ----------- interpolate_boundary_values for ConstraintMatrix --------------



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>                   &mapping,
   const DoFHandlerType                                                                        &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
   ConstraintMatrix                                                                            &constraints,
   const ComponentMask                                                                         &component_mask_)
  {
    std::map<types::global_dof_index,number> boundary_values;
    interpolate_boundary_values (mapping, dof, function_map,
                                 boundary_values, component_mask_);
    typename std::map<types::global_dof_index,number>::const_iterator boundary_value =
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



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &mapping,
   const DoFHandlerType                                                      &dof,
   const types::boundary_id                                                   boundary_component,
   const Function<DoFHandlerType::space_dimension,number>                    &boundary_function,
   ConstraintMatrix                                                          &constraints,
   const ComponentMask                                                       &component_mask)
  {
    std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> function_map;
    function_map[boundary_component] = &boundary_function;
    interpolate_boundary_values (mapping, dof, function_map, constraints,
                                 component_mask);
  }



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const DoFHandlerType                                   &dof,
   const types::boundary_id                                boundary_component,
   const Function<DoFHandlerType::space_dimension,number> &boundary_function,
   ConstraintMatrix                                       &constraints,
   const ComponentMask                                    &component_mask)
  {
    interpolate_boundary_values
    (StaticMappingQ1<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::mapping,
     dof, boundary_component,
     boundary_function, constraints, component_mask);
  }



  template <typename DoFHandlerType, typename number>
  void
  interpolate_boundary_values
  (const DoFHandlerType                                                                        &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &function_map,
   ConstraintMatrix                                                                            &constraints,
   const ComponentMask                                                                         &component_mask)
  {
    interpolate_boundary_values
    (StaticMappingQ1<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::mapping,
     dof, function_map,
     constraints, component_mask);
  }




// -------- implementation for project_boundary_values with std::map --------


  namespace
  {
    // keep the first argument non-reference since we use it
    // with 1e-8 * number
    template <typename number1, typename number2>
    bool real_part_bigger_than(const number1 a,
                               const number2 &b)
    {
      return a > b;
    }

    template <typename number1, typename number2>
    bool real_part_bigger_than(const std::complex<number1> a,
                               const std::complex<number2> &b)
    {
      Assert(std::abs(a.imag()) <= 1e-15*std::abs(a) , ExcInternalError());
      Assert(std::abs(b.imag()) <= 1e-15*std::abs(b) , ExcInternalError());
      return a.real() > b.real();
    }

    // this function is needed to get an idea where
    // rhs.norm_sqr()  is too small for a given type.
    template <typename number>
    number min_number(const number &/*dummy*/)
    {
      return std::numeric_limits<number>::min();
    }

    // Sine rhs.norm_sqr() is non-negative real, in complex case we
    // take the numeric limits of the underlying type used in std::complex<>.
    template <typename number>
    number min_number(const std::complex<number> &/*dummy*/)
    {
      return std::numeric_limits<number>::min();
    }


    template <typename number>
    void invert_mass_matrix(const SparseMatrix<number> &mass_matrix,
                            const FilteredMatrix<Vector<number> > &filtered_mass_matrix,
                            FilteredMatrix<Vector<number> > &filtered_preconditioner,
                            const Vector<number> &rhs,
                            Vector<number> &boundary_projection)
    {
      // Allow for a maximum of 5*n steps to reduce the residual by 10^-12. n
      // steps may not be sufficient, since roundoff errors may accumulate for
      // badly conditioned matrices
      ReductionControl        control(5*rhs.size(), 0., 1.e-12, false, false);
      GrowingVectorMemory<Vector<number> > memory;
      SolverCG<Vector<number> >            cg(control,memory);

      PreconditionSSOR<SparseMatrix<number> > prec;
      prec.initialize(mass_matrix, 1.2);
      filtered_preconditioner.initialize(prec, true);
      // solve
      cg.solve (filtered_mass_matrix, boundary_projection, rhs, filtered_preconditioner);
      filtered_preconditioner.apply_constraints(boundary_projection, true);
      filtered_preconditioner.clear();
    }

    template <typename number>
    void invert_mass_matrix(const SparseMatrix<std::complex<number> > &/*mass_matrix*/,
                            const FilteredMatrix<Vector<std::complex<number> > > &/*filtered_mass_matrix*/,
                            FilteredMatrix<Vector<std::complex<number> > > &/*filtered_preconditioner*/,
                            const Vector<std::complex<number> > &/*rhs*/,
                            Vector<std::complex<number> > &/*boundary_projection*/)
    {
      Assert(false, ExcNotImplemented());
    }


    template <int dim, int spacedim, template <int, int> class DoFHandlerType,
              template <int,int> class M_or_MC, template <int> class Q_or_QC, typename number>
    void
    do_project_boundary_values (const M_or_MC<dim, spacedim>                                         &mapping,
                                const DoFHandlerType<dim, spacedim>                                  &dof,
                                const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                                const Q_or_QC<dim-1>                                                 &q,
                                std::map<types::global_dof_index,number>                             &boundary_values,
                                std::vector<unsigned int>                                             component_mapping)
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
      for (typename std::map<types::boundary_id, const Function<spacedim,number>*>::const_iterator i=boundary_functions.begin();
           i!=boundary_functions.end(); ++i)
        selected_boundary_components.insert (i->first);

      DoFTools::map_dof_to_boundary_indices (dof, selected_boundary_components,
                                             dof_to_boundary_mapping);

      // Done if no degrees of freedom on the boundary
      if (dof.n_boundary_dofs (boundary_functions) == 0)
        return;

      // set up sparsity structure
      DynamicSparsityPattern dsp(dof.n_boundary_dofs(boundary_functions),
                                 dof.n_boundary_dofs(boundary_functions));
      DoFTools::make_boundary_sparsity_pattern (dof,
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
      if (dim>=3)
        {
#ifdef DEBUG
          // Assert that there are no hanging nodes at the boundary
          int level = -1;
          for (typename DoFHandlerType<dim,spacedim>::active_cell_iterator cell = dof.begin_active();
               cell != dof.end(); ++cell)
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (cell->at_boundary(f))
                  {
                    if (level == -1)
                      level = cell->level();
                    else
                      {
                        Assert (level == cell->level(),
                                ExcMessage("The mesh you use in projecting boundary values "
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


      MatrixCreator::create_boundary_mass_matrix (mapping, dof, q,
                                                  mass_matrix, boundary_functions,
                                                  rhs, dof_to_boundary_mapping, (const Function<spacedim,number> *) 0,
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

      FilteredMatrix<Vector<number> > filtered_mass_matrix(mass_matrix, true);
      FilteredMatrix<Vector<number> > filtered_precondition;
      std::vector<bool> excluded_dofs(mass_matrix.m(), false);

      // we assemble mass matrix with unit weight,
      // thus it will be real-valued irrespectively of the underlying algebra
      // with positive elements on diagonal.
      // Thus in order to extend this filtering to complex-algebra simply take
      // the real-part of element.
      number max_element = 0.;
      for (unsigned int i=0; i<mass_matrix.m(); ++i)
        if (real_part_bigger_than(mass_matrix.diag_element(i),max_element))
          max_element = mass_matrix.diag_element(i);

      for (unsigned int i=0; i<mass_matrix.m(); ++i)
        if (real_part_bigger_than(1.e-8 * max_element,mass_matrix.diag_element(i)))
          {
            filtered_mass_matrix.add_constraint(i, 0.);
            filtered_precondition.add_constraint(i, 0.);
            mass_matrix.diag_element(i) = 1.;
            excluded_dofs[i] = true;
          }

      Vector<number> boundary_projection (rhs.size());

      // cannot reduce residual in a useful way if we are close to the square
      // root of the minimal double value
      if (rhs.norm_sqr() < 1e28 * min_number(number()))
        boundary_projection = 0;
      else
        {
          invert_mass_matrix(mass_matrix,
                             filtered_mass_matrix,
                             filtered_precondition,
                             rhs,
                             boundary_projection);
        }
      // fill in boundary values
      for (unsigned int i=0; i<dof_to_boundary_mapping.size(); ++i)
        if (dof_to_boundary_mapping[i] != DoFHandler<dim,spacedim>::invalid_dof_index
            && ! excluded_dofs[dof_to_boundary_mapping[i]])
          {
            AssertIsFinite(boundary_projection(dof_to_boundary_mapping[i]));

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

  template <int dim, int spacedim, typename number>
  void
  project_boundary_values (const Mapping<dim, spacedim>                                         &mapping,
                           const DoFHandler<dim, spacedim>                                      &dof,
                           const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                           const Quadrature<dim-1>                                              &q,
                           std::map<types::global_dof_index,number>                             &boundary_values,
                           std::vector<unsigned int>                                             component_mapping)
  {
    do_project_boundary_values(mapping, dof, boundary_functions, q,
                               boundary_values, component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values (const DoFHandler<dim,spacedim>                                       &dof,
                           const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                           const Quadrature<dim-1>                                              &q,
                           std::map<types::global_dof_index,number>                             &boundary_values,
                           std::vector<unsigned int>                                             component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
                            boundary_values, component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void project_boundary_values (const hp::MappingCollection<dim, spacedim>                           &mapping,
                                const hp::DoFHandler<dim,spacedim>                                   &dof,
                                const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                                const hp::QCollection<dim-1>                                         &q,
                                std::map<types::global_dof_index,number>                             &boundary_values,
                                std::vector<unsigned int>                                             component_mapping)
  {
    do_project_boundary_values (mapping, dof,
                                boundary_functions,
                                q, boundary_values,
                                component_mapping);
  }



  template <int dim, int spacedim, typename number>
  void project_boundary_values (const hp::DoFHandler<dim,spacedim>                                   &dof,
                                const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_function,
                                const hp::QCollection<dim-1>                                         &q,
                                std::map<types::global_dof_index,number>                             &boundary_values,
                                std::vector<unsigned int>                                             component_mapping)
  {
    project_boundary_values (hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof,
                             boundary_function,
                             q, boundary_values,
                             component_mapping);
  }


  // ----- implementation for project_boundary_values with ConstraintMatrix -----



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values (const Mapping<dim, spacedim>                                         &mapping,
                           const DoFHandler<dim,spacedim>                                       &dof,
                           const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                           const Quadrature<dim-1>                                              &q,
                           ConstraintMatrix                                                     &constraints,
                           std::vector<unsigned int>                                             component_mapping)
  {
    std::map<types::global_dof_index,number> boundary_values;
    project_boundary_values (mapping, dof, boundary_functions, q,
                             boundary_values, component_mapping);
    typename std::map<types::global_dof_index,number>::const_iterator boundary_value =
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



  template <int dim, int spacedim, typename number>
  void
  project_boundary_values (const DoFHandler<dim,spacedim>                                       &dof,
                           const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_functions,
                           const Quadrature<dim-1>                                              &q,
                           ConstraintMatrix                                                     &constraints,
                           std::vector<unsigned int>                                             component_mapping)
  {
    project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
                            constraints, component_mapping);
  }




  namespace internal
  {
    /**
     * A structure that stores the dim DoF indices that correspond to a
     * vector-valued quantity at a single support point.
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
     * Add the constraint $\vec n \cdot \vec u = inhom$ to the list of
     * constraints.
     *
     * Here, $\vec u$ is represented by the set of given DoF indices, and
     * $\vec n$ by the vector specified as the second argument.
     *
     * The function does not add constraints if a degree of freedom is already
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
     * Add the constraint $(\vec u-\vec u_\Gamma) \| \vec t$ to the list of
     * constraints. In 2d, this is a single constraint, in 3d these are two
     * constraints
     *
     * Here, $\vec u$ is represented by the set of given DoF indices, and
     * $\vec t$ by the vector specified as the second argument.
     *
     * The function does not add constraints if a degree of freedom is already
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
     * Given a vector, compute a set of dim-1 vectors that are orthogonal to
     * the first one and mutually orthonormal as well.
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
          orthonormals[0] = cross_product_3d (vector, tmp);
          orthonormals[1] = cross_product_3d (vector, orthonormals[0]);

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

      std::vector<Tensor<1,dim> > tangentials (fe_values.n_quadrature_points);
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
          tangentials[q_point] /= tangentials[q_point].norm();

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
                const double tangential_solution_component
                  = (values[q_point] (first_vector_component) * tangentials[q_point][0]
                     + values[q_point] (first_vector_component + 1) * tangentials[q_point][1]
                     + values[q_point] (first_vector_component + 2) * tangentials[q_point][2]);
                dof_values[i]
                += (fe_values.JxW (q_point)
                    * tangential_solution_component
                    * (fe_values[vec].value (fe.face_to_cell_index (i, face), q_point) *
                       tangentials[q_point])
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
          std::vector<Tensor<1,dim> > tangentials (fe_values.n_quadrature_points);

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
              tangentials[q_point] /= tangentials[q_point].norm();

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
                          * tangentials[q_point][0]
                          + values[q_point] (first_vector_component + 1)
                          * tangentials[q_point][1])
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
    template<typename cell_iterator>
    void
    compute_edge_projection_l2 (const cell_iterator &cell,
                                const unsigned int face,
                                const unsigned int line,
                                hp::FEValues<3> &hp_fe_values,
                                const Function<3> &boundary_function,
                                const unsigned int first_vector_component,
                                std::vector<double> &dof_values,
                                std::vector<bool> &dofs_processed)
    {
      // This function computes the L2-projection of the given
      // boundary function on 3D edges and returns the constraints
      // associated with the edge functions for the given cell.
      //
      // In the context of this function, by associated DoFs we mean:
      // the DoFs corresponding to the group of components making up the vector
      // with first component first_vector_component (length dim).
      const unsigned int dim = 3;
      const FiniteElement<dim> &fe = cell->get_fe ();

      // reinit for this cell, face and line.
      hp_fe_values.reinit
      (cell,
       (cell->active_fe_index () * GeometryInfo<dim>::faces_per_cell + face)
       * GeometryInfo<dim>::lines_per_face + line);

      // Initialize the required objects.
      const FEValues<dim> &
      fe_values = hp_fe_values.get_present_fe_values ();
      // Store degree as fe.degree-1
      // For nedelec elements FE_Nedelec<dim> (0) returns fe.degree = 1.
      const unsigned int degree = fe.degree - 1;

      const std::vector<Point<dim> > &
      quadrature_points = fe_values.get_quadrature_points ();
      std::vector<Vector<double> > values (fe_values.n_quadrature_points,
                                           Vector<double> (fe.n_components ()));

      // Get boundary function values
      // at quadrature points.
      boundary_function.vector_value_list (quadrature_points, values);

      // Find the group of vector components (dim of them,
      // starting at first_vector_component) are within an FESystem.
      //
      // If not using FESystem then must be using FE_Nedelec,
      // which has one base element and one copy of it (with 3 components).
      std::pair<unsigned int, unsigned int> base_indices (0, 0);
      if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) != 0)
        {
          unsigned int fe_index = 0;
          unsigned int fe_index_old = 0;
          unsigned int i = 0;

          // Find base element:
          // base_indices.first
          //
          // Then select which copy of that base element
          // [ each copy is of length fe.base_element(base_indices.first).n_components() ]
          // corresponds to first_vector_component:
          // base_index.second
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

      // Find DoFs we want to constrain:
      // There are fe.dofs_per_line DoFs associated with the
      // given line on the given face on the given cell.
      //
      // Want to know which of these DoFs (there are degree+1 of interest)
      // are associated with the components given by first_vector_component.
      // Then we can make a map from the associated line DoFs to the face DoFs.
      //
      // For a single FE_Nedelec<3> element this is simple:
      //    We know the ordering of local DoFs goes
      //    lines -> faces -> cells
      //
      // For a set of FESystem<3> elements we need to pick out the matching base element and
      // the index within this ordering.
      //
      // We call the map associated_edge_dof_to_face_dof
      std::vector<unsigned int> associated_edge_dof_to_face_dof (degree + 1);

      // Lowest DoF in the base element allowed for this edge:
      const unsigned int lower_bound =
        fe.base_element(base_indices.first).face_to_cell_index(line * (degree + 1), face);
      // Highest DoF in the base element allowed for this edge:
      const unsigned int upper_bound =
        fe.base_element(base_indices.first).face_to_cell_index((line + 1) * (degree + 1) - 1, face);

      unsigned int associated_edge_dof_index = 0;
      //       for (unsigned int face_idx = 0; face_idx < fe.dofs_per_face; ++face_idx)
      for (unsigned int line_idx = 0; line_idx < fe.dofs_per_line; ++line_idx)
        {
          // Assuming DoFs on a face are numbered in order by lines then faces.
          // i.e. line 0 has degree+1 dofs numbered 0,..,degree
          //      line 1 has degree+1 dofs numbered (degree+1),..,2*(degree+1) and so on.
          const unsigned int face_idx = line*fe.dofs_per_line + line_idx;
          // Note, assuming that the edge orientations are "standard"
          //       i.e. cell->line_orientation(line) = true.
          const unsigned int cell_idx = fe.face_to_cell_index (face_idx, face);

          // Check this cell_idx belongs to the correct base_element, component and line:
          if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
               && (fe.system_to_base_index (cell_idx).first == base_indices)
               && (lower_bound <= fe.system_to_base_index (cell_idx).second)
               && (fe.system_to_base_index (cell_idx).second <= upper_bound))
              || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                  && (line * (degree + 1) <= face_idx)
                  && (face_idx <= (line + 1) * (degree + 1) - 1)))
            {
              associated_edge_dof_to_face_dof[associated_edge_dof_index] = face_idx;
              ++associated_edge_dof_index;
            }
        }
      // Sanity check:
      const unsigned int associated_edge_dofs = associated_edge_dof_index;
      Assert (associated_edge_dofs == degree + 1,
              ExcMessage ("Error: Unexpected number of 3D edge DoFs"));

      // Matrix and RHS vectors to store linear system:
      // We have (degree+1) basis functions for an edge
      FullMatrix<double> edge_matrix (degree + 1,degree + 1);
      FullMatrix<double> edge_matrix_inv (degree + 1,degree + 1);
      Vector<double> edge_rhs (degree + 1);
      Vector<double> edge_solution (degree + 1);

      const FEValuesExtractors::Vector vec (first_vector_component);

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

      const double tol = 0.5 * cell->face (face)->line (line)->diameter () / fe.degree;
      const std::vector<Point<dim> > &
      reference_quadrature_points = fe_values.get_quadrature ().get_points ();

      // Project the boundary function onto the shape functions for this edge
      // and set up a linear system of equations to get the values for the DoFs
      // associated with this edge.
      for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; ++q_point)
        {
          // Compute the tangential
          // of the edge at
          // the quadrature point.
          Point<dim> shifted_reference_point_1 = reference_quadrature_points[q_point];
          Point<dim> shifted_reference_point_2 = reference_quadrature_points[q_point];

          shifted_reference_point_1 (edge_coordinate_direction[face][line]) += tol;
          shifted_reference_point_2 (edge_coordinate_direction[face][line]) -= tol;
          Tensor<1,dim> tangential
            = (0.5 *
               (fe_values.get_mapping ()
                .transform_unit_to_real_cell (cell,
                                              shifted_reference_point_1)
                -
                fe_values.get_mapping ()
                .transform_unit_to_real_cell (cell,
                                              shifted_reference_point_2))
               / tol);
          tangential
          /= tangential.norm ();

          // Compute the entires of the linear system
          // Note the system is symmetric so we could only compute the lower/upper triangle.
          //
          // The matrix entries are
          // \int_{edge} (tangential*edge_shape_function_i)*(tangential*edge_shape_function_j) dS
          //
          // The RHS entries are:
          // \int_{edge} (tangential*boundary_value)*(tangential*edge_shape_function_i) dS.
          for (unsigned int j = 0; j < associated_edge_dofs; ++j)
            {
              const unsigned int j_face_idx = associated_edge_dof_to_face_dof[j];
              const unsigned int j_cell_idx = fe.face_to_cell_index (j_face_idx, face);
              for (unsigned int i = 0; i < associated_edge_dofs; ++i)
                {
                  const unsigned int i_face_idx = associated_edge_dof_to_face_dof[i];
                  const unsigned int i_cell_idx = fe.face_to_cell_index (i_face_idx, face);

                  edge_matrix(i,j)
                  += fe_values.JxW (q_point)
                     * (fe_values[vec].value (i_cell_idx, q_point) * tangential)
                     * (fe_values[vec].value (j_cell_idx, q_point) * tangential);
                }
              // Compute the RHS entries:
              edge_rhs(j) += fe_values.JxW (q_point)
                             * (values[q_point] (first_vector_component) * tangential [0]
                                + values[q_point] (first_vector_component + 1) * tangential [1]
                                + values[q_point] (first_vector_component + 2) * tangential [2])
                             * (fe_values[vec].value (j_cell_idx, q_point) * tangential);
            }
        }

      // Invert linear system
      edge_matrix_inv.invert(edge_matrix);
      edge_matrix_inv.vmult(edge_solution,edge_rhs);

      // Store computed DoFs
      for (unsigned int i = 0; i < associated_edge_dofs; ++i)
        {
          dof_values[associated_edge_dof_to_face_dof[i]] = edge_solution(i);
          dofs_processed[associated_edge_dof_to_face_dof[i]] = true;
        }
    }


    template<int dim, typename cell_iterator>
    void
    compute_edge_projection_l2 (const cell_iterator &,
                                const unsigned int,
                                const unsigned int,
                                hp::FEValues<dim> &,
                                const Function<dim> &,
                                const unsigned int,
                                std::vector<double> &,
                                std::vector<bool> &)
    {
      // dummy implementation of above function
      // for all other dimensions
      Assert (false, ExcInternalError ());
    }

    template<int dim, typename cell_iterator>
    void
    compute_face_projection_curl_conforming_l2 (const cell_iterator &cell,
                                                const unsigned int face,
                                                hp::FEFaceValues<dim> &hp_fe_face_values,
                                                const Function<dim> &boundary_function,
                                                const unsigned int first_vector_component,
                                                std::vector<double> &dof_values,
                                                std::vector<bool> &dofs_processed)
    {
      // This function computes the L2-projection of the boundary
      // function on the interior of faces only. In 3D, this should only be
      // called after first calling compute_edge_projection_l2, as it relies on
      // edge constraints which are found.

      // In the context of this function, by associated DoFs we mean:
      // the DoFs corresponding to the group of components making up the vector
      // with first component first_vector_component (with total components dim).

      // Copy to the standard FEFaceValues object:
      hp_fe_face_values.reinit (cell, face);
      const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();

      // Initialize the required objects.
      const FiniteElement<dim> &fe = cell->get_fe ();
      const std::vector<Point<dim> > &
      quadrature_points = fe_face_values.get_quadrature_points ();
      const unsigned int degree = fe.degree - 1;

      std::vector<Vector<double> >
      values (fe_face_values.n_quadrature_points, Vector<double> (fe.n_components ()));

      // Get boundary function values at quadrature points.
      boundary_function.vector_value_list (quadrature_points, values);

      // Find where the group of vector components (dim of them,
      // starting at first_vector_component) are within an FESystem.
      //
      // If not using FESystem then must be using FE_Nedelec,
      // which has one base element and one copy of it (with 3 components).
      std::pair<unsigned int, unsigned int> base_indices (0, 0);
      if (dynamic_cast<const FESystem<dim>*> (&cell->get_fe ()) != 0)
        {
          unsigned int fe_index = 0;
          unsigned int fe_index_old = 0;
          unsigned int i = 0;

          // Find base element:
          // base_indices.first
          //
          // Then select which copy of that base element
          // [ each copy is of length fe.base_element(base_indices.first).n_components() ]
          // corresponds to first_vector_component:
          // base_index.second
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

      switch (dim)
        {
        case 2:
          // NOTE: This is very similar to compute_edge_projection as used in 3D,
          //       and contains a lot of overlap with that function.
        {

          // Find the DoFs we want to constrain. There are degree+1 in total.
          // Create a map from these to the face index
          // Note:
          //    - for a single FE_Nedelec<2> element this is
          //      simply 0 to fe.dofs_per_face
          //    - for FESystem<2> this just requires matching the
          //      base element, fe.system_to_base_index.first.first
          //      and the copy of the base element we're interested
          //      in, fe.system_to_base_index.first.second
          std::vector<unsigned int> associated_edge_dof_to_face_dof (degree + 1);

          unsigned int associated_edge_dof_index = 0;
          for (unsigned int face_idx = 0; face_idx < fe.dofs_per_face; ++face_idx)
            {
              const unsigned int cell_idx = fe.face_to_cell_index (face_idx, face);
              if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                   && (fe.system_to_base_index (cell_idx).first == base_indices))
                  || (dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0))
                {
                  associated_edge_dof_to_face_dof[associated_edge_dof_index] = face_idx;
                  ++associated_edge_dof_index;
                }
            }
          // Sanity check:
          const unsigned int associated_edge_dofs = associated_edge_dof_index;
          Assert (associated_edge_dofs == degree + 1,
                  ExcMessage ("Error: Unexpected number of 2D edge DoFs"));

          // Matrix and RHS vectors to store:
          // We have (degree+1) edge basis functions
          FullMatrix<double> edge_matrix (degree + 1,degree + 1);
          FullMatrix<double> edge_matrix_inv (degree + 1,degree + 1);
          Vector<double> edge_rhs (degree + 1);
          Vector<double> edge_solution (degree + 1);

          const FEValuesExtractors::Vector vec (first_vector_component);

          // Project the boundary function onto the shape functions for this edge
          // and set up a linear system of equations to get the values for the DoFs
          // associated with this edge.
          for (unsigned int q_point = 0;
               q_point < fe_face_values.n_quadrature_points; ++q_point)
            {
              // Compute the entires of the linear system
              // Note the system is symmetric so we could only compute the lower/upper triangle.
              //
              // The matrix entries are
              // \int_{edge} (tangential * edge_shape_function_i) * (tangential * edge_shape_function_j) dS
              //
              // The RHS entries are:
              // \int_{edge} (tangential* boundary_value) * (tangential * edge_shape_function_i) dS.
              //
              // In 2D, tangential*vector is equivalent to cross_product_3d(normal, vector), so we use this instead.
              // This avoids possible issues with the computation of the tangent.

              // Store the normal at this quad point:
              Tensor<1,dim> normal_at_q_point = fe_face_values.normal_vector(q_point);
              for (unsigned int j = 0; j < associated_edge_dofs; ++j)
                {
                  const unsigned int j_face_idx = associated_edge_dof_to_face_dof[j];
                  const unsigned int j_cell_idx = fe.face_to_cell_index (j_face_idx, face);

                  Tensor<1,dim> phi_j = fe_face_values[vec].value (j_cell_idx, q_point);
                  for (unsigned int i = 0; i < associated_edge_dofs; ++i)
                    {
                      const unsigned int i_face_idx = associated_edge_dof_to_face_dof[i];
                      const unsigned int i_cell_idx = fe.face_to_cell_index (i_face_idx, face);

                      Tensor<1,dim> phi_i = fe_face_values[vec].value (i_cell_idx, q_point);

                      // Using n cross phi
                      edge_matrix(i,j)
                      += fe_face_values.JxW (q_point)
                         * ((phi_i[1]*normal_at_q_point[0] - phi_i[0]*normal_at_q_point[1])
                            * (phi_j[1]*normal_at_q_point[0] - phi_j[0]*normal_at_q_point[1]));
                    }
                  // Using n cross phi
                  edge_rhs(j)
                  += fe_face_values.JxW (q_point)
                     * ((values[q_point] (first_vector_component+1) * normal_at_q_point[0]
                         - values[q_point] (first_vector_component) * normal_at_q_point[1])
                        * (phi_j[1]*normal_at_q_point[0] - phi_j[0]*normal_at_q_point[1]));
                }
            }

          // Invert linear system
          edge_matrix_inv.invert(edge_matrix);
          edge_matrix_inv.vmult(edge_solution,edge_rhs);

          // Store computed DoFs
          for (unsigned int associated_edge_dof_index = 0; associated_edge_dof_index < associated_edge_dofs; ++associated_edge_dof_index)
            {
              dof_values[associated_edge_dof_to_face_dof[associated_edge_dof_index]] = edge_solution (associated_edge_dof_index);
              dofs_processed[associated_edge_dof_to_face_dof[associated_edge_dof_index]] = true;
            }
          break;
        }

        case 3:
        {
          const FEValuesExtractors::Vector vec (first_vector_component);

          // First group DoFs associated with edges which we already know.
          // Sort these into groups of dofs (0 -> degree+1 of them) by each edge.
          // This will help when computing the residual for the face projections.
          //
          // This matches with the search done in compute_edge_projection.
          const unsigned int lines_per_face = GeometryInfo<dim>::lines_per_face;
          std::vector<std::vector<unsigned int> >
          associated_edge_dof_to_face_dof(lines_per_face, std::vector<unsigned int> (degree + 1));
          std::vector<unsigned int> associated_edge_dofs (lines_per_face);

          for (unsigned int line = 0; line < lines_per_face; ++line)
            {
              // Lowest DoF in the base element allowed for this edge:
              const unsigned int lower_bound =
                fe.base_element(base_indices.first).face_to_cell_index(line * (degree + 1), face);
              // Highest DoF in the base element allowed for this edge:
              const unsigned int upper_bound =
                fe.base_element(base_indices.first).face_to_cell_index((line + 1) * (degree + 1) - 1, face);
              unsigned int associated_edge_dof_index = 0;
              for (unsigned int line_idx = 0; line_idx < fe.dofs_per_line; ++line_idx)
                {
                  const unsigned int face_idx = line*fe.dofs_per_line + line_idx;
                  const unsigned int cell_idx = fe.face_to_cell_index(face_idx, face);
                  // Check this cell_idx belongs to the correct base_element, component and line:
                  if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                       && (fe.system_to_base_index (cell_idx).first == base_indices)
                       && (lower_bound <= fe.system_to_base_index (cell_idx).second)
                       && (fe.system_to_base_index (cell_idx).second <= upper_bound))
                      || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)
                          && (line * (degree + 1) <= face_idx)
                          && (face_idx <= (line + 1) * (degree + 1) - 1)))
                    {
                      associated_edge_dof_to_face_dof[line][associated_edge_dof_index] = face_idx;
                      ++associated_edge_dof_index;
                    }
                }
              // Sanity check:
              associated_edge_dofs[line] = associated_edge_dof_index;
              Assert (associated_edge_dofs[line] == degree + 1,
                      ExcMessage ("Error: Unexpected number of 3D edge DoFs"));
            }

          // Next find the face DoFs associated with the vector components
          // we're interested in. There are 2*degree*(degree+1) DoFs associated
          // with each face (not including edges!).
          //
          // Create a map mapping from the consecutively numbered associated_dofs
          // to the face DoF (which can be transferred to a local cell index).
          //
          // For FE_Nedelec<3> we just need to have a face numbering greater than
          // the number of edge DoFs (=lines_per_face*(degree+1).
          //
          // For FE_System<3> we need to base the base_indices (base element and
          // copy within that base element) and ensure we're above the number of
          // edge DoFs within that base element.
          std::vector<unsigned int> associated_face_dof_to_face_dof (2*degree*(degree + 1));

          // Skip the edge DoFs, so we start at lines_per_face*(fe.dofs_per_line).
          unsigned int associated_face_dof_index = 0;
          for (unsigned int face_idx = lines_per_face*(fe.dofs_per_line); face_idx < fe.dofs_per_face; ++face_idx)
            {
              const unsigned int cell_idx = fe.face_to_cell_index (face_idx, face);
              if (((dynamic_cast<const FESystem<dim>*> (&fe) != 0)
                   && (fe.system_to_base_index (cell_idx).first == base_indices))
                  || ((dynamic_cast<const FE_Nedelec<dim>*> (&fe) != 0)))
                {
                  associated_face_dof_to_face_dof[associated_face_dof_index] = face_idx;
                  ++associated_face_dof_index;
                }
            }
          // Sanity check:
          const unsigned int associated_face_dofs = associated_face_dof_index;
          Assert (associated_face_dofs == 2*degree*(degree + 1),
                  ExcMessage ("Error: Unexpected number of 3D face DoFs"));

          // Storage for the linear system.
          // There are 2*degree*(degree+1) DoFs associated with a face in 3D.
          // Note this doesn't include the DoFs associated with edges on that face.
          FullMatrix<double> face_matrix (2*degree*(degree + 1));
          FullMatrix<double> face_matrix_inv (2*degree*(degree + 1));
          Vector<double> face_rhs (2*degree*(degree + 1));
          Vector<double> face_solution (2*degree*(degree + 1));

          // Project the boundary function onto the shape functions for this face
          // and set up a linear system of equations to get the values for the DoFs
          // associated with this face. We also must include the residuals from the
          // shape funcations associated with edges.
          Tensor<1, dim> tmp;
          Tensor<1, dim> cross_product_i;
          Tensor<1, dim> cross_product_j;
          Tensor<1, dim> cross_product_rhs;

          // Loop to construct face linear system.
          for (unsigned int q_point = 0;
               q_point < fe_face_values.n_quadrature_points; ++q_point)
            {
              // First calculate the residual from the edge functions
              // store the result in tmp.
              //
              // Edge_residual =
              //        boundary_value - (
              //            \sum_(edges on face)
              //                 \sum_(DoFs on edge) edge_dof_value*edge_shape_function
              //                   )
              for (unsigned int d = 0; d < dim; ++d)
                {
                  tmp[d]=0.0;
                }
              for (unsigned int line = 0; line < lines_per_face; ++line)
                {
                  for (unsigned int associated_edge_dof = 0; associated_edge_dof < associated_edge_dofs[line]; ++associated_edge_dof)
                    {
                      const unsigned int face_idx = associated_edge_dof_to_face_dof[line][associated_edge_dof];
                      const unsigned int cell_idx = fe.face_to_cell_index(face_idx, face);
                      tmp -= dof_values[face_idx]*fe_face_values[vec].value(cell_idx, q_point);
                    }
                }

              for (unsigned int d=0; d<dim; ++d)
                {
                  tmp[d] += values[q_point] (first_vector_component + d);
                }

              // Tensor of normal vector on the face at q_point;
              const Tensor<1,dim> normal_vector = fe_face_values.normal_vector(q_point);

              // Now compute the linear system:
              // On a face:
              // The matrix entries are:
              // \int_{face} (n x face_shape_function_i) \cdot ( n x face_shape_function_j) dS
              //
              // The RHS entries are:
              // \int_{face} (n x (Edge_residual) \cdot (n x face_shape_function_i) dS

              for (unsigned int j = 0; j < associated_face_dofs; ++j)
                {
                  const unsigned int j_face_idx = associated_face_dof_to_face_dof[j];
                  const unsigned int cell_j = fe.face_to_cell_index (j_face_idx, face);

                  cross_product_j =
                    cross_product_3d(normal_vector,
                                     fe_face_values[vec].value(cell_j, q_point));

                  for (unsigned int i = 0; i < associated_face_dofs; ++i)
                    {
                      const unsigned int i_face_idx = associated_face_dof_to_face_dof[i];
                      const unsigned int cell_i = fe.face_to_cell_index (i_face_idx, face);
                      cross_product_i =
                        cross_product_3d(normal_vector,
                                         fe_face_values[vec].value(cell_i, q_point));

                      face_matrix(i, j) += fe_face_values.JxW(q_point) *
                                           cross_product_i * cross_product_j;
                    }
                  // compute rhs
                  cross_product_rhs = cross_product_3d(normal_vector, tmp);
                  face_rhs(j) += fe_face_values.JxW(q_point) *
                                 cross_product_rhs * cross_product_j;
                }
            }

          // Solve lienar system:
          face_matrix_inv.invert(face_matrix);
          face_matrix_inv.vmult(face_solution, face_rhs);


          // Store computed DoFs:
          for (unsigned int associated_face_dof = 0; associated_face_dof < associated_face_dofs; ++associated_face_dof)
            {
              dof_values[associated_face_dof_to_face_dof[associated_face_dof]] = face_solution(associated_face_dof);
              dofs_processed[associated_face_dof_to_face_dof[associated_face_dof]] = true;
            }
          break;
        }
        default:
          Assert (false, ExcNotImplemented ());
        }
    }

    template <int dim, typename DoFHandlerType>
    void
    compute_project_boundary_values_curl_conforming_l2
    (const DoFHandlerType                  &dof_handler,
     const unsigned int                     first_vector_component,
     const Function<dim>                   &boundary_function,
     const types::boundary_id               boundary_component,
     ConstraintMatrix                      &constraints,
     const hp::MappingCollection<dim, dim> &mapping_collection)
    {
      // L2-projection based interpolation formed in one (in 2D) or two (in 3D) steps.
      //
      // In 2D we only need to constrain edge DoFs.
      //
      // In 3D we need to constrain both edge and face DoFs. This is done in two parts.
      //
      // For edges, since the face shape functions are zero here ("bubble functions"),
      // we project the tangential component of the boundary function and compute
      // the L2-projection. This returns the values for the DoFs associated with each edge shape
      // function. In 3D, this is computed by internals::compute_edge_projection_l2, in 2D,
      // it is handled by compute_face_projection_curl_conforming_l2.
      //
      // For faces we compute the residual of the boundary function which is satisfied
      // by the edge shape functions alone. Which can then be used to calculate the
      // remaining face DoF values via a projection which leads to a linear system to
      // solve. This is handled by compute_face_projection_curl_conforming_l2
      //
      // For details see (for example) section 4.2:
      // Electromagnetic scattering simulation using an H (curl) conforming hp finite element
      // method in three dimensions, PD Ledger, K Morgan, O Hassan,
      // Int. J.  Num. Meth. Fluids, Volume 53, Issue 8, pages 1267–1296, 20 March 2007:
      // http://onlinelibrary.wiley.com/doi/10.1002/fld.1223/abstract

      // Create hp FEcollection, dof_handler can be either hp or standard type.
      // From here on we can treat it like a hp-namespace object.
      const hp::FECollection<dim> fe_collection (dof_handler.get_fe ());

      // Create face quadrature collection
      hp::QCollection<dim - 1> face_quadrature_collection;
      for (unsigned int i = 0; i < fe_collection.size (); ++i)
        {
          const QGauss<dim - 1>  reference_face_quadrature (2 * fe_collection[i].degree + 1);
          face_quadrature_collection.push_back(reference_face_quadrature);
        }

      hp::FEFaceValues<dim> fe_face_values (mapping_collection, fe_collection,
                                            face_quadrature_collection,
                                            update_values |
                                            update_quadrature_points |
                                            update_normal_vectors |
                                            update_JxW_values);

      // Storage for dof values found and whether they have been processed:
      std::vector<bool> dofs_processed;
      std::vector<double> dof_values;
      std::vector<types::global_dof_index> face_dof_indices;
      typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active ();

      switch (dim)
        {
        case 2:
        {
          for (; cell != dof_handler.end (); ++cell)
            {
              if (cell->at_boundary () && cell->is_locally_owned ())
                {
                  for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                    {
                      if (cell->face (face)->boundary_id () == boundary_component)
                        {
                          // If the FE is an FE_Nothing object there is no work to do
                          if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                            {
                              return;
                            }

                          // This is only implemented for FE_Nedelec elements.
                          // If the FE is a FESystem we cannot check this.
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

                          // Compute the projection of the boundary function on the edge.
                          // In 2D this is all that's required.
                          compute_face_projection_curl_conforming_l2 (cell, face, fe_face_values,
                                                                      boundary_function,
                                                                      first_vector_component,
                                                                      dof_values, dofs_processed);

                          // store the local->global map:
                          face_dof_indices.resize(dofs_per_face);
                          cell->face (face)->get_dof_indices (face_dof_indices,
                                                              cell->active_fe_index ());

                          // Add the computed constraints to the constraint matrix,
                          // assuming the degree of freedom is not already constrained.
                          for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                            {
                              if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                                  && !(constraints.is_constrained (face_dof_indices[dof])))
                                {
                                  constraints.add_line (face_dof_indices[dof]);
                                  if (std::abs (dof_values[dof]) > 1e-13)
                                    {
                                      constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
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
          for (unsigned int i = 0; i < fe_collection.size (); ++i)
            {
              const QGauss<dim-2> reference_edge_quadrature (2 * fe_collection[i].degree + 1);
              for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                {
                  for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_face; ++line)
                    {
                      edge_quadrature_collection.push_back
                      (QProjector<dim>::project_to_face
                       (QProjector<dim - 1>::project_to_face
                        (reference_edge_quadrature, line), face));
                    }
                }
            }

          hp::FEValues<dim> fe_edge_values (mapping_collection, fe_collection,
                                            edge_quadrature_collection,
                                            update_jacobians |
                                            update_JxW_values |
                                            update_quadrature_points |
                                            update_values);

          for (; cell != dof_handler.end (); ++cell)
            {
              if (cell->at_boundary () && cell->is_locally_owned ())
                {
                  for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                    {
                      if (cell->face (face)->boundary_id () == boundary_component)
                        {
                          // If the FE is an FE_Nothing object there is no work to do
                          if (dynamic_cast<const FE_Nothing<dim>*> (&cell->get_fe ()) != 0)
                            {
                              return;
                            }

                          // This is only implemented for FE_Nedelec elements.
                          // If the FE is a FESystem we cannot check this.
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

                          // First compute the projection on the edges.
                          for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_face; ++line)
                            {
                              compute_edge_projection_l2 (cell, face, line, fe_edge_values,
                                                          boundary_function,
                                                          first_vector_component,
                                                          dof_values, dofs_processed);
                            }

                          // If there are higher order shape functions, then we
                          // still need to compute the face projection
                          if (degree > 0)
                            {
                              compute_face_projection_curl_conforming_l2 (cell, face, fe_face_values,
                                                                          boundary_function,
                                                                          first_vector_component,
                                                                          dof_values,
                                                                          dofs_processed);
                            }

                          // Store the computed values in the global vector.
                          face_dof_indices.resize(dofs_per_face);
                          cell->face (face)->get_dof_indices (face_dof_indices,
                                                              cell->active_fe_index ());

                          for (unsigned int dof = 0; dof < dofs_per_face; ++dof)
                            {
                              if (dofs_processed[dof] && constraints.can_store_line (face_dof_indices[dof])
                                  && !(constraints.is_constrained (face_dof_indices[dof])))
                                {
                                  constraints.add_line (face_dof_indices[dof]);

                                  if (std::abs (dof_values[dof]) > 1e-13)
                                    {
                                      constraints.set_inhomogeneity (face_dof_indices[dof], dof_values[dof]);
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
          Assert (false, ExcNotImplemented ());
        }
    }

  }


  template <int dim>
  void
  project_boundary_values_curl_conforming_l2 (const DoFHandler<dim> &dof_handler,
                                              const unsigned int first_vector_component,
                                              const Function<dim> &boundary_function,
                                              const types::boundary_id boundary_component,
                                              ConstraintMatrix &constraints,
                                              const Mapping<dim> &mapping)
  {
    // non-hp version - calls the internal
    // compute_project_boundary_values_curl_conforming_l2() function
    // above after recasting the mapping.

    hp::MappingCollection<dim> mapping_collection (mapping);
    internals::
    compute_project_boundary_values_curl_conforming_l2(dof_handler,
                                                       first_vector_component,
                                                       boundary_function,
                                                       boundary_component,
                                                       constraints,
                                                       mapping_collection);
  }

  template <int dim>
  void
  project_boundary_values_curl_conforming_l2 (const hp::DoFHandler<dim> &dof_handler,
                                              const unsigned int first_vector_component,
                                              const Function<dim> &boundary_function,
                                              const types::boundary_id boundary_component,
                                              ConstraintMatrix &constraints,
                                              const hp::MappingCollection<dim, dim> &mapping_collection)
  {
    // hp version - calls the internal
    // compute_project_boundary_values_curl_conforming_l2() function above.
    internals::
    compute_project_boundary_values_curl_conforming_l2(dof_handler,
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
      const std::vector<Tensor<1,2> > &normals = fe_values.get_all_normal_vectors ();
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
      const std::vector<Tensor<1,3> > &normals = fe_values.get_all_normal_vectors ();
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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
              if (cell->face (face)->boundary_id () == boundary_component)
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



  template <int dim, template <int, int> class DoFHandlerType, int spacedim>
  void
  compute_no_normal_flux_constraints (const DoFHandlerType<dim,spacedim> &dof_handler,
                                      const unsigned int                  first_vector_component,
                                      const std::set<types::boundary_id> &boundary_ids,
                                      ConstraintMatrix                   &constraints,
                                      const Mapping<dim, spacedim>       &mapping)
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

  template <int dim, template <int, int> class DoFHandlerType, int spacedim>
  void
  compute_nonzero_normal_flux_constraints
  (const DoFHandlerType<dim,spacedim>   &dof_handler,
   const unsigned int                    first_vector_component,
   const std::set<types::boundary_id>   &boundary_ids,
   typename FunctionMap<spacedim>::type &function_map,
   ConstraintMatrix                     &constraints,
   const Mapping<dim, spacedim>         &mapping)
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
        std::pair<Tensor<1,dim>, typename DoFHandlerType<dim,spacedim>::active_cell_iterator> >
        DoFToNormalsMap;
    std::map<internal::VectorDoFTuple<dim>, Vector<double> >
    dof_vector_to_b_values;

    DoFToNormalsMap dof_to_normals_map;

    // now loop over all cells and all faces
    typename DoFHandlerType<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    std::set<types::boundary_id>::iterator b_id;
    for (; cell!=endc; ++cell)
      if (!cell->is_artificial())
        for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if ((b_id=boundary_ids.find(cell->face(face_no)->boundary_id()))
              != boundary_ids.end())
            {
              const FiniteElement<dim> &fe = cell->get_fe ();
              typename DoFHandlerType<dim,spacedim>::face_iterator face = cell->face(face_no);

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
                    Tensor<1,dim> normal_vector
                      = (cell->face(face_no)->get_boundary().normal_vector
                         (cell->face(face_no),
                          fe_values.quadrature_point(i)));
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
        typedef std::map<typename DoFHandlerType<dim,spacedim>::active_cell_iterator,
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
            typedef std::map<typename DoFHandlerType<dim,spacedim>::active_cell_iterator,
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
                    tangent = cross_product_3d (normals[0], normals[dim-2]);
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



  template <int dim, template <int, int> class DoFHandlerType, int spacedim>
  void
  compute_normal_flux_constraints (const DoFHandlerType<dim,spacedim> &dof_handler,
                                   const unsigned int                  first_vector_component,
                                   const std::set<types::boundary_id> &boundary_ids,
                                   ConstraintMatrix                   &constraints,
                                   const Mapping<dim, spacedim>       &mapping)
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

  template <int dim, template <int, int> class DoFHandlerType, int spacedim>
  void
  compute_nonzero_tangential_flux_constraints
  (const DoFHandlerType<dim,spacedim>   &dof_handler,
   const unsigned int                    first_vector_component,
   const std::set<types::boundary_id>   &boundary_ids,
   typename FunctionMap<spacedim>::type &function_map,
   ConstraintMatrix                     &constraints,
   const Mapping<dim, spacedim>         &mapping)
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
    for (typename DoFHandlerType<dim,spacedim>::active_cell_iterator cell =
           dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      if (!cell->is_artificial())
        for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if ((b_id=boundary_ids.find(cell->face(face_no)->boundary_id()))
              != boundary_ids.end())
            {
              const FiniteElement<dim> &fe = cell->get_fe();
              typename DoFHandlerType<dim,spacedim>::face_iterator face=cell->face(face_no);

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
    template <int dim, int spacedim, typename Number>
    struct IDScratchData
    {
      IDScratchData (const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                     const dealii::hp::FECollection<dim,spacedim> &fe,
                     const dealii::hp::QCollection<dim> &q,
                     const UpdateFlags update_flags);

      IDScratchData (const IDScratchData &data);

      void resize_vectors (const unsigned int n_q_points,
                           const unsigned int n_components);

      std::vector<Vector<Number> > function_values;
      std::vector<std::vector<Tensor<1,spacedim,Number> > > function_grads;
      std::vector<double> weight_values;
      std::vector<Vector<double> > weight_vectors;

      std::vector<Vector<Number> > psi_values;
      std::vector<std::vector<Tensor<1,spacedim,Number> > > psi_grads;
      std::vector<Number> psi_scalar;

      std::vector<double>         tmp_values;
      std::vector<Vector<double> > tmp_vector_values;
      std::vector<Tensor<1,spacedim> > tmp_gradients;
      std::vector<std::vector<Tensor<1,spacedim> > > tmp_vector_gradients;

      dealii::hp::FEValues<dim,spacedim> x_fe_values;
    };


    template <int dim, int spacedim, typename Number>
    IDScratchData<dim,spacedim,Number>
    ::IDScratchData(const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                    const dealii::hp::FECollection<dim,spacedim> &fe,
                    const dealii::hp::QCollection<dim> &q,
                    const UpdateFlags update_flags)
      :
      x_fe_values(mapping, fe, q, update_flags)
    {}

    template <int dim, int spacedim, typename Number>
    IDScratchData<dim,spacedim,Number>::IDScratchData (const IDScratchData &data)
      :
      x_fe_values(data.x_fe_values.get_mapping_collection(),
                  data.x_fe_values.get_fe_collection(),
                  data.x_fe_values.get_quadrature_collection(),
                  data.x_fe_values.get_update_flags())
    {}

    template <int dim, int spacedim, typename Number>
    void
    IDScratchData<dim,spacedim,Number>::resize_vectors (const unsigned int n_q_points,
                                                        const unsigned int n_components)
    {
      function_values.resize (n_q_points,
                              Vector<Number>(n_components));
      function_grads.resize (n_q_points,
                             std::vector<Tensor<1,spacedim,Number> >(n_components));

      weight_values.resize (n_q_points);
      weight_vectors.resize (n_q_points,
                             Vector<double>(n_components));

      psi_values.resize (n_q_points,
                         Vector<Number>(n_components));
      psi_grads.resize (n_q_points,
                        std::vector<Tensor<1,spacedim,Number> >(n_components));
      psi_scalar.resize (n_q_points);

      tmp_values.resize (n_q_points);
      tmp_vector_values.resize (n_q_points,
                                Vector<double>(n_components));
      tmp_gradients.resize (n_q_points);
      tmp_vector_gradients.resize (n_q_points,
                                   std::vector<Tensor<1,spacedim> >(n_components));
    }

    namespace
    {
      template<typename number>
      double mean_to_double(const number &mean_value)
      {
        return mean_value;
      }

      template<typename number>
      double mean_to_double(const std::complex<number> &mean_value)
      {
        // we need to return double as a norm, but mean value is a complex
        // number. Panick and return real-part while warning the user that
        // he shall never do that.
        Assert ( false, ExcMessage("Mean value norm is not implemented for complex-valued vectors") );
        return mean_value.real();
      }
    }


    // avoid compiling inner function for many vector types when we always
    // really do the same thing by putting the main work into this helper
    // function
    template <int dim, int spacedim, typename Number>
    double
    integrate_difference_inner (const Function<spacedim>   &exact_solution,
                                const NormType              &norm,
                                const Function<spacedim>    *weight,
                                const UpdateFlags            update_flags,
                                const double                 exponent,
                                const unsigned int           n_components,
                                IDScratchData<dim,spacedim,Number> &data)
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
          // points. try to do this as efficient as possible by avoiding a
          // second virtual function call in case the function really has only
          // one component
          //
          // TODO: we have to work a bit here because the Function<dim,double>
          //   interface of the argument denoting the exact function only
          //   provides us with double/Tensor<1,dim> values, rather than
          //   with the correct data type. so evaluate into a temp
          //   object, then copy around
          if (fe_is_system)
            {
              exact_solution.vector_value_list (fe_values.get_quadrature_points(),
                                                data.tmp_vector_values);
              for (unsigned int i=0; i<n_q_points; ++i)
                data.psi_values[i] = data.tmp_vector_values[i];
            }
          else
            {
              exact_solution.value_list (fe_values.get_quadrature_points(),
                                         data.tmp_values);
              for (unsigned int i=0; i<n_q_points; ++i)
                data.psi_values[i](0) = data.tmp_values[i];
            }

          // then subtract finite element fe_function
          for (unsigned int q=0; q<n_q_points; ++q)
            for (unsigned int i=0; i<data.psi_values[q].size(); ++i)
              data.psi_values[q][i] -= data.function_values[q][i];
        }

      // Do the same for gradients, if required
      if (update_flags & update_gradients)
        {
          // try to be a little clever to avoid recursive virtual function
          // calls when calling gradient_list for functions that are really
          // scalar functions
          if (fe_is_system)
            {
              exact_solution.vector_gradient_list (fe_values.get_quadrature_points(),
                                                   data.tmp_vector_gradients);
              for (unsigned int i=0; i<n_q_points; ++i)
                for (unsigned int comp=0; comp<data.psi_grads[i].size(); ++comp)
                  data.psi_grads[i][comp] = data.tmp_vector_gradients[i][comp];
            }
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
                {
                  // compute (f.n) n
                  const Number f_dot_n
                    = (data.psi_grads[q][k] * fe_values.normal_vector(q));
                  const Tensor<1,spacedim,Number> f_dot_n_times_n (f_dot_n * fe_values.normal_vector(q));

                  data.psi_grads[q][k] -= (data.function_grads[q][k] + f_dot_n_times_n);
                }
          else
            for (unsigned int k=0; k<n_components; ++k)
              for (unsigned int q=0; q<n_q_points; ++q)
                for (unsigned int d=0; d<spacedim; ++d)
                  data.psi_grads[q][k][d] -= data.function_grads[q][k][d];
        }

      double diff = 0;
      Number diff_mean = 0;

      // First work on function values:
      switch (norm)
        {
        case mean:
          // Compute values in quadrature points and integrate
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              Number sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += data.psi_values[q](k) * data.weight_vectors[q](k);
              diff_mean += sum * fe_values.JxW(q);
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
                sum += std::pow(
                         static_cast<double>(numbers::NumberTraits<Number>::abs_square(data.psi_values[q](k))),
                         exponent/2.) * data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }

          // Compute the root only if no derivative values are added later
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
                sum += numbers::NumberTraits<Number>::abs_square(data.psi_values[q](k)) *
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
              diff = std::max (diff, double(std::abs(data.psi_values[q](k)*
                                                     data.weight_vectors[q](k))));
          break;

        case H1_seminorm:
        case Hdiv_seminorm:
        case W1p_seminorm:
        case W1infty_seminorm:
          // function values are not used for these norms
          break;

        default:
          Assert (false, ExcNotImplemented());
          break;
        }

      // Now compute terms depending on derivatives:
      switch (norm)
        {
        case W1p_seminorm:
        case W1p_norm:
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              double sum = 0;
              for (unsigned int k=0; k<n_components; ++k)
                sum += std::pow(
                         data.psi_grads[q][k].norm_square(),
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
                sum += data.psi_grads[q][k].norm_square() *
                       data.weight_vectors[q](k);
              diff += sum * fe_values.JxW(q);
            }
          diff = std::sqrt(diff);
          break;

        case Hdiv_seminorm:
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              Assert (n_components >= dim,
                      ExcMessage ("You can only ask for the Hdiv norm for a finite element "
                                  "with at least 'dim' components. In that case, this function "
                                  "will take the divergence of the first 'dim' components."));
              Number sum = 0;
              // take the trace of the derivatives scaled by the weight and square it
              for (unsigned int k=0; k<dim; ++k)
                sum += data.psi_grads[q][k][k] * std::sqrt(data.weight_vectors[q](k));
              diff += numbers::NumberTraits<Number>::abs_square(sum) * fe_values.JxW(q);
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
                t = std::max(t,
                             double(std::abs(data.psi_grads[q][k][d]) *
                                    data.weight_vectors[q](k)));

          // then add seminorm to norm if that had previously been computed
          diff += t;
        }
        break;
        default:
          break;
        }

      if (norm == mean)
        diff = mean_to_double(diff_mean);

      // append result of this cell to the end of the vector
      AssertIsFinite(diff);
      return diff;
    }



    template <int dim, class InVector, class OutVector, typename DoFHandlerType, int spacedim>
    static
    void
    do_integrate_difference (const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                             const DoFHandlerType                              &dof,
                             const InVector                                    &fe_function,
                             const Function<spacedim>                          &exact_solution,
                             OutVector                                         &difference,
                             const dealii::hp::QCollection<dim>                &q,
                             const NormType                                    &norm,
                             const Function<spacedim>                          *weight,
                             const double                                       exponent_1)
    {
      typedef typename InVector::value_type Number;
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

      difference.reinit (dof.get_triangulation().n_active_cells());

      switch (norm)
        {
        case L2_norm:
        case H1_seminorm:
        case H1_norm:
        case Hdiv_seminorm:
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
        case Hdiv_seminorm:
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
      IDScratchData<dim,spacedim, Number> data(mapping, fe_collection, q, update_flags);

      // loop over all cells
      for (typename DoFHandlerType::active_cell_iterator cell = dof.begin_active();
           cell != dof.end(); ++cell)
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
              fe_values.get_function_gradients (fe_function, data.function_grads);

            difference(cell->active_cell_index()) =
              integrate_difference_inner<dim,spacedim, Number> (exact_solution, norm, weight,
                                                                update_flags, exponent,
                                                                n_components, data);
          }
        else
          // the cell is a ghost cell or is artificial. write a zero into the
          // corresponding value of the returned vector
          difference(cell->active_cell_index()) = 0;
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

  template <int dim, int spacedim, class InVector>
  double compute_global_error(const Triangulation<dim,spacedim> &tria,
                              const InVector &cellwise_error,
                              const NormType &norm,
                              const double exponent)
  {
    Assert( cellwise_error.size() == tria.n_active_cells(),
            ExcMessage("input vector cell_error has invalid size!"));
#ifdef DEBUG
    {
      // check that off-processor entries are zero. Otherwise we will compute
      // wrong results below!
      typename InVector::size_type i = 0;
      typename Triangulation<dim,spacedim>::cell_iterator it = tria.begin_active();
      for (; i<cellwise_error.size(); ++i, ++it)
        if (!it->is_locally_owned())
          Assert(std::fabs(cellwise_error[i]) <  1e-20,
                 ExcMessage("cellwise_error of cells that are not locally owned need to be zero!"));
    }
#endif

    MPI_Comm comm = MPI_COMM_SELF;
#ifdef DEAL_II_WITH_MPI
    if (const parallel::Triangulation<dim,spacedim> *ptria =
          dynamic_cast<const parallel::Triangulation<dim,spacedim>*>(&tria))
      comm = ptria->get_communicator();
#endif

    switch (norm)
      {
      case L2_norm:
      case H1_seminorm:
      case H1_norm:
      case Hdiv_seminorm:
      {
        const double local = cellwise_error.l2_norm();
        return std::sqrt (Utilities::MPI::sum (local * local, comm));
      }

      case L1_norm:
      {
        const double local = cellwise_error.l1_norm();
        return Utilities::MPI::sum (local, comm);
      }

      case Linfty_norm:
      case W1infty_seminorm:
      {
        const double local = cellwise_error.linfty_norm();
        return Utilities::MPI::max (local, comm);
      }

      case W1infty_norm:
      {
        AssertThrow(false, ExcMessage(
                      "compute_global_error() is impossible for "
                      "the W1infty_norm. See the documentation for "
                      "NormType::W1infty_norm for more information."));
        return std::numeric_limits<double>::infinity();
      }

      case mean:
      {
        // Note: mean is defined as int_\Omega f = sum_K \int_K f, so we need
        // the sum of the cellwise errors not the Euclidean mean value that
        // is returned by Vector<>::mean_value().
        const double local = cellwise_error.mean_value()
                             * cellwise_error.size();
        return Utilities::MPI::sum (local, comm);
      }

      case Lp_norm:
      case W1p_norm:
      case W1p_seminorm:
      {
        double local = 0;
        typename InVector::size_type i;
        typename Triangulation<dim,spacedim>::cell_iterator it = tria.begin_active();
        for (i = 0; i<cellwise_error.size(); ++i, ++it)
          if (it->is_locally_owned())
            local += std::pow(cellwise_error[i], exponent);

        return std::pow (Utilities::MPI::sum (local, comm), 1./exponent);
      }

      default:
        AssertThrow(false, ExcNotImplemented());
        break;
      }
    return 0.0;
  }

  template <int dim, typename VectorType, int spacedim>
  void
  point_difference (const DoFHandler<dim,spacedim> &dof,
                    const VectorType               &fe_function,
                    const Function<spacedim, typename VectorType::value_type>       &exact_function,
                    Vector<typename VectorType::value_type>                 &difference,
                    const Point<spacedim>          &point)
  {
    point_difference(StaticMappingQ1<dim>::mapping,
                     dof,
                     fe_function,
                     exact_function,
                     difference,
                     point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_difference (const Mapping<dim, spacedim>   &mapping,
                    const DoFHandler<dim,spacedim> &dof,
                    const VectorType               &fe_function,
                    const Function<spacedim, typename VectorType::value_type>       &exact_function,
                    Vector<typename VectorType::value_type>                 &difference,
                    const Point<spacedim>          &point)
  {
    typedef typename VectorType::value_type Number;
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
    std::vector<Vector<Number> > u_value(1, Vector<Number> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    if (fe.n_components() == 1)
      difference(0) = exact_function.value(point);
    else
      exact_function.vector_value(point, difference);

    for (unsigned int i=0; i<difference.size(); ++i)
      difference(i) -= u_value[0](i);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value (const DoFHandler<dim,spacedim> &dof,
               const VectorType               &fe_function,
               const Point<spacedim>          &point,
               Vector<typename VectorType::value_type>                 &value)
  {

    point_value (StaticMappingQ1<dim,spacedim>::mapping,
                 dof,
                 fe_function,
                 point,
                 value);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value (const hp::DoFHandler<dim,spacedim> &dof,
               const VectorType                   &fe_function,
               const Point<spacedim>              &point,
               Vector<typename VectorType::value_type>                     &value)
  {
    point_value(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                dof,
                fe_function,
                point,
                value);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value (const DoFHandler<dim,spacedim> &dof,
               const VectorType               &fe_function,
               const Point<spacedim>          &point)
  {
    return point_value (StaticMappingQ1<dim,spacedim>::mapping,
                        dof,
                        fe_function,
                        point);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value (const hp::DoFHandler<dim,spacedim> &dof,
               const VectorType                   &fe_function,
               const Point<spacedim>              &point)
  {
    return point_value(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                       dof,
                       fe_function,
                       point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value (const Mapping<dim, spacedim>   &mapping,
               const DoFHandler<dim,spacedim> &dof,
               const VectorType               &fe_function,
               const Point<spacedim>          &point,
               Vector<typename VectorType::value_type>                 &value)
  {
    typedef typename VectorType::value_type Number;
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
    std::vector<Vector<Number> > u_value(1, Vector<Number> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value (const hp::MappingCollection<dim, spacedim> &mapping,
               const hp::DoFHandler<dim,spacedim>         &dof,
               const VectorType                           &fe_function,
               const Point<spacedim>                      &point,
               Vector<typename VectorType::value_type>                             &value)
  {
    typedef typename VectorType::value_type Number;
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
    std::vector<Vector<Number> > u_value(1, Vector<Number> (fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value (const Mapping<dim, spacedim>   &mapping,
               const DoFHandler<dim,spacedim> &dof,
               const VectorType               &fe_function,
               const Point<spacedim>          &point)
  {
    Assert(dof.get_fe().n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    Vector<typename VectorType::value_type> value(1);
    point_value(mapping, dof, fe_function, point, value);

    return value(0);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value (const hp::MappingCollection<dim, spacedim> &mapping,
               const hp::DoFHandler<dim,spacedim>         &dof,
               const VectorType                           &fe_function,
               const Point<spacedim>                      &point)
  {
    Assert(dof.get_fe().n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    Vector<typename VectorType::value_type> value(1);
    point_value(mapping, dof, fe_function, point, value);

    return value(0);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient (const DoFHandler<dim,spacedim> &dof,
                  const VectorType               &fe_function,
                  const Point<spacedim>          &point,
                  std::vector<Tensor<1, spacedim, typename VectorType::value_type> > &gradients)
  {

    point_gradient (StaticMappingQ1<dim,spacedim>::mapping,
                    dof,
                    fe_function,
                    point,
                    gradients);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient (const hp::DoFHandler<dim,spacedim> &dof,
                  const VectorType                   &fe_function,
                  const Point<spacedim>              &point,
                  std::vector<Tensor<1, spacedim, typename VectorType::value_type> > &gradients)
  {
    point_gradient(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                   dof,
                   fe_function,
                   point,
                   gradients);
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient (const DoFHandler<dim,spacedim> &dof,
                  const VectorType               &fe_function,
                  const Point<spacedim>          &point)
  {
    return point_gradient (StaticMappingQ1<dim,spacedim>::mapping,
                           dof,
                           fe_function,
                           point);
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient (const hp::DoFHandler<dim,spacedim> &dof,
                  const VectorType                   &fe_function,
                  const Point<spacedim>              &point)
  {
    return point_gradient(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                          dof,
                          fe_function,
                          point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient (const Mapping<dim, spacedim>   &mapping,
                  const DoFHandler<dim,spacedim> &dof,
                  const VectorType               &fe_function,
                  const Point<spacedim>          &point,
                  std::vector<Tensor<1, spacedim, typename VectorType::value_type> > &gradient)
  {
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(gradient.size() == fe.n_components(),
           ExcDimensionMismatch(gradient.size(), fe.n_components()));

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

    FEValues<dim> fe_values(mapping, fe, quadrature, update_gradients);
    fe_values.reinit(cell_point.first);

    // then use this to get the gradients of
    // the given fe_function at this point
    typedef typename VectorType::value_type Number;
    std::vector<std::vector<Tensor<1, dim, Number> > >
    u_gradient(1, std::vector<Tensor<1, dim, Number> > (fe.n_components()));
    fe_values.get_function_gradients(fe_function, u_gradient);

    gradient = u_gradient[0];
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient (const hp::MappingCollection<dim, spacedim> &mapping,
                  const hp::DoFHandler<dim,spacedim>         &dof,
                  const VectorType                           &fe_function,
                  const Point<spacedim>                      &point,
                  std::vector<Tensor<1, spacedim, typename VectorType::value_type> > &gradient)
  {
    typedef typename VectorType::value_type Number;
    const hp::FECollection<dim, spacedim> &fe = dof.get_fe();

    Assert(gradient.size() == fe.n_components(),
           ExcDimensionMismatch(gradient.size(), fe.n_components()));

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
    hp::FEValues<dim, spacedim> hp_fe_values(mapping, fe, hp::QCollection<dim>(quadrature), update_gradients);
    hp_fe_values.reinit(cell_point.first);
    const FEValues<dim, spacedim> &fe_values = hp_fe_values.get_present_fe_values();

    std::vector<std::vector<Tensor<1, dim, Number> > >
    u_gradient(1, std::vector<Tensor<1, dim, Number> > (fe.n_components()));
    fe_values.get_function_gradients(fe_function, u_gradient);

    gradient = u_gradient[0];
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient (const Mapping<dim, spacedim>   &mapping,
                  const DoFHandler<dim,spacedim> &dof,
                  const VectorType               &fe_function,
                  const Point<spacedim>          &point)
  {
    Assert(dof.get_fe().n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    std::vector<Tensor<1, dim, typename VectorType::value_type> > gradient(1);
    point_gradient (mapping, dof, fe_function, point, gradient);

    return gradient[0];
  }



  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient (const hp::MappingCollection<dim, spacedim> &mapping,
                  const hp::DoFHandler<dim,spacedim>         &dof,
                  const VectorType                           &fe_function,
                  const Point<spacedim>                      &point)
  {
    Assert(dof.get_fe().n_components() == 1,
           ExcMessage ("Finite element is not scalar as is necessary for this function"));

    std::vector<Tensor<1, dim, typename VectorType::value_type> > gradient(1);
    point_gradient (mapping, dof, fe_function, point, gradient);

    return gradient[0];
  }



  template <typename VectorType>
  void
  subtract_mean_value(VectorType              &v,
                      const std::vector<bool> &p_select)
  {
    if (p_select.size() == 0)
      {
        // In case of an empty boolean mask operate on the whole vector:
        v.add( - v.mean_value() );
      }
    else
      {
        // This function is not implemented for distributed vectors.
        Assert(!v.supports_distributed_data, ExcNotImplemented());

        const unsigned int n = v.size();

        Assert(p_select.size() == n,
               ExcDimensionMismatch(p_select.size(), n));

        typename VectorType::value_type s = 0.;
        unsigned int counter = 0;
        for (unsigned int i=0; i<n; ++i)
          if (p_select[i])
            {
              typename VectorType::value_type vi = v(i);
              s += vi;
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

  namespace
  {
    template <typename Number>
    void set_possibly_complex_number(const double &r,
                                     const double &,
                                     Number &n)
    {
      n = r;
    }



    template <typename Type>
    void set_possibly_complex_number(const double &r,
                                     const double &i,
                                     std::complex<Type> &n)
    {
      n = std::complex<Type>(r,i);
    }
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value (const Mapping<dim, spacedim>   &mapping,
                      const DoFHandler<dim,spacedim> &dof,
                      const Quadrature<dim>          &quadrature,
                      const VectorType               &v,
                      const unsigned int             component)
  {
    typedef typename VectorType::value_type Number;
    Assert (v.size() == dof.n_dofs(),
            ExcDimensionMismatch (v.size(), dof.n_dofs()));
    Assert (component < dof.get_fe().n_components(),
            ExcIndexRange(component, 0, dof.get_fe().n_components()));

    FEValues<dim,spacedim> fe(mapping, dof.get_fe(), quadrature,
                              UpdateFlags(update_JxW_values
                                          | update_values));

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell;
    std::vector<Vector<Number> > values(quadrature.size(),
                                        Vector<Number> (dof.get_fe().n_components()));

    Number mean = Number();
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

#ifdef DEAL_II_WITH_MPI
    // if this was a distributed DoFHandler, we need to do the reduction
    // over the entire domain
    if (const parallel::Triangulation<dim,spacedim> *
        p_triangulation
        = dynamic_cast<const parallel::Triangulation<dim,spacedim> *>(&dof.get_triangulation()))
      {
        // The type used to store the elements of the global vector may be a
        // real or a complex number. Do the global reduction always with real
        // and imaginary types so that we don't have to distinguish, and to this
        // end just copy everything into a complex number and, later, back into
        // the original data type.
        std::complex<double> mean_double = mean;
        double my_values[3] = { mean_double.real(), mean_double.imag(), area };
        double global_values[3];

        MPI_Allreduce (my_values, global_values, 3, MPI_DOUBLE,
                       MPI_SUM,
                       p_triangulation->get_communicator());

        set_possibly_complex_number(global_values[0], global_values[1],
                                    mean);
        area = global_values[2];
      }
#endif

    return (mean/area);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  compute_mean_value (const DoFHandler<dim,spacedim> &dof,
                      const Quadrature<dim>          &quadrature,
                      const VectorType               &v,
                      const unsigned int             component)
  {
    return compute_mean_value(StaticMappingQ1<dim,spacedim>::mapping, dof, quadrature, v, component);
  }


  template<typename DoFHandlerType, typename VectorType>
  void get_position_vector(const DoFHandlerType &dh,
                           VectorType           &vector,
                           const ComponentMask  &mask)
  {
    AssertDimension(vector.size(), dh.n_dofs());

    const unsigned int dim=DoFHandlerType::dimension;
    const unsigned int spacedim=DoFHandlerType::space_dimension;
    const FiniteElement<dim, spacedim> &fe = dh.get_fe();


    // Construct default fe_mask;
    const ComponentMask fe_mask(mask.size() ? mask :
                                ComponentMask(fe.get_nonzero_components(0).size(), true));

    AssertDimension(fe_mask.size(), fe.get_nonzero_components(0).size());

    std::vector<unsigned int> fe_to_real(fe_mask.size(), numbers::invalid_unsigned_int);
    unsigned int size = 0;
    for (unsigned int i=0; i<fe_mask.size(); ++i)
      {
        if (fe_mask[i])
          fe_to_real[i] = size++;
      }
    Assert(size == spacedim,
           ExcMessage("The Component Mask you provided is invalid. It has to select exactly spacedim entries."));


    if ( fe.has_support_points() )
      {
        typename DoFHandlerType::active_cell_iterator cell;
        const Quadrature<dim> quad(fe.get_unit_support_points());

        MappingQ<dim,spacedim> map_q(fe.degree);
        FEValues<dim,spacedim> fe_v(map_q, fe, quad, update_quadrature_points);
        std::vector<types::global_dof_index> dofs(fe.dofs_per_cell);

        AssertDimension(fe.dofs_per_cell, fe.get_unit_support_points().size());
        Assert(fe.is_primitive(), ExcMessage("FE is not Primitive! This won't work."));

        for (cell = dh.begin_active(); cell != dh.end(); ++cell)
          {
            fe_v.reinit(cell);
            cell->get_dof_indices(dofs);
            const std::vector<Point<spacedim> > &points = fe_v.get_quadrature_points();
            for (unsigned int q = 0; q < points.size(); ++q)
              {
                unsigned int comp = fe.system_to_component_index(q).first;
                if (fe_mask[comp])
                  vector(dofs[q]) = points[q][fe_to_real[comp]];
              }
          }
      }
    else
      {
        // Construct a FiniteElement with FE_Q^spacedim, and call this
        // function again.
        //
        // Once we have this, interpolate with the given finite element
        // to get a Mapping which is interpolatory at the support points
        // of FE_Q(fe.degree())
        const FESystem<dim,spacedim> *fe_system = dynamic_cast<const FESystem<dim, spacedim> *>(&fe);
        Assert(fe_system, ExcNotImplemented());
        unsigned int degree = numbers::invalid_unsigned_int;

        // Get information about the blocks
        for (unsigned int i=0; i<fe_mask.size(); ++i)
          if (fe_mask[i])
            {
              const unsigned int base_i = fe_system->component_to_base_index(i).first;
              Assert(degree == numbers::invalid_unsigned_int ||
                     degree == fe_system->base_element(base_i).degree,
                     ExcNotImplemented());
              Assert(fe_system->base_element(base_i).is_primitive(),
                     ExcNotImplemented());
              degree = fe_system->base_element(base_i).degree;
            }

        // We create an intermediate FE_Q vector space, and then
        // interpolate from that vector space to this one, by
        // carefully selecting the right components.

        FESystem<dim,spacedim> feq(FE_Q<dim,spacedim>(degree), spacedim);
        DoFHandlerType dhq(dh.get_triangulation());
        dhq.distribute_dofs(feq);
        Vector<double> eulerq(dhq.n_dofs());
        const ComponentMask maskq(spacedim, true);
        get_position_vector(dhq, eulerq);

        FullMatrix<double> transfer(fe.dofs_per_cell, feq.dofs_per_cell);
        FullMatrix<double> local_transfer(feq.dofs_per_cell);
        const std::vector<Point<dim> > &points = feq.get_unit_support_points();

        // Here we construct the interpolation matrix from
        // FE_Q^spacedim to the FiniteElement used by
        // euler_dof_handler.
        //
        // In order to construct such interpolation matrix, we have to
        // solve the following system:
        //
        // v_j phi_j(q_i) = w_k psi_k(q_i) = w_k delta_ki = w_i
        //
        // where psi_k are the basis functions for fe_q, and phi_i are
        // the basis functions of the target space while q_i are the
        // support points for the fe_q space. With this choice of
        // interpolation points, on the matrices is the identity
        // matrix, and we have to invert only one matrix. The
        // resulting vector will be interpolatory at the support
        // points of fe_q, even if the finite element does not have
        // support points.
        //
        // Morally, we should invert the matrix T_ij = phi_i(q_j),
        // however in general this matrix is not invertible, since
        // there may be components which do not contribute to the
        // displacement vector. Since we are not interested in those
        // components, we construct a square matrix with the same
        // number of components of the FE_Q system. The FE_Q system
        // was constructed above in such a way that the polynomial
        // degree of the FE_Q system and that of the given FE are the
        // same on the cell, which should guarantee that, for the
        // displacement components only, the interpolation matrix is
        // invertible. We construct a mapping between indices first,
        // and check that this is the case. If not, we bail out, not
        // knowing what to do in this case.

        std::vector<unsigned int> fe_to_feq(fe.dofs_per_cell, numbers::invalid_unsigned_int);
        unsigned int index=0;
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          if (fe_mask[fe.system_to_component_index(i).first])
            fe_to_feq[i] = index++;

        // If index is not the same as feq.dofs_per_cell, we won't
        // know how to invert the resulting matrix. Bail out.
        Assert(index == feq.dofs_per_cell, ExcNotImplemented());

        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          {
            const unsigned int comp_j = fe.system_to_component_index(j).first;
            if (fe_mask[comp_j])
              for (unsigned int i=0; i<points.size(); ++i)
                {
                  if ( fe_to_real[comp_j] == feq.system_to_component_index(i).first)
                    local_transfer(i, fe_to_feq[j]) = fe.shape_value(j, points[i]);
                }
          }

        // Now we construct the rectangular interpolation matrix. This
        // one is filled only with the information from the components
        // of the displacement. The rest is set to zero.
        local_transfer.invert(local_transfer);
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          if (fe_to_feq[i] != numbers::invalid_unsigned_int)
            for (unsigned int j=0; j<feq.dofs_per_cell; ++j)
              transfer(i, j) = local_transfer(fe_to_feq[i], j);

        // The interpolation matrix is then passed to the
        // VectorTools::interpolate() function to generate the correct
        // interpolation.
        interpolate(dhq, dh, transfer, eulerq, vector);
      }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
