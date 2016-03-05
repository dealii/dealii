// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#include <deal.II/base/numbers.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/base/std_cxx11/bind.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace
{
  template <typename CellIterator>
  inline
  void advance_by_n (CellIterator &cell,
                     const unsigned int n)
  {
    // store a pointer to the end iterator, since we can't get at it any more
    // once cell is already the end iterator (in that case dereferencing
    // cell-> triggers an assertion)
    const CellIterator endc = cell->get_dof_handler().end();
    for (unsigned int t=0; ((t<n) && (cell!=endc)); ++t, ++cell)
      ;
  }
}


namespace internal
{
  namespace
  {
    /**
     * All small temporary data objects that are needed once per thread by the
     * several functions of the error estimator are gathered in this struct.
     * The reason for this structure is mainly that we have a number of
     * functions that operate on cells or faces and need a number of small
     * temporary data objects. Since these functions may run in parallel, we
     * cannot make these objects member variables of the enclosing class. On
     * the other hand, declaring them locally in each of these functions would
     * require their reallocating every time we visit the next cell or face,
     * which we found can take a significant amount of time if it happens
     * often even in the single threaded case (10-20 per cent in our
     * measurements); however, most importantly, memory allocation requires
     * synchronisation in multithreaded mode. While that is done by the C++
     * library and has not to be handcoded, it nevertheless seriously damages
     * the ability to efficiently run the functions of this class in parallel,
     * since they are quite often blocked by these synchronisation points,
     * slowing everything down by a factor of two or three.
     *
     * Thus, every thread gets an instance of this class to work with and
     * needs not allocate memory itself, or synchronise with other threads.
     *
     * The sizes of the arrays are initialized with the maximal number of
     * entries necessary for the hp case. Within the loop over individual
     * cells, we then resize the arrays as necessary. Since for std::vector
     * resizing to a smaller size doesn't imply memory allocation, this is
     * fast.
     */
    template <typename DoFHandlerType,typename number>
    struct ParallelData
    {
      static const unsigned int dim      = DoFHandlerType::dimension;
      static const unsigned int spacedim = DoFHandlerType::space_dimension;

      /**
       * The finite element to be used.
       */
      const dealii::hp::FECollection<dim,spacedim> finite_element;

      /**
       * The quadrature formulas to be used for the faces.
       */
      const dealii::hp::QCollection<dim-1> face_quadratures;

      /**
       * FEFaceValues objects to integrate over the faces of the current and
       * potentially of neighbor cells.
       */
      dealii::hp::FEFaceValues<dim,spacedim>    fe_face_values_cell;
      dealii::hp::FEFaceValues<dim,spacedim>    fe_face_values_neighbor;
      dealii::hp::FESubfaceValues<dim,spacedim> fe_subface_values;

      /**
       * A vector to store the jump of the normal vectors in the quadrature
       * points for each of the solution vectors (i.e. a temporary value).
       * This vector is not allocated inside the functions that use it, but
       * rather globally, since memory allocation is slow, in particular in
       * presence of multiple threads where synchronisation makes things even
       * slower.
       */
      std::vector<std::vector<std::vector<number> > > phi;

      /**
       * A vector for the gradients of the finite element function on one cell
       *
       * Let psi be a short name for <tt>a grad u_h</tt>, where the third
       * index be the component of the finite element, and the second index
       * the number of the quadrature point. The first index denotes the index
       * of the solution vector.
       */
      std::vector<std::vector<std::vector<Tensor<1,spacedim,number> > > > psi;

      /**
       * The same vector for a neighbor cell
       */
      std::vector<std::vector<std::vector<Tensor<1,spacedim,number> > > > neighbor_psi;

      /**
       * The normal vectors of the finite element function on one face
       */
      std::vector<Tensor<1,spacedim> > normal_vectors;

      /**
       * Normal vectors of the opposing face.
       */
      std::vector<Tensor<1,spacedim> > neighbor_normal_vectors;

      /**
       * Two arrays needed for the values of coefficients in the jumps, if
       * they are given.
       */
      std::vector<double>                  coefficient_values1;
      std::vector<dealii::Vector<double> > coefficient_values;

      /**
       * Array for the products of Jacobian determinants and weights of
       * quadraturs points.
       */
      std::vector<double>          JxW_values;

      /**
       * The subdomain id we are to care for.
       */
      const types::subdomain_id subdomain_id;
      /**
       * The material id we are to care for.
       */
      const types::material_id material_id;

      /**
       * Some more references to input data to the
       * KellyErrorEstimator::estimate() function.
       */
      const typename FunctionMap<spacedim>::type *neumann_bc;
      const ComponentMask                component_mask;
      const Function<spacedim>                   *coefficients;

      /**
       * Constructor.
       */
      template <class FE>
      ParallelData (const FE                                           &fe,
                    const dealii::hp::QCollection<dim-1>               &face_quadratures,
                    const dealii::hp::MappingCollection<dim,spacedim> &mapping,
                    const bool         need_quadrature_points,
                    const unsigned int n_solution_vectors,
                    const types::subdomain_id subdomain_id,
                    const types::material_id material_id,
                    const typename FunctionMap<spacedim>::type *neumann_bc,
                    const ComponentMask                component_mask,
                    const Function<spacedim>                   *coefficients);

      /**
       * Resize the arrays so that they fit the number of quadrature points
       * associated with the given finite element index into the hp
       * collections.
       */
      void resize (const unsigned int active_fe_index);
    };


    template <typename DoFHandlerType,typename number>
    template <class FE>
    ParallelData<DoFHandlerType,number>::
    ParallelData
    (const FE                                           &fe,
     const dealii::hp::QCollection<dim-1>               &face_quadratures,
     const dealii::hp::MappingCollection<dim, spacedim> &mapping,
     const bool                                          need_quadrature_points,
     const unsigned int                                  n_solution_vectors,
     const types::subdomain_id                           subdomain_id,
     const types::material_id                            material_id,
     const typename FunctionMap<spacedim>::type         *neumann_bc,
     const ComponentMask                                 component_mask,
     const Function<spacedim>                           *coefficients)
      :
      finite_element (fe),
      face_quadratures (face_quadratures),
      fe_face_values_cell (mapping,
                           finite_element,
                           face_quadratures,
                           update_gradients      |
                           update_JxW_values     |
                           (need_quadrature_points  ?
                            update_quadrature_points :
                            UpdateFlags()) |
                           update_normal_vectors),
      fe_face_values_neighbor (mapping,
                               finite_element,
                               face_quadratures,
                               update_gradients|
                               update_normal_vectors),
      fe_subface_values (mapping,
                         finite_element,
                         face_quadratures,
                         update_gradients|
                         update_normal_vectors),
      phi (n_solution_vectors,
           std::vector<std::vector<number> >
           (face_quadratures.max_n_quadrature_points(),
            std::vector<number> (fe.n_components()))),
      psi (n_solution_vectors,
           std::vector<std::vector<Tensor<1,spacedim,number> > >
           (face_quadratures.max_n_quadrature_points(),
            std::vector<Tensor<1,spacedim,number> > (fe.n_components()))),
      neighbor_psi (n_solution_vectors,
                    std::vector<std::vector<Tensor<1,spacedim,number> > >
                    (face_quadratures.max_n_quadrature_points(),
                     std::vector<Tensor<1,spacedim,number> > (fe.n_components()))),
      normal_vectors (face_quadratures.max_n_quadrature_points()),
      neighbor_normal_vectors (face_quadratures.max_n_quadrature_points()),
      coefficient_values1 (face_quadratures.max_n_quadrature_points()),
      coefficient_values (face_quadratures.max_n_quadrature_points(),
                          dealii::Vector<double> (fe.n_components())),
      JxW_values (face_quadratures.max_n_quadrature_points()),
      subdomain_id (subdomain_id),
      material_id (material_id),
      neumann_bc (neumann_bc),
      component_mask (component_mask),
      coefficients (coefficients)
    {}



    template <typename DoFHandlerType, typename number>
    void
    ParallelData<DoFHandlerType,number>::resize (const unsigned int active_fe_index)
    {
      const unsigned int n_q_points   = face_quadratures[active_fe_index].size();
      const unsigned int n_components = finite_element.n_components();

      normal_vectors.resize(n_q_points);
      neighbor_normal_vectors.resize(n_q_points);
      coefficient_values1.resize(n_q_points);
      coefficient_values.resize(n_q_points);
      JxW_values.resize(n_q_points);

      for (unsigned int i=0; i<phi.size(); ++i)
        {
          phi[i].resize(n_q_points);
          psi[i].resize(n_q_points);
          neighbor_psi[i].resize(n_q_points);

          for (unsigned int qp=0; qp<n_q_points; ++qp)
            {
              phi[i][qp].resize(n_components);
              psi[i][qp].resize(n_components);
              neighbor_psi[i][qp].resize(n_components);
            }
        }

      for (unsigned int qp=0; qp<n_q_points; ++qp)
        coefficient_values[qp].reinit(n_components);
    }



    /**
     * Copy data from the local_face_integrals map of a single ParallelData
     * object into a global such map. This is the copier stage of a WorkStream
     * pipeline.
     */
    template <typename DoFHandlerType>
    void
    copy_local_to_global
    (const std::map<typename DoFHandlerType::face_iterator,std::vector<double> > &local_face_integrals,
     std::map<typename DoFHandlerType::face_iterator,std::vector<double> > &face_integrals)
    {

      // now copy locally computed elements into the global map
      for (typename std::map<typename DoFHandlerType::face_iterator,std::vector<double> >::const_iterator
           p=local_face_integrals.begin();
           p!=local_face_integrals.end();
           ++p)
        {
          // double check that the element does not already exists in the
          // global map
          Assert (face_integrals.find (p->first) == face_integrals.end(),
                  ExcInternalError());

          for (unsigned int i=0; i<p->second.size(); ++i)
            {
              Assert (numbers::is_finite(p->second[i]), ExcInternalError());
              Assert (p->second[i] >= 0, ExcInternalError());
            }

          face_integrals[p->first] = p->second;
        }
    }


    /**
     * Actually do the computation based on the evaluated gradients in
     * ParallelData.
     */
    template <typename DoFHandlerType, typename number>
    std::vector<double>
    integrate_over_face
    (ParallelData<DoFHandlerType,number>                 &parallel_data,
     const typename DoFHandlerType::face_iterator        &face,
     dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_cell)
    {
      const unsigned int n_q_points         = parallel_data.psi[0].size(),
                         n_components       = parallel_data.finite_element.n_components(),
                         n_solution_vectors = parallel_data.psi.size();

      // now psi contains the following:
      // - for an internal face, psi=[grad u]
      // - for a neumann boundary face, psi=grad u
      // each component being the mentioned value at one of the quadrature
      // points

      // next we have to multiply this with the normal vector. Since we have
      // taken the difference of gradients for internal faces, we may chose
      // the normal vector of one cell, taking that of the neighbor would only
      // change the sign. We take the outward normal.

      parallel_data.normal_vectors =
        fe_face_values_cell.get_present_fe_values().get_all_normal_vectors();

      for (unsigned int n=0; n<n_solution_vectors; ++n)
        for (unsigned int component=0; component<n_components; ++component)
          for (unsigned int point=0; point<n_q_points; ++point)
            parallel_data.phi[n][point][component]
              = (parallel_data.psi[n][point][component] *
                 parallel_data.normal_vectors[point]);

      if (face->at_boundary() == false)
        {
          // compute the jump in the gradients

          for (unsigned int n=0; n<n_solution_vectors; ++n)
            for (unsigned int component=0; component<n_components; ++component)
              for (unsigned int p=0; p<n_q_points; ++p)
                parallel_data.phi[n][p][component]
                += (parallel_data.neighbor_psi[n][p][component] *
                    parallel_data.neighbor_normal_vectors[p]);
        }

      // if a coefficient was given: use that to scale the jump in the
      // gradient
      if (parallel_data.coefficients != 0)
        {
          // scalar coefficient
          if (parallel_data.coefficients->n_components == 1)
            {
              parallel_data.coefficients
              ->value_list (fe_face_values_cell.get_present_fe_values()
                            .get_quadrature_points(),
                            parallel_data.coefficient_values1);
              for (unsigned int n=0; n<n_solution_vectors; ++n)
                for (unsigned int component=0; component<n_components; ++component)
                  for (unsigned int point=0; point<n_q_points; ++point)
                    parallel_data.phi[n][point][component] *=
                      parallel_data.coefficient_values1[point];
            }
          else
            // vector-valued coefficient
            {
              parallel_data.coefficients
              ->vector_value_list (fe_face_values_cell.get_present_fe_values()
                                   .get_quadrature_points(),
                                   parallel_data.coefficient_values);
              for (unsigned int n=0; n<n_solution_vectors; ++n)
                for (unsigned int component=0; component<n_components; ++component)
                  for (unsigned int point=0; point<n_q_points; ++point)
                    parallel_data.phi[n][point][component] *=
                      parallel_data.coefficient_values[point](component);
            }
        }


      if (face->at_boundary() == true)
        // neumann boundary face. compute difference between normal derivative
        // and boundary function
        {
          const types::boundary_id boundary_id = face->boundary_id();

          Assert (parallel_data.neumann_bc->find(boundary_id) !=
                  parallel_data.neumann_bc->end(),
                  ExcInternalError ());
          // get the values of the boundary function at the quadrature points
          if (n_components == 1)
            {
              std::vector<double> g(n_q_points);
              parallel_data.neumann_bc->find(boundary_id)->second
              ->value_list (fe_face_values_cell.get_present_fe_values()
                            .get_quadrature_points(), g);

              for (unsigned int n=0; n<n_solution_vectors; ++n)
                for (unsigned int point=0; point<n_q_points; ++point)
                  parallel_data.phi[n][point][0] -= g[point];
            }
          else
            {
              std::vector<dealii::Vector<double> >
              g(n_q_points, dealii::Vector<double>(n_components));
              parallel_data.neumann_bc->find(boundary_id)->second
              ->vector_value_list (fe_face_values_cell.get_present_fe_values()
                                   .get_quadrature_points(),
                                   g);

              for (unsigned int n=0; n<n_solution_vectors; ++n)
                for (unsigned int component=0; component<n_components; ++component)
                  for (unsigned int point=0; point<n_q_points; ++point)
                    parallel_data.phi[n][point][component] -= g[point](component);
            }
        }




      // now phi contains the following:
      // - for an internal face, phi=[a du/dn]
      // - for a neumann boundary face, phi=a du/dn-g
      // each component being the mentioned value at one of the quadrature
      // points

      parallel_data.JxW_values
        = fe_face_values_cell.get_present_fe_values().get_JxW_values();

      // take the square of the phi[i] for integration, and sum up
      std::vector<double> face_integral (n_solution_vectors, 0);
      for (unsigned int n=0; n<n_solution_vectors; ++n)
        for (unsigned int component=0; component<n_components; ++component)
          if (parallel_data.component_mask[component] == true)
            for (unsigned int p=0; p<n_q_points; ++p)
              face_integral[n] += numbers::NumberTraits<number>::abs_square(parallel_data.phi[n][p][component]) *
                                  parallel_data.JxW_values[p];

      return face_integral;
    }

    /**
     * A factor to scale the integral for the face at the boundary. Used for
     * Neumann BC.
     */
    template <typename DoFHandlerType>
    double boundary_face_factor(const typename DoFHandlerType::active_cell_iterator &cell,
                                const unsigned int                       face_no,
                                const dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_cell,
                                const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      switch (strategy)
        {
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter_over_24 :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::face_diameter_over_twice_max_degree :
        {
          const double cell_degree = fe_face_values_cell.get_fe_collection()[cell->active_fe_index()].degree;
          return cell->face(face_no)->diameter() / cell_degree;
        }
        default:
        {
          Assert (false, ExcNotImplemented());
          return -std::numeric_limits<double>::max();
        }
        }
    }


    /**
     * A factor to scale the integral for the regular face.
     */
    template <typename DoFHandlerType>
    double regular_face_factor(const typename DoFHandlerType::active_cell_iterator &cell,
                               const unsigned int                       face_no,
                               const dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_cell,
                               const dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_neighbor,
                               const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      switch (strategy)
        {
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter_over_24 :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::face_diameter_over_twice_max_degree :
        {
          const double cell_degree     = fe_face_values_cell.get_fe_collection()[cell->active_fe_index()].degree;
          const double neighbor_degree = fe_face_values_neighbor.get_fe_collection()[cell->neighbor(face_no)->active_fe_index()].degree;
          return cell->face(face_no)->diameter() / std::max(cell_degree,neighbor_degree) / 2.0;
        }
        default:
        {
          Assert (false, ExcNotImplemented());
          return -std::numeric_limits<double>::max();
        }
        }
    }

    /**
     * A factor to scale the integral for the irregular face.
     */
    template <typename DoFHandlerType>
    double irregular_face_factor(const typename DoFHandlerType::active_cell_iterator &cell,
                                 const typename DoFHandlerType::active_cell_iterator &neighbor_child,
                                 const unsigned int                       face_no,
                                 const unsigned int                       subface_no,
                                 const dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values,
                                 dealii::hp::FESubfaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension>    &fe_subface_values,
                                 const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      switch (strategy)
        {
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter_over_24 :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter :
        {
          return 1.0;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::face_diameter_over_twice_max_degree :
        {
          const double cell_degree = fe_face_values.get_fe_collection()[cell->active_fe_index()].degree;
          const double neighbor_child_degree = fe_subface_values.get_fe_collection()[neighbor_child->active_fe_index()].degree;
          return cell->face(face_no)->child(subface_no)->diameter()/std::max(neighbor_child_degree,cell_degree)/2.0;
        }
        default:
        {
          Assert (false, ExcNotImplemented());
          return -std::numeric_limits<double>::max();
        }
        }
    }

    /**
     * A factor used when summing up all the contribution from different faces
     * of each cell.
     */
    template <typename DoFHandlerType>
    double cell_factor(const typename DoFHandlerType::active_cell_iterator &cell,
                       const unsigned int                       /*face_no*/,
                       const DoFHandlerType                    &/*dof_handler*/,
                       const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      switch (strategy)
        {
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter_over_24 :
        {
          return cell->diameter()/24;
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::cell_diameter :
        {
          return cell->diameter();
        }
        case KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::face_diameter_over_twice_max_degree :
        {
          return 1.0;
        }
        default:
        {
          Assert (false, ExcNotImplemented());
          return -std::numeric_limits<double>::max();
        }
        }
    }



    /**
     * Actually do the computation on a face which has no hanging nodes (it is
     * regular), i.e. either on the other side there is nirvana (face is at
     * boundary), or the other side's refinement level is the same as that of
     * this side, then handle the integration of these both cases together.
     */
    template <typename InputVector, typename DoFHandlerType>
    void
    integrate_over_regular_face (const std::vector<const InputVector *>   &solutions,
                                 ParallelData<DoFHandlerType, typename InputVector::value_type> &parallel_data,
                                 std::map<typename DoFHandlerType::face_iterator,std::vector<double> > &local_face_integrals,
                                 const typename DoFHandlerType::active_cell_iterator &cell,
                                 const unsigned int                       face_no,
                                 dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_cell,
                                 dealii::hp::FEFaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_face_values_neighbor,
                                 const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      const unsigned int dim = DoFHandlerType::dimension;
      (void)dim;

      const typename DoFHandlerType::face_iterator face = cell->face(face_no);
      const unsigned int n_solution_vectors = solutions.size();


      // initialize data of the restriction
      // of this cell to the present face
      fe_face_values_cell.reinit (cell, face_no,
                                  cell->active_fe_index());

      // get gradients of the finite element
      // function on this cell
      for (unsigned int n=0; n<n_solution_vectors; ++n)
        fe_face_values_cell.get_present_fe_values()
        .get_function_gradients (*solutions[n], parallel_data.psi[n]);

      double factor;
      // now compute over the other side of the face
      if (face->at_boundary() == false)
        // internal face; integrate jump of gradient across this face
        {
          Assert (cell->neighbor(face_no).state() == IteratorState::valid,
                  ExcInternalError());

          const typename DoFHandlerType::active_cell_iterator neighbor = cell->neighbor(face_no);

          // find which number the current face has relative to the
          // neighboring cell
          const unsigned int neighbor_neighbor
            = cell->neighbor_of_neighbor (face_no);
          Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell,
                  ExcInternalError());

          // get restriction of finite element function of @p{neighbor} to the
          // common face. in the hp case, use the quadrature formula that
          // matches the one we would use for the present cell
          fe_face_values_neighbor.reinit (neighbor, neighbor_neighbor,
                                          cell->active_fe_index());

          factor = regular_face_factor<DoFHandlerType>(cell,face_no,
                                                       fe_face_values_cell,fe_face_values_neighbor,
                                                       strategy);

          // get gradients on neighbor cell
          for (unsigned int n=0; n<n_solution_vectors; ++n)
            {
              fe_face_values_neighbor.get_present_fe_values()
              .get_function_gradients (*solutions[n],
                                       parallel_data.neighbor_psi[n]);
            }

          parallel_data.neighbor_normal_vectors =
            fe_face_values_neighbor.get_present_fe_values().get_all_normal_vectors();

        }
      else
        {
          factor = boundary_face_factor<DoFHandlerType>(cell,face_no,
                                                        fe_face_values_cell,
                                                        strategy);
        }

      // now go to the generic function that does all the other things
      local_face_integrals[face] =
        integrate_over_face (parallel_data, face,
                             fe_face_values_cell);

      for (unsigned int i = 0; i < local_face_integrals[face].size(); i++)
        local_face_integrals[face][i] *= factor;
    }




    /**
     * The same applies as for the function above, except that integration is
     * over face @p face_no of @p cell, where the respective neighbor is
     * refined, so that the integration is a bit more complex.
     */
    template <typename InputVector, typename DoFHandlerType>
    void
    integrate_over_irregular_face (const std::vector<const InputVector *>   &solutions,
                                   ParallelData<DoFHandlerType, typename InputVector::value_type> &parallel_data,
                                   std::map<typename DoFHandlerType::face_iterator,std::vector<double> > &local_face_integrals,
                                   const typename DoFHandlerType::active_cell_iterator    &cell,
                                   const unsigned int                          face_no,
                                   dealii::hp::FEFaceValues<DoFHandlerType::dimension,DoFHandlerType::space_dimension>    &fe_face_values,
                                   dealii::hp::FESubfaceValues<DoFHandlerType::dimension, DoFHandlerType::space_dimension> &fe_subface_values,
                                   const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      const unsigned int dim = DoFHandlerType::dimension;
      (void)dim;

      const typename DoFHandlerType::cell_iterator neighbor = cell->neighbor(face_no);
      (void)neighbor;
      const unsigned int n_solution_vectors = solutions.size();
      const typename DoFHandlerType::face_iterator
      face=cell->face(face_no);

      Assert (neighbor.state() == IteratorState::valid, ExcInternalError());
      Assert (face->has_children(), ExcInternalError());

      // set up a vector of the gradients of the finite element function on
      // this cell at the quadrature points
      //
      // let psi be a short name for [a grad u_h], where the second index be
      // the component of the finite element, and the first index the number
      // of the quadrature point

      // store which number @p{cell} has in the list of neighbors of
      // @p{neighbor}
      const unsigned int neighbor_neighbor
        = cell->neighbor_of_neighbor (face_no);
      Assert (neighbor_neighbor<GeometryInfo<dim>::faces_per_cell,
              ExcInternalError());

      // loop over all subfaces
      for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
        {
          // get an iterator pointing to the cell behind the present subface
          const typename DoFHandlerType::active_cell_iterator neighbor_child
            = cell->neighbor_child_on_subface (face_no, subface_no);
          Assert (!neighbor_child->has_children(),
                  ExcInternalError());

          // restrict the finite element on the present cell to the subface
          fe_subface_values.reinit (cell, face_no, subface_no,
                                    cell->active_fe_index());

          // restrict the finite element on the neighbor cell to the common
          // @p{subface}.
          fe_face_values.reinit (neighbor_child, neighbor_neighbor,
                                 cell->active_fe_index());

          const double factor = irregular_face_factor<DoFHandlerType>(cell,
                                                                      neighbor_child,
                                                                      face_no,
                                                                      subface_no,
                                                                      fe_face_values,
                                                                      fe_subface_values,
                                                                      strategy);

          // store the gradient of the solution in psi
          for (unsigned int n=0; n<n_solution_vectors; ++n)
            fe_subface_values.get_present_fe_values()
            .get_function_gradients (*solutions[n], parallel_data.psi[n]);

          // store the gradient from the neighbor's side in @p{neighbor_psi}
          for (unsigned int n=0; n<n_solution_vectors; ++n)
            fe_face_values.get_present_fe_values()
            .get_function_gradients (*solutions[n], parallel_data.neighbor_psi[n]);

          // call generic evaluate function
          parallel_data.neighbor_normal_vectors =
            fe_subface_values.get_present_fe_values().get_all_normal_vectors();

          local_face_integrals[neighbor_child->face(neighbor_neighbor)] =
            integrate_over_face (parallel_data, face, fe_face_values);
          for (unsigned int i = 0; i < local_face_integrals[neighbor_child->face(neighbor_neighbor)].size(); i++)
            local_face_integrals[neighbor_child->face(neighbor_neighbor)][i] *= factor;
        }

      // finally loop over all subfaces to collect the contributions of the
      // subfaces and store them with the mother face
      std::vector<double> sum (n_solution_vectors, 0);
      for (unsigned int subface_no=0; subface_no<face->n_children(); ++subface_no)
        {
          Assert (local_face_integrals.find(face->child(subface_no)) !=
                  local_face_integrals.end(),
                  ExcInternalError());
          Assert (local_face_integrals[face->child(subface_no)][0] >= 0,
                  ExcInternalError());

          for (unsigned int n=0; n<n_solution_vectors; ++n)
            sum[n] += local_face_integrals[face->child(subface_no)][n];
        }

      local_face_integrals[face] = sum;
    }


    /**
     * Computate the error on the faces of a single cell.
     *
     * This function is only needed in two or three dimensions.  The error
     * estimator in one dimension is implemented separately.
     */
    template <typename InputVector, typename DoFHandlerType>
    void
    estimate_one_cell (const typename DoFHandlerType::active_cell_iterator &cell,
                       ParallelData<DoFHandlerType, typename InputVector::value_type> &parallel_data,
                       std::map<typename DoFHandlerType::face_iterator,std::vector<double> > &local_face_integrals,
                       const std::vector<const InputVector *> &solutions,
                       const typename KellyErrorEstimator<DoFHandlerType::dimension,DoFHandlerType::space_dimension>::Strategy strategy)
    {
      const unsigned int dim = DoFHandlerType::dimension;
      const unsigned int n_solution_vectors = solutions.size();

      const types::subdomain_id subdomain_id = parallel_data.subdomain_id;
      const unsigned int material_id  = parallel_data.material_id;

      // empty our own copy of the local face integrals
      local_face_integrals.clear();

      // loop over all faces of this cell
      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          const typename DoFHandlerType::face_iterator
          face=cell->face(face_no);

          // make sure we do work only once: this face may either be regular
          // or irregular. if it is regular and has a neighbor, then we visit
          // the face twice, once from every side. let the one with the lower
          // index do the work. if it is at the boundary, or if the face is
          // irregular, then do the work below
          if ((face->has_children() == false) &&
              !cell->at_boundary(face_no) &&
              (!cell->neighbor_is_coarser(face_no) &&
               (cell->neighbor(face_no)->index() < cell->index() ||
                (cell->neighbor(face_no)->index() == cell->index() &&
                 cell->neighbor(face_no)->level() < cell->level()))))
            continue;

          // if the neighboring cell is less refined than the present one,
          // then do nothing since we integrate over the subfaces when we
          // visit the coarse cells.
          if (face->at_boundary() == false)
            if (cell->neighbor_is_coarser(face_no))
              continue;

          // if this face is part of the boundary but not of the neumann
          // boundary -> nothing to do. However, to make things easier when
          // summing up the contributions of the faces of cells, we enter this
          // face into the list of faces with contribution zero.
          if (face->at_boundary()
              &&
              (parallel_data.neumann_bc->find(face->boundary_id()) ==
               parallel_data.neumann_bc->end()))
            {
              local_face_integrals[face]
                = std::vector<double> (n_solution_vectors, 0.);
              continue;
            }

          // finally: note that we only have to do something if either the
          // present cell is on the subdomain we care for (and the same for
          // material_id), or if one of the neighbors behind the face is on
          // the subdomain we care for
          if ( ! ( ((subdomain_id == numbers::invalid_subdomain_id)
                    ||
                    (cell->subdomain_id() == subdomain_id))
                   &&
                   ((material_id == numbers::invalid_material_id)
                    ||
                    (cell->material_id() == material_id))) )
            {
              // ok, cell is unwanted, but maybe its neighbor behind the face
              // we presently work on? oh is there a face at all?
              if (face->at_boundary())
                continue;

              bool care_for_cell = false;
              if (face->has_children() == false)
                care_for_cell |= ((cell->neighbor(face_no)->subdomain_id()
                                   == subdomain_id) ||
                                  (subdomain_id == numbers::invalid_subdomain_id))
                                 &&
                                 ((cell->neighbor(face_no)->material_id()
                                   == material_id) ||
                                  (material_id == numbers::invalid_material_id));
              else
                {
                  for (unsigned int sf=0; sf<face->n_children(); ++sf)
                    if (((cell->neighbor_child_on_subface(face_no,sf)
                          ->subdomain_id() == subdomain_id)
                         &&
                         (material_id ==
                          numbers::invalid_material_id))
                        ||
                        ((cell->neighbor_child_on_subface(face_no,sf)
                          ->material_id() == material_id)
                         &&
                         (subdomain_id ==
                          numbers::invalid_subdomain_id)))
                      {
                        care_for_cell = true;
                        break;
                      }
                }

              // so if none of the neighbors cares for this subdomain or
              // material either, then try next face
              if (care_for_cell == false)
                continue;
            }

          // so now we know that we care for this face, let's do something
          // about it. first re-size the arrays we may use to the correct
          // size:
          parallel_data.resize (cell->active_fe_index());


          // then do the actual integration
          if (face->has_children() == false)
            // if the face is a regular one, i.e.  either on the other side
            // there is nirvana (face is at boundary), or the other side's
            // refinement level is the same as that of this side, then handle
            // the integration of these both cases together
            integrate_over_regular_face (solutions,
                                         parallel_data,
                                         local_face_integrals,
                                         cell, face_no,
                                         parallel_data.fe_face_values_cell,
                                         parallel_data.fe_face_values_neighbor,
                                         strategy);

          else
            // otherwise we need to do some special computations which do not
            // fit into the framework of the above function
            integrate_over_irregular_face (solutions,
                                           parallel_data,
                                           local_face_integrals,
                                           cell, face_no,
                                           parallel_data.fe_face_values_cell,
                                           parallel_data.fe_subface_values,
                                           strategy);
        }
    }
  }
}





// the following function is still independent of dimension, but it
// calls dimension dependent functions
template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim, spacedim>::
estimate (const Mapping<dim, spacedim>               &mapping,
          const DoFHandlerType                       &dof_handler,
          const Quadrature<dim-1>                    &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
            component_mask, coefficients, n_threads, subdomain_id, material_id, strategy);
}


template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim,spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const Quadrature<dim-1>                    &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
  estimate(StaticMappingQ1<dim,spacedim>::mapping, dof_handler, quadrature, neumann_bc, solution,
           error, component_mask, coefficients, n_threads,
           subdomain_id, material_id, strategy);
}


template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim, spacedim>::
estimate (const Mapping<dim, spacedim>               &mapping,
          const DoFHandlerType                       &dof_handler,
          const hp::QCollection<dim-1>               &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
  // just pass on to the other function
  const std::vector<const InputVector *> solutions (1, &solution);
  std::vector<Vector<float>*>              errors (1, &error);
  estimate (mapping, dof_handler, quadrature, neumann_bc, solutions, errors,
            component_mask, coefficients, n_threads, subdomain_id, material_id, strategy);
}


template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim, spacedim>::
estimate (const DoFHandlerType                       &dof_handler,
          const hp::QCollection<dim-1>               &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const InputVector                          &solution,
          Vector<float>                              &error,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
  estimate(StaticMappingQ1<dim, spacedim>::mapping, dof_handler, quadrature, neumann_bc, solution,
           error, component_mask, coefficients, n_threads,
           subdomain_id, material_id, strategy);
}




template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim, spacedim>::
estimate (const Mapping<dim, spacedim>               &mapping,
          const DoFHandlerType                       &dof_handler,
          const hp::QCollection<dim-1>               &face_quadratures,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const std::vector<const InputVector *>     &solutions,
          std::vector<Vector<float>*>                &errors,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int,
          const types::subdomain_id                   subdomain_id_,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
#ifdef DEAL_II_WITH_P4EST
  if (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
      (&dof_handler.get_triangulation())
      != 0)
    Assert ((subdomain_id_ == numbers::invalid_subdomain_id)
            ||
            (subdomain_id_ ==
             dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>&>
             (dof_handler.get_triangulation()).locally_owned_subdomain()),
            ExcMessage ("For parallel distributed triangulations, the only "
                        "valid subdomain_id that can be passed here is the "
                        "one that corresponds to the locally owned subdomain id."));

  const types::subdomain_id subdomain_id
    = ((dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
        (&dof_handler.get_triangulation())
        != 0)
       ?
       dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>&>
       (dof_handler.get_triangulation()).locally_owned_subdomain()
       :
       subdomain_id_);
#else
  const types::subdomain_id subdomain_id
    = subdomain_id_;
#endif

  const unsigned int n_components = dof_handler.get_fe().n_components();
  (void)n_components;

  // sanity checks
  Assert (solutions.size() > 0,
          ExcNoSolutions());
  Assert (solutions.size() == errors.size(),
          ExcIncompatibleNumberOfElements(solutions.size(), errors.size()));

  for (typename FunctionMap<spacedim>::type::const_iterator i=neumann_bc.begin();
       i!=neumann_bc.end(); ++i)
    Assert (i->second->n_components == n_components,
            ExcInvalidBoundaryFunction(i->first,
                                       i->second->n_components,
                                       n_components));

  Assert (component_mask.represents_n_components(n_components),
          ExcInvalidComponentMask());
  Assert (component_mask.n_selected_components(n_components) > 0,
          ExcInvalidComponentMask());

  Assert ((coefficients == 0) ||
          (coefficients->n_components == n_components) ||
          (coefficients->n_components == 1),
          ExcInvalidCoefficient());

  for (unsigned int n=0; n<solutions.size(); ++n)
    Assert (solutions[n]->size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(solutions[n]->size(),
                                 dof_handler.n_dofs()));

  const unsigned int n_solution_vectors = solutions.size();

  // Map of integrals indexed by the corresponding face. In this map we store
  // the integrated jump of the gradient for each face.  At the end of the
  // function, we again loop over the cells and collect the contributions of
  // the different faces of the cell.
  std::map<typename DoFHandlerType::face_iterator,std::vector<double> > face_integrals;

  // all the data needed in the error estimator by each of the threads is
  // gathered in the following structures
  const hp::MappingCollection<dim,spacedim> mapping_collection(mapping);
  const internal::ParallelData<DoFHandlerType,typename InputVector::value_type>
  parallel_data (dof_handler.get_fe(),
                 face_quadratures,
                 mapping_collection,
                 (!neumann_bc.empty() || (coefficients != 0)),
                 solutions.size(),
                 subdomain_id,
                 material_id,
                 &neumann_bc,
                 component_mask,
                 coefficients);
  std::map<typename DoFHandlerType::face_iterator,std::vector<double> > sample_local_face_integrals;

  // now let's work on all those cells:
  WorkStream::run (dof_handler.begin_active(),
                   static_cast<typename DoFHandlerType::active_cell_iterator>(dof_handler.end()),
                   std_cxx11::bind (&internal::estimate_one_cell<InputVector,DoFHandlerType>,
                                    std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::ref(solutions),strategy),
                   std_cxx11::bind (&internal::copy_local_to_global<DoFHandlerType>,
                                    std_cxx11::_1, std_cxx11::ref(face_integrals)),
                   parallel_data,
                   sample_local_face_integrals);

  // finally add up the contributions of the faces for each cell

  // reserve one slot for each cell and set it to zero
  for (unsigned int n=0; n<n_solution_vectors; ++n)
    {
      (*errors[n]).reinit (dof_handler.get_triangulation().n_active_cells());
      for (unsigned int i=0; i<dof_handler.get_triangulation().n_active_cells(); ++i)
        (*errors[n])(i)=0;
    }

  // now walk over all cells and collect information from the faces. only do
  // something if this is a cell we care for based on the subdomain id
  unsigned int present_cell=0;
  for (typename DoFHandlerType::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end();
       ++cell, ++present_cell)
    if ( ((subdomain_id == numbers::invalid_subdomain_id)
          ||
          (cell->subdomain_id() == subdomain_id))
         &&
         ((material_id == numbers::invalid_material_id)
          ||
          (cell->material_id() == material_id)))
      {
        // loop over all faces of this cell
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          {
            Assert(face_integrals.find(cell->face(face_no))
                   != face_integrals.end(),
                   ExcInternalError());
            const double factor = internal::cell_factor<DoFHandlerType>(cell,
                                                                        face_no,
                                                                        dof_handler,
                                                                        strategy);

            for (unsigned int n=0; n<n_solution_vectors; ++n)
              {
                // make sure that we have written a meaningful value into this
                // slot
                Assert (face_integrals[cell->face(face_no)][n] >= 0,
                        ExcInternalError());

                (*errors[n])(present_cell)
                += (face_integrals[cell->face(face_no)][n] * factor);
              }
          }

        for (unsigned int n=0; n<n_solution_vectors; ++n)
          (*errors[n])(present_cell) = std::sqrt((*errors[n])(present_cell));
      }
}



template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void
KellyErrorEstimator<dim, spacedim>::
estimate (const Mapping<dim, spacedim>               &mapping,
          const DoFHandlerType                       &dof_handler,
          const Quadrature<dim-1>                    &quadrature,
          const typename FunctionMap<spacedim>::type &neumann_bc,
          const std::vector<const InputVector *>     &solutions,
          std::vector<Vector<float>*>                &errors,
          const ComponentMask                        &component_mask,
          const Function<spacedim>                   *coefficients,
          const unsigned int                          n_threads,
          const types::subdomain_id                   subdomain_id,
          const types::material_id                    material_id,
          const Strategy                              strategy)
{
  // forward to the function with the QCollection
  estimate (mapping, dof_handler,
            hp::QCollection<dim-1>(quadrature),
            neumann_bc, solutions,
            errors, component_mask, coefficients,
            n_threads, subdomain_id, material_id, strategy);
}


template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void KellyErrorEstimator<dim, spacedim>::estimate
(const DoFHandlerType                       &dof_handler,
 const Quadrature<dim-1>                    &quadrature,
 const typename FunctionMap<spacedim>::type &neumann_bc,
 const std::vector<const InputVector *>     &solutions,
 std::vector<Vector<float>*>                &errors,
 const ComponentMask                        &component_mask,
 const Function<spacedim>                   *coefficients,
 const unsigned int                          n_threads,
 const types::subdomain_id                   subdomain_id,
 const types::material_id                    material_id,
 const Strategy                              strategy)
{
  estimate(StaticMappingQ1<dim, spacedim>::mapping, dof_handler, quadrature, neumann_bc, solutions,
           errors, component_mask, coefficients, n_threads,
           subdomain_id, material_id, strategy);
}



template <int dim, int spacedim>
template <typename InputVector, typename DoFHandlerType>
void KellyErrorEstimator<dim, spacedim>::estimate
(const DoFHandlerType                       &dof_handler,
 const hp::QCollection<dim-1>               &quadrature,
 const typename FunctionMap<spacedim>::type &neumann_bc,
 const std::vector<const InputVector *>     &solutions,
 std::vector<Vector<float>*>                &errors,
 const ComponentMask                        &component_mask,
 const Function<spacedim>                   *coefficients,
 const unsigned int                          n_threads,
 const types::subdomain_id                   subdomain_id,
 const types::material_id                    material_id,
 const Strategy                              strategy)
{
  estimate(StaticMappingQ1<dim, spacedim>::mapping, dof_handler, quadrature, neumann_bc, solutions,
           errors, component_mask, coefficients, n_threads,
           subdomain_id, material_id, strategy);
}

DEAL_II_NAMESPACE_CLOSE
