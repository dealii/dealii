// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2017 by the deal.II authors
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

#ifndef dealii__fe_tools_interpolate_templates_H
#define dealii__fe_tools_interpolate_templates_H


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/hp/dof_handler.h>


#include <deal.II/base/index_set.h>

#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace FETools
{
  template <int dim, int spacedim,
            template <int, int> class DoFHandlerType1,
            template <int, int> class DoFHandlerType2,
            class InVector, class OutVector>
  void
  interpolate(const DoFHandlerType1<dim, spacedim> &dof1,
              const InVector                       &u1,
              const DoFHandlerType2<dim, spacedim> &dof2,
              OutVector                            &u2)
  {
    ConstraintMatrix dummy;
    dummy.close();
    interpolate(dof1, u1, dof2, dummy, u2);
  }



  template <int dim, int spacedim,
            template <int, int> class DoFHandlerType1,
            template <int, int> class DoFHandlerType2,
            class InVector, class OutVector>
  void
  interpolate (const DoFHandlerType1<dim, spacedim> &dof1,
               const InVector                       &u1,
               const DoFHandlerType2<dim, spacedim> &dof2,
               const ConstraintMatrix               &constraints,
               OutVector                            &u2)
  {
    Assert(&dof1.get_triangulation()==&dof2.get_triangulation(), ExcTriangulationMismatch());

    Assert(u1.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));


    const IndexSet u2_elements = u2.locally_owned_elements();
#ifdef DEBUG
    const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
    const IndexSet &dof2_local_dofs = dof2.locally_owned_dofs();
    const IndexSet u1_elements = u1.locally_owned_elements();
    Assert(u1_elements == dof1_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
    Assert(u2_elements == dof2_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
#endif

    // allocate vectors at maximal
    // size. will be reinited in inner
    // cell, but Vector makes sure that
    // this does not lead to
    // reallocation of memory
    Vector<typename OutVector::value_type> u1_local(DoFTools::max_dofs_per_cell(dof1));
    Vector<typename OutVector::value_type> u2_local(DoFTools::max_dofs_per_cell(dof2));

    // have a map for interpolation
    // matrices. shared_ptr make sure
    // that memory is released again
    std::map<const FiniteElement<dim,spacedim> *,
        std::map<const FiniteElement<dim,spacedim> *,
        std::shared_ptr<FullMatrix<double> > > >
        interpolation_matrices;

    typename DoFHandlerType1<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active(),
                                                                 endc1 = dof1.end();
    typename DoFHandlerType2<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active(),
                                                                 endc2 = dof2.end();
    (void)endc2;

    std::vector<types::global_dof_index> dofs;
    dofs.reserve (DoFTools::max_dofs_per_cell (dof2));

    u2 = 0;
    OutVector touch_count(u2);
    touch_count = 0;

    // for distributed triangulations,
    // we can only interpolate u1 on
    // a cell, which this processor owns,
    // so we have to know the subdomain_id
    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    for (; cell1!=endc1; ++cell1, ++cell2)
      if ((cell1->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          Assert(cell1->get_fe().n_components() == cell2->get_fe().n_components(),
                 ExcDimensionMismatch (cell1->get_fe().n_components(),
                                       cell2->get_fe().n_components()));

          // for continuous elements on
          // grids with hanging nodes we
          // need hanging node
          // constraints. Consequently,
          // if there are no constraints
          // then hanging nodes are not
          // allowed.
          const bool hanging_nodes_not_allowed
            = ((cell2->get_fe().dofs_per_vertex != 0) &&
               (constraints.n_constraints() == 0));

          if (hanging_nodes_not_allowed)
            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              Assert (cell1->at_boundary(face) ||
                      cell1->neighbor(face)->level() == cell1->level(),
                      ExcHangingNodesNotAllowed(0));


          const unsigned int dofs_per_cell1 = cell1->get_fe().dofs_per_cell;
          const unsigned int dofs_per_cell2 = cell2->get_fe().dofs_per_cell;
          u1_local.reinit (dofs_per_cell1);
          u2_local.reinit (dofs_per_cell2);

          // check if interpolation
          // matrix for this particular
          // pair of elements is already
          // there
          if (interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()].get() == nullptr)
            {
              std::shared_ptr<FullMatrix<double> >
              interpolation_matrix (new FullMatrix<double> (dofs_per_cell2,
                                                            dofs_per_cell1));
              interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]
                = interpolation_matrix;

              get_interpolation_matrix(cell1->get_fe(),
                                       cell2->get_fe(),
                                       *interpolation_matrix);
            }

          cell1->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]
          ->vmult(u2_local, u1_local);

          dofs.resize (dofs_per_cell2);
          cell2->get_dof_indices(dofs);

          for (unsigned int i=0; i<dofs_per_cell2; ++i)
            {
              // if dof is locally_owned
              const types::global_dof_index gdi = dofs[i];
              if (u2_elements.is_element(gdi))
                {
                  ::dealii::internal::ElementAccess<OutVector>::add(u2_local(i),
                                                                    dofs[i], u2);
                  ::dealii::internal::ElementAccess<OutVector>::add(1,
                                                                    dofs[i], touch_count);
                }
            }
        }
    // cell1 is at the end, so should
    // be cell2
    Assert (cell2 == endc2, ExcInternalError());

    u2.compress(VectorOperation::add);
    touch_count.compress(VectorOperation::add);

    // if we work on parallel distributed
    // vectors, we have to ensure, that we only
    // work on dofs this processor owns.
    IndexSet  locally_owned_dofs = dof2.locally_owned_dofs();

    // when a discontinuous element is
    // interpolated to a continuous
    // one, we take the mean values.
    // for parallel vectors check,
    // if this component is owned by
    // this processor.
    for (types::global_dof_index i=0; i<dof2.n_dofs(); ++i)
      if (locally_owned_dofs.is_element(i))
        {
          Assert(static_cast<typename OutVector::value_type>(
                   ::dealii::internal::ElementAccess<OutVector>::get(
                     touch_count, i)) != typename OutVector::value_type(0),
                 ExcInternalError());


          const double val = ::dealii::internal::ElementAccess<OutVector>::get(
                               u2, i);
          ::dealii::internal::ElementAccess<OutVector>::set(
            val/::dealii::internal::ElementAccess<OutVector>::get(touch_count,i), i, u2);
        }

    // finish the work on parallel vectors
    u2.compress(VectorOperation::insert);
    // Apply hanging node constraints.
    constraints.distribute(u2);
  }



  template <int dim,
            template <int, int> class DoFHandlerType,
            class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DoFHandlerType<dim, spacedim> &dof1,
                   const InVector                      &u1,
                   const FiniteElement<dim,spacedim>   &fe2,
                   OutVector                           &u1_interpolated)
  {
    Assert(dof1.get_fe().n_components() == fe2.n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

#ifdef DEBUG
    const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
    const IndexSet u1_elements = u1.locally_owned_elements();
    const IndexSet u1_interpolated_elements = u1_interpolated.locally_owned_elements();
    Assert(u1_elements == dof1_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
    Assert(u1_interpolated_elements == dof1_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
#endif

    Vector<typename OutVector::value_type> u1_local(DoFTools::max_dofs_per_cell(dof1));
    Vector<typename OutVector::value_type> u1_int_local(DoFTools::max_dofs_per_cell(dof1));

    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    typename DoFHandlerType<dim, spacedim>::active_cell_iterator cell = dof1.begin_active(),
                                                                 endc = dof1.end();

    // map from possible fe objects in
    // dof1 to the back_interpolation
    // matrices
    std::map<const FiniteElement<dim> *,
        std::shared_ptr<FullMatrix<double> > > interpolation_matrices;

    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          // For continuous elements on
          // grids with hanging nodes we
          // need hanging node
          // constraints. Consequently,
          // when the elements are
          // continuous no hanging node
          // constraints are allowed.
          const bool hanging_nodes_not_allowed=
            (cell->get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

          if (hanging_nodes_not_allowed)
            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              Assert (cell->at_boundary(face) ||
                      cell->neighbor(face)->level() == cell->level(),
                      ExcHangingNodesNotAllowed(0));

          const unsigned int dofs_per_cell1 = cell->get_fe().dofs_per_cell;

          // make sure back_interpolation
          // matrix is available
          if (interpolation_matrices[&cell->get_fe()] == nullptr)
            {
              interpolation_matrices[&cell->get_fe()] =
                std::shared_ptr<FullMatrix<double> >
                (new FullMatrix<double>(dofs_per_cell1, dofs_per_cell1));
              get_back_interpolation_matrix(cell->get_fe(), fe2,
                                            *interpolation_matrices[&cell->get_fe()]);
            }

          u1_local.reinit (dofs_per_cell1);
          u1_int_local.reinit (dofs_per_cell1);

          cell->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell->get_fe()]->vmult(u1_int_local, u1_local);
          cell->set_dof_values(u1_int_local, u1_interpolated);
        };

    // if we work on a parallel vector, we have to finish the work
    u1_interpolated.compress(VectorOperation::insert);
  }



  namespace internal
  {
    namespace
    {
      template <int dim, int spacedim, class InVector>
      void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                             const ConstraintMatrix &constraints1,
                             const InVector &u1,
                             const DoFHandler<dim,spacedim> &dof2,
                             const ConstraintMatrix &constraints2,
                             InVector &u1_interpolated)
      {
        Vector<typename InVector::value_type> u2(dof2.n_dofs());
        interpolate(dof1, u1, dof2, constraints2, u2);
        interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
      }

      // special version for PETSc
#ifdef DEAL_II_WITH_PETSC
      template <int dim, int spacedim>
      void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                             const ConstraintMatrix &constraints1,
                             const PETScWrappers::MPI::Vector &u1,
                             const DoFHandler<dim,spacedim> &dof2,
                             const ConstraintMatrix &constraints2,
                             PETScWrappers::MPI::Vector &u1_interpolated)
      {
        // if u1 is a parallel distributed PETSc vector, we create a
        // vector u2 with based on the sets of locally owned and relevant
        // dofs of dof2
        IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet  dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        PETScWrappers::MPI::Vector  u2_out (dof2_locally_owned_dofs,
                                            u1.get_mpi_communicator());
        interpolate(dof1, u1, dof2, constraints2, u2_out);
        PETScWrappers::MPI::Vector  u2 (dof2_locally_owned_dofs,
                                        dof2_locally_relevant_dofs,
                                        u1.get_mpi_communicator());
        u2 = u2_out;
        interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
      }
#endif

      // special version for Trilinos
#ifdef DEAL_II_WITH_TRILINOS
      template <int dim, int spacedim>
      void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                             const ConstraintMatrix &constraints1,
                             const TrilinosWrappers::MPI::Vector &u1,
                             const DoFHandler<dim,spacedim> &dof2,
                             const ConstraintMatrix &constraints2,
                             TrilinosWrappers::MPI::Vector &u1_interpolated)
      {
        // if u1 is a parallel distributed Trilinos vector, we create a
        // vector u2 with based on the sets of locally owned and relevant
        // dofs of dof2
        IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet  dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        TrilinosWrappers::MPI::Vector  u2_out (dof2_locally_owned_dofs,
                                               u1.get_mpi_communicator());
        interpolate(dof1, u1, dof2, constraints2, u2_out);
        TrilinosWrappers::MPI::Vector  u2 (dof2_locally_owned_dofs,
                                           dof2_locally_relevant_dofs,
                                           u1.get_mpi_communicator());
        u2 = u2_out;
        interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
      }
#endif

      // special version for LinearAlgebra::distributed::Vector
      template <int dim, int spacedim, typename Number>
      void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                             const ConstraintMatrix &constraints1,
                             const LinearAlgebra::distributed::Vector<Number> &u1,
                             const DoFHandler<dim,spacedim> &dof2,
                             const ConstraintMatrix &constraints2,
                             LinearAlgebra::distributed::Vector<Number> &u1_interpolated)
      {
        IndexSet dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        LinearAlgebra::distributed::Vector<Number>
        u2 (dof2_locally_owned_dofs,
            dof2_locally_relevant_dofs,
            u1.get_mpi_communicator());

        interpolate(dof1, u1, dof2, constraints2, u2);
        u2.update_ghost_values ();
        interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
      }
    }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void back_interpolate(const DoFHandler<dim,spacedim> &dof1,
                        const ConstraintMatrix &constraints1,
                        const InVector &u1,
                        const DoFHandler<dim,spacedim> &dof2,
                        const ConstraintMatrix &constraints2,
                        OutVector &u1_interpolated)
  {
    // For discontinuous elements without constraints take the simpler version
    // of the back_interpolate function.
    if (dof1.get_fe().dofs_per_vertex==0 && dof2.get_fe().dofs_per_vertex==0
        && constraints1.n_constraints()==0 && constraints2.n_constraints()==0)
      back_interpolate(dof1, u1, dof2.get_fe(), u1_interpolated);
    else
      {
        Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
               ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
        Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
        Assert(u1_interpolated.size()==dof1.n_dofs(),
               ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

        // For continuous elements first interpolate to dof2, taking into
        // account constraints2, and then interpolate back to dof1 taking into
        // account constraints1
        internal::back_interpolate(dof1, constraints1, u1, dof2, constraints2,
                                   u1_interpolated);
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference (const DoFHandler<dim,spacedim> &dof1,
                                 const InVector &u1,
                                 const FiniteElement<dim,spacedim> &fe2,
                                 OutVector &u1_difference)
  {
    Assert(dof1.get_fe().n_components() == fe2.n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_difference.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));

#ifdef DEBUG
    const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
    const IndexSet u1_elements = u1.locally_owned_elements();
    const IndexSet u1_difference_elements = u1_difference.locally_owned_elements();
    Assert(u1_elements == dof1_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
    Assert(u1_difference_elements == dof1_local_dofs,
           ExcMessage("The provided vector and DoF handler should have the same"
                      " index sets."));
#endif

    // For continuous elements on grids
    // with hanging nodes we need
    // hanging node
    // constraints. Consequently, when
    // the elements are continuous no
    // hanging node constraints are
    // allowed.
    const bool hanging_nodes_not_allowed=
      (dof1.get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

    const unsigned int dofs_per_cell=dof1.get_fe().dofs_per_cell;

    Vector<typename OutVector::value_type> u1_local(dofs_per_cell);
    Vector<typename OutVector::value_type> u1_diff_local(dofs_per_cell);

    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    FullMatrix<double> difference_matrix(dofs_per_cell, dofs_per_cell);
    get_interpolation_difference_matrix(dof1.get_fe(), fe2,
                                        difference_matrix);

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof1.begin_active(),
                                                            endc = dof1.end();

    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          if (hanging_nodes_not_allowed)
            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              Assert (cell->at_boundary(face) ||
                      cell->neighbor(face)->level() == cell->level(),
                      ExcHangingNodesNotAllowed(0));

          cell->get_dof_values(u1, u1_local);
          difference_matrix.vmult(u1_diff_local, u1_local);
          cell->set_dof_values(u1_diff_local, u1_difference);
        }

    // if we work on a parallel PETSc vector
    // we have to finish the work and
    // update ghost values
    u1_difference.compress(VectorOperation::insert);
  }



  namespace internal
  {
    namespace
    {
      template <int dim, class InVector, class OutVector, int spacedim>
      void interpolation_difference (const DoFHandler<dim,spacedim> &dof1,
                                     const ConstraintMatrix &constraints1,
                                     const InVector &u1,
                                     const DoFHandler<dim,spacedim> &dof2,
                                     const ConstraintMatrix &constraints2,
                                     OutVector &u1_difference)
      {
        back_interpolate(dof1, constraints1, u1, dof2, constraints2, u1_difference);
        u1_difference.sadd(-1., 1., u1);
      }

      // special version for Trilinos
#ifdef DEAL_II_WITH_TRILINOS
      template <int dim, int spacedim>
      void interpolation_difference (const DoFHandler<dim,spacedim> &dof1,
                                     const ConstraintMatrix &constraints1,
                                     const TrilinosWrappers::MPI::Vector &u1,
                                     const DoFHandler<dim,spacedim> &dof2,
                                     const ConstraintMatrix &constraints2,
                                     TrilinosWrappers::MPI::Vector &u1_difference)
      {
        back_interpolate(dof1, constraints1, u1, dof2, constraints2, u1_difference);

        // Trilinos vectors with and without ghost entries are very different
        // and we cannot use the sadd function directly, so we have to create
        // a completely distributed vector first and copy the local entries
        // from the vector with ghost entries
        TrilinosWrappers::MPI::Vector  u1_completely_distributed (u1_difference.vector_partitioner ());

        u1_completely_distributed = u1;

        u1_difference.sadd(-1, u1_completely_distributed);
      }
#endif
    }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference(const DoFHandler<dim,spacedim> &dof1,
                                const ConstraintMatrix &constraints1,
                                const InVector &u1,
                                const DoFHandler<dim,spacedim> &dof2,
                                const ConstraintMatrix &constraints2,
                                OutVector &u1_difference)
  {
    // For discontinuous elements
    // without constraints take the
    // cheaper version of the
    // interpolation_difference function.
    if (dof1.get_fe().dofs_per_vertex==0 && dof2.get_fe().dofs_per_vertex==0
        && constraints1.n_constraints()==0 && constraints2.n_constraints()==0)
      interpolation_difference(dof1, u1, dof2.get_fe(), u1_difference);
    else
      {
        internal::interpolation_difference(dof1, constraints1, u1, dof2, constraints2, u1_difference);
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void project_dg(const DoFHandler<dim,spacedim> &dof1,
                  const InVector &u1,
                  const DoFHandler<dim,spacedim> &dof2,
                  OutVector &u2)
  {
    Assert(&dof1.get_triangulation()==&dof2.get_triangulation(), ExcTriangulationMismatch());
    Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator end = dof2.end();

    const unsigned int n1 = dof1.get_fe().dofs_per_cell;
    const unsigned int n2 = dof2.get_fe().dofs_per_cell;

    Vector<typename OutVector::value_type> u1_local(n1);
    Vector<typename OutVector::value_type> u2_local(n2);
    std::vector<types::global_dof_index> dofs(n2);

    FullMatrix<double> matrix(n2,n1);
    get_projection_matrix(dof1.get_fe(), dof2.get_fe(), matrix);

    u2 = 0;
    while (cell2 != end)
      {
        cell1->get_dof_values(u1, u1_local);
        matrix.vmult(u2_local, u1_local);
        cell2->get_dof_indices(dofs);
        for (unsigned int i=0; i<n2; ++i)
          {
            ::dealii::internal::ElementAccess<OutVector>::add(u2_local(i),
                                                              dofs[i], u2);
          }

        ++cell1;
        ++cell2;
      }
  }
} // end of namespace FETools

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   fe_tools_interpolate_templates.h     ---------------------------*/
/* end of #ifndef dealii__fe_tools_interpolate_templates_H */
#endif
