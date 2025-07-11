// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_tools_interpolate_templates_H
#define dealii_fe_tools_interpolate_templates_H


#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <iostream>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace FETools
{
  template <int dim, int spacedim, class InVector, class OutVector>
  void
  interpolate(const DoFHandler<dim, spacedim> &dof1,
              const InVector                  &u1,
              const DoFHandler<dim, spacedim> &dof2,
              OutVector                       &u2)
  {
    AffineConstraints<typename OutVector::value_type> dummy;
    dummy.close();
    interpolate(dof1, u1, dof2, dummy, u2);
  }



  template <int dim, int spacedim, class InVector, class OutVector>
  void
  interpolate(
    const DoFHandler<dim, spacedim>                         &dof1,
    const InVector                                          &u1,
    const DoFHandler<dim, spacedim>                         &dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints,
    OutVector                                               &u2)
  {
    Assert(&dof1.get_triangulation() == &dof2.get_triangulation(),
           ExcTriangulationMismatch());

    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size() == dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));


    const IndexSet u2_elements = u2.locally_owned_elements();
    if constexpr (running_in_debug_mode())
      {
        const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
        const IndexSet &dof2_local_dofs = dof2.locally_owned_dofs();
        const IndexSet  u1_elements     = u1.locally_owned_elements();
        Assert(u1_elements == dof1_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
        Assert(u2_elements == dof2_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
      }

    // allocate vectors at maximal
    // size. will be reinited in inner
    // cell, but Vector makes sure that
    // this does not lead to
    // reallocation of memory
    Vector<typename OutVector::value_type> u1_local(
      dof1.get_fe_collection().max_dofs_per_cell());
    Vector<typename OutVector::value_type> u2_local(
      dof2.get_fe_collection().max_dofs_per_cell());

    // have a map for interpolation matrices.
    // Using a unique_ptr makes sure that the
    // memory is released again automatically.
    std::map<const FiniteElement<dim, spacedim> *,
             std::map<const FiniteElement<dim, spacedim> *,
                      std::unique_ptr<FullMatrix<double>>>>
      interpolation_matrices;

    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell1 = dof1.begin_active(),
      endc1 = dof1.end();
    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell2 = dof2.begin_active(),
      endc2 = dof2.end();
    (void)endc2;

    std::vector<types::global_dof_index> dofs;
    dofs.reserve(dof2.get_fe_collection().max_dofs_per_cell());

    u2 = typename OutVector::value_type(0.);
    OutVector touch_count(u2);
    touch_count = typename OutVector::value_type(0.);

    // for distributed triangulations,
    // we can only interpolate u1 on
    // a cell, which this processor owns,
    // so we have to know the subdomain_id
    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    for (; cell1 != endc1; ++cell1, ++cell2)
      if ((cell1->subdomain_id() == subdomain_id) ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          Assert(cell1->get_fe().n_components() ==
                   cell2->get_fe().n_components(),
                 ExcDimensionMismatch(cell1->get_fe().n_components(),
                                      cell2->get_fe().n_components()));

          if constexpr (running_in_debug_mode())
            {
              // For continuous elements on grids with hanging nodes we need
              // hanging node constraints. Consequently, when the elements are
              // continuous no hanging node constraints are allowed.
              const bool hanging_nodes_not_allowed =
                ((cell2->get_fe().n_dofs_per_vertex() != 0) &&
                 (constraints.n_constraints() == 0));

              if (hanging_nodes_not_allowed)
                for (const unsigned int face : cell1->face_indices())
                  Assert(cell1->at_boundary(face) ||
                           cell1->neighbor(face)->level() == cell1->level(),
                         ExcHangingNodesNotAllowed());
            }

          const unsigned int dofs_per_cell1 = cell1->get_fe().n_dofs_per_cell();
          const unsigned int dofs_per_cell2 = cell2->get_fe().n_dofs_per_cell();
          u1_local.reinit(dofs_per_cell1);
          u2_local.reinit(dofs_per_cell2);

          // check if interpolation
          // matrix for this particular
          // pair of elements is already
          // there
          if (interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]
                .get() == nullptr)
            {
              auto interpolation_matrix =
                std::make_unique<FullMatrix<double>>(dofs_per_cell2,
                                                     dofs_per_cell1);

              get_interpolation_matrix(cell1->get_fe(),
                                       cell2->get_fe(),
                                       *interpolation_matrix);

              interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()] =
                std::move(interpolation_matrix);
            }

          cell1->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]->vmult(
            u2_local, u1_local);

          dofs.resize(dofs_per_cell2);
          cell2->get_dof_indices(dofs);

          for (unsigned int i = 0; i < dofs_per_cell2; ++i)
            {
              // if dof is locally_owned
              const types::global_dof_index gdi = dofs[i];
              if (u2_elements.is_element(gdi))
                {
                  ::dealii::internal::ElementAccess<OutVector>::add(u2_local(i),
                                                                    dofs[i],
                                                                    u2);
                  ::dealii::internal::ElementAccess<OutVector>::add(
                    1, dofs[i], touch_count);
                }
            }
        }
    // cell1 is at the end, so should
    // be cell2
    Assert(cell2 == endc2, ExcInternalError());

    u2.compress(VectorOperation::add);
    touch_count.compress(VectorOperation::add);

    // if we work on parallel distributed
    // vectors, we have to ensure, that we only
    // work on dofs this processor owns.
    const IndexSet &locally_owned_dofs = dof2.locally_owned_dofs();

    // when a discontinuous element is
    // interpolated to a continuous
    // one, we take the mean values.
    // for parallel vectors check,
    // if this component is owned by
    // this processor.
    for (types::global_dof_index i = 0; i < dof2.n_dofs(); ++i)
      if (locally_owned_dofs.is_element(i))
        {
          Assert(static_cast<typename OutVector::value_type>(
                   ::dealii::internal::ElementAccess<OutVector>::get(
                     touch_count, i)) != typename OutVector::value_type(0),
                 ExcInternalError());


          const typename OutVector::value_type val =
            ::dealii::internal::ElementAccess<OutVector>::get(u2, i);
          ::dealii::internal::ElementAccess<OutVector>::set(
            val /
              ::dealii::internal::ElementAccess<OutVector>::get(touch_count, i),
            i,
            u2);
        }

    // finish the work on parallel vectors
    u2.compress(VectorOperation::insert);
    // Apply hanging node constraints.
    constraints.distribute(u2);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DoFHandler<dim, spacedim>    &dof1,
                   const InVector                     &u1,
                   const FiniteElement<dim, spacedim> &fe2,
                   OutVector                          &u1_interpolated)
  {
    Assert(dof1.get_fe(0).n_components() == fe2.n_components(),
           ExcDimensionMismatch(dof1.get_fe(0).n_components(),
                                fe2.n_components()));
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

    if constexpr (running_in_debug_mode())
      {
        const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
        const IndexSet  u1_elements     = u1.locally_owned_elements();
        const IndexSet  u1_interpolated_elements =
          u1_interpolated.locally_owned_elements();
        Assert(u1_elements == dof1_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
        Assert(u1_interpolated_elements == dof1_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
      }

    Vector<typename OutVector::value_type> u1_local(
      dof1.get_fe_collection().max_dofs_per_cell());
    Vector<typename OutVector::value_type> u1_int_local(
      dof1.get_fe_collection().max_dofs_per_cell());

    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof1.begin_active(),
      endc = dof1.end();

    // map from possible FE objects in
    // dof1 to the back_interpolation
    // matrices
    std::map<const FiniteElement<dim> *, std::unique_ptr<FullMatrix<double>>>
      interpolation_matrices;

    for (; cell != endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id) ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          if constexpr (running_in_debug_mode())
            {
              // For continuous elements on grids with hanging nodes we need
              // hanging node constraints. Consequently, when the elements are
              // continuous no hanging node constraints are allowed.
              const bool hanging_nodes_not_allowed =
                (cell->get_fe().n_dofs_per_vertex() != 0) ||
                (fe2.n_dofs_per_vertex() != 0);

              if (hanging_nodes_not_allowed)
                for (const unsigned int face : cell->face_indices())
                  Assert(cell->at_boundary(face) ||
                           cell->neighbor(face)->level() == cell->level(),
                         ExcHangingNodesNotAllowed());
            }

          const unsigned int dofs_per_cell1 = cell->get_fe().n_dofs_per_cell();

          // make sure back_interpolation matrix is available
          if (interpolation_matrices[&cell->get_fe()] == nullptr)
            {
              interpolation_matrices[&cell->get_fe()] =
                std::make_unique<FullMatrix<double>>(dofs_per_cell1,
                                                     dofs_per_cell1);
              get_back_interpolation_matrix(
                cell->get_fe(), fe2, *interpolation_matrices[&cell->get_fe()]);
            }

          u1_local.reinit(dofs_per_cell1);
          u1_int_local.reinit(dofs_per_cell1);

          cell->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell->get_fe()]->vmult(u1_int_local,
                                                         u1_local);
          cell->set_dof_values(u1_int_local, u1_interpolated);
        };

    // if we work on a parallel vector, we have to finish the work
    u1_interpolated.compress(VectorOperation::insert);
  }



  namespace internal
  {
    template <int dim, int spacedim, class InVector>
    std::enable_if_t<is_serial_vector<InVector>::value>
    back_interpolate(
      const DoFHandler<dim, spacedim>                        &dof1,
      const AffineConstraints<typename InVector::value_type> &constraints1,
      const InVector                                         &u1,
      const DoFHandler<dim, spacedim>                        &dof2,
      const AffineConstraints<typename InVector::value_type> &constraints2,
      InVector                                               &u1_interpolated)
    {
      Vector<typename InVector::value_type> u2(dof2.n_dofs());
      interpolate(dof1, u1, dof2, constraints2, u2);
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }



    // special version for PETSc
#ifdef DEAL_II_WITH_PETSC
    template <int dim, int spacedim>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &dof1,
      const AffineConstraints<PETScWrappers::MPI::Vector::value_type>
                                       &constraints1,
      const PETScWrappers::MPI::Vector &u1,
      const DoFHandler<dim, spacedim>  &dof2,
      const AffineConstraints<PETScWrappers::MPI::Vector::value_type>
                                 &constraints2,
      PETScWrappers::MPI::Vector &u1_interpolated)
    {
      // if u1 is a parallel distributed PETSc vector, we create a
      // vector u2 with based on the sets of locally owned and relevant
      // dofs of dof2
      const IndexSet &dof2_locally_owned_dofs = dof2.locally_owned_dofs();
      const IndexSet  dof2_locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof2);

      PETScWrappers::MPI::Vector u2_out(dof2_locally_owned_dofs,
                                        u1.get_mpi_communicator());
      interpolate(dof1, u1, dof2, constraints2, u2_out);
      PETScWrappers::MPI::Vector u2(dof2_locally_owned_dofs,
                                    dof2_locally_relevant_dofs,
                                    u1.get_mpi_communicator());
      u2 = u2_out;
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }



    template <int dim, int spacedim>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<PETScWrappers::MPI::BlockVector::value_type> &,
      const PETScWrappers::MPI::BlockVector &,
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<PETScWrappers::MPI::BlockVector::value_type> &,
      PETScWrappers::MPI::BlockVector &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }
#endif



    // special version for Trilinos
#ifdef DEAL_II_WITH_TRILINOS
    template <int dim, int spacedim>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &dof1,
      const AffineConstraints<
        typename TrilinosWrappers::MPI::Vector::value_type> &constraints1,
      const TrilinosWrappers::MPI::Vector                   &u1,
      const DoFHandler<dim, spacedim>                       &dof2,
      const AffineConstraints<
        typename TrilinosWrappers::MPI::Vector::value_type> &constraints2,
      TrilinosWrappers::MPI::Vector                         &u1_interpolated)
    {
      // if u1 is a parallel distributed Trilinos vector, we create a
      // vector u2 with based on the sets of locally owned and relevant
      // dofs of dof2
      const IndexSet &dof2_locally_owned_dofs = dof2.locally_owned_dofs();
      const IndexSet  dof2_locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof2);

      TrilinosWrappers::MPI::Vector u2_out(dof2_locally_owned_dofs,
                                           u1.get_mpi_communicator());
      interpolate(dof1, u1, dof2, constraints2, u2_out);
      TrilinosWrappers::MPI::Vector u2(dof2_locally_owned_dofs,
                                       dof2_locally_relevant_dofs,
                                       u1.get_mpi_communicator());
      u2 = u2_out;
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }



    template <int dim, int spacedim>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &dof1,
      const AffineConstraints<
        typename TrilinosWrappers::MPI::BlockVector::value_type> &constraints1,
      const TrilinosWrappers::MPI::BlockVector                   &u1,
      const DoFHandler<dim, spacedim>                            &dof2,
      const AffineConstraints<
        typename TrilinosWrappers::MPI::BlockVector::value_type> &constraints2,
      TrilinosWrappers::MPI::BlockVector &u1_interpolated)
    {
      if (u1.n_blocks() == 0)
        return;
      const MPI_Comm  mpi_communicator = u1.block(0).get_mpi_communicator();
      const IndexSet &dof2_locally_owned_dofs = dof2.locally_owned_dofs();
      const IndexSet  dof2_locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof2);

      TrilinosWrappers::MPI::Vector u2_out(dof2_locally_owned_dofs,
                                           mpi_communicator);
      interpolate(dof1, u1, dof2, constraints2, u2_out);
      TrilinosWrappers::MPI::Vector u2(dof2_locally_owned_dofs,
                                       dof2_locally_relevant_dofs,
                                       mpi_communicator);
      u2 = u2_out;
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }



    template <int dim, int spacedim>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<
        typename LinearAlgebra::EpetraWrappers::Vector::value_type> &,
      const LinearAlgebra::EpetraWrappers::Vector &,
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<
        typename LinearAlgebra::EpetraWrappers::Vector::value_type> &,
      LinearAlgebra::EpetraWrappers::Vector &)
    {
      AssertThrow(false, ExcNotImplemented());
    }

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    template <int dim, int spacedim, typename Number, typename MemorySpace>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<typename LinearAlgebra::TpetraWrappers::
                                Vector<Number, MemorySpace>::value_type> &,
      const LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &,
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<typename LinearAlgebra::TpetraWrappers::
                                Vector<Number, MemorySpace>::value_type> &,
      LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &)
    {
      AssertThrow(false, ExcNotImplemented());
    }

    template <int dim, int spacedim, typename Number, typename MemorySpace>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<typename LinearAlgebra::TpetraWrappers::
                                BlockVector<Number, MemorySpace>::value_type> &,
      const LinearAlgebra::TpetraWrappers::BlockVector<Number, MemorySpace> &,
      const DoFHandler<dim, spacedim> &,
      const AffineConstraints<typename LinearAlgebra::TpetraWrappers::
                                BlockVector<Number, MemorySpace>::value_type> &,
      LinearAlgebra::TpetraWrappers::BlockVector<Number, MemorySpace> &)
    {
      AssertThrow(false, ExcNotImplemented());
    }
#  endif
#endif



    // special version for LinearAlgebra::distributed::Vector
    template <int dim, int spacedim, typename Number>
    void
    back_interpolate(
      const DoFHandler<dim, spacedim>                  &dof1,
      const AffineConstraints<Number>                  &constraints1,
      const LinearAlgebra::distributed::Vector<Number> &u1,
      const DoFHandler<dim, spacedim>                  &dof2,
      const AffineConstraints<Number>                  &constraints2,
      LinearAlgebra::distributed::Vector<Number>       &u1_interpolated)
    {
      const IndexSet &dof2_locally_owned_dofs = dof2.locally_owned_dofs();
      const IndexSet  dof2_locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof2);

      LinearAlgebra::distributed::Vector<Number> u2(dof2_locally_owned_dofs,
                                                    dof2_locally_relevant_dofs,
                                                    u1.get_mpi_communicator());

      interpolate(dof1, u1, dof2, constraints2, u2);
      u2.update_ghost_values();
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }



    // special version for LinearAlgebra::distributed::BlockVector
    template <int dim, int spacedim, typename Number>
    void
    back_interpolate(const DoFHandler<dim, spacedim> &,
                     const AffineConstraints<Number> &,
                     const LinearAlgebra::distributed::BlockVector<Number> &,
                     const DoFHandler<dim, spacedim> &,
                     const AffineConstraints<Number> &,
                     LinearAlgebra::distributed::BlockVector<Number> &)
    {
      AssertThrow(false, ExcNotImplemented());
    }
  } // namespace internal



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(
    const DoFHandler<dim, spacedim>                         &dof1,
    const AffineConstraints<typename OutVector::value_type> &constraints1,
    const InVector                                          &u1,
    const DoFHandler<dim, spacedim>                         &dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints2,
    OutVector                                               &u1_interpolated)
  {
    // For discontinuous elements without constraints take the simpler version
    // of the back_interpolate function.
    if (dof1.get_fe().n_dofs_per_vertex() == 0 &&
        dof2.get_fe().n_dofs_per_vertex() == 0 &&
        constraints1.n_constraints() == 0 && constraints2.n_constraints() == 0)
      back_interpolate(dof1, u1, dof2.get_fe(), u1_interpolated);
    else
      {
        Assert(dof1.get_fe(0).n_components() == dof2.get_fe(0).n_components(),
               ExcDimensionMismatch(dof1.get_fe(0).n_components(),
                                    dof2.get_fe(0).n_components()));
        Assert(u1.size() == dof1.n_dofs(),
               ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
        Assert(u1_interpolated.size() == dof1.n_dofs(),
               ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

        // For continuous elements first interpolate to dof2, taking into
        // account constraints2, and then interpolate back to dof1 taking into
        // account constraints1
        internal::back_interpolate(
          dof1, constraints1, u1, dof2, constraints2, u1_interpolated);
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolation_difference(const DoFHandler<dim, spacedim>    &dof1,
                           const InVector                     &u1,
                           const FiniteElement<dim, spacedim> &fe2,
                           OutVector                          &u1_difference)
  {
    Assert(dof1.get_fe(0).n_components() == fe2.n_components(),
           ExcDimensionMismatch(dof1.get_fe(0).n_components(),
                                fe2.n_components()));
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_difference.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));

    if constexpr (running_in_debug_mode())
      {
        const IndexSet &dof1_local_dofs = dof1.locally_owned_dofs();
        const IndexSet  u1_elements     = u1.locally_owned_elements();
        const IndexSet  u1_difference_elements =
          u1_difference.locally_owned_elements();
        Assert(u1_elements == dof1_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
        Assert(u1_difference_elements == dof1_local_dofs,
               ExcMessage(
                 "The provided vector and DoF handler should have the same"
                 " index sets."));
      }

    const unsigned int dofs_per_cell = dof1.get_fe().n_dofs_per_cell();

    Vector<typename OutVector::value_type> u1_local(dofs_per_cell);
    Vector<typename OutVector::value_type> u1_diff_local(dofs_per_cell);

    const types::subdomain_id subdomain_id =
      dof1.get_triangulation().locally_owned_subdomain();

    FullMatrix<double> difference_matrix(dofs_per_cell, dofs_per_cell);
    get_interpolation_difference_matrix(dof1.get_fe(), fe2, difference_matrix);

    typename DoFHandler<dim, spacedim>::active_cell_iterator
      cell = dof1.begin_active(),
      endc = dof1.end();

    for (; cell != endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id) ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          if constexpr (running_in_debug_mode())
            {
              // For continuous elements on grids with hanging nodes we need
              // hanging node constraints. Consequently, when the elements are
              // continuous no hanging node constraints are allowed.
              const bool hanging_nodes_not_allowed =
                (dof1.get_fe().n_dofs_per_vertex() != 0) ||
                (fe2.n_dofs_per_vertex() != 0);

              if (hanging_nodes_not_allowed)
                for (const unsigned int face : cell->face_indices())
                  Assert(cell->at_boundary(face) ||
                           cell->neighbor(face)->level() == cell->level(),
                         ExcHangingNodesNotAllowed());
            }

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
    template <int dim, class InVector, class OutVector, int spacedim>
    void
    interpolation_difference(
      const DoFHandler<dim, spacedim>                         &dof1,
      const AffineConstraints<typename OutVector::value_type> &constraints1,
      const InVector                                          &u1,
      const DoFHandler<dim, spacedim>                         &dof2,
      const AffineConstraints<typename OutVector::value_type> &constraints2,
      OutVector                                               &u1_difference)
    {
      back_interpolate(
        dof1, constraints1, u1, dof2, constraints2, u1_difference);
      u1_difference.sadd(-1., 1., u1);
    }

    // special version for Trilinos
#ifdef DEAL_II_WITH_TRILINOS
    template <int dim, int spacedim>
    void
    interpolation_difference(
      const DoFHandler<dim, spacedim> &dof1,
      const AffineConstraints<TrilinosWrappers::MPI::Vector::value_type>
                                          &constraints1,
      const TrilinosWrappers::MPI::Vector &u1,
      const DoFHandler<dim, spacedim>     &dof2,
      const AffineConstraints<TrilinosWrappers::MPI::Vector::value_type>
                                    &constraints2,
      TrilinosWrappers::MPI::Vector &u1_difference)
    {
      back_interpolate(
        dof1, constraints1, u1, dof2, constraints2, u1_difference);

      // Trilinos vectors with and without ghost entries are very different
      // and we cannot use the sadd function directly, so we have to create
      // a completely distributed vector first and copy the local entries
      // from the vector with ghost entries
      TrilinosWrappers::MPI::Vector u1_completely_distributed;
      u1_completely_distributed.reinit(u1_difference, true);
      u1_completely_distributed = u1;

      u1_difference.sadd(-1, u1_completely_distributed);
    }
#endif
  } // namespace internal



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolation_difference(
    const DoFHandler<dim, spacedim>                         &dof1,
    const AffineConstraints<typename OutVector::value_type> &constraints1,
    const InVector                                          &u1,
    const DoFHandler<dim, spacedim>                         &dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints2,
    OutVector                                               &u1_difference)
  {
    // For discontinuous elements
    // without constraints take the
    // cheaper version of the
    // interpolation_difference function.
    if (dof1.get_fe().n_dofs_per_vertex() == 0 &&
        dof2.get_fe().n_dofs_per_vertex() == 0 &&
        constraints1.n_constraints() == 0 && constraints2.n_constraints() == 0)
      interpolation_difference(dof1, u1, dof2.get_fe(), u1_difference);
    else
      {
        internal::interpolation_difference<dim, InVector, OutVector, spacedim>(
          dof1, constraints1, u1, dof2, constraints2, u1_difference);
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  project_dg(const DoFHandler<dim, spacedim> &dof1,
             const InVector                  &u1,
             const DoFHandler<dim, spacedim> &dof2,
             OutVector                       &u2)
  {
    Assert(&dof1.get_triangulation() == &dof2.get_triangulation(),
           ExcTriangulationMismatch());
    Assert(dof1.get_fe(0).n_components() == dof2.get_fe(0).n_components(),
           ExcDimensionMismatch(dof1.get_fe(0).n_components(),
                                dof2.get_fe(0).n_components()));
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size() == dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

    typename DoFHandler<dim, spacedim>::active_cell_iterator cell1 =
      dof1.begin_active();
    typename DoFHandler<dim, spacedim>::active_cell_iterator cell2 =
      dof2.begin_active();
    typename DoFHandler<dim, spacedim>::active_cell_iterator end = dof2.end();

    const unsigned int n1 = dof1.get_fe().n_dofs_per_cell();
    const unsigned int n2 = dof2.get_fe().n_dofs_per_cell();

    Vector<typename OutVector::value_type> u1_local(n1);
    Vector<typename OutVector::value_type> u2_local(n2);
    std::vector<types::global_dof_index>   dofs(n2);

    FullMatrix<double> matrix(n2, n1);
    get_projection_matrix(dof1.get_fe(), dof2.get_fe(), matrix);

    u2 = typename OutVector::value_type(0.);
    while (cell2 != end)
      {
        cell1->get_dof_values(u1, u1_local);
        matrix.vmult(u2_local, u1_local);
        cell2->get_dof_indices(dofs);
        for (unsigned int i = 0; i < n2; ++i)
          {
            ::dealii::internal::ElementAccess<OutVector>::add(u2_local(i),
                                                              dofs[i],
                                                              u2);
          }

        ++cell1;
        ++cell2;
      }
  }
} // end of namespace FETools

DEAL_II_NAMESPACE_CLOSE

#endif
