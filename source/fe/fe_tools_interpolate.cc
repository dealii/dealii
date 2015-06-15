// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
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

#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <deal.II/base/index_set.h>

#ifdef DEAL_II_WITH_P4EST
#  include <p4est_bits.h>
#  include <p4est_extended.h>
#  include <p4est_vtk.h>
#  include <p4est_ghost.h>
#  include <p4est_communication.h>
#  include <p4est_iterate.h>

#  include <p8est_bits.h>
#  include <p8est_extended.h>
#  include <p8est_vtk.h>
#  include <p8est_ghost.h>
#  include <p8est_communication.h>
#  include <p8est_iterate.h>
#endif

#include <iostream>


DEAL_II_NAMESPACE_OPEN

namespace FETools
{
  template <int dim, int spacedim,
            template <int, int> class DH1,
            template <int, int> class DH2,
            class InVector, class OutVector>
  void
  interpolate(const DH1<dim, spacedim> &dof1,
              const InVector           &u1,
              const DH2<dim, spacedim> &dof2,
              OutVector                &u2)
  {
    ConstraintMatrix dummy;
    dummy.close();
    interpolate(dof1, u1, dof2, dummy, u2);
  }



  template <int dim, int spacedim,
            template <int, int> class DH1,
            template <int, int> class DH2,
            class InVector, class OutVector>
  void
  interpolate (const DH1<dim, spacedim> &dof1,
               const InVector           &u1,
               const DH2<dim, spacedim> &dof2,
               const ConstraintMatrix   &constraints,
               OutVector                &u2)
  {
    Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());

    Assert(u1.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

#ifdef DEAL_II_WITH_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          // if u1 is a parallel distributed
          // PETSc vector, we check the local
          // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.n_locally_owned_dofs()));
        }

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u2) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof2) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u2)->local_size() == dof2.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u2)->local_size(), dof2.n_locally_owned_dofs()));
        }
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
        std_cxx11::shared_ptr<FullMatrix<double> > > >
        interpolation_matrices;

    typename DH1<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active(),
                                                     endc1 = dof1.end();
    typename DH2<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active(),
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
      dof1.get_tria().locally_owned_subdomain();

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
          // constraints. Consequentely,
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
          if (interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()].get() == 0)
            {
              std_cxx11::shared_ptr<FullMatrix<double> >
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
              u2(dofs[i])+=u2_local(i);
              touch_count(dofs[i]) += 1;
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
          Assert(touch_count(i)!=0, ExcInternalError());
          u2(i) /= touch_count(i);
        }

    // finish the work on parallel vectors
    u2.compress(VectorOperation::insert);
    // Apply hanging node constraints.
    constraints.distribute(u2);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DoFHandler<dim,spacedim>    &dof1,
                   const InVector           &u1,
                   const FiniteElement<dim,spacedim> &fe2,
                   OutVector                &u1_interpolated)
  {
    Assert(dof1.get_fe().n_components() == fe2.n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

#ifdef DEAL_II_WITH_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          // if u1 is a parallel distributed
          // PETSc vector, we check the local
          // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.n_locally_owned_dofs()));
        }

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size(), dof1.n_locally_owned_dofs()));
        }
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

    const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

    Vector<typename OutVector::value_type> u1_local(dofs_per_cell1);
    Vector<typename OutVector::value_type> u1_int_local(dofs_per_cell1);

    const types::subdomain_id subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof1.begin_active(),
                                                            endc = dof1.end();

    FullMatrix<double> interpolation_matrix(dofs_per_cell1, dofs_per_cell1);
    get_back_interpolation_matrix(dof1.get_fe(), fe2,
                                  interpolation_matrix);
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
          interpolation_matrix.vmult(u1_int_local, u1_local);
          cell->set_dof_values(u1_int_local, u1_interpolated);
        }

    // if we work on a parallel PETSc vector
    // we have to finish the work
    u1_interpolated.compress(VectorOperation::insert);
  }



  template <int dim,
            template <int> class DH,
            class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DH<dim>            &dof1,
                   const InVector           &u1,
                   const FiniteElement<dim,spacedim> &fe2,
                   OutVector                &u1_interpolated)
  {
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

#ifdef DEAL_II_WITH_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          // if u1 is a parallel distributed
          // PETSc vector, we check the local
          // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.n_locally_owned_dofs()));
        }

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size(), dof1.n_locally_owned_dofs()));
        }
#endif

    Vector<typename OutVector::value_type> u1_local(DoFTools::max_dofs_per_cell(dof1));
    Vector<typename OutVector::value_type> u1_int_local(DoFTools::max_dofs_per_cell(dof1));

    const types::subdomain_id subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    typename DH<dim>::active_cell_iterator cell = dof1.begin_active(),
                                           endc = dof1.end();

    // map from possible fe objects in
    // dof1 to the back_interpolation
    // matrices
    std::map<const FiniteElement<dim> *,
        std_cxx11::shared_ptr<FullMatrix<double> > > interpolation_matrices;

    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == numbers::invalid_subdomain_id))
        {
          Assert(cell->get_fe().n_components() == fe2.n_components(),
                 ExcDimensionMismatch(cell->get_fe().n_components(),
                                      fe2.n_components()));

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
          if (interpolation_matrices[&cell->get_fe()] != 0)
            {
              interpolation_matrices[&cell->get_fe()] =
                std_cxx11::shared_ptr<FullMatrix<double> >
                (new FullMatrix<double>(dofs_per_cell1, dofs_per_cell1));
              get_back_interpolation_matrix(dof1.get_fe(), fe2,
                                            *interpolation_matrices[&cell->get_fe()]);
            }

          u1_local.reinit (dofs_per_cell1);
          u1_int_local.reinit (dofs_per_cell1);

          cell->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell->get_fe()]->vmult(u1_int_local, u1_local);
          cell->set_dof_values(u1_int_local, u1_interpolated);
        };

    // if we work on a parallel PETSc vector
    // we have to finish the work
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

      // special version for parallel::distributed::Vector
      template <int dim, int spacedim, typename Number>
      void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                             const ConstraintMatrix &constraints1,
                             const parallel::distributed::Vector<Number> &u1,
                             const DoFHandler<dim,spacedim> &dof2,
                             const ConstraintMatrix &constraints2,
                             parallel::distributed::Vector<Number> &u1_interpolated)
      {
        IndexSet dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        parallel::distributed::Vector<Number>
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

#ifdef DEAL_II_WITH_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          // if u1 is a parallel distributed
          // PETSc vector, we check the local
          // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.n_locally_owned_dofs()));
        }

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference)->local_size() == dof1.n_locally_owned_dofs(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference)->local_size(), dof1.n_locally_owned_dofs()));
        }
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
      dof1.get_tria().locally_owned_subdomain();

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
        u1_difference.sadd(-1, u1);
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
    Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
    Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator end = dof2.end();

    const unsigned int n1 = dof1.get_fe().dofs_per_cell;
    const unsigned int n2 = dof2.get_fe().dofs_per_cell;

    Vector<double> u1_local(n1);
    Vector<double> u2_local(n2);
    std::vector<types::global_dof_index> dofs(n2);

    FullMatrix<double> matrix(n2,n1);
    get_projection_matrix(dof1.get_fe(), dof2.get_fe(), matrix);

    while (cell2 != end)
      {
        cell1->get_dof_values(u1, u1_local);
        matrix.vmult(u2_local, u1_local);
        cell2->get_dof_indices(dofs);
        for (unsigned int i=0; i<n2; ++i)
          {
            u2(dofs[i])+=u2_local(i);
          }

        ++cell1;
        ++cell2;
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                   const InVector &u1,
                   const DoFHandler<dim,spacedim> &dof2,
                   OutVector &u2)
  {
    ConstraintMatrix dummy;
    dummy.close();
    extrapolate(dof1, u1, dof2, dummy, u2);
  }



  namespace internal
  {
    namespace p4est
    {
      /**
       * A structure whose explicit specializations contain pointers to the
       * relevant p4est_* and p8est_* functions. Using this structure, for
       * example by saying functions<dim>::quadrant_compare, we can write code
       * in a dimension independent way, either calling p4est_quadrant_compare
       * or p8est_quadrant_compare, depending on template argument.
       */
      template <int dim> struct functions;

      template <>
      struct functions<2>
      {
        static
        int (&quadrant_compare) (const void *v1, const void *v2);

        static
        int (&quadrant_overlaps_tree) (dealii::internal::p4est::types<2>::tree *tree,
                                       const dealii::internal::p4est::types<2>::quadrant *q);

        static
        int (&comm_find_owner) (dealii::internal::p4est::types<2>::forest *p4est,
                                const dealii::internal::p4est::types<2>::locidx which_tree,
                                const dealii::internal::p4est::types<2>::quadrant *q,
                                const int guess);
      };

      int (&functions<2>::quadrant_compare) (const void *v1, const void *v2)
        = p4est_quadrant_compare;

      int (&functions<2>::quadrant_overlaps_tree) (dealii::internal::p4est::types<2>::tree *tree,
                                                   const dealii::internal::p4est::types<2>::quadrant *q)
        = p4est_quadrant_overlaps_tree;

      int (&functions<2>::comm_find_owner) (dealii::internal::p4est::types<2>::forest *p4est,
                                            const dealii::internal::p4est::types<2>::locidx which_tree,
                                            const dealii::internal::p4est::types<2>::quadrant *q,
                                            const int guess)
        = p4est_comm_find_owner;

      template <>
      struct functions<3>
      {
        static
        int (&quadrant_compare) (const void *v1, const void *v2);

        static
        int (&quadrant_overlaps_tree) (dealii::internal::p4est::types<3>::tree *tree,
                                       const dealii::internal::p4est::types<3>::quadrant *q);

        static
        int (&comm_find_owner) (dealii::internal::p4est::types<3>::forest *p4est,
                                const dealii::internal::p4est::types<3>::locidx which_tree,
                                const dealii::internal::p4est::types<3>::quadrant *q,
                                const int guess);
      };

      int (&functions<3>::quadrant_compare) (const void *v1, const void *v2)
        = p8est_quadrant_compare;

      int (&functions<3>::quadrant_overlaps_tree) (dealii::internal::p4est::types<3>::tree *tree,
                                                   const dealii::internal::p4est::types<3>::quadrant *q)
        = p8est_quadrant_overlaps_tree;

      int (&functions<3>::comm_find_owner) (dealii::internal::p4est::types<3>::forest *p4est,
                                            const dealii::internal::p4est::types<3>::locidx which_tree,
                                            const dealii::internal::p4est::types<3>::quadrant *q,
                                            const int guess)
        = p8est_comm_find_owner;

      template <int dim, int spacedim>
      bool
      tree_exists_locally (const typename dealii::internal::p4est::types<dim>::forest *parallel_forest,
                           const typename dealii::internal::p4est::types<dim>::topidx coarse_grid_cell)
      {
        Assert (coarse_grid_cell < parallel_forest->connectivity->num_trees,
                ExcInternalError());
        return ((coarse_grid_cell >= parallel_forest->first_local_tree)
                &&
                (coarse_grid_cell <= parallel_forest->last_local_tree));
      }
    }

    /**
     * Implementation of the @p extrapolate function
     * on parallel distributed grids.
     */
    template <int dim,int spacedim>
    class ExtrapolateImplementation
    {
    public:

      ExtrapolateImplementation ();

      ~ExtrapolateImplementation ();

      template <class InVector,class OutVector>
      void extrapolate_parallel (const InVector &u2_relevant,
                                 const DoFHandler<dim,spacedim> &dof2,
                                 OutVector &u2);

    private:

      /**
       * A structure holding all data
       * of cells needed from other processes
       * for the extrapolate algorithm.
       */
      struct CellData
      {
        CellData ();

        CellData (const unsigned int dofs_per_cell);

        Vector<double>  dof_values;

        unsigned int tree_index;

        typename dealii::internal::p4est::types<dim>::quadrant  quadrant;

        int receiver;

        bool operator < (const CellData &rhs) const
        {
          if (p4est::functions<dim>::quadrant_compare (&quadrant, &rhs.quadrant) < 0)
            return true;

          return false;
        }

        unsigned int bytes_for_buffer () const
        {
          return (sizeof(unsigned int) +                                              // dofs_per_cell
                  dof_values.size() * sizeof(double) +                                // dof_values
                  sizeof(unsigned int) +                                              // tree_index
                  sizeof(typename dealii::internal::p4est::types<dim>::quadrant));    // quadrant
        }

        void pack_data (std::vector<char> &buffer) const
        {
          buffer.resize(bytes_for_buffer());

          char *ptr = &buffer[0];

          unsigned int n_dofs = dof_values.size ();
          std::memcpy(ptr, &n_dofs, sizeof(unsigned int));
          ptr += sizeof(unsigned int);

          std::memcpy(ptr,dof_values.begin(),n_dofs*sizeof(double));
          ptr += n_dofs*sizeof(double);

          std::memcpy(ptr,&tree_index,sizeof(unsigned int));
          ptr += sizeof(unsigned int);

          std::memcpy(ptr,&quadrant,sizeof(typename dealii::internal::p4est::types<dim>::quadrant));
          ptr += sizeof(typename dealii::internal::p4est::types<dim>::quadrant);

          Assert (ptr == &buffer[0]+buffer.size(),
                  ExcInternalError());
        }

        void unpack_data (const std::vector<char> &buffer)
        {
          const char *ptr = &buffer[0];
          unsigned int n_dofs;
          memcpy(&n_dofs, ptr, sizeof(unsigned int));
          ptr += sizeof(unsigned int);

          dof_values.reinit(n_dofs);
          std::memcpy(dof_values.begin(),ptr,n_dofs * sizeof(double));
          ptr += n_dofs * sizeof(double);

          std::memcpy(&tree_index,ptr,sizeof(unsigned int));
          ptr += sizeof(unsigned int);

          std::memcpy(&quadrant,ptr,sizeof(typename dealii::internal::p4est::types<dim>::quadrant));
          ptr += sizeof(typename dealii::internal::p4est::types<dim>::quadrant);

          Assert (ptr == &buffer[0]+buffer.size(),
                  ExcInternalError());
        }
      };

      // Problem: The function extrapolates a polynomial
      // function from a finer mesh of size $h$ to a polynmial
      // function of higher degree but on a coarser mesh of
      // size $2h$. Therefor the mesh has to consist of patches
      // of four (2d) or eight (3d) cells and the function requires
      // that the mesh is refined globally at least once.
      // The algorithm starts on the coarsest level of the grid,
      // loops over all cells and if a cell has at least one active child,
      // dof values are set via
      //   cell->get_interpolated_dof_values(u_input, dof_values)
      //   cell->set_dof_values_by_interpolation(dof_values, u_output)
      // both *_interpolation_* functions traverse recursively
      // over all children down to the active cells
      // and get/set dof values by interpolation.
      //
      // On distributed parallel grids the problem arises, that
      // if a cell has at least one active child, there is no guarantee
      // that all children of this cell belong to this process.
      // There might be children which are owned by other processes
      // and the algorithm needs to find and has to get the
      // dof values from these processes and so on.
      //
      // Algorithm:
      // 1) Loop over all coarse cells
      //      From each coarse cell traverse down the tree
      //      of refined cells and search for active children
      //      If there is an active child, check, if all
      //      other children down to the finest level are part
      //      of this process, if not, add the cell to the list
      //      of needs.
      // 2) Send the list of needs to other processes
      //      Each process has a list of needs and after
      //      the loop over all coarse cells (all trees)
      //      is finished, this list can be send to other processes.
      // 3) Compute the needs required by other processes
      //      While computing the values required from other
      //      processes there can arise new needs and they are
      //      stored in a list again.
      // 4) Send the computed values and the list of new needs around
      //
      // This procedure has to be repeated until no process needs any
      // new need and all needs are computed, but there are at most the
      // "number of grid levels" rounds of sending/receiving cell data.
      //
      // Then each process has all data needed for extrapolation.

      // driver function sending/receiving all values from
      // different processes
      template <class InVector>
      void compute_all_non_local_data (const DoFHandler<dim,spacedim> &dof2,
                                       const InVector                 &u);

      // traverse recursively over
      // the cells of this tree and
      // interpolate over patches which
      // are part of this process
      template <class InVector,class OutVector>
      void interpolate_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                    const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                    const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                    const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                    const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                    const InVector                                                &u1,
                                    OutVector                                                     &u2);

      // get dof values for this
      // cell by interpolation
      // if a child is reached which
      // is not part of this process
      // a new need is created
      template <class InVector,typename number>
      void get_interpolated_dof_values (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                        const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                        const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                        const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                        const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                        const InVector                                                &u,
                                        Vector<number>                                                &interpolated_values,
                                        std::vector<CellData>                                         &new_needs);

      // set dof values for this
      // cell by interpolation
      template <class OutVector,typename number>
      void set_dof_values_by_interpolation (const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                            const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                            const Vector<number>                                          &interpolated_values,
                                            OutVector                                                     &u);

      // compute all cell_data
      // needed from other processes
      // to interpolate the part of
      // this process
      void compute_needs (const DoFHandler<dim,spacedim> &dof2,
                          std::vector<CellData>          &new_needs);

      // traverse over the tree
      // and look for patches this
      // process has to work on
      void traverse_tree_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                      const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                      const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                      const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                      const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                      std::vector<CellData>                                         &new_needs);

      // traverse recursively
      // over a patch and look
      // for cells needed from
      // other processes for interpolation
      void traverse_patch_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                       const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                       const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                       const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                       const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                       std::vector<CellData>                                         &new_needs);

      // compute dof values of all
      // cells collected in cells_to_compute
      // computed cells are deleted
      // from cells_to_compute and
      // stored in computed_cells
      // during computation there can
      // be new cells needed from
      // other processes, these cells
      // are stored in new_needs
      template <class InVector>
      void compute_cells (const DoFHandler<dim,spacedim> &dof2,
                          const InVector                 &u,
                          std::vector<CellData>          &cells_to_compute,
                          std::vector<CellData>          &computed_cells,
                          std::vector<CellData>          &new_needs);

      // traverse over the tree
      // and compute cells
      template <class InVector>
      void compute_cells_in_tree_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                              const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                              const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                              const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                              const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                              const InVector                                                &u,
                                              std::vector<CellData>                                         &cells_to_compute,
                                              std::vector<CellData>                                         &computed_cells,
                                              std::vector<CellData>                                         &new_needs);

      // send all cell_data in the vector
      // cells_to_send to their receivers
      // and receives a vector of cell_data
      void send_cells (const std::vector<CellData>  &cells_to_send,
                       std::vector<CellData>  &received_cells) const;

      // add new cell_data to
      // the ordered list new_needs
      // uses cell_data_insert
      void add_new_need (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                         const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                         const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                         const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                         std::vector<CellData>                                         &new_needs) const;

      // binary search in cells_list
      // assume that cells_list
      // is ordered
      static int cell_data_search (const CellData               &cell_data,
                                   const std::vector<CellData>  &cells_list);

      // insert cell_data into a sorted
      // vector cells_list at the correct
      // position if cell_data
      // not exists already in cells_list
      static void cell_data_insert (const CellData         &cell_data,
                                    std::vector<CellData>  &cells_list);

      MPI_Comm  communicator;

      // a vector of all cells this process has
      // computed or received data
      std::vector<CellData>  available_cells;
    };

    template <int dim,int spacedim>
    ExtrapolateImplementation<dim,spacedim>::
    ExtrapolateImplementation ()
    {}



    template <int dim,int spacedim>
    ExtrapolateImplementation<dim,spacedim>::
    ~ExtrapolateImplementation ()
    {}



    template <int dim,int spacedim>
    ExtrapolateImplementation<dim,spacedim>::
    CellData::CellData ()
      : tree_index (0),
        receiver (0)
    {}



    template <int dim,int spacedim>
    ExtrapolateImplementation<dim,spacedim>::
    CellData::CellData (const unsigned int dofs_per_cell)
      : tree_index (0),
        receiver (0)
    {
      dof_values.reinit (dofs_per_cell);
    }



    template <int dim,int spacedim>
    template <class InVector,class OutVector>
    void
    ExtrapolateImplementation<dim,spacedim>::
    interpolate_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                             const typename dealii::internal::p4est::types<dim>::tree      &tree,
                             const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                             const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                             const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                             const InVector                                                &u1,
                             OutVector                                                     &u2)
    {
      // check if this cell exists in the local p4est
      int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                 &p4est_cell,
                                 p4est::functions<dim>::quadrant_compare);

      // this cell and none of it's children belongs to us
      if (idx == -1 && (p4est::functions<dim>::
                        quadrant_overlaps_tree (const_cast<typename dealii::internal::p4est::types<dim>::tree *>(&tree),
                                                &p4est_cell)
                        == false))
        return;

      bool p4est_has_children = (idx == -1);

      bool locally_owned_children = false;
      if (p4est_has_children)
        {
          Assert (dealii_cell->has_children (), ExcInternalError ());

          // check if at least one
          // child is locally owned
          // on our process
          for (unsigned int child_n=0; child_n<dealii_cell->n_children(); ++child_n)
            if (dealii_cell->child(child_n)->active())
              if (dealii_cell->child(child_n)->is_locally_owned())
                {
                  locally_owned_children=true;
                  break;
                }
        }

      if (locally_owned_children)
        {
          const FiniteElement<dim,spacedim> &fe            = dealii_cell->get_dof_handler().get_fe();
          const unsigned int                 dofs_per_cell = fe.dofs_per_cell;

          Vector<typename OutVector::value_type> interpolated_values(dofs_per_cell);

          std::vector<CellData>  new_needs;
          get_interpolated_dof_values (forest,
                                       tree,
                                       tree_index,
                                       dealii_cell,
                                       p4est_cell,
                                       u1,
                                       interpolated_values,
                                       new_needs);

          // at this point of
          // the procedure no new
          // needs should come up
          Assert (new_needs.size () == 0,
                  ExcInternalError ());

          set_dof_values_by_interpolation (dealii_cell,
                                           p4est_cell,
                                           interpolated_values,
                                           u2);
        }

      // traverse recursively over this tree
      if (p4est_has_children)
        {
          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              interpolate_recursively (forest,
                                       tree,
                                       tree_index,
                                       dealii_cell->child(c),
                                       p4est_child[c],
                                       u1,
                                       u2);
            }
        }
    }



    template <int dim,int spacedim>
    template <class InVector,typename number>
    void
    ExtrapolateImplementation<dim,spacedim>::
    get_interpolated_dof_values (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                 const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                 const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                 const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                 const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                 const InVector                                                &u,
                                 Vector<number>                                                &interpolated_values,
                                 std::vector<CellData>                                         &new_needs)
    {
      if (!dealii_cell->has_children ())
        {
          if (dealii_cell->is_locally_owned ())
            {
              // if this is one of our cells,
              // get dof values from input vector
              dealii_cell->get_dof_values (u, interpolated_values);
            }
          else
            {
              add_new_need (forest,
                            tree_index,
                            dealii_cell,
                            p4est_cell,
                            new_needs);
            }
        }
      else
        {
          const FiniteElement<dim,spacedim> &fe            = dealii_cell->get_dof_handler().get_fe();
          const unsigned int                 dofs_per_cell = fe.dofs_per_cell;

          Assert (&fe != 0,
                  ExcNotInitialized());
          Assert (interpolated_values.size() == dofs_per_cell,
                  ExcDimensionMismatch(interpolated_values.size(), dofs_per_cell));
          Assert (u.size() == dealii_cell->get_dof_handler().n_dofs(),
                  ExcDimensionMismatch(u.size(), dealii_cell->get_dof_handler().n_dofs()));

          Vector<number> tmp1(dofs_per_cell);
          Vector<number> tmp2(dofs_per_cell);

          interpolated_values = 0;
          std::vector<bool> restriction_is_additive (dofs_per_cell);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            restriction_is_additive[i] = fe.restriction_is_additive(i);

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          bool found_child = true;
          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              if (p4est::functions<dim>::
                  quadrant_overlaps_tree (const_cast<typename dealii::internal::p4est::types<dim>::tree *>(&tree),
                                          &p4est_child[c])
                  == false)
                {
                  // this is a cell this process needs
                  // data from another process

                  // check if this cell
                  // available from other
                  // computed patches
                  CellData  cell_data;
                  cell_data.quadrant = p4est_child[c];
                  int pos = cell_data_search (cell_data, available_cells);

                  if (pos == -1)
                    {
                      // data is not available
                      // create a new need
                      found_child = false;

                      add_new_need (forest,
                                    tree_index,
                                    dealii_cell->child(c),
                                    p4est_child[c],
                                    new_needs);
                    }
                  else
                    {
                      Assert (available_cells[pos].dof_values.size() == dofs_per_cell,
                              ExcDimensionMismatch(available_cells[pos].dof_values.size(), dofs_per_cell));

                      tmp1 = available_cells[pos].dof_values;
                    }
                }
              else
                {
                  // get the values from the present child, if necessary by
                  // interpolation itself either from its own children or
                  // by interpolating from the finite element on an active
                  // child to the finite element space requested here
                  get_interpolated_dof_values (forest,
                                               tree,
                                               tree_index,
                                               dealii_cell->child(c),
                                               p4est_child[c],
                                               u,
                                               tmp1,
                                               new_needs);
                }

              if (found_child)
                {
                  // interpolate these to the mother cell
                  fe.get_restriction_matrix(c, dealii_cell->refinement_case()).vmult (tmp2, tmp1);

                  // and add up or set them in the output vector
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    if (restriction_is_additive[i])
                      interpolated_values(i) += tmp2(i);
                    else if (tmp2(i) != number())
                      interpolated_values(i) = tmp2(i);
                }
            }

          if (found_child == false)
            interpolated_values = 0;
        }
    }



    template <int dim,int spacedim>
    template <class OutVector,typename number>
    void
    ExtrapolateImplementation<dim,spacedim>::
    set_dof_values_by_interpolation (const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                     const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                     const Vector<number>                                          &local_values,
                                     OutVector                                                     &u)
    {
      if (!dealii_cell->has_children ())
        {
          if (dealii_cell->is_locally_owned ())
            {
              // if this is one of our cells,
              // set dof values in output vector
              dealii_cell->set_dof_values (local_values, u);
            }
        }
      else
        {
          const FiniteElement<dim,spacedim> &fe            = dealii_cell->get_dof_handler().get_fe();
          const unsigned int                 dofs_per_cell = fe.dofs_per_cell;

          Assert (&fe != 0,
                  ExcNotInitialized());
          Assert (local_values.size() == dofs_per_cell,
                  ExcDimensionMismatch(local_values.size(), dofs_per_cell));
          Assert (u.size() == dealii_cell->get_dof_handler().n_dofs(),
                  ExcDimensionMismatch(u.size(), dealii_cell->get_dof_handler().n_dofs()));

          Vector<number> tmp(dofs_per_cell);

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              if (tmp.size() > 0)
                fe.get_prolongation_matrix(c, dealii_cell->refinement_case())
                .vmult (tmp, local_values);

              set_dof_values_by_interpolation (dealii_cell->child(c),
                                               p4est_child[c],
                                               tmp,
                                               u);
            }
        }
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    compute_needs (const DoFHandler<dim,spacedim> &dof2,
                   std::vector<CellData>          &new_needs)
    {
      const parallel::distributed::Triangulation< dim, spacedim > *tr
        = (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
           (&dof2.get_tria()));

      Assert (tr != 0, ExcInternalError());

      typename DoFHandler<dim,spacedim>::cell_iterator
      cell=dof2.begin(0),
      endc=dof2.end(0);
      for (; cell!=endc; ++cell)
        {
          if (p4est::tree_exists_locally<dim,spacedim> (tr->parallel_forest,
                                                        tr->coarse_cell_to_p4est_tree_permutation[cell->index()])
              == false)
            continue;

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          const unsigned int tree_index = tr->coarse_cell_to_p4est_tree_permutation[cell->index()];
          typename dealii::internal::p4est::types<dim>::tree *tree = tr->init_tree(cell->index());

          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          // make sure that each cell on the
          // coarsest level is at least once
          // refined, otherwise, these cells
          // can't be treated and would
          // generate a bogus result
          {
            int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree->quadrants),
                                       &p4est_coarse_cell,
                                       p4est::functions<dim>::quadrant_compare);

            AssertThrow (idx == -1, ExcGridNotRefinedAtLeastOnce ());
          }

          traverse_tree_recursively (*tr->parallel_forest,
                                     *tree,
                                     tree_index,
                                     cell,
                                     p4est_coarse_cell,
                                     new_needs);
        }
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    traverse_tree_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                               const typename dealii::internal::p4est::types<dim>::tree      &tree,
                               const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                               const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                               const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                               std::vector<CellData>                                         &new_needs)
    {
      // check if this cell exists in the local p4est
      int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                 &p4est_cell,
                                 p4est::functions<dim>::quadrant_compare);

      // this cell and none of it's children belongs to us
      if (idx == -1 && (p4est::functions<dim>::
                        quadrant_overlaps_tree (const_cast<typename dealii::internal::p4est::types<dim>::tree *>(&tree),
                                                &p4est_cell)
                        == false))
        return;

      bool p4est_has_children = (idx == -1);

      // this cell is part of a patch
      // this process has to interpolate on
      // if there is at least one locally
      // owned child
      bool locally_owned_children = false;
      if (p4est_has_children)
        {
          Assert (dealii_cell->has_children (), ExcInternalError ());

          for (unsigned int child_n=0; child_n<dealii_cell->n_children(); ++child_n)
            if (dealii_cell->child(child_n)->active())
              if (dealii_cell->child(child_n)->is_locally_owned())
                {
                  locally_owned_children=true;
                  break;
                }
        }

      // traverse the patch recursively
      // to find cells needed from
      // other processes
      if (locally_owned_children)
        {
          traverse_patch_recursively (forest,
                                      tree,
                                      tree_index,
                                      dealii_cell,
                                      p4est_cell,
                                      new_needs);
        }

      // traverse recursively over this tree
      if (p4est_has_children)
        {
          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              traverse_tree_recursively (forest,
                                         tree,
                                         tree_index,
                                         dealii_cell->child(c),
                                         p4est_child[c],
                                         new_needs);
            }
        }
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    traverse_patch_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                std::vector<CellData>                                         &new_needs)
    {
      if (dealii_cell->has_children ())
        {
          Assert (dealii_cell->has_children (), ExcInternalError ());

          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              // check if this child
              // is locally available
              // in the p4est
              if (p4est::functions<dim>::
                  quadrant_overlaps_tree (const_cast<typename dealii::internal::p4est::types<dim>::tree *>(&tree),
                                          &p4est_child[c])
                  == false)
                {
                  // this is a cell this process needs
                  // data from another process
                  // so add the cell to the list
                  add_new_need (forest,
                                tree_index,
                                dealii_cell->child(c),
                                p4est_child[c],
                                new_needs);
                }
              else
                {
                  // at least some part of
                  // the tree rooted in this
                  // child is locally available
                  traverse_patch_recursively (forest,
                                              tree,
                                              tree_index,
                                              dealii_cell->child(c),
                                              p4est_child[c],
                                              new_needs);
                }
            }
        }
    }




    template <int dim,int spacedim>
    template <class InVector>
    void
    ExtrapolateImplementation<dim,spacedim>::
    compute_cells (const DoFHandler<dim,spacedim> &dof2,
                   const InVector                 &u,
                   std::vector<CellData>          &cells_to_compute,
                   std::vector<CellData>          &computed_cells,
                   std::vector<CellData>          &new_needs)
    {
      const parallel::distributed::Triangulation< dim, spacedim > *tr
        = (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
           (&dof2.get_tria()));

      Assert (tr != 0, ExcInternalError());

      // collect in a set all trees this
      // process has to compute cells on
      std::set<unsigned int>  trees;
      for (typename std::vector<CellData>::const_iterator it=cells_to_compute.begin(); it!=cells_to_compute.end(); ++it)
        trees.insert (it->tree_index);

      typename DoFHandler<dim,spacedim>::cell_iterator
      cell=dof2.begin(0),
      endc=dof2.end(0);
      for (; cell!=endc; ++cell)
        {
          // check if this is a tree this process has to
          // work on and that this tree is in the p4est
          const unsigned int tree_index = tr->coarse_cell_to_p4est_tree_permutation[cell->index()];

          if ((trees.find (tree_index) == trees.end ())
              ||
              (p4est::tree_exists_locally<dim,spacedim> (tr->parallel_forest,
                                                         tree_index)
               == false))
            continue;

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          typename dealii::internal::p4est::types<dim>::tree *tree = tr->init_tree(cell->index());

          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          compute_cells_in_tree_recursively (*tr->parallel_forest,
                                             *tree,
                                             tree_index,
                                             cell,
                                             p4est_coarse_cell,
                                             u,
                                             cells_to_compute,
                                             computed_cells,
                                             new_needs);
        }
    }



    template <int dim,int spacedim>
    template <class InVector>
    void
    ExtrapolateImplementation<dim,spacedim>::
    compute_cells_in_tree_recursively (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                                       const typename dealii::internal::p4est::types<dim>::tree      &tree,
                                       const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                                       const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                                       const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                                       const InVector                                                &u,
                                       std::vector<CellData>                                         &cells_to_compute,
                                       std::vector<CellData>                                         &computed_cells,
                                       std::vector<CellData>                                         &new_needs)
    {
      if (cells_to_compute.size () == 0)
        return;

      // check if this cell exists in the local p4est
      int idx = sc_array_bsearch(const_cast<sc_array_t *>(&tree.quadrants),
                                 &p4est_cell,
                                 p4est::functions<dim>::quadrant_compare);

      // this cell and none of it's children belongs to us
      if (idx == -1 && (p4est::functions<dim>::
                        quadrant_overlaps_tree (const_cast<typename dealii::internal::p4est::types<dim>::tree *>(&tree),
                                                &p4est_cell)
                        == false))
        return;

      bool p4est_has_children = (idx == -1);

      // check if this quadrant is in the list
      CellData  cell_data;
      cell_data.quadrant = p4est_cell;
      int pos = cell_data_search (cell_data, cells_to_compute);
      if (pos != -1)
        {
          std::vector<CellData>  tmp;
          // compute dof values
          get_interpolated_dof_values (forest,
                                       tree,
                                       tree_index,
                                       dealii_cell,
                                       p4est_cell,
                                       u,
                                       cells_to_compute[pos].dof_values,
                                       tmp);
          // there is no new cell_data
          // this process needs for computing this cell,
          // store cell_data in the list of
          // computed cells and erase this cell
          // from the list of cells to compute
          if (tmp.size () == 0)
            {
              cell_data_insert (cells_to_compute[pos], computed_cells);
              cells_to_compute.erase (cells_to_compute.begin () + pos);
            }
          else
            {
              for (unsigned int i=0; i<tmp.size (); ++i)
                cell_data_insert (tmp[i], new_needs);
            }
        }

      // search recursively over this tree
      if (p4est_has_children)
        {
          typename dealii::internal::p4est::types<dim>::quadrant
          p4est_child[GeometryInfo<dim>::max_children_per_cell];

          dealii::internal::p4est::init_quadrant_children<dim> (p4est_cell,
                                                                p4est_child);

          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              compute_cells_in_tree_recursively (forest,
                                                 tree,
                                                 tree_index,
                                                 dealii_cell->child(c),
                                                 p4est_child[c],
                                                 u,
                                                 cells_to_compute,
                                                 computed_cells,
                                                 new_needs);
            }
        }
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    send_cells (const std::vector<CellData>  &cells_to_send,
                std::vector<CellData>        &received_cells) const
    {
      std::vector<std::vector<char> > sendbuffers (cells_to_send.size());
      std::vector<std::vector<char> >::iterator buffer = sendbuffers.begin();
      std::vector<MPI_Request> requests (cells_to_send.size());
      std::vector<unsigned int> destinations;

      // send data
      unsigned int idx=0;
      for (typename std::vector<CellData>::const_iterator it=cells_to_send.begin();
           it!=cells_to_send.end();
           ++it, ++idx)
        {
          destinations.push_back (it->receiver);

          it->pack_data (*buffer);
          MPI_Isend (&(*buffer)[0], buffer->size(),
                     MPI_BYTE,
                     it->receiver,
                     123,
                     communicator,
                     &requests[idx]);
        }

      Assert(destinations.size()==cells_to_send.size(), ExcInternalError());

      std::vector<unsigned int> senders
        = Utilities::MPI::compute_point_to_point_communication_pattern(communicator, destinations);

      // receive data
      std::vector<char> receive;
      CellData  cell_data;
      for (unsigned int index=0; index<senders.size (); ++index)
        {
          MPI_Status status;
          int len;
          MPI_Probe(MPI_ANY_SOURCE, 123, communicator, &status);
          MPI_Get_count(&status, MPI_BYTE, &len);
          receive.resize (len);

          char *buf = &receive[0];
          MPI_Recv (buf, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, communicator, &status);

          cell_data.unpack_data (receive);

          // this process has to send this
          // cell back to the sender
          // the receiver is the old sender
          cell_data.receiver = status.MPI_SOURCE;

          received_cells.push_back (cell_data);
        }

      if (requests.size () > 0)
        MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);

      // finally sort the list of cells
      std::sort (received_cells.begin (), received_cells.end ());
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    add_new_need (const typename dealii::internal::p4est::types<dim>::forest    &forest,
                  const typename dealii::internal::p4est::types<dim>::locidx    &tree_index,
                  const typename DoFHandler<dim,spacedim>::cell_iterator        &dealii_cell,
                  const typename dealii::internal::p4est::types<dim>::quadrant  &p4est_cell,
                  std::vector<CellData>                                         &new_needs) const
    {
      const FiniteElement<dim,spacedim> &fe            = dealii_cell->get_dof_handler().get_fe();
      const unsigned int                 dofs_per_cell = fe.dofs_per_cell;

      Assert (&fe != 0, ExcNotInitialized());

      CellData  cell_data (dofs_per_cell);
      cell_data.quadrant = p4est_cell;
      cell_data.tree_index = tree_index;
      cell_data.receiver = p4est::functions<dim>::
                           comm_find_owner (const_cast<typename dealii::internal::p4est::types<dim>::forest *> (&forest),
                                            tree_index,
                                            &p4est_cell,
                                            dealii_cell->level_subdomain_id ());

      cell_data_insert (cell_data, new_needs);
    }



    template <int dim,int spacedim>
    int
    ExtrapolateImplementation<dim,spacedim>::
    cell_data_search (const CellData               &cell_data,
                      const std::vector<CellData>  &cells_list)
    {
      typename std::vector<CellData>::const_iterator  bound
        = std::lower_bound (cells_list.begin(), cells_list.end(), cell_data);

      if ((bound != cells_list.end ())
          &&
          !(cell_data < *bound))
        return (int)(bound - cells_list.begin ());

      return (-1);
    }



    template <int dim,int spacedim>
    void
    ExtrapolateImplementation<dim,spacedim>::
    cell_data_insert (const CellData         &cell_data,
                      std::vector<CellData>  &cells_list)
    {
      typename std::vector<CellData>::iterator  bound
        = std::lower_bound (cells_list.begin(), cells_list.end(), cell_data);

      if ((bound == cells_list.end ())
          ||
          (cell_data < *bound))
        cells_list.insert (bound, 1, cell_data);
    }



    template <int dim,int spacedim>
    template <class InVector>
    void
    ExtrapolateImplementation<dim,spacedim>::
    compute_all_non_local_data (const DoFHandler<dim,spacedim> &dof2,
                                const InVector                 &u)
    {
      std::vector<CellData>  cells_we_need,
          cells_to_compute,
          received_cells,
          received_needs,
          new_needs,
          computed_cells,
          cells_to_send;

      // compute all the cells needed
      // from other processes
      compute_needs (dof2, cells_we_need);

      // send the cells needed to there
      // owners and receive a list other
      // processes need from us
      send_cells (cells_we_need, received_needs);

      // the list of received needs can contain
      // some cells more than ones because different
      // processes may need data from the same cell
      // to compute data only ones, generate a vector
      // with unique entries and distribute computed
      // data afterwards back to a vector with correct
      // receivers to send the data back
      // computing cell_data can cause some new cells
      // needed for this ones
      // if a cell is computed send it back to
      // their senders, maybe receive new needs and
      // compute again, do not wait that all cells
      // are computed or all needs are collected,
      // otherwise we can run into a deadlock if
      // a cell needed from another process,
      // itself needs some data from us
      unsigned int ready = 0;
      do
        {
          for (unsigned int i=0; i<received_needs.size(); ++i)
            cell_data_insert (received_needs[i], cells_to_compute);

          compute_cells (dof2, u, cells_to_compute, computed_cells, new_needs);

          // if there are no cells to compute and no new needs, stop
          ready = Utilities::MPI::sum (new_needs.size () + cells_to_compute.size (), communicator);

          for (typename std::vector<CellData>::const_iterator comp=computed_cells.begin ();
               comp != computed_cells.end ();
               ++comp)
            {
              // store computed cells
              cell_data_insert (*comp, available_cells);

              // and generate a vector
              // of computed cells with
              // correct receivers
              // then delete this received
              // need from the list
              typename std::vector<CellData>::iterator recv=std::begin (received_needs);
              while (recv != std::end (received_needs))
                {
                  if (dealii::internal::p4est::quadrant_is_equal<dim>(recv->quadrant, comp->quadrant))
                    {
                      recv->dof_values = comp->dof_values;
                      cells_to_send.push_back (*recv);
                      received_needs.erase (recv);
                      recv = std::begin (received_needs);
                    }
                  else
                    ++recv;
                }
            }

          send_cells (cells_to_send, received_cells);

          // strore received cell_data
          for (typename std::vector<CellData>::const_iterator recv=received_cells.begin ();
               recv != received_cells.end ();
               ++recv)
            {
              cell_data_insert (*recv, available_cells);
            }

          // finally send and receive new
          // needs and start a new round
          send_cells (new_needs, received_needs);
        }
      while (ready != 0);
    }



    template <int dim,int spacedim>
    template <class InVector,class OutVector>
    void
    ExtrapolateImplementation<dim,spacedim>::
    extrapolate_parallel (const InVector &u2_relevant,
                          const DoFHandler<dim,spacedim> &dof2,
                          OutVector &u2)
    {
      const parallel::distributed::Triangulation< dim, spacedim > *tr
        = (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
           (&dof2.get_tria()));

      Assert (tr != 0, ExcMessage ("Exrapolate in parallel only works for parallel distributed triangulations!"));

      communicator = tr->get_communicator ();

      compute_all_non_local_data (dof2, u2_relevant);

      // after all dof values on
      // not owned patch cells
      // are computed, start
      // the interpolation
      u2 = 0;

      typename DoFHandler<dim,spacedim>::cell_iterator
      cell=dof2.begin(0),
      endc=dof2.end(0);
      for (; cell!=endc; ++cell)
        {
          if (p4est::tree_exists_locally<dim,spacedim> (tr->parallel_forest,
                                                        tr->coarse_cell_to_p4est_tree_permutation[cell->index()])
              == false)
            continue;

          typename dealii::internal::p4est::types<dim>::quadrant p4est_coarse_cell;
          const unsigned int tree_index = tr->coarse_cell_to_p4est_tree_permutation[cell->index()];
          typename dealii::internal::p4est::types<dim>::tree *tree = tr->init_tree(cell->index());

          dealii::internal::p4est::init_coarse_quadrant<dim>(p4est_coarse_cell);

          interpolate_recursively (*tr->parallel_forest,
                                   *tree,
                                   tree_index,
                                   cell,
                                   p4est_coarse_cell,
                                   u2_relevant,
                                   u2);
        }

      u2.compress(VectorOperation::insert);
    }



    namespace
    {
      template <int dim, class InVector, class OutVector, int spacedim>
      void extrapolate_serial(const InVector &u3,
                              const DoFHandler<dim,spacedim> &dof2,
                              OutVector &u2)
      {
        const unsigned int dofs_per_cell  = dof2.get_fe().dofs_per_cell;
        Vector<typename OutVector::value_type> dof_values(dofs_per_cell);

        // then traverse grid bottom up
        for (unsigned int level=0; level<dof2.get_tria().n_levels()-1; ++level)
          {
            typename DoFHandler<dim,spacedim>::cell_iterator cell=dof2.begin(level),
                                                             endc=dof2.end(level);

            for (; cell!=endc; ++cell)
              if (!cell->active())
                {
                  // check whether this
                  // cell has active
                  // children
                  bool active_children=false;
                  for (unsigned int child_n=0; child_n<cell->n_children(); ++child_n)
                    if (cell->child(child_n)->active())
                      {
                        active_children=true;
                        break;
                      }

                  // if there are active
                  // children, this process
                  // has to work on this
                  // cell. get the data
                  // from the one vector
                  // and set it on the
                  // other
                  if (active_children)
                    {
                      cell->get_interpolated_dof_values(u3, dof_values);
                      cell->set_dof_values_by_interpolation(dof_values, u2);
                    }
                }
          }
      }



      template <int dim, class InVector, class OutVector, int spacedim>
      void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                       const InVector &u1,
                       const DoFHandler<dim,spacedim> &dof2,
                       const ConstraintMatrix &constraints,
                       OutVector &u2)
      {
        // make sure that each cell on the
        // coarsest level is at least once
        // refined, otherwise, these cells
        // can't be treated and would
        // generate a bogus result
        {
          typename DoFHandler<dim,spacedim>::cell_iterator cell = dof2.begin(0),
                                                           endc = dof2.end(0);
          for (; cell!=endc; ++cell)
            Assert (cell->has_children(), ExcGridNotRefinedAtLeastOnce());
        }

        OutVector u3;
        u3.reinit(u2);
        interpolate(dof1, u1, dof2, constraints, u3);

        extrapolate_serial (u3, dof2, u2);
      }



      template <int dim, class InVector, class OutVector, int spacedim>
      void extrapolate_parallel(const InVector &u2_relevant,
                                const DoFHandler<dim,spacedim> &dof2,
                                OutVector &u2)
      {
        if (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>(&dof2.get_tria()) != 0)
          {
            internal::ExtrapolateImplementation<dim,spacedim>  implementation;
            implementation.extrapolate_parallel (u2_relevant, dof2, u2);
          }
        else
          {
            extrapolate_serial (u2_relevant, dof2, u2);
          }
      }



      // special version for PETSc
#ifdef DEAL_II_WITH_PETSC
      template <int dim,int spacedim>
      void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                       const PETScWrappers::MPI::Vector &u1,
                       const DoFHandler<dim,spacedim> &dof2,
                       const ConstraintMatrix &constraints2,
                       PETScWrappers::MPI::Vector &u2)
      {
        IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet  dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        PETScWrappers::MPI::Vector u3 (dof2_locally_owned_dofs,
                                       u1.get_mpi_communicator ());
        interpolate (dof1, u1, dof2, constraints2, u3);
        PETScWrappers::MPI::Vector u3_relevant (dof2_locally_owned_dofs,
                                                dof2_locally_relevant_dofs,
                                                u1.get_mpi_communicator ());
        u3_relevant = u3;

        extrapolate_parallel (u3_relevant, dof2, u2);
      }
#endif



      // special version for Trilinos
#ifdef DEAL_II_WITH_TRILINOS
      template <int dim,int spacedim>
      void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                       const TrilinosWrappers::MPI::Vector &u1,
                       const DoFHandler<dim,spacedim> &dof2,
                       const ConstraintMatrix &constraints2,
                       TrilinosWrappers::MPI::Vector &u2)
      {
        IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet  dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        TrilinosWrappers::MPI::Vector u3 (dof2_locally_owned_dofs,
                                          u1.get_mpi_communicator ());
        interpolate (dof1, u1, dof2, constraints2, u3);
        TrilinosWrappers::MPI::Vector u3_relevant (dof2_locally_relevant_dofs,
                                                   u1.get_mpi_communicator ());
        u3_relevant = u3;

        extrapolate_parallel (u3_relevant, dof2, u2);
      }
#endif


      // special version for parallel::distributed::Vector
      template <int dim,int spacedim,typename Number>
      void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                       const parallel::distributed::Vector<Number> &u1,
                       const DoFHandler<dim,spacedim> &dof2,
                       const ConstraintMatrix &constraints2,
                       parallel::distributed::Vector<Number> &u2)
      {
        IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
        IndexSet  dof2_locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs (dof2,
                                                 dof2_locally_relevant_dofs);

        parallel::distributed::Vector<Number>  u3 (dof2_locally_owned_dofs,
                                                   dof2_locally_relevant_dofs,
                                                   u2.get_mpi_communicator ());

        interpolate (dof1, u1, dof2, constraints2, u3);
        u3.update_ghost_values ();

        extrapolate_parallel (u3, dof2, u2);
      }
    }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate(const DoFHandler<dim,spacedim> &dof1,
                   const InVector &u1,
                   const DoFHandler<dim,spacedim> &dof2,
                   const ConstraintMatrix &constraints,
                   OutVector &u2)
  {
    Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
           ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
    Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

    internal::extrapolate (dof1, u1, dof2, constraints, u2);

    // Apply hanging node constraints.
    constraints.distribute(u2);
  }

} // end of namespace FETools



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_tools_interpolate.inst"


/*----------------------------   fe_tools.cc     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE
