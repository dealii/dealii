// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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

#ifndef dealii_distributed_solution_transfer_h
#define dealii_distributed_solution_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    /**
     * Transfer a discrete FE function (like a solution vector) by
     * interpolation while refining and/or coarsening a distributed grid and
     * handles the necessary communication.
     *
     * @note It is important to note, that if you use more than one
     * SolutionTransfer object at the same time, that the calls to prepare_*()
     * and interpolate()/deserialize() need to be in the same order.
     *
     * <h3>Note on ghost elements</h3> In a parallel computation PETSc or
     * Trilinos vector may contain ghost elements or not. For reading in
     * information with prepare_for_coarsening_and_refinement() or
     * prepare_serialization() you need to supply vectors with ghost elements,
     * so that all locally_active elements can be read. On the other hand,
     * ghosted vectors are generally not writable, so for calls to
     * interpolate() or deserialize() you need to supply distributed vectors
     * without ghost elements. More precisely, during interpolation the
     * current algorithm writes into all locally active degrees of freedom.
     *
     * <h3>Transferring a solution</h3> Here VectorType is your favorite
     * vector type, e.g. PETScWrappers::MPI::Vector,
     * TrilinosWrappers::MPI::Vector, or corresponding block vectors.
     * @code
     * SolutionTransfer<dim, VectorType> soltrans(dof_handler);
     * // flag some cells for refinement and coarsening, e.g.
     * GridRefinement::refine_and_coarsen_fixed_fraction(tria,
     *                                                   error_indicators,
     *                                                   0.3,
     *                                                   0.05);
     *
     * // prepare the triangulation,
     * tria.prepare_coarsening_and_refinement();
     *
     * // prepare the SolutionTransfer object for coarsening and refinement
     * // and give the solution vector that we intend to interpolate later,
     * soltrans.prepare_for_coarsening_and_refinement(solution);
     *
     * // actually execute the refinement,
     * tria.execute_coarsening_and_refinement ();
     *
     * // redistribute dofs,
     * dof_handler.distribute_dofs (fe);
     *
     * // and interpolate the solution
     * VectorType interpolated_solution;
     *
     * //create VectorType in the right size here
     * soltrans.interpolate(interpolated_solution);
     * @endcode
     *
     * Different from PETSc and Trilinos vectors,
     * LinearAlgebra::distributed::Vector allows writing into ghost elements.
     * For a ghosted vector the interpolation step can be accomplished via
     * @code
     * interpolated_solution.zero_out_ghosts();
     * soltrans.interpolate(interpolated_solution);
     * interpolated_solution.update_ghost_values();
     * @endcode
     *
     * As the grid is distributed, it is important to note that the old
     * solution(s) must be copied to one that also provides access to the
     * locally relevant DoF values (these values required for the interpolation
     * process):
     * @code
     * // Create initial indexsets pertaining to the grid before refinement
     * IndexSet locally_owned_dofs, locally_relevant_dofs;
     * locally_owned_dofs = dof_handler.locally_owned_dofs();
     * DoFTools::extract_locally_relevant_dofs(dof_handler,
     * locally_relevant_dofs);
     *
     * // The solution vector only knows about locally owned DoFs
     * TrilinosWrappers::MPI::Vector solution;
     * solution.reinit(locally_owned_dofs,
     *                 mpi_communicator);
     * ...
     * // Transfer solution to vector that provides access to
     * // locally relevant DoFs
     * TrilinosWrappers::MPI::Vector old_solution;
     * old_solution.reinit(locally_owned_dofs,
     *                     locally_relevant_dofs,
     *                     mpi_communicator);
     * old_solution = solution;
     * ...
     * // Refine grid
     * // Recreate locally_owned_dofs and locally_relevant_dofs index sets
     * ...
     * solution.reinit(locally_owned_dofs, mpi_communicator);
     * soltrans.refine_interpolate(old_solution, solution);
     * @endcode
     *
     * <h3>Use for Serialization</h3>
     *
     * This class can be used to serialize and later deserialize a distributed
     * mesh with solution vectors to a file. If you use more than one
     * DoFHandler and therefore more than one SolutionTransfer object, they
     * need to be serialized and deserialized in the same order.
     *
     * If vector has the locally relevant DoFs, serialization works as
     * follows:
     * @code
     * parallel::distributed::SolutionTransfer<dim,VectorType>
     *   sol_trans(dof_handler);
     * sol_trans.prepare_serialization (vector);
     *
     * triangulation.save(filename);
     * @endcode
     * For deserialization the vector needs to be a distributed vector
     * (without ghost elements):
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * parallel::distributed::SolutionTransfer<dim,VectorType>
     *   sol_trans(dof_handler);
     * sol_trans.deserialize (distributed_vector);
     * @endcode
     *
     *
     * <h3>Note on usage with hp::DoFHandler</h3>
     *
     * If an object of the hp::DoFHandler class is registered with an
     * instantiation of this parallel::distributed::SolutionTransfer
     * class, the following requirements have to be met:
     * <ul>
     *   <li>
     *     The hp::DoFHandler needs to be explicitly mentioned
     *     in the parallel::distributed::SolutionTransfer type, i.e.:
     *     @code
     *     parallel::distributed::SolutionsTransfer<dim, VectorType,
     *       hp::DoFHandler<dim, spacedim>> sol_trans(hp_dof_handler);
     *     @endcode
     *   </li>
     *   <li>
     *     The transfer of the <tt>active_fe_index</tt> of each cell
     *     has to be scheduled as well in the parallel::distributed case,
     *     since ownerships of cells change during mesh repartitioning.
     *     This can be achieved with the
     *     parallel::distributed::ActiveFEIndicesTransfer class. The order in
     *     which both objects append data during the packing process does not
     *     matter. However, the unpacking process after refinement or
     *     deserialization requires the `active_fe_index` to be distributed
     *     before interpolating the data of the
     *     parallel::distributed::SolutionTransfer object, so that the correct
     *     FiniteElement object is associated with each cell. See the
     *     documentation on parallel::distributed::ActiveFEIndicesTransfer
     *     for further instructions.
     *   </li>
     * </ul>
     *
     * Since data on hp::DoFHandler objects is associated with many different
     * FiniteElement objects, each cell's data has to be processed with its
     * correpsonding `active_fe_index`. Further, if refinement is involved,
     * data will be packed on the parent cell with its `active_fe_index` and
     * unpacked later with the same index on its children. If cells get
     * coarsened into one, data will be packed on the latter with the least
     * dominating `active_fe_index` amongst its children, as determined by the
     * function hp::FECollection::find_least_face_dominating_fe_in_collection(),
     * and unpacked on the same cell with the same index.
     *
     * Code snippets to demonstrate the usage of the
     * parallel::distributed::SolutionTransfer class with hp::DoFHandler
     * objects are provided in the following. Here VectorType
     * is your favorite vector type, e.g. PETScWrappers::MPI::Vector,
     * TrilinosWrappers::MPI::Vector, or corresponding block vectors.
     *
     * After refinement, the order in which to unpack the transferred data
     * is important:
     * @code
     * //[prepare triangulation for refinement ...]
     *
     * parallel::distributed::ActiveFEIndicesTransfer<dim,spacedim>
     *   feidx_trans(hp_dof_handler);
     * parallel::distributed::
     *   SolutionTransfer<dim,VectorType,hp::DoFHandler<dim,spacedim>>
     *     sol_trans(hp_dof_handler);
     *
     * feidx_trans.prepare_for_transfer();
     * sol_trans.prepare_for_coarsening_and_refinement(solution);
     * triangulation.execute_coarsening_and_refinement();
     *
     * feidx_trans.unpack();
     * hp_dof_handler.distribute_dofs(fe_collection);
     *
     * VectorType interpolated_solution;
     * //[create VectorType in the right size here ...]
     * soltrans.interpolate(interpolated_solution);
     * @endcode
     *
     * If vector has the locally relevant DoFs, serialization works as follows:
     * @code
     * parallel::distributed::ActiveFEIndicesTransfer<dim,spacedim>
     *   feidx_trans(hp_dof_handler);
     * parallel::distributed::
     *   SolutionTransfer<dim, VectorType, hp::DoFHandler<dim,spacedim>>
     *     sol_trans(hp_dof_handler);
     *
     * feidx_trans.prepare_for_transfer();
     * sol_trans.prepare_serialization(vector);
     *
     * triangulation.save(filename);
     * @endcode
     *
     * For deserialization the vector needs to be a distributed vector
     * (without ghost elements):
     * @code
     * //[create coarse mesh...]
     * triangulation.load(filename);
     *
     * hp::DoFHandler<dim,spacedim>   hp_dof_handler(triangulation);
     * hp::FECollection<dim,spacedim> fe_collection;
     * //[prepare identical fe_collection...]
     * hp_dof_handler.distribute_dofs(fe_collection);
     *
     * parallel::distributed::ActiveFEIndicesTransfer<dim,spacedim>
     *   feidx_trans(hp_dof_handler);
     * feidx_trans.deserialize();
     * hp_dof_handler.distribute_dofs(fe_collection);
     *
     * parallel::distributed::
     *   SolutionTransfer<dim,VectorType,hp::DoFHandler<dim,spacedim>>
     *     sol_trans(dof_handler);
     * sol_trans.deserialize(distributed_vector);
     * @endcode
     *
     *
     * <h3>Interaction with hanging nodes</h3>
     *
     * In essence, this class implements the same steps as does
     * dealii::SolutionTransfer (though the implementation is entirely
     * separate). Consequently, the same issue with hanging nodes and
     * coarsening can happen with this class as happens with
     * dealii::SolutionTransfer. See there for an extended discussion.
     *
     * @ingroup distributed
     * @author Timo Heister, 2009-2011
     */
    template <int dim,
              typename VectorType,
              typename DoFHandlerType = DoFHandler<dim>>
    class SolutionTransfer
    {
#ifndef DEAL_II_MSVC
      static_assert(dim == DoFHandlerType::dimension,
                    "The dimension explicitly provided as a template "
                    "argument, and the dimension of the DoFHandlerType "
                    "template argument must match.");
#endif
    public:
      /**
       * Constructor.
       *
       * @param[in] dof The DoFHandler or hp::DoFHandler on which all
       *   operations will happen. At the time when this constructor
       *   is called, the DoFHandler still points to the triangulation
       *   before the refinement in question happens.
       */
      SolutionTransfer(const DoFHandlerType &dof);

      /**
       * Destructor.
       */
      ~SolutionTransfer() = default;

      /**
       * Prepare the current object for coarsening and refinement. It
       * stores the dof indices of each cell and stores the dof values of the
       * vectors in @p all_in in each cell that'll be coarsened. @p all_in
       * includes all vectors that are to be interpolated onto the new
       * (refined and/or coarsened) grid.
       */
      void
      prepare_for_coarsening_and_refinement(
        const std::vector<const VectorType *> &all_in);

      /**
       * Same as the previous function but for only one discrete function to be
       * interpolated.
       */
      void
      prepare_for_coarsening_and_refinement(const VectorType &in);

      /**
       * Interpolate the data previously stored in this object before the mesh
       * was refined or coarsened onto the current set of cells. Do so for
       * each of the vectors provided to
       * prepare_for_coarsening_and_refinement() and write the result into the
       * given set of vectors.
       */
      void
      interpolate(std::vector<VectorType *> &all_out);

      /**
       * Same as the previous function. It interpolates only one function. It
       * assumes the vectors having the right sizes (i.e.
       * <tt>in.size()==n_dofs_old</tt>, <tt>out.size()==n_dofs_refined</tt>)
       *
       * Multiple calling of this function is NOT allowed. Interpolating
       * several functions can be performed in one step by using
       * <tt>interpolate (all_in, all_out)</tt>
       */
      void
      interpolate(VectorType &out);


      /**
       * Prepare the serialization of the given vector. The serialization is
       * done by Triangulation::save(). The given vector needs all information
       * on the locally active DoFs (it must be ghosted). See documentation of
       * this class for more information.
       */
      void
      prepare_serialization(const VectorType &in);


      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      prepare_serialization(const std::vector<const VectorType *> &all_in);


      /**
       * Execute the deserialization of the given vector. This needs to be
       * done after calling Triangulation::load(). The given vector must be a
       * fully distributed vector without ghost elements. See documentation of
       * this class for more information.
       */
      void
      deserialize(VectorType &in);


      /**
       * Same as the function above, only for a list of vectors.
       */
      void
      deserialize(std::vector<VectorType *> &all_in);

    private:
      /**
       * Pointer to the degree of freedom handler to work with.
       */
      SmartPointer<const DoFHandlerType,
                   SolutionTransfer<dim, VectorType, DoFHandlerType>>
        dof_handler;

      /**
       * A vector that stores pointers to all the vectors we are supposed to
       * copy over from the old to the new mesh.
       */
      std::vector<const VectorType *> input_vectors;

      /**
       * The handle that the Triangulation has assigned to this object
       * with which we can access our memory offset and our pack function.
       */
      unsigned int handle;

      /**
       * A callback function used to pack the data on the current mesh into
       * objects that can later be retrieved after refinement, coarsening and
       * repartitioning.
       */
      std::vector<char>
      pack_callback(
        const typename Triangulation<dim, DoFHandlerType::space_dimension>::
          cell_iterator &cell,
        const typename Triangulation<dim, DoFHandlerType::space_dimension>::
          CellStatus status);

      /**
       * A callback function used to unpack the data on the current mesh that
       * has been packed up previously on the mesh before refinement,
       * coarsening and repartitioning.
       */
      void
      unpack_callback(
        const typename Triangulation<dim, DoFHandlerType::space_dimension>::
          cell_iterator &cell,
        const typename Triangulation<dim, DoFHandlerType::space_dimension>::
          CellStatus status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &                        data_range,
        std::vector<VectorType *> &all_out);


      /**
       * Registers the pack_callback() function to the
       * parallel::distributed::Triangulation that has been assigned to the
       * DoFHandler class member and stores the returning handle.
       */
      void
      register_data_attach();
    };


  } // namespace distributed
} // namespace parallel



DEAL_II_NAMESPACE_CLOSE

#endif
