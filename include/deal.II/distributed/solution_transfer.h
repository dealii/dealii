// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2014 by the deal.II authors
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

#ifndef __deal2__distributed_solution_transfer_h
#define __deal2__distributed_solution_transfer_h

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
     * Transfers a discrete FE function (like a solution vector) by
     * interpolation while refining and/or
     * coarsening a distributed grid and
     * handles the necessary communication.
     *
     * @note It is important to note, that if you use more than one
     * SolutionTransfer object at the same time, that the calls to
     * prepare_*() and interpolate()/deserialize() need to be in the same order.
     *
     * <h3>Note on ghost elements</h3>
     * In a parallel computation PETSc or Trilinos vector may contain ghost
     * elements or not. For reading in information with
     * prepare_for_coarsening_and_refinement() or prepare_serialization() you need
     * to supply vectors with ghost elements, so that all locally_active elements
     * can be read. On the other hand, ghosted vectors are generally not writable,
     * so for calls to interpolate() or deserialize() you need to supply
     * distributed vectors without ghost elements.
     *
     * <h3>Transfering a solution</h3>
     * Here VECTOR is your favorite vector type, e.g. PETScWrappers::MPI::Vector,
     * TrilinosWrappers::MPI::Vector, or corresponding blockvectors.
     * @code
     SolutionTransfer<dim, VECTOR> soltrans(dof_handler);
                                         // flag some cells for refinement
                                         // and coarsening, e.g.
     GridRefinement::refine_and_coarsen_fixed_fraction(
       tria, error_indicators, 0.3, 0.05);
                                         // prepare the triangulation,
     tria.prepare_coarsening_and_refinement();
                                         // prepare the SolutionTransfer object
                                         // for coarsening and refinement and give
                                         // the solution vector that we intend to
                                         // interpolate later,
     soltrans.prepare_for_coarsening_and_refinement(solution);
                                         // actually execute the refinement,
     tria.execute_coarsening_and_refinement ();
                                         // redistribute dofs,
     dof_handler.distribute_dofs (fe);
                                         // and interpolate the solution
     VECTOR interpolated_solution;
     //create VECTOR in the right size here
     soltrans.interpolate(interpolated_solution);
    @endcode
     *
     * <h3>Use for Serialization</h3>
     *
     * This class can be used to serialize and later deserialize a distributed mesh with solution vectors
     * to a file. If you use more than one DoFHandler and therefore more than one SolutionTransfer object, they need to be serialized and deserialized in the same order.
     *
     * If vector has the locally relevant DoFs, serialization works as follows:
     *@code

    parallel::distributed::SolutionTransfer<dim,VECTOR> sol_trans(dof_handler);
    sol_trans.prepare_serialization (vector);

    triangulation.save(filename);
    @endcode
    * For deserialization the vector needs to be a distributed vector (without ghost elements):
    * @code
    //[create coarse mesh...]
    triangulation.load(filename);

    parallel::distributed::SolutionTransfer<dim,VECTOR> sol_trans(dof_handler);
    sol_trans.deserialize (distributed_vector);
    @endcode
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
    template<int dim, typename VECTOR, class DH=DoFHandler<dim> >
    class SolutionTransfer
    {
    public:
      /**
       * Constructor, takes the current DoFHandler
       * as argument.
       */
      SolutionTransfer(const DH &dof);
      /**
       * Destructor.
       */
      ~SolutionTransfer();

      /**
       * Prepares the @p SolutionTransfer for
       * coarsening and refinement. It
       * stores the dof indices of each cell and
       * stores the dof values of the vectors in
       * @p all_in in each cell that'll be coarsened.
       * @p all_in includes all vectors
       * that are to be interpolated
       * onto the new (refined and/or
       * coarsenend) grid.
       */
      void prepare_for_coarsening_and_refinement (const std::vector<const VECTOR *> &all_in);

      /**
       * Same as previous function
       * but for only one discrete function
       * to be interpolated.
       */
      void prepare_for_coarsening_and_refinement (const VECTOR &in);

      /**
       * Interpolate the data previously stored in this object before
       * the mesh was refined or coarsened onto the current set of
       * cells. Do so for each of the vectors provided to
       * prepare_for_coarsening_and_refinement() and write the result into
       * the given set of vectors.
       */
      void interpolate (std::vector<VECTOR *> &all_out);

      /**
       * Same as the previous function.
       * It interpolates only one function.
       * It assumes the vectors having the
       * right sizes (i.e. <tt>in.size()==n_dofs_old</tt>,
       * <tt>out.size()==n_dofs_refined</tt>)
       *
       * Multiple calling of this function is
       * NOT allowed. Interpolating
       * several functions can be performed
       * in one step by using
       * <tt>interpolate (all_in, all_out)</tt>
       */
      void interpolate (VECTOR &out);


      /**
       * Return the size in bytes that need
       * to be stored per cell.
       */
      unsigned int get_data_size() const;


      /**
      * Prepare the serialization of the
      * given vector. The serialization is
      * done by Triangulation::save(). The
      * given vector needs all information
      * on the locally active DoFs (it must
      * be ghosted). See documentation of
      * this class for more information.
      */
      void prepare_serialization(const VECTOR &in);


      /**
       * Same as the function above, only
       * for a list of vectors.
       */
      void prepare_serialization(const std::vector<const VECTOR *> &all_in);


      /**
       * Execute the deserialization of the
       * given vector. This needs to be
       * done after calling
       * Triangulation::load(). The given
       * vector must be a fully distributed
       * vector without ghost elements. See
       * documentation of this class for
       * more information.
       */
      void deserialize(VECTOR &in);


      /**
       * Same as the function above, only
       * for a list of vectors.
       */
      void deserialize(std::vector<VECTOR *> &all_in);

    private:
      /**
       * Pointer to the degree of
       * freedom handler to work
       * with.
       */
      SmartPointer<const DH,SolutionTransfer<dim,VECTOR,DH> > dof_handler;

      /**
       * A vector that stores
       * pointers to all the
       * vectors we are supposed to
       * copy over from the old to
       * the new mesh.
       */
      std::vector<const VECTOR *> input_vectors;

      /**
       * The offset that the
       * Triangulation has assigned
       * to this object starting at
       * which we are allowed to
       * write.
       */
      unsigned int offset;

      /**
       * A callback function used
       * to pack the data on the
       * current mesh into objects
       * that can later be
       * retrieved after
       * refinement, coarsening and
       * repartitioning.
       */
      void pack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell,
                         const typename Triangulation<dim,dim>::CellStatus status,
                         void *data);

      /**
       * A callback function used
       * to unpack the data on the
       * current mesh that has been
       * packed up previously on
       * the mesh before
       * refinement, coarsening and
       * repartitioning.
       */
      void unpack_callback(const typename Triangulation<dim,dim>::cell_iterator &cell,
                           const typename Triangulation<dim,dim>::CellStatus status,
                           const void *data,
                           std::vector<VECTOR *> &all_out);


      /**
       *
       */
      void register_data_attach(const std::size_t size);

    };


  }
}



DEAL_II_NAMESPACE_CLOSE

#endif
