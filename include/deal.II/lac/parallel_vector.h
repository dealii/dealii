// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2016 by the deal.II authors
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

#ifndef dealii__parallel_vector_h
#define dealii__parallel_vector_h

#include <deal.II/base/config.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <cstring>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace distributed
  {
    /*! @addtogroup Vectors
     *@{
     */


    /**
     * Implementation of a parallel vector class. The design of this class is
     * similar to the standard ::dealii::Vector class in deal.II, with the
     * exception that storage is distributed with MPI.
     *
     * The vector is designed for the following scheme of parallel
     * partitioning:
     * <ul>
     * <li> The indices held by individual processes (locally owned part) in
     * the MPI parallelization form a contiguous range
     * <code>[my_first_index,my_last_index)</code>.
     * <li> Ghost indices residing on arbitrary positions of other processors
     * are allowed. It is in general more efficient if ghost indices are
     * clustered, since they are stored as a set of intervals. The
     * communication pattern of the ghost indices is determined when calling
     * the function <code>reinit (locally_owned, ghost_indices,
     * communicator)</code>, and retained until the partitioning is changed.
     * This allows for efficient parallel communication of indices. In
     * particular, it stores the communication pattern, rather than having to
     * compute it again for every communication. For more information on ghost
     * vectors, see also the
     * @ref GlossGhostedVector "glossary entry on vectors with ghost elements".
     * <li> Besides the usual global access operator() it is also possible to
     * access vector entries in the local index space with the function @p
     * local_element(). Locally owned indices are placed first, [0,
     * local_size()), and then all ghost indices follow after them
     * contiguously, [local_size(), local_size()+n_ghost_entries()).
     * </ul>
     *
     * Functions related to parallel functionality:
     * <ul>
     * <li> The function <code>compress()</code> goes through the data
     * associated with ghost indices and communicates it to the owner process,
     * which can then add it to the correct position. This can be used e.g.
     * after having run an assembly routine involving ghosts that fill this
     * vector. Note that the @p insert mode of @p compress() does not set the
     * elements included in ghost entries but simply discards them, assuming
     * that the owning processor has set them to the desired value already
     * (See also the
     * @ref GlossCompress "glossary entry on compress").
     * <li> The <code>update_ghost_values()</code> function imports the data
     * from the owning processor to the ghost indices in order to provide read
     * access to the data associated with ghosts.
     * <li> It is possible to split the above functions into two phases, where
     * the first initiates the communication and the second one finishes it.
     * These functions can be used to overlap communication with computations
     * in other parts of the code.
     * <li> Of course, reduction operations (like norms) make use of
     * collective all-to-all MPI communications.
     * </ul>
     *
     * This vector can take two different states with respect to ghost
     * elements:
     * <ul>
     * <li> After creation and whenever zero_out_ghosts() is called (or
     * <code>operator= (0.)</code>), the vector does only allow writing into
     * ghost elements but not reading from ghost elements.
     * <li> After a call to update_ghost_values(), the vector does not allow
     * writing into ghost elements but only reading from them. This is to
     * avoid undesired ghost data artifacts when calling compress() after
     * modifying some vector entries. The current status of the ghost entries
     * (read mode or write mode) can be queried by the method
     * has_ghost_elements(), which returns <code>true</code> exactly when
     * ghost elements have been updated and <code>false</code> otherwise,
     * irrespective of the actual number of ghost entries in the vector layout
     * (for that information, use n_ghost_entries() instead).
     * </ul>
     *
     * This vector uses the facilities of the class dealii::Vector<Number> for
     * implementing the operations on the local range of the vector. In
     * particular, it also inherits thread parallelism that splits most
     * vector-vector operations into smaller chunks if the program uses
     * multiple threads. This may or may not be desired when working also with
     * MPI.
     *
     * <h4>Limitations regarding the vector size</h4>
     *
     * This vector class is based on two different number types for indexing.
     * The so-called global index type encodes the overall size of the vector.
     * Its type is types::global_dof_index. The largest possible value is
     * <code>2^32-1</code> or approximately 4 billion in case 64 bit integers
     * are disabled at configuration of deal.II (default case) or
     * <code>2^64-1</code> or approximately <code>10^19</code> if 64 bit
     * integers are enabled (see the glossary entry on
     * @ref GlobalDoFIndex
     * for further information).
     *
     * The second relevant index type is the local index used within one MPI
     * rank. As opposed to the global index, the implementation assumes 32-bit
     * unsigned integers unconditionally. In other words, to actually use a
     * vector with more than four billion entries, you need to use MPI with
     * more than one rank (which in general is a safe assumption since four
     * billion entries consume at least 16 GB of memory for floats or 32 GB of
     * memory for doubles) and enable 64-bit indices. If more than 4 billion
     * local elements are present, the implementation tries to detect that,
     * which triggers an exception and aborts the code. Note, however, that
     * the detection of overflow is tricky and the detection mechanism might
     * fail in some circumstances. Therefore, it is strongly recommended to
     * not rely on this class to automatically detect the unsupported case.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
     */
    using LinearAlgebra::distributed::Vector;

    /*@}*/
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
