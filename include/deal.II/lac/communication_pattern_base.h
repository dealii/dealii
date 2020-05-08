// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_communication_pattern_base_h
#define dealii_communication_pattern_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class IndexSet;
#endif

namespace LinearAlgebra
{
  /**
   * CommunicationPattern is an abstract class that is used to define a
   * communication plan that can be called repeatedly to efficiently obtain
   * off-processor elements. The idea is to decouple the communication pattern
   * from the data that needs to be communicated. The goal is to reuse the same
   * communication pattern for different containers.
   * This is similar to the way SparseMatrix and SparsityPattern works.
   *
   * @author Bruno Turcksin, 2015.
   */
  class CommunicationPatternBase
  {
  public:
    /**
     * Destructor.
     */
    virtual ~CommunicationPatternBase() = default;

    /**
     * Reinitialize the communication pattern. The first argument
     * `vector_space_vector_index_set` is the index set associated to a
     * VectorSpaceVector object. The second argument
     * `read_write_vector_index_set` is the index set associated to a
     * ReadWriteVector object.
     */
    virtual void
    reinit(const IndexSet &vector_space_vector_index_set,
           const IndexSet &read_write_vector_index_set,
           const MPI_Comm &communicator) = 0;

    /**
     * Return a constant reference to the underlying MPI communicator.
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const = 0;
  };

} // end of namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
