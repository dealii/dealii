// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#ifndef dealii_vector_repartitioner_h
#define dealii_vector_repartitioner_h


#include <deal.II/base/config.h>

#include <deal.II/base/partitioner.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>


DEAL_II_NAMESPACE_OPEN

/**
 * Class for repartitioning data living on the same triangulation, which
 * has been partitioned differently.
 */
template <typename VectorType>
class VectorRepartitioner
{
public:
  /**
   * Transfer data.
   */
  void
  update_forwards(VectorType &dst, const VectorType &src) const;

  /**
   * Transfer data back.
   */
  void
  update_backwards(VectorType &dst, const VectorType &src) const;
};

/**
 * Class for repartitioning data living on the same triangulation, which
 * has been partitioned differently.
 *
 * Specialization for LinearAlgebra::distributed::Vector.
 */
template <typename Number>
class VectorRepartitioner<LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * Transfer data.
   */
  void
  update_forwards(LinearAlgebra::distributed::Vector<Number> &      dst,
                  const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Transfer data back.
   */
  void
  update_backwards(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Initialize transfer based on a destination and a source DoFHandler.
   */
  template <int dim, int spacedim>
  void
  reinit(const DoFHandler<dim, spacedim> &dof_handler_dst,
         const DoFHandler<dim, spacedim> &dof_handler_src);

private:
  /**
   * Partitioner needed by an intermediate vector, which is needed for
   * collecting all degrees of freedom of the destination cells.
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> extended_partitioner;

  /**
   * Indices for copying the data from/to the intermediate vector.
   */
  std::vector<unsigned int> indices;
};


DEAL_II_NAMESPACE_CLOSE

#endif
