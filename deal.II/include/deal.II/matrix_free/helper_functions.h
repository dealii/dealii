//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__matrix_free_helper_functions_h
#define __deal2__matrix_free_helper_functions_h


#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <boost/functional/hash.hpp>

DEAL_II_NAMESPACE_OPEN



namespace internal
{
namespace MatrixFreeFunctions
{
  // forward declaration of internal data structure
  template <typename Number> struct ConstraintValues;

                                /**
                                 * A struct that collects all information
                                 * related to parallelization with threads:
                                 * The work is subdivided into tasks that can
                                 * be done independently.
                                 */
  struct TaskInfo
  {
                                /**
                                 * Constructor.
                                 */
    TaskInfo ();

                                /**
                                 * Clears all the data fields and resets them
                                 * to zero.
                                 */
    void clear ();

                                /**
                                 * Returns the memory consumption of
                                 * the class.
                                 */
    std::size_t memory_consumption () const;

    unsigned int block_size;
    unsigned int n_blocks;
    unsigned int block_size_last;
    unsigned int position_short_block;
    bool use_multithreading;
    bool use_partition_partition;
    bool use_coloring_only;

    std::vector<unsigned int> partition_color_blocks_row_index;
    std::vector<unsigned int> partition_color_blocks_data;
    unsigned int evens;
    unsigned int odds;
    unsigned int n_blocked_workers;
    unsigned int n_workers;

    std::vector<unsigned int> partition_evens;
    std::vector<unsigned int> partition_odds;
    std::vector<unsigned int> partition_n_blocked_workers;
    std::vector<unsigned int> partition_n_workers;
  };



                                /**
                                 * A struct that collects all information
                                 * related to the size of the problem and MPI
                                 * parallelization.
                                 */
  struct SizeInfo
  {
                                /**
                                 * Constructor.
                                 */
    SizeInfo ();

                                /**
                                 * Clears all data fields and resets the sizes
                                 * to zero.
                                 */
    void clear();
    
                                /**
                                 * Prints minimum, average, and
                                 * maximal memory consumption over the
                                 * MPI processes.
                                 */
    template <typename STREAM>
    void print_memory_statistics (STREAM     &out,
                                  std::size_t data_length) const;

                                /**
                                 * Determines the position of cells
                                 * with ghosts for distributed-memory
                                 * calculations.
                                 */
    void make_layout (const unsigned int n_active_cells_in,
                      const unsigned int vectorization_length_in,
                      std::vector<unsigned int> &boundary_cells,
                      std::vector<unsigned int> &irregular_cells);

    unsigned int n_active_cells;
    unsigned int n_macro_cells;
    unsigned int boundary_cells_start;
    unsigned int boundary_cells_end;
    unsigned int vectorization_length;

                                /**
                                 * index sets to describe the layout of cells:
                                 * locally owned cells and locally active
                                 * cells
                                 */
    IndexSet locally_owned_cells;
    IndexSet ghost_cells;

                                /**
                                 * MPI communicator
                                 */
    MPI_Comm communicator;
    unsigned int my_pid;
    unsigned int n_procs;
  };

                                /**
                                 * Data type to identify cell type.
                                 */ 
  enum CellType {cartesian=0, affine=1, general=2, undefined=3};

  // ----------------- hash structure --------------------------------

                                /**
                                 * A class that is
                                 * used to quickly find out whether two
                                 * vectors of floating point numbers are the
                                 * same without going through all the
                                 * elements: store a hash value for each
                                 * vector. Generate the
                                 * hash value by a sum of all values
                                 * multiplied by random numbers (cast to
                                 * int). Of course, this is not a sure
                                 * criterion and one must manually check for
                                 * equality before this hash is telling
                                 * something useful. However, inequalities are
                                 * easily detected (unless roundoff spoils the
                                 * hash function)
                                 */
  struct HashValue
  {
                                // Constructor: sets the size of Number values
                                // with the typical magnitude that is to be
                                // expected.
    HashValue (const double element_size = 1.);
    
                                // get hash value for a vector of floating
                                // point numbers (which are assumed to be of
                                // order of magnitude one). Do this by first
                                // truncating everything that is smaller than
                                // the scaling (in order to eliminate noise
                                // from roundoff errors) and then calling the
                                // boost hash function
    unsigned int operator ()(const std::vector<double> &vec);

                                // get hash value for a tensor of rank
                                // two where the magnitude of the
                                // entries is given by the parameter
                                // weight
    template <int dim, typename number>
    unsigned int operator ()(const Tensor<2,dim,VectorizedArray<number> > 
                             &input,
                             const bool     is_diagonal);
 
    
    const double scaling;
  };

  // Note: Implementation in matrix_free.templates.h

} // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
