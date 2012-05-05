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
  namespace internal
  {
    template <typename Number> struct ConstraintValues;
  }


                                // set minimum grain size for parallel
                                // computations
  namespace internal
  {
    const unsigned int minimum_parallel_grain_size = 500;
  }


                                    /*
                                    * Compressed data type to store a two
                                    * dimensional array. The data is stored in
                                    * a single standard vector. In a second
                                    * vector, the first element belonging to
                                    * each row is stored.
                                    */
  template<typename T>
  struct CompressedMatrix
  {
    AlignedVector<T> data;
    std::vector<unsigned int> row_index;
    T* operator[] (const unsigned int row) {
      return begin(row);
    };
    const T* operator[] (const unsigned int row) const {
      return begin(row);
    };
    const T* begin(const unsigned int row) const {
      AssertIndexRange (row, row_index.size()-1);
      return data.begin() + row_index[row];
    };
    const T* end(const unsigned int row) const {
      AssertIndexRange (row, row_index.size()-1);
      return data.begin() + row_index[row+1];
    };
    unsigned int row_length (const unsigned int row) const {
      AssertIndexRange (row, row_index.size()-1);
      return row_index[row+1] - row_index[row];
    };
    T* begin(const unsigned int row) {
      AssertIndexRange (row, row_index.size()-1);
      return data.begin() + row_index[row];
    };
    T* end(const unsigned int row) {
      AssertIndexRange (row, row_index.size()-1);
      return data.begin() + row_index[row+1];
    };
    void complete_last_row() {
      row_index.push_back (data.size());
    }
    void swap (CompressedMatrix<T> &other) {
      data.swap (other.data);
      row_index.swap (other.row_index);
    }
    void print (std::ostream &out) const
    {
      for (unsigned int row=0; row<row_index.size(); ++row)
        {
          for (const T* iterator=begin(row); iterator != end(row); ++iterator)
            out << *iterator << " ";
          out << std::endl;
        }
    };
    void clear()
    {
      data.clear();
      row_index.clear();
    }
    unsigned int memory_consumption() const
    {
      return MemoryConsumption::memory_consumption(data)+
        MemoryConsumption::memory_consumption(row_index);
    };
  };

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
    TaskInfo ()
    {
      clear();
    }

                                /**
                                 * Clears all the data fields and resets them
                                 * to zero.
                                 */
    void clear ()
    {
      block_size = 0;
      n_blocks = 0;
      block_size_last = 0;
      position_short_block = 0;
      use_multithreading = false;
      use_partition_partition = false;
      use_coloring_only = false;
      partition_color_blocks.clear();
      evens = 0;
      odds = 0;
      n_blocked_workers = 0;
      n_workers = 0;
      partition_evens.clear();
      partition_odds.clear();
      partition_n_blocked_workers.clear();
      partition_n_workers.clear();
    }

    std::size_t memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (partition_color_blocks) +
              MemoryConsumption::memory_consumption (partition_evens) +
              MemoryConsumption::memory_consumption (partition_odds) +
              MemoryConsumption::memory_consumption (partition_n_blocked_workers) +
              MemoryConsumption::memory_consumption (partition_n_workers));
    }

    unsigned int block_size;
    unsigned int n_blocks;
    unsigned int block_size_last;
    unsigned int position_short_block;
    bool use_multithreading;
    bool use_partition_partition;
    bool use_coloring_only;

    CompressedMatrix<unsigned int> partition_color_blocks;
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
    SizeInfo ()
    {
      clear();
    }

                                /**
                                 * Clears all data fields and resets the sizes
                                 * to zero.
                                 */
    void clear()
    {
      n_active_cells = 0;
      n_macro_cells  = 0;
      boundary_cells_start = 0;
      boundary_cells_end   = 0;
      n_vectors = 0;
      locally_owned_cells = IndexSet();
      ghost_cells = IndexSet();
      communicator = MPI_COMM_SELF;
      my_pid = 0;
      n_procs = 0;
    }

    template <typename STREAM>
    void print_mem (STREAM     &out,
                    std::size_t data_length) const
    {
      Utilities::MPI::MinMaxAvg memory_c;
      if (Utilities::System::job_supports_mpi() == true)
        {
          memory_c = Utilities::MPI::min_max_avg (1e-6*data_length,
                                                  communicator);
        }
      else
        {
          memory_c.sum = 1e-6*data_length;
          memory_c.min = memory_c.sum;
          memory_c.max = memory_c.sum;
          memory_c.avg = memory_c.sum;
          memory_c.min_index = 0;
          memory_c.max_index = 0;
        }
      if (n_procs < 2)
        out << memory_c.min;
      else
        out << memory_c.min << "/" << memory_c.avg << "/" << memory_c.max;
      out << " MB" << std::endl;
    }

    void make_layout (const unsigned int n_active_cells_in,
                      const unsigned int n_boundary_cells,
                      const unsigned int n_vectors_in,
                      std::vector<unsigned int> &irregular_cells)
    {
      n_vectors = n_vectors_in;
      n_active_cells = n_active_cells_in;

                                // check that number of boundary cells is
                                // divisible by n_vectors or that it contains
                                // all cells
      Assert (n_boundary_cells % n_vectors == 0 ||
              n_boundary_cells == n_active_cells, ExcInternalError());
      n_macro_cells = (n_active_cells+n_vectors-1)/n_vectors;
      irregular_cells.resize (n_macro_cells);
      if (n_macro_cells*n_vectors > n_active_cells)
        {
          irregular_cells[n_macro_cells-1] =
            n_vectors - (n_macro_cells*n_vectors - n_active_cells);
        }
      if (n_procs > 1)
        {
          const unsigned int n_macro_boundary_cells =
            (n_boundary_cells+n_vectors-1)/n_vectors;
          boundary_cells_start = (n_macro_cells-n_macro_boundary_cells)/2;
          boundary_cells_end   = boundary_cells_start + n_macro_boundary_cells;
        }
      else
        boundary_cells_start = boundary_cells_end = n_macro_cells;
    }

    unsigned int n_active_cells;
    unsigned int n_macro_cells;
    unsigned int boundary_cells_start;
    unsigned int boundary_cells_end;
    unsigned int n_vectors;

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



  namespace internal
  {
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
      HashValue (const double element_size = 1.)
        :
        scaling (element_size * std::numeric_limits<double>::epsilon() *
                 1024.)
      {};

                                // get hash value for a vector of floating
                                // point numbers (which are assumed to be of
                                // order of magnitude one). Do this by first
                                // truncating everything that is smaller than
                                // the scaling (in order to eliminate noise
                                // from roundoff errors) and then calling the
                                // boost hash function
      unsigned int operator ()(const std::vector<double> &vec)
      {
        std::vector<double> mod_vec(vec);
        for (unsigned int i=0; i<mod_vec.size(); ++i)
          mod_vec[i] -= fmod (mod_vec[i], scaling);
        return static_cast<unsigned int>(boost::hash_range (mod_vec.begin(), mod_vec.end()));
      };

                                // get hash value for a tensor of rank
                                // two where the magnitude of the
                                // entries is given by the parameter
                                // weight
      template <int dim, typename number>
      unsigned int operator ()(const Tensor<2,dim,VectorizedArray<number> > &input,
                               const bool     is_diagonal)
      {
        const unsigned int n_vectors = VectorizedArray<number>::n_array_elements;

        if (is_diagonal)
          {
            number mod_tensor [dim][n_vectors];
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<n_vectors; ++j)
                mod_tensor[i][j] = input[i][i][j] - fmod (input[i][i][j],
                                                          number(scaling));
            return static_cast<unsigned int>(boost::hash_range
                                             (&mod_tensor[0][0],
                                              &mod_tensor[0][0]+dim*n_vectors));
          }
        else
          {
            number mod_tensor [dim][dim][n_vectors];
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int d=0; d<dim; ++d)
                for (unsigned int j=0; j<n_vectors; ++j)
                  mod_tensor[i][d][j] = input[i][d][j] - fmod (input[i][d][j],
                                                               number(scaling));
            return static_cast<unsigned int>(boost::hash_range
                                             (&mod_tensor[0][0][0],
                                              &mod_tensor[0][0][0]+
                                              dim*dim*n_vectors));
          }
      };

      const double scaling;
    };

  } // end of namespace internal

} // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
