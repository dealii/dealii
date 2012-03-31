//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011, 2012 by deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mpi_h
#define __deal2__mpi_h

#include <deal.II/base/config.h>
#include <vector>

#if defined(DEAL_II_COMPILER_SUPPORTS_MPI) || defined(DEAL_II_USE_PETSC)
#  include <mpi.h>
#else
                                 // without MPI, we would still like to use
                                 // some constructs with MPI data
                                 // types. Therefore, create some dummies
typedef int MPI_Comm;
const int MPI_COMM_SELF = 0;
typedef int MPI_Datatype;
typedef int MPI_Op;
namespace MPI
{
  static const unsigned int UNSIGNED = 0;
  static const unsigned int LONG_DOUBLE = 0;
  static const unsigned int LONG_DOUBLE_COMPLEX = 0;
  static const unsigned int MAX = 0;
  static const unsigned int MIN = 0;
  static const unsigned int SUM = 0;
}
#endif

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
                                   /**
                                    * A namespace for utility functions that
                                    * abstract certain operations using the
                                    * Message Passing Interface (MPI) or
                                    * provide fallback operations in
                                    * case deal.II is configured not to use
                                    * MPI at all.
                                    *
                                    * @ingroup utilities
                                    */
  namespace MPI
  {
                                     /**
                                      * Return the number of MPI processes
                                      * there exist in the given communicator
                                      * object. If this is a sequential job,
                                      * it returns 1.
                                      */
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator);

                                     /**
                                      * Return the number of the present MPI
                                      * process in the space of processes
                                      * described by the given
                                      * communicator. This will be a unique
                                      * value for each process between zero
                                      * and (less than) the number of all
                                      * processes (given by
                                      * get_n_mpi_processes()).
                                      */
    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator);

                                     /**
                                      * Consider an unstructured
                                      * communication pattern where
                                      * every process in an MPI
                                      * universe wants to send some
                                      * data to a subset of the other
                                      * processors. To do that, the
                                      * other processors need to know
                                      * who to expect messages
                                      * from. This function computes
                                      * this information.
                                      *
                                      * @param mpi_comm A communicator
                                      * that describes the processors
                                      * that are going to communicate
                                      * with each other.
                                      *
                                      * @param destinations The list
                                      * of processors the current
                                      * process wants to send
                                      * information to. This list need
                                      * not be sorted in any way. If
                                      * it contains duplicate entries
                                      * that means that multiple
                                      * messages are intended for a
                                      * given destination.
                                      *
                                      * @return A list of processors
                                      * that have indicated that they
                                      * want to send something to the
                                      * current processor. The
                                      * resulting list is not
                                      * sorted. It may contain
                                      * duplicate entries if
                                      * processors enter the same
                                      * destination more than once in
                                      * their destinations list.
                                      */
    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm & mpi_comm,
                                                  const std::vector<unsigned int> & destinations);

                                     /**
                                      * Given a communicator, generate a new
                                      * communicator that contains the same
                                      * set of processors but that has a
                                      * different, unique identifier.
                                      *
                                      * This functionality can be used to
                                      * ensure that different objects, such as
                                      * distributed matrices, each have unique
                                      * communicators over which they can
                                      * interact without interfering with each
                                      * other.
                                      *
                                      * When no longer needed, the
                                      * communicator created here needs to
                                      * be destroyed using
                                      * <code>MPI_Comm_free</code>.
                                      */
    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator);

                                     /**
                                      * Return the sum over all processors of the value @p t. This function
                                      * is collective over all processors given in the communicator. If
                                      * deal.II is not configured for use of MPI, this function simply
                                      * returns the value of @p t. This function corresponds to the
                                      * <code>MPI_Allreduce</code> function, i.e. all processors receive
                                      * the result of this operation.
                                      *
                                      * @note This function is only implemented for certain template
                                      * arguments <code>T</code>, namely <code>float, double, int,
                                      * unsigned int</code>.
                                      */
    template <typename T>
    T sum (const T &t,
           const MPI_Comm &mpi_communicator);

                                     /**
                                      * Like the previous function,
                                      * but take the sums over the
                                      * elements of an array
                                      * of length N. In other words,
                                      * the i-th element of the
                                      * results array is the sum over
                                      * the i-th entries of the input
                                      * arrays from each processor.
                                      */
    template <typename T, unsigned int N>
    inline
    void sum (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&sums)[N]);

                                     /**
                                      * Like the previous function,
                                      * but take the sums over the
                                      * elements of a std::vector. In other words,
                                      * the i-th element of the
                                      * results array is the sum over
                                      * the i-th entries of the input
                                      * arrays from each processor.
                                      */
    template <typename T>
    inline
    void sum (const std::vector<T> &values,
              const MPI_Comm &mpi_communicator,
              std::vector<T> &sums);

                                     /**
                                      * Return the maximum over all processors of the value @p t. This function
                                      * is collective over all processors given in the communicator. If
                                      * deal.II is not configured for use of MPI, this function simply
                                      * returns the value of @p t. This function corresponds to the
                                      * <code>MPI_Allreduce</code> function, i.e. all processors receive
                                      * the result of this operation.
                                      *
                                      * @note This function is only implemented for certain template
                                      * arguments <code>T</code>, namely <code>float, double, int,
                                      * unsigned int</code>.
                                      */
    template <typename T>
    T max (const T &t,
           const MPI_Comm &mpi_communicator);

                                     /**
                                      * Like the previous function,
                                      * but take the maxima over the
                                      * elements of an array
                                      * of length N. In other words,
                                      * the i-th element of the
                                      * results array is the maximum of
                                      * the i-th entries of the input
                                      * arrays from each processor.
                                      */
    template <typename T, unsigned int N>
    inline
    void max (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&maxima)[N]);

                                     /**
                                      * Like the previous function,
                                      * but take the maximum over the
                                      * elements of a std::vector. In other words,
                                      * the i-th element of the
                                      * results array is the maximum over
                                      * the i-th entries of the input
                                      * arrays from each processor.
                                      */
    template <typename T>
    inline
    void max (const std::vector<T> &values,
              const MPI_Comm &mpi_communicator,
              std::vector<T> &maxima);

                                     /**
                                      * Data structure to store the result of
                                      * min_max_avg().
                                      */
    struct MinMaxAvg
    {
        double sum;
        double min;
        double max;
        unsigned int min_index;
        unsigned int max_index;
        double avg;
    };

                                     /**
                                      * Returns sum, average, minimum,
                                      * maximum, processor id of minimum and
                                      * maximum as a collective operation of
                                      * on the given MPI communicator @param
                                      * mpi_communicator . Each processor's
                                      * value is given in @param my_value and
                                      * the result will be returned in @p
                                      * result . The result is available on all
                                      * machines.
                                      */
    MinMaxAvg
    min_max_avg (const double my_value,
                 const MPI_Comm &mpi_communicator);



                                     /**
                                      * A class that is used to initialize the
                                      * MPI system at the beginning of a
                                      * program and to shut it down again at
                                      * the end.
                                      *
                                      * If a program uses MPI one would
                                      * typically just create an object of
                                      * this type at the beginning of
                                      * <code>main()</code>. The constructor
                                      * of this class then runs
                                      * <code>MPI_Init()</code> with the given
                                      * arguments. At the end of the program,
                                      * the compiler will invoke the
                                      * destructor of this object which in
                                      * turns calls <code>MPI_Finalize</code>
                                      * to shut down the MPI system.
                                      *
                                      * This class is used in step-32, for example.
                                      */
    class MPI_InitFinalize
    {
      public:
                                         /**
                                          * Constructor. Takes the arguments
                                          * from the command line (in case of
                                          * MPI, the number of processes is
                                          * specified there), and sets up a
                                          * respective communicator by calling
                                          * <tt>MPI_Init()</tt>. This
                                          * constructor can only be called once
                                          * in a program, since MPI cannot be
                                          * initialized twice.
                                          */
        MPI_InitFinalize (int    &argc,
                          char** &argv);

                                         /**
                                          * Destructor. Calls
                                          * <tt>MPI_Finalize()</tt> in
                                          * case this class owns the MPI
                                          * process.
                                          */
        ~MPI_InitFinalize();

      private:
                                         /**
                                          * This flag tells the class
                                          * whether it owns the MPI
                                          * process (i.e., it has been
                                          * constructed using the
                                          * argc/argv input, or it has
                                          * been copied). In the former
                                          * case, the command
                                          * <tt>MPI_Finalize()</tt> will
                                          * be called at destruction.
                                          */
        const bool owns_mpi;
    };

    namespace internal
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                                       /**
                                        * Return the corresponding MPI data
                                        * type id for the argument given.
                                        */
      inline MPI_Datatype mpi_type_id (const int *)
      {
        return MPI_INT;
      }


      inline MPI_Datatype mpi_type_id (const long int *)
      {
        return MPI_LONG;
      }


      inline MPI_Datatype mpi_type_id (const unsigned int *)
      {
        return MPI_UNSIGNED;
      }


      inline MPI_Datatype mpi_type_id (const unsigned long int *)
      {
        return MPI_UNSIGNED_LONG;
      }


      inline MPI_Datatype mpi_type_id (const float *)
      {
        return MPI_FLOAT;
      }


      inline MPI_Datatype mpi_type_id (const double *)
      {
        return MPI_DOUBLE;
      }


      inline MPI_Datatype mpi_type_id (const long double *)
      {
        return MPI_LONG_DOUBLE;
      }
#endif
    }


    template <typename T>
    inline
    T sum (const T &t,
           const MPI_Comm &mpi_communicator)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      T sum;
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&t)),
                     &sum, 1, internal::mpi_type_id(&t), MPI_SUM,
                     mpi_communicator);
      return sum;
#else
      (void)mpi_communicator;
      return t;
#endif
    }


    template <typename T, unsigned int N>
    inline
    void sum (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&sums)[N])
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&values[0])),
                     &sums[0], N, internal::mpi_type_id(values), MPI_SUM,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      for (unsigned int i=0; i<N; ++i)
        sums[i] = values[i];
#endif
    }


    template <typename T>
    inline
    void sum (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &sums)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      sums.resize (values.size());
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&values[0])),
                     &sums[0], values.size(), internal::mpi_type_id((T*)0), MPI_SUM,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      sums = values;
#endif
    }


    template <typename T>
    inline
    T max (const T &t,
           const MPI_Comm &mpi_communicator)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      T sum;
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&t)),
                     &sum, 1, internal::mpi_type_id(&t), MPI_MAX,
                     mpi_communicator);
      return sum;
#else
      (void)mpi_communicator;
      return t;
#endif
    }


    template <typename T, unsigned int N>
    inline
    void max (const T (&values)[N],
              const MPI_Comm &mpi_communicator,
              T (&maxima)[N])
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&values[0])),
                     &maxima[0], N, internal::mpi_type_id(values), MPI_MAX,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      for (unsigned int i=0; i<N; ++i)
        maxima[i] = values[i];
#endif
    }


    template <typename T>
    inline
    void max (const std::vector<T> &values,
              const MPI_Comm       &mpi_communicator,
              std::vector<T>       &maxima)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      maxima.resize (values.size());
      MPI_Allreduce (const_cast<void*>(static_cast<const void*>(&values[0])),
                     &maxima[0], values.size(), internal::mpi_type_id((T*)0), MPI_MAX,
                     mpi_communicator);
#else
      (void)mpi_communicator;
      maxima = values;
#endif
    }
  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
