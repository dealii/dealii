//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/lac/sparse_direct.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <cerrno>
#include <iostream>
#include <list>
#include <typeinfo>


DEAL_II_NAMESPACE_OPEN


// if we know that at least one of the HSL functions are there,
// include the respective include file. Otherwise save some CPU cycles
// in the compiler
#if defined(HAVE_HSL_MA27) || defined(HAVE_HSL_MA47)
#  include <hsl/hsl.h>
#endif

// include UMFPACK file.
#ifdef HAVE_LIBUMFPACK
extern "C" {
#  include <umfpack.h>
}
#endif

// if the HSL functions are not there, define them empty and throw an
// exception
#ifndef HAVE_HSL_MA27
namespace HSL
{
  namespace MA27
  {
    typedef types::global_dof_index size_type;

    extern "C"
    void ma27ad_ (const size_type *,
                  const size_type *,
                  const size_type *,
                  const size_type *,
                  size_type *,
                  const size_type *,
                  size_type *,
                  size_type *,
                  size_type *,
                  int *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C"
    void ma27bd_ (const size_type *,
                  const size_type *,
                  const size_type *,
                  const size_type *,
                  double *,
                  const size_type *,
                  size_type *,
                  const size_type *,
                  const size_type *,
                  const size_type *,
                  size_type *,
                  size_type *,
                  int *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C"
    void ma27cd_ (const size_type *,
                  const double *,
                  const size_type *,
                  const size_type *,
                  const size_type *,
                  double *,
                  const size_type *,
                  double *,
                  const size_type *,
                  const size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C" void ma27x1_ (size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C" void ma27x2_ (size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C" void ma27x3_ (const size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }
  }
}
#endif  // ifndef HAVE_HSL_MA27


#ifndef HAVE_HSL_MA47
namespace HSL
{
  namespace MA47
  {
    typedef types::global_dof_index size_type;

    extern "C"
    void ma47id_ (double *,
                  size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C"
    void ma47ad_ (const size_type *,
                  const size_type *,
                  size_type *,
                  size_type *,
                  size_type *,
                  const size_type *,
                  size_type *,
                  const size_type *,
                  double *,
                  int *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C"
    void ma47bd_ (const size_type *,
                  const size_type *,
                  const size_type *,
                  double *,
                  const size_type *,
                  size_type *,
                  const size_type *,
                  const size_type *,
                  const double *,
                  const size_type *,
                  size_type *,
                  double *,
                  int *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }


    extern "C"
    void ma47cd_ (const size_type *,
                  const double *,
                  const size_type *,
                  const size_type *,
                  const size_type *,
                  double *,
                  double *,
                  size_type *,
                  const size_type *)
    {
      AssertThrow (false,
                   ExcMessage("You can only use the HSL functions after putting "
                              "the respective files in the right place, "
                              "re-configuring the library and re-building it!"));
    }
  }
}
#endif   // ifndef HAVE_HSL_MA47




namespace
{
  /**
   * Output an error message and terminate the program.
   */
  void die (const std::string &text,
            const pid_t child)
  {
    std::cerr << "+++++ detached_ma27 driver(" << child << "): " << text
              << std::endl;
    std::abort ();
  }


  /**
   * Output an error message and terminate the program. Write two error
   * codes.
   */
  template <typename T1, typename T2>
  void die (const std::string &text,
            const T1 t1,
            const T2 t2,
            const pid_t child)
  {
    std::cerr << "+++++ detached_ma27 driver(" << child << "): " << text
              << " code1=" << t1 << ", code2=" << t2
              << std::endl;
    std::abort ();
  }
}



/* -------------------------- MA27 ---------------------------- */

Threads::Mutex SparseDirectMA27::static_synchronisation_lock;


struct SparseDirectMA27::DetachedModeData
{
  /**
   * Mutex to assure that only one
   * thread is currently talking
   * through the pipe.
   */
  Threads::Mutex mutex;

  /**
   * File handles for the pipe
   * between server (computing
   * process) and client (display
   * process).
   */
  int server_client_pipe[2];
  int client_server_pipe[2];

  /**
   * PID of the forked child
   * process.
   */
  pid_t child_pid;

  /**
   * Put a message from the server
   * to the client program. Obey
   * all the rules the operating
   * system sets, and create a log
   * entry for this communication
   */
  template <typename T>
  void put (const T *t,
            const std::size_t N,
            const char * /*debug_info*/) const
  {
    size_type count = 0;
    while (count < sizeof(T)*N)
      {
        // repeat writing until
        // syscall is not
        // interrupted
        int ret = -1;
#ifndef DEAL_II_MSVC
        do
          ret = write (server_client_pipe[1],
                       reinterpret_cast<const char *> (t) + count,
                       sizeof(T) * N - count);
        while ((ret<0) && (errno==EINTR));
#else
        Assert (false,
                ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif
        if (ret < 0)
          die ("error on client side in 'put'", ret, errno, child_pid);

        count += ret;
      };

    std::fflush (NULL);
  }


  /**
   * Get a message from the client
   * program. Obey all the rules
   * the operating system sets, and
   * create a log entry for this
   * communication
   */
  template <typename T>
  void get (T *t,
            const std::size_t N,
            const char * /*debug_info*/) const
  {
    size_type count = 0;
    while (count < sizeof(T)*N)
      {
        int ret = -1;
#ifndef DEAL_II_MSVC
        do
          ret = write (server_client_pipe[1],
                       reinterpret_cast<const char *> (t) + count,
                       sizeof(T) * N - count);
        while ((ret<0) && (errno==EINTR));
#else
        Assert (false,
                ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif

        if (ret < 0)
          die ("error on client side in 'get'", ret, errno, child_pid);

        count += ret;
      }
  }
};



SparseDirectMA27::SparseDirectMA27 (const double LIW_factor_1,
                                    const double LIW_factor_2,
                                    const double LA_factor,
                                    const double LIW_increase_factor_1,
                                    const double LIW_increase_factor_2,
                                    const double LA_increase_factor,
                                    const bool   suppress_output)
  :
  suppress_output (suppress_output),
  detached_mode (false),
  detached_mode_data (0),
  LIW_factor_1 (LIW_factor_1),
  LIW_factor_2 (LIW_factor_2),
  LA_factor (LA_factor),
  LIW_increase_factor_1 (LIW_increase_factor_1),
  LIW_increase_factor_2 (LIW_increase_factor_2),
  LA_increase_factor (LA_increase_factor),
  initialize_called (false),
  factorize_called (false),
  sparsity_pattern (0, typeid(*this).name())
{}



SparseDirectMA27::~SparseDirectMA27()
{
  if (detached_mode)
    if (detached_mode_data != 0)
      {
        // close down client
        Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
        // Assign the result of write
        // and reset the variable to
        // avoid compiler warnings
#ifndef DEAL_II_MSVC
//TODO:[WB] Shouldn't t be used to trace errors?
        ssize_t t = write (detached_mode_data->server_client_pipe[1], "7", 1);
        (void)t;
#else
        Assert (false,
                ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif
        // then also delete data
        delete detached_mode_data;
        detached_mode_data = 0;
      }
}



void
SparseDirectMA27::set_detached_mode ()
{
  Assert (initialize_called == false,
          ExcInitializeAlreadyCalled());
  detached_mode = true;
}



bool
SparseDirectMA27::detached_mode_set () const
{
  return detached_mode;
}



void
SparseDirectMA27::initialize (const SparsityPattern &sp)
{
  Assert (initialize_called == false,
          ExcInitializeAlreadyCalled());


  // first thing is: if detached mode
  // is requested, then we need to
  // spawn an instance of the
  // detached solver and open
  // communication channels with it
  if (detached_mode_set())
    {
      Assert (detached_mode_data == 0, ExcInternalError());
      detached_mode_data = new DetachedModeData();

      // create pipes to which we can
      // write and from which the
      // slave process will read its
      // stdin

      // Assign the return value to a
      // variable to avoid compiler
      // warnings
#ifndef DEAL_II_MSVC
//TODO:[WB] Use t to trace errors?
      int t = pipe(detached_mode_data->server_client_pipe);
      (void)t;
#else
      Assert (false,
              ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif
      // fflush(NULL) is said to be a
      // good idea before fork()
      std::fflush(NULL);

      // now fork and create child
      // process
#ifndef DEAL_II_MSVC
      // BG comment out until pipes are implemented in MSVC
      detached_mode_data->child_pid = fork();
#else
      Assert (false,
              ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif

      if (detached_mode_data->child_pid == 0)
        // child process starts here
        {
          // copy read end of input
          // pipe to stdin, and
          // likewise with write end
          // of pipe to stdout
#ifndef DEAL_II_MSVC
          dup2(detached_mode_data->server_client_pipe[0], 0);
          close(detached_mode_data->server_client_pipe[0]);

          dup2(detached_mode_data->client_server_pipe[1], 1);
          close(detached_mode_data->client_server_pipe[1]);

          // then dispose of this
          // copy of the program, and
          // run the detached solver
          // slave instead
          const char *const program_name = DEAL_II_PATH"/lib/bin/detached_ma27";
          const char *const child_argv[] = { program_name, NULL };
          execv(program_name, const_cast<char *const *>(child_argv));


          // usually execv does not
          // return. if it does, then an
          // error happened and we report it
          // herewith:
          AssertThrow (false,
                       ExcMessage ("execv returned, which it is not supposed to do!"));
          std::exit(1);

#else
          Assert (false,
                  ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif
        };
      // parent process continues
      // here.  first thing is to
      // send the process id of the
      // present process. this is
      // used to make sure that the
      // client can end itself when
      // it finds that the master
      // process was somehow
      // terminated without sending
      // him this information
#ifndef DEAL_II_MSVC
      const pid_t parent_pid = getpid();
      detached_mode_data->put (&parent_pid, 1, "parent_pid");
#else
      Assert (false,
              ExcMessage ("Detached mode isn't currently implemented on Windows"));
#endif
    };


  // suppress error output if
  // requested
  if (suppress_output)
    {
      const size_type LP = 0;
      call_ma27x3 (&LP);
    };

  sparsity_pattern = &sp;

  const size_type 
  n_rows           = sparsity_pattern->n_rows();

  // first count number of nonzero elements in the upper right part. the
  // matrix is symmetric, so this suffices
  n_nonzero_elements = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (SparsityPattern::iterator col = sparsity_pattern->begin(row);
         col < sparsity_pattern->end(row); ++col)
      if (row <= col->column())
        ++n_nonzero_elements;


  // fill the row numbers and column numbers arrays from the sparsity
  // pattern. note that we have Fortran convention, i.e. indices need to be
  // 1-base, as opposed to C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  size_type global_index = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (SparsityPattern::iterator col = sparsity_pattern->begin(row);
         col < sparsity_pattern->end(row); ++col)
      // note that the matrix must be
      // symmetric, so only treat the
      // upper right part
      if (row <= col->column())
        {
          Assert (global_index < n_nonzero_elements, ExcInternalError());

          row_numbers[global_index] = row+1;
          column_numbers[global_index] = col->column()+1;
          ++global_index;
        };
  Assert (global_index == n_nonzero_elements, ExcInternalError());

  // initialize scratch arrays and
  // variables
  LIW = static_cast<size_type>((2*n_nonzero_elements + 3*n_rows + 1) *
                                  LIW_factor_1);
  IW.resize    (detached_mode_set() ? 0 : LIW);
  IKEEP.resize (detached_mode_set() ? 0 : 3*n_rows);
  IW1.resize   (detached_mode_set() ? 0 : 2*n_rows);

  // no output please
  IFLAG = 0;

  // loop until memory requirements
  // are satisfied or we are not
  // allowed to allocate more memory
  // no more
  bool call_succeeded = true;
  do
    {
      call_ma27ad (&n_rows, &n_nonzero_elements,
                   &row_numbers[0], &column_numbers[0],
                   &IW[0], &LIW, &IKEEP[0],
                   &IW1[0], &NSTEPS, &IFLAG);
      call_succeeded = (IFLAG==0);

      // if enough memory or no
      // increase allowed: exit loop
      if (call_succeeded || (LIW_increase_factor_1 <= 1))
        break;

      // otherwise: increase LIW and retry
      LIW = static_cast<size_type>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);
    }
  while (true);

  // if we were not allowed to
  // allocate more memory, then throw
  // an exception
  AssertThrow (call_succeeded, ExcMA27AFailed(IFLAG));

  // catch returned values from the
  // COMMON block. we need these
  // values in order to set array
  // sizes in the next function
  call_ma27x1 (&NRLNEC);
  call_ma27x2 (&NIRNEC);

  // note that we have already been
  // in this function
  initialize_called = true;
}



template <typename number>
void
SparseDirectMA27::factorize (const SparseMatrix<number> &matrix)
{
  // if necessary, initialize process
  if (initialize_called == false)
    initialize (matrix.get_sparsity_pattern());

  // make sure the sparsity patterns
  // are the same
  Assert (sparsity_pattern == &matrix.get_sparsity_pattern(),
          ExcDifferentSparsityPatterns());


  // set LA and fill the A array of
  // values
  LA = std::max (static_cast<int>(NRLNEC * LA_factor),
                 static_cast<int>(n_nonzero_elements));
  A.resize (LA);
  fill_A (matrix);

  // if necessary extend IW
  if (LIW < NIRNEC * LIW_factor_2)
    {
      LIW = static_cast<size_type>(NIRNEC * LIW_factor_2);
      IW.resize (LIW);
    };

  const size_type n_rows = matrix.get_sparsity_pattern().n_rows();

  // loop until memory requirements
  // are satisfied or we are not
  // allowed to allocate more memory
  // no more
  bool call_succeeded = true;
  do
    {
      call_ma27bd (&n_rows, &n_nonzero_elements,
                   &row_numbers[0], &column_numbers[0],
                   &A[0], &LA,
                   &IW[0], &LIW, &IKEEP[0], &NSTEPS, &MAXFRT,
                   &IW1[0], &IFLAG);
      call_succeeded = (IFLAG==0);

      // if enough memory or no
      // increase allowed: exit
      // loop. delete data that is no
      // more used
      if (call_succeeded)
        {
          std::vector<size_type> tmp1, tmp2, tmp3;
          row_numbers.swap (tmp1);
          column_numbers.swap (tmp2);
          IKEEP.swap (tmp3);

          break;
        };


      // otherwise: increase LIW or
      // LA if that is allowed and
      // retry
      switch (IFLAG)
        {
        case -3:
        {
          if (LIW_increase_factor_2 <= 1)
            goto exit_loop;

          LIW = static_cast<size_type>(LIW * LIW_increase_factor_2);
          IW.resize (LIW);
          break;
        };

        case -4:
        {
          if (LA_increase_factor <= 1)
            goto exit_loop;
          // increase A. note that
          // since the function has
          // already part of the
          // array @p{A}, we have
          // to re-fill it with the
          // original values. minor
          // clue: since the old
          // entries are no more
          // needed, we can discard
          // them; we use this to
          // first release all
          // memory (through the
          // call to @p{swap} and
          // the subsequent call to
          // the destructor of the
          // @p{tmp} object) and
          // only then re-allocate
          // it. If we called
          // @p{resize} directly,
          // this would first
          // allocate more memory,
          // then copy the old
          // contents, and only
          // then release the old
          // memory, but keeping
          // both memory regions at
          // the same time could
          // sometimes be more than
          // we can do, leading to
          // an exception on the
          // allocation.
          std::cout << "<*>" << std::flush;

          LA  = static_cast<size_type>(LA * LA_increase_factor);
          if (true)
            {
              std::vector<double> tmp;
              A.swap (tmp);
            };

          A.resize (LA);
          fill_A (matrix);

          break;
        };

        // ups, other return
        // value, don't know
        // what to do here
        default:
          AssertThrow (false, ExcMA27BFailed(IFLAG));
        };
      continue;

exit_loop:
      break;
    }
  while (true);

  AssertThrow (call_succeeded, ExcMA27BFailed(IFLAG));

  // note that we have been here
  // already and release the sparsity
  // pattern object, since we won't
  // need it any more
  factorize_called = true;
  sparsity_pattern = 0;
}



template <>
void
SparseDirectMA27::solve (Vector<double> &rhs_and_solution) const
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());

  const size_type n_rows = rhs_and_solution.size();
  call_ma27cd (&n_rows, &A[0], &LA,
               &IW[0], &LIW, &MAXFRT,
               &rhs_and_solution(0), &IW1[0], &NSTEPS);
}



template <>
void
SparseDirectMA27::solve (Vector<float> &rhs_and_solution) const
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());

  // first have to convert data type to
  // doubles
  Vector<double> tmp (rhs_and_solution.size());
  tmp = rhs_and_solution;

  const size_type n_rows = rhs_and_solution.size();
  call_ma27cd (&n_rows, &A[0], &LA,
               &IW[0], &LIW, &MAXFRT,
               &tmp(0), &IW1[0], &NSTEPS);

  // then copy result back
  rhs_and_solution = tmp;
}



template <typename number>
void
SparseDirectMA27::solve (const SparseMatrix<number> &matrix,
                         Vector<double>             &rhs_and_solution)
{
  initialize (matrix.get_sparsity_pattern());
  factorize (matrix);
  solve (rhs_and_solution);
}



std::size_t
SparseDirectMA27::memory_consumption () const
{
  return (sizeof(*this) +
          MemoryConsumption::memory_consumption (row_numbers) +
          MemoryConsumption::memory_consumption (column_numbers) +
          MemoryConsumption::memory_consumption (A) +
          MemoryConsumption::memory_consumption (IW) +
          MemoryConsumption::memory_consumption (IKEEP) +
          MemoryConsumption::memory_consumption (IW1));
}



Threads::Mutex &
SparseDirectMA27::get_synchronisation_lock () const
{
  if (detached_mode)
    return non_static_synchronisation_lock;
  else
    return static_synchronisation_lock;
}



template <typename number>
void
SparseDirectMA27::fill_A (const SparseMatrix<number> &matrix)
{
  Assert (n_nonzero_elements <= A.size(), ExcInternalError());

  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();

  const size_type n_rows = sparsity_pattern.n_rows();

  size_type global_index = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (typename SparseMatrix<number>::const_iterator col=matrix.begin(row);
         col < matrix.end(row); ++col)
      // note that the matrix must be
      // symmetric, so only treat the
      // upper right part
      if (row <= col->column())
        {
          Assert (global_index < n_nonzero_elements, ExcInternalError());

          A[global_index] = col->value();
          ++global_index;

          // make sure that the symmetric
          // entry exists and has the same
          // value, unless this one is zero
          Assert ((col->value() == 0)
                  ||
                  (std::fabs(col->value() - matrix(col->column(),row))
                   <= 1e-15 * std::fabs (col->value())),
                  ExcMatrixNotSymmetric());
        }
      else
        // lower left part. just check
        // symmetry
        Assert ((col->value() == 0)
                ||
                (std::fabs(col->value() - matrix(col->column(),row))
                 <= 1e-15 * std::fabs (col->value())),
                ExcMatrixNotSymmetric());

  Assert (global_index == n_nonzero_elements, ExcInternalError());
}




void SparseDirectMA27::call_ma27ad (const size_type *N,
                                    const size_type *NZ,
                                    const size_type *IRN,
                                    const size_type *ICN,
                                    size_type       *IW,
                                    const size_type *LIW,
                                    size_type       *IKEEP,
                                    size_type       *IW1,
                                    size_type       *NSTEPS,
                                    int             *IFLAG)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27ad_ (N, NZ, IRN, ICN, IW, LIW,
                        IKEEP, IW1, NSTEPS, IFLAG);
  else
    {
      Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
      // first write the data we have
      // to push over, i.e. first
      // function index, then array
      // sizes, then arrays
      detached_mode_data->put ("1", 1, "ACTION 1");

      detached_mode_data->put (N,   1, "N");
      detached_mode_data->put (NZ,  1, "NZ");
      detached_mode_data->put (IRN, *NZ, "IRN");
      detached_mode_data->put (ICN, *NZ, "ICN");
      detached_mode_data->put (LIW, 1, "LIW");
      detached_mode_data->put (IFLAG, 1, "IFLAG");

      // all other fields are kept at
      // the client. array should not
      // be in used on this side
      Assert (this->IKEEP.size() == 0, ExcInternalError());
      Assert (this->IW1.size() == 0, ExcInternalError());

      // next get back what we need
      // to know
      detached_mode_data->get (IFLAG, 1, "IFLAG");
    };
}



void SparseDirectMA27::call_ma27bd (const size_type *N,
                                    const size_type *NZ,
                                    const size_type *IRN,
                                    const size_type *ICN,
                                    double          *A,
                                    const size_type *LA,
                                    size_type       *IW,
                                    const size_type *LIW,
                                    const size_type *IKEEP,
                                    const size_type *NSTEPS,
                                    size_type       *MAXFRT,
                                    size_type       *IW1,
                                    int             *IFLAG)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27bd_ (N, NZ, IRN, ICN, A, LA, IW, LIW,
                        IKEEP, NSTEPS, MAXFRT, IW1, IFLAG);
  else
    {
      // basically, everything is
      // already over the line,
      // except for A and LA
      Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
      detached_mode_data->put ("2", 1, "ACTION 2");

      detached_mode_data->put (LA, 1, "LA");
      detached_mode_data->put (A,  *LA, "A");

      // next get back what we need
      // to know
      detached_mode_data->get (IFLAG, 1, "IFLAG");
    };
}



void SparseDirectMA27::call_ma27cd (const size_type *N,
                                    const double    *A,
                                    const size_type *LA,
                                    const size_type *IW,
                                    const size_type *LIW,
                                    const size_type *MAXFRT,
                                    double          *RHS,
                                    const size_type *IW1,
                                    const size_type *NSTEPS) const
{
  if (detached_mode_set() == false)
    {
      std::vector<double> W(*MAXFRT);
      HSL::MA27::ma27cd_ (N, A, LA, IW, LIW, &W[0], MAXFRT, RHS, IW1, NSTEPS);
    }
  else
    {
      detached_mode_data->put ("3", 1, "ACTION 3");

      // we only have to push and get
      // the rhs vector
      detached_mode_data->put (RHS, *N, "RHS");
      detached_mode_data->get (RHS, *N, "RHS");
    };
}



void SparseDirectMA27::call_ma27x1 (size_type *NRLNEC)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x1_ (NRLNEC);
  else
    {
      Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
      // ma27x1 only reads data, so
      // don't send anything except
      // for the id
      detached_mode_data->put ("4", 1, "ACTION 4");
      detached_mode_data->get (NRLNEC, 1, "NRLNEC");
    };
}



void SparseDirectMA27::call_ma27x2 (size_type *NIRNEC)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x2_ (NIRNEC);
  else
    {
      Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
      // ma27x2 only reads data, so
      // don't send anything except
      // for the id
      detached_mode_data->put ("5", 1, "ACTION 5");
      detached_mode_data->get (NIRNEC, 1, "NIRNEC");
    };
}



void SparseDirectMA27::call_ma27x3 (const size_type *LP)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x3_ (LP);
  else
    {
      Threads::Mutex::ScopedLock lock (detached_mode_data->mutex);
      // ma27x2 only reads data, so
      // don't send anything except
      // for the id
      detached_mode_data->put ("6", 1, "ACTION 6");
      detached_mode_data->put (LP, 1, "LP");
    };
}





/* -------------------------- MA47 ---------------------------- */

Threads::Mutex SparseDirectMA47::synchronisation_lock;


SparseDirectMA47::SparseDirectMA47 (const double LIW_factor_1,
                                    const double LIW_factor_2,
                                    const double LA_factor,
                                    const double LIW_increase_factor_1,
                                    const double LIW_increase_factor_2,
                                    const double LA_increase_factor,
                                    const bool   suppress_output)
  :
  suppress_output (suppress_output),
  LIW_factor_1 (LIW_factor_1),
  LIW_factor_2 (LIW_factor_2),
  LA_factor (LA_factor),
  LIW_increase_factor_1 (LIW_increase_factor_1),
  LIW_increase_factor_2 (LIW_increase_factor_2),
  LA_increase_factor (LA_increase_factor),
  initialize_called (false),
  factorize_called (false),
  matrix (0, typeid(*this).name())
{}



void
SparseDirectMA47::initialize (const SparseMatrix<double> &m)
{
  Assert (initialize_called == false,
          ExcInitializeAlreadyCalled());

  // some initialization stuff
  call_ma47id (CNTL, ICNTL);
  if (suppress_output)
    ICNTL[0] = 0;

  // then start with work
  matrix = &m;
  const SparsityPattern &sparsity_pattern = matrix->get_sparsity_pattern();

  const size_type 
  n_rows           = sparsity_pattern.n_rows();

  // first count number of nonzero
  // elements in the upper right
  // part. the matrix is symmetric,
  // so this suffices
  n_nonzero_elements = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (SparseMatrix<double>::const_iterator col = m.begin(row);
         col < m.end(row); ++col)
      // skip zero elements, as required by the docs of MA47
      if (row <= col->column() && col->value() != 0)
        ++n_nonzero_elements;


  // fill the row numbers and column
  // numbers arrays from the sparsity
  // pattern. note that we have
  // Fortran convention, i.e. indices
  // need to be 1-base, as opposed to
  // C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  size_type global_index = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (SparseMatrix<double>::const_iterator col = m.begin(row);
         col < m.end(row); ++col)
      // note that the matrix must be
      // symmetric, so only treat the
      // upper right part
      if ((row <= col->column()) && (col->value() != 0))
        {
          Assert (global_index < n_nonzero_elements, ExcInternalError());

          row_numbers[global_index] = row+1;
          column_numbers[global_index] = col->column()+1;
          ++global_index;
        };
  Assert (global_index == n_nonzero_elements, ExcInternalError());

  // initialize scratch arrays and
  // variables
  LIW = static_cast<size_type>((2*n_nonzero_elements + 5*n_rows + 4) *
                                  LIW_factor_1);
  IW.resize (LIW);
  KEEP.resize (n_nonzero_elements + 5*n_rows + 2);

  // declare output info fields
  bool call_succeeded;
  do
    {
      call_ma47ad(&n_rows, &n_nonzero_elements,
                  &row_numbers[0], &column_numbers[0],
                  &IW[0], &LIW, &KEEP[0],
                  &ICNTL[0], &INFO[0]);
      call_succeeded = (INFO[0] == 0);

      // if enough memory or no
      // increase allowed: exit loop
      if (call_succeeded || (LIW_increase_factor_1 <= 1))
        break;

      // otherwise: increase LIW and retry
      LIW = static_cast<size_type>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);
    }
  while (true);

  AssertThrow (call_succeeded, ExcMA47AFailed(INFO[0]));

  // note that we have already been
  // in this function
  initialize_called = true;
}



void
SparseDirectMA47::factorize (const SparseMatrix<double> &m)
{
  Assert (factorize_called == false,
          ExcCantFactorizeAgain());

  // if necessary, initialize process
  if (initialize_called == false)
    initialize (m);

  // make sure the matrices
  // are the same
  Assert (matrix == &m, ExcDifferentMatrices());


  // set LA and fill the A array of
  // values
  LA = std::max (static_cast<int>(INFO[5] * LA_factor),
                 static_cast<int>(n_nonzero_elements));
  A.resize (LA);
  fill_A (m);

  // if necessary extend IW
  if (LIW < INFO[6] * LIW_factor_2)
    {
      LIW = static_cast<size_type>(INFO[6] * LIW_factor_2);
      IW.resize (LIW);
    };

  const size_type n_rows = m.get_sparsity_pattern().n_rows();
  IW1.resize (2*n_rows+2);

  // output info flags
  bool call_succeeded;
  do
    {
      call_ma47bd (&n_rows, &n_nonzero_elements, &column_numbers[0],
                   &A[0], &LA,
                   &IW[0], &LIW, &KEEP[0], &CNTL[0], &ICNTL[0],
                   &IW1[0], &INFO[0]);
      call_succeeded = (INFO[0] == 0);

      // if enough memory or no
      // increase allowed: exit loop
      if (call_succeeded)
        break;

      // otherwise: increase LIW or
      // LA if that is allowed and
      // retry
      switch (INFO[0])
        {
        case -3:
        {
          if (LIW_increase_factor_2 <= 1)
            goto exit_loop;

          LIW = static_cast<size_type>(LIW * LIW_increase_factor_2);
          IW.resize (LIW);
          break;
        };

        case -4:
        {
          if (LA_increase_factor <= 1)
            goto exit_loop;
          // increase A. note that
          // since the function has
          // already part of the
          // array @p{A}, we have
          // to re-fill it with the
          // original values. minor
          // clue: since the old
          // entries are no more
          // needed, we can discard
          // them; we use this to
          // first release all
          // memory (through the
          // call to @p{swap} and
          // the subsequent call to
          // the destructor of the
          // @p{tmp} object) and
          // only then re-allocate
          // it. If we called
          // @p{resize} directly,
          // this would first
          // allocate more memory,
          // then copy the old
          // contents, and only
          // then release the old
          // memory, but keeping
          // both memory regions at
          // the same time could
          // sometimes be more than
          // we can do, leading to
          // an exception on the
          // allocation.
          std::cout << "<*>" << std::flush;

          LA  = static_cast<size_type>(LA * LA_increase_factor);
          if (true)
            {
              std::vector<double> tmp;
              A.swap (tmp);
            };

          A.resize (LA);
          fill_A (m);

          break;
        };

        // ups, other return
        // value, don't know
        // what to do here
        default:
          AssertThrow (false, ExcMA47BFailed(INFO[0]));
        };
      continue;

exit_loop:
      break;
    }
  while (true);

  AssertThrow (call_succeeded, ExcMA47BFailed(INFO[0]));

  // note that we have been here
  // already
  factorize_called = true;
}



void
SparseDirectMA47::solve (Vector<double> &rhs_and_solution)
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());

  const size_type n_rows = rhs_and_solution.size();
  call_ma47cd (&n_rows, &A[0], &LA,
               &IW[0], &LIW,
               &rhs_and_solution(0), &IW1[0], &ICNTL[0]);
}



void
SparseDirectMA47::solve (const SparseMatrix<double> &matrix,
                         Vector<double>             &rhs_and_solution)
{
  initialize (matrix);
  factorize (matrix);
  solve (rhs_and_solution);
}



std::size_t
SparseDirectMA47::memory_consumption () const
{
  return (sizeof(*this) +
          MemoryConsumption::memory_consumption (row_numbers) +
          MemoryConsumption::memory_consumption (column_numbers) +
          MemoryConsumption::memory_consumption (A) +
          MemoryConsumption::memory_consumption (IW) +
          MemoryConsumption::memory_consumption (KEEP) +
          MemoryConsumption::memory_consumption (IW1));
}



Threads::Mutex &
SparseDirectMA47::get_synchronisation_lock () const
{
  return synchronisation_lock;
}



void
SparseDirectMA47::fill_A (const SparseMatrix<double> &matrix)
{
  Assert (n_nonzero_elements <= A.size(), ExcInternalError());

  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();

  const size_type n_rows = sparsity_pattern.n_rows();

  size_type global_index = 0;
  for (size_type row=0; row<n_rows; ++row)
    for (SparseMatrix<double>::const_iterator col=matrix.begin(row);
         col < matrix.end(row); ++col)
      // note that the matrix must be
      // symmetric, so only treat the
      // upper right part
      if ((row <= col->column()) && (col->value() != 0))
        {
          Assert (global_index < n_nonzero_elements, ExcInternalError());

          A[global_index] = col->value();
          ++global_index;

          // make sure that the symmetric
          // entry exists and has the same
          // value, unless this one is zero
          Assert ((col->value() == 0)
                  ||
                  (col->value() == matrix(col->column(),row)),
                  ExcMatrixNotSymmetric());
        }
      else
        // lower left part. just check
        // symmetry
        Assert ((col->value() == 0)
                ||
                (col->value() == matrix(col->column(),row)),
                ExcMatrixNotSymmetric());

  Assert (global_index == n_nonzero_elements, ExcInternalError());
}



void
SparseDirectMA47::call_ma47id (double    *CNTL,   // length 2
                               size_type *ICNTL)  // length 7
{
  HSL::MA47::ma47id_ (CNTL, ICNTL);
}



void
SparseDirectMA47::
call_ma47ad (const size_type    *n_rows,             //scalar
             const size_type    *n_nonzero_elements, //scalar
             size_type          *row_numbers,        //length n_nonzero
             size_type          *column_numbers,     //length n_nonzero
             size_type          *IW,                 //length LIW
             const size_type    *LIW,                //scalar
             size_type          *KEEP,               //n_nonzero+5*n_rows+2
             const unsigned int *ICNTL,              //length 7
             int                *INFO)               //length 24
{
  double RINFO[4];
  HSL::MA47::ma47ad_(n_rows, n_nonzero_elements,
                     row_numbers, column_numbers,
                     IW, LIW, KEEP,
                     ICNTL, &RINFO[0], INFO);
}



void
SparseDirectMA47::
call_ma47bd (const size_type    *n_rows,             //scalar
             const size_type    *n_nonzero_elements, //scalar
             const size_type    *column_numbers,     //length n_nonzero
             double             *A,                  //length LA
             const size_type    *LA,                 //scalar
             size_type          *IW,                 //length LIW
             const size_type    *LIW,                //scalar
             const size_type    *KEEP,               //n_nonzero+5*n_rows+2
             const double       *CNTL,               //length 2
             const unsigned int *ICNTL,              //length 7
             size_type          *IW1,                //2*n_rows+2
             int                *INFO)               //length 24
{
  double RINFO[4];
  HSL::MA47::ma47bd_(n_rows, n_nonzero_elements, column_numbers,
                     A, LA,
                     IW, LIW, KEEP, CNTL, ICNTL,
                     IW1, &RINFO[0], INFO);
}



void
SparseDirectMA47::
call_ma47cd (const size_type    *n_rows,           //scalar
             const double       *A,                //length LA
             const size_type    *LA,               //scalar
             const size_type    *IW,               //length LIW
             const size_type    *LIW,              //scalar
             double             *rhs_and_solution, //length n_rows
             size_type          *IW1,              //length 2*n_rows+2
             const unsigned int *ICNTL)            //length 7
{
  std::vector<double> W(*n_rows);
  HSL::MA47::ma47cd_(n_rows, A, LA,
                     IW, LIW, &W[0],
                     rhs_and_solution, IW1, ICNTL);
}



SparseDirectUMFPACK::~SparseDirectUMFPACK ()
{
  clear ();
}


void
SparseDirectUMFPACK::
initialize (const SparsityPattern &)
{}


#ifdef HAVE_LIBUMFPACK

SparseDirectUMFPACK::SparseDirectUMFPACK ()
  :
  symbolic_decomposition (0),
  numeric_decomposition (0),
  control (UMFPACK_CONTROL)
{
  umfpack_dl_defaults (&control[0]);
}



void
SparseDirectUMFPACK::clear ()
{
  // delete objects that haven't been deleted
  // yet
  if (symbolic_decomposition != 0)
    {
      umfpack_dl_free_symbolic (&symbolic_decomposition);
      symbolic_decomposition = 0;
    }

  if (numeric_decomposition != 0)
    {
      umfpack_dl_free_numeric (&numeric_decomposition);
      numeric_decomposition = 0;
    }

  {
    std::vector<long int> tmp;
    tmp.swap (Ap);
  }

  {
    std::vector<long int> tmp;
    tmp.swap (Ai);
  }

  {
    std::vector<double> tmp;
    tmp.swap (Ax);
  }

  umfpack_dl_defaults (&control[0]);
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const SparseMatrix<number> &matrix)
{
  // do the copying around of entries
  // so that the diagonal entry is in the
  // right place. note that this is easy to
  // detect: since all entries apart from the
  // diagonal entry are sorted, we know that
  // the diagonal entry is in the wrong place
  // if and only if its column index is
  // larger than the column index of the
  // second entry in a row
  //
  // ignore rows with only one or no entry
  for (size_type row=0; row<matrix.m(); ++row)
    {
      // we may have to move some elements
      // that are left of the diagonal but
      // presently after the diagonal entry
      // to the left, whereas the diagonal
      // entry has to move to the right. we
      // could first figure out where to
      // move everything to, but for
      // simplicity we just make a series
      // of swaps instead (this is kind of
      // a single run of bubble-sort, which
      // gives us the desired result since
      // the array is already "almost"
      // sorted)
      //
      // in the first loop, the condition
      // in the while-header also checks
      // that the row has at least two
      // entries and that the diagonal
      // entry is really in the wrong place
      long int cursor = Ap[row];
      while ((cursor < Ap[row+1]-1) &&
             (Ai[cursor] > Ai[cursor+1]))
        {
          std::swap (Ai[cursor], Ai[cursor+1]);
          std::swap (Ax[cursor], Ax[cursor+1]);
          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const SparseMatrixEZ<number> &matrix)
{
  //same thing for SparseMatrixEZ
  for (size_type row=0; row<matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      while ((cursor < Ap[row+1]-1) &&
             (Ai[cursor] > Ai[cursor+1]))
        {
          std::swap (Ai[cursor], Ai[cursor+1]);
          std::swap (Ax[cursor], Ax[cursor+1]);
          ++cursor;
        }
    }
}



template <typename number>
void
SparseDirectUMFPACK::
sort_arrays (const BlockSparseMatrix<number> &matrix)
{
  // the case for block matrices is a
  // bit more difficult, since all we
  // know is that *within each
  // block*, the diagonal of that
  // block may come first. however,
  // that means that there may be as
  // many entries per row in the
  // wrong place as there are block
  // columns. we can do the same
  // thing as above, but we have to
  // do it multiple times
  for (size_type row=0; row<matrix.m(); ++row)
    {
      long int cursor = Ap[row];
      for (size_type block=0; block<matrix.n_block_cols(); ++block)
        {

          // find the next
          // out-of-order element
          while ((cursor < Ap[row+1]-1) &&
                 (Ai[cursor] < Ai[cursor+1]))
            ++cursor;

          // if there is none, then
          // just go on
          if (cursor == Ap[row+1]-1)
            break;

          // otherwise swap this entry
          // with successive ones as
          // long as necessary
          long int element = cursor;
          while ((element < Ap[row+1]-1) &&
                 (Ai[element] > Ai[element+1]))
            {
              std::swap (Ai[element], Ai[element+1]);
              std::swap (Ax[element], Ax[element+1]);
              ++element;
            }
        }
    }
}



template <class Matrix>
void
SparseDirectUMFPACK::
factorize (const Matrix &matrix)
{
  Assert (matrix.m() == matrix.n(), ExcNotQuadratic())

  clear ();

  const size_type N = matrix.m();

  // copy over the data from the matrix to
  // the data structures UMFPACK wants. note
  // two things: first, UMFPACK wants
  // compressed column storage whereas we
  // always do compressed row storage; we
  // work around this by, rather than
  // shuffling things around, copy over the
  // data we have, but then call the
  // umfpack_dl_solve function with the
  // UMFPACK_At argument, meaning that we
  // want to solve for the transpose system
  //
  // second: the data we have in the sparse
  // matrices is "almost" right already;
  // UMFPACK wants the entries in each row
  // (i.e. really: column) to be sorted in
  // ascending order. we almost have that,
  // except that we usually store the
  // diagonal first in each row to allow for
  // some optimizations. thus, we have to
  // resort things a little bit, but only
  // within each row
  //
  // final note: if the matrix has entries in
  // the sparsity pattern that are actually
  // occupied by entries that have a zero
  // numerical value, then we keep them
  // anyway. people are supposed to provide
  // accurate sparsity patterns.
  Ap.resize (N+1);
  Ai.resize (matrix.n_nonzero_elements());
  Ax.resize (matrix.n_nonzero_elements());

  // first fill row lengths array
  Ap[0] = 0;
  for (size_type row=1; row<=N; ++row)
    Ap[row] = Ap[row-1] + matrix.get_row_length(row-1);
  Assert (static_cast<size_type>(Ap.back()) == Ai.size(),
          ExcInternalError());

  // then copy over matrix
  // elements. note that for sparse
  // matrices, iterators are sorted
  // so that they traverse each row
  // from start to end before moving
  // on to the next row. however,
  // this isn't true for block
  // matrices, so we have to do a bit
  // of book keeping
  {
    // have an array that for each
    // row points to the first entry
    // not yet written to
    std::vector<long int> row_pointers = Ap;

    // loop over the elements of the matrix row by row, as suggested
    // in the documentation of the sparse matrix iterator class
    for (size_type row = 0; row < matrix.m(); ++row)
      {
        for (typename Matrix::const_iterator p=matrix.begin(row);
            p!=matrix.end(row); ++p)
          {
            // write entry into the first
            // free one for this row
            Ai[row_pointers[row]] = p->column();
            Ax[row_pointers[row]] = p->value();

            // then move pointer ahead
            ++row_pointers[row];
          }
      }

    // at the end, we should have
    // written all rows completely
    for (size_type i=0; i<Ap.size()-1; ++i)
      Assert (row_pointers[i] == Ap[i+1], ExcInternalError());
  }

  // make sure that the elements in
  // each row are sorted. we have to
  // be more careful for block sparse
  // matrices, so ship this task out
  // to a different function
  sort_arrays (matrix);

  int status;
  status = umfpack_dl_symbolic (N, N,
                                &Ap[0], &Ai[0], &Ax[0],
                                &symbolic_decomposition,
                                &control[0], 0);
  AssertThrow (status == UMFPACK_OK,
               ExcUMFPACKError("umfpack_dl_symbolic", status));

  status = umfpack_dl_numeric (&Ap[0], &Ai[0], &Ax[0],
                               symbolic_decomposition,
                               &numeric_decomposition,
                               &control[0], 0);
  AssertThrow (status == UMFPACK_OK,
               ExcUMFPACKError("umfpack_dl_numeric", status));

  umfpack_dl_free_symbolic (&symbolic_decomposition) ;
}



void
SparseDirectUMFPACK::solve (Vector<double> &rhs_and_solution) const
{
  // make sure that some kind of factorize()
  // call has happened before
  Assert (Ap.size() != 0, ExcNotInitialized());
  Assert (Ai.size() != 0, ExcNotInitialized());
  Assert (Ai.size() == Ax.size(), ExcNotInitialized());

  Vector<double> rhs (rhs_and_solution.size());
  rhs = rhs_and_solution;

  // solve the system. note that since
  // UMFPACK wants compressed column storage
  // instead of the compressed row storage
  // format we use in deal.II's
  // SparsityPattern classes, we solve for
  // UMFPACK's A^T instead
  const int status
    = umfpack_dl_solve (UMFPACK_At,
                        &Ap[0], &Ai[0], &Ax[0],
                        rhs_and_solution.begin(), rhs.begin(),
                        numeric_decomposition,
                        &control[0], 0);
  AssertThrow (status == UMFPACK_OK, ExcUMFPACKError("umfpack_dl_solve", status));
}



template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix   &matrix,
                            Vector<double> &rhs_and_solution)
{
  factorize (matrix);
  solve (rhs_and_solution);
}


#else


SparseDirectUMFPACK::SparseDirectUMFPACK ()
  :
  symbolic_decomposition (0),
  numeric_decomposition (0),
  control (0)
{}


void
SparseDirectUMFPACK::clear ()
{}


template <class Matrix>
void SparseDirectUMFPACK::factorize (const Matrix &)
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but you called ./configure without the necessary --with-umfpack switch."));
}


void
SparseDirectUMFPACK::solve (Vector<double> &) const
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but you called ./configure without the necessary --with-umfpack switch."));
}


template <class Matrix>
void
SparseDirectUMFPACK::solve (const Matrix &,
                            Vector<double> &)
{
  AssertThrow(false, ExcMessage("To call this function you need UMFPACK, but you called ./configure without the necessary --with-umfpack switch."));
}


#endif


template <class Matrix>
void
SparseDirectUMFPACK::initialize (const Matrix        &M,
                                 const AdditionalData)
{
  this->factorize(M);
}


void
SparseDirectUMFPACK::vmult (
  Vector<double>       &dst,
  const Vector<double> &src) const
{
  dst = src;
  this->solve(dst);
}


void
SparseDirectUMFPACK::Tvmult (
  Vector<double> &,
  const Vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


void
SparseDirectUMFPACK::vmult_add (
  Vector<double> &,
  const Vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}


void
SparseDirectUMFPACK::Tvmult_add (
  Vector<double> &,
  const Vector<double> &) const
{
  Assert(false, ExcNotImplemented());
}

#ifdef DEAL_II_USE_MUMPS
SparseDirectMUMPS::SparseDirectMUMPS ()
  :
  initialize_called (false)
{}

SparseDirectMUMPS::~SparseDirectMUMPS ()
{
  id.job = -2;
  dmumps_c (&id);

  // Do some cleaning
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      delete[] a;
      delete[] irn;
      delete[] jcn;
    }
}

template <class Matrix>
void SparseDirectMUMPS::initialize_matrix (const Matrix &matrix)
{
  // Check we haven't been here before:
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());

  // Initialize MUMPS instance:
  id.job = -1;
  id.par =  1;
  id.sym =  0;

  // Use MPI_COMM_WORLD as communicator
  id.comm_fortran = -987654;
  dmumps_c (&id);

  // Hand over matrix and right-hand side
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      // Objects denoting a MUMPS data structure:
      //
      // Set number of unknowns
      n   = matrix.n ();

      // number of nonzero elements in matrix
      nz  = matrix.n_actually_nonzero_elements ();

      // representation of the matrix
      a   = new double[nz];

      // matrix indices pointing to the row and
      // column dimensions respectively of the
      // matrix representation above (a): ie. a[k]
      // is the matrix element (irn[k], jcn[k])
      irn = new int[nz];
      jcn = new int[nz];

      size_type index = 0;

      // loop over the elements of the matrix row by row, as suggested
      // in the documentation of the sparse matrix iterator class
      for (size_type row = 0; row < matrix.m(); ++row)
        {
          for (typename Matrix::const_iterator ptr = matrix.begin (row);
              ptr != matrix.end (row); ++ptr)
            if (std::abs (ptr->value ()) > 0.0)
              {
                a[index]   = ptr->value ();
                irn[index] = row + 1;
                jcn[index] = ptr->column () + 1;
                ++index;
              }
        }

      id.n   = n;
      id.nz  = nz;
      id.irn = irn;
      id.jcn = jcn;
      id.a   = a;
    }

  // No outputs
  id.icntl[0] = -1;
  id.icntl[1] = -1;
  id.icntl[2] = -1;
  id.icntl[3] =  0;

  // Exit by setting this flag:
  initialize_called = true;
}

template <class Matrix>
void SparseDirectMUMPS::initialize (const Matrix &matrix,
                                    const Vector<double>       &vector)
{
  // Hand over matrix and right-hand side
  initialize_matrix (matrix);

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      // Object denoting a MUMPS data structure
      rhs = new double[n];

      for (size_type i = 0; i < n; ++i)
        rhs[i] = vector (i);

      id.rhs = rhs;
    }
}

void SparseDirectMUMPS::copy_solution (Vector<double> &vector)
{
  // Copy solution into the given vector
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      for (size_type i=0; i<n; ++i)
        vector(i) = rhs[i];

      delete[] rhs;
    }
}

template <class Matrix>
void SparseDirectMUMPS::initialize (const Matrix &matrix)
{
  // Initialize MUMPS instance:
  initialize_matrix (matrix);
  // Start factorization
  id.job = 4;
  dmumps_c (&id);
}

void SparseDirectMUMPS::solve (Vector<double> &vector)
{
  // Check that the solver has been initialized
  // by the routine above:
  Assert (initialize_called == true, ExcNotInitialized());

  // and that the matrix has at least one
  // nonzero element:
  Assert (nz != 0, ExcNotInitialized());

  // Start solver
  id.job = 6;
  dmumps_c (&id);
  copy_solution (vector);
}

void SparseDirectMUMPS::vmult (Vector<double>       &dst,
                               const Vector<double> &src)
{
  // Check that the solver has been initialized
  // by the routine above:
  Assert (initialize_called == true, ExcNotInitialized());

  // and that the matrix has at least one
  // nonzero element:
  Assert (nz != 0, ExcNotInitialized());

  // Hand over right-hand side
  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    {
      // Object denoting a MUMPS data structure:
      rhs = new double[n];

      for (size_type i = 0; i < n; ++i)
        rhs[i] = src (i);

      id.rhs = rhs;
    }

  // Start solver
  id.job = 3;
  dmumps_c (&id);
  copy_solution (dst);
}

#endif // DEAL_II_USE_MUMPS

// explicit instantiations for SparseMatrixMA27
template
void SparseDirectMA27::factorize (const SparseMatrix<double> &matrix);

template
void SparseDirectMA27::factorize (const SparseMatrix<float> &matrix);

template
void SparseDirectMA27::solve (const SparseMatrix<double> &matrix,
                              Vector<double>             &rhs_and_solution);

template
void SparseDirectMA27::solve (const SparseMatrix<float> &matrix,
                              Vector<double>             &rhs_and_solution);


// explicit instantiations for SparseMatrixUMFPACK
#define InstantiateUMFPACK(MATRIX) \
  template    \
  void SparseDirectUMFPACK::factorize (const MATRIX &);    \
  template    \
  void SparseDirectUMFPACK::solve (const MATRIX   &,    \
                                   Vector<double> &);    \
  template    \
  void SparseDirectUMFPACK::initialize (const MATRIX &,    \
                                        const AdditionalData);

InstantiateUMFPACK(SparseMatrix<double>)
InstantiateUMFPACK(SparseMatrix<float>)
InstantiateUMFPACK(SparseMatrixEZ<double>)
InstantiateUMFPACK(SparseMatrixEZ<float>)
InstantiateUMFPACK(BlockSparseMatrix<double>)
InstantiateUMFPACK(BlockSparseMatrix<float>)

// explicit instantiations for SparseDirectMUMPS
#ifdef DEAL_II_USE_MUMPS
#define InstantiateMUMPS(MATRIX) \
  template \
  void SparseDirectMUMPS::initialize (const MATRIX &, const Vector<double> &);

InstantiateMUMPS(SparseMatrix<double>)
InstantiateMUMPS(SparseMatrix<float>)
// InstantiateMUMPS(SparseMatrixEZ<double>)
// InstantiateMUMPS(SparseMatrixEZ<float>)
InstantiateMUMPS(BlockSparseMatrix<double>)
InstantiateMUMPS(BlockSparseMatrix<float>)
#endif

DEAL_II_NAMESPACE_CLOSE
