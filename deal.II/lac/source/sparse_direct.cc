//----------------------------  sparse_direct.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_direct.cc  ---------------------------


#include <lac/sparse_direct.h>
#include <base/memory_consumption.h>
#include <base/thread_management.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>

#include <iostream>
#include <list>
#include <typeinfo>

#include <sys/wait.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#ifndef DEAL_II_USE_DIRECT_ERRNO_H
#  include <errno.h>
#else
#  include </usr/include/errno.h>
#endif
#include <sys/errno.h>


// if we know that at least one of the HSL functions are there,
// include the respective include file. Otherwise save some CPU cycles
// in the compiler
#if defined(HAVE_HSL_MA27) || defined(HAVE_HSL_MA47)
#  include <hsl/hsl.h>
#endif

// if the HSL functions are not there, define them empty and throw an
// exception
#ifndef HAVE_HSL_MA27
namespace HSL
{
  namespace MA27
  {
    extern "C"
    void ma27ad_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  unsigned int       *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

    
    extern "C"
    void ma27bd_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };


    extern "C"
    void ma27cd_ (const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };


    extern "C" void ma27x1_ (unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
      

    extern "C" void ma27x2_ (unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
      
    
    extern "C" void ma27x3_ (const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
  }
};
#endif  // ifndef HAVE_HSL_MA27


#ifndef HAVE_HSL_MA47
namespace HSL
{
  namespace MA47
  {
    extern "C"
    void ma47id_ (double       *,
		  unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
    

    extern "C"
    void ma47ad_ (const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  unsigned int       *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  double             *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

      
    extern "C"
    void ma47bd_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  unsigned int       *,
		  double             *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

    
    extern "C"
    void ma47cd_ (const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  double             *,
		  unsigned int       *,
		  const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };    
  }
};
#endif   // ifndef HAVE_HSL_MA47



namespace CommunicationsLog
{
  enum Direction { put, get };
  
  struct Record 
  {
      Direction             direction;
      const std::type_info* type;
      unsigned int          count;
      unsigned int          scheduled_bytes;
      unsigned int          completed_bytes;
      std::string           description;
  };

  std::list<Record> communication_log;

  
  template <typename T>
  void record_communication (const Direction    direction,
                             const unsigned int count,
                             const unsigned int completed_bytes,
                             const std::string &descr)
  {
    Record record = {direction, &typeid(T), count,
                     sizeof(T)*count, completed_bytes, descr};
    communication_log.push_back (record);
  };


  void list_communication () 
  {
                                     // make sure only one thread is
                                     // writing out at a time
    static Threads::ThreadMutex write_lock;
    write_lock.acquire ();
    
    std::cerr << "++++++++++++++++++++++++++++++" << std::endl
              << "Communiction log history:" << std::endl;
    
    for (std::list<Record>::const_iterator i=communication_log.begin();
         i!=communication_log.end(); ++i)
      std::cerr << "++ "
                << (i->direction == put ? "put" : "get")
                << " "
                << i->count << " objects of type "
                << i->type->name()
                << ", " << i->completed_bytes
                << " of " << i->scheduled_bytes
                << " bytes completed, description="
                << i->description
                << std::endl;
    std::cerr << "++++++++++++++++++++++++++++++" << std::endl;

    write_lock.release ();
  };  
};



/**
 * Output an error message and terminate the program.
 */
void die (const std::string &text)
{
  std::cerr << "+++++ detached_ma27 driver: " << text
            << std::endl;
  CommunicationsLog::list_communication ();
  abort ();
};


/**
 * Output an error message and terminate the program. Write two error
 * codes.
 */
template <typename T1, typename T2>
void die (const std::string &text, const T1 t1, const T2 t2)
{
  std::cerr << "+++++ detached_ma27 driver: " << text
            << " code1=" << t1 << ", code2=" << t2
            << std::endl;
  CommunicationsLog::list_communication ();
  abort ();
};



/**
 * loop and check every once in a while whether the mother process is
 * still existing or has died without giving us due notice. if the
 * latter is the case, then also exit this process
 *
 * check by calling "kill(PID,0)", where PID is the pid of the
 * child. if the return value is non-null, then kill couldn't find
 * out about the process, so it is apparently gone
 */
extern "C"
void monitor_child_liveness (const pid_t child_pid) 
{
  while (true)
    {
      std::cerr << "+++ monitoring child " << child_pid << std::endl;
      
      int ret = kill (child_pid, 0);
      if (ret != 0)
        if ((ret == -1) && (errno == ESRCH))
          die ("Child process seems to have died!");
        else
          die ("Unspecified error while checking for other process!", ret, errno);

                                       // ok, master still running,
                                       // take a little rest and then
                                       // ask again
      sleep (10);
    };
};


Threads::ThreadManager child_liveness_watchers;



/* -------------------------- MA27 ---------------------------- */

Threads::ThreadMutex SparseDirectMA27::static_synchronisation_lock;


struct SparseDirectMA27::DetachedModeData 
{
                                     /**
                                      * Mutex to assure that only one
                                      * thread is currently talking
                                      * through the pipe.
                                      */
    Threads::ThreadMutex mutex;
    
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

    template <typename T>
    void put (const T *t, const size_t N, const char *debug_info) const
      {
        unsigned int count = 0;
        while (count < sizeof(T)*N)
          {
                                             // repeat writing until
                                             // syscall is not
                                             // interrupted
            int ret;
            do
              ret = write (server_client_pipe[1],
                           reinterpret_cast<const char *> (t),
                           sizeof(T) * N);
            while ((ret<0) && (errno==EINTR));
            if (ret < 0)
              die ("error on client side in 'put'", ret, errno);

            count += ret;
          };
        
        fflush (NULL);
        CommunicationsLog::
          record_communication<T> (CommunicationsLog::put, N, count, debug_info);
      };

    
    template <typename T>
    void get (T *t, const size_t N, const char *debug_info) const
      {
        unsigned int count = 0;
        while (count < sizeof(T)*N)
          {
            int ret;
            do
              ret = read (client_server_pipe[0],
                          reinterpret_cast<char *> (t) + count,
                          sizeof(T) * N - count);
            while ((ret<0) && (errno==EINTR));
            
            if (ret < 0)
              die ("error on client side in 'get'", ret, errno);

            count += ret;
          };
        
        CommunicationsLog::
          record_communication<T> (CommunicationsLog::get, N, count, debug_info);
      };    
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
		sparsity_pattern (0)
{};



SparseDirectMA27::~SparseDirectMA27() 
{
  if (detached_mode)
    if (detached_mode_data != 0)
      {

                                         // first close down client
        detached_mode_data->mutex.acquire ();
        write (detached_mode_data->server_client_pipe[1], "7", 1);
        detached_mode_data->mutex.release ();
        
                                         // then also delete data
        delete detached_mode_data;
        detached_mode_data = 0;
      };
};



void
SparseDirectMA27::set_detached_mode () 
{
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());
  detached_mode = true;
};



bool
SparseDirectMA27::detached_mode_set () const
{
  return detached_mode;
};



void
SparseDirectMA27::initialize (const SparsityPattern &sp)
{
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());


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
      pipe(detached_mode_data->server_client_pipe);
      pipe(detached_mode_data->client_server_pipe);      
                                       // fflush(NULL) is said to be a
                                       // good idea before fork()
      std::fflush(NULL);

                                       // now fork and create child
                                       // process
      detached_mode_data->child_pid = fork();
      if (detached_mode_data->child_pid == 0)
                                         // child process starts here
        {
                                           // copy read end of input
                                           // pipe to stdin, and
                                           // likewise with write end
                                           // of pipe to stdout
          dup2(detached_mode_data->server_client_pipe[0], 0);
          close(detached_mode_data->server_client_pipe[0]);

          dup2(detached_mode_data->client_server_pipe[1], 1);
          close(detached_mode_data->client_server_pipe[1]);

                                           // then dispose of this
                                           // copy of the program, and
                                           // run the detached solver
                                           // slave instead
          const char * const program_name = DEAL_II_PATH"/lib/bin/detached_ma27";
          const char * const child_argv[] = { program_name, 0 };
          execv(program_name, const_cast<char * const *>(child_argv));

                                           // usually execv does not
                                           // return. if it does, then an
                                           // error happened and we report it
                                           // herewith:
          AssertThrow (false,
                       ExcMessage ("execv returned, which it is not supposed to do!"));
          std::exit(1);
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
      const pid_t parent_pid = getpid();
      detached_mode_data->put (&parent_pid, 1, "parent_pid");

      Threads::spawn (child_liveness_watchers,
                      Threads::encapsulate(&monitor_child_liveness)
                      .collect_args (detached_mode_data->child_pid));
    };
  
  
				   // suppress error output if
				   // requested
  if (suppress_output)
    {
      const unsigned int LP = 0;
      call_ma27x3 (&LP);
    };
  
  sparsity_pattern = &sp;
  
  const unsigned int
    n_rows           = sparsity_pattern->n_rows();
  const unsigned int * const
    rowstart_indices = sparsity_pattern->get_rowstart_indices();
  const unsigned int * const
    col_nums         = sparsity_pattern->get_column_numbers();

				   // first count number of nonzero
				   // elements in the upper right
				   // part. the matrix is symmetric,
				   // so this suffices
  n_nonzero_elements = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
      if (row <= *col)
	++n_nonzero_elements;
  

				   // fill the row numbers and column
				   // numbers arrays from the sparsity
				   // pattern. note that we have
				   // Fortran convention, i.e. indices
				   // need to be 1-base, as opposed to
				   // C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if (row <= *col)
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  row_numbers[global_index] = row+1;
	  column_numbers[global_index] = *col+1;
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());
  
				   // initialize scratch arrays and
				   // variables
  LIW = static_cast<unsigned int>((2*n_nonzero_elements + 3*n_rows + 1) *
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
      LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);
    }
  while (true);

				   // if we were not allowed to
				   // allocate more memory, then throw
				   // an exception
  if (!call_succeeded)
    throw ExcMA27AFailed(IFLAG);

				   // catch returned values from the
				   // COMMON block. we need these
				   // values in order to set array
				   // sizes in the next function
  call_ma27x1 (&NRLNEC);
  call_ma27x2 (&NIRNEC);

				   // note that we have already been
				   // in this function
  initialize_called = true;
};



void
SparseDirectMA27::factorize (const SparseMatrix<double> &matrix)
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
      LIW = static_cast<unsigned int>(NIRNEC * LIW_factor_2);
      IW.resize (LIW);
    };
  
  const unsigned int n_rows = matrix.get_sparsity_pattern().n_rows();
  
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
          std::vector<unsigned int> tmp1, tmp2, tmp3;
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
	    
	    LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_2);
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
	    
	    LA  = static_cast<unsigned int>(LA * LA_increase_factor);
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
		throw ExcMA27BFailed(IFLAG);
	};
      continue;

      exit_loop:
      break;
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA27BFailed(IFLAG);

				   // note that we have been here
				   // already
  factorize_called = true;
};



void
SparseDirectMA27::solve (Vector<double> &rhs_and_solution) const
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());
  
  const unsigned int n_rows = rhs_and_solution.size();
  call_ma27cd (&n_rows, &A[0], &LA,
               &IW[0], &LIW, &MAXFRT,
               &rhs_and_solution(0), &IW1[0], &NSTEPS);
};



void
SparseDirectMA27::solve (const SparseMatrix<double> &matrix,
			 Vector<double>             &rhs_and_solution)
{
  initialize (matrix.get_sparsity_pattern());
  factorize (matrix);
  solve (rhs_and_solution);
};



unsigned int
SparseDirectMA27::memory_consumption () const
{
  return (sizeof(*this) +
	  MemoryConsumption::memory_consumption (row_numbers) +
	  MemoryConsumption::memory_consumption (column_numbers) +
	  MemoryConsumption::memory_consumption (A) +
	  MemoryConsumption::memory_consumption (IW) +
	  MemoryConsumption::memory_consumption (IKEEP) +
	  MemoryConsumption::memory_consumption (IW1));
};



Threads::ThreadMutex &
SparseDirectMA27::get_synchronisation_lock () const
{
  if (detached_mode)
    return non_static_synchronisation_lock;
  else
    return static_synchronisation_lock;    
};



void
SparseDirectMA27::fill_A (const SparseMatrix<double> &matrix)
{
  Assert (n_nonzero_elements <= A.size(), ExcInternalError());

  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();
  
  const unsigned int n_rows = sparsity_pattern.n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern.get_column_numbers();

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if (row <= *col)
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  A[global_index] = matrix(row,*col);
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());  
};



    
void SparseDirectMA27::call_ma27ad (const unsigned int *N,
                                    const unsigned int *NZ,
                                    const unsigned int *IRN,
                                    const unsigned int *ICN,
                                    unsigned int       *IW,
                                    const unsigned int *LIW,
                                    unsigned int       *IKEEP,
                                    unsigned int       *IW1,
                                    unsigned int       *NSTEPS,
                                    int                *IFLAG)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27ad_ (N, NZ, IRN, ICN, IW, LIW,
                        IKEEP, IW1, NSTEPS, IFLAG);
  else
    {
      detached_mode_data->mutex.acquire ();
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
      
      detached_mode_data->mutex.release ();
    };
};



void SparseDirectMA27::call_ma27bd (const unsigned int *N,
                                    const unsigned int *NZ,
                                    const unsigned int *IRN,
                                    const unsigned int *ICN,
                                    double             *A,
                                    const unsigned int *LA,
                                    unsigned int       *IW,
                                    const unsigned int *LIW,
                                    const unsigned int *IKEEP,
                                    const unsigned int *NSTEPS,
                                    unsigned int       *MAXFRT,
                                    unsigned int       *IW1,
                                    int                *IFLAG)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27bd_ (N, NZ, IRN, ICN, A, LA, IW, LIW,
                        IKEEP, NSTEPS, MAXFRT, IW1, IFLAG);
  else
    {
                                       // basically, everything is
                                       // already over the line,
                                       // except for A and LA
      detached_mode_data->mutex.acquire ();
      detached_mode_data->put ("2", 1, "ACTION 2");
      
      detached_mode_data->put (LA, 1, "LA");
      detached_mode_data->put (A,  *LA, "A");

                                       // next get back what we need
                                       // to know
      detached_mode_data->get (IFLAG, 1, "IFLAG");
      
      detached_mode_data->mutex.release ();
    };
};



void SparseDirectMA27::call_ma27cd (const unsigned int *N,
                                    const double       *A,
                                    const unsigned int *LA,
                                    const unsigned int *IW,
                                    const unsigned int *LIW,
                                    const unsigned int *MAXFRT,
                                    double             *RHS,
                                    const unsigned int *IW1,
                                    const unsigned int *NSTEPS) const
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
};



void SparseDirectMA27::call_ma27x1 (unsigned int *NRLNEC)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x1_ (NRLNEC);
  else
    {
      detached_mode_data->mutex.acquire ();
                                       // ma27x1 only reads data, so
                                       // don't send anything except
                                       // for the id
      detached_mode_data->put ("4", 1, "ACTION 4");
      detached_mode_data->get (NRLNEC, 1, "NRLNEC");
      detached_mode_data->mutex.release ();
    };
};



void SparseDirectMA27::call_ma27x2 (unsigned int *NIRNEC)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x2_ (NIRNEC);
  else
    {
      detached_mode_data->mutex.acquire ();
                                       // ma27x2 only reads data, so
                                       // don't send anything except
                                       // for the id
      detached_mode_data->put ("5", 1, "ACTION 5");
      detached_mode_data->get (NIRNEC, 1, "NIRNEC");
      detached_mode_data->mutex.release ();
    };
};



void SparseDirectMA27::call_ma27x3 (const unsigned int *LP)
{
  if (detached_mode_set() == false)
    HSL::MA27::ma27x3_ (LP);
  else
    {
      detached_mode_data->mutex.acquire ();
                                       // ma27x2 only reads data, so
                                       // don't send anything except
                                       // for the id
      detached_mode_data->put ("6", 1, "ACTION 6");
      detached_mode_data->put (LP, 1, "LP");
      detached_mode_data->mutex.release ();
    };
};

  



/* -------------------------- MA47 ---------------------------- */

Threads::ThreadMutex SparseDirectMA47::synchronisation_lock;


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
		matrix (0)
{};



void
SparseDirectMA47::initialize (const SparseMatrix<double> &m)
{
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());

                                   // some initialization stuff
  call_ma47id (CNTL, ICNTL);
  if (suppress_output)
    ICNTL[0] = 0;

                                   // then start with work
  matrix = &m;
  const SparsityPattern &sparsity_pattern = matrix->get_sparsity_pattern();
  
  const unsigned int
    n_rows           = sparsity_pattern.n_rows();
  const unsigned int * const
    rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int * const
    col_nums         = sparsity_pattern.get_column_numbers();

				   // first count number of nonzero
				   // elements in the upper right
				   // part. the matrix is symmetric,
				   // so this suffices
  n_nonzero_elements = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // skip zero elements, as
				       // required by the docs of MA47
      if ((row <= *col) && (m(row,*col) != 0))
	++n_nonzero_elements;
  

				   // fill the row numbers and column
				   // numbers arrays from the sparsity
				   // pattern. note that we have
				   // Fortran convention, i.e. indices
				   // need to be 1-base, as opposed to
				   // C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if ((row <= *col) && (m(row,*col) != 0))
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  row_numbers[global_index] = row+1;
	  column_numbers[global_index] = *col+1;
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());
  
				   // initialize scratch arrays and
				   // variables
  LIW = static_cast<unsigned int>((2*n_nonzero_elements + 5*n_rows + 4) *
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
      LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);      
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA47AFailed(INFO[0]);

				   // note that we have already been
				   // in this function
  initialize_called = true;
};



void
SparseDirectMA47::factorize (const SparseMatrix<double> &m)
{
  Assert (factorize_called == false, ExcCantFactorizeAgain());
  
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
      LIW = static_cast<unsigned int>(INFO[6] * LIW_factor_2);
      IW.resize (LIW);
    };

  const unsigned int n_rows = m.get_sparsity_pattern().n_rows();
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
	    
	    LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_2);
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
	    
	    LA  = static_cast<unsigned int>(LA * LA_increase_factor);
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
		throw ExcMA47BFailed(INFO[0]);
	};
      continue;

      exit_loop:
      break;
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA47BFailed(INFO[0]);

				   // note that we have been here
				   // already
  factorize_called = true;
};



void
SparseDirectMA47::solve (Vector<double> &rhs_and_solution)
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());
  
  const unsigned int n_rows = rhs_and_solution.size();
  call_ma47cd (&n_rows, &A[0], &LA,
               &IW[0], &LIW,
               &rhs_and_solution(0), &IW1[0], &ICNTL[0]);
};



void
SparseDirectMA47::solve (const SparseMatrix<double> &matrix,
			 Vector<double>             &rhs_and_solution)
{
  initialize (matrix);
  factorize (matrix);
  solve (rhs_and_solution);
};



unsigned int
SparseDirectMA47::memory_consumption () const
{
  return (sizeof(*this) +
	  MemoryConsumption::memory_consumption (row_numbers) +
	  MemoryConsumption::memory_consumption (column_numbers) +
	  MemoryConsumption::memory_consumption (A) +
	  MemoryConsumption::memory_consumption (IW) +
	  MemoryConsumption::memory_consumption (KEEP) +
	  MemoryConsumption::memory_consumption (IW1));
};



Threads::ThreadMutex &
SparseDirectMA47::get_synchronisation_lock () const
{
  return synchronisation_lock;
};



void
SparseDirectMA47::fill_A (const SparseMatrix<double> &matrix)
{
  Assert (n_nonzero_elements <= A.size(), ExcInternalError());

  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();
  
  const unsigned int n_rows = sparsity_pattern.n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern.get_column_numbers();

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if ((row <= *col) && (matrix(row,*col) != 0))
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  A[global_index] = matrix(row,*col);
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());  
};



void
SparseDirectMA47::call_ma47id (double       *CNTL,   // length 2
                               unsigned int *ICNTL)  // length 7
{
  HSL::MA47::ma47id_ (CNTL, ICNTL);
};



void
SparseDirectMA47::
call_ma47ad (const unsigned int *n_rows,             //scalar
             const unsigned int *n_nonzero_elements, //scalar
             unsigned int       *row_numbers,        //length n_nonzero
             unsigned int       *column_numbers,     //length n_nonzero
             unsigned int       *IW,                 //length LIW
             const unsigned int *LIW,                //scalar
             unsigned int       *KEEP,               //n_nonzero+5*n_rows+2
             const unsigned int *ICNTL,              //length 7
             int                *INFO)               //length 24
{
  double RINFO[4];  
  HSL::MA47::ma47ad_(n_rows, n_nonzero_elements,
                     row_numbers, column_numbers,
                     IW, LIW, KEEP,
                     ICNTL, &RINFO[0], INFO);
};



void
SparseDirectMA47::
call_ma47bd (const unsigned int *n_rows,             //scalar
             const unsigned int *n_nonzero_elements, //scalar
             const unsigned int *column_numbers,     //length n_nonzero
             double             *A,                  //length LA
             const unsigned int *LA,                 //scalar
             unsigned int       *IW,                 //length LIW
             const unsigned int *LIW,                //scalar
             const unsigned int *KEEP,               //n_nonzero+5*n_rows+2
             const double       *CNTL,               //length 2
             const unsigned int *ICNTL,              //length 7
             unsigned int       *IW1,                //2*n_rows+2
             int                *INFO)               //length 24
{
  double RINFO[4];  
  HSL::MA47::ma47bd_(n_rows, n_nonzero_elements, column_numbers,
                     A, LA,
                     IW, LIW, KEEP, CNTL, ICNTL,
                     IW1, &RINFO[0], INFO);
};



void
SparseDirectMA47::
call_ma47cd (const unsigned int *n_rows,           //scalar
             const double       *A,                //length LA
             const unsigned int *LA,               //scalar
             const unsigned int *IW,               //length LIW
             const unsigned int *LIW,              //scalar
             double             *rhs_and_solution, //length n_rows
             unsigned int       *IW1,              //length 2*n_rows+2
             const unsigned int *ICNTL)            //length 7
{
  std::vector<double> W(*n_rows);
  HSL::MA47::ma47cd_(n_rows, A, LA,
		     IW, LIW, &W[0],
		     rhs_and_solution, IW1, ICNTL);  
};
