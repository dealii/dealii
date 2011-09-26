//---------------------------------------------------------------------------
//      $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cerrno>
#include <cmath>
#include <unistd.h>
#include <sys/types.h>
#include <sstream>
#include <iostream>

#ifdef HAVE_STD_NUMERIC_LIMITS
#  include <limits>
#else
#  include <limits.h>
#endif

#ifdef DEAL_II_USE_TRILINOS
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#    include <Epetra_MpiComm.h>
#    include <deal.II/lac/vector_memory.h>
#    include <deal.II/lac/trilinos_vector.h>
#    include <deal.II/lac/trilinos_block_vector.h>
#  endif
#  include "Teuchos_RCP.hpp"
#  include "Epetra_SerialComm.h"
#endif



DEAL_II_NAMESPACE_OPEN


namespace Utilities
{


  DeclException2 (ExcInvalidNumber2StringConversersion,
		  unsigned int, unsigned int,
		  << "When trying to convert " << arg1
		  << " to a string with " << arg2 << " digits");
  DeclException1 (ExcInvalidNumber,
		  unsigned int,
		  << "Invalid number " << arg1);
  DeclException1 (ExcCantConvertString,
		  std::string,
		  << "Can't convert the string " << arg1
                  << " to the desired type");


  std::string
  int_to_string (const unsigned int i,
		 const unsigned int digits)
  {
				     // if second argument is invalid, then do
				     // not pad the resulting string at all
    if (digits == numbers::invalid_unsigned_int)
      return int_to_string (i, needed_digits(i));


    AssertThrow ( ! ((digits==1 && i>=10)   ||
		     (digits==2 && i>=100)  ||
		     (digits==3 && i>=1000) ||
		     (digits==4 && i>=10000)||
		     (digits==5 && i>=100000)||
		     (i>=1000000)),
		  ExcInvalidNumber2StringConversersion(i, digits));

    std::string s;
    switch (digits)
      {
	case 6:
	      s += '0' + i/100000;
	case 5:
	      s += '0' + (i%100000)/10000;
	case 4:
	      s += '0' + (i%10000)/1000;
	case 3:
	      s += '0' + (i%1000)/100;
	case 2:
	      s += '0' + (i%100)/10;
	case 1:
	      s += '0' + i%10;
	      break;
	default:
	      s += "invalid digits information";
      };
    return s;
  }



  unsigned int
  needed_digits (const unsigned int max_number)
  {
    if (max_number < 10)
      return 1;
    if (max_number < 100)
      return 2;
    if (max_number < 1000)
      return 3;
    if (max_number < 10000)
      return 4;
    if (max_number < 100000)
      return 5;
    if (max_number < 1000000)
      return 6;
    AssertThrow (false, ExcInvalidNumber(max_number));
    return 0;
  }



  int
  string_to_int (const std::string &s)
  {
    std::istringstream ss(s);

#ifdef HAVE_STD_NUMERIC_LIMITS
    static const int max_int = std::numeric_limits<int>::max();
#else
    static const int max_int = INT_MAX;
#endif

    int i = max_int;
    ss >> i;
                                     // check for errors
    AssertThrow (i != max_int, ExcCantConvertString (s));

//TODO: The test for errors above doesn't work, as can easily be
//verified. furthermore, it doesn't catch cases like when calling
//string_to_int("1.23.4") since it just reads in however much it can, without
//realizing that there is more

    return i;
  }



  std::vector<int>
  string_to_int (const std::vector<std::string> &s)
  {
    std::vector<int> tmp (s.size());
    for (unsigned int i=0; i<s.size(); ++i)
      tmp[i] = string_to_int (s[i]);
    return tmp;
  }



  double
  string_to_double (const std::string &s)
  {
    std::istringstream ss(s);

#ifdef HAVE_STD_NUMERIC_LIMITS
    static const double max_double = std::numeric_limits<double>::max();
#else
    static const double max_double = DBL_MAX;
#endif

    double i = max_double;
    ss >> i;

                                     // check for errors
    AssertThrow (i != max_double, ExcCantConvertString (s));

//TODO: The test for errors above doesn't work, as can easily be
//verified. furthermore, it doesn't catch cases like when calling
//string_to_int("1.23.4") since it just reads in however much it can, without
//realizing that there is more

    return i;
  }



  std::vector<double>
  string_to_double (const std::vector<std::string> &s)
  {
    std::vector<double> tmp (s.size());
    for (unsigned int i=0; i<s.size(); ++i)
      tmp[i] = string_to_double (s[i]);
    return tmp;
  }



  std::vector<std::string>
  split_string_list (const std::string &s,
                     const char         delimiter)
  {
    std::string tmp = s;
    std::vector<std::string> split_list;
    split_list.reserve (std::count (tmp.begin(), tmp.end(), delimiter)+1);

				     // split the input list
    while (tmp.length() != 0)
      {
        std::string name;
	name = tmp;

	if (name.find(delimiter) != std::string::npos)
	  {
	    name.erase (name.find(delimiter), std::string::npos);
	    tmp.erase (0, tmp.find(delimiter)+1);
	  }
	else
	  tmp = "";

	while ((name.length() != 0) &&
	       (name[0] == ' '))
	  name.erase (0,1);

	while (name[name.length()-1] == ' ')
	  name.erase (name.length()-1, 1);

	split_list.push_back (name);
      }

    return split_list;
  }



  std::vector<std::string>
  break_text_into_lines (const std::string &original_text,
                         const unsigned int width,
                         const char delimiter)
  {
    std::string              text = original_text;
    std::vector<std::string> lines;

                                     // remove trailing spaces
    while ((text.length() != 0) && (text[text.length()-1] == delimiter))
      text.erase(text.length()-1,1);

                                     // then split the text into lines
    while (text.length() != 0)
      {
                                         // in each iteration, first remove
                                         // leading spaces
        while ((text.length() != 0) && (text[0] == delimiter))
          text.erase(0, 1);

                                         // if we can fit everything into one
                                         // line, then do so. otherwise, we have
                                         // to keep breaking
        if (text.length() < width)
          {
            lines.push_back (text);
            text = "";
          }
        else
          {
                                             // starting at position width, find the
                                             // location of the previous space, so
                                             // that we can break around there
            int location = std::min<int>(width,text.length()-1);
            for (; location>0; --location)
              if (text[location] == delimiter)
                break;

                                             // if there are no spaces, then try if
                                             // there are spaces coming up
            if (location == 0)
              for (location = std::min<int>(width,text.length()-1);
                   location<static_cast<int>(text.length());
                   ++location)
                if (text[location] == delimiter)
                  break;

                                             // now take the text up to the found
                                             // location and put it into a single
                                             // line, and remove it from 'text'
            lines.push_back (std::string (text, 0, location));
            text.erase (0, location);
          }
      }

    return lines;
  }



  bool
  match_at_string_start (const std::string &name,
			 const std::string &pattern)
  {
    if (pattern.size() > name.size())
      return false;

    for (unsigned int i=0; i<pattern.size(); ++i)
      if (pattern[i] != name[i])
	return false;

    return true;
  }



  std::pair<int, unsigned int>
  get_integer_at_position (const std::string &name,
			   const unsigned int position)
  {
    Assert (position < name.size(), ExcInternalError());

    const std::string test_string (name.begin()+position,
				   name.end());

    std::istringstream str(test_string);

    int i;
    if (str >> i)
      {
					 // compute the number of
					 // digits of i. assuming it
					 // is less than 8 is likely
					 // ok
	if (i<10)
	  return std::make_pair (i, 1U);
	else if (i<100)
	  return std::make_pair (i, 2U);
	else if (i<1000)
	  return std::make_pair (i, 3U);
	else if (i<10000)
	  return std::make_pair (i, 4U);
	else if (i<100000)
	  return std::make_pair (i, 5U);
	else if (i<1000000)
	  return std::make_pair (i, 6U);
	else if (i<10000000)
	  return std::make_pair (i, 7U);
	else
	  {
	    Assert (false, ExcNotImplemented());
	    return std::make_pair (-1, numbers::invalid_unsigned_int);
	  }
      }
    else
      return std::make_pair (-1, numbers::invalid_unsigned_int);
  }



  double
  generate_normal_random_number (const double a,
				 const double sigma)
  {
				     // if no noise: return now
    if (sigma == 0)
      return a;

#ifdef HAVE_RAND_R
    static unsigned int seed = 0xabcd1234;
    const double y = 1.0*rand_r(&seed)/RAND_MAX;
#else
    const double y = 1.0*rand()/RAND_MAX;
#endif

				     // find x such that y=erf(x). do so
				     // using a Newton method to find
				     // the zero of F(x)=erf(x)-y. start
				     // at x=0
    double x = 0;
    unsigned int iteration = 0;
    while (true)
      {
	const double residual = 0.5+erf(x/std::sqrt(2.)/sigma)/2-y;
	if (std::fabs(residual) < 1e-7)
	  break;
	const double F_prime = 1./std::sqrt(2*3.1415926536)/sigma *
			       std::exp(-x*x/sigma/sigma/2);
	x += -residual / F_prime;

					 // make sure that we don't
					 // recurse endlessly
	++iteration;
	Assert (iteration < 20, ExcInternalError());
      };
    return x+a;
  }



  std::vector<unsigned int>
  reverse_permutation (const std::vector<unsigned int> &permutation)
  {
    const unsigned int n = permutation.size();

    std::vector<unsigned int> out (n);
    for (unsigned int i=0; i<n; ++i)
      out[i] = n - 1 - permutation[i];

    return out;
  }



  std::vector<unsigned int>
  invert_permutation (const std::vector<unsigned int> &permutation)
  {
    const unsigned int n = permutation.size();

    std::vector<unsigned int> out (n, numbers::invalid_unsigned_int);

    for (unsigned int i=0; i<n; ++i)
      {
	Assert (permutation[i] < n, ExcIndexRange (permutation[i], 0, n));
	out[permutation[i]] = i;
      }

				     // check that we have actually reached
				     // all indices
    for (unsigned int i=0; i<n; ++i)
      Assert (out[i] != numbers::invalid_unsigned_int,
	      ExcMessage ("The given input permutation had duplicate entries!"));

    return out;
  }


  namespace MPI
  {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
                                // Unfortunately, we have to work
                                // around an oddity in the way PETSc
                                // and some gcc versions interact. If
                                // we use PETSc's MPI dummy
                                // implementation, it expands the
                                // calls to the two MPI functions
                                // basically as ``(n_jobs=1, 0)'',
                                // i.e. it assigns the number one to
                                // the variable holding the number of
                                // jobs, and then uses the comma
                                // operator to let the entire
                                // expression have the value zero. The
                                // latter is important, since
                                // ``MPI_Comm_size'' returns an error
                                // code that we may want to check (we
                                // don't here, but one could in
                                // principle), and the trick with the
                                // comma operator makes sure that both
                                // the number of jobs is correctly
                                // assigned, and the return value is
                                // zero. Unfortunately, if some recent
                                // versions of gcc detect that the
                                // comma expression just stands by
                                // itself, i.e. the result is not
                                // assigned to another variable, then
                                // they warn ``right-hand operand of
                                // comma has no effect''. This
                                // unwanted side effect can be
                                // suppressed by casting the result of
                                // the entire expression to type
                                // ``void'' -- not beautiful, but
                                // helps calming down unwarranted
                                // compiler warnings...
    unsigned int n_mpi_processes (const MPI_Comm &mpi_communicator)
    {
      int n_jobs=1;
      (void) MPI_Comm_size (mpi_communicator, &n_jobs);

      return n_jobs;
    }


    unsigned int this_mpi_process (const MPI_Comm &mpi_communicator)
    {
      int rank=0;
      (void) MPI_Comm_rank (mpi_communicator, &rank);

      return rank;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      MPI_Comm new_communicator;
      MPI_Comm_dup (mpi_communicator, &new_communicator);
      return new_communicator;
    }


    std::vector<unsigned int>
    compute_point_to_point_communication_pattern (const MPI_Comm & mpi_comm,
						  const std::vector<unsigned int> & destinations)
    {
      unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
      unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_comm);

      for (unsigned int i=0; i<destinations.size(); ++i)
	{
	  Assert (destinations[i] < n_procs,
		  ExcIndexRange (destinations[i], 0, n_procs));
	  Assert (destinations[i] != myid,
		  ExcMessage ("There is no point in communicating with ourselves."));
	}


				       // let all processors
				       // communicate the maximal
				       // number of destinations they
				       // have
      unsigned int my_n_destinations = destinations.size();
      unsigned int max_n_destinations = 0;

      MPI_Allreduce (&my_n_destinations, &max_n_destinations, 1, MPI_UNSIGNED,
		    MPI_MAX, mpi_comm);

				       // now that we know the number
				       // of data packets every
				       // processor wants to send, set
				       // up a buffer with the maximal
				       // size and copy our
				       // destinations in there,
				       // padded with -1's
      std::vector<unsigned int> my_destinations(max_n_destinations,
						numbers::invalid_unsigned_int);
      std::copy (destinations.begin(), destinations.end(),
		 my_destinations.begin());

				       // now exchange these (we could
				       // communicate less data if we
				       // used MPI_Allgatherv, but
				       // we'd have to communicate
				       // my_n_destinations to all
				       // processors in this case,
				       // which is more expensive than
				       // the reduction operation
				       // above in MPI_Allreduce)
      std::vector<unsigned int> all_destinations (max_n_destinations * n_procs);
      MPI_Allgather (&my_destinations[0], max_n_destinations, MPI_UNSIGNED,
		     &all_destinations[0], max_n_destinations, MPI_UNSIGNED,
		     mpi_comm);

				       // now we know who is going to
				       // communicate with
				       // whom. collect who is going
				       // to communicate with us!
      std::vector<unsigned int> origins;
      for (unsigned int i=0; i<n_procs; ++i)
	for (unsigned int j=0; j<max_n_destinations; ++j)
	  if (all_destinations[i*max_n_destinations + j] == myid)
	    origins.push_back (i);
	  else if (all_destinations[i*max_n_destinations + j] ==
		   numbers::invalid_unsigned_int)
	    break;

      return origins;
    }


    namespace
    {
				       // custom MIP_Op for
				       // calculate_collective_mpi_min_max_avg
      void max_reduce ( const void * in_lhs_,
			void * inout_rhs_,
			int * len,
			MPI_Datatype * )
      {
	const MinMaxAvg * in_lhs = static_cast<const MinMaxAvg*>(in_lhs_);
	MinMaxAvg * inout_rhs = static_cast<MinMaxAvg*>(inout_rhs_);

	Assert(*len==1, ExcInternalError());

	inout_rhs->sum += in_lhs->sum;
	if (inout_rhs->min>in_lhs->min)
	  {
	    inout_rhs->min = in_lhs->min;
	    inout_rhs->min_index = in_lhs->min_index;
	  }
	else if (inout_rhs->min == in_lhs->min)
	  { // choose lower cpu index when tied to make operator cumutative
	    if (inout_rhs->min_index > in_lhs->min_index)
	      inout_rhs->min_index = in_lhs->min_index;
	  }

	if (inout_rhs->max < in_lhs->max)
	  {
	  inout_rhs->max = in_lhs->max;
	  inout_rhs->max_index = in_lhs->max_index;
	  }
	else if (inout_rhs->max == in_lhs->max)
	  { // choose lower cpu index when tied to make operator cumutative
	    if (inout_rhs->max_index > in_lhs->max_index)
	      inout_rhs->max_index = in_lhs->max_index;
	  }
      }
    }



    MinMaxAvg
    min_max_avg(const double my_value,
		const MPI_Comm &mpi_communicator)
    {
      MinMaxAvg result;

      const unsigned int my_id
	= dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
      const unsigned int numproc
	= dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

      MPI_Op op;
      int ierr = MPI_Op_create((MPI_User_function *)&max_reduce, true, &op);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      MinMaxAvg in;
      in.sum = in.min = in.max = my_value;
      in.min_index = in.max_index = my_id;

      MPI_Datatype type;
      int lengths[]={3,2};
      MPI_Aint displacements[]={0,offsetof(MinMaxAvg, min_index)};
      MPI_Datatype types[]={MPI_DOUBLE, MPI_INT};

      ierr = MPI_Type_struct(2, lengths, displacements, types, &type);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Type_commit(&type);
      ierr = MPI_Allreduce ( &in, &result, 1, type, op, mpi_communicator );
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Type_free (&type);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      ierr = MPI_Op_free(&op);
      AssertThrow(ierr == MPI_SUCCESS, ExcInternalError());

      result.avg = result.sum / numproc;

      return result;
    }

#else

    unsigned int get_n_mpi_processes (const MPI_Comm &)
    {
      return 1;
    }



    unsigned int get_this_mpi_process (const MPI_Comm &)
    {
      return 0;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      return mpi_communicator;
    }




    void min_max_avg(const double my_value,
		     const MPI_Comm &)
    {
      MinMaxAvg result;

      result.sum = my_value;
      result.avg = my_value;
      result.min = my_value;
      result.max = my_value;
      result.min_index = 0;
      result.max_index = 0;

      return result;
    }

#endif
  }


  namespace System
  {
#if defined(__linux__)

    double get_cpu_load ()
    {
      std::ifstream cpuinfo;
      cpuinfo.open("/proc/loadavg");

      AssertThrow(cpuinfo, ExcIO());

      double load;
      cpuinfo >> load;

      return load;
    }

#else

    double get_cpu_load ()
    {
      return 0.;
    }

#endif



    void get_memory_stats (MemoryStats & stats)
    {
      stats.VmPeak = stats.VmSize = stats.VmHWM = stats.VmRSS = 0;

				       // parsing /proc/self/stat would be a
				       // lot easier, but it does not contain
				       // VmHWM, so we use /status instead.
#if defined(__linux__)
      std::ifstream file("/proc/self/status");
      std::string line;
      std::string name;
      while (!file.eof())
	{
	  file >> name;
	  if (name == "VmPeak:")
	    file >> stats.VmPeak;
	  else if (name == "VmSize:")
	    file >> stats.VmSize;
	  else if (name == "VmHWM:")
	    file >> stats.VmHWM;
	  else if (name == "VmRSS:")
	    {
	      file >> stats.VmRSS;
	      break; //this is always the last entry
	    }

	  getline(file, line);
	}
#endif
    }



    std::string get_hostname ()
    {
      const unsigned int N=1024;
      char hostname[N];
      gethostname (&(hostname[0]), N-1);
      return hostname;
    }



    std::string get_time ()
    {
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1);

      std::ostringstream o;
      o << time->tm_hour << ":"
        << (time->tm_min < 10 ? "0" : "") << time->tm_min << ":"
        << (time->tm_sec < 10 ? "0" : "") << time->tm_sec;

      return o.str();
    }


#ifdef DEAL_II_COMPILER_SUPPORTS_MPI

    bool job_supports_mpi ()
    {
      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);

      return (MPI_has_been_started > 0);
    }


#else

    bool job_supports_mpi ()
    {
      return false;
    }

#endif



    MPI_InitFinalize::MPI_InitFinalize (int    &argc,
					char** &argv)
		    :
		    owns_mpi (true)
    {
      static bool constructor_has_already_run = false;
      Assert (constructor_has_already_run == false,
	      ExcMessage ("You can only create a single object of this class "
			  "in a program since it initializes the MPI system."));

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);
      AssertThrow (MPI_has_been_started == 0,
		   ExcMessage ("MPI error. You can only start MPI once!"));

      int mpi_err;
      mpi_err = MPI_Init (&argc, &argv);
      AssertThrow (mpi_err == 0,
		   ExcMessage ("MPI could not be initialized."));
#else
				       // make sure the compiler doesn't warn
				       // about these variables
      (void)argc;
      (void)argv;
#endif

      constructor_has_already_run = true;
    }


    MPI_InitFinalize::~MPI_InitFinalize()
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI

#  if defined(DEAL_II_USE_TRILINOS) && !defined(__APPLE__)
				       // make memory pool release all
				       // vectors that are no longer
				       // used at this point. this is
				       // relevant because the static
				       // object destructors run for
				       // these vectors at the end of
				       // the program would run after
				       // MPI_Finalize is called,
				       // leading to errors
				       //
				       // TODO: On Mac OS X, shared libs can
 				       // only depend on other libs listed
 				       // later on the command line. This
 				       // means that libbase can't depend on
 				       // liblac, and we can't destroy the
 				       // memory pool here as long as we have
 				       // separate libraries. Consequently,
 				       // the #ifdef above. Deal will then
 				       // just continue to seg fault upon
 				       // completion of main()
      GrowingVectorMemory<TrilinosWrappers::MPI::Vector>
	::release_unused_memory ();
      GrowingVectorMemory<TrilinosWrappers::MPI::BlockVector>
	::release_unused_memory ();
#  endif

      int mpi_err = 0;

      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);
      if (program_uses_mpi() == true && owns_mpi == true &&
	  MPI_has_been_started != 0)
	{
	  if (std::uncaught_exception())
	    {
	      std::cerr << "ERROR: Uncaught exception in MPI_InitFinalize on proc "
			<< get_this_mpi_process(MPI_COMM_WORLD)
			<< ". Skipping MPI_Finalize() to avoid a deadlock."
			<< std::endl;
	    }
	  else
	    mpi_err = MPI_Finalize();
	}


      AssertThrow (mpi_err == 0,
		   ExcMessage ("An error occurred while calling MPI_Finalize()"));
#endif
    }



    bool
    program_uses_mpi ()
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      int MPI_has_been_started = 0;
      MPI_Initialized(&MPI_has_been_started);

      return true && (MPI_has_been_started > 0);
#else
      return false;
#endif
    }


    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator)
    {
      return MPI::n_mpi_processes (mpi_communicator);
    }

    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator)
    {
      return MPI::this_mpi_process (mpi_communicator);
    }



    void calculate_collective_mpi_min_max_avg(const MPI_Comm &mpi_communicator,
					      double my_value,
					      MinMaxAvg & result)
    {
      result = Utilities::MPI::min_max_avg (my_value,
					    mpi_communicator);
    }


  }


#ifdef DEAL_II_USE_TRILINOS

  namespace Trilinos
  {
    const Epetra_Comm&
    comm_world()
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      static Teuchos::RCP<Epetra_MpiComm>
	communicator = Teuchos::rcp (new Epetra_MpiComm (MPI_COMM_WORLD), true);
#else
      static Teuchos::RCP<Epetra_SerialComm>
	communicator = Teuchos::rcp (new Epetra_SerialComm (), true);
#endif

      return *communicator;
    }



    const Epetra_Comm&
    comm_self()
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      static Teuchos::RCP<Epetra_MpiComm>
	communicator = Teuchos::rcp (new Epetra_MpiComm (MPI_COMM_SELF), true);
#else
      static Teuchos::RCP<Epetra_SerialComm>
	communicator = Teuchos::rcp (new Epetra_SerialComm (), true);
#endif

      return *communicator;
    }



    Epetra_Comm *
    duplicate_communicator (const Epetra_Comm &communicator)
    {
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI

				       // see if the communicator is in fact a
				       // parallel MPI communicator; if so,
				       // return a duplicate of it
      const Epetra_MpiComm
	*mpi_comm = dynamic_cast<const Epetra_MpiComm *>(&communicator);
      if (mpi_comm != 0)
	return new Epetra_MpiComm(Utilities::System::
				  duplicate_communicator(mpi_comm->GetMpiComm()));
#endif

				       // if we don't support MPI, or if the
				       // communicator in question was in fact
				       // not an MPI communicator, return a
				       // copy of the same object again
      Assert (dynamic_cast<const Epetra_SerialComm*>(&communicator)
	      != 0,
	      ExcInternalError());
      return new Epetra_SerialComm(dynamic_cast<const Epetra_SerialComm&>(communicator));
    }



    void destroy_communicator (Epetra_Comm &communicator)
    {
      Assert (&communicator != 0, ExcInternalError());

				       // save the communicator, reset
				       // the map, and delete the
				       // communicator if this whole
				       // thing was created as an MPI
				       // communicator
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
      Epetra_MpiComm
	*mpi_comm = dynamic_cast<Epetra_MpiComm *>(&communicator);
      if (mpi_comm != 0)
	{
	  MPI_Comm comm = mpi_comm->GetMpiComm();
	  *mpi_comm = Epetra_MpiComm(MPI_COMM_SELF);
	  MPI_Comm_free (&comm);
	}
#endif
    }



    unsigned int get_n_mpi_processes (const Epetra_Comm &mpi_communicator)
    {
      return mpi_communicator.NumProc();
    }


    unsigned int get_this_mpi_process (const Epetra_Comm &mpi_communicator)
    {
      return (unsigned int)mpi_communicator.MyPID();
    }



    Epetra_Map
    duplicate_map (const Epetra_BlockMap &map,
		   const Epetra_Comm     &comm)
    {
      if (map.LinearMap() == true)
	{
					   // each processor stores a
					   // contiguous range of
					   // elements in the
					   // following constructor
					   // call
	  return Epetra_Map (map.NumGlobalElements(),
			     map.NumMyElements(),
			     map.IndexBase(),
			     comm);
	}
      else
	{
					   // the range is not
					   // contiguous
	  return Epetra_Map (map.NumGlobalElements(),
			     map.NumMyElements(),
			     map.MyGlobalElements (),
			     0,
			     comm);
	}
    }
  }

#endif

}

DEAL_II_NAMESPACE_CLOSE
