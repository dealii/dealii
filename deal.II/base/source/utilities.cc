//---------------------------------------------------------------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/utilities.h>
#include <base/exceptions.h>

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

#ifdef HAVE_STD_NUMERIC_LIMITS
#  include <limits>
#else
#  include <limits.h>
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
    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator)
    {
      int n_jobs=1;
      (void) MPI_Comm_size (mpi_communicator, &n_jobs);
      
      return n_jobs;
    }

    
    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator)
    {
      int rank=0;
      (void) MPI_Comm_rank (mpi_communicator, &rank);
      
      return rank;
    }


    MPI_Comm duplicate_communicator (const MPI_Comm &mpi_communicator)
    {
      MPI_Comm new_communicator;
      MPI_Comm_Dup (mpi_communicator, &new_communicator);
      return new_communicator;
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
#endif



    MPI_InitFinalize::MPI_InitFinalize (int    &argc,
					char** &argv)
		    :
		    owns_mpi (true)
    {
      static bool constructor_has_already_run = false;
      Assert (constructor_has_already_run == false,
	      ExcMessage ("You can only create a single object of this class "
			  "in a program sinc it initialize the MPI system."));
      
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
      int mpi_err = 0;

      if (program_uses_mpi() == true && owns_mpi == true)
	mpi_err = MPI_Finalize();

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


    unsigned int get_n_mpi_processes (const Epetra_Comm &mpi_communicator)
    {
      return mpi_communicator.NumProc();
    }
    

    unsigned int get_this_mpi_process (const Epetra_Comm &mpi_communicator)
    {
      return (unsigned int)mpi_communicator.MyPID();
    }
  }
  
#endif

}

DEAL_II_NAMESPACE_CLOSE
