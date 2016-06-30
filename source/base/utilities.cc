// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_local_storage.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/math/special_functions/erf.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#ifndef DEAL_II_MSVC
#  include <stdlib.h>
#endif

#ifdef DEAL_II_MSVC
#  include <winsock2.h>
#endif


#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
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
  int_to_string (const unsigned int value, const unsigned int digits)
  {
    return to_string(value,digits);
  }

  template <typename number>
  std::string
  to_string (const number value, const unsigned int digits)
  {
    std::string lc_string = boost::lexical_cast<std::string>(value);

    if (digits == numbers::invalid_unsigned_int)
      return lc_string;
    else if (lc_string.size() < digits)
      {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                              ?
                                              1
                                              :
                                              0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
      }
    return lc_string;
  }


  std::string
  replace_in_string(const std::string &input,
                    const std::string &from,
                    const std::string &to)
  {
    if (from.empty())
      return input;

    std::string out = input;
    std::string::size_type pos = out.find(from);

    while (pos != std::string::npos)
      {
        out.replace(pos, from.size(), to);
        pos = out.find(from, pos + to.size());
      }
    return out;
  }

  std::string
  trim(const std::string &input)
  {
    std::string::size_type left = 0;
    std::string::size_type right = input.size() - 1;

    for (; left < input.size(); ++left)
      {
        if (!std::isspace(input[left]))
          {
            break;
          }
      }

    for (; right >= left; --right)
      {
        if (!std::isspace(input[right]))
          {
            break;
          }
      }

    return std::string(input, left, right - left + 1);
  }



  std::string
  dim_string(const int dim, const int spacedim)
  {
    if (dim == spacedim)
      return int_to_string(dim);
    else
      return int_to_string(dim)+","+int_to_string(spacedim);
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
  string_to_int (const std::string &s_)
  {
    // trim whitespace on either side of the text if necessary
    std::string s = s_;
    while ((s.size() > 0) && (s[0] == ' '))
      s.erase (s.begin());
    while ((s.size() > 0) && (s[s.size()-1] == ' '))
      s.erase (s.end()-1);

    // now convert and see whether we succeed. note that strtol only
    // touches errno if an error occurred, so if we want to check
    // whether an error happened, we need to make sure that errno==0
    // before calling strtol since otherwise it may be that the
    // conversion succeeds and that errno remains at the value it
    // was before, whatever that was
    char *p;
    errno = 0;
    const int i = std::strtol(s.c_str(), &p, 10);
    AssertThrow ( !((errno != 0) || (s.size() == 0) || ((s.size()>0) && (*p != '\0'))),
                  ExcMessage ("Can't convert <" + s + "> to an integer."));

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
  string_to_double (const std::string &s_)
  {
    // trim whitespace on either side of the text if necessary
    std::string s = s_;
    while ((s.size() > 0) && (s[0] == ' '))
      s.erase (s.begin());
    while ((s.size() > 0) && (s[s.size()-1] == ' '))
      s.erase (s.end()-1);

    // now convert and see whether we succeed. note that strtol only
    // touches errno if an error occurred, so if we want to check
    // whether an error happened, we need to make sure that errno==0
    // before calling strtol since otherwise it may be that the
    // conversion succeeds and that errno remains at the value it
    // was before, whatever that was
    char *p;
    errno = 0;
    const double d = std::strtod(s.c_str(), &p);
    AssertThrow ( !((errno != 0) || (s.size() == 0) || ((s.size()>0) && (*p != '\0'))),
                  ExcMessage ("Can't convert <" + s + "> to a double."));

    return d;
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
    // keep the currently remaining part of the input string in 'tmp' and
    // keep chopping elements of the list off the front
    std::string tmp = s;

    // as discussed in the documentation, eat whitespace from the end
    // of the string
    while (tmp.length() != 0 && tmp[tmp.length()-1] == ' ')
      tmp.erase (tmp.length()-1, 1);

    // split the input list until it is empty. since in every iteration
    // 'tmp' is what's left of the string after the next delimiter,
    // and since we've stripped trailing space already, 'tmp' will
    // be empty at one point if 's' ended in a delimiter, even if
    // there was space after the last delimiter. this matches what's
    // discussed in the documentation
    std::vector<std::string> split_list;
    split_list.reserve (std::count (tmp.begin(), tmp.end(), delimiter)+1);
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

        // strip spaces from this element's front and end
        while ((name.length() != 0) && (name[0] == ' '))
          name.erase (0,1);
        while (name.length() != 0 && name[name.length()-1] == ' ')
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

        std::size_t pos_newline = text.find_first_of("\n", 0);
        if (pos_newline != std::string::npos && pos_newline <= width)
          {
            std::string line (text, 0, pos_newline);
            while ((line.length() != 0) && (line[line.length()-1] == delimiter))
              line.erase(line.length()-1,1);
            lines.push_back (line);
            text.erase (0, pos_newline+1);
            continue;
          }

        // if we can fit everything into one
        // line, then do so. otherwise, we have
        // to keep breaking
        if (text.length() < width)
          {
            // remove trailing spaces
            while ((text.length() != 0) && (text[text.length()-1] == delimiter))
              text.erase(text.length()-1,1);
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
            std::string line (text, 0, location);
            while ((line.length() != 0) && (line[line.length()-1] == delimiter))
              line.erase(line.length()-1,1);
            lines.push_back (line);
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

    // we would want to use rand(), but that function is not reentrant
    // in a thread context. one could use rand_r, but this does not
    // produce reproducible results between threads either (though at
    // least it is reentrant). these two approaches being
    // non-workable, use a thread-local random number generator here
    static Threads::ThreadLocalStorage<boost::mt19937> random_number_generator;
    return boost::normal_distribution<>(a,sigma)(random_number_generator.get());
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

  std::vector<unsigned long long int>
  reverse_permutation (const std::vector<unsigned long long int> &permutation)
  {
    const unsigned long long int n = permutation.size();

    std::vector<unsigned long long int> out (n);
    for (unsigned long long int i=0; i<n; ++i)
      out[i] = n - 1 - permutation[i];

    return out;
  }



  std::vector<unsigned long long int>
  invert_permutation (const std::vector<unsigned long long int> &permutation)
  {
    const unsigned long long int n = permutation.size();

    std::vector<unsigned long long int> out (n, numbers::invalid_unsigned_int);

    for (unsigned long long int i=0; i<n; ++i)
      {
        Assert (permutation[i] < n, ExcIndexRange (permutation[i], 0, n));
        out[permutation[i]] = i;
      }

    // check that we have actually reached
    // all indices
    for (unsigned long long int i=0; i<n; ++i)
      Assert (out[i] != numbers::invalid_unsigned_int,
              ExcMessage ("The given input permutation had duplicate entries!"));

    return out;
  }


  template <typename Integer>
  std::vector<Integer>
  reverse_permutation (const std::vector<Integer> &permutation)
  {
    const unsigned int n = permutation.size();

    std::vector<Integer> out (n);
    for (unsigned int i=0; i<n; ++i)
      out[i] = n - 1 - permutation[i];

    return out;
  }



  template <typename Integer>
  std::vector<Integer>
  invert_permutation (const std::vector<Integer> &permutation)
  {
    const unsigned int n = permutation.size();

    std::vector<Integer> out (n, numbers::invalid_unsigned_int);

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



    void get_memory_stats (MemoryStats &stats)
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
#if defined(DEAL_II_HAVE_UNISTD_H) && defined(DEAL_II_HAVE_GETHOSTNAME)
      const unsigned int N=1024;
      char hostname[N];
      gethostname (&(hostname[0]), N-1);
#else
      std::string hostname("unknown");
#endif
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



    std::string get_date ()
    {
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1);

      std::ostringstream o;
      o << time->tm_year + 1900 << "/"
        << time->tm_mon + 1 << "/"
        << time->tm_mday;

      return o.str();
    }



    void posix_memalign (void **memptr, size_t alignment, size_t size)
    {
#ifndef DEAL_II_MSVC
      const int ierr = ::posix_memalign (memptr, alignment, size);

      AssertThrow (ierr == 0, ExcOutOfMemory());
      AssertThrow (*memptr != 0, ExcOutOfMemory());
#else
      // Windows does not appear to have posix_memalign. just use the
      // regular malloc in that case
      *memptr = malloc (size);
      (void)alignment;
      AssertThrow (*memptr != 0, ExcOutOfMemory());
#endif
    }



    bool job_supports_mpi ()
    {
      return Utilities::MPI::job_supports_mpi();
    }
  }


#ifdef DEAL_II_WITH_TRILINOS

  namespace Trilinos
  {
    const Epetra_Comm &
    comm_world()
    {
#ifdef DEAL_II_WITH_MPI
      static Teuchos::RCP<Epetra_MpiComm>
      communicator = Teuchos::rcp (new Epetra_MpiComm (MPI_COMM_WORLD), true);
#else
      static Teuchos::RCP<Epetra_SerialComm>
      communicator = Teuchos::rcp (new Epetra_SerialComm (), true);
#endif

      return *communicator;
    }



    const Epetra_Comm &
    comm_self()
    {
#ifdef DEAL_II_WITH_MPI
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
#ifdef DEAL_II_WITH_MPI

      // see if the communicator is in fact a
      // parallel MPI communicator; if so,
      // return a duplicate of it
      const Epetra_MpiComm
      *mpi_comm = dynamic_cast<const Epetra_MpiComm *>(&communicator);
      if (mpi_comm != 0)
        return new Epetra_MpiComm(Utilities::MPI::
                                  duplicate_communicator(mpi_comm->GetMpiComm()));
#endif

      // if we don't support MPI, or if the
      // communicator in question was in fact
      // not an MPI communicator, return a
      // copy of the same object again
      Assert (dynamic_cast<const Epetra_SerialComm *>(&communicator)
              != 0,
              ExcInternalError());
      return new Epetra_SerialComm(dynamic_cast<const Epetra_SerialComm &>(communicator));
    }



    void destroy_communicator (Epetra_Comm &communicator)
    {
      // save the communicator, reset the map, and delete the communicator if
      // this whole thing was created as an MPI communicator
#ifdef DEAL_II_WITH_MPI
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

  template std::string to_string<int> (int, unsigned int);
  template std::string to_string<long int> (long int, unsigned int);
  template std::string to_string<long long int> (long long int, unsigned int);
  template std::string to_string<unsigned int> (unsigned int, unsigned int);
  template std::string to_string<unsigned long int> (unsigned long int, unsigned int);
  template std::string to_string<unsigned long long int> (unsigned long long int, unsigned int);
  template std::string to_string<float> (float, unsigned int);
  template std::string to_string<double> (double, unsigned int);
  template std::string to_string<long double> (long double, unsigned int);

}

DEAL_II_NAMESPACE_CLOSE
