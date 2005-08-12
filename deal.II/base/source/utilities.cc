//---------------------------------------------------------------------------
//      $Id$   
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
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
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cerrno>
#include <cmath>
#include <unistd.h>
#include <sys/types.h>
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif
#ifdef HAVE_STD_NUMERIC_LIMITS
# include <limits>
#endif


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
    AssertThrow ( ! ((digits==1 && i>=10)   ||
		     (digits==2 && i>=100)  ||
		     (digits==3 && i>=1000) ||
		     (digits==4 && i>=10000)||
		     (i>=100000)),
		  ExcInvalidNumber2StringConversersion(i, digits));
  
    std::string s;
    switch (digits) 
      {
	case 4:
	      s += '0' + i/1000;
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
#ifdef HAVE_STD_STRINGSTREAM
    std::istringstream ss(s);
#else
    std::ostrstream ss(s.c_str());
#endif

    int i = std::numeric_limits<int>::max();
    ss >> i;

                                     // check for errors
    AssertThrow (i != std::numeric_limits<int>::max(),
                 ExcCantConvertString (s));
    
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

#ifdef HAVE_STD_STRINGSTREAM
      std::ostringstream o;
#else
      std::ostrstream o;
#endif
      o << time->tm_hour << ":"
        << (time->tm_min < 10 ? "0" : "") << time->tm_min << ":"
        << (time->tm_sec < 10 ? "0" : "") << time->tm_sec;
#ifndef HAVE_STD_STRINGSTREAM
      o << std::ends;
#endif
      return o.str();
    }
    
  }
  
}
