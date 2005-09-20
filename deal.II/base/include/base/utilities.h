//---------------------------------------------------------------------------
//    $Id$
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
#ifndef __deal2__utilities_h
#define __deal2__utilities_h

#include <base/config.h>

#include <vector>
#include <string>

#ifdef DEAL_II_USE_PETSC
#include <petscvec.h>
#endif    

/**
 * A namespace for utility functions that are not particularly specific to
 * finite element computing or numerical programs, but nevertheless are needed
 * in various contexts when writing applications.
 *
 * @author Wolfgang Bangerth, 2005
 */
namespace Utilities
{

                                   /**
                                    * Convert a number @p i to a string, with
                                    * as many digits as given to fill with
                                    * leading zeros.
                                    */
  std::string
  int_to_string (const unsigned int i,
		 const unsigned int digits);

                                   /**
                                    * Determine how many digits are needed to
                                    * represent numbers at most as large as
                                    * the given number.
                                    */
  unsigned int
  needed_digits (const unsigned int max_number);

                                   /**
                                    * Given a string, convert it to an
                                    * integer. Throw an assertion if that is
                                    * not possible.
                                    */
  int
  string_to_int (const std::string &s);


                                   /**
                                    * Given a list of strings, convert it to a
                                    * list of integers. Throw an assertion if
                                    * that is not possible.
                                    */
  std::vector<int>
  string_to_int (const std::vector<std::string> &s);

                                   /**
                                    * Given a string that contains text
                                    * separated by a @p delimiter, split it into
                                    * its components; for each component,
                                    * remove leading and trailing spaces.
                                    *
                                    * The default value of the delimiter is a
                                    * comma, so that the function splits comma
                                    * separated lists of strings.
                                    */
  std::vector<std::string>
  split_string_list (const std::string &s,
                     const char         delimiter = ',');

                                   /**
                                    * Take a text, usually a documentation or
                                    * something, and try to break it into
                                    * individual lines of text at most @p
                                    * width characters wide, by breaking at
                                    * spaces in the text. If this is not
                                    * possible, return the shortest lines than
                                    * are longer than @p width.
                                    */
  std::vector<std::string>
  break_text_into_lines (const std::string &original_text,
                         const unsigned int width);

  
                                   /**
                                    * Generate a random number from a
                                    * normalized Gaussian probability
                                    * distribution centered around @p a and
                                    * with standard deviation @p sigma.
                                    */
  double
  generate_normal_random_number (const double a,
				 const double sigma);


                                   /**
                                    * A namespace for utility functions that
                                    * probe system properties.
                                    */
  namespace System
  {
    
                                     /**
                                      * Return the CPU load as returned by
                                      * "uptime". Note that the interpretation
                                      * of this number depends on the actual
                                      * number of processors in the
                                      * machine. This is presently only
                                      * implemented on Linux, using the
                                      * /proc/loadavg pseudo-file, on other
                                      * systems we simply return zero.
                                      */
    double get_cpu_load ();


                                     /**
                                      * Return the name of the host this
                                      * process runs on.
                                      */
    std::string get_hostname ();


                                     /**
                                      * Return the present time as HH:MM:SS.
                                      */
    std::string get_time ();

#ifdef DEAL_II_USE_PETSC
				     /**
				      * Return the number of MPI
				      * processes there exist.
				      */
    unsigned int get_n_mpi_processes (const MPI_Comm &mpi_communicator);

				     /**
				      * Return the number of the
				      * present MPI process in the
				      * space of processes. This will
				      * be a unique value for each
				      * process between zero and (less
				      * than) the number of all
				      * processes (given by
				      * get_n_mpi_processes()).
				      */
    unsigned int get_this_mpi_process (const MPI_Comm &mpi_communicator);
#endif // DEAL_II_USE_PETSC    
  }
}


#endif

