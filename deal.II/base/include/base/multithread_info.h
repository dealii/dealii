//----------------------------   multithread_info.h     ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------   multithread_info.h     ---------------------------
#ifndef __deal2__multithread_info_h
#define __deal2__multithread_info_h
//----------------------------   multithread_info.h     ---------------------------


#include <base/exceptions.h>


/**
 * This class provides information about the system.
 * At the moment this is just the number of cpus. If deal is
 * compiled with multithreading support,  all functions with
 * use #n_default_threads# as the default number of threads
 * to start all functions.
 * This variable #n_default_threads# is set to the number of
 * cpus per default, but can be adjusted by the user to fit
 * the requirements.
 */
class MultithreadInfo {
  public:
				     /**
				      * The constructor determines the
				      * number of cpus in the system.
				      * At the moment detection of cpus
				      * is only impelented on Linux
				      * computers with the proc filesystem
				      * and on suns.
				      * Other platforms will hopefully follow.
				      * The number of cpus present is set to
				      * one if detection failed or if
				      * detection is not supported.
				      */
    MultithreadInfo();

				     /**
				      * The number of cpus in the system.
				      * It is one, if detection is not implemented.
				      */
    const unsigned int n_cpus;

				     /**
				      * The number of threads to use as
				      * a default value for all functions
				      * that support multithreading.
				      * At start time this is #n_cpus# or
				      * one, if detection of the number
				      * of cpus is not possibly
				      */
    unsigned int n_default_threads;

				     /**
				      * Exception
				      */
    DeclException0(ExcProcNotPresent);
    
    
  private:

				     /**
				      * Private function to determine
				      * the number of cpus.
				      * Implementation for Linux and suns.
				      */
    static unsigned int get_n_cpus();
};


				 
extern MultithreadInfo multithread_info;




//----------------------------   multithread_info.h     ---------------------------
// end of #ifndef __deal2__multithread_info_h
#endif
//----------------------------   multithread_info.h     ---------------------------



