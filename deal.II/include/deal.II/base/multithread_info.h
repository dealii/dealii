//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__multithread_info_h
#define __deal2__multithread_info_h
//---------------------------------------------------------------------------


#include <base/config.h>
#include <base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides information about the system which may be of
 * use in multithreaded programs.  At the moment this is just the
 * number of cpus. If deal.II is compiled with multithreading support,
 * some functions will use multiple threads for their action. Currently
 * the library supports both thread-based and task-based parallelism. @ref threads
 * describes the different uses of each. Thread-based parallel methods will
 * use the member variable n_default_threads of this class as a guide to the
 * number of threads to start.  This variable
 * n_default_threads is set to the number of CPUs by default, but
 * can be adjusted by the user to fit the requirements. The default number of 
 * threads used for task-based parallel methods is selected automatically
 * by the Threading Building
 * Blocks library. See @ref threads for more information on this.
 *
 * @ingroup threads
 * @author Thomas Richter, Wolfgang Bangerth, 2000
 */
class MultithreadInfo
{
  public:
				     /**
				      * The constructor determines the
				      * number of CPUs in the system.
				      * At the moment detection of
				      * CPUs is only implemented on
				      * Linux computers with the /proc
				      * filesystem and on Sun
				      * machines.  The number of CPUs
				      * present is set to one if
				      * detection failed or if
				      * detection is not supported.
				      */
    MultithreadInfo ();

				     /**
				      * The number of CPUs in the
				      * system.  It is one if
				      * detection is not implemented
				      * or failed.
				      *
				      * If it is one, although you
				      * are on a multi-processor
				      * machine, please refer to the
				      * documentation in
				      * <tt>multithread_info.cc</tt>
				      * near to the <tt>error</tt> directive.
				      */
    const unsigned int n_cpus;

				     /**
				      * The number of threads to use as
				      * a default value for all functions
				      * that support multithreading.
				      * At start time this is <tt>n_cpus</tt> or
				      * one, if detection of the number
				      * of CPUs is not possible.
				      */
    unsigned int n_default_threads;

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    static std::size_t memory_consumption ();

				     /**
				      * Exception
				      */
    DeclException0(ExcProcNotPresent);
    
  private:

				     /**
				      * Private function to determine
				      * the number of CPUs.
				      * Implementation for Linux, OSF,
				      * SGI, and Sun machines; if no
				      * detection of the number of CPUs is
				      * supported, or if detection
				      * fails, this function returns
				      * one.
				      */
    static unsigned int get_n_cpus();
};



/**
 * Global variable of type <tt>MultithreadInfo</tt> which you may ask for the
 * number of CPUs in your system, as well as for the default number of
 * threads that multithreaded functions shall use.
 *
 * @ingroup threads
 */
extern MultithreadInfo multithread_info;




//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef __deal2__multithread_info_h
#endif
//---------------------------------------------------------------------------
