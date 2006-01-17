//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup threads Multithreading
 *
 * On machines with more than one processor (or multicore processors),
 * it is often profitable to run several parts of the computations in
 * parallel. For example, one could have several threads running in
 * parallel, each of which assembles the cell matrices of a subset of
 * the triangulation and then writes them into the global
 * matrix. Since assembling matrices is often an expensive operation,
 * this frequently leads to significant savings in compute time on
 * multiprocessor machines.
 *
 * In a similar way, it is often also profitable to use multiple
 * threads on a single-CPU system if a significant amount of input or
 * output tasks has to be performed. In such cases, the program
 * usually has to wait for disks or network storages to provide the
 * requested data, or to flush buffers. If this is done on a separate
 * thread, other threads of the program can continue to do other, more
 * interesting things at the same time, using the CPU downtime.
 *
 * Support for this model of computations, i.e. using multiple threads
 * on a shared-memory machine (SMP machine) is provided mainly through
 * the Threads namespace that offers functions to create new threads
 * as well as synchronisation primitives. The MultithreadInfo class
 * allows to query certain properties of the system, such as the
 * number of CPUs.  The use of these classes is explained in the
 * step-9, step-13 and step-14 tutorial programs.
 *
 * On the other hand, programs running on distributed memory machines
 * (i.e. clusters) need a different programming model built on top of
 * MPI and PETSc that is described in the step-17 and later example
 * programs.
 */
