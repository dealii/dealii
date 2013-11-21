// ---------------------------------------------------------------------
// $Id$
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __tests_tests_h
#define __tests_tests_h

// common definitions used in all the tests

#include <deal.II/base/config.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/thread_management.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>


// Cygwin has a different implementation for rand() which causes many tests to fail.
// This here is a reimplementation that gives the same sequence of numbers as a program
// that uses rand() on a typical linux machine.
int myrand(bool reseed=false, int seed=1)
{
  static int r[32];
  static int k;
  static bool inited=false;
  if (!inited || reseed)
    {
      //srand treats a seed 0 as 1 for some reason
      r[0]=(seed==0)?1:seed;

      for (int i=1; i<31; i++)
        {
          r[i] = (16807LL * r[i-1]) % 2147483647;
          if (r[i] < 0)
            r[i] += 2147483647;
        }
      k=31;
      for (int i=31; i<34; i++)
        {
          r[k%32] = r[(k+32-31)%32];
          k=(k+1)%32;
        }

      for (int i=34; i<344; i++)
        {
          r[k%32] = r[(k+32-31)%32] + r[(k+32-3)%32];
          k=(k+1)%32;
        }
      inited=true;
      if (reseed==true)
        return 0;// do not generate new no
    }

  r[k%32] = r[(k+32-31)%32] + r[(k+32-3)%32];
  int ret = r[k%32];
  k=(k+1)%32;
  return (unsigned int)ret >> 1;
}

// reseed our random number generator
void mysrand(int seed)
{
  myrand(true, seed);
}


// I hope this is safe to just override the standard
// rand() functions here.
int rand()
{
  return myrand();
}

void srand(int seed)
{
  mysrand(seed);
}



// given the name of a file, copy it to deallog
// and then delete it
inline void cat_file(const char *filename)
{
  std::ifstream in(filename);
  Assert (in, dealii::ExcIO());

  while (in)
    {
      std::string s;
      std::getline(in, s);
      dealii::deallog.get_file_stream() << s << "\n";
    }
  in.close();

  std::remove (filename);
}


#ifdef DEAL_II_WITH_PETSC
#include <petscsys.h>

inline void check_petsc_allocations()
{
  PetscStageLog stageLog;
  PetscLogGetStageLog(&stageLog);

  // I don't quite understand petsc and it looks like stageLog->stageInfo->classLog->classInfo[i].id
  // is always -1, so we look it up in stageLog->classLog, make sure it has the same number of entries:
  Assert(stageLog->stageInfo->classLog->numClasses == stageLog->classLog->numClasses, dealii::ExcInternalError());

  bool errors = false;
  for (int i=0;i<stageLog->stageInfo->classLog->numClasses;++i)
    {
      if (stageLog->stageInfo->classLog->classInfo[i].destructions
	  != stageLog->stageInfo->classLog->classInfo[i].creations)
	{
	  errors = true;
	  std::cerr << "ERROR: PETSc objects leaking of type '"
	       <<  stageLog->classLog->classInfo[i].name << "'"
	       << " with " << stageLog->stageInfo->classLog->classInfo[i].creations
	       << " creations and only "
	       << stageLog->stageInfo->classLog->classInfo[i].destructions
	       << " destructions." << std::endl;
	  }
    }

  if (errors)
    throw dealii::ExcMessage("PETSc memory leak");
}
#endif


// implicitly use the deal.II namespace everywhere, without us having to say
// so in each and every testcase
using namespace dealii;

// Function for initialize deallog. Normally, it should be called at
// the beginning of main() like
//
// initlog();
//
// This will open the correct output file, divert log output there and
// switch off screen output. If screen output is desired, provide the
// optional second argument as 'true'.
std::string deallogname;
std::ofstream deallogfile;

inline
void
initlog(bool console=false)
{
  deallogname = "output";
  deallogfile.open(deallogname.c_str());
  deallog.attach(deallogfile);
  if (!console)
    deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
  deallog.threshold_float(1.e-8);
}


inline
void
mpi_initlog(bool console=false)
{
#ifdef DEAL_II_WITH_MPI
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  if (myid == 0)
    {
      deallogname = "output";
      deallogfile.open(deallogname.c_str());
      deallog.attach(deallogfile);
      if (!console)
        deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
      deallog.threshold_float(1.e-8);
    }
#else
  // can't use this function if not using MPI
  Assert (false, ExcInternalError());
#endif
}



/* helper class to include the deallogs of all processors
   on proc 0 */
struct MPILogInitAll
{
  MPILogInitAll(bool console=false)
  {
#ifdef DEAL_II_WITH_MPI
    unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    deallogname = "output";
    if (myid != 0)
      deallogname = deallogname + Utilities::int_to_string(myid);
    deallogfile.open(deallogname.c_str());
    deallog.attach(deallogfile);
    if (!console)
      deallog.depth_console(0);

//TODO: Remove this line and replace by test_mode()
    deallog.threshold_float(1.e-8);
    deallog.push(Utilities::int_to_string(myid));
#else
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }

  ~MPILogInitAll()
  {
#ifdef DEAL_II_WITH_MPI
    unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
    unsigned int nproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

    deallog.pop();

    if (myid!=0)
      {
        deallog.detach();
        deallogfile.close();
      }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEAL_II_WITH_PETSC
    check_petsc_allocations();
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (myid==0)
      {
        for (unsigned int i=1; i<nproc; ++i)
          {
            std::string filename = "output" + Utilities::int_to_string(i);
            std::ifstream in(filename.c_str());
            Assert (in, ExcIO());

            while (in)
              {
                std::string s;
                std::getline(in, s);
                deallog.get_file_stream() << s << "\n";
              }
            in.close();
            std::remove (filename.c_str());
          }
      }
#else
    // can't use this function if not using MPI
    Assert (false, ExcInternalError());
#endif
  }

};



/*
 * Some tests (notably base/thread*, base/task*) create output that
 * comes out in random order. To make the output of these tests comparable,
 * we need to sort them.
 *
 * This function does just that with the file given. All streams writing
 * to this should be closed when calling this function.
 */
void sort_file_contents (const std::string &filename)
{
  int error = std::system ((std::string ("LC_ALL=C sort ") + filename + " -o " + filename).c_str());
  Assert (error == 0,
	  ExcInternalError());
}


/*
 * Replace all occurences of ' &' by '& ' from the given file to hopefully be more compiler
 * independent with respect to __PRETTY_FUNCTION__
 */
void unify_pretty_function (const std::string &filename)
{
  int error = std::system ((std::string ("sed -i -e 's/ \\&/ \\& /g' -e 's/ & ,/\\&,/g' -e 's/ \\& )/\\&)/g' -e 's/ \\& /\\& /g' ") + filename).c_str());
  Assert (error == 0,
	  ExcInternalError());
}





#ifndef DEAL_II_STACKTRACE_SWITCH
#define DEAL_II_STACKTRACE_SWITCH

// A structure and a variable that are used to make sure that we do
// not show a stacktrace in out testcases, since this would lead to
// trouble with automatic comparison of output against a fixed
// baseline
struct SwitchOffStacktrace
{
  SwitchOffStacktrace ()
  {
    deal_II_exceptions::suppress_stacktrace_in_exceptions ();
  }
} deal_II_stacktrace_dummy;

#endif

DEAL_II_NAMESPACE_OPEN
namespace internal
{

  namespace Vector
  {

    extern unsigned int minimum_parallel_grain_size;

  }

  namespace SparseMatrix
  {

    extern unsigned int minimum_parallel_grain_size;

  }

}

// A structure and a variable that are used to set grainsizes for parallel
// mode smaller than they would otherwise be. this is used to test that the
// parallel algorithms in lac/ work alright.
struct SetGrainSizes
{

  SetGrainSizes ()
  {

    internal::Vector::minimum_parallel_grain_size = 2;

    internal::SparseMatrix::minimum_parallel_grain_size = 2;

  }

}
set_grain_sizes;

DEAL_II_NAMESPACE_CLOSE

#endif // __tests_tests_h
