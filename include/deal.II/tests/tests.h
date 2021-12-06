//-----------------------------------------------------------
//
//    Copyright (C) 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_tests_tests_h
#define dealii_tests_tests_h

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN


namespace Testing
{
  // given the name of a file, copy it to deallog
  // and then delete it
  void
  cat_file(const char *filename)
  {
    {
      std::ifstream in(filename);
      Assert(in, dealii::ExcIO());
      deallog.get_file_stream() << in.rdbuf() << "\n";
    }

    std::remove(filename);
  }



  // Function to initialize deallog. Normally, it should be called at
  // the beginning of main() like
  //
  // initlog();
  //
  // This will open the correct output file, divert log output there and
  // switch off screen output. If screen output is desired, provide the
  // optional first argument as 'true'.
  std::string   deallogname;
  std::ofstream deallogfile;

  void
  initlog(const bool                    console = false,
          const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                                std::ios::left)
  {
    deallogname = "output";
    deallogfile.open(deallogname);
    deallog.attach(deallogfile, true, flags);
    deallog.depth_console(console ? 10 : 0);
  }



  inline void
  mpi_initlog(const bool                    console = false,
              const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                                    std::ios::left)
  {
#ifdef DEAL_II_WITH_MPI
    unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    if (myid == 0)
      {
        deallogname = "output";
        deallogfile.open(deallogname.c_str());
        deallog.attach(deallogfile, true, flags);
        deallog.depth_console(console ? 10 : 0);
      }
#else
    (void)console;
    (void)flags;
    // can't use this function if not using MPI
    Assert(false, ExcInternalError());
#endif
  }



  /**
   * A helper class that gives each MPI process its own output file
   * for the `deallog` stream, and at the end of the program (or,
   * more correctly, the end of the current object), concatenates them
   * all into the output file used on processor 0.
   */
  struct MPILogInitAll
  {
    MPILogInitAll(const bool                    console = false,
                  const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                                        std::ios::left)
    {
#ifdef DEAL_II_WITH_MPI
      const unsigned int myid =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
#else
      constexpr unsigned int myid = 0;
#endif
      if (myid == 0)
        {
          if (!deallog.has_file())
            {
              deallogfile.open("output");
              deallog.attach(deallogfile, true, flags);
            }
        }
      else
        {
          deallogname = "output" + Utilities::int_to_string(myid);
          deallogfile.open(deallogname.c_str());
          deallog.attach(deallogfile, true, flags);
        }

      deallog.depth_console(console ? 10 : 0);

      deallog.push(Utilities::int_to_string(myid));
    }

    ~MPILogInitAll()
    {
      // pop the prefix for the MPI rank of the current process
      deallog.pop();

#ifdef DEAL_II_WITH_MPI
      const unsigned int myid =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      const unsigned int nproc =
        Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

      if (myid != 0)
        {
          deallog.detach();
          deallogfile.close();
        }

      MPI_Barrier(MPI_COMM_WORLD);

#  ifdef DEAL_II_WITH_PETSC
      check_petsc_allocations();
      MPI_Barrier(MPI_COMM_WORLD);
#  endif

      if (myid == 0)
        {
          for (unsigned int i = 1; i < nproc; ++i)
            {
              std::string filename = "output" + Utilities::int_to_string(i);
              cat_file(filename.c_str());
            }
        }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  };

} // namespace Testing


DEAL_II_NAMESPACE_CLOSE


#endif
