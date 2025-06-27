// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

//
// The content of this file is a slightly modified version of the
// header-only library BigMPICompat - a tiny MPI 4.x compatibility
// library released under MIT license at
// https://github.com/tjhei/BigMPICompat and relicensed and included
// here with permission.
//

#ifndef BIG_MPI_COMPAT_H
#define BIG_MPI_COMPAT_H

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_MPI

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <mpi.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

// required for std::numeric_limits used below.
#  include <limits>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    /**
     * This namespace contains symbols to support MPI routines with
     * large "counts" on MPI implementations that implement version
     * 3.x of the standard, where @p count is a signed integer.
     */
    namespace LargeCount
    {
      /**
       * This is the largest @p count supported when it is represented
       * with a signed integer (old MPI routines).
       */
      inline constexpr MPI_Count mpi_max_int_count =
        std::numeric_limits<int>::max();

      /**
       * Create a contiguous type of (possibly large) @p count.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      Type_contiguous_c(MPI_Count     count,
                        MPI_Datatype  oldtype,
                        MPI_Datatype *newtype)
      {
#  if MPI_VERSION >= 4
        return MPI_Type_contiguous_c(count, oldtype, newtype);
#  else
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_Type_contiguous(count, oldtype, newtype);
        else
          {
            int             ierr;
            const MPI_Count max_signed_int = LargeCount::mpi_max_int_count;

            MPI_Count size_old;
            ierr = MPI_Type_size_x(oldtype, &size_old);

            MPI_Count n_chunks             = count / max_signed_int;
            MPI_Count n_remaining_elements = count % max_signed_int;

            MPI_Datatype chunks;
            ierr = MPI_Type_vector(
              n_chunks, max_signed_int, max_signed_int, oldtype, &chunks);
            if (ierr != MPI_SUCCESS)
              return ierr;

            MPI_Datatype remainder;
            ierr =
              MPI_Type_contiguous(n_remaining_elements, oldtype, &remainder);
            if (ierr != MPI_SUCCESS)
              return ierr;

            int          blocklengths[2]  = {1, 1};
            MPI_Aint     displacements[2] = {0,
                                             static_cast<MPI_Aint>(n_chunks) *
                                               size_old * max_signed_int};
            MPI_Datatype types[2]         = {chunks, remainder};
            ierr                          = MPI_Type_create_struct(
              2, blocklengths, displacements, types, newtype);
            if (ierr != MPI_SUCCESS)
              return ierr;

            ierr = MPI_Type_commit(newtype);
            if (ierr != MPI_SUCCESS)
              return ierr;

            ierr = MPI_Type_free(&chunks);
            if (ierr != MPI_SUCCESS)
              return ierr;

            ierr = MPI_Type_free(&remainder);
            if (ierr != MPI_SUCCESS)
              return ierr;

#    ifndef MPI_COMPAT_SKIP_SIZE_CHECK
            MPI_Count size_new;
            ierr = MPI_Type_size_x(*newtype, &size_new);
            if (ierr != MPI_SUCCESS)
              return ierr;

            if (size_old * count != size_new)
              {
                // This error can happen when you are using a very old and
                // buggy MPI implementation. There is nothing we can do
                // here, unfortunately. Please update your installation.
                //	  std::cerr
                //<< "MPI_Type_contiguous_c() produced an invalid result.
                // Expected =
                //"
                //<< size_old << " * " << count << " = " << size_old * count
                //<< " but received " << size_new << std::endl;
                return MPI_ERR_INTERN;
              }
#    endif

            return MPI_SUCCESS;
          }
#  endif
      }


      /**
       * Send a package to rank @p dest with a (possibly large) @p count.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      Send_c(const void  *buf,
             MPI_Count    count,
             MPI_Datatype datatype,
             int          dest,
             int          tag,
             MPI_Comm     comm)
      {
#  if MPI_VERSION >= 4
        return MPI_Send_c(buf, count, datatype, dest, tag, comm);
#  else
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_Send(buf, count, datatype, dest, tag, comm);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Send(buf, 1, bigtype, dest, tag, comm);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
#  endif
      }

      /**
       * Receive a package from rank @p source with a (possibly large) @p count.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      Recv_c(void        *buf,
             MPI_Count    count,
             MPI_Datatype datatype,
             int          source,
             int          tag,
             MPI_Comm     comm,
             MPI_Status  *status)
      {
#  if MPI_VERSION >= 4
        return MPI_Recv_c(buf, count, datatype, source, tag, comm, status);
#  else
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_Recv(buf, count, datatype, source, tag, comm, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Recv(buf, 1, bigtype, source, tag, comm, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
#  endif
      }

      /**
       * Broadcast a message of possibly large @p count of data from
       * the process with rank "root" to all other processes.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      Bcast_c(void        *buf,
              MPI_Count    count,
              MPI_Datatype datatype,
              unsigned int root_mpi_rank,
              MPI_Comm     comm)
      {
#  if MPI_VERSION >= 4
        return MPI_Bcast_c(buf, count, datatype, root_mpi_rank, comm);
#  else
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_Bcast(buf, count, datatype, root_mpi_rank, comm);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Bcast(buf, 1, bigtype, root_mpi_rank, comm);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
#  endif
      }

      /**
       * Write a possibly large @p count of data at the location @p offset.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      File_write_at_c(MPI_File     fh,
                      MPI_Offset   offset,
                      const void  *buf,
                      MPI_Count    count,
                      MPI_Datatype datatype,
                      MPI_Status  *status)
      {
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_File_write_at(fh, offset, buf, count, datatype, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_File_write_at(fh, offset, buf, 1, bigtype, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
      }

      /**
       * Collectively write a possibly large @p count of data at the
       * location @p offset.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      File_write_at_all_c(MPI_File     fh,
                          MPI_Offset   offset,
                          const void  *buf,
                          MPI_Count    count,
                          MPI_Datatype datatype,
                          MPI_Status  *status)
      {
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_File_write_at_all(
            fh, offset, buf, count, datatype, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_File_write_at_all(fh, offset, buf, 1, bigtype, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
      }

      /**
       * Collectively write a possibly large @p count of data in order.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      File_write_ordered_c(MPI_File     fh,
                           const void  *buf,
                           MPI_Count    count,
                           MPI_Datatype datatype,
                           MPI_Status  *status)
      {
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_File_write_ordered(fh, buf, count, datatype, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_File_write_ordered(fh, buf, 1, bigtype, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
      }

      /**
       * Read a possibly large @p count of data at the
       * location @p offset.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      File_read_at_c(MPI_File     fh,
                     MPI_Offset   offset,
                     void        *buf,
                     MPI_Count    count,
                     MPI_Datatype datatype,
                     MPI_Status  *status)
      {
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_File_read_at(fh, offset, buf, count, datatype, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_File_read_at(fh, offset, buf, 1, bigtype, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
      }

      /**
       * Collectively read a possibly large @p count of data at the
       * location @p offset.
       *
       * See the MPI 4.x standard for details.
       */
      inline int
      File_read_at_all_c(MPI_File     fh,
                         MPI_Offset   offset,
                         void        *buf,
                         MPI_Count    count,
                         MPI_Datatype datatype,
                         MPI_Status  *status)
      {
        if (count <= LargeCount::mpi_max_int_count)
          return MPI_File_read_at_all(fh, offset, buf, count, datatype, status);

        MPI_Datatype bigtype;
        int          ierr;
        ierr = Type_contiguous_c(count, datatype, &bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        ierr = MPI_Type_commit(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_File_read_at_all(fh, offset, buf, 1, bigtype, status);
        if (ierr != MPI_SUCCESS)
          return ierr;

        ierr = MPI_Type_free(&bigtype);
        if (ierr != MPI_SUCCESS)
          return ierr;
        return MPI_SUCCESS;
      }

    } // namespace LargeCount
  }   // namespace MPI
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
