// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/base/hdf5.h>

#  include <hdf5.h>

#  include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace HDF5
{
  namespace internal
  {
    namespace HDF5ObjectImplementation
    {
      /**
       * Abort, should there be an exception being processed (see the error
       * message).
       */
      void
      check_exception(const std::string &type, const std::string &name)
      {
#  ifdef DEAL_II_WITH_MPI
        if (std::uncaught_exceptions() > 0)
          {
            std::cerr
              << "---------------------------------------------------------\n"
              << "HDF5 " + type + " objects call " + name +
                   " to end access to\n"
              << "them and release any resources they acquired. This call\n"
              << "requires MPI synchronization. Since an exception is\n"
              << "currently uncaught, this synchronization would likely\n"
              << "deadlock because only the current process is trying to\n"
              << "destroy the object. As a consequence, the program will be\n"
              << "aborted and the HDF5 might be corrupted.\n"
              << "---------------------------------------------------------"
              << std::endl;

            MPI_Abort(MPI_COMM_WORLD, 1);
          }
#  else
        (void)type;
        (void)name;
#  endif
      }
    } // namespace HDF5ObjectImplementation
  }   // namespace internal



  HDF5Object::HDF5Object(const std::string &name, const bool mpi)
    : name(name)
    , mpi(mpi)
  {}



  std::string
  HDF5Object::get_name() const
  {
    return name;
  }



  DataSet::DataSet(const std::string &name,
                   const hid_t       &parent_group_id,
                   const bool         mpi)
    : HDF5Object(name, mpi)
    , query_io_mode(false)
    , io_mode(H5D_MPIO_NO_COLLECTIVE)
    , local_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
    , global_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("DataSet",
                                                          "H5Dclose");
      // Release the HDF5 resource
      const herr_t ret = H5Dclose(*pointer);
      AssertNothrow(ret >= 0, ExcInternalError());
      delete pointer;
    });
    dataspace      = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("DataSpace",
                                                          "H5Sclose");
      // Release the HDF5 resource
      const herr_t ret = H5Sclose(*pointer);
      AssertNothrow(ret >= 0, ExcInternalError());
      delete pointer;
    });

    // rank_ret can take a negative value if the functions
    // H5Sget_simple_extent_ndims and H5Sget_simple_extent_dims fail. rank is
    // unsigned int therefore it can not be used to store the return value of
    // H5Sget_simple_extent_ndims and H5Sget_simple_extent_dims
    int rank_ret;

    *hdf5_reference = H5Dopen2(parent_group_id, name.data(), H5P_DEFAULT);
    Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Dopen2"));
    *dataspace = H5Dget_space(*hdf5_reference);
    Assert(*dataspace >= 0, ExcMessage("Error at H5Dget_space"));
    rank_ret = H5Sget_simple_extent_ndims(*dataspace);
    Assert(rank_ret >= 0, ExcInternalError());
    rank = rank_ret;
    dimensions.resize(rank);
    rank_ret =
      H5Sget_simple_extent_dims(*dataspace, dimensions.data(), nullptr);
    AssertDimension(rank_ret, static_cast<int>(rank));

    size = 1;
    for (const auto &dimension : dimensions)
      {
        size *= dimension;
      }
  }



  DataSet::DataSet(const std::string            &name,
                   const hid_t                  &parent_group_id,
                   const std::vector<hsize_t>   &dimensions,
                   const std::shared_ptr<hid_t> &t_type,
                   const bool                    mpi)
    : HDF5Object(name, mpi)
    , rank(dimensions.size())
    , dimensions(dimensions)
    , query_io_mode(false)
    , io_mode(H5D_MPIO_NO_COLLECTIVE)
    , local_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
    , global_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("DataSet",
                                                          "H5Dclose");
      // Release the HDF5 resource
      const herr_t ret = H5Dclose(*pointer);
      AssertNothrow(ret >= 0, ExcInternalError());
      delete pointer;
    });
    dataspace      = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("DataSpace",
                                                          "H5Sclose");
      // Release the HDF5 resource
      const herr_t ret = H5Sclose(*pointer);
      AssertNothrow(ret >= 0, ExcInternalError());
      delete pointer;
    });

    *dataspace = H5Screate_simple(rank, dimensions.data(), nullptr);
    Assert(*dataspace >= 0, ExcMessage("Error at H5Screate_simple"));

    *hdf5_reference = H5Dcreate2(parent_group_id,
                                 name.data(),
                                 *t_type,
                                 *dataspace,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
    Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Dcreate2"));

    size = 1;
    for (const auto &dimension : dimensions)
      {
        size *= dimension;
      }
  }



  void
  DataSet::set_query_io_mode(const bool new_query_io_mode)
  {
    query_io_mode = new_query_io_mode;
  }



  std::vector<hsize_t>
  DataSet::get_dimensions() const
  {
    return dimensions;
  }



  std::string
  DataSet::get_io_mode()
  {
    Assert(query_io_mode,
           ExcMessage("query_io_mode should be true in order to use io_mode"));
    switch (io_mode)
      {
        case (H5D_MPIO_NO_COLLECTIVE):
          return "H5D_MPIO_NO_COLLECTIVE";
          break;
        case (H5D_MPIO_CHUNK_INDEPENDENT):
          return "H5D_MPIO_CHUNK_INDEPENDENT";
          break;
        case (H5D_MPIO_CHUNK_COLLECTIVE):
          return "H5D_MPIO_CHUNK_COLLECTIVE";
          break;
        case (H5D_MPIO_CHUNK_MIXED):
          return "H5D_MPIO_CHUNK_MIXED";
          break;
        case (H5D_MPIO_CONTIGUOUS_COLLECTIVE):
          return "H5D_MPIO_CONTIGUOUS_COLLECTIVE";
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          return "Internal error";
          break;
      }
    // The function should not reach this line.
    DEAL_II_ASSERT_UNREACHABLE();
    return "Internal error";
  }



  H5D_mpio_actual_io_mode_t
  DataSet::get_io_mode_as_hdf5_type()
  {
    Assert(query_io_mode,
           ExcMessage(
             "query_io_mode should be true in order to use get_io_mode"));
    return io_mode;
  }



  std::string
  DataSet::get_local_no_collective_cause()
  {
    Assert(
      query_io_mode,
      ExcMessage(
        "query_io_mode should be true in order to use get_local_no_collective_cause()"));
    return internal::no_collective_cause_to_string(local_no_collective_cause);
  }



  std::uint32_t
  DataSet::get_local_no_collective_cause_as_hdf5_type()
  {
    Assert(
      query_io_mode,
      ExcMessage(
        "query_io_mode should be true in order to use get_local_no_collective_cause()"));
    return local_no_collective_cause;
  }



  std::string
  DataSet::get_global_no_collective_cause()
  {
    Assert(
      query_io_mode,
      ExcMessage(
        "query_io_mode should be true in order to use get_global_no_collective_cause()"));
    return internal::no_collective_cause_to_string(global_no_collective_cause);
  }



  std::uint32_t
  DataSet::get_global_no_collective_cause_as_hdf5_type()
  {
    Assert(
      query_io_mode,
      ExcMessage(
        "query_io_mode should be true in order to use get_global_no_collective_cause()"));
    return global_no_collective_cause;
  }


  bool
  DataSet::get_query_io_mode() const
  {
    return query_io_mode;
  }



  unsigned int
  DataSet::get_size() const
  {
    return size;
  }



  unsigned int
  DataSet::get_rank() const
  {
    return rank;
  }



  Group::Group(const std::string    &name,
               const Group          &parentGroup,
               const bool            mpi,
               const GroupAccessMode mode)
    : HDF5Object(name, mpi)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("Group", "H5Gclose");
      // Release the HDF5 resource
      const herr_t ret = H5Gclose(*pointer);
      AssertNothrow(ret >= 0, ExcInternalError());
      delete pointer;
    });
    switch (mode)
      {
        case (GroupAccessMode::open):
          *hdf5_reference =
            H5Gopen2(*(parentGroup.hdf5_reference), name.data(), H5P_DEFAULT);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Gopen2"));
          break;
        case (GroupAccessMode::create):
          *hdf5_reference = H5Gcreate2(*(parentGroup.hdf5_reference),
                                       name.data(),
                                       H5P_DEFAULT,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Gcreate2"));
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }
  }



  Group::Group(const std::string &name, const bool mpi)
    : HDF5Object(name, mpi)
  {}



  Group
  Group::open_group(const std::string &name) const
  {
    return {name, *this, mpi, GroupAccessMode::open};
  }



  Group
  Group::create_group(const std::string &name) const
  {
    return {name, *this, mpi, GroupAccessMode::create};
  }



  DataSet
  Group::open_dataset(const std::string &name) const
  {
    return {name, *hdf5_reference, mpi};
  }



  File::File(const std::string &name, const FileAccessMode mode)
    : File(name,
           mode,
           false,
#  ifdef DEAL_II_WITH_MPI
           MPI_COMM_NULL
#  else
           0
#  endif
      )
  {}



  File::File(const std::string   &name,
             const FileAccessMode mode,
             const MPI_Comm       mpi_communicator)
    : File(name, mode, true, mpi_communicator)
  {}



  File::File(const std::string   &name,
             const FileAccessMode mode,
             const bool           mpi,
             const MPI_Comm       mpi_communicator)
    : Group(name, mpi)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
      internal::HDF5ObjectImplementation::check_exception("File", "H5Fclose");
      // Release the HDF5 resource
      const herr_t err = H5Fclose(*pointer);
      AssertNothrow(err >= 0, ExcInternalError());
      delete pointer;
    });

    hid_t  plist;
    herr_t ret;

    if (mpi)
      {
#  ifdef DEAL_II_WITH_MPI
#    ifdef H5_HAVE_PARALLEL
        const MPI_Info info = MPI_INFO_NULL;

        plist = H5Pcreate(H5P_FILE_ACCESS);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_fapl_mpio(plist, mpi_communicator, info);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_fapl_mpio"));
#    else
        AssertThrow(false, ExcMessage("HDF5 parallel support is disabled."));
#    endif // H5_HAVE_PARALLEL
#  else
        AssertThrow(false, ExcMessage("MPI support is disabled."));
#  endif // DEAL_II_WITH_MPI
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    switch (mode)
      {
        case (FileAccessMode::open):
          *hdf5_reference = H5Fopen(name.data(), H5F_ACC_RDWR, plist);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Fopen"));
          break;
        case (FileAccessMode::create):
          *hdf5_reference =
            H5Fcreate(name.data(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Fcreate"));
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
          break;
      }

    if (mpi)
      {
        // Release the HDF5 resource
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }

    (void)ret;
    (void)mpi_communicator;
  }
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_HDF5
