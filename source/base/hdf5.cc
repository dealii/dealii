// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/base/hdf5.h>

#  include <deal.II/lac/full_matrix.h>
#  include <deal.II/lac/vector.h>

#  include <hdf5.h>

#  include <memory>
#  include <numeric>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace HDF5
{
  namespace internal
  {
    // This function gives the HDF5 datatype corresponding to the C++ type. In
    // the case of std::complex types the HDF5 handlers are automatically freed
    // using the destructor of std::shared_ptr.
    // std::shared_ptr is used instead of std::unique_ptr because the destructor
    // of std::shared_ptr doesn't have to be defined in the template argument.
    // In the other hand, the destructor of std::unique has to be defined in the
    // template argument. Native types such as H5T_NATIVE_DOUBLE does not
    // require a constructor, but compound types such as std::complex<double>
    // require a destructor to free the HDF5 resources.
    template <typename number>
    std::shared_ptr<hid_t>
    get_hdf5_datatype()
    {
      std::shared_ptr<hid_t> t_type;
      if (std::is_same<number, float>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_FLOAT);
        }
      else if (std::is_same<number, double>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_DOUBLE);
        }
      else if (std::is_same<number, int>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_INT);
        }
      else if (std::is_same<number, unsigned int>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_UINT);
        }
      else if (std::is_same<number, std::complex<float>>::value)
        {
          t_type  = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
            // Relase the HDF5 resource
            H5Tclose(*pointer);
            delete pointer;
          });
          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          H5Tinsert(*t_type, "r", 0, H5T_NATIVE_FLOAT);
          H5Tinsert(*t_type, "i", sizeof(float), H5T_NATIVE_FLOAT);
          return t_type;
        }
      else if (std::is_same<number, std::complex<double>>::value)
        {
          t_type  = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
            // Relase the HDF5 resource
            H5Tclose(*pointer);
            delete pointer;
          });
          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          H5Tinsert(*t_type, "r", 0, H5T_NATIVE_DOUBLE);
          H5Tinsert(*t_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);
          return t_type;
        }

      // The function should not reach this point
      Assert(false, ExcInternalError());
      return t_type;
    }



    // Several HDF5 functions such as H5Screate_simple() require a
    // one-dimensional array that specifies the size of each dimension of the
    // container, see:
    // https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-CreateSimple
    // The function get_container_dimensions returns a vector with the
    // dimensions of the container.
    // For a std::vector the function returns std::vector<hsize_t>{vector_size}
    // For a Vector the function returns std::vector<hsize_t>{vector_size}
    // For a FullMatrix the function returns std::vector<hsize_t>{rows, columns}
    //
    // Instead of using functions with std::enable_if, "constexpr if" from C++17
    // could be used
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      std::vector<hsize_t>>::type
    get_container_dimensions(const Container &data)
    {
      return std::vector<hsize_t>{data.size()};
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      std::vector<hsize_t>>::type
    get_container_dimensions(const Container &data)
    {
      return std::vector<hsize_t>{data.size()};
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      std::vector<hsize_t>>::type
    get_container_dimensions(const Container &data)
    {
      return std::vector<hsize_t>{data.m(), data.n()};
    }



    // This function returns the total size of the container.
    // For a std::vector the function returns int(vector_size)
    // For a Vector the function returns int(vector_size)
    // For a FullMatrix the function returns int(rows*columns)
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      unsigned int>::type
    get_container_size(const Container &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      unsigned int>::type
    get_container_size(const Container &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      unsigned int>::type
    get_container_size(const Container &data)
    {
      return static_cast<unsigned int>(data.m() * data.n());
    }



    // This function returns the pointer to the raw data of a container
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      void *>::type
    get_container_pointer(Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong
      return data.data();
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      void *>::type
    get_container_pointer(Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong.
      // Use the first element of FullMatrix to get the pointer to the raw data
      return &data[0];
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      void *>::type
    get_container_pointer(Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong.
      // Use the first element of FullMatrix to get the pointer to the raw data
      return &data[0][0];
    }



    // This function returns the pointer to the raw data of a container
    // The returned pointer is const, which means that it can be used only to
    // read the data
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      const void *>::type
    get_container_const_pointer(const Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong
      return data.data();
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      const void *>::type
    get_container_const_pointer(const Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong.
      // Use the first element of FullMatrix to get the pointer to the raw data
      return &*data.begin();
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      const void *>::type
    get_container_const_pointer(const Container &data)
    {
      // It is very important to pass the variable "data" by reference otherwise
      // the pointer will be wrong.
      // Use the first element of FullMatrix to get the pointer to the raw data
      return &data[0][0];
    }



    // This function initializes and returns a container of type std::vector,
    // Vector or FullMatrix. The function does not set the values of the
    // elements of the container. The container can store data of a HDF5 dataset
    // or a HDF5 selection. The dimensions parameter holds the dimensions of the
    // HDF5 dataset or selection.
    //
    // In the case of a std::vector, the size of the vector will be the total
    // size given by dimensions. For example in the case of a dataset of rank 3,
    // the dimensions are std::vector<hsize_t>{dim_0,dim_1,dim_2}. The size of
    // the returned std::vector will be dim_0*dim_1*dim_2
    //
    // In the case of a dealii::Vector, the size of the returned dealii::Vector
    // will be as well dim_0*dim_1*dim_2
    //
    // A FullMatrix can store only data of HDF5 datasets with rank 2. The size
    // of the FullMatrix will be FullMatrix(dim_0,dim_2)
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(std::vector<hsize_t> dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(std::vector<hsize_t> dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      Container>::type
    initialize_container(std::vector<hsize_t> dimensions)
    {
      AssertDimension(dimensions.size(), 2);
      return Container(dimensions[0], dimensions[1]);
    }



    // Convert a HDF5 no_collective_cause code to a human readable string
    std::string
    no_collective_cause_to_string(uint32_t no_collective_cause)
    {
      std::string message;

      auto append_to_message = [&message](const char *p) {
        if (message.length() > 0)
          message += ", ";
        message += p;
      };

      // The first is not a bitmask comparison, the rest are bitmask
      // comparisons.
      // https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause
      // See H5Ppublic.h
      // Hex codes are used because the HDF5 Group can deprecate some of the
      // enum codes. For example the enum code H5D_MPIO_FILTERS is not defined
      // in 1.10.2 because it is possible to use compressed datasets with the
      // MPI/IO driver.

      // H5D_MPIO_COLLECTIVE
      if (no_collective_cause == 0x00)
        {
          append_to_message("H5D_MPIO_COLLECTIVE");
        }
      // H5D_MPIO_SET_INDEPENDENT
      if ((no_collective_cause & 0x01) == 0x01)
        {
          append_to_message("H5D_MPIO_SET_INDEPENDENT");
        }
      // H5D_MPIO_DATATYPE_CONVERSION
      if ((no_collective_cause & 0x02) == 0x02)
        {
          append_to_message("H5D_MPIO_DATATYPE_CONVERSION");
        }
      // H5D_MPIO_DATA_TRANSFORMS
      if ((no_collective_cause & 0x04) == 0x04)
        {
          append_to_message("H5D_MPIO_DATA_TRANSFORMS");
        }
      // H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES
      if ((no_collective_cause & 0x10) == 0x10)
        {
          append_to_message("H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES");
        }
      // H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET
      if ((no_collective_cause & 0x20) == 0x20)
        {
          append_to_message("H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET");
        }
      // H5D_MPIO_FILTERS
      if ((no_collective_cause & 0x40) == 0x40)
        {
          append_to_message("H5D_MPIO_FILTERS");
        }
      return message;
    }
  } // namespace internal



  HDF5Object::HDF5Object(const std::string name, const bool mpi)
    : name(name)
    , mpi(mpi)
  {}



  template <typename T>
  T
  HDF5Object::get_attribute(const std::string &attr_name)
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();
    T                      value;
    hid_t                  attr;


    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Aopen"));
    ret = H5Aread(attr, *t_type, &value);
    Assert(ret >= 0, ExcMessage("Error at H5Aread"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));
    return value;
  }



  template <>
  bool
  HDF5Object::get_attribute(const std::string &attr_name)
  {
    // The enum field generated by h5py can be casted to int
    int   int_value;
    hid_t attr;
    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Aopen"));
    ret = H5Aread(attr, H5T_NATIVE_INT, &int_value);
    Assert(ret >= 0, ExcMessage("Error at H5Aread"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));
    // The int can be casted to a bool
    bool bool_value = (bool)int_value;
    return bool_value;
  }



  template <>
  std::string
  HDF5Object::get_attribute(const std::string &attr_name)
  {
    // Reads a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c
    //
    // In the case of a variable length string the user does not have to reserve
    // memory for string_out. The call HAread will reserve the memory and the
    // user has to free the memory.
    //
    // Todo:
    // - Use H5Dvlen_reclaim instead of free

    char *string_out;
    hid_t attr;
    hid_t type;

    /* Create a datatype to refer to. */
    type = H5Tcopy(H5T_C_S1);
    Assert(type >= 0, ExcInternalError());

    // Python strings are encoded in UTF8
    ret = H5Tset_cset(type, H5T_CSET_UTF8);
    Assert(type >= 0, ExcInternalError());

    ret = H5Tset_size(type, H5T_VARIABLE);
    Assert(ret >= 0, ExcInternalError());

    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcInternalError());

    ret = H5Aread(attr, type, &string_out);
    Assert(ret >= 0, ExcInternalError());

    std::string string_value(string_out);
    // The memory of the variable length string has to be freed.
    // H5Dvlen_reclaim could be also used
    free(string_out);
    H5Tclose(type);
    H5Aclose(attr);
    return string_value;
  }



  template <typename T>
  void
  HDF5Object::set_attribute(const std::string &attr_name, const T value)
  {
    hid_t                  attr;
    hid_t                  aid;
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();


    /*
     * Create scalar attribute.
     */
    aid = H5Screate(H5S_SCALAR);
    Assert(aid >= 0, ExcMessage("Error at H5Screate"));
    attr = H5Acreate2(*hdf5_reference,
                      attr_name.data(),
                      *t_type,
                      aid,
                      H5P_DEFAULT,
                      H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Acreate2"));

    /*
     * Write scalar attribute.
     */
    ret = H5Awrite(attr, *t_type, &value);
    Assert(ret >= 0, ExcMessage("Error at H5Awrite"));

    ret = H5Sclose(aid);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));
  }



  template <>
  void
  HDF5Object::set_attribute(const std::string &attr_name,
                            const std::string  value)
  {
    // Writes a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c
    //
    // In the case of a variable length string, H5Awrite needs the address of a
    // (char *). For this reason the std::string value has been copied to a C
    // string.

    hid_t attr;
    hid_t aid;
    hid_t t_type;

    // Reserve space for the string and the null terminator
    char *c_string_value = (char *)malloc(sizeof(char) * (value.size() + 1));
    strcpy(c_string_value, value.data());

    /* Create a datatype to refer to. */
    t_type = H5Tcopy(H5T_C_S1);
    Assert(t_type >= 0, ExcInternalError());

    // Python strings are encoded in UTF8
    ret = H5Tset_cset(t_type, H5T_CSET_UTF8);
    Assert(t_type >= 0, ExcInternalError());

    ret = H5Tset_size(t_type, H5T_VARIABLE);
    Assert(ret >= 0, ExcInternalError());

    /*
     * Create scalar attribute.
     */
    aid = H5Screate(H5S_SCALAR);
    Assert(aid >= 0, ExcMessage("Error at H5Screate"));
    attr = H5Acreate2(
      *hdf5_reference, attr_name.data(), t_type, aid, H5P_DEFAULT, H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Acreate2"));

    /*
     * Write scalar attribute.
     */
    ret = H5Awrite(attr, t_type, &c_string_value);
    Assert(ret >= 0, ExcInternalError());

    free(c_string_value);
    ret = H5Sclose(aid);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));
  }



  DataSet::DataSet(const std::string &name,
                   const hid_t &      parent_group_id,
                   const bool         mpi)
    : HDF5Object(name, mpi)
    , query_io_mode(false)
    , io_mode(H5D_MPIO_NO_COLLECTIVE)
    , local_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
    , global_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Dclose(*pointer);
      delete pointer;
    });
    dataspace      = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Sclose(*pointer);
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
    rank          = rank_ret;
    hsize_t *dims = (hsize_t *)malloc(rank * sizeof(hsize_t));
    rank_ret      = H5Sget_simple_extent_dims(*dataspace, dims, NULL);
    AssertDimension(rank_ret, static_cast<int>(rank));
    dimensions.assign(dims, dims + rank);
    free(dims);

    size = 1;
    for (auto &&dimension : dimensions)
      {
        size *= dimension;
      }
  }



  DataSet::DataSet(const std::string &         name,
                   const hid_t &               parent_group_id,
                   const std::vector<hsize_t> &dimensions,
                   std::shared_ptr<hid_t>      t_type,
                   const bool                  mpi)
    : HDF5Object(name, mpi)
    , rank(dimensions.size())
    , dimensions(dimensions)
    , query_io_mode(false)
    , io_mode(H5D_MPIO_NO_COLLECTIVE)
    , local_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
    , global_no_collective_cause(H5D_MPIO_SET_INDEPENDENT)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Dclose(*pointer);
      delete pointer;
    });
    dataspace      = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Sclose(*pointer);
      delete pointer;
    });

    *dataspace = H5Screate_simple(rank, dimensions.data(), NULL);
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
    for (auto &&dimension : dimensions)
      {
        size *= dimension;
      }
  }



  template <typename Container>
  Container
  DataSet::read()
  {
    std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t plist;

    Container data = internal::initialize_container<Container>(dimensions);

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcInternalError());
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcInternalError());
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  H5S_ALL,
                  H5S_ALL,
                  plist,
                  internal::get_container_pointer(data));
    Assert(ret >= 0, ExcInternalError());

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcInternalError());
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0, ExcInternalError());
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcInternalError());
      }

    return data;
  }



  template <typename number>
  std::vector<number>
  DataSet::read_selection(const std::vector<hsize_t> &coordinates)
  {
    Assert(coordinates.size() % rank == 0,
           ExcMessage(
             "The dimension of coordinates has to be divisible by the rank"));
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    hid_t                  plist;
    hid_t                  memory_dataspace;

    std::vector<hsize_t> data_dimensions{
      static_cast<hsize_t>(coordinates.size() / rank)};

    std::vector<number> data(
      static_cast<unsigned int>(coordinates.size() / rank));

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  data.data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    return data;
  }



  template <typename Container>
  Container
  DataSet::read_hyperslab(const std::vector<hsize_t> &offset,
                          const std::vector<hsize_t> &count)
  {
    std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t plist;
    hid_t memory_dataspace;

    // In this particular overload of read_hyperslab the data_dimensions are
    // the same as count
    std::vector<hsize_t> data_dimensions = count;

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(
      *dataspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  internal::get_container_pointer(data));
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    return data;
  }



  template <typename number>
  void
  DataSet::read_none()
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    std::vector<hsize_t>   data_dimensions = {0};

    hid_t memory_dataspace;
    hid_t plist;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("H5Sselect_none"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    // The pointer of data can safely be NULL, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dread(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, NULL);
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));
  }



  template <typename Container>
  void
  DataSet::write(const Container &data)
  {
    AssertDimension(size, internal::get_container_size(data));
    std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t plist;

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   H5S_ALL,
                   H5S_ALL,
                   plist,
                   internal::get_container_const_pointer(data));
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcInternalError());
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0, ExcInternalError());
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
  }



  template <typename number>
  void
  DataSet::write_selection(const std::vector<number> & data,
                           const std::vector<hsize_t> &coordinates)
  {
    AssertDimension(coordinates.size(), data.size() * rank);
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    std::vector<hsize_t>   data_dimensions = {data.size()};

    hid_t memory_dataspace;
    hid_t plist;


    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   data.data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container &           data,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &count)
  {
    AssertDimension(std::accumulate(count.begin(),
                                    count.end(),
                                    1,
                                    std::multiplies<unsigned int>()),
                    internal::get_container_size(data));
    std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    // In this particular overload of write_hyperslab the data_dimensions are
    // the same as count
    std::vector<hsize_t> data_dimensions = count;

    hid_t memory_dataspace;
    hid_t plist;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(
      *dataspace, H5S_SELECT_SET, offset.data(), NULL, count.data(), NULL);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }
    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   internal::get_container_const_pointer(data));
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container &           data,
                           const std::vector<hsize_t> &data_dimensions,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &stride,
                           const std::vector<hsize_t> &count,
                           const std::vector<hsize_t> &block)
  {
    std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();

    hid_t memory_dataspace;
    hid_t plist;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              stride.data(),
                              count.data(),
                              block.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_dxpl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }
    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   internal::get_container_const_pointer(data));
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcMessage("Error at H5Pget_mpio_actual_io_mode"));
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0,
                   ExcMessage("Error at H5Pget_mpio_no_collective_cause"));
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
  }



  template <typename number>
  void
  DataSet::write_none()
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    std::vector<hsize_t>   data_dimensions = {0};

    hid_t memory_dataspace;
    hid_t plist;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), NULL);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5PSselect_none"));

    if (mpi)
      {
        plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    // The pointer of data can safely be NULL, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dwrite(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, NULL);
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    if (mpi)
      {
        if (query_io_mode)
          {
            ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
            Assert(ret >= 0, ExcInternalError());
            ret = H5Pget_mpio_no_collective_cause(plist,
                                                  &local_no_collective_cause,
                                                  &global_no_collective_cause);
            Assert(ret >= 0, ExcInternalError());
          }
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
  }



  std::string
  DataSet::get_io_mode()
  {
    Assert(query_io_mode,
           ExcMessage("query_io_mode should be true in order to use io_mode"));
    switch (io_mode)
      {
        case (H5D_MPIO_NO_COLLECTIVE):
          return std::string("H5D_MPIO_NO_COLLECTIVE");
          break;
        case (H5D_MPIO_CHUNK_INDEPENDENT):
          return std::string("H5D_MPIO_CHUNK_INDEPENDENT");
          break;
        case (H5D_MPIO_CHUNK_COLLECTIVE):
          return std::string("H5D_MPIO_CHUNK_COLLECTIVE");
          break;
        case (H5D_MPIO_CHUNK_MIXED):
          return std::string("H5D_MPIO_CHUNK_MIXED");
          break;
        case (H5D_MPIO_CONTIGUOUS_COLLECTIVE):
          return std::string("H5D_MPIO_CONTIGUOUS_COLLECTIVE");
          break;
        default:
          Assert(false, ExcInternalError());
          return std::string("Internal error");
          break;
      }
    // The function should not reach this line.
    Assert(false, ExcInternalError());
    return std::string("Internal error");
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



  uint32_t
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



  uint32_t
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



  void
  DataSet::set_query_io_mode(bool new_query_io_mode)
  {
    query_io_mode = new_query_io_mode;
  }



  std::vector<hsize_t>
  DataSet::get_dimensions() const
  {
    return dimensions;
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



  Group::Group(const std::string &name,
               const Group &      parentGroup,
               const bool         mpi,
               const Mode         mode)
    : HDF5Object(name, mpi)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Gclose(*pointer);
      delete pointer;
    });
    switch (mode)
      {
        case (Mode::create):
          *hdf5_reference = H5Gcreate2(*(parentGroup.hdf5_reference),
                                       name.data(),
                                       H5P_DEFAULT,
                                       H5P_DEFAULT,
                                       H5P_DEFAULT);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Gcreate2"));
          break;
        case (Mode::open):
          *hdf5_reference =
            H5Gopen2(*(parentGroup.hdf5_reference), name.data(), H5P_DEFAULT);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Gopen2"));
          break;
        default:
          Assert(false, ExcInternalError());
          break;
      }
  }



  Group::Group(const std::string &name, const bool mpi)
    : HDF5Object(name, mpi)
  {}



  Group
  Group::group(const std::string &name)
  {
    return Group(name, *this, mpi, Mode::open);
  }



  Group
  Group::create_group(const std::string &name)
  {
    return Group(name, *this, mpi, Mode::create);
  }



  DataSet
  Group::dataset(const std::string &name)
  {
    return DataSet(name, *hdf5_reference, mpi);
  }



  template <typename number>
  DataSet
  Group::create_dataset(const std::string &         name,
                        const std::vector<hsize_t> &dimensions) const
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    return DataSet(name, *hdf5_reference, dimensions, t_type, mpi);
  }



  template <typename Container>
  void
  Group::write_dataset(const std::string &name, const Container &data) const
  {
    std::vector<hsize_t> dimensions = internal::get_container_dimensions(data);
    auto                 dataset =
      create_dataset<typename Container::value_type>(name, dimensions);
    dataset.write(data);
  }

  File::File(const std::string &name,
             const bool         mpi,
             const MPI_Comm     mpi_communicator,
             const Mode         mode)
    : Group(name, mpi)
  {
    hdf5_reference = std::shared_ptr<hid_t>(new hid_t, [](auto pointer) {
      // Relase the HDF5 resource
      H5Fclose(*pointer);
      delete pointer;
    });

    hid_t          plist;
    const MPI_Info info = MPI_INFO_NULL;

    if (mpi)
      {
#  ifndef DEAL_II_WITH_MPI
        AssertThrow(false, ExcMessage("MPI support is disabled."));
#  endif // DEAL_II_WITH_MPI
#  ifndef H5_HAVE_PARALLEL
        AssertThrow(false, ExcMessage("HDF5 parallel support is disabled."));
#  endif // H5_HAVE_PARALLEL
        plist = H5Pcreate(H5P_FILE_ACCESS);
        Assert(plist >= 0, ExcMessage("Error at H5Pcreate"));
        ret = H5Pset_fapl_mpio(plist, mpi_communicator, info);
        Assert(ret >= 0, ExcMessage("Error at H5Pset_fapl_mpio"));
      }
    else
      {
        plist = H5P_DEFAULT;
      }

    switch (mode)
      {
        case (Mode::create):
          *hdf5_reference =
            H5Fcreate(name.data(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Fcreate"));
          break;
        case (Mode::open):
          *hdf5_reference = H5Fopen(name.data(), H5F_ACC_RDWR, plist);
          Assert(*hdf5_reference >= 0, ExcMessage("Error at H5Fopen"));
          break;
        default:
          Assert(false, ExcInternalError());
          break;
      }

    if (mpi)
      {
        // Relase the HDF5 resource
        ret = H5Pclose(plist);
        Assert(ret >= 0, ExcMessage("Error at H5Pclose"));
      }
  }



  File::File(const std::string &name,
             const MPI_Comm     mpi_communicator,
             const Mode         mode)
    : File(name, true, mpi_communicator, mode)
  {}



  File::File(const std::string &name, const Mode mode)
    : File(name, false, MPI_COMM_NULL, mode)
  {}



#  ifndef DOXYGEN

  // instantiations of functions

#    include "hdf5.inst"

  // Instantiations for int and unsigned int can be found below. Only
  // std::vector is provided with instantiations with int and unsigned int.
  // Instantiations of FullMatrix and Vector are provided with float, double,
  // std::complex<float> and std::complex<double> in hdf5.inst

  template int
  HDF5Object::get_attribute<int>(const std::string &attr_name);
  template unsigned int
  HDF5Object::get_attribute<unsigned int>(const std::string &attr_name);
  // The specializations of HDF5Object::get_attribute<std::string>
  // and HDF5Object::get_attribute<bool> have been defined above

  template void
  HDF5Object::set_attribute<int>(const std::string &attr_name, int value);
  template void
  HDF5Object::set_attribute<unsigned int>(const std::string &attr_name,
                                          unsigned int       value);

  template std::vector<int>
  DataSet::read<std::vector<int>>();
  template std::vector<unsigned int>
  DataSet::read<std::vector<unsigned int>>();

  template std::vector<int>
  DataSet::read_selection<int>(const std::vector<hsize_t> &coordinates);
  template std::vector<unsigned int>
  DataSet::read_selection<unsigned int>(
    const std::vector<hsize_t> &coordinates);

  template std::vector<int>
  DataSet::read_hyperslab<std::vector<int>>(const std::vector<hsize_t> &offset,
                                            const std::vector<hsize_t> &count);
  template std::vector<unsigned int>
  DataSet::read_hyperslab<std::vector<unsigned int>>(
    const std::vector<hsize_t> &offset,
    const std::vector<hsize_t> &count);

  template void
  DataSet::read_none<int>();
  template void
  DataSet::read_none<unsigned int>();

  template void
  DataSet::write<std::vector<int>>(const std::vector<int> &data);
  template void
  DataSet::write<std::vector<unsigned int>>(
    const std::vector<unsigned int> &data);

  template void
  DataSet::write_selection<int>(const std::vector<int> &    data,
                                const std::vector<hsize_t> &coordinates);
  template void
  DataSet::write_selection<unsigned int>(
    const std::vector<unsigned int> &data,
    const std::vector<hsize_t> &     coordinates);

  template void
  DataSet::write_hyperslab<std::vector<int>>(const std::vector<int> &    data,
                                             const std::vector<hsize_t> &offset,
                                             const std::vector<hsize_t> &count);
  template void
  DataSet::write_hyperslab<std::vector<unsigned int>>(
    const std::vector<unsigned int> &data,
    const std::vector<hsize_t> &     offset,
    const std::vector<hsize_t> &     count);

  template void
  DataSet::write_hyperslab<std::vector<int>>(
    const std::vector<int> &    data,
    const std::vector<hsize_t> &data_dimensions,
    const std::vector<hsize_t> &offset,
    const std::vector<hsize_t> &stride,
    const std::vector<hsize_t> &count,
    const std::vector<hsize_t> &block);
  template void
  DataSet::write_hyperslab<std::vector<unsigned int>>(
    const std::vector<unsigned int> &data,
    const std::vector<hsize_t> &     data_dimensions,
    const std::vector<hsize_t> &     offset,
    const std::vector<hsize_t> &     stride,
    const std::vector<hsize_t> &     count,
    const std::vector<hsize_t> &     block);

  template void
  DataSet::write_none<int>();
  template void
  DataSet::write_none<unsigned int>();

  template DataSet
  Group::create_dataset<int>(const std::string &         name,
                             const std::vector<hsize_t> &dimensions) const;
  template DataSet
  Group::create_dataset<unsigned int>(
    const std::string &         name,
    const std::vector<hsize_t> &dimensions) const;

  template void
  Group::write_dataset<std::vector<int>>(const std::string &     name,
                                         const std::vector<int> &data) const;
  template void
  Group::write_dataset<std::vector<unsigned int>>(
    const std::string &              name,
    const std::vector<unsigned int> &data) const;

#  endif // DOXYGEN

} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_HDF5
