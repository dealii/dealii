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

#ifndef dealii_hdf5_h
#define dealii_hdf5_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/base/array_view.h>

#  include <deal.II/lac/full_matrix.h>

#  include <hdf5.h>

#  include <numeric>

DEAL_II_NAMESPACE_OPEN

// It is necessary to turn clang-format off in order to maintain the Doxygen
// links because they are longer than 80 characters
// clang-format off
/**
 * Namespace containing deal.II's HDF5 interface.
 *
 * The [Hierarchical Data Format (HDF)](https://www.hdfgroup.org/) is a cross
 * platform and a high I/O performance format designed to store large amounts of
 * data. It supports serial and MPI I/O access. This set of classes provides an
 * interface to the [HDF5 library](https://www.hdfgroup.org/downloads/hdf5/).
 *
 * The tutorial step-62 shows how to use deal.II's HDF5 interface.
 *
 * # Groups, Datasets and attributes
 * An HDF5 file is organized in
 * [groups](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Groups/HDF5_Groups.htm)
 * and
 * [datasets](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Datasets/HDF5_Datasets.htm).
 * Groups can contain datasets and other groups. Datasets are objects composed by
 * a collection of data elements. Datasets are equivalent to tensors and matrices.
 * In addition, attributes can be attached to the root file, a group or a
 * dataset. An [HDF5
 * attribute](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Attributes/HDF5_Attributes.htm)
 * is a small meta data. The methods HDF5Object::get_attribute() and
 * HDF5Object::set_attribute() can be used to get and set attributes.
 *
 * An example is shown below
 * @code
 * HDF5::File data_file(filename, HDF5::File::FileAccessMode::create);
 * double double_attribute = 2.2;
 * data_file.set_attribute("double_attribute", double_attribute);
 * auto group = data_file.create_group("group");
 * group.set_attribute("simulation_type", "elastic_equation");
 * auto dataset = group.create_dataset<double>("dataset_name", dimensions);
 * dataset.set_attribute("complex_double_attribute",
 *                       std::complex<double>(2,2.3));
 * @endcode
 *
 * # MPI I/O
 * An HDF5 file can be opened/created with serial (one single process) or
 * MPI support (several processes access the same HDF5 file).
 * File::File(const std::string &, const FileAccessMode)
 * opens/creates an HDF5 file for serial operations.
 * File::File(const std::string &, const FileAccessMode, const MPI_Comm )
 * creates or opens an HDF5 file in parallel using MPI. The HDF5 calls that
 * modify the structure of the file are always collective, whereas writing
 * and reading raw data in a dataset can be done independently or collectively.
 * [Collective access is usually faster](https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/)
 * since it allows MPI to do optimizations. In the deal.II's HDF5 interface all
 * the calls are set to collective in order to maximize the performance. This
 * means that all the MPI processes have to contribute to every single call, even
 * if they don't have data to write. MPI HDF5 requires that deal.II and HDF5 have
 * been compiled with MPI support.
 *
 * ## Write a hyperslab in parallel
 * Hyperslabs are portions of datasets. A hyperslab can be a contiguous
 * collection of points in a dataset, or it can be a regular pattern of points
 * or blocks in a datataset. Hyperslabs are equivalent to python numpy and h5py
 * [slices](http://docs.h5py.org/en/latest/high/dataset.html#reading-writing-data).
 * See the <a
 * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
 * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
 * <a
 * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
 * definition</a>.
 *
 * The example below shows how to write a simple rectangular hyperslab. The
 * offset defines the origin of the hyperslab in the original dataset. The
 * dimensions of the hyperslab are `hyperslab_dimensions = {2, 5}`. Note that
 * each process can write a hyperslab with a different size. If a process does
 * not write any data at all, the process should call the function
 * DataSet::write_none() because the operation is *collective* and all the MPI
 * processes have to contribute to the call, even if they don't have data to
 * write.
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     // hyperslab_data can be std::vector, FullMatrix or Vector
 *     FullMatrix<double> hyperslab_data = {...};
 *     std::vector<hsize_t> hyperslab_offset     = {1, 2};
 *     std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *     dataset.write_hyperslab(hyperslab_data,
 *                             hyperslab_offset,
 *                             hyperslab_dimensions);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 * @endcode
 *
 * The function
 * DataSet::write_hyperslab(const Container &,const std::vector<hsize_t> &, const std::vector<hsize_t> &)
 * is used to write simple hyperslabs and the function
 * DataSet::write_hyperslab(const Container &,const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &)
 * is used to write complex hyperslabs.
 *
 * ## Write unordered data in parallel
 * The example below shows how to write a selection of data. Note that each
 * process can write a different amount of data. If a process does not write
 * any data at all, the process should call the function
 * DataSet::write_none() because the operation is *collective* and all the MPI
 * processes have to contribute to the call, even if they don't have data to
 * write. A more detailed example can be found in step-62.
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 *
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     std::vector<hsize_t> coordinates = {0,
 *                                         0, // first point
 *                                         0,
 *                                         2, // second point
 *                                         3,
 *                                         4, // third point
 *                                         25,
 *                                         12}; // fourth point
 *     std::vector<double>  data        = {2, 3, 5, 6};
 *     dataset.write_selection(data, coordinates);
 *   }
 * else if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1)
 *   {
 *     std::vector<hsize_t> coordinates = {5,
 *                                         0, // first point
 *                                         0,
 *                                         4, // second point
 *                                         5,
 *                                         4, // third point
 *                                         26,
 *                                         12}; // fourth point
 *     std::vector<double>  data        = {9, 4, 7, 6};
 *     dataset.write_selection(data, coordinates);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 * @endcode
 *
 * ## Query the I/O mode that HDF5 used in the last parallel I/O call
 * The default access mode in the deal.II's HDF5 C++ interface  is collective
 * which is typically faster since it allows MPI to do more optimizations. In
 * some cases, such as when there is type conversion, the HDF5 library can
 * decide to do independent I/O instead of collective I/O, even if the user asks
 * for collective I/O. See the following
 * [article](https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/).
 * In cases where maximum performance is a requirement, it is important to
 * make sure that all MPI read/write operations are collective. The HDF5 library
 * provides API routines that can be used after the read/write I/O operations to
 * query the I/O mode. If DataSet::query_io_mode is True, then after
 * every read/write operation the deal.II's HDF5 interface calls the routines
 * [H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)
 * and
 * [H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
 * The results are stored in DataSet::io_mode, DataSet::local_no_collective_cause
 * and DataSet::get_global_no_collective_cause. We suggest to query the I/O mode
 * only in Debug mode because it requires calling additional HDF5 routines.
 *
 * The following code can be used to query the I/O method.
 * @code
 * auto dataset = group.create_dataset<double>("name", dimensions);
 * if constexpr (running_in_debug_mode())
 *   dataset.set_query_io_mode(true);
 *
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     dataset.write(data);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 *
 * if(dataset.get_query_io_mode()){
 *   pcout << "IO mode: " << dataset.io_mode() << std::endl;
 *   pcout << "Local no collective cause: "
 *         << dataset.local_no_collective_cause() << std::endl;
 *   pcout << "Global no collective cause: "
 *         << dataset.get_global_no_collective_cause() <<
 * std::endl;
 * }
 * @endcode
 *
 * If the write operation was collective then the output should be
 * @code
 * IO mode: H5D_MPIO_CONTIGUOUS_COLLECTIVE
 * Local no collective cause: H5D_MPIO_COLLECTIVE
 * Global no collective cause: H5D_MPIO_COLLECTIVE
 * @endcode
 * See DataSet::get_io_mode(), DataSet::get_local_no_collective_cause() and
 * DataSet::get_global_no_collective_cause() for all the possible return
 * codes.
 *
 * # Rank of HDF5 datasets and hyperslabs
 * The deal.II's HDF5 interface can be used to write/read data to datasets and
 * hyperslabs of any particular rank. `FullMatrix` can only be used to
 * write/read data to datasets and hyperslabs of rank 2. In the other hand,
 * `std::vector` and `Vector` can be used to write/read data to datasets and
 * hyperslabs of rank 1, 2, 3 and higher, the data is organized in
 * [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
 * which is commonly used in C and C++ matrices. We can re-write the code from
 * the previous section using std::vector
 * @code
 * // Dataset of rank 2. dim_0 = 50, dim_1 = 30
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     // hyperslab_data can be std::vector, FullMatrix or Vector
 *     std::vector<double> hyperslab_data = {0,1,2,3,4,5};
 *     // hyperslab of rank 2. dim_0 = 2 and dim_1 = 3
 *     std::vector<hsize_t> hyperslab_offset     = {1, 2};
 *     std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *     dataset.write_hyperslab(hyperslab_data,
 *                             hyperslab_offset,
 *                             hyperslab_dimensions);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 * @endcode
 * The previous code writes the following hyperslab matrix
 * @code
 * 0 1
 * 2 3
 * 4 5
 * @endcode
 *
 * # Datatypes
 * Attribute datatypes can be float, `double`, `std::complex<float>`,
 * `std::complex<double>`, `int`, `unsigned int`, `bool` and `std::string`.
 * HDF5Object::get_attribute() and HDF5Object::set_attribute() can be used with
 * all of these datatypes.
 *
 * Dataset datatypes can be `float`, `double`, `std::complex<float>`,
 * `std::complex<double>`, `int` and `unsigned int`. DataSet::read(),
 * DataSet::write(), DataSet::read_selection(), etc. can be used with all of
 * these datatypes. Note that the dataset datatype can not be `bool`, the
 * reason is that it can not be assumed that `std::vector<bool>` stores the
 * elements in a contiguous way.
 *
 *
 * ## Complex numbers and HDF5
 * There is no official HDF5 format to store `std::complex` numbers in a HDF5
 * file. But the *de facto* standard is to store the `std::complex` number in a
 * compound type in which `r` corresponds to the real part and `i` corresponds
 * to the imaginary part. In this interface we define two compound types one for
 * `std::complex<double>` which corresponds to `(double,double)` and another
 * one for `std::complex<float>` which corresponds to `(float,float)`. These two
 * types correspond respectively to the types of python/numpy/h5py:
 * `complex128` and `complex64`. This means that the files generated by this
 * interface will be read correctly by python/numpy/h5py and at the same time
 * this interface is able to read the files generated by python/numpy/h5py.
 *
 * # Data exchange with python scripts
 * The HDF5 format can be used to exchange data with python scripts. The strings
 * are stored as HDF5 variable-length UTF-8 strings and the complex numbers, as
 * explained above, are stored as HDF5 compound datatypes compatible with
 * [h5py](https://www.h5py.org/) and [numpy](http://www.numpy.org/).
 *
 * The following python script writes the parameters for a deal.II simulation:
 * ~~~~~~~~~~~~~{.py}
 * h5_file = h5py.File('simulation.hdf5','w')
 * data = h5_file.create_group('data')
 * data.attrs['nb_frequency_points'] = 50 # int
 * data.attrs['rho'] = 2300.5 # double
 * data.attrs['save_vtk_files'] = True # bool
 * data.attrs['simulation_type'] = 'elastic_equation' # utf8 string
 * ~~~~~~~~~~~~~
 *
 * C++ deal.II simulation with MPI HDF5:
 * @code
 * HDF5::File data_file("simulation.hdf5",
 *                      HDF5::File::FileAccessMode::open,
 *                      MPI_COMM_WORLD);
 * HDF5::Group data = data_file.open_group("data");
 *
 * auto nb_frequency_points = data.get_attribute<int>("nb_frequency_points");
 * auto rho = data.get_attribute<double>("rho");
 * auto save_vtk_files = data.get_attribute<bool>("save_vtk_files");
 * auto simulation_type = data.get_attribute<std::string>("simulation_type");
 *
 * std::vector<std::complex<double>> displacement = {...};
 *
 * data.write_dataset("displacement", displacement);
 *
 * // Write the simulation metadata
 * data.set_attribute("active_cells", triangulation.n_active_cells());
 * @endcode
 *
 * Read the simulation results with python:
 * ~~~~~~~~~~~~~{.py}
 * h5_file = h5py.File('simulation.hdf5','r+')
 * data = h5_file['data']
 * displacement = data['displacement'] # complex128 dtype
 * active_cells = data.attrs['degrees_of_freedom'])
 * ~~~~~~~~~~~~~
 *
 * # HDF5 and thread safety
 * By default HDF5 is not thread-safe. The HDF5 library can be configured to be
 * thread-safe, see [the HDF5
 * documentation](https://support.hdfgroup.org/HDF5/faq/threadsafe.html). The
 * thread-safe HDF5 version serializes the API but does not provide any level of
 * concurrency. To achieve high parallel performance with HDF5, we advice to use
 * HDF5 with MPI.
 */
// clang-format on
namespace HDF5
{
  /**
   * Base class for the HDF5 objects.
   */
  class HDF5Object
  {
  protected:
    /**
     * Constructor. @p string_name is the name of the HDF5 Object. If @p mpi is
     * True then MPI I/O is used.
     */
    HDF5Object(const std::string &name, const bool mpi);

  public:
    /**
     * Reads an attribute. @p T can be `float`, `double`, `std::complex<float>`,
     * `std::complex<double>`, `int`, `unsigned int`, `bool` or `std::string`.
     * Note that the encoding of `std::string` is UTF8 in order to be compatible
     * with python3.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    T
    get_attribute(const std::string &attr_name) const;

    /**
     * Writes an attribute. @p T can be `float`, `double`, `std::complex<float>`,
     * `std::complex<double>`, `int`, `unsigned int`, `bool` or `std::string`.
     * Note that the encoding of `std::string` is UTF8 in order to be compatible
     * with python3.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    set_attribute(const std::string &attr_name, const T value);

    /**
     * Returns the #name of the object. In the case of File, #name corresponds
     * to the file name. In the case of Group and DataSet, #name corresponds to
     * the name of the object in the HDF5 file.
     */
    std::string
    get_name() const;

  protected:
    /**
     * Name of the HDF5Object. In the case of File, @p name corresponds to the
     * file name. In the case of Group and DataSet @p name corresponds to the
     * name of the object in the HDF5 file.
     */
    const std::string name;

    /**
     * HDF5 identifier for the objects File, Group and DataSet. The
     * `std::shared_ptr<>` pointer allows the object to be copied. For example
     * several parts of the program can share and access the same group; when
     * all the functions that access the group are closed, the HDF5 resources of
     * the group will be automatically released.
     */
    std::shared_ptr<hid_t> hdf5_reference;

    /**
     * If true use parallel HDF5, if false use serial HDF5.
     */
    const bool mpi;
  };

  /**
   * This class implements an HDF5 DataSet.
   */
  class DataSet : public HDF5Object
  {
    friend class Group;

  protected:
    /**
     * Open dataset. This is an internal constructor. The function
     * Group::open_dataset() should be used to open a dataset.
     */
    DataSet(const std::string &name, const hid_t &parent_group_id, bool mpi);

    /**
     * Create dataset. This is an internal constructor. The function
     * Group::create_dataset() should be used to create a dataset.
     */
    DataSet(const std::string            &name,
            const hid_t                  &parent_group_id,
            const std::vector<hsize_t>   &dimensions,
            const std::shared_ptr<hid_t> &t_type,
            const bool                    mpi);

  public:
    /**
     * Reads all the data of the dataset.
     *
     * Datatype conversion takes place at the time of the read operation and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    template <typename Container>
    Container
    read();

    /**
     * Reads data of a subset of the dataset.
     *
     * Datatype conversion takes place at the time of the read operation and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * The selected elements can be scattered and take any shape in the dataset.
     * For example, in the case of a dataset with rank 4 a selection of 3 points
     * will be described by a 3-by-4 array. Note the indexing is zero-based. To
     * select the points (1,1,1,1), (14,6,12,18), and (8,22,30,22), the point
     * selection array would be as follows:
     *
     * @code
     *    0  0  0  0
     *
     *   13  5 11 17
     *
     *    7 21 29 21
     * @endcode
     *
     * <a
     * href="https://support.hdfgroup.org/newsletters/newsletter140.html">Parallel
     * HDF5 supports collective I/O on point selections.</a>
     *
     * Datatype conversion takes place at the time of the read operation and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename Container>
    Container
    read_selection(const std::vector<hsize_t> &coordinates);

    // clang-format off
    /**
     * Reads a hyperslab from the dataset. The parameters are summarized
     * below:
     *  - @p offset: The starting location for the hyperslab.
     *  - @p count: The number of elements to select along each dimension.
     *
     * When reading a hyperslab, HDF5 also allows to provide "stride" and
     * "block" parameters (see the [HDF5 documentation](https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab)).
     * These are not used by the current function and set to `nullptr`. However
     * these parameters can be used with the function
     * read_hyperslab(const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &)
     *
     * See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
     * <a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    // clang-format on
    template <typename Container>
    Container
    read_hyperslab(const std::vector<hsize_t> &offset,
                   const std::vector<hsize_t> &count);

    /**
     * Writes a data hyperslab to the dataset. The parameters are summarized
     * below:
     *  - @p dataset_dimensions: the dimensions of the data memory block.
     *  - @p offset: The starting location for the hyperslab.
     *  - @p stride: The number of elements to separate each element or block to
     *               be selected.
     *  - @p count: The number of elements or blocks to select along each
     *              dimension.
     *  - @p block: The size of the block selected from the dataspace.
     *
     * See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
     * <a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    template <typename Container>
    Container
    read_hyperslab(const std::vector<hsize_t> &data_dimensions,
                   const std::vector<hsize_t> &offset,
                   const std::vector<hsize_t> &stride,
                   const std::vector<hsize_t> &count,
                   const std::vector<hsize_t> &block);

    /**
     * This function does not read any data, but it can contribute to a
     * collective read call. @p number can be `float`, `double`,
     * `std::complex<float>`, `std::complex<double>`, `int` or `unsigned int`.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename number>
    void
    read_none();

    /**
     * Writes data in the dataset. @p number can be `float`, `double`,
     * `std::complex<float>`, `std::complex<double>`, `int` or `unsigned int`.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    template <typename Container>
    void
    write(const Container &data);

    /**
     * Writes data to a subset of the dataset. @p number can be `float`, `double`,
     * `std::complex<float>`, `std::complex<double>`, `int` or `unsigned int`.
     *
     * The selected elements can be scattered and take any shape in the dataset.
     * For example, in the case of a dataset with rank 4 a selection of 3 points
     * will be described by a 3-by-4 array. Note the indexing is zero-based. To
     * select the points (1,1,1,1), (14,6,12,18), and (8,22,30,22), the point
     * selection array would be as follows:
     *
     * @code
     *    0  0  0  0
     *
     *   13  5 11 17
     *
     *    7 21 29 21
     * @endcode
     *
     * <a
     * href="https://support.hdfgroup.org/newsletters/newsletter140.html">Parallel
     * HDF5 supports collective I/O on point selections.</a>
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename Container>
    void
    write_selection(const Container            &data,
                    const std::vector<hsize_t> &coordinates);

    // clang-format off
    /**
     * Writes a data hyperslab to the dataset. The parameters are summarized
     * below:
     *  - @p offset: The starting location for the hyperslab.
     *  - @p count: The number of elements to select along each dimension.
     *
     * When writing a hyperslab, HDF5 also allows to provide "stride" and
     * "block" parameters (see the [HDF5 documentation](https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab)).
     * These are not used by the current function and set to `nullptr`. However
     * these parameters can be used with the function
     * write_hyperslab(const Container &data, const std::vector<hsize_t> &data_dimensions, const std::vector<hsize_t> &offset, const std::vector<hsize_t> &stride, const std::vector<hsize_t> &count, const std::vector<hsize_t> &block).
     *
     * See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
     * <a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    // clang-format on
    template <typename Container>
    void
    write_hyperslab(const Container            &data,
                    const std::vector<hsize_t> &offset,
                    const std::vector<hsize_t> &count);

    /**
     * Writes a data hyperslab to the dataset. The parameters are summarized
     * below:
     *  - @p dataset_dimensions: the dimensions of the data memory block.
     *  - @p offset: The starting location for the hyperslab.
     *  - @p stride: The number of elements to separate each element or block to be
     *               selected.
     *  - @p count: The number of elements or blocks to select along each
     *              dimension.
     *  - @p block: The size of the block selected from the dataspace.
     *
     * See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
     * <a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    template <typename Container>
    void
    write_hyperslab(const Container            &data,
                    const std::vector<hsize_t> &data_dimensions,
                    const std::vector<hsize_t> &offset,
                    const std::vector<hsize_t> &stride,
                    const std::vector<hsize_t> &count,
                    const std::vector<hsize_t> &block);

    /**
     * This function does not write any data, but it can contribute to a
     * collective write call. In the context of a collective MPI write call,
     * if a process does not write any data at all, the process should call
     * this function because the operation is *collective* and all the MPI
     * processes have to contribute to the call, even if they don't have data
     * to write. @p number can be `float`, `double`, `std::complex<float>`,
     * `std::complex<double>`, `int` or `unsigned int`.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * An example of how to use this function can be found in step-62.
     */
    template <typename number>
    void
    write_none();

    /**
     * This function returns the boolean query_io_mode.
     *
     * In cases where maximum performance has to be achieved, it is important to
     * make sure that all MPI read/write operations are collective. The HDF5
     * library provides API routines that can be used after the read/write I/O
     * operations to query the I/O mode. If query_io_mode is set to true, then
     * after every read/write operation the deal.II's HDF5 interface calls the
     * routines
     * [H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)
     * and
     * [H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
     * The results are stored in io_mode, local_no_collective_cause and
     * global_no_collective_cause. We suggest to query the I/O mode only in
     * Debug mode because it requires calling additional HDF5 routines.
     */
    bool
    get_query_io_mode() const;

    /**
     * This function sets the boolean query_io_mode.
     */
    void
    set_query_io_mode(const bool new_query_io_mode);

    /**
     * This function returns the I/O mode that was used on the last
     * parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">H5Pget_mpio_actual_io_mode</a>.
     *
     * The return value is a `std::string` and can be
     * Value                          | Meaning
     * ------------------------------ | -------
     * H5D_MPIO_NO_COLLECTIVE         | No collective I/O was performed. Collective I/O was not requested or collective I/O isn't possible on this dataset.
     * H5D_MPIO_CHUNK_INDEPENDENT     | HDF5 performed chunk collective optimization schemes and each chunk was accessed independently.
     * H5D_MPIO_CHUNK_COLLECTIVE      | HDF5 performed chunk collective optimization and each chunk was accessed collectively.
     * H5D_MPIO_CHUNK_MIXED           | HDF5 performed chunk collective optimization and some chunks were accessed independently, some collectively.
     * H5D_MPIO_CONTIGUOUS_COLLECTIVE | Collective I/O was performed on a contiguous dataset.
     */
    std::string
    get_io_mode();

    /**
     * This function returns the I/O mode that was used on the last
     * parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">H5Pget_mpio_actual_io_mode</a>.
     * The return type is `H5D_mpio_actual_io_mode_t` which corresponds to the
     * value returned by H5Pget_mpio_actual_io_mode.
     *
     * The return value can be
     * Value                          | Meaning
     * ------------------------------ | -------
     * H5D_MPIO_NO_COLLECTIVE         | No collective I/O was performed. Collective I/O was not requested or collective I/O isn't possible on this dataset.
     * H5D_MPIO_CHUNK_INDEPENDENT     | HDF5 performed chunk collective optimization and each chunk was accessed independently.
     * H5D_MPIO_CHUNK_COLLECTIVE      | HDF5 performed chunk collective optimization and each chunk was accessed collectively.
     * H5D_MPIO_CHUNK_MIXED           | HDF5 performed chunk collective optimization and some chunks were accessed independently, some collectively.
     * H5D_MPIO_CONTIGUOUS_COLLECTIVE | Collective I/O was performed on a contiguous dataset.
     */
    H5D_mpio_actual_io_mode_t
    get_io_mode_as_hdf5_type();

    /**
     * This function returns the local causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     *
     * The return value is a string and can be
     * Value                                      | Meaning
     * ------------------------------------------ | -------
     * H5D_MPIO_COLLECTIVE                        | Collective I/O was performed successfully.
     * H5D_MPIO_SET_INDEPENDENT                   | Collective I/O was not performed because independent I/O was requested.
     * H5D_MPIO_DATATYPE_CONVERSION               | Collective I/O was not performed because datatype conversions were required.
     * H5D_MPIO_DATA_TRANSFORMS                   | Collective I/O was not performed because data transforms needed to be applied.
     * H5D_MPIO_SET_MPIPOSIX                      | Collective I/O was not performed because the selected file driver was MPI-POSIX.
     * H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES   | Collective I/O was not performed because one of the dataspaces was neither simple nor scalar.
     * H5D_MPIO_POINT_SELECTIONS                  | Collective I/O was not performed because there were point selections in one of the dataspaces.
     * H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | Collective I/O was not performed because the dataset was neither contiguous nor chunked.
     * H5D_MPIO_FILTERS                           | Collective I/O was not performed because filters needed to be applied.
     */
    std::string
    get_local_no_collective_cause();

    /**
     * This function returns the local causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     * The return type is `std::uint32_t` and corresponds to the value returned
     * by
     * [H5Pget_mpio_no_collective_cause](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
     *
     * The return value can be
     * Value                                      | Meaning
     * ------------------------------------------ | -------
     * H5D_MPIO_COLLECTIVE                        | Collective I/O was performed successfully.
     * H5D_MPIO_SET_INDEPENDENT                   | Collective I/O was not performed because independent I/O was requested.
     * H5D_MPIO_DATATYPE_CONVERSION               | Collective I/O was not performed because datatype conversions were required.
     * H5D_MPIO_DATA_TRANSFORMS                   | Collective I/O was not performed because data transforms needed to be applied.
     * H5D_MPIO_SET_MPIPOSIX                      | Collective I/O was not performed because the selected file driver was MPI-POSIX.
     * H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES   | Collective I/O was not performed because one of the dataspaces was neither simple nor scalar.
     * H5D_MPIO_POINT_SELECTIONS                  | Collective I/O was not performed because there were point selections in one of the dataspaces.
     * H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | Collective I/O was not performed because the dataset was neither contiguous nor chunked.
     * H5D_MPIO_FILTERS                           | Collective I/O was not performed because filters needed to be applied.
     */
    std::uint32_t
    get_local_no_collective_cause_as_hdf5_type();

    /**
     * This function retrieves the global causes that broke collective I/O on
     * the last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     *
     * The return value is a std::string and can be
     * Value                                      | Meaning
     * ------------------------------------------ | -------
     * H5D_MPIO_COLLECTIVE                        | Collective I/O was performed successfully.
     * H5D_MPIO_SET_INDEPENDENT                   | Collective I/O was not performed because independent I/O was requested.
     * H5D_MPIO_DATATYPE_CONVERSION               | Collective I/O was not performed because datatype conversions were required.
     * H5D_MPIO_DATA_TRANSFORMS                   | Collective I/O was not performed because data transforms needed to be applied.
     * H5D_MPIO_SET_MPIPOSIX                      | Collective I/O was not performed because the selected file driver was MPI-POSIX.
     * H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES   | Collective I/O was not performed because one of the dataspaces was neither simple nor scalar.
     * H5D_MPIO_POINT_SELECTIONS                  | Collective I/O was not performed because there were point selections in one of the dataspaces.
     * H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | Collective I/O was not performed because the dataset was neither contiguous nor chunked.
     * H5D_MPIO_FILTERS                           | Collective I/O was not performed because filters needed to be applied.
     */
    std::string
    get_global_no_collective_cause();

    /**
     * This function returns the global causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     * The return type is `std::uint32_t` and corresponds to the value returned
     * by H5Pget_mpio_no_collective_cause.
     *
     * The return value can be
     * Value                                      | Meaning
     * ------------------------------------------ | -------
     * H5D_MPIO_COLLECTIVE                        | Collective I/O was performed successfully.
     * H5D_MPIO_SET_INDEPENDENT                   | Collective I/O was not performed because independent I/O was requested.
     * H5D_MPIO_DATATYPE_CONVERSION               | Collective I/O was not performed because datatype conversions were required.
     * H5D_MPIO_DATA_TRANSFORMS                   | Collective I/O was not performed because data transforms needed to be applied.
     * H5D_MPIO_SET_MPIPOSIX                      | Collective I/O was not performed because the selected file driver was MPI-POSIX.
     * H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES   | Collective I/O was not performed because one of the dataspaces was neither simple nor scalar.
     * H5D_MPIO_POINT_SELECTIONS                  | Collective I/O was not performed because there were point selections in one of the dataspaces.
     * H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | Collective I/O was not performed because the dataset was neither contiguous nor chunked.
     * H5D_MPIO_FILTERS                           | Collective I/O was not performed because filters needed to be applied.
     */
    std::uint32_t
    get_global_no_collective_cause_as_hdf5_type();

    /**
     * This function returns the dimensions of the dataset. The vector
     * dimensions is a one-dimensional array of size rank specifying the size of
     * each dimension of the dataset.
     */
    std::vector<hsize_t>
    get_dimensions() const;

    /**
     * This function returns the total number of elements in the dataset.
     */
    unsigned int
    get_size() const;

    /**
     * This function returns the rank of the dataset.
     */
    unsigned int
    get_rank() const;

  private:
    /**
     * Rank of the DataSet
     */
    unsigned int rank;

    /**
     * The vector `dimensions` is a one-dimensional array of size rank
     * specifying the size of each dimension of the dataset.
     */
    std::vector<hsize_t> dimensions;

    /**
     * HDF5 dataspace identifier.
     */
    std::shared_ptr<hid_t> dataspace;

    /**
     * Total number of elements in the dataset.
     */
    unsigned int size;

    /**
     * If query_io_mode is set to true, then after every read/write operation
     * the deal.II's HDF5 interface calls the routines
     * [H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)
     * and
     * [H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
     * The results are stored in io_mode, local_no_collective_cause and
     * global_no_collective_cause.
     */
    bool query_io_mode;

    /**
     * I/O mode that was performed on the last parallel I/O call.
     */
    H5D_mpio_actual_io_mode_t io_mode;

    /**
     * Local causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     */
    std::uint32_t local_no_collective_cause;

    /**
     * Global causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     */
    std::uint32_t global_no_collective_cause;
  };

  /**
   * This class implements an HDF5 Group
   */
  class Group : public HDF5Object
  {
  protected:
    /**
     * Group access mode
     */
    enum class GroupAccessMode
    {
      /**
       * Opens an existing group
       */
      open,
      /**
       * Creates a new group
       */
      create
    };
    /**
     * This constructor creates or opens a group depending on the value of
     * @p mode. The group will be placed inside the group @p parent_group. The
     * parameter @p mpi defines if the I/O operations are serial or
     * parallel. This is an internal constructor, the functions open_group() and
     * create_group() of the current class should be used to open or create a
     * group.
     */
    Group(const std::string    &name,
          const Group          &parent_group,
          const bool            mpi,
          const GroupAccessMode mode);

    /**
     * Internal constructor used by File. The constructor sets the protected
     * const members of HDF5Group: @p name and @p mpi. It does not create or
     * open a Group.
     */
    Group(const std::string &name, const bool mpi);

  public:
    /**
     * Opens a sub-group of the current Group or File.
     */
    Group
    open_group(const std::string &name) const;

    /**
     * Creates a sub-group in the current Group or File.
     */
    Group
    create_group(const std::string &name) const;

    /**
     * Opens a dataset.
     */
    DataSet
    open_dataset(const std::string &name) const;

    /**
     * Creates a dataset. @p number can be `float`, `double`,
     * `std::complex<float>`, `std::complex<double>`, `int` or `unsigned int`.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     */
    template <typename number>
    DataSet
    create_dataset(const std::string          &name,
                   const std::vector<hsize_t> &dimensions) const;

    /**
     * Create and write data to a dataset. @p number can be `float`, `double`,
     * `std::complex<float>`, `std::complex<double>`, `int` or `unsigned int`.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>  section in the HDF5
     * User's Guide.
     *
     * `Container` can be `std::vector<float>`, `std::vector<double>`,
     * `std::vector<std::complex<float>>`, `std::vector<std::complex<double>>`,
     * `std::vector<int>`, `std::vector<unsigned int>`, `Vector<float>`,
     * `Vector<double>`, `Vector<std::complex<float>>`,
     * `Vector<std::complex<double>>`, `FullMatrix<float>`,
     * `FullMatrix<double>`, `FullMatrix<std::complex<float>>` or
     * `FullMatrix<std::complex<double>>`.
     */
    template <typename Container>
    void
    write_dataset(const std::string &name, const Container &data) const;
  };

  /**
   * This class implements an HDF5 File
   */
  class File : public Group
  {
  public:
    /**
     * File access mode
     */
    enum class FileAccessMode
    {
      /**
       * Read/write, file must exist
       */
      open,
      /**
       * Create file, truncate if exists
       */
      create
    };

    /**
     * Creates or opens an HDF5 file for serial operations. This call does not
     * require MPI support. It creates or opens an HDF5 file depending on the
     * value of @p mode.
     */
    File(const std::string &name, const FileAccessMode mode);

    /**
     * Creates or opens an HDF5 file in parallel using MPI. This requires that
     * deal.II and HDF5 were compiled with MPI support. It creates or opens a
     * HDF5 file depending on the value of @p mode. @p mpi_communicator
     * defines the processes that participate in this call; `MPI_COMM_WORLD` is
     * a common value for the MPI communicator.
     */
    File(const std::string   &name,
         const FileAccessMode mode,
         const MPI_Comm       mpi_communicator);

  private:
    /**
     * Delegation internal constructor.
     * File(const std::string &, const MPI_Comm , const Mode);
     * and
     * File(const std::string &, const Mode)
     * should be used to open or create HDF5 files.
     */
    File(const std::string   &name,
         const FileAccessMode mode,
         const bool           mpi,
         const MPI_Comm       mpi_communicator);
  };

  namespace internal
  {
    /**
     * This function returns the HDF5 datatype corresponding to the C++ type.
     * In the case of std::complex types the HDF5 handlers are automatically
     * freed using the destructor of `std::shared_ptr`. `std::shared_ptr` is
     * used instead of `std::unique_ptr` because the destructor of
     * `std::shared_ptr` doesn't have to be defined in the template argument. In
     * the other hand, the destructor of `std::unique` has to be defined in the
     * template argument. Native types such as `H5T_NATIVE_DOUBLE` do not
     * require a destructor, but compound types such as std::complex<double>
     * require a destructor to free the HDF5 resources.
     */
    template <typename number>
    std::shared_ptr<hid_t>
    get_hdf5_datatype();

    /**
     * Return the dimensions of `data`. For a std::vector this function returns
     * `std::vector<hsize_t>{vector_size}`.
     *
     * Several HDF5 functions such as H5Screate_simple() require a
     * one-dimensional array that specifies the size of each dimension of the
     * container, see:
     * https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-CreateSimple
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const std::vector<number> &data);

    /**
     * Return the dimensions of `data`. For a Vector this function returns
     * `std::vector<hsize_t>{vector_size}`.
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const Vector<number> &data);

    /**
     * Return the dimensions of `data`. For a FullMatrix the function returns
     * `std::vector<hsize_t>{rows, columns}`.
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const FullMatrix<number> &data);

    /**
     * This function returns the total size of the container. For a std::vector
     * the function returns `int(vector_size)`.
     */
    template <typename number>
    unsigned int
    get_container_size(const std::vector<number> &data);

    /**
     * This function returns the total size of the container. For a Vector the
     * function returns `int(vector_size)`.
     */
    template <typename number>
    unsigned int
    get_container_size(const Vector<number> &data);

    /**
     * This function returns the total size of the container. For a FullMatrix
     * the function returns `int(rows*columns)`.
     */
    template <typename number>
    unsigned int
    get_container_size(const FullMatrix<number> &data);

    /**
     * This function initializes and returns a container of type std::vector,
     * Vector or FullMatrix. The function does not set the values of the
     * elements of the container. The container can store data of a HDF5 dataset
     * or a HDF5 selection. The dimensions parameter holds the dimensions of the
     * HDF5 dataset or selection.
     *
     * In the case of a std::vector, the size of the vector will be the total
     * size given by dimensions. For example in the case of a dataset of rank 3,
     * the dimensions are `std::vector<hsize_t>{dim_0,dim_1,dim_2}`. The size of
     * the returned std::vector will be `dim_0*dim_1*dim_2`.
     *
     * In the case of a dealii::Vector, the size of the returned dealii::Vector
     * will be as well `dim_0*dim_1*dim_2`.
     *
     * A FullMatrix can store only data of HDF5 datasets with rank 2. The size
     * of the FullMatrix will be FullMatrix(dim_0,dim_2)
     */
    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, std::vector<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * Same as above.
     */
    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, Vector<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * Same as above.
     */
    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, FullMatrix<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * This helper function sets the property list of the read and write
     * operations of DataSet. A property list has to be created for the MPI
     * driver. For the serial driver the default H5P_DEFAULT can be used. In
     * addition H5Pset_dxpl_mpio is used to set the MPI mode to collective.
     */
    inline void
    set_plist(hid_t &plist, const bool mpi);

    /**
     * This helper function releases the property list handler of the read and
     * write operations of DataSet. For the serial version there is no need to
     * release the property list handler because H5P_DEFAULT has been used. If
     * query_io_mode is True then H5Pget_mpio_actual_io_mode and
     * H5Pget_mpio_no_collective_cause are used to check if the operation has
     * been collective.
     */
    inline void
    release_plist(hid_t                     &plist,
                  H5D_mpio_actual_io_mode_t &io_mode,
                  std::uint32_t             &local_no_collective_cause,
                  std::uint32_t             &global_no_collective_cause,
                  const bool                 mpi,
                  const bool                 query_io_mode);

    /**
     * Convert a HDF5 no_collective_cause code to a human readable string.
     */
    inline std::string
    no_collective_cause_to_string(const std::uint32_t no_collective_cause);
  } // namespace internal



  // definitions

  namespace internal
  {
    template <typename number>
    std::shared_ptr<hid_t>
    get_hdf5_datatype()
    {
      static_assert(std::is_same_v<number, float> ||
                      std::is_same_v<number, double> ||
                      std::is_same_v<number, int> ||
                      std::is_same_v<number, bool> ||
                      std::is_same_v<number, unsigned int> ||
                      std::is_same_v<number, std::complex<float>> ||
                      std::is_same_v<number, std::complex<double>>,
                    "The data type you are trying to get the HDF5 tag for "
                    "is not supported by this function.");

      if (std::is_same_v<number, float>)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_FLOAT);
        }
      else if (std::is_same_v<number, double>)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_DOUBLE);
        }
      else if (std::is_same_v<number, int>)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_INT);
        }
      else if (std::is_same_v<number, unsigned int>)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_UINT);
        }
      else if (std::is_same_v<number, bool>)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_HBOOL);
        }
      else if (std::is_same_v<number, std::complex<float>>)
        {
          std::shared_ptr<hid_t> t_type =
            std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
              // Release the HDF5 resource
              const herr_t ret = H5Tclose(*pointer);
              AssertNothrow(ret >= 0, ExcInternalError());
              delete pointer;
            });

          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          herr_t ret = H5Tinsert(*t_type, "r", 0, H5T_NATIVE_FLOAT);
          Assert(ret >= 0, ExcInternalError());
          ret = H5Tinsert(*t_type, "i", sizeof(float), H5T_NATIVE_FLOAT);
          Assert(ret >= 0, ExcInternalError());
          (void)ret;
          return t_type;
        }
      else if (std::is_same_v<number, std::complex<double>>)
        {
          std::shared_ptr<hid_t> t_type =
            std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
              // Release the HDF5 resource
              const herr_t ret = H5Tclose(*pointer);
              AssertNothrow(ret >= 0, ExcInternalError());
              delete pointer;
            });
          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          herr_t ret = H5Tinsert(*t_type, "r", 0, H5T_NATIVE_DOUBLE);
          Assert(ret >= 0, ExcInternalError());
          ret = H5Tinsert(*t_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);
          Assert(ret >= 0, ExcInternalError());
          (void)ret;
          return t_type;
        }

      // The function should not reach this point
      DEAL_II_ASSERT_UNREACHABLE();
      return {};
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const std::vector<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.size()};
      return dimensions;
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const Vector<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.size()};
      return dimensions;
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const FullMatrix<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.m(), data.n()};
      return dimensions;
    }



    template <typename number>
    unsigned int
    get_container_size(const std::vector<number> &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename number>
    unsigned int
    get_container_size(const Vector<number> &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename number>
    unsigned int
    get_container_size(const FullMatrix<number> &data)
    {
      return static_cast<unsigned int>(data.m() * data.n());
    }



    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, std::vector<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, Vector<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    std::enable_if_t<
      std::is_same_v<Container, FullMatrix<typename Container::value_type>>,
      Container>
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      // If the rank is higher than 2, then remove single-dimensional entries
      // from the shape defined by dimensions. This is equivalent to the squeeze
      // function of python/numpy. For example the following code would convert
      // the vector {1,3,1,2} to {3,2}
      std::vector<hsize_t> squeezed_dimensions;

      if (dimensions.size() > 2)
        {
          for (const auto &dimension : dimensions)
            {
              if (dimension > 1)
                squeezed_dimensions.push_back(dimension);
            }
        }
      else
        {
          squeezed_dimensions = dimensions;
        }

      AssertDimension(squeezed_dimensions.size(), 2);
      return Container(squeezed_dimensions[0], squeezed_dimensions[1]);
    }


    inline void
    set_plist(hid_t &plist, const bool mpi)
    {
      if (mpi)
        {
#  ifdef DEAL_II_WITH_MPI
          plist = H5Pcreate(H5P_DATASET_XFER);
          Assert(plist >= 0, ExcInternalError());
          const herr_t ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
          (void)ret;
          Assert(ret >= 0, ExcInternalError());
#  else
          AssertThrow(false, ExcNotImplemented());
#  endif
        }
      else
        {
          plist = H5P_DEFAULT;
        }

      (void)plist;
      (void)mpi;
    }


    inline void
    release_plist(hid_t                     &plist,
                  H5D_mpio_actual_io_mode_t &io_mode,
                  std::uint32_t             &local_no_collective_cause,
                  std::uint32_t             &global_no_collective_cause,
                  const bool                 mpi,
                  const bool                 query_io_mode)
    {
      if (mpi)
        {
#  ifdef DEAL_II_WITH_MPI
          herr_t ret;
          (void)ret;
          if (query_io_mode)
            {
              ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
              Assert(ret >= 0, ExcInternalError());
              ret =
                H5Pget_mpio_no_collective_cause(plist,
                                                &local_no_collective_cause,
                                                &global_no_collective_cause);
              Assert(ret >= 0, ExcInternalError());
            }
          ret = H5Pclose(plist);
          Assert(ret >= 0, ExcInternalError());
#  else
          AssertThrow(false, ExcNotImplemented());
#  endif
        }

      (void)plist;
      (void)io_mode;
      (void)local_no_collective_cause;
      (void)global_no_collective_cause;
      (void)mpi;
      (void)query_io_mode;
    }


    inline std::string
    no_collective_cause_to_string(const std::uint32_t no_collective_cause)
    {
      std::string message;

      auto append_to_message = [&message](const char *p) {
        if (message.size() > 0)
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


  template <typename T>
  T
  HDF5Object::get_attribute(const std::string &attr_name) const
  {
    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();
    T                            value;
    hid_t                        attr;
    herr_t                       ret;

    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Aopen"));
    (void)ret;
    ret = H5Aread(attr, *t_type, &value);
    Assert(ret >= 0, ExcMessage("Error at H5Aread"));
    (void)ret;
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));

    return value;
  }



  template <>
  inline std::string
  HDF5Object::get_attribute(const std::string &attr_name) const
  {
    // Reads a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c
    //
    // In the case of a variable length string the user does not have to reserve
    // memory for string_out. H5Aread will reserve the memory and the
    // user has to free the memory.
    //
    // Todo:
    // - Use H5Dvlen_reclaim instead of free

    char  *string_out;
    hid_t  attr;
    hid_t  type;
    herr_t ret;

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
    std::free(string_out);
    ret = H5Tclose(type);
    Assert(ret >= 0, ExcInternalError());

    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcInternalError());

    (void)ret;
    return string_value;
  }



  template <typename T>
  void
  HDF5Object::set_attribute(const std::string &attr_name, const T value)
  {
    hid_t  attr;
    hid_t  aid;
    herr_t ret;

    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();


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

    (void)ret;
  }



  template <>
  inline void
  HDF5Object::set_attribute(const std::string &attr_name,
                            const std::string  value) // NOLINT
  {
    // Writes a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c

    hid_t  attr;
    hid_t  aid;
    hid_t  t_type;
    herr_t ret;

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
     * In most of the cases H5Awrite and H5Dwrite take a pointer to the data.
     * But in the particular case of a variable length string, H5Awrite takes
     * the address of the pointer of the string.
     */
    const char *c_string_value = value.c_str();
    ret                        = H5Awrite(attr, t_type, &c_string_value);
    Assert(ret >= 0, ExcInternalError());

    ret = H5Sclose(aid);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));

    (void)ret;
  }



  template <typename Container>
  Container
  DataSet::read()
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    herr_t ret;

    Container data = internal::initialize_container<Container>(dimensions);

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  H5S_ALL,
                  H5S_ALL,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcInternalError());


    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_selection(const std::vector<hsize_t> &coordinates)
  {
    Assert(coordinates.size() % rank == 0,
           ExcMessage(
             "The dimension of coordinates has to be divisible by the rank"));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    std::vector<hsize_t> data_dimensions{
      static_cast<hsize_t>(coordinates.size() / rank)};

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_hyperslab(const std::vector<hsize_t> &offset,
                          const std::vector<hsize_t> &count)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    // In this particular overload of read_hyperslab the data_dimensions are
    // the same as count
    std::vector<hsize_t> data_dimensions = count;

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              nullptr,
                              count.data(),
                              nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_hyperslab(const std::vector<hsize_t> &data_dimensions,
                          const std::vector<hsize_t> &offset,
                          const std::vector<hsize_t> &stride,
                          const std::vector<hsize_t> &count,
                          const std::vector<hsize_t> &block)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              stride.data(),
                              count.data(),
                              block.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename number>
  void
  DataSet::read_none()
  {
    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    const std::vector<hsize_t>   data_dimensions = {0};

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("H5Sselect_none"));

    internal::set_plist(plist, mpi);

    // The pointer of data can safely be nullptr, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dread(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write(const Container &data)
  {
    AssertDimension(size, internal::get_container_size(data));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    herr_t ret;

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   H5S_ALL,
                   H5S_ALL,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_selection(const Container            &data,
                           const std::vector<hsize_t> &coordinates)
  {
    AssertDimension(coordinates.size(), data.size() * rank);
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    const std::vector<hsize_t> data_dimensions =
      internal::get_container_dimensions(data);

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;


    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container            &data,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &count)
  {
    AssertDimension(std::accumulate(count.begin(),
                                    count.end(),
                                    1,
                                    std::multiplies<unsigned int>()),
                    internal::get_container_size(data));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    // In this particular overload of write_hyperslab the data_dimensions are
    // the same as count
    const std::vector<hsize_t> &data_dimensions = count;

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              nullptr,
                              count.data(),
                              nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container            &data,
                           const std::vector<hsize_t> &data_dimensions,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &stride,
                           const std::vector<hsize_t> &count,
                           const std::vector<hsize_t> &block)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              stride.data(),
                              count.data(),
                              block.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename number>
  void
  DataSet::write_none()
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    std::vector<hsize_t>   data_dimensions = {0};

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5PSselect_none"));

    internal::set_plist(plist, mpi);

    // The pointer of data can safely be nullptr, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dwrite(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename number>
  DataSet
  Group::create_dataset(const std::string          &name,
                        const std::vector<hsize_t> &dimensions) const
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    return {name, *hdf5_reference, dimensions, t_type, mpi};
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
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE


#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_HDF5

#endif // dealii_hdf5_h
