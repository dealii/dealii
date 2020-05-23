// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_hdf5_h
#define dealii_hdf5_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/lac/full_matrix.h>

#  include <hdf5.h>

#  include <vector>

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
 * group.set_attribute("simulation_type", std::string("elastic_equation"));
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
 * File::File(const std::string &, const FileAccessMode, const MPI_Comm)
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
 * #ifdef DEBUG
 * dataset.set_query_io_mode(true);
 * #endif
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
 * # Complex numbers and HDF5
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
 *
 * @author Daniel Garcia-Sanchez, 2018, 2019.
 */
// clang-format on
namespace HDF5
{
  /**
   * Base class for the HDF5 objects.
   *
   * @author Daniel Garcia-Sanchez, 2018, 2019.
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
     * Name of the HDF5Oject. In the case of File, @p name corresponds to the
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
   *
   * @author Daniel Garcia-Sanchez, 2018, 2019.
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
    DataSet(const std::string &           name,
            const hid_t &                 parent_group_id,
            const std::vector<hsize_t> &  dimensions,
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
    write_selection(const Container &           data,
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
    write_hyperslab(const Container &           data,
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
    write_hyperslab(const Container &           data,
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
     * The return type is `uint32_t` and corresponds to the value returned by
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
    uint32_t
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
     * The return type is `uint32_t` and corresponds to the value returned by
     * H5Pget_mpio_no_collective_cause.
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
    uint32_t
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
    uint32_t local_no_collective_cause;

    /**
     * Global causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     */
    uint32_t global_no_collective_cause;
  };

  /**
   * This class implements an HDF5 Group
   *
   * @author Daniel Garcia-Sanchez, 2018,2019.
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
    Group(const std::string &   name,
          const Group &         parent_group,
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
    create_dataset(const std::string &         name,
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
   *
   * @author Daniel Garcia-Sanchez, 2018, 2019.
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
    File(const std::string &  name,
         const FileAccessMode mode,
         const MPI_Comm       mpi_communicator);

  private:
    /**
     * Delegation internal constructor.
     * File(const std::string &, const MPI_Comm, const Mode);
     * and
     * File(const std::string &, const Mode)
     * should be used to open or create HDF5 files.
     */
    File(const std::string &  name,
         const FileAccessMode mode,
         const bool           mpi,
         const MPI_Comm       mpi_communicator);
  };
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_HDF5

#endif // dealii_hdf5_h
