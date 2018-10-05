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

#ifndef dealii_hdf5_h
#define dealii_hdf5_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/lac/full_matrix.h>

#  include <hdf5.h>

#  include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Namespace containing the HDF5 interface.
 *
 * The [Hierarchical Data Format (HDF)](https://www.hdfgroup.org/) is a cross
 * platform and high I/O performance format designed to store large amounts of
 * data. It supports serial and MPI I/O access. This set of classes provide an
 * interface to the [C HDF5 library](https://www.hdfgroup.org/downloads/hdf5/).
 *
 * # Data exchange with python scripts
 * The HDF5 format can be used to exchange data with python scripts. The strings
 * are stored as HDF5 variable-length UTF-8 strings and the complex numbers are
 * stored as HDF5 compound datatypes compatible with
 * [h5py](https://www.h5py.org/) and [numpy](http://www.numpy.org/).
 *
 * This python script writes the parameters for a deal.ii simulation:
 * ~~~~~~~~~~~~~{.py}
 * h5_file = h5py.File('simulation.hdf5','w')
 * data = h5_file.create_group('data')
 * data.attrs['nb_frequency_points'] = 50 # int
 * data.attrs['rho'] = 2300.5 # double
 * data.attrs['save_vtk_files'] = True # bool
 * data.attrs['simulation_type'] = 'elastic_equation' # utf8 string
 * ~~~~~~~~~~~~~
 *
 * C++ deal.ii simulation with MPI HDF5:
 * @code
 * hdf5::File data_file("simulation.hdf5",
 *                      MPI_COMM_WORLD,
 *                      HDF5::File::Mode::open);
 * hdf5::Group data = data_file.group("data");
 *
 * auto nb_frequency_points = data.attr<int>("nb_frequency_points");
 * auto rho = data.attr<double>("rho");
 * auto save_vtk_files = data.attr<bool>("save_vtk_files");
 * auto simulation_type = data.attr<std::string>("simulation_type");
 *
 * std::vector<std::complex<double>> displacement = {...};
 *
 * auto some_data = data.write_dataset("displacement", displacement);
 *
 * // Write the simulation metadata
 * data.write_attr("active_cells", triangulation.n_active_cells());
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
 * # Groups, Datasets and attributes
 * The HDF5 file is organized in
 * [groups](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Groups/HDF5_Groups.htm)
 * and
 * [datasets](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Datasets/HDF5_Datasets.htm).
 * In the most comon case the file structure is a tree. Groups can contain
 * datasets and other groups. Datasets are objects composed of a collection of
 * data elements which can be seen as tensors or a matrices. The methods of the
 * DataSet class have been instantiated for the types: `float`, `double`,
 * `std::complex<float>`, `std::complex<double>`, `int` and `unsigned int`.
 *
 * In addition, attributes can be attached to the root file, a group or a
 * dataset. An [HDF5
 * attribute](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Attributes/HDF5_Attributes.htm)
 * is a small meta data. The methods HDF5Object::attr(const std::string) and
 * HDF5Object::write_attr(const std::string, const T) have been instantiated for
 * the types: `float`, `double`, `std::complex<float>`, `std::complex<double>`,
 * `int`, `unsigned int`, `bool`, and `std::string`.
 *
 * Below an example code can be found. Note that, if the group already exists
 * the method Group::group(std::string) should be used instead of
 * Group::create_group(std::string).
 * @code
 * HDF5::File data_file(filename, HDF5::File::Mode::create);
 * double double_attribute = 2.2;
 * data_file.write_attr("double_attribute", double attribute);
 * auto group = data_file.create_group("group");
 * group.write_attr("simulation_type", std::string("elastic_equation"));
 * auto dataset = group.create_dataset<double>("dataset_name", dimensions);
 * dataset.write_attr("complex_double_attribute", std::complex<double>(2,2.3));
 * @endcode
 *
 * # MPI I/O
 * The HDF5 calls that modify the structure of the file are
 * always collective, whereas writing and reading raw data in a dataset can be
 * done independently or collectively. [Collective access is usually
 * faster](https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/) since it
 * allows MPI to do optimizations. In these set of classes all the calls are set
 * to collective in order to maximize the performance. This means that all the
 * MPI processes have to contribute to every single call, even if they don't
 * have data to write.
 *
 * ## Write a hyperslab in parallel
 * Hyperslabs are portions of datasets. A hyperslab can be a contiguous
 * collection of points in a dataset, or it can be a regular pattern of points
 * or blocks in a datataset.
 *
 * See the <a
 * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
 * and Data Transfer</a>  section in the HDF5 User's Guide. See as well the
 * <a
 * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
 * definition</a>.
 *
 * The example below shows how to write a simple rectangular hyperslab. The
 * offset defines the origin of the hyperslab in the original dataset. The
 * dimensions of the hyperslab are `hyperslab_dimensions = {2, 5}`.
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     // data can be std::vector, FullMatrix or Vector
 *     FullMatrix<double> data = {...};
 *     std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *     std::vector<hsize_t> hyperslab_offset     = {1, 2};
 *     dataset.write_hyperslab(hyperslab_data,
 *                             hyperslab_offset,
 *                             hyperslab_dimensions);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 * @endcode
 */
// It is necessary to turn clang-format off in order to maintain the Doxygen
// link
// clang-format off
/**
 * The function
 * DataSet::write_hyperslab(const Container<number> &,const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &, const std::vector<hsize_t> &)
 * can be used to write complex hyperslabs.
 */
// clang-format on
/**
 *
 * ## Write unordered data in parallel
 * The example below shows how to write a selection of data.
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 *
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     std::vector<hsize_t> coordinates_a = {0,
 *                                           0, // first point
 *                                           0,
 *                                           2, // second point
 *                                           3,
 *                                           4, // third point
 *                                           25,
 *                                           12}; // fourth point
 *     std::vector<double>  data_a        = {2, 3, 5, 6};
 *     dataset.write_selection(data_a, coordinates_a);
 *   }
 * else if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1)
 *   {
 *     std::vector<hsize_t> coordinates_b = {5,
 *                                           0, // first point
 *                                           0,
 *                                           4, // second point
 *                                           5,
 *                                           4, // third point
 *                                           26,
 *                                           12}; // fourth point
 *     std::vector<double>  data_b        = {9, 4, 7, 6};
 *     dataset.write_selection(data_b, coordinates_b);
 *   }
 * else
 *   {
 *     dataset.write_none<double>();
 *   }
 * @endcode
 *
 * ## Query the I/O mode that HDF5 used on the last parallel I/O call
 * The default access mode in the HDF5 C++ interface of deal.ii is collective
 * which is typically faster since it allows MPI to do more optimizations. In
 * some cases, such as when there is type conversion, the HDF5 library can
 * decide to do independent I/O instead of collective I/O, even if the user asks
 * for collective I/O. See the following
 * [article](https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/)
 *
 * In cases where maximum performance is a requirement, it is important to
 * make sure that all MPI read/write operations are collective. The HDF5 library
 * provides API routines that can be used after the read/write I/O operations to
 * query the I/O mode. If DataSet::query_io_mode() is set to True, then after
 * every read/write operation the HDF5 deal.ii interface calls the routines
 * [H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)
 * and
 * [H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
 * The results are stored in DataSet::io_mode(),
 * DataSet::local_no_collective_cause() and
 * DataSet::global_no_collective_cause(). We suggest to query the I/O mode only
 * in Debug mode because it requires calling additional HDF5 routines.
 *
 * The following code can be used to query the I/O method.
 * @code
 * auto dataset = group.create_dataset<double>("name", dimensions);
 * #ifdef DEBUG
 * dataset.query_io_mode(true);
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
 * if(dataset.query_io_mode()){
 *   pcout << "IO mode: " << dataset.io_mode<std::string>() << std::endl;
 *   pcout << "Local no collective cause: "
 *         << dataset.local_no_collective_cause<std::string>() << std::endl;
 *   pcout << "Global no collective cause: "
 *         << dataset.global_no_collective_cause<std::string>() << std::endl;
 * }
 * @endcode
 *
 * If the write operation was collective then the output should be
 * @code
 * IO mode: H5D_MPIO_CONTIGUOUS_COLLECTIVE
 * Local no collective cause: H5D_MPIO_COLLECTIVE
 * Global no collective cause: H5D_MPIO_COLLECTIVE
 * @endcode
 * See DataSet::io_mode(), DataSet::local_no_collective_cause() and
 * DataSet::global_no_collective_cause() for all the possible returned codes.
 *
 * # Rank of the HDF5 datasets and hyperslabs
 * The deal.ii HDF5 interface can be used to write/read data to datasets and
 * hyperslabs of any particular rank. `FullMatrix` can only be used to
 * write/read data to datasets and hyperslabs of rank 2. `std::vector` and
 * `Vector` can be used to write/read data to datasets and hyperslabs of
 * rank 1, 2, 3 and higher, the data is organized in
 * [row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order)
 * which is commonly used in C and C++ matrices. We can re-write the code from
 * the previous section using std::vector
 * @code
 * // Dataset of rank 2. dim_0 = 50, dim_1 = 30
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 *   {
 *     // data can be std::vector, FullMatrix or Vector
 *     std::vector<double> data = {0,1,2,3,4,5};
 *     // hyperslab of rank 2. dim_0 = 2 and dim_1 = 3
 *     std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *     std::vector<hsize_t> hyperslab_offset     = {1, 2};
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
 * @author Daniel Garcia-Sanchez, 2018
 */
namespace HDF5
{
  /**
   * Base class for the HDF5 objects.
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class HDF5Object
  {
  protected:
    /**
     * Constructor. @p string_name is the name of the HDF5 Object. If @p mpi is
     * True then MPI I/O is used.
     */
    HDF5Object(const std::string name, const bool mpi);

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
    attr(const std::string &attr_name);

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
    write_attr(const std::string &attr_name, const T value);


    const std::string name;

  protected:
    /**
     * HDF5 identifier for the objects File, Group and DataSet. The
     * `std::share_ptr<>` pointer allows the object to be copied. For example
     * several parts of the program can share and access the same group; when
     * all the functions that access the group are closed, the HDF5 resources of
     * the group will be automatically released.
     */
    std::shared_ptr<hid_t> hdf5_reference;

    /**
     * If true use parallel HDF5, if false use serial HDF5.
     */
    const bool mpi;

    /**
     * Returned value of the last HDF5 C library call.
     *
     * The HDF5 calls return a non-negative value if successful; otherwise they
     * return a negative value.
     *
     * Calls such as
     * [H5DWrite](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html#Dataset-Write)
     * use this variable to store the returned value. Other calls such as
     * [H5Pcreate](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-Create)
     * don't use this variable to store the returned value.
     */
    herr_t ret;
  };

  /**
   * This class implements a HDF5 DataSet.
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class DataSet : public HDF5Object
  {
    friend class Group;

  protected:
    /**
     * Open dataset. This is an internal constructor. The function
     * Group::dataset() should be used to open a dataset.
     */
    DataSet(const std::string &name, const hid_t &parent_group_id, bool mpi);

    /**
     * Create dataset. This is an internal constructor. The function
     * Group::create_dataset() should be used to create a dataset.
     */
    DataSet(const std::string &         name,
            const hid_t &               parent_group_id,
            const std::vector<hsize_t> &dimensions,
            std::shared_ptr<hid_t>      t_type,
            const bool                  mpi);

  public:
    /**
     * Reads data of the dataset.
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
    read();

    /**
     * Reads data of a subset of the dataset.
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
    template <typename number>
    std::vector<number>
    read_selection(const std::vector<hsize_t> &coordinates);

    /**
     * Reads a hyperslab from the dataset. The parameters are summarized
     * below:
     *  - Offset: The starting location for the hyperslab.
     *  - Count: The number of elements to select along each dimension.
     *
     * Stride and block are set to NULL.
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
    read_hyperslab(const std::vector<hsize_t> &offset,
                   const std::vector<hsize_t> &count);

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
    template <typename number>
    void
    write_selection(const std::vector<number> & data,
                    const std::vector<hsize_t> &coordinates);

    /**
     * Writes a data hyperslab to the dataset. The parameters are summarized
     * below:
     *  - Offset: The starting location for the hyperslab.
     *  - Count: The number of elements to select along each dimension.
     */
    // It is necessary to turn clang-format off in order to maintain the Doxygen
    // link
    // clang-format off
    /**
     *
     * Stride and block are set to NULL. For complex hyperslabs see
     * write_hyperslab(const Container<number> &data, const std::vector<hsize_t> &data_dimensions, const std::vector<hsize_t> &offset, const std::vector<hsize_t> &stride, const std::vector<hsize_t> &count, const std::vector<hsize_t> &block).
     */
    // clang-format on
    /**
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
    template <typename Container>
    void
    write_hyperslab(const Container &           data,
                    const std::vector<hsize_t> &offset,
                    const std::vector<hsize_t> &count);

    /**
     * Writes a data hyperslab to the dataset. The parameters are summarized
     * below:
     *  - Dataset_dimensions: the dimensions of the data memory block.
     *  - Offset: The starting location for the hyperslab.
     *  - Stride: The number of elements to separate each element or block to be
     * selected.
     *  - Count: The number of elements or blocks to select along each
     * dimension.
     *  - Block: The size of the block selected from the dataspace.
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
     * collective write call. @p number can be `float`, `double`,
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
    write_none();

    /**
     * This funcion retrieves the type of I/O that was performed on the last
     * parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">H5Pget_mpio_actual_io_mode</a>.
     * The return type @p T can be `H5D_mpio_actual_io_mode_t` or `std::string`.
     * The type `H5D_mpio_actual_io_mode_t` corresponds to the value returned by
     * H5Pget_mpio_actual_io_mode and `std::string` is a human readable
     * conversion.
     *
     * The returned value value can be
     * Value                          | Meaning
     * ------------------------------ | -------
     * H5D_MPIO_NO_COLLECTIVE         | No collective I/O was performed. Collective I/O was not requested or collective I/O isn't possible on this dataset.
     * H5D_MPIO_CHUNK_INDEPENDENT     | HDF5 performed one the chunk collective optimization schemes and each chunk was accessed independently.
     * H5D_MPIO_CHUNK_COLLECTIVE      | HDF5 performed one the chunk collective optimization schemes and each chunk was accessed collectively.
     * H5D_MPIO_CHUNK_MIXED           | HDF5 performed one the chunk collective optimization schemes and some chunks were accessed independently, some collectively.
     * H5D_MPIO_CONTIGUOUS_COLLECTIVE | Collective I/O was performed on a contiguous dataset.
     */
    template <typename T>
    T
    io_mode();

    /**
     * This funcion retrieves the local causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     * The return type @p T can be `uint32_t` or `std::string`. The type
     * `uint32_t` corresponds to the value returned by
     * [H5Pget_mpio_no_collective_cause](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause)
     * and `std::string` is a human readable conversion.
     *
     * The returned value value can be
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
    template <typename T>
    T
    local_no_collective_cause();

    /**
     * This funcion retrieves the global causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>.
     * The return type @p T can be `uint32_t` or `std::string`. The type
     * `uint32_t` corresponds to the value returned by
     * H5Pget_mpio_no_collective_cause and `std::string` is a human readable
     * conversion.
     *
     * The returned value value can be
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
    template <typename T>
    T
    global_no_collective_cause();

    /**
     * This function retrieves the boolean query_io_mode().
     *
     * In cases where maximum performance has to be achieved, it is important to
     * make sure that all MPI read/write operations are collective. The HDF5
     * library provides API routines that can be used after the read/write I/O
     * operations to query the I/O mode. If query_io_mode() is set to True, then
     * after every read/write operation the HDF5 deal.ii interface calls the
     * routines
     * [H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)
     * and
     * [H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause).
     * The results are stored in io_mode(), local_no_collective_cause() and
     * global_no_collective_cause(). We suggest to query the I/O mode only in
     * Debug mode because it requires calling additional HDF5 routines.
     */
    bool
    query_io_mode() const;

    /**
     * This function sets the boolean query_io_mode().
     */
    void
    query_io_mode(bool query_io_mode);

    /**
     * This funcion returns the dimensions of the dataset.
     */
    std::vector<hsize_t>
    dimensions() const;

    /**
     * This funcion returns the total size of the dataset.
     */
    unsigned int
    size() const;

    /**
     * This funcion returns the rank of the dataset.
     */
    unsigned int
    rank() const;

  private:
    unsigned int              _rank;
    std::vector<hsize_t>      _dimensions;
    std::shared_ptr<hid_t>    dataspace;
    unsigned int              _size;
    bool                      _query_io_mode;
    H5D_mpio_actual_io_mode_t _io_mode;
    uint32_t                  _local_no_collective_cause;
    uint32_t                  _global_no_collective_cause;
  };

  /**
   * This class implements a HDF5 Group
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class Group : public HDF5Object
  {
  protected:
    /**
     * Group access mode
     */
    enum class Mode
    {
      /**
       * Creates a new group
       */
      create,
      /**
       * Opens an existing group
       */
      open
    };

    Group(const std::string &name,
          const Group &      parent_group,
          const bool         mpi,
          const Mode         mode);

    // This constructor is used by File
    Group(const std::string &name, const bool mpi);

  public:
    /**
     * Opens a group.
     */
    Group
    group(const std::string &name);

    /**
     * Creates a group.
     */
    Group
    create_group(const std::string &name);

    /**
     * Opens a dataset.
     */
    DataSet
    dataset(const std::string &name);

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
     * Creates and writes data to a dataset. @p number can be `float`, `double`,
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
   * This class implements a HDF5 Group
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class File : public Group
  {
  public:
    /**
     * File access mode
     */
    enum class Mode
    {
      /**
       * Create file, truncate if exists
       */
      create,
      /**
       * Read/write, file must exist
       */
      open
    };

  private:
    File(const std::string &name,
         const bool         mpi,
         const MPI_Comm     mpi_communicator,
         const Mode         mode);

  public:
    /**
     * Creates or opens a hdf5 file in parallel.
     */
    File(const std::string &name,
         const MPI_Comm     mpi_communicator,
         const Mode         mode = Mode::create);

    /**
     * Creates or opens a hdf5 file.
     */
    File(const std::string &name, const Mode mode = Mode::create);
  };
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_HDF5

#endif // dealii_hdf5_h
