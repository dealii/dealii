// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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
 * Namespace containing the HDF5 wrappers.
 *
 * This set of classes can be used to store the data in the hdf5 file format.
 *
 * The classes can be used to improve the interaction with python scripts.
 */
// clang-format off
/**
 * This python script writes the parameters for a deal.ii simulation:
 * @code
 * h5_file = h5py.File('simulation.hdf5','w')
 * data = h5_file.create_group('data')
 * data.attrs['nb_frequency_points'] = 50 # int
 * data.attrs['nb_slice_points'] = 10 # int
 * data.attrs['rho'] = 2300.5 # double
 * data.attrs['simulation_type'] = 'elastic_equation' # utf8 string
 * @endcode
 *
 * C++ deal.ii simulation with MPI HDF5:
 * @code
 * hdf5::File data_file("simulation.hdf5", MPI_COMM_WORLD, hdf5::File::Mode::open);
 * hdf5::Group data = data_file.group("data");
 *
 * // Read some parameters
 * const int nb_frequency_points = data.attr<int>("nb_frequency_points");
 * const int nb_slice_points = data.attr<int>("nb_cut_points");
 * const double rho = data.attr<double>("rho");
 * const std::string simulation_type =
 * data.attr<std::string>("simulation_type");
 *
 * // Create some distributed datasets.
 * // Note that the keyword template is used between data and create because
 * // the function create_dataset returns a dependent template.
 * auto frequency_dataset = data.create_dataset<double>("frequency", std::vector<hsize_t>{nb_frequency_points});
 *
 * auto displacement = data.create_dataset<std::complex<double>>("displacement",
 *                                                           std::vector<hsize_t>{nb_slice_points, nb_frequency_points}));
 *
 * // ... do some calculations ...
 * // The displacement data is distributed across the MPI processes.
 * // The displacement_data vector is local. The coordinates vector corresponds
 * // to the coordinates of the data in the hdf5 dataset.
 * if (coordinates.size() > 0) {
 *   displacement.write_data_selection(displacement_data, coordinates);
 * } else {
 *   // Even if this MPI process has no displacement data to write,
 *   // it has to contribute to the collective call
 *   displacement.write_data_none();
 * }
 *
 * // ... do some calculations ...
 * // The whole frequency dataset can be written by one single MPI process,
 * // but the rest of the processes have to contribute to the collective call
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) {
 *   frequency_dataset.write_data(frequency);
 * } else {
 *   frequency_dataset.write_data_none();
 * }
 *
 * // Write the simulation metadata
 * data.write_attr("active_cells", triangulation.n_active_cells());
 * data.write_attr("degrees_of_freedom", dof_handler.n_dofs());
 * @endcode
 *
 * Read the simulation data with python:
 * @code
 * h5_file = h5py.File('simulation.hdf5','r+')
 * data = h5_file['data']
 * frequency = data['frequency'] # float64 dtype
 * displacement = data['displacement'] # complex128 dtype
 * print('Degrees of freedom: ' + repr(data.attrs['degrees_of_freedom']))
 * @endcode
 */
// clang-format on
/**
 *
 * @author Daniel Garcia-Sanchez, 2018
 */
namespace HDF5

{
  /**
   * General class for the HDF5 objects.
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class HDF5Object
  {
  protected:
    HDF5Object(const std::string name, bool mpi);

  public:
    enum class Mode
    {
      create,
      open
    };

    /**
     * Reads an attribute. T can be double, int, unsigned int, bool,
     * std::complex<double> or std::string. Note that the encoding of
     * std::string is UTF8 in order to be compatible with python3.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    T
    attr(const std::string attr_name) const;

    /**
     * Writes an attribute. T can be double, int, unsigned int, bool,
     * std::complex<double> or std::string. Note that the encoding of
     * std::string is UTF8 in order to be compatible with python3.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    write_attr(const std::string attr_name, const T value) const;


    const std::string name;

  protected:
    std::shared_ptr<hid_t> hdf5_reference;
    // If true use mpi hdf5. If false use serial hdf5.
    const bool mpi;
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
    // Open dataset
    DataSet(const std::string name, const hid_t &parent_group_id, bool mpi);

    // Create dataset
    DataSet(const std::string      name,
            const hid_t &          parent_group_id,
            std::vector<hsize_t>   dimensions,
            std::shared_ptr<hid_t> t_type,
            bool                   mpi);

  public:
    /**
     * Reads data of the dataset. T can be double, int, unsigned int, bool
     * or std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <template <class...> class Container, typename T>
    Container<T>
    read_data();

    /**
     * Writes data in the dataset. T can be double, int, unsigned int,
     * or std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <template <class...> class Container, typename T>
    void
    write_data(const Container<T> &data);

    /**
     * Writes data to a subset of the dataset. T can be double, int, unsigned
     * int, bool or std::complex<double>.
     *
     * The selected elements can be scattered and take any shape in the dataset.
     * For examples, in the 4D case, we will be selecting three points and a 4D
     * dataspace has rank 4, so the selection will be described in a 3-by-4
     * array. To select the points (1,1,1,1), (14,6,12,18), and (8,22,30,22),
     * the point selection array would be as follows:
     *
     *    0  0  0  0
     *
     *   13  5 11 17
     *
     *    7 21 29 21
     *
     * <a
     * href="https://support.hdfgroup.org/newsletters/newsletter140.html">Parallel
     * HDF5 supports collective I/O on point selections.</a>
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    write_data_selection(const std::vector<T> &     data,
                         const std::vector<hsize_t> coordinates);

    /**
     * Writes data to a subset of the dataset. T can be double, int, unsigned
     * int, bool or std::complex<double>.
     *
     * The selected elements form a hyperslab in the dataset.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <template <class...> class Container, typename T>
    void
    write_data_hyperslab(const Container<T> &       data,
                         const std::vector<hsize_t> offset,
                         const std::vector<hsize_t> count);

    /**
     * This function does not write any data, but can contribute to a collective
     * write call. T can be double, int, unsigned int, bool or
     * std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    write_data_none();

    /**
     * This funcion retrieves the type of I/O that was performed on the last
     * parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">"H5Pget_mpio_actual_io_mode"</a>.
     * The return type T can be H5D_mpio_actual_io_mode_t or std::string. The
     * type H5D_mpio_actual_io_mode_t corresponds to the value returned by
     * H5Pget_mpio_actual_io_mode and std::string is a human readable
     * conversion.
     */
    template <typename T>
    T
    io_mode();

    /**
     * This funcion retrieves the local causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">"H5Pget_mpio_no_collective_cause"</a>.
     * The return type T can be uint32_t or std::string. The type uint32_t
     * corresponds to the value returned by H5Pget_mpio_no_collective_cause and
     * std::string is a human readable conversion.
     */
    template <typename T>
    T
    local_no_collective_cause();

    /**
     * This funcion retrieves the global causes that broke collective I/O on the
     * last parallel I/O call. See <a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">"H5Pget_mpio_no_collective_cause"</a>.
     * The return type T can be uint32_t or std::string. The type uint32_t
     * corresponds to the value returned by H5Pget_mpio_no_collective_cause and
     * std::string is a human readable conversion.
     */
    template <typename T>
    T
    global_no_collective_cause();

    /**
     * This function retrieves the IO mode checking. If check_io_mode is true,
     * then after every read and write operation in the dataset, it will be
     * retrieved the type of I/O that was performed on the last parallel I/O
     * call If check_io_mode is false then no checking will be performed.
     */
    bool
    check_io_mode() const;

    /**
     * This funcion sets the IO mode checking. If check_io_mode is true, then
     * after every read and write operation in the dataset, it will be retrieved
     * the type of I/O that was performed on the last parallel I/O call If
     * check_io_mode is false then no checking will be performed.
     */
    void
    check_io_mode(bool check_io_mode);

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
    bool                      _check_io_mode;
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
    Group(const std::string name,
          const Group &     parent_group,
          const bool        mpi,
          const Mode        mode);
    Group(const std::string name, const bool mpi);

  public:
    /**
     * Opens a group.
     */
    Group
    group(const std::string name);

    /**
     * Creates a group.
     */
    Group
    create_group(const std::string name);

    /**
     * Opens a dataset.
     */
    DataSet
    dataset(const std::string name);

    /**
     * Creates a dataset. T can be double, int, unsigned int, bool or
     * std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    DataSet
    create_dataset(const std::string          name,
                   const std::vector<hsize_t> dimensions) const;

    /**
     * Creates and writes data to a dataset. T can be double, int, unsigned int,
     * bool or std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    write_dataset(const std::string name, const std::vector<T> &data) const;

    /**
     * Creates and writes data to a dataset. T can be double, int, unsigned int,
     * bool or std::complex<double>.
     *
     * Datatype conversion takes place at the time of a read or write and is
     * automatic. See the <a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">"Data
     * Transfer: Datatype Conversion and Selection"</a>  section in the HDF5
     * User's Guide.
     */
    template <typename T>
    void
    write_dataset(const std::string name, const FullMatrix<T> &data) const;
  };

  /**
   * This class implements a HDF5 Group
   *
   * @author Daniel Garcia-Sanchez, 2018
   */
  class File : public Group
  {
  private:
    File(const std::string name,
         const bool        mpi,
         const MPI_Comm    mpi_communicator,
         const Mode        mode);

  public:
    /**
     * Creates or opens a hdf5 file.
     */
    File(const std::string name,
         const MPI_Comm    mpi_communicator,
         const Mode        mode = Mode::create);

    /**
     * Creates or opens a hdf5 file in parallel.
     */
    File(const std::string name, const Mode mode = Mode::create);
  };
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_HDF5

#endif // dealii_hdf5_h
