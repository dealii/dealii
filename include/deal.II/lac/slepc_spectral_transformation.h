// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_slepc_spectral_transformation_h
#  define dealii_slepc_spectral_transformation_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_SLEPC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_solver.h>

#    include <petscksp.h>

#    include <slepceps.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#    ifndef DOXYGEN
namespace PETScWrappers
{
  // forward declarations
  class SolverBase;
} // namespace PETScWrappers
#    endif

namespace SLEPcWrappers
{
  // forward declaration
  class SolverBase;

  /**
   * Base class for spectral transformation classes using the SLEPc solvers
   * which are selected based on flags passed to the spectral transformation.
   *
   * <code>SLEPcWrappers::TransformationXXX</code>, where <code>XXX</code> is
   * your favourite transformation type, can then be implemented in
   * application codes in the following way for <code>XXX=INVERT</code> with
   * the solver object <code>eigensolver</code>:
   * @code
   * // Set a transformation, this one shifts the eigenspectrum by 3.142..
   * SLEPcWrappers::TransformationShift::AdditionalData
   *   additional_data(3.142);
   * SLEPcWrappers::TransformationShift shift(mpi_communicator,additional_data);
   * eigensolver.set_transformation(shift);
   * @endcode
   * and later calling the <code>solve()</code> function as usual:
   * @code
   * SolverControl solver_control (1000, 1e-9);
   * SolverArnoldi system (solver_control, mpi_communicator);
   * eigensolver.solve (A, B, lambda, x, size_of_spectrum);
   * @endcode
   *
   * @note These options can also be set at the command line.
   *
   * @ingroup SLEPcWrappers
   */
  class TransformationBase
  {
  protected:
    /**
     * Constructor.
     */
    explicit TransformationBase(const MPI_Comm mpi_communicator);

  public:
    /**
     * Destructor.
     */
    virtual ~TransformationBase();

    /**
     * Set a flag to indicate how the transformed matrices are being stored in
     * the spectral transformations.
     *
     * The possible values are given by the enumerator STMatMode in the SLEPc
     * library
     * https://slepc.upv.es/documentation/current/docs/manualpages/ST/STMatMode.html
     */
    void
    set_matrix_mode(const STMatMode mode);

    /**
     * Set solver to be used when solving a system of linear algebraic
     * equations inside the eigensolver.
     */
    void
    set_solver(const PETScWrappers::SolverBase &solver);

  protected:
    /**
     * SLEPc spectral transformation object.
     */
    ST st;

    // Make the solver class a friend, since it needs to set spectral
    // transformation object.
    friend class SolverBase;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Shift.
   *
   * @ingroup SLEPcWrappers
   */
  class TransformationShift : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the shift parameter to zero.
       */
      explicit AdditionalData(const double shift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;
    };


    /**
     * Constructor.
     */
    explicit TransformationShift(const MPI_Comm        mpi_communicator,
                                 const AdditionalData &data = AdditionalData());


  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Shift
   * and Invert.
   *
   * @ingroup SLEPcWrappers
   */
  class TransformationShiftInvert : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. By default, set the shift parameter to zero.
       */
      explicit AdditionalData(const double shift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;
    };


    /**
     * Constructor.
     */
    explicit TransformationShiftInvert(
      const MPI_Comm        mpi_communicator,
      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;

    // Make the solver class a friend, since it may need to set target
    // equal the provided shift value.
    friend class SolverBase;
  };

  /**
   * An implementation of the transformation interface using the SLEPc Cayley.
   *
   * @ingroup SLEPcWrappers
   */
  class TransformationCayley : public TransformationBase
  {
  public:
    /**
     * Standardized data struct to pipe additional data to the solver.
     */
    struct AdditionalData
    {
      /**
       * Constructor. Requires two shift parameters
       */
      explicit AdditionalData(const double shift_parameter     = 0,
                              const double antishift_parameter = 0);

      /**
       * Shift parameter.
       */
      const double shift_parameter;

      /**
       * Antishift parameter.
       */
      const double antishift_parameter;
    };


    /**
     * Constructor.
     */
    explicit TransformationCayley(
      const MPI_Comm        mpi_communicator,
      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * Store a copy of the flags for this particular solver.
     */
    const AdditionalData additional_data;
  };

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#  else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_SLEPC

/*--------------------   slepc_spectral_transformation.h   ------------------*/

#endif

/*--------------------   slepc_spectral_transformation.h   ------------------*/
