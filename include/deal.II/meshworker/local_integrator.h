// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_mesh_worker_local_integrator_h
#define dealii_mesh_worker_local_integrator_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>

#include <functional>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  // Forward declarations
#ifndef DOXYGEN
  template <int dim, int spacedim, typename number>
  class DoFInfo;
  template <int dim, int spacedim>
  class IntegrationInfo;
#endif

  /**
   * A local integrator object, which can be used to simplify the call of
   * loop(). Instead of providing the three local integration functions
   * separately, we bundle them as virtual functions in this class.
   *
   * Additionally, since we cannot have a virtual null function, we provide
   * flags, which allow us to indicate, whether we want to integrate on
   * boundary and interior faces. These flags are true by default, but can be
   * modified by applications to speed up the loop.
   *
   * If a function is not overloaded in a derived class, but its usage flag is
   * true, the function will cause an exception ExcPureFunction.
   *
   * @deprecated This class is deprecated. It used to be the basis for
   *   integration via the MeshWorker::integration_loop() function, but the
   *   same functionality is available via MeshWorker::loop().
   *
   * @ingroup MeshWorker
   */
  template <int dim, int spacedim = dim, typename number = double>
  class DEAL_II_DEPRECATED LocalIntegrator : public EnableObserverPointer
  {
  public:
    /**
     * The constructor setting default values, namely all integration flags to
     * true.
     */
    LocalIntegrator();

    /**
     * The constructor setting integration flags to specified values.
     */
    LocalIntegrator(bool use_cell, bool use_boundary, bool use_face);

    /**
     * The empty virtual destructor.
     */
    virtual ~LocalIntegrator() override = default;

    /**
     * Virtual function for integrating on cells. Throws exception
     * PureFunctionCalled if not overloaded by a derived class.
     */
    virtual void
    cell(DoFInfo<dim, spacedim, number> &dinfo,
         IntegrationInfo<dim, spacedim> &info) const;
    /**
     * Virtual function for integrating on boundary faces. Throws exception
     * PureFunctionCalled if not overloaded by a derived class.
     */
    virtual void
    boundary(DoFInfo<dim, spacedim, number> &dinfo,
             IntegrationInfo<dim, spacedim> &info) const;
    /**
     * Virtual function for integrating on interior faces. Throws exception
     * PureFunctionCalled if not overloaded by a derived class.
     */
    virtual void
    face(DoFInfo<dim, spacedim, number> &dinfo1,
         DoFInfo<dim, spacedim, number> &dinfo2,
         IntegrationInfo<dim, spacedim> &info1,
         IntegrationInfo<dim, spacedim> &info2) const;

    /**
     * The flag indicating whether the cell integrator cell() is to be used in
     * the loop. Defaults to <tt>true</tt>.
     */
    bool use_cell;

    /**
     * The flag indicating whether the boundary integrator boundary() is to be
     * used in the loop. Defaults to <tt>true</tt>.
     */
    bool use_boundary;

    /**
     * The flag indicating whether the interior face integrator face() is to
     * be used in the loop. Defaults to <tt>true</tt>.
     */
    bool use_face;

    /**
     * The names of the input vectors. If this vector is nonempty, it can be
     * used by application programs to automatically select and verify the
     * input vectors used for integration.
     *
     * @note This variable is currently not used by the library, but it is
     * provided to help develop application programs.
     *
     * @deprecated Because the library itself does not use this field, it is
     *   better placed in derived classes.
     */
    DEAL_II_DEPRECATED
    std::vector<std::string> input_vector_names;

    /**
     * The names of the results produced. If this vector is nonempty, it can
     * be used by application programs to automatically assign names to output
     * values and/or verify the names of vectors.
     *
     * @note This variable is currently not used by the library, but it is
     * provided to help develop application programs.
     *
     * @deprecated Because the library itself does not use this field, it is
     *   better placed in derived classes.
     */
    DEAL_II_DEPRECATED
    std::vector<std::string> output_names;

    /**
     * This error is thrown if one of the virtual functions cell(),
     * boundary(), or face() is called without being overloaded in a derived
     * class. Consider setting #use_cell, #use_boundary, and #use_face to
     * false, respectively.
     *
     * @ingroup Exceptions
     */
    DeclException0(ExcPureFunction);
  };
} // namespace MeshWorker



DEAL_II_NAMESPACE_CLOSE

#endif
