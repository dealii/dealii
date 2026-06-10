// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


#ifndef dealii_dof_output_operator_h
#define dealii_dof_output_operator_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/base/event.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  /**
   * An output operator writing a separate file in each step and writing the
   * vectors as finite element functions with respect to a given DoFHandler.
   *
   * @deprecated The caller should set up and manage a DataOut object instead.
   */
  template <typename VectorType, int dim, int spacedim = dim>
  class DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "The caller should set up and manage a DataOut object instead.")
    DoFOutputOperator : public OutputOperator<VectorType>
  {
  public:
    /*
     * Constructor. The <tt>filename</tt> is the common base name of
     * all files and the argument <tt>digits</tt> should be the number
     * of digits of the highest number in the sequence. File names by
     * default have the form "outputNNN" with NNN the number set by the
     * last step command. Numbers with less digits are filled with
     * zeros from the left.
     */
    DoFOutputOperator(const std::string &filename_base = "output",
                      const unsigned int digits        = 3);

    void
    parse_parameters(ParameterHandler &param);
    void
    initialize(const DoFHandler<dim, spacedim> &dof_handler);

    virtual OutputOperator<VectorType> &
    operator<<(const AnyData &vectors) override;

  private:
    ObserverPointer<const DoFHandler<dim, spacedim>,
                    DoFOutputOperator<VectorType, dim, spacedim>>
      dof;

    const std::string  filename_base;
    const unsigned int digits;

    DataOut<dim> out;
  };

  /*-------------------- Functions: DoFOutputOperator ------------------------*/


  template <typename VectorType, int dim, int spacedim>
  inline void
  DoFOutputOperator<VectorType, dim, spacedim>::initialize(
    const DoFHandler<dim, spacedim> &dof_handler)
  {
    dof = &dof_handler;
  }



  template <typename VectorType, int dim, int spacedim>
  DoFOutputOperator<VectorType, dim, spacedim>::DoFOutputOperator(
    const std::string &filename_base,
    const unsigned int digits)
    : filename_base(filename_base)
    , digits(digits)
  {
    out.set_default_format(DataOutBase::gnuplot);
  }



  template <typename VectorType, int dim, int spacedim>
  inline void
  DoFOutputOperator<VectorType, dim, spacedim>::parse_parameters(
    ParameterHandler &param)
  {
    out.parse_parameters(param);
  }



  template <typename VectorType, int dim, int spacedim>
  inline OutputOperator<VectorType> &
  DoFOutputOperator<VectorType, dim, spacedim>::operator<<(const AnyData &data)
  {
    Assert((dof != nullptr), ExcNotInitialized());
    out.attach_dof_handler(*dof);
    for (unsigned int i = 0; i < data.size(); ++i)
      {
        const VectorType *p = data.try_read_ptr<VectorType>(i);
        if (p != nullptr)
          {
            out.add_data_vector(*p, data.name(i));
          }
      }
    std::ostringstream streamOut;
    streamOut << filename_base << std::setw(digits) << std::setfill('0')
              << this->step << out.default_suffix();
    std::ofstream out_filename(streamOut.str());
    out.build_patches();
    out.write(out_filename);
    out.clear();
    return *this;
  }
} // namespace Algorithms


DEAL_II_NAMESPACE_CLOSE

#endif
