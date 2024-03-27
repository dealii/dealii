// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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


#ifndef dealii_function_unit_to_real_h
#define dealii_function_unit_to_real_h


#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>


DEAL_II_NAMESPACE_OPEN
namespace Functions
{
  /**
   * Wrapper class for transforming a function from real to unit space.
   *
   * @ingroup functions
   */
  template <int dim>
  class FunctionUnitToReal : public Function<dim>
  {
  public:
    FunctionUnitToReal(const Function<dim> &function);

    void
    set_active_cell(
      const typename Triangulation<dim>::active_cell_iterator &cell);

    double
    value(const Point<dim> & point,
          const unsigned int component = 0) const override;

    Tensor<1, dim>
    gradient(const Point<dim> & point,
             const unsigned int component = 0) const override;

    SymmetricTensor<2, dim>
    hessian(const Point<dim> & point,
            const unsigned int component = 0) const override;

  private:
    const SmartPointer<const Function<dim>> function;

    unsigned int cell_level, cell_index;

    SmartPointer<const Triangulation<dim>> triangulation;

    const MappingQ1<dim> mapping;

    /**
     * We use a FEValues object to compute the Jacobian and its derivatives. To
     * use it we need an element.
     */
    const FE_Nothing<dim> dummy_element;
  };

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif
