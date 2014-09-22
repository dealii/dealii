// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__integrators_patches_h
#define __deal2__integrators_patches_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/meshworker/dof_info.h>

DEAL_II_NAMESPACE_OPEN

namespace LocalIntegrators
{
  /**
   * @brief Integrators writing patches with values in quadrature points
   *
   * @author Guido Kanschat
   * @date 2011
   */
  namespace Patches
  {
    template <int dim>
    inline
    void
    points_and_values(
      Table<2, double> &result,
      const FEValuesBase<dim> &fe,
      const VectorSlice<const std::vector<std::vector<double> > > &input)
    {
      const unsigned int n_comp = fe.get_fe().n_components();
      AssertVectorVectorDimension(input, n_comp, fe.n_quadrature_points);
      AssertDimension(result.n_rows(), fe.n_quadrature_points);
      AssertDimension(result.n_cols(), n_comp+dim);

      for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
        {
          for (unsigned int d=0; d<dim; ++d)
            result(k,d) = fe.quadrature_point(k)[d];
          for (unsigned int i=0; i<n_comp; ++i)
            result(k,dim+i) = input[i][k];
        }
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
