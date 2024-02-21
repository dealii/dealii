// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FE::subface_interpolation



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  if (dim == 1)
    return;

  // check all combinations of fe1 and fe2
  for (unsigned int subface = 0;
       subface < GeometryInfo<dim>::max_children_per_face;
       ++subface)
    {
      FullMatrix<double> face_constraints;
      try
        {
          face_constraints.reinit(fe1.dofs_per_face, fe1.dofs_per_face);
          fe1.get_subface_interpolation_matrix(fe1, subface, face_constraints);

          deallog << fe1.get_name() << "  vs.  " << fe1.get_name() << std::endl;
          output_matrix(face_constraints);
        }
      catch (...)
        {}

      try
        {
          face_constraints.reinit(fe2.dofs_per_face, fe2.dofs_per_face);
          fe2.get_subface_interpolation_matrix(fe2, subface, face_constraints);

          deallog << fe2.get_name() << "  vs.  " << fe2.get_name() << std::endl;
          output_matrix(face_constraints);
        }
      catch (...)
        {}

      if (fe1.dofs_per_face <= fe2.dofs_per_face)
        try
          {
            face_constraints.reinit(fe2.dofs_per_face, fe1.dofs_per_face);
            fe1.get_subface_interpolation_matrix(fe2,
                                                 subface,
                                                 face_constraints);

            deallog << fe1.get_name() << "  vs.  " << fe2.get_name()
                    << std::endl;
            output_matrix(face_constraints);
          }
        catch (...)
          {}

      if (fe2.dofs_per_face <= fe1.dofs_per_face)
        try
          {
            face_constraints.reinit(fe1.dofs_per_face, fe2.dofs_per_face);
            fe2.get_subface_interpolation_matrix(fe1,
                                                 subface,
                                                 face_constraints);

            deallog << fe2.get_name() << "  vs.  " << fe1.get_name()
                    << std::endl;
            output_matrix(face_constraints);
          }
        catch (...)
          {}
    }
}
