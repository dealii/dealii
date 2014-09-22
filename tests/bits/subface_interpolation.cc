// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


#include "../tests.h"
#include "fe_tools_common.h"

// check
//   FE::subface_interpolation


std::string output_file_name = "output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  if (dim == 1)
    return;

  // check all combinations of fe1 and fe2
  for (unsigned int subface=0; subface<GeometryInfo<dim>::max_children_per_face; ++subface)
    {
      FullMatrix<double> face_constraints;
      try
        {
          face_constraints.reinit (fe1.dofs_per_face,
                                   fe1.dofs_per_face);
          fe1.get_subface_interpolation_matrix (fe1, subface, face_constraints);

          deallog << fe1.get_name()
                  << "  vs.  "
                  << fe1.get_name()
                  << std::endl;
          output_matrix (face_constraints);
        }
      catch (...)
        {
        }

      try
        {
          face_constraints.reinit (fe2.dofs_per_face,
                                   fe2.dofs_per_face);
          fe2.get_subface_interpolation_matrix (fe2, subface, face_constraints);

          deallog << fe2.get_name()
                  << "  vs.  "
                  << fe2.get_name()
                  << std::endl;
          output_matrix (face_constraints);
        }
      catch (...)
        {
        }

      if (fe1.dofs_per_face <= fe2.dofs_per_face)
        try
          {
            face_constraints.reinit (fe2.dofs_per_face,
                                     fe1.dofs_per_face);
            fe1.get_subface_interpolation_matrix (fe2, subface, face_constraints);

            deallog << fe1.get_name()
                    << "  vs.  "
                    << fe2.get_name()
                    << std::endl;
            output_matrix (face_constraints);
          }
        catch (...)
          {
          }

      if (fe2.dofs_per_face <= fe1.dofs_per_face)
        try
          {
            face_constraints.reinit (fe1.dofs_per_face,
                                     fe2.dofs_per_face);
            fe2.get_subface_interpolation_matrix (fe1, subface, face_constraints);

            deallog << fe2.get_name()
                    << "  vs.  "
                    << fe1.get_name()
                    << std::endl;
            output_matrix (face_constraints);
          }
        catch (...)
          {
          }
    }
}

