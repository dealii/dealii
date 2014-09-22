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
#include <deal.II/fe/fe_face.h>

// check
//   FE_FaceQ/P::subface_interpolation


std::string output_file_name = "output";


void
output_matrix (const FullMatrix<double> &m)
{
  if ((m.m() == 0) || (m.n() == 0))
    {
      deallog << "(Empty matrix)" << std::endl;
      return;
    }

  deallog << m.l1_norm() << ' ' << m.linfty_norm()
          << std::endl;
  if (m.m() == m.n())
    deallog << m.frobenius_norm() << std::endl;

  for (unsigned int i=0; i<std::min(m.m(),m.n()); ++i)
    deallog << m(i,i) << ' ' << m(i,std::min(m.m(),m.n())-i-1) << ' ';
  deallog << std::endl;
}


template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
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


template <int dim>
void
check (const unsigned int degree1,
       const unsigned int degree2)
{
  {
    deallog << "Checking FE_FaceQ in " << dim << "d:"
            << std::endl;

    FE_FaceQ<dim> fe1 (degree1), fe2 (degree2);

    check_this (fe1, fe2);
  }
  {
    deallog << "Checking FE_FaceP in " << dim << "d:"
            << std::endl;

    FE_FaceP<dim> fe1 (degree1), fe2 (degree2);

    check_this (fe1, fe2);
  }
}


template <int dim>
void
check ()
{
  check<dim>(0,0);
  check<dim>(0,1);
  check<dim>(0,2);
  check<dim>(0,3);
  check<dim>(1,2);
  check<dim>(1,3);
  check<dim>(2,3);
}



int
main()
{
  try
    {
      std::ofstream logfile(output_file_name.c_str());
      deallog << std::setprecision (6);
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      check<2>();
      check<3>();

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
