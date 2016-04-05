/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */

//Check whether the restricted prolongation is the identity for nested FE spaces
#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/fe_q_bubbles.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim, int spacedim>
void
check(const FiniteElement<dim, spacedim> &fe,
      const bool isotropic_only = false,
      unsigned int nested_size = 0)
{
  deallog<<fe.get_name()<<std::endl;
  const unsigned int dpc = fe.dofs_per_cell;

  if (nested_size==0)
    nested_size=dpc;

  // loop over all possible refinement cases
  unsigned int ref_case = (isotropic_only)
                          ? RefinementCase<dim>::isotropic_refinement
                          : RefinementCase<dim>::cut_x;
  for (; ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
    {
      deallog << "RefinementCase " << ref_case << std::endl;
      // create a respective refinement on the triangulation
      Triangulation<dim, spacedim> tr;
      GridGenerator::hyper_cube (tr, 0, 1);
      tr.begin_active()->set_refine_flag(RefinementCase<dim>(ref_case));
      tr.execute_coarsening_and_refinement();

      DoFHandler<dim, spacedim> dh(tr);
      dh.distribute_dofs(fe);

      const unsigned int n_dofs = dh.n_dofs();

      FullMatrix<double> restriction_global(dpc,n_dofs);
      FullMatrix<double> prolongation_global(n_dofs,dpc);

      std::vector<types::global_dof_index> ldi (dpc);

      //now create the matrix coarse to fine (prolongation)
      //and fine to coarse (restriction) with respect to all dofs
      unsigned int child_no = 0;
      typename dealii::DoFHandler<dim, spacedim>::active_cell_iterator cell
        = dh.begin_active();
      for (; cell!=dh.end(); ++cell, ++child_no)
        {
          FullMatrix<double> restriction_local
            = fe.get_restriction_matrix(child_no, RefinementCase<dim>(ref_case));
          FullMatrix<double> prolongation_local
            = fe.get_prolongation_matrix(child_no, RefinementCase<dim>(ref_case));

          cell->get_dof_indices(ldi);

          for (unsigned int j=0; j<dpc; ++j)
            {
              const bool add = fe.restriction_is_additive(j);
              for (unsigned int i=0; i<dpc; ++i)
                {
                  prolongation_global(ldi[i],j)=prolongation_local(i,j);
                  if (add)
                    restriction_global(j,ldi[i])+=restriction_local(j,i);
                  else if (restriction_local(j,i) != 0)
                    restriction_global(j,ldi[i]) =restriction_local(j,i);
                }
            }
        }

      FullMatrix<double> result(dpc);
      restriction_global.mmult(result, prolongation_global);

//     deallog << std::endl
//             << "Restriction"
//             << std::endl;
//     for (unsigned int i=0; i<dpc; ++i)
//     {
//       for (unsigned int j=0; j<n_dofs; ++j)
//       {
//         deallog << std::setw(8) << std::setprecision(8)
//                 << restriction_global(i,j);
//         deallog << ' ';
//       }
//       deallog<<std::endl;
//     }
//     deallog<<std::endl;
//
//     deallog << std::endl
//             << "Prolongation"
//             << std::endl;
//     for (unsigned int i=0; i<n_dofs; ++i)
//     {
//       for (unsigned int j=0; j<dpc; ++j)
//       {
//         deallog << std::setw(8) << std::setprecision(8)
//                 << prolongation_global(i,j);
//         deallog << ' ';
//       }
//       deallog<<std::endl;
//     }
//     deallog<<std::endl;

      bool is_identity = true;
      for (unsigned int i=0; i<nested_size; ++i)
        for (unsigned int j=0; j<nested_size; ++j)
          {
            const double expected = (i==j)?1.:0.;
            if (std::fabs(result(i,j)-expected)>1.e-12)
              {
                deallog << i << " " << j << " " << result(i,j) << std::endl;
                is_identity = false;
              }
          }

      if (is_identity)
        deallog << "OK" << std::endl;
    }
}



int main ()
{
  initlog();
  deallog.depth_file (1);
  deallog.threshold_double(1.e-10);
  for (unsigned int i=1; i<=3; ++i)
    {
      {
        FE_Q<2> fe(i);
        check(fe);
      }
      {
        FE_Q<3> fe(i);
        check(fe);
      }
      {
        FE_Q_DG0<2> fe(i);
        check(fe);
      }
      {
        FE_Q_DG0<3> fe(i);
        check(fe);
      }
      {
        FE_Q_Bubbles<2> fe(i);
        const unsigned int n_q_dofs = FE_Q<2>(i).dofs_per_cell;
        check(fe, false, n_q_dofs);
      }
      {
        FE_Q_Bubbles<3> fe(i);
        const unsigned int n_q_dofs = FE_Q<3>(i).dofs_per_cell;
        check(fe, false, n_q_dofs);
      }
    }

  return 0;
}
