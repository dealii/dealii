//----------------------------  anna_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors and Anna Schneebeli
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  anna_1.cc  ---------------------------


// check has_support_on_face for some elements
//
// this program is a modified version of one by Roy Stogner, 
// University of Texas at Austin

#include <base/logstream.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_system.h>
#include <fstream>


template <int dim>
void check (const FiniteElement<dim> &fe) {
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
         face++)
      for (unsigned int i = 0; i < fe.dofs_per_cell; i++)
        if (fe.has_support_on_face(i, face))
          deallog << "Basis function " << i
                  << " has support on face " << face << std::endl;
}


#define check_el(fe) { deallog << #fe << std::endl; check(fe); }

template <int dim>
void check () 
{
  deallog << "************ dim = " << dim << std::endl;
  check_el (FE_Q<dim>(1));
  check_el (FE_Q<dim>(2));
  check_el (FE_Q<dim>(3));

  check_el (FE_DGQ<dim>(0));
  check_el (FE_DGQ<dim>(1));
  check_el (FE_DGQ<dim>(2));
  check_el (FE_DGQ<dim>(3));

  check_el (FE_DGP<dim>(0));
  check_el (FE_DGP<dim>(1));
  check_el (FE_DGP<dim>(2));

  if (dim > 1)
    check_el (FE_Nedelec<dim>(1));

  check_el (FESystem<dim> (FE_Q<dim>(1), 2));
  check_el (FESystem<dim> (FE_Q<dim>(1), 1,
                           FE_Q<dim>(1), 1));

  if (dim > 1)
    check_el (FESystem<dim> (FE_Nedelec<dim>(1), 2));
};


int main () 
{
  std::ofstream logfile("roy_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

//  check<1> ();
  check<2> ();
//  check<3> ();  
  
  return 0;
};
