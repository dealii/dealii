//----------------------------  unit_support_points.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  unit_support_points.cc  ---------------------------


// check equivalence of get_unit_support_points() and
// unit_support_point(), as well as the case where an element does not
// have support points

#include <base/logstream.h>
#include <fe/fe_system.h>		
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_nedelec.h>
#include <fstream>


template <int dim>
void check_cell1 (const FiniteElement<dim> &fe) 
{
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    Assert (fe.get_unit_support_points()[i] ==
            fe.unit_support_point(i),
            ExcInternalError());
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void check_face1 (const FiniteElement<dim> &fe) 
{
  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
    Assert (fe.get_unit_face_support_points()[i] ==
            fe.unit_face_support_point(i),
            ExcInternalError());
  deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void check_face1 (const FiniteElement<1> &) 
{}



template <int dim>
void check1 (const FiniteElement<dim> &fe) 
{
  check_cell1 (fe);
  check_face1 (fe);
}


template <int dim>
void check_cell2 (const FiniteElement<dim> &fe,
                  const unsigned int        comp) 
{
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    if (fe.system_to_component_index(i).first == comp)
      deallog << i << " " << fe.unit_support_point(i)
              << std::endl;
  deallog << "dim=" << dim << ", cell=ok" << std::endl;
}


template <int dim>
void check_face2 (const FiniteElement<dim> &fe,
                  const unsigned int        comp) 
{
      for (unsigned int i=0; i<fe.dofs_per_face; ++i)
        if (fe.face_system_to_component_index(i).first == comp)
          deallog << i << " " << fe.unit_face_support_point(i)
                  << std::endl;
      deallog << "dim=" << dim << ", face=ok" << std::endl;
}


void check_face2 (const FiniteElement<1> &,
                  const unsigned int) 
{}



template <int dim>
void check2 (const FiniteElement<dim> &fe,
             const unsigned int        comp) 
{
  check_cell2 (fe, comp);
  check_face2 (fe, comp);
}



template <int dim>
void check () 
{
  check1 (FE_Q<dim>(2));
  check1 (FE_DGQ<dim>(2));
  check1 (FESystem<dim> (FE_Q<dim> (2), 2,
                         FE_DGQ<dim> (2),1));

  check1 (FESystem<dim> (FE_Q<dim> (2), 2,
                         FE_DGQ<dim> (2),1));

  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGQ<dim> (2), 1),
          1);
  check2 (FESystem<dim> (FE_Q<dim> (2), 1,
                         FE_DGP<dim> (2), 1),
          0);
}

    

int main () 
{
  std::ofstream logfile("unit_support_points.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check<1> ();
  check<2> ();
  check<3> ();
  return 0;
}
