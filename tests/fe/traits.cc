//----------------------------  fe_traits_test.cc  ---------------------------
//    traits.cc,v 1.2 2004/01/23 16:34:24 wolf Exp
//    Version: 
//
//    Copyright (C) 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_traits_test.cc  ---------------------------


#include "../tests.h"
#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_system.h>


template <int dim>
void check (const FiniteElement<dim> &fe)
{
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

				   // first check whether shape
				   // functions are primitive:
  deallog << "  Primitivity: ";
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    deallog << (fe.is_primitive(i) ? 1 : 0);
  deallog << std::endl;

  deallog << "  Overall primitivity: " << fe.is_primitive() << std::endl;

				   // then check n_nonzero_components
  deallog << "  n_nonzero_components: ";
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    deallog << fe.n_nonzero_components(i);
  deallog << std::endl;

				   // finally check component pattern
				   // for each shape function
  deallog << "  component pattern for each shape function:" << std::endl;
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      deallog << "    shape function i=" << i
	      << " [" << fe.system_to_component_index(i).first
	      << ','
	      << fe.system_to_component_index(i).second
	      << "]: ";
		       
      for (unsigned int j=0; j<fe.n_components(); ++j)
	{
	  Assert (fe.get_nonzero_components(i).size() == fe.n_components(),
		  ExcInternalError());
	  deallog << (fe.get_nonzero_components(i)[j] ? 1 : 0);
	};
      deallog << std::endl;
    };
}



template <int dim>
void check ()
{
				   // check usual Lagrange elements
  for (unsigned int p=1; p<3; ++p)
    {
      deallog << "Checking FE_Q<" << dim << ">(" << p << "): "
	      << std::endl;
      check (FE_Q<dim>(p));
    };

				   // check DG Lagrange elements
  for (unsigned int p=0; p<3; ++p)
    {
      deallog << "Checking FE_DGQ<" << dim << ">(" << p << "): "
	      << std::endl;
      check (FE_DGQ<dim>(p));
    };

				   // check DG-P elements
  for (unsigned int p=0; p<3; ++p)
    {
      deallog << "Checking FE_DGP<" << dim << ">(" << p << "): "
	      << std::endl;
      check (FE_DGP<dim>(p));
    };

				   // check systems of Q-elements
  for (unsigned int p=1; p<3; ++p)
    {
      deallog << "Checking FE_Q<" << dim << ">(" << p << ")^2: "
	      << std::endl;
      check (FESystem<dim> (FE_Q<dim>(p),2));
    };

				   // check systems of systems of
				   // Q-elements
  for (unsigned int p=1; p<3; ++p)
    {
      deallog << "Checking FE_Q<" << dim << ">(" << p << ")^2^2: "
	      << std::endl;
      check (FESystem<dim> (FESystem<dim> (FE_Q<dim>(p),2), 2));
    };
}

  


int main ()
{
  std::ofstream logfile("traits.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}

