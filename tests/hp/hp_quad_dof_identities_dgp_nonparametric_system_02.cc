//----------------------------  hp_quad_dof_identities_dgp_nonparametric_system_02.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_quad_dof_identities_dgp_nonparametric_system_02.cc  ---------------------------


// check FESystem(FE_DGPNonparametric)::hp_quad_dof_identities, but with a different
// arrangement of base elements and multiplicities than in the _01 test


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_system.h>

#include <fstream>


template <int dim>
void test ()
{
  hp::FECollection<dim> fe_collection;
  for (unsigned int i=1; i<8-dim; ++i)
    {
				       // add the system three times, with
				       // different numbers of base elements
				       // and multiplicities
      fe_collection.push_back (FESystem<dim>(FE_DGPNonparametric<dim>(i),3));
      fe_collection.push_back (FESystem<dim>(FE_DGPNonparametric<dim>(i),2,
					     FE_DGPNonparametric<dim>(i),1));
      fe_collection.push_back (FESystem<dim>(FE_DGPNonparametric<dim>(i),1,
					     FE_DGPNonparametric<dim>(i),2));
    }

  for (unsigned int i=0; i<fe_collection.size(); ++i)
    for (unsigned int j=0; j<fe_collection.size(); ++j)
      {
	const std::vector<std::pair<unsigned int, unsigned int> >
	  identities = fe_collection[i].hp_quad_dof_identities (fe_collection[j]);

	deallog << "Identities for "
		<< fe_collection[i].get_name() << " and "
		<< fe_collection[j].get_name() << ": "
		<< identities.size()
		<< std::endl;
	
	for (unsigned int k=0; k<identities.size(); ++k)
	  {
	    Assert (identities[k].first < fe_collection[i].dofs_per_quad,
		    ExcInternalError());
	    Assert (identities[k].second < fe_collection[j].dofs_per_quad,
		    ExcInternalError());
	    
	    deallog << identities[k].first << ' '
		    << identities[k].second
		    << std::endl;
	  }
      }
}



int main ()
{
  std::ofstream logfile("hp_quad_dof_identities_dgp_nonparametric_system_02/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
