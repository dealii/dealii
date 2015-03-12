// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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



// check FESystem(FE_Q)::hp_quad_dof_identities, but with a different
// arrangement of base elements and multiplicities than in the _01 test


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/fe/fe_q.h>
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
      fe_collection.push_back (FESystem<dim>(FE_Q<dim>(i),3));
      fe_collection.push_back (FESystem<dim>(FE_Q<dim>(i),2,
                                             FE_Q<dim>(i),1));
      fe_collection.push_back (FESystem<dim>(FE_Q<dim>(i),1,
                                             FE_Q<dim>(i),2));
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
  std::ofstream logfile("output");
  logfile.precision(2);

  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2> ();
  test<3> ();

  deallog << "OK" << std::endl;
}
