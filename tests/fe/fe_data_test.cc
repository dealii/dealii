//----------------------------  fe_data_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_data_test.cc  ---------------------------


#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>


template<int dim>
void test_fe_datas()
{
  vector<FiniteElement<dim> *> fe_datas;
  fe_datas.push_back(new FE_Q<dim> (1));
  fe_datas.push_back(new FE_Q<dim> (4));
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (2), 2));
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (1), 2,
				       FE_Q<dim> (2), 1));
				   // for this an assertion thrown in
				   // FESystem::build_interface_constraints
				   // for the following calls for dim=3
//  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (3), 2));
//  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (1), 2,
//				       FE_Q<dim> (3), 1));
//  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (4), 2));
  
  deallog << endl << "dim=" << dim << endl;
  for (unsigned int n=0; n<fe_datas.size(); ++n)
    {
      FiniteElementData<dim> *fe_data=fe_datas[n];
      deallog << "fe_data[" << n <<"]:" << endl;
      deallog << "dofs_per_vertex=" << fe_data->dofs_per_vertex << endl;
      deallog << "dofs_per_line=" << fe_data->dofs_per_line << endl;
      deallog << "dofs_per_quad=" << fe_data->dofs_per_quad << endl;
      deallog << "dofs_per_hex=" << fe_data->dofs_per_hex << endl;
      deallog << "first_line_index=" << fe_data->first_line_index << endl;
      deallog << "first_quad_index=" << fe_data->first_quad_index << endl;
      deallog << "first_hex_index=" << fe_data->first_hex_index << endl;
      deallog << "first_face_line_index=" << fe_data->first_face_line_index << endl;
      deallog << "first_face_quad_index=" << fe_data->first_face_quad_index << endl;
      deallog << "dofs_per_face=" << fe_data->dofs_per_face << endl;
      deallog << "dofs_per_cell=" << fe_data->dofs_per_cell << endl;
      deallog << "components=" << fe_data->components << endl;
    }

				   // delete all FiniteElementDatas
  for (unsigned int i=0; i<fe_datas.size(); ++i)
    delete fe_datas[i];
}

int main(int,char)
{
  ofstream logfile("fe_data_test.dat");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test_fe_datas<1>();
  test_fe_datas<2>();
  test_fe_datas<3>();
}



