//----------------------------  fe_data_test.cc  ---------------------------
//    fe_data_test.cc,v 1.14 2003/11/28 11:52:35 guido Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_data_test.cc  ---------------------------


#include "../tests.h"
#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

//TODO: Find support_on_face problems for test-no. > 7
//TODO: Check support_on_face in 3D

template<int dim>
void test_fe_datas()
{
  std::vector<FiniteElement<dim> *> fe_datas;
  fe_datas.push_back(new FE_Q<dim> (1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_Q<dim> (2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_Q<dim> (4));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGQ<dim> (1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGQ<dim> (2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGQ<dim> (4));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGP<dim> (1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGP<dim> (2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_DGP<dim> (4));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (2), 2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (1), 2,
				       FE_Q<dim> (2), 1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  
				   // Check Raviart-Thomas in 2d only
  if (dim==2)
    {
      fe_datas.push_back(new FE_RaviartThomas<dim>(0));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FE_RaviartThomas<dim>(1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FE_RaviartThomas<dim>(2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
    }
  

				   // for dim==3 the constraints are
				   // only hardcoded for Q1-Q2
  if (dim!=3)
    {
      fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (3), 2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (1), 2,
					   FE_Q<dim> (3), 1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (4), 2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
    }

				   // have systems of systems, and
				   // construct hierarchies of
				   // subsequently weirder elements by
				   // taking each of them in turn as
				   // basis of others
  fe_datas.push_back (new FESystem<dim> (FESystem<dim> (FE_Q<dim>(1),2),2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back (new FESystem<dim> (FESystem<dim> (FE_Q<dim>(1),2),1,
					 FESystem<dim> (FE_DGQ<dim>(1),2),1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back (new FESystem<dim> (FESystem<dim> (FE_Q<dim>(1),1,
							FE_Q<dim>(2),1),1,
					 FESystem<dim> (FE_Q<dim>(2),2),1,
					 FESystem<dim> (FE_DGQ<dim>(2),2),1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back (new FESystem<dim> (*fe_datas[fe_datas.size()-3], 2,
					 *fe_datas[fe_datas.size()-2], 1,
					 *fe_datas[fe_datas.size()-1], 2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  
  deallog << std::endl << "dim=" << dim << std::endl;
  for (unsigned int n=0; n<fe_datas.size(); ++n)
    {
      FiniteElement<dim> *fe_data=fe_datas[n];
      deallog << "fe_data[" << n <<"]:" << fe_data->get_name() << std::endl;
      deallog << "dofs_per_vertex=" << fe_data->dofs_per_vertex << std::endl;
      deallog << "dofs_per_line=" << fe_data->dofs_per_line << std::endl;
      deallog << "dofs_per_quad=" << fe_data->dofs_per_quad << std::endl;
      deallog << "dofs_per_hex=" << fe_data->dofs_per_hex << std::endl;
      deallog << "first_line_index=" << fe_data->first_line_index << std::endl;
      deallog << "first_quad_index=" << fe_data->first_quad_index << std::endl;
      deallog << "first_hex_index=" << fe_data->first_hex_index << std::endl;
      deallog << "first_face_line_index=" << fe_data->first_face_line_index << std::endl;
      deallog << "first_face_quad_index=" << fe_data->first_face_quad_index << std::endl;
      deallog << "dofs_per_face=" << fe_data->dofs_per_face << std::endl;
      deallog << "dofs_per_cell=" << fe_data->dofs_per_cell << std::endl;
      deallog << "components=" << fe_data->components << std::endl
	      << "degree=" << fe_data->tensor_degree() << std::endl;

      for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
	{
	  deallog << "support on face " << f << ':';
	  for (unsigned int s=0;s<fe_data->dofs_per_cell;++s)
	    if (fe_datas[n]->has_support_on_face(s, f))
	      deallog << '\t' << s;
	  deallog << std::endl;
	}
      deallog << std::endl;
    }

				   // delete all FiniteElementDatas
  for (unsigned int i=0; i<fe_datas.size(); ++i)
    delete fe_datas[i];
}

int main(int,char)
{
  std::ofstream logfile("fe_data_test.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
//  deallog.log_execution_time(true);
  
  test_fe_datas<1>();
  test_fe_datas<2>();
  test_fe_datas<3>();
}



