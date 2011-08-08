//----------------------------  fe_data_test.cc  ---------------------------
//    fe_data_test.cc,v 1.14 2003/11/28 11:52:35 guido Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2010, 2011 by the deal.II authors
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

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomials_raviart_thomas.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_dg_vector.h>
#include <deal.II/fe/fe_face.h>

#include <deal.II/fe/fe_system.h>

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
  fe_datas.push_back(new FE_Q_Hierarchical<dim> (1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_Q_Hierarchical<dim> (2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FE_Q_Hierarchical<dim> (4));
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
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (2), 2));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  fe_datas.push_back(new FESystem<dim>(FE_Q<dim> (1), 2,
				       FE_Q<dim> (2), 1));
  deallog << (*fe_datas.rbegin())->get_name() << std::endl;
  

				   // Check vector elements in 2d and higher only
  if (dim>1)
    {
				       // Face elements
      fe_datas.push_back(new FE_FaceQ<dim> (0));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FE_FaceQ<dim> (1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FE_FaceQ<dim> (3));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
				   // Vector DG elements
      fe_datas.push_back(
	new FE_DGRaviartThomas<dim>(0));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(
	new FE_DGRaviartThomas<dim>(1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(
	new FE_DGNedelec<dim>(0));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(
	new FE_DGNedelec<dim>(1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;

				       // Hdiv elements
      FE_RaviartThomas<dim>* rt0 = new FE_RaviartThomas<dim>(0);
      fe_datas.push_back(rt0);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;

      FE_RaviartThomas<dim>* rt1 = new FE_RaviartThomas<dim>(1);
      fe_datas.push_back(rt1);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;

      fe_datas.push_back(new FE_RaviartThomas<dim>(2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FESystem<dim>(*rt1, 1,
					   FE_DGQ<dim> (1), 1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;

				       // Hcurl elements
      FE_Nedelec<dim>* ned0 = new FE_Nedelec<dim>(0);
      fe_datas.push_back(ned0);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      FE_Nedelec<dim>* ned1 = new FE_Nedelec<dim>(1);
      fe_datas.push_back(ned1);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
    }
  if (dim==2)
    {
      fe_datas.push_back(
	new FE_DGBDM<dim>(1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(
	new FE_DGBDM<dim>(2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      
      fe_datas.push_back(new FE_BDM<dim>(1));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FE_BDM<dim>(2));
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;     
    }
  if (dim>1)
    {
      FE_RaviartThomasNodal<dim>* rt0 = new FE_RaviartThomasNodal<dim>(0);
      FE_RaviartThomasNodal<dim>* rt1 = new FE_RaviartThomasNodal<dim>(1);
      fe_datas.push_back(rt0);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(rt1);
      deallog << (*fe_datas.rbegin())->get_name() << std::endl;
      fe_datas.push_back(new FESystem<dim>(*rt1, 1,
					   FE_DGQ<dim> (1), 1));
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
      deallog << "primitive=" << (fe_data->is_primitive() ? "yes" : "no") << std::endl
	      << "components=" << fe_data->components << std::endl
	      << "blocks=" << fe_data->n_blocks() << std::endl
	      << "degree=" << fe_data->tensor_degree() << std::endl
	      << "conformity=";
      if (fe_data->conforms(FiniteElementData<dim>::L2)) deallog << " L2";
      if (fe_data->conforms(FiniteElementData<dim>::Hcurl)) deallog << " Hcurl";
      if (fe_data->conforms(FiniteElementData<dim>::Hdiv)) deallog << " Hdiv";
      if (fe_data->conforms(FiniteElementData<dim>::H1)) deallog << " H1";
      if (fe_data->conforms(FiniteElementData<dim>::H2)) deallog << " H2";
      deallog << std::endl;
      deallog << "unit_support_points=" << fe_data->get_unit_support_points().size()
	      << std::endl;
      deallog << "unit_face_support_points=" << fe_data->get_unit_face_support_points().size()
	      << std::endl;
      deallog << "generalized_support_points=" << fe_data->get_generalized_support_points().size()
	      << std::endl;
      deallog << "generalized_face_support_points=" << fe_data->get_generalized_face_support_points().size()
	      << std::endl;
      
      deallog << "face_to_equivalent_cell_index:";
      for (unsigned int i=0;i<fe_data->dofs_per_face;++i)
	deallog << ' ' << fe_data->face_to_equivalent_cell_index(i);
      deallog << std::endl;
      for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
	{
	  deallog << "face_to_cell_index:";
	  for (unsigned int i=0;i<fe_data->dofs_per_face;++i)
	    deallog << ' ' << fe_data->face_to_cell_index(i, f);
	  deallog << std::endl;
	}
      
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

int main()
{
  std::ofstream logfile("fe_data_test/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
//  deallog.log_execution_time(true);
//  deallog.log_time_differences(true);
  
  test_fe_datas<1>();
  test_fe_datas<2>();
  test_fe_datas<3>();
}



