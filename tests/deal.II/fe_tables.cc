// $Id$

#include <fe/fe_lib.lagrange.h>
#include <fe/fe_lib.dg.h>
#include <fe/fe_system.h>
#include <iomanip>
#include <base/logstream.h>

template<int dim>
inline void
print_fe_statistics(const FiniteElement<dim>& fe)
{
  deallog << "total_dofs" << " " << fe.total_dofs;
  deallog << ": vertex" << " " << fe.dofs_per_vertex;
  deallog << "  line" << " " << fe.dofs_per_line;
  deallog << "  quad" <<  " " <<fe.dofs_per_quad << endl;
  deallog << "n_transform_fct " << fe.n_transform_functions << endl;
  deallog << "n_components " << fe.n_components << endl;
  deallog.push("components");
  for (unsigned i=0;i<fe.total_dofs;++i)
    {
      pair<unsigned,unsigned> p = fe.system_to_component_index(i);
      deallog << "Index " << i << " ("
	      << p.first << "," << p.second << ")" << endl;
    }
  deallog.pop();
}

main()
{
  deallog.push("GeometryInfo");

  deallog.push("1D");
  deallog << " vertices: " << GeometryInfo<1>::vertices_per_cell
	  << " lines: " << GeometryInfo<1>::lines_per_cell
	  << " quads: " << GeometryInfo<1>::quads_per_cell
	  << " hexes: " << GeometryInfo<1>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("2D");
  deallog << " vertices: " << GeometryInfo<2>::vertices_per_cell
	  << " lines: " << GeometryInfo<2>::lines_per_cell
	  << " quads: " << GeometryInfo<2>::quads_per_cell
	  << " hexes: " << GeometryInfo<2>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("3D");
  deallog << " vertices: " << GeometryInfo<3>::vertices_per_cell
	  << " lines: " << GeometryInfo<3>::lines_per_cell
	  << " quads: " << GeometryInfo<3>::quads_per_cell
	  << " hexes: " << GeometryInfo<3>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.push("4D");
  deallog << " vertices: " << GeometryInfo<4>::vertices_per_cell
	  << " lines: " << GeometryInfo<4>::lines_per_cell
	  << " quads: " << GeometryInfo<4>::quads_per_cell
	  << " hexes: " << GeometryInfo<4>::hexes_per_cell
	  << endl;
  deallog.pop();
  
  deallog.pop();
  
  FEDGConstant<2>   dg0;
  FEDGLinear<2>   dg1;
  FELinear<2>       q1;
  FEQuadraticSub<2> q2;
  FECubicSub<2>     q3;

  FESystem<2> q2_3(q2,3);
  FESystem<2> q3_3(q3,3);
  
  
  deallog.push("Q0  ");
  print_fe_statistics(dg0);
  deallog.pop();
  deallog << endl;
  
  deallog.push("Q1  ");
  print_fe_statistics(q1);
  deallog.pop();
  deallog << endl;
  
  deallog.push("Q2  ");
  print_fe_statistics(q2);
  deallog.pop();
  deallog << endl;

  deallog.push("Q3  ");
  print_fe_statistics(q3);
  deallog.pop();
  deallog << endl;

  deallog.push("DG1 ");
  print_fe_statistics(dg1);
  deallog.pop();
  deallog << endl;

  {
    FESystem<2> q1_3(q1,3);
    deallog.push("Q1^3");
    print_fe_statistics(q1_3);
    deallog.pop();
  }
  deallog << endl;

  deallog.push("Q2^3");
  print_fe_statistics(q2_3);
  deallog.pop();
  deallog << endl;

  deallog.push("Q3^3");
  print_fe_statistics(q3_3);
  deallog.pop();
  deallog << endl;
  
  {
    FESystem<2> q1_2q0(q1,2,dg0,1);
    deallog.push("Q1^2 Q0");
    print_fe_statistics(q1_2q0);
    deallog.pop();
  }
  deallog << endl;
  
  {
    FESystem<2> q2_2q1(q2,2,q1,1);
    deallog.push("Q2^2 Q1");
    print_fe_statistics(q2_2q1);
    deallog.pop();
  }
  deallog << endl;

  {
    FESystem<2> q2_2dg1(q2,2,dg1,1);
    deallog.push("Q2^2 DG1");
    print_fe_statistics(q2_2dg1);
    deallog.pop();
  }
  deallog << endl;
}
