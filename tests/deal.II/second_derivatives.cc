// 	$Id$	
// test for correctness of gradients on a given cell

// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <grid/tria.h>
#include <grid/tria_boundary.h>
#include <grid/dof.h>
#include <fe/fe_values.h>
#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <lac/vector.h>




int main () {
  Triangulation<2> tria;
  tria.create_hypercube (0,1);

  FELinear<2> fe;
  DoFHandler<2> dof(&tria);
  dof.distribute_dofs(fe);

  StraightBoundary<2> b;
  QTrapez<2> q;
  FEValues<2> fevalues(fe,q,update_second_derivatives);
  
  
  Vector<double> val(4);

  cout << "------------------------------------------------------------" << endl
       << "Testing transformation of 2nd derivatives of shape function:" << endl;
  
				   // test for each of the four
				   // shape functions. first loop:
				   // unit cell, second loop:
				   // one vertex moved
  bool testcase_succeeded = true;
  for (unsigned int loop=0; loop<2; ++loop)
    {
      cout << "Test loop: " << loop << endl
	   << "-----------" << endl;
	  
      				   // move one vertex of the only cell
      if (loop==1)
	tria.begin_active()->vertex(2)(0) = 2;
      fevalues.reinit (dof.begin_active(),b);
      
      for (unsigned int vertex=0; vertex<4; ++vertex)
	{
	  val.clear ();
	  val(vertex) = 1;
	  
	  vector<Tensor<2,2> > derivs(4);
	  fevalues.get_function_2nd_derivatives (val, derivs);
	  
	  cout << "  Vertex " << vertex << ": ";
	  switch (loop)
	    {
	      case 0:
		    for (unsigned int point=0; point<4; ++point)
		      {
							 // on the unit cell,
			bool ok=true;

			if ((derivs[point][0][0] != 0) ||
			    (derivs[point][1][1] != 0))
			  ok = false;

			switch (vertex)
			  {
			    case 0:
			    case 2:
				  if ((derivs[point][0][1] != 1) ||
				      (derivs[point][1][0] != 1))
				    ok = false;
				  break;
			    case 1:
			    case 3:
				  if ((derivs[point][0][1] != -1) ||
				      (derivs[point][1][0] != -1))
				    ok = false;
				  break;
			  };  
			
			if (!ok)
			  {
			    testcase_succeeded = false;
			    cout << "WRONG! ";
			  }
			else
			  cout << "OK ";
		      };
		    break;

	      case 1:
						     // are these numbers ok?
						     // I have never checked
						     // them, exept for the
						     // symmetry of the tensors
						     // and some of the signs of
						     // the entries! Someone
						     // should really try to
						     // verify at least
						     // some of the numbers
		    static double _true_value[4][4][2][2]
		      = { { {{0,1},
			     {1,0}},
			    {{0,0},
			     {0,0}},
			    {{0,1},
			     {1,-2}},
			    {{0,0},
			     {0,0}}  },
			  { {{0,-1},
			     {-1,0}},
			    {{0,0},
			     {0,0}},
			    {{0,-1},
			     {-1,2}},
			    {{0,0},
			     {0,0}}  },
			  { {{0,0},
			     {0,0}},
			    {{0,-.25},
			     {-.25,0}},
			    {{0,0},
			     {0,0}},
			    {{0,-.25},
			     {-.25,0.5}}  },
			  { {{0,0},
			     {0,0}},
			    {{0,.25},
			     {.25,0}},
			    {{0,0},
			     {0,0}},
			    {{0,.25},
			     {.25,-.5}}  }  };
		    static Tensor<2,2> true_value[4][4]
		      = {{ Tensor<2,2>(_true_value[0][0]),
			   Tensor<2,2>(_true_value[0][1]),
			   Tensor<2,2>(_true_value[0][2]),
			   Tensor<2,2>(_true_value[0][3])  },
			 { Tensor<2,2>(_true_value[1][0]),
			   Tensor<2,2>(_true_value[1][1]),
			   Tensor<2,2>(_true_value[1][2]),
			   Tensor<2,2>(_true_value[1][3])  },
			 { Tensor<2,2>(_true_value[2][0]),
			   Tensor<2,2>(_true_value[2][1]),
			   Tensor<2,2>(_true_value[2][2]),
			   Tensor<2,2>(_true_value[2][3])  },
			 { Tensor<2,2>(_true_value[3][0]),
			   Tensor<2,2>(_true_value[3][1]),
			   Tensor<2,2>(_true_value[3][2]),
			   Tensor<2,2>(_true_value[3][3])  }};
		    
		    for (unsigned int point=0; point<4; ++point)
		      {
			if (derivs[point] != true_value[vertex][point])
			  {
			    testcase_succeeded = false;
			    cout << "WRONG! ";
			  }
			else
			  cout << "OK ";
		      };
		    break;
	    };

	  cout << endl;
	};
    };
  

  cout << "------------------------------------------------------" << endl;

  if (testcase_succeeded)
    return 0;
  else
    return 1;
};
