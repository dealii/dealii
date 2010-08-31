#include "../include/local_assemble_elastic_matrix.h"
#include <base/quadrature_lib.h>

template <int dim>
LocalAssembleElasticMatrix<dim>::LocalAssembleElasticMatrix() : 
fe(0, "Local Assemble Fe Pointer"),
  fe_v(0, "Local Assemble FeValues Pointer")
{
}


template <int dim>
void LocalAssembleElasticMatrix<dim>::reinit (FiniteElement<dim>& myfe) 
{ 
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  fe = &myfe;
  QGauss<dim> quadrature(2*fe->degree + 1);
  QGauss<dim-1> face_quadrature(2*fe->degree + 1);
  UpdateFlags flags (update_values |
		     update_gradients |
		     update_q_points  |
		     update_JxW_values);
  fe_v = new FEValues<dim>(*fe, quadrature, flags);
    
  this->flags = assemble_cell;
    
}

template <int dim>
LocalAssembleElasticMatrix<dim>::~LocalAssembleElasticMatrix()
{
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  fe = 0;
}

template <int dim>
void LocalAssembleElasticMatrix<dim>::parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Elastic Moduli");
  C.parse_parameters(prm);
  prm.leave_subsection();
}

template <int dim>
void LocalAssembleElasticMatrix<dim>::assemble_cell_term
(const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
 FullMatrix<double> &cell_m)
{
  cell_m = 0;
  Assert(fe, ExcNotInitialized());
  
  fe_v->reinit(cell);

  std::vector<Point<dim> > points(fe_v->n_quadrature_points);
  points = fe_v->get_quadrature_points();
  unsigned int comp_i, comp_j;  

  for(unsigned int point =0; point<fe_v->n_quadrature_points; ++point) {
    SymmetricTensor<4,dim> C_qp = C(points[point]);
    
    for (unsigned int i=0; i<fe_v->dofs_per_cell; ++i) {
      comp_i = fe->system_to_component_index(i).first;
      for (unsigned int j=0; j<fe_v->dofs_per_cell; ++j) {
	comp_j = fe->system_to_component_index(j).first;
	
	for(unsigned int b=0; b<dim; ++b) {
	  for(unsigned int n=0; n<dim; ++n) {
	    cell_m(i,j) +=  C_qp[comp_i][b][comp_j][n] * 
	      fe_v->shape_grad(j, point)[n] *
	      fe_v->shape_grad(i, point)[b] *
	      fe_v->JxW(point);
	  } //n
	}//b
      } //j
    }//i
  }//qp
}
