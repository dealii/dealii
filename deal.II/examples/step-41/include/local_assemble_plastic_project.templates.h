#include "../include/local_assemble_plastic_project.h"
#include <base/quadrature_lib.h>

template <int dim>
LocalAssemblePlasticProject<dim>::LocalAssemblePlasticProject() : 
fe(0, "Local Assemble Fe Pointer"),
  fe_v(0, "Local Assemble FeValues Pointer")
{
}


template <int dim>
void LocalAssemblePlasticProject<dim>::reinit (FiniteElement<dim> &myfe,
					       Table<2, SymmetricTensor<2,dim> > &my_plastic_strain) 
{ 
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  fe = &myfe;
  QGauss<dim> quadrature(2*fe->degree + 1);
  QGauss<dim-1> face_quadrature(2*fe->degree + 1);
  UpdateFlags flags (update_gradients |
		     update_q_points  |
		     update_JxW_values);
  fe_v = new FEValues<dim>(*fe, quadrature, flags);
    
  this->flags = assemble_rhs_cell;
    
  plastic_strain = my_plastic_strain;
}

template <int dim>
LocalAssemblePlasticProject<dim>::~LocalAssemblePlasticProject()
{
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  fe = 0;
}

template <int dim>
void LocalAssemblePlasticProject<dim>::parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Elastic Moduli");
  C.parse_parameters(prm);
  prm.leave_subsection();
}


template <int dim>
void LocalAssemblePlasticProject<dim>::assemble_rhs_term
(const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
 Vector<double> &cell_rhs)
{
  cell_rhs = 0;
  Assert(fe, ExcNotInitialized());

  fe_v->reinit(cell);
  
  unsigned int size = fe->n_components();
  unsigned int n_q_points = fe_v->n_quadrature_points;
  
  std::vector<Vector<double> > load_vector (n_q_points, Vector<double>(size) );
  std::vector<Point<dim> > points = fe_v->get_quadrature_points();  

  for (unsigned int qp=0; qp<n_q_points; ++qp) {
    
    SymmetricTensor<4,dim> C_qp = C(points[qp]);

    for (unsigned int i=0; i<fe_v->dofs_per_cell; ++i) {
      unsigned int comp_i = fe->system_to_component_index(i).first;
      
      for(unsigned int b=0; b<dim; ++b) {
	for(unsigned int m=0; m<dim; ++m) {
	  for(unsigned int n=0; n<dim; ++n) {
	   
	    cell_rhs(i) += ( C_qp[comp_i][b][m][n] *
			     plastic_strain(cell->index(), qp)[m][n] *
			     fe_v->shape_grad(i, qp)[b] *
			     fe_v->JxW(qp) );
	    
	  } //n
	} //m 
      } //b
    } //i
  } //qp
}


