#include "../include/local_assemble_scalar_project.h"
#include <base/quadrature_lib.h>

template <int dim>
LocalAssembleScalarProject<dim>::LocalAssembleScalarProject() : 
fe(0, "Local Assemble Fe Pointer"),
  fe_v(0, "Local Assemble FeValues Pointer")
{
}


template <int dim>
void LocalAssembleScalarProject<dim>::reinit (FiniteElement<dim> &myfe,
					      Table<2, double> &my_hard,
					      Table<2, double> &my_iter) 
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
		     update_q_points  |
		     update_JxW_values);
  fe_v = new FEValues<dim>(*fe, quadrature, flags);
    
  this->flags = assemble_rhs_cell;
    
  hard_table = my_hard;
  iter_table = my_iter;
}

template <int dim>
LocalAssembleScalarProject<dim>::~LocalAssembleScalarProject()
{
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  fe = 0;
}


template <int dim>
void LocalAssembleScalarProject<dim>::assemble_rhs_term
(const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
 Vector<double> &cell_rhs)
{
  cell_rhs = 0;
  Assert(fe, ExcNotInitialized());

  fe_v->reinit(cell);
  
  unsigned int n_q_points = fe_v->n_quadrature_points; 

  for (unsigned int qp=0; qp<n_q_points; ++qp) {

    for (unsigned int i=0; i<fe_v->dofs_per_cell; ++i) {
      unsigned int comp_i = fe->system_to_component_index(i).first;
      
      if(comp_i == 0) {
	
	cell_rhs(i) += ( hard_table(cell->index(), qp) *
			 fe_v->shape_value(i, qp) *
			 fe_v->JxW(qp) );
      }

      if(comp_i == 1) {
	
	cell_rhs(i) += ( iter_table(cell->index(), qp) *
			 fe_v->shape_value(i, qp) *
			 fe_v->JxW(qp) );
      }


    } //i
  } //qp
}


