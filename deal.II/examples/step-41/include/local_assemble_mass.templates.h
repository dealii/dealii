#include "../include/local_assemble_mass.h"
#include "utilities.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

template <int dim, typename DH>
LocalAssembleMass<dim, DH>::LocalAssembleMass() : 
    fe(0, "Local Assemble Fe Pointer"),
    fe_v(0, "Local Assemble FeValues Pointer")
{
}


template <int dim, typename DH>
void LocalAssembleMass<dim, DH>::reinit (FiniteElement<dim> &myfe,
					   Function<dim> &f) 
{ 
    smart_delete(fe_v);
    fe = &myfe;
    QGauss<dim> quadrature(2*fe->degree + 1);
    UpdateFlags flags (update_values |
		       update_q_points  |
		       update_JxW_values);
    fe_v = new FEValues<dim>(*fe, quadrature, flags);
    
    this->flags = 
	assemble_cell|
	assemble_rhs_cell; 
    
    rhs = &f;
}

template <int dim, typename DH>
LocalAssembleMass<dim, DH>::~LocalAssembleMass()
{
    smart_delete(fe_v);
    fe = 0;
}

template <int dim, typename DH>
void LocalAssembleMass<dim, DH>::assemble_cell_term
(const typename DH::active_cell_iterator& cell, 
 FullMatrix<double> &cell_m)
{
  cell_m = 0;
  Assert(fe, ExcNotInitialized());
  
  fe_v->reinit(cell);
  unsigned int comp_i = 0, comp_j = 0;
  for (unsigned int i=0; i<fe_v->dofs_per_cell; ++i) {
      comp_i = fe->system_to_component_index(i).first;
      for (unsigned int j=0; j<fe_v->dofs_per_cell; ++j) {      
	  comp_j = fe->system_to_component_index(j).first;
	  if(comp_i == comp_j) 
	      for(unsigned int q_point =0; q_point<fe_v->n_quadrature_points; ++q_point) {
		  cell_m(i,j) +=  ( fe_v->shape_value(j,q_point) *
				    fe_v->shape_value(i,q_point) * 
				    fe_v->JxW(q_point) );
	      }
      }
  }
}


template <int dim, typename DH>
void LocalAssembleMass<dim, DH>::assemble_rhs_term
(const typename DH::active_cell_iterator& cell, 
 Vector<double> &cell_rhs)
{
  cell_rhs = 0;
  Assert(fe, ExcNotInitialized());

  fe_v->reinit(cell);
  
  unsigned int size = fe->n_components();
  unsigned int n_q_points = fe_v->n_quadrature_points;
  
  std::vector<Vector<double> > load_vector (n_q_points, Vector<double>(size) );
  
  /* Evaluate rhs and solution on the quadrature points. */
  rhs->vector_value_list (fe_v->get_quadrature_points(), load_vector);

  unsigned int comp_i;
  for (unsigned int i=0; i<fe_v->dofs_per_cell; ++i) {
    comp_i = fe->system_to_component_index(i).first;
    for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
      cell_rhs(i) +=
	( load_vector[q_point](comp_i) *
	  fe_v->shape_value(i,q_point) *
	  fe_v->JxW(q_point) );
    }
  }
}
