#include "../include/local_assemble_elastic_rhs.h"
#include <deal.II/base/quadrature_lib.h>

template <int dim>
LocalAssembleElasticRHS<dim>::LocalAssembleElasticRHS() : 
fe(0, "Local Assemble Fe Pointer"),
  fe_v(0, "Local Assemble FeValues Pointer"),
  fe_face_v(0, "Local Assemble FeFaceValues Pointer")
{
}


template <int dim>
void LocalAssembleElasticRHS<dim>::reinit (FiniteElement<dim> &myfe,
					Function<dim> &bf,
					Function<dim> &nbc,
					std::map<char, std::vector<bool> > & n_map) 
{ 
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  if(fe_face_v) {
    FEFaceValues<dim> * p = fe_face_v;
    fe_face_v = 0;
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
    
  fe_face_v = new FEFaceValues<dim>(*fe, face_quadrature, flags);
    
  this->flags = 
    assemble_rhs_cell|
    assemble_rhs_boundary;
    
  rhs = &bf;
  neumann = &nbc;
  neumann_map = n_map;
}

template <int dim>
LocalAssembleElasticRHS<dim>::~LocalAssembleElasticRHS()
{
  if(fe_v) {
    FEValues<dim> * p = fe_v;
    fe_v = 0;
    delete p;
  }
  if(fe_face_v) {
    FEFaceValues<dim> * p = fe_face_v;
    fe_face_v = 0;
    delete p;
  }
  fe = 0;
}


template <int dim>
void LocalAssembleElasticRHS<dim>::assemble_rhs_term
(const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
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

template<int dim>
void LocalAssembleElasticRHS<dim>::assemble_rhs_boundary_term
(const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
 const unsigned int face_no,
 Vector<double> &cell_rhs)
{
  cell_rhs = 0;
  // See if we need to do anything here
  char id = cell->face(face_no)->boundary_indicator();
  if(neumann_map.find(id) == neumann_map.end()) 
    return;
  
  std::vector<bool> & filter = neumann_map[id];
  
  Assert(fe, ExcNotInitialized());
  fe_face_v->reinit(cell,face_no);
  
  unsigned int size = fe->n_components();
  unsigned int n_q_points = fe_face_v->n_quadrature_points;
  
  /** Vector of boundary values.*/
  std::vector<Vector<double> > neumann_vector (n_q_points, Vector<double>(size) );
  
  /* Evaluate rhs and solution on the quadrature points. */
  neumann->vector_value_list (fe_face_v->get_quadrature_points(), neumann_vector);
  
  unsigned int comp_i;
  for (unsigned int i=0; i<fe_face_v->dofs_per_cell; ++i) {
    comp_i = fe->system_to_component_index(i).first;
    if(filter[comp_i]) 
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point) {
	cell_rhs(i) +=
	  ( neumann_vector[q_point](comp_i) *
	    fe_face_v->shape_value(i,q_point) *
	    fe_face_v->JxW(q_point) );
      }
  }
}
