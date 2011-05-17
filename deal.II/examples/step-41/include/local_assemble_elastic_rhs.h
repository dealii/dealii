#ifndef LOCAL_ASSEMBLE_ELASTIC_RHS
#define LOCAL_ASSEMBLE_ELASTIC_RHS

#include "local_assemble.h"
#include "parsed_symmetric_tensor_function.h"
#include <deal.II/base/parsed_function.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe.h>

#include <fstream>
#include <iostream>
#include <deal.II/base/parameter_handler.h>


template <int dim>
class LocalAssembleElasticRHS : public LocalAssembleBase<dim>
{
 public:
  
  LocalAssembleElasticRHS(); 
  
  ~LocalAssembleElasticRHS(); 
  
  virtual void assemble_rhs_term
    (const typename MGDoFHandler<dim>::active_cell_iterator&,
     Vector<double> &);
  
  /** This object will be called for each boundary face of a triangulation.*/
  virtual void assemble_rhs_boundary_term
    (const typename MGDoFHandler<dim>::active_cell_iterator&, const unsigned int, 
     Vector<double> &);


  void reinit(FiniteElement<dim>&,
	      Function<dim> &,
	      Function<dim> &,
	      std::map<char, std::vector<bool> > &);
			
 private:
  /** A pointer to the finite element.*/
  SmartPointer<FiniteElement<dim> > fe;

  /** A pointer to the rhs function.*/
  SmartPointer<Function<dim> > rhs;

  /** A pointer to the Neumann function.*/
  SmartPointer<Function<dim> > neumann;
  
  /** A map of ids and components for neumann bc. */
  std::map<char, std::vector<bool> > neumann_map;

  /** A pointer to fe_values objects. */
  SmartPointer<FEValues<dim> > fe_v;
  
  /** A pointer to fe_face_values objects. */
  SmartPointer<FEFaceValues<dim> > fe_face_v;
  
};

#endif
