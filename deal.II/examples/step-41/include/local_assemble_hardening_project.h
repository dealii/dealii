#ifndef LOCAL_ASSEMBLE_HARDENING_PROJECT
#define LOCAL_ASSEMBLE_HARDENING_PROJECT

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
class LocalAssemblePlasticProject : public LocalAssembleBase<dim>
{
 public:
  
  LocalAssemblePlasticProject(); 
  
  ~LocalAssemblePlasticProject(); 
  
  virtual void assemble_rhs_term
    (const typename MGDoFHandler<dim>::active_cell_iterator&,
     Vector<double> &);
  
  void parameters(ParameterHandler &prm);


  void reinit(FiniteElement<dim>&,
	      Table<2, SymmetricTensor<2,dim> > &);
			
 private:
  /** A pointer to the finite element.*/
  SmartPointer<FiniteElement<dim> > fe;

  /** A pointer to the plastic strain */
  Table<dim,SymmetricTensor<2,dim> >  plastic_strain;

  /** A pointer to fe_values objects. */
  SmartPointer<FEValues<dim> > fe_v;

  ParsedSymmetricTensorFunction<4, dim> C;
  
};

#endif
