#ifndef LOCAL_ASSEMBLE_ELASTIC_MATRIX
#define LOCAL_ASSEMBLE_ELASTIC_MATRIX

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
class LocalAssembleElasticMatrix : public LocalAssembleBase<dim>
{
 public:
  
  LocalAssembleElasticMatrix(); 
  
  ~LocalAssembleElasticMatrix(); 
  
  virtual void  assemble_cell_term
    (const typename MGDoFHandler<dim>::active_cell_iterator& cell, 
     FullMatrix<double> &);


  void reinit(FiniteElement<dim>&);

  void parameters(ParameterHandler &prm);
			
 private:
  /** A pointer to the finite element.*/
  SmartPointer<FiniteElement<dim> > fe;

  /** A pointer to the rhs function.*/
  ParsedSymmetricTensorFunction<4, dim>  C;

  /** A pointer to fe_values objects. */
  SmartPointer<FEValues<dim> > fe_v;
  
};

#endif
