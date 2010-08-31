#include "local_assemble_base.h"
                
template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_cell_term (const typename DH::active_cell_iterator&, 
						 FullMatrix<double> &)
{
  AssertThrow(false, ExcPureFunctionCalled());
}

template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_boundary_term
(const typename DH::active_cell_iterator&,const  unsigned int, 
 FullMatrix<double> &) 
{
  AssertThrow(false, ExcPureFunctionCalled());
}

template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator&,const unsigned int, 
 const typename DH::active_cell_iterator&,const unsigned int, 
 FullMatrix<double> &,
 FullMatrix<double> &) 
{
  AssertThrow(false, ExcPureFunctionCalled());
}

template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator&,
 const unsigned int, const unsigned int,
 const typename DH::active_cell_iterator&,const unsigned int, 
 FullMatrix<double> &,
 FullMatrix<double> &) 
{
  AssertThrow(false, ExcPureFunctionCalled());
}

template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator&,const unsigned int, 
 const typename DH::active_cell_iterator&,
 const unsigned int, 
 const unsigned int,
 FullMatrix<double> &,
 FullMatrix<double> &) 
{
  AssertThrow(false, ExcPureFunctionCalled());
}


template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_rhs_term
(const typename DH::active_cell_iterator&, 
 Vector<double> &)
{
  AssertThrow(false, ExcPureFunctionCalled());
}

template <int dim, typename DH>  
void LocalAssembleBase<dim,DH>::assemble_rhs_boundary_term
(const typename DH::active_cell_iterator&, 
 const unsigned, 
 Vector<double> &)
{
  AssertThrow(false, ExcPureFunctionCalled());
}

