#include "local_assemble.h"
#include <ostream>

#if deal_II_dimension != 1

template <int dim, typename DH>
LocalAssemble<dim,DH>::LocalAssemble() : 
  fe_face_v(0, "Local Assemble FeFaceValues Pointer"),
  fe_face_n_v(0, "Local Assemble FeFaceValues on Neighbor Pointer"),
  fe_sub_face_v(0, "Local Assemble FesubfaceValues Pointer")
{
}

template <int dim, typename DH>
void  LocalAssemble<dim,DH>::assemble_face_terms (FEFaceValuesBase<dim> &,
						  FEFaceValuesBase<dim> &,
						  FullMatrix<double> &,
						  FullMatrix<double> &)
{
  Assert(false, ExcPureFunctionCalled());
}

  
template <int dim, typename DH>
void LocalAssemble<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator& cell, 
 const unsigned int face_no,
 const typename DH::active_cell_iterator& n_cell, 
 const unsigned int n_face_no,
 FullMatrix<double> &cell_m,
 FullMatrix<double> &cell_n_m)
{
  this->fe_face_v->reinit(cell, face_no);
  this->fe_face_n_v->reinit(n_cell, n_face_no);
  assemble_face_terms(*fe_face_v, *fe_face_n_v, cell_m, cell_n_m);
}

template <int dim, typename DH>
void LocalAssemble<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator& cell, 
 const unsigned int face_no,
 const typename DH::active_cell_iterator& n_cell, 
 const unsigned int n_face_no,
 const unsigned int n_sub_face_no,
 FullMatrix<double> &cell_m,
 FullMatrix<double> &cell_n_m)
{
  fe_face_v->reinit(cell, face_no);
  fe_sub_face_v->reinit(n_cell, n_face_no, n_sub_face_no);
  assemble_face_terms(*fe_face_v, *fe_sub_face_v, cell_m, cell_n_m);
}

template <int dim, typename DH>
void LocalAssemble<dim,DH>::assemble_face_term
(const typename DH::active_cell_iterator& cell, 
 const unsigned int face_no,
 const unsigned int subface_no,
 const typename DH::active_cell_iterator& n_cell, 
 const unsigned int n_face_no,
 FullMatrix<double> &cell_m,
 FullMatrix<double> &cell_n_m)
{
  fe_sub_face_v->reinit(cell, face_no, subface_no);
  fe_face_v->reinit(n_cell, n_face_no);
  assemble_face_terms(*fe_sub_face_v, *fe_face_v, cell_m, cell_n_m);
}
#endif
