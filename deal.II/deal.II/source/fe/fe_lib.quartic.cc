/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <fe/fe_lib.lagrange.h>
#include <grid/tria_iterator.h>
#include <grid/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>



// declare explicit specializations before use:
template <> void FEQ4<deal_II_dimension>::initialize_matrices ();



#if deal_II_dimension == 1

template <>
FEQ4<1>::FEQ4 () :
		FEQ1Mapping<1> (1, 3) {
  initialize_matrices ();
};


template <>
FEQ4<1>::FEQ4 (const int) :
		FEQ1Mapping<1> (0, 5) {
  initialize_matrices ();
};




template <>
void FEQ4<1>::initialize_matrices () {
  prolongation[0](0,0) = 1.0;
  prolongation[0](1,3) = 1.0;
  prolongation[0](2,0) = 35.0/128.0;
  prolongation[0](2,1) = -5.0/128.0;
  prolongation[0](2,2) = 35.0/32.0;
  prolongation[0](2,3) = -35.0/64.0;
  prolongation[0](2,4) = 7.0/32.0;
  prolongation[0](3,2) = 1.0;
  prolongation[0](4,0) = -5.0/128.0;
  prolongation[0](4,1) = 3.0/128.0;
  prolongation[0](4,2) = 15.0/32.0;
  prolongation[0](4,3) = 45.0/64.0;
  prolongation[0](4,4) = -5.0/32.0;
  prolongation[1](0,3) = 1.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](2,0) = 3.0/128.0;
  prolongation[1](2,1) = -5.0/128.0;
  prolongation[1](2,2) = -5.0/32.0;
  prolongation[1](2,3) = 45.0/64.0;
  prolongation[1](2,4) = 15.0/32.0;
  prolongation[1](3,4) = 1.0;
  prolongation[1](4,0) = -5.0/128.0;
  prolongation[1](4,1) = 35.0/128.0;
  prolongation[1](4,2) = 7.0/32.0;
  prolongation[1](4,3) = -35.0/64.0;
  prolongation[1](4,4) = 35.0/32.0;

  restriction[0](0,0) = 1.0;
  restriction[0](2,3) = 1.0;
  restriction[0](3,1) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](3,0) = 1.0;
  restriction[1](4,3) = 1.0;
};



template <>
double
FEQ4<1>::shape_value(const unsigned int i,
			     const Point<1>     &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
  switch (i)
    {
      case 0: return 32.0/3.0*xi*xi*xi*xi-80.0/3.0*xi*xi*xi+70.0/3.0*xi*xi-25.0/3.0*xi+1.0;
      case 1: return 32.0/3.0*xi*xi*xi*xi-16.0*xi*xi*xi+22.0/3.0*xi*xi-xi;
      case 2: return -128.0/3.0*xi*xi*xi*xi+96.0*xi*xi*xi-208.0/3.0*xi*xi+16.0*xi;
      case 3: return 64.0*xi*xi*xi*xi-128.0*xi*xi*xi+76.0*xi*xi-12.0*xi;
      case 4: return -128.0/3.0*xi*xi*xi*xi+224.0/3.0*xi*xi*xi-112.0/3.0*xi*xi+16.0/3.0*xi;
    };
  return 0.;
};



template <>
Tensor<1,1>
FEQ4<1>::shape_grad(const unsigned int i,
			    const Point<1>    &p) const
{
  Assert((i<total_dofs), ExcInvalidIndex(i));
  const double xi = p(0);
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<1>(128.0/3.0*xi*xi*xi-80.0*xi*xi+140.0/3.0*xi-25.0/3.0);
      case 1: return Point<1>(128.0/3.0*xi*xi*xi-48.0*xi*xi+44.0/3.0*xi-1.0);
      case 2: return Point<1>(-512.0/3.0*xi*xi*xi+288.0*xi*xi-416.0/3.0*xi+16.0);
      case 3: return Point<1>(256.0*xi*xi*xi-384.0*xi*xi+152.0*xi-12.0);
      case 4: return Point<1>(-512.0/3.0*xi*xi*xi+224.0*xi*xi-224.0/3.0*xi+16.0/3.0);
    };
  return Point<1>();
};



template <>
Tensor<2,1>
FEQ4<1>::shape_grad_grad (const unsigned int i,
				  const Point<1>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0);
  Tensor<2,1> return_value;
  switch (i) 
    {
      case 0: return_value[0][0] = 128.0*xi*xi-160.0*xi+140.0/3.0;
      case 1: return_value[0][0] = 128.0*xi*xi-96.0*xi+44.0/3.0;
      case 2: return_value[0][0] = -512.0*xi*xi+576.0*xi-416.0/3.0;
      case 3: return_value[0][0] = 768.0*xi*xi-768.0*xi+152.0;
      case 4: return_value[0][0] = -512.0*xi*xi+448.0*xi-224.0/3.0;
    };

  return return_value;
};



template <>
void FEQ4<1>::get_unit_support_points (vector<Point<1> > &unit_points) const {
  FiniteElement<1>::get_unit_support_points (unit_points);
};



template <>
void FEQ4<1>::get_support_points (const DoFHandler<1>::cell_iterator &cell,
				  vector<Point<1> >  &support_points) const {
  FiniteElement<1>::get_support_points (cell, support_points);
};



template <>
void FEQ4<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
				       vector<Point<1> >  &) const {
  Assert (false, ExcInternalError());
};



template <>
void FEQ4<1>::get_local_mass_matrix (const DoFHandler<1>::cell_iterator &cell,
					     FullMatrix<double> &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double h = cell->vertex(1)(0) - cell->vertex(0)(0);
  Assert (h>0, ExcJacobiDeterminantHasWrongSign());

  const double t1 = 146.0/2835.0*h;
  const double t2 = 29.0/5670.0*h;
  const double t3 = 148.0/2835.0*h;
  const double t4 = 29.0/945.0*h;
  const double t5 = 4.0/405.0*h;
  const double t6 = 128.0/405.0*h;
  const double t7 = 64.0/945.0*h;
  const double t8 = 128.0/2835.0*h;
  local_mass_matrix(0,0) = t1;
  local_mass_matrix(0,1) = -t2;
  local_mass_matrix(0,2) = t3;
  local_mass_matrix(0,3) = -t4;
  local_mass_matrix(0,4) = t5;
  local_mass_matrix(1,0) = -t2;
  local_mass_matrix(1,1) = t1;
  local_mass_matrix(1,2) = t5;
  local_mass_matrix(1,3) = -t4;
  local_mass_matrix(1,4) = t3;
  local_mass_matrix(2,0) = t3;
  local_mass_matrix(2,1) = t5;
  local_mass_matrix(2,2) = t6;
  local_mass_matrix(2,3) = -t7;
  local_mass_matrix(2,4) = t8;
  local_mass_matrix(3,0) = -t4;
  local_mass_matrix(3,1) = -t4;
  local_mass_matrix(3,2) = -t7;
  local_mass_matrix(3,3) = 104.0/315.0*h;
  local_mass_matrix(3,4) = -t7;
  local_mass_matrix(4,0) = t5;
  local_mass_matrix(4,1) = t3;
  local_mass_matrix(4,2) = t8;
  local_mass_matrix(4,3) = -t7;
  local_mass_matrix(4,4) = t6;
};

#endif




#if deal_II_dimension == 2

template <>
FEQ4<2>::FEQ4 () :
		FEQ1Mapping<2> (1, 3, 9)
{
  interface_constraints(0,3) = 1.0;
  interface_constraints(1,0) = 35.0/128.0;
  interface_constraints(1,1) = -5.0/128.0;
  interface_constraints(1,2) = 35.0/32.0;
  interface_constraints(1,3) = -35.0/64.0;
  interface_constraints(1,4) = 7.0/32.0;
  interface_constraints(2,2) = 1.0;
  interface_constraints(3,0) = -5.0/128.0;
  interface_constraints(3,1) = 3.0/128.0;
  interface_constraints(3,2) = 15.0/32.0;
  interface_constraints(3,3) = 45.0/64.0;
  interface_constraints(3,4) = -5.0/32.0;
  interface_constraints(4,0) = 3.0/128.0;
  interface_constraints(4,1) = -5.0/128.0;
  interface_constraints(4,2) = -5.0/32.0;
  interface_constraints(4,3) = 45.0/64.0;
  interface_constraints(4,4) = 15.0/32.0;
  interface_constraints(5,4) = 1.0;
  interface_constraints(6,0) = -5.0/128.0;
  interface_constraints(6,1) = 35.0/128.0;
  interface_constraints(6,2) = 7.0/32.0;
  interface_constraints(6,3) = -35.0/64.0;
  interface_constraints(6,4) = 35.0/32.0;

  initialize_matrices ();
};



template <>
FEQ4<2>::FEQ4 (const int) :
		FEQ1Mapping<2> (0, 0, 25)
{
  initialize_matrices ();
};



template <>
void FEQ4<2>::initialize_matrices () {
  prolongation[0](0,0) = 1.0;
  prolongation[0](1,5) = 1.0;
  prolongation[0](2,24) = 1.0;
  prolongation[0](3,14) = 1.0;
  prolongation[0](4,0) = 35.0/128.0;
  prolongation[0](4,1) = -5.0/128.0;
  prolongation[0](4,4) = 35.0/32.0;
  prolongation[0](4,5) = -35.0/64.0;
  prolongation[0](4,6) = 7.0/32.0;
  prolongation[0](5,4) = 1.0;
  prolongation[0](6,0) = -5.0/128.0;
  prolongation[0](6,1) = 3.0/128.0;
  prolongation[0](6,4) = 15.0/32.0;
  prolongation[0](6,5) = 45.0/64.0;
  prolongation[0](6,6) = -5.0/32.0;
  prolongation[0](7,5) = 35.0/128.0;
  prolongation[0](7,11) = -5.0/128.0;
  prolongation[0](7,20) = 35.0/32.0;
  prolongation[0](7,22) = 7.0/32.0;
  prolongation[0](7,24) = -35.0/64.0;
  prolongation[0](8,20) = 1.0;
  prolongation[0](9,5) = -5.0/128.0;
  prolongation[0](9,11) = 3.0/128.0;
  prolongation[0](9,20) = 15.0/32.0;
  prolongation[0](9,22) = -5.0/32.0;
  prolongation[0](9,24) = 45.0/64.0;
  prolongation[0](10,8) = -5.0/128.0;
  prolongation[0](10,14) = 35.0/128.0;
  prolongation[0](10,21) = 7.0/32.0;
  prolongation[0](10,23) = 35.0/32.0;
  prolongation[0](10,24) = -35.0/64.0;
  prolongation[0](11,23) = 1.0;
  prolongation[0](12,8) = 3.0/128.0;
  prolongation[0](12,14) = -5.0/128.0;
  prolongation[0](12,21) = -5.0/32.0;
  prolongation[0](12,23) = 15.0/32.0;
  prolongation[0](12,24) = 45.0/64.0;
  prolongation[0](13,0) = 35.0/128.0;
  prolongation[0](13,3) = -5.0/128.0;
  prolongation[0](13,13) = 35.0/32.0;
  prolongation[0](13,14) = -35.0/64.0;
  prolongation[0](13,15) = 7.0/32.0;
  prolongation[0](14,13) = 1.0;
  prolongation[0](15,0) = -5.0/128.0;
  prolongation[0](15,3) = 3.0/128.0;
  prolongation[0](15,13) = 15.0/32.0;
  prolongation[0](15,14) = 45.0/64.0;
  prolongation[0](15,15) = -5.0/32.0;
  prolongation[0](16,0) = 1225.0/16384.0;
  prolongation[0](16,1) = -175.0/16384.0;
  prolongation[0](16,2) = 25.0/16384.0;
  prolongation[0](16,3) = -175.0/16384.0;
  prolongation[0](16,4) = 1225.0/4096.0;
  prolongation[0](16,5) = -1225.0/8192.0;
  prolongation[0](16,6) = 245.0/4096.0;
  prolongation[0](16,7) = -175.0/4096.0;
  prolongation[0](16,8) = 175.0/8192.0;
  prolongation[0](16,9) = -35.0/4096.0;
  prolongation[0](16,10) = -175.0/4096.0;
  prolongation[0](16,11) = 175.0/8192.0;
  prolongation[0](16,12) = -35.0/4096.0;
  prolongation[0](16,13) = 1225.0/4096.0;
  prolongation[0](16,14) = -1225.0/8192.0;
  prolongation[0](16,15) = 245.0/4096.0;
  prolongation[0](16,16) = 1225.0/1024.0;
  prolongation[0](16,17) = 245.0/1024.0;
  prolongation[0](16,18) = 49.0/1024.0;
  prolongation[0](16,19) = 245.0/1024.0;
  prolongation[0](16,20) = -1225.0/2048.0;
  prolongation[0](16,21) = -245.0/2048.0;
  prolongation[0](16,22) = -245.0/2048.0;
  prolongation[0](16,23) = -1225.0/2048.0;
  prolongation[0](16,24) = 1225.0/4096.0;
  prolongation[0](17,0) = -175.0/16384.0;
  prolongation[0](17,1) = 105.0/16384.0;
  prolongation[0](17,2) = -15.0/16384.0;
  prolongation[0](17,3) = 25.0/16384.0;
  prolongation[0](17,4) = 525.0/4096.0;
  prolongation[0](17,5) = 1575.0/8192.0;
  prolongation[0](17,6) = -175.0/4096.0;
  prolongation[0](17,7) = 105.0/4096.0;
  prolongation[0](17,8) = -105.0/8192.0;
  prolongation[0](17,9) = 21.0/4096.0;
  prolongation[0](17,10) = -75.0/4096.0;
  prolongation[0](17,11) = -225.0/8192.0;
  prolongation[0](17,12) = 25.0/4096.0;
  prolongation[0](17,13) = -175.0/4096.0;
  prolongation[0](17,14) = 175.0/8192.0;
  prolongation[0](17,15) = -35.0/4096.0;
  prolongation[0](17,16) = 525.0/1024.0;
  prolongation[0](17,17) = -175.0/1024.0;
  prolongation[0](17,18) = -35.0/1024.0;
  prolongation[0](17,19) = 105.0/1024.0;
  prolongation[0](17,20) = 1575.0/2048.0;
  prolongation[0](17,21) = 175.0/2048.0;
  prolongation[0](17,22) = 315.0/2048.0;
  prolongation[0](17,23) = -525.0/2048.0;
  prolongation[0](17,24) = -1575.0/4096.0;
  prolongation[0](18,0) = 25.0/16384.0;
  prolongation[0](18,1) = -15.0/16384.0;
  prolongation[0](18,2) = 9.0/16384.0;
  prolongation[0](18,3) = -15.0/16384.0;
  prolongation[0](18,4) = -75.0/4096.0;
  prolongation[0](18,5) = -225.0/8192.0;
  prolongation[0](18,6) = 25.0/4096.0;
  prolongation[0](18,7) = 45.0/4096.0;
  prolongation[0](18,8) = 135.0/8192.0;
  prolongation[0](18,9) = -15.0/4096.0;
  prolongation[0](18,10) = 45.0/4096.0;
  prolongation[0](18,11) = 135.0/8192.0;
  prolongation[0](18,12) = -15.0/4096.0;
  prolongation[0](18,13) = -75.0/4096.0;
  prolongation[0](18,14) = -225.0/8192.0;
  prolongation[0](18,15) = 25.0/4096.0;
  prolongation[0](18,16) = 225.0/1024.0;
  prolongation[0](18,17) = -75.0/1024.0;
  prolongation[0](18,18) = 25.0/1024.0;
  prolongation[0](18,19) = -75.0/1024.0;
  prolongation[0](18,20) = 675.0/2048.0;
  prolongation[0](18,21) = -225.0/2048.0;
  prolongation[0](18,22) = -225.0/2048.0;
  prolongation[0](18,23) = 675.0/2048.0;
  prolongation[0](18,24) = 2025.0/4096.0;
  prolongation[0](19,0) = -175.0/16384.0;
  prolongation[0](19,1) = 25.0/16384.0;
  prolongation[0](19,2) = -15.0/16384.0;
  prolongation[0](19,3) = 105.0/16384.0;
  prolongation[0](19,4) = -175.0/4096.0;
  prolongation[0](19,5) = 175.0/8192.0;
  prolongation[0](19,6) = -35.0/4096.0;
  prolongation[0](19,7) = -75.0/4096.0;
  prolongation[0](19,8) = -225.0/8192.0;
  prolongation[0](19,9) = 25.0/4096.0;
  prolongation[0](19,10) = 105.0/4096.0;
  prolongation[0](19,11) = -105.0/8192.0;
  prolongation[0](19,12) = 21.0/4096.0;
  prolongation[0](19,13) = 525.0/4096.0;
  prolongation[0](19,14) = 1575.0/8192.0;
  prolongation[0](19,15) = -175.0/4096.0;
  prolongation[0](19,16) = 525.0/1024.0;
  prolongation[0](19,17) = 105.0/1024.0;
  prolongation[0](19,18) = -35.0/1024.0;
  prolongation[0](19,19) = -175.0/1024.0;
  prolongation[0](19,20) = -525.0/2048.0;
  prolongation[0](19,21) = 315.0/2048.0;
  prolongation[0](19,22) = 175.0/2048.0;
  prolongation[0](19,23) = 1575.0/2048.0;
  prolongation[0](19,24) = -1575.0/4096.0;
  prolongation[0](20,4) = 35.0/128.0;
  prolongation[0](20,10) = -5.0/128.0;
  prolongation[0](20,16) = 35.0/32.0;
  prolongation[0](20,19) = 7.0/32.0;
  prolongation[0](20,23) = -35.0/64.0;
  prolongation[0](21,7) = 3.0/128.0;
  prolongation[0](21,13) = -5.0/128.0;
  prolongation[0](21,16) = 15.0/32.0;
  prolongation[0](21,17) = -5.0/32.0;
  prolongation[0](21,20) = 45.0/64.0;
  prolongation[0](22,4) = -5.0/128.0;
  prolongation[0](22,10) = 3.0/128.0;
  prolongation[0](22,16) = 15.0/32.0;
  prolongation[0](22,19) = -5.0/32.0;
  prolongation[0](22,23) = 45.0/64.0;
  prolongation[0](23,7) = -5.0/128.0;
  prolongation[0](23,13) = 35.0/128.0;
  prolongation[0](23,16) = 35.0/32.0;
  prolongation[0](23,17) = 7.0/32.0;
  prolongation[0](23,20) = -35.0/64.0;
  prolongation[0](24,16) = 1.0;
  prolongation[1](0,5) = 1.0;
  prolongation[1](1,1) = 1.0;
  prolongation[1](2,8) = 1.0;
  prolongation[1](3,24) = 1.0;
  prolongation[1](4,0) = 3.0/128.0;
  prolongation[1](4,1) = -5.0/128.0;
  prolongation[1](4,4) = -5.0/32.0;
  prolongation[1](4,5) = 45.0/64.0;
  prolongation[1](4,6) = 15.0/32.0;
  prolongation[1](5,6) = 1.0;
  prolongation[1](6,0) = -5.0/128.0;
  prolongation[1](6,1) = 35.0/128.0;
  prolongation[1](6,4) = 7.0/32.0;
  prolongation[1](6,5) = -35.0/64.0;
  prolongation[1](6,6) = 35.0/32.0;
  prolongation[1](7,1) = 35.0/128.0;
  prolongation[1](7,2) = -5.0/128.0;
  prolongation[1](7,7) = 35.0/32.0;
  prolongation[1](7,8) = -35.0/64.0;
  prolongation[1](7,9) = 7.0/32.0;
  prolongation[1](8,7) = 1.0;
  prolongation[1](9,1) = -5.0/128.0;
  prolongation[1](9,2) = 3.0/128.0;
  prolongation[1](9,7) = 15.0/32.0;
  prolongation[1](9,8) = 45.0/64.0;
  prolongation[1](9,9) = -5.0/32.0;
  prolongation[1](10,8) = -5.0/128.0;
  prolongation[1](10,14) = 3.0/128.0;
  prolongation[1](10,21) = 15.0/32.0;
  prolongation[1](10,23) = -5.0/32.0;
  prolongation[1](10,24) = 45.0/64.0;
  prolongation[1](11,21) = 1.0;
  prolongation[1](12,8) = 35.0/128.0;
  prolongation[1](12,14) = -5.0/128.0;
  prolongation[1](12,21) = 35.0/32.0;
  prolongation[1](12,23) = 7.0/32.0;
  prolongation[1](12,24) = -35.0/64.0;
  prolongation[1](13,5) = 35.0/128.0;
  prolongation[1](13,11) = -5.0/128.0;
  prolongation[1](13,20) = 35.0/32.0;
  prolongation[1](13,22) = 7.0/32.0;
  prolongation[1](13,24) = -35.0/64.0;
  prolongation[1](14,20) = 1.0;
  prolongation[1](15,5) = -5.0/128.0;
  prolongation[1](15,11) = 3.0/128.0;
  prolongation[1](15,20) = 15.0/32.0;
  prolongation[1](15,22) = -5.0/32.0;
  prolongation[1](15,24) = 45.0/64.0;
  prolongation[1](16,0) = 105.0/16384.0;
  prolongation[1](16,1) = -175.0/16384.0;
  prolongation[1](16,2) = 25.0/16384.0;
  prolongation[1](16,3) = -15.0/16384.0;
  prolongation[1](16,4) = -175.0/4096.0;
  prolongation[1](16,5) = 1575.0/8192.0;
  prolongation[1](16,6) = 525.0/4096.0;
  prolongation[1](16,7) = -175.0/4096.0;
  prolongation[1](16,8) = 175.0/8192.0;
  prolongation[1](16,9) = -35.0/4096.0;
  prolongation[1](16,10) = 25.0/4096.0;
  prolongation[1](16,11) = -225.0/8192.0;
  prolongation[1](16,12) = -75.0/4096.0;
  prolongation[1](16,13) = 105.0/4096.0;
  prolongation[1](16,14) = -105.0/8192.0;
  prolongation[1](16,15) = 21.0/4096.0;
  prolongation[1](16,16) = -175.0/1024.0;
  prolongation[1](16,17) = 525.0/1024.0;
  prolongation[1](16,18) = 105.0/1024.0;
  prolongation[1](16,19) = -35.0/1024.0;
  prolongation[1](16,20) = 1575.0/2048.0;
  prolongation[1](16,21) = -525.0/2048.0;
  prolongation[1](16,22) = 315.0/2048.0;
  prolongation[1](16,23) = 175.0/2048.0;
  prolongation[1](16,24) = -1575.0/4096.0;
  prolongation[1](17,0) = -175.0/16384.0;
  prolongation[1](17,1) = 1225.0/16384.0;
  prolongation[1](17,2) = -175.0/16384.0;
  prolongation[1](17,3) = 25.0/16384.0;
  prolongation[1](17,4) = 245.0/4096.0;
  prolongation[1](17,5) = -1225.0/8192.0;
  prolongation[1](17,6) = 1225.0/4096.0;
  prolongation[1](17,7) = 1225.0/4096.0;
  prolongation[1](17,8) = -1225.0/8192.0;
  prolongation[1](17,9) = 245.0/4096.0;
  prolongation[1](17,10) = -35.0/4096.0;
  prolongation[1](17,11) = 175.0/8192.0;
  prolongation[1](17,12) = -175.0/4096.0;
  prolongation[1](17,13) = -175.0/4096.0;
  prolongation[1](17,14) = 175.0/8192.0;
  prolongation[1](17,15) = -35.0/4096.0;
  prolongation[1](17,16) = 245.0/1024.0;
  prolongation[1](17,17) = 1225.0/1024.0;
  prolongation[1](17,18) = 245.0/1024.0;
  prolongation[1](17,19) = 49.0/1024.0;
  prolongation[1](17,20) = -1225.0/2048.0;
  prolongation[1](17,21) = -1225.0/2048.0;
  prolongation[1](17,22) = -245.0/2048.0;
  prolongation[1](17,23) = -245.0/2048.0;
  prolongation[1](17,24) = 1225.0/4096.0;
  prolongation[1](18,0) = 25.0/16384.0;
  prolongation[1](18,1) = -175.0/16384.0;
  prolongation[1](18,2) = 105.0/16384.0;
  prolongation[1](18,3) = -15.0/16384.0;
  prolongation[1](18,4) = -35.0/4096.0;
  prolongation[1](18,5) = 175.0/8192.0;
  prolongation[1](18,6) = -175.0/4096.0;
  prolongation[1](18,7) = 525.0/4096.0;
  prolongation[1](18,8) = 1575.0/8192.0;
  prolongation[1](18,9) = -175.0/4096.0;
  prolongation[1](18,10) = 21.0/4096.0;
  prolongation[1](18,11) = -105.0/8192.0;
  prolongation[1](18,12) = 105.0/4096.0;
  prolongation[1](18,13) = -75.0/4096.0;
  prolongation[1](18,14) = -225.0/8192.0;
  prolongation[1](18,15) = 25.0/4096.0;
  prolongation[1](18,16) = 105.0/1024.0;
  prolongation[1](18,17) = 525.0/1024.0;
  prolongation[1](18,18) = -175.0/1024.0;
  prolongation[1](18,19) = -35.0/1024.0;
  prolongation[1](18,20) = -525.0/2048.0;
  prolongation[1](18,21) = 1575.0/2048.0;
  prolongation[1](18,22) = 175.0/2048.0;
  prolongation[1](18,23) = 315.0/2048.0;
  prolongation[1](18,24) = -1575.0/4096.0;
  prolongation[1](19,0) = -15.0/16384.0;
  prolongation[1](19,1) = 25.0/16384.0;
  prolongation[1](19,2) = -15.0/16384.0;
  prolongation[1](19,3) = 9.0/16384.0;
  prolongation[1](19,4) = 25.0/4096.0;
  prolongation[1](19,5) = -225.0/8192.0;
  prolongation[1](19,6) = -75.0/4096.0;
  prolongation[1](19,7) = -75.0/4096.0;
  prolongation[1](19,8) = -225.0/8192.0;
  prolongation[1](19,9) = 25.0/4096.0;
  prolongation[1](19,10) = -15.0/4096.0;
  prolongation[1](19,11) = 135.0/8192.0;
  prolongation[1](19,12) = 45.0/4096.0;
  prolongation[1](19,13) = 45.0/4096.0;
  prolongation[1](19,14) = 135.0/8192.0;
  prolongation[1](19,15) = -15.0/4096.0;
  prolongation[1](19,16) = -75.0/1024.0;
  prolongation[1](19,17) = 225.0/1024.0;
  prolongation[1](19,18) = -75.0/1024.0;
  prolongation[1](19,19) = 25.0/1024.0;
  prolongation[1](19,20) = 675.0/2048.0;
  prolongation[1](19,21) = 675.0/2048.0;
  prolongation[1](19,22) = -225.0/2048.0;
  prolongation[1](19,23) = -225.0/2048.0;
  prolongation[1](19,24) = 2025.0/4096.0;
  prolongation[1](20,6) = 35.0/128.0;
  prolongation[1](20,12) = -5.0/128.0;
  prolongation[1](20,17) = 35.0/32.0;
  prolongation[1](20,18) = 7.0/32.0;
  prolongation[1](20,21) = -35.0/64.0;
  prolongation[1](21,7) = 35.0/128.0;
  prolongation[1](21,13) = -5.0/128.0;
  prolongation[1](21,16) = 7.0/32.0;
  prolongation[1](21,17) = 35.0/32.0;
  prolongation[1](21,20) = -35.0/64.0;
  prolongation[1](22,6) = -5.0/128.0;
  prolongation[1](22,12) = 3.0/128.0;
  prolongation[1](22,17) = 15.0/32.0;
  prolongation[1](22,18) = -5.0/32.0;
  prolongation[1](22,21) = 45.0/64.0;
  prolongation[1](23,7) = -5.0/128.0;
  prolongation[1](23,13) = 3.0/128.0;
  prolongation[1](23,16) = -5.0/32.0;
  prolongation[1](23,17) = 15.0/32.0;
  prolongation[1](23,20) = 45.0/64.0;
  prolongation[1](24,17) = 1.0;
  prolongation[2](0,24) = 1.0;
  prolongation[2](1,8) = 1.0;
  prolongation[2](2,2) = 1.0;
  prolongation[2](3,11) = 1.0;
  prolongation[2](4,8) = -5.0/128.0;
  prolongation[2](4,14) = 3.0/128.0;
  prolongation[2](4,21) = 15.0/32.0;
  prolongation[2](4,23) = -5.0/32.0;
  prolongation[2](4,24) = 45.0/64.0;
  prolongation[2](5,21) = 1.0;
  prolongation[2](6,8) = 35.0/128.0;
  prolongation[2](6,14) = -5.0/128.0;
  prolongation[2](6,21) = 35.0/32.0;
  prolongation[2](6,23) = 7.0/32.0;
  prolongation[2](6,24) = -35.0/64.0;
  prolongation[2](7,1) = 3.0/128.0;
  prolongation[2](7,2) = -5.0/128.0;
  prolongation[2](7,7) = -5.0/32.0;
  prolongation[2](7,8) = 45.0/64.0;
  prolongation[2](7,9) = 15.0/32.0;
  prolongation[2](8,9) = 1.0;
  prolongation[2](9,1) = -5.0/128.0;
  prolongation[2](9,2) = 35.0/128.0;
  prolongation[2](9,7) = 7.0/32.0;
  prolongation[2](9,8) = -35.0/64.0;
  prolongation[2](9,9) = 35.0/32.0;
  prolongation[2](10,2) = -5.0/128.0;
  prolongation[2](10,3) = 3.0/128.0;
  prolongation[2](10,10) = -5.0/32.0;
  prolongation[2](10,11) = 45.0/64.0;
  prolongation[2](10,12) = 15.0/32.0;
  prolongation[2](11,12) = 1.0;
  prolongation[2](12,2) = 35.0/128.0;
  prolongation[2](12,3) = -5.0/128.0;
  prolongation[2](12,10) = 7.0/32.0;
  prolongation[2](12,11) = -35.0/64.0;
  prolongation[2](12,12) = 35.0/32.0;
  prolongation[2](13,5) = 3.0/128.0;
  prolongation[2](13,11) = -5.0/128.0;
  prolongation[2](13,20) = -5.0/32.0;
  prolongation[2](13,22) = 15.0/32.0;
  prolongation[2](13,24) = 45.0/64.0;
  prolongation[2](14,22) = 1.0;
  prolongation[2](15,5) = -5.0/128.0;
  prolongation[2](15,11) = 35.0/128.0;
  prolongation[2](15,20) = 7.0/32.0;
  prolongation[2](15,22) = 35.0/32.0;
  prolongation[2](15,24) = -35.0/64.0;
  prolongation[2](16,0) = 9.0/16384.0;
  prolongation[2](16,1) = -15.0/16384.0;
  prolongation[2](16,2) = 25.0/16384.0;
  prolongation[2](16,3) = -15.0/16384.0;
  prolongation[2](16,4) = -15.0/4096.0;
  prolongation[2](16,5) = 135.0/8192.0;
  prolongation[2](16,6) = 45.0/4096.0;
  prolongation[2](16,7) = 25.0/4096.0;
  prolongation[2](16,8) = -225.0/8192.0;
  prolongation[2](16,9) = -75.0/4096.0;
  prolongation[2](16,10) = 25.0/4096.0;
  prolongation[2](16,11) = -225.0/8192.0;
  prolongation[2](16,12) = -75.0/4096.0;
  prolongation[2](16,13) = -15.0/4096.0;
  prolongation[2](16,14) = 135.0/8192.0;
  prolongation[2](16,15) = 45.0/4096.0;
  prolongation[2](16,16) = 25.0/1024.0;
  prolongation[2](16,17) = -75.0/1024.0;
  prolongation[2](16,18) = 225.0/1024.0;
  prolongation[2](16,19) = -75.0/1024.0;
  prolongation[2](16,20) = -225.0/2048.0;
  prolongation[2](16,21) = 675.0/2048.0;
  prolongation[2](16,22) = 675.0/2048.0;
  prolongation[2](16,23) = -225.0/2048.0;
  prolongation[2](16,24) = 2025.0/4096.0;
  prolongation[2](17,0) = -15.0/16384.0;
  prolongation[2](17,1) = 105.0/16384.0;
  prolongation[2](17,2) = -175.0/16384.0;
  prolongation[2](17,3) = 25.0/16384.0;
  prolongation[2](17,4) = 21.0/4096.0;
  prolongation[2](17,5) = -105.0/8192.0;
  prolongation[2](17,6) = 105.0/4096.0;
  prolongation[2](17,7) = -175.0/4096.0;
  prolongation[2](17,8) = 1575.0/8192.0;
  prolongation[2](17,9) = 525.0/4096.0;
  prolongation[2](17,10) = -35.0/4096.0;
  prolongation[2](17,11) = 175.0/8192.0;
  prolongation[2](17,12) = -175.0/4096.0;
  prolongation[2](17,13) = 25.0/4096.0;
  prolongation[2](17,14) = -225.0/8192.0;
  prolongation[2](17,15) = -75.0/4096.0;
  prolongation[2](17,16) = -35.0/1024.0;
  prolongation[2](17,17) = -175.0/1024.0;
  prolongation[2](17,18) = 525.0/1024.0;
  prolongation[2](17,19) = 105.0/1024.0;
  prolongation[2](17,20) = 175.0/2048.0;
  prolongation[2](17,21) = 1575.0/2048.0;
  prolongation[2](17,22) = -525.0/2048.0;
  prolongation[2](17,23) = 315.0/2048.0;
  prolongation[2](17,24) = -1575.0/4096.0;
  prolongation[2](18,0) = 25.0/16384.0;
  prolongation[2](18,1) = -175.0/16384.0;
  prolongation[2](18,2) = 1225.0/16384.0;
  prolongation[2](18,3) = -175.0/16384.0;
  prolongation[2](18,4) = -35.0/4096.0;
  prolongation[2](18,5) = 175.0/8192.0;
  prolongation[2](18,6) = -175.0/4096.0;
  prolongation[2](18,7) = 245.0/4096.0;
  prolongation[2](18,8) = -1225.0/8192.0;
  prolongation[2](18,9) = 1225.0/4096.0;
  prolongation[2](18,10) = 245.0/4096.0;
  prolongation[2](18,11) = -1225.0/8192.0;
  prolongation[2](18,12) = 1225.0/4096.0;
  prolongation[2](18,13) = -35.0/4096.0;
  prolongation[2](18,14) = 175.0/8192.0;
  prolongation[2](18,15) = -175.0/4096.0;
  prolongation[2](18,16) = 49.0/1024.0;
  prolongation[2](18,17) = 245.0/1024.0;
  prolongation[2](18,18) = 1225.0/1024.0;
  prolongation[2](18,19) = 245.0/1024.0;
  prolongation[2](18,20) = -245.0/2048.0;
  prolongation[2](18,21) = -1225.0/2048.0;
  prolongation[2](18,22) = -1225.0/2048.0;
  prolongation[2](18,23) = -245.0/2048.0;
  prolongation[2](18,24) = 1225.0/4096.0;
  prolongation[2](19,0) = -15.0/16384.0;
  prolongation[2](19,1) = 25.0/16384.0;
  prolongation[2](19,2) = -175.0/16384.0;
  prolongation[2](19,3) = 105.0/16384.0;
  prolongation[2](19,4) = 25.0/4096.0;
  prolongation[2](19,5) = -225.0/8192.0;
  prolongation[2](19,6) = -75.0/4096.0;
  prolongation[2](19,7) = -35.0/4096.0;
  prolongation[2](19,8) = 175.0/8192.0;
  prolongation[2](19,9) = -175.0/4096.0;
  prolongation[2](19,10) = -175.0/4096.0;
  prolongation[2](19,11) = 1575.0/8192.0;
  prolongation[2](19,12) = 525.0/4096.0;
  prolongation[2](19,13) = 21.0/4096.0;
  prolongation[2](19,14) = -105.0/8192.0;
  prolongation[2](19,15) = 105.0/4096.0;
  prolongation[2](19,16) = -35.0/1024.0;
  prolongation[2](19,17) = 105.0/1024.0;
  prolongation[2](19,18) = 525.0/1024.0;
  prolongation[2](19,19) = -175.0/1024.0;
  prolongation[2](19,20) = 315.0/2048.0;
  prolongation[2](19,21) = -525.0/2048.0;
  prolongation[2](19,22) = 1575.0/2048.0;
  prolongation[2](19,23) = 175.0/2048.0;
  prolongation[2](19,24) = -1575.0/4096.0;
  prolongation[2](20,6) = 3.0/128.0;
  prolongation[2](20,12) = -5.0/128.0;
  prolongation[2](20,17) = -5.0/32.0;
  prolongation[2](20,18) = 15.0/32.0;
  prolongation[2](20,21) = 45.0/64.0;
  prolongation[2](21,9) = 35.0/128.0;
  prolongation[2](21,15) = -5.0/128.0;
  prolongation[2](21,18) = 35.0/32.0;
  prolongation[2](21,19) = 7.0/32.0;
  prolongation[2](21,22) = -35.0/64.0;
  prolongation[2](22,6) = -5.0/128.0;
  prolongation[2](22,12) = 35.0/128.0;
  prolongation[2](22,17) = 7.0/32.0;
  prolongation[2](22,18) = 35.0/32.0;
  prolongation[2](22,21) = -35.0/64.0;
  prolongation[2](23,9) = -5.0/128.0;
  prolongation[2](23,15) = 3.0/128.0;
  prolongation[2](23,18) = 15.0/32.0;
  prolongation[2](23,19) = -5.0/32.0;
  prolongation[2](23,22) = 45.0/64.0;
  prolongation[2](24,18) = 1.0;
  prolongation[3](0,14) = 1.0;
  prolongation[3](1,24) = 1.0;
  prolongation[3](2,11) = 1.0;
  prolongation[3](3,3) = 1.0;
  prolongation[3](4,8) = -5.0/128.0;
  prolongation[3](4,14) = 35.0/128.0;
  prolongation[3](4,21) = 7.0/32.0;
  prolongation[3](4,23) = 35.0/32.0;
  prolongation[3](4,24) = -35.0/64.0;
  prolongation[3](5,23) = 1.0;
  prolongation[3](6,8) = 3.0/128.0;
  prolongation[3](6,14) = -5.0/128.0;
  prolongation[3](6,21) = -5.0/32.0;
  prolongation[3](6,23) = 15.0/32.0;
  prolongation[3](6,24) = 45.0/64.0;
  prolongation[3](7,5) = 3.0/128.0;
  prolongation[3](7,11) = -5.0/128.0;
  prolongation[3](7,20) = -5.0/32.0;
  prolongation[3](7,22) = 15.0/32.0;
  prolongation[3](7,24) = 45.0/64.0;
  prolongation[3](8,22) = 1.0;
  prolongation[3](9,5) = -5.0/128.0;
  prolongation[3](9,11) = 35.0/128.0;
  prolongation[3](9,20) = 7.0/32.0;
  prolongation[3](9,22) = 35.0/32.0;
  prolongation[3](9,24) = -35.0/64.0;
  prolongation[3](10,2) = -5.0/128.0;
  prolongation[3](10,3) = 35.0/128.0;
  prolongation[3](10,10) = 35.0/32.0;
  prolongation[3](10,11) = -35.0/64.0;
  prolongation[3](10,12) = 7.0/32.0;
  prolongation[3](11,10) = 1.0;
  prolongation[3](12,2) = 3.0/128.0;
  prolongation[3](12,3) = -5.0/128.0;
  prolongation[3](12,10) = 15.0/32.0;
  prolongation[3](12,11) = 45.0/64.0;
  prolongation[3](12,12) = -5.0/32.0;
  prolongation[3](13,0) = 3.0/128.0;
  prolongation[3](13,3) = -5.0/128.0;
  prolongation[3](13,13) = -5.0/32.0;
  prolongation[3](13,14) = 45.0/64.0;
  prolongation[3](13,15) = 15.0/32.0;
  prolongation[3](14,15) = 1.0;
  prolongation[3](15,0) = -5.0/128.0;
  prolongation[3](15,3) = 35.0/128.0;
  prolongation[3](15,13) = 7.0/32.0;
  prolongation[3](15,14) = -35.0/64.0;
  prolongation[3](15,15) = 35.0/32.0;
  prolongation[3](16,0) = 105.0/16384.0;
  prolongation[3](16,1) = -15.0/16384.0;
  prolongation[3](16,2) = 25.0/16384.0;
  prolongation[3](16,3) = -175.0/16384.0;
  prolongation[3](16,4) = 105.0/4096.0;
  prolongation[3](16,5) = -105.0/8192.0;
  prolongation[3](16,6) = 21.0/4096.0;
  prolongation[3](16,7) = 25.0/4096.0;
  prolongation[3](16,8) = -225.0/8192.0;
  prolongation[3](16,9) = -75.0/4096.0;
  prolongation[3](16,10) = -175.0/4096.0;
  prolongation[3](16,11) = 175.0/8192.0;
  prolongation[3](16,12) = -35.0/4096.0;
  prolongation[3](16,13) = -175.0/4096.0;
  prolongation[3](16,14) = 1575.0/8192.0;
  prolongation[3](16,15) = 525.0/4096.0;
  prolongation[3](16,16) = -175.0/1024.0;
  prolongation[3](16,17) = -35.0/1024.0;
  prolongation[3](16,18) = 105.0/1024.0;
  prolongation[3](16,19) = 525.0/1024.0;
  prolongation[3](16,20) = 175.0/2048.0;
  prolongation[3](16,21) = 315.0/2048.0;
  prolongation[3](16,22) = -525.0/2048.0;
  prolongation[3](16,23) = 1575.0/2048.0;
  prolongation[3](16,24) = -1575.0/4096.0;
  prolongation[3](17,0) = -15.0/16384.0;
  prolongation[3](17,1) = 9.0/16384.0;
  prolongation[3](17,2) = -15.0/16384.0;
  prolongation[3](17,3) = 25.0/16384.0;
  prolongation[3](17,4) = 45.0/4096.0;
  prolongation[3](17,5) = 135.0/8192.0;
  prolongation[3](17,6) = -15.0/4096.0;
  prolongation[3](17,7) = -15.0/4096.0;
  prolongation[3](17,8) = 135.0/8192.0;
  prolongation[3](17,9) = 45.0/4096.0;
  prolongation[3](17,10) = -75.0/4096.0;
  prolongation[3](17,11) = -225.0/8192.0;
  prolongation[3](17,12) = 25.0/4096.0;
  prolongation[3](17,13) = 25.0/4096.0;
  prolongation[3](17,14) = -225.0/8192.0;
  prolongation[3](17,15) = -75.0/4096.0;
  prolongation[3](17,16) = -75.0/1024.0;
  prolongation[3](17,17) = 25.0/1024.0;
  prolongation[3](17,18) = -75.0/1024.0;
  prolongation[3](17,19) = 225.0/1024.0;
  prolongation[3](17,20) = -225.0/2048.0;
  prolongation[3](17,21) = -225.0/2048.0;
  prolongation[3](17,22) = 675.0/2048.0;
  prolongation[3](17,23) = 675.0/2048.0;
  prolongation[3](17,24) = 2025.0/4096.0;
  prolongation[3](18,0) = 25.0/16384.0;
  prolongation[3](18,1) = -15.0/16384.0;
  prolongation[3](18,2) = 105.0/16384.0;
  prolongation[3](18,3) = -175.0/16384.0;
  prolongation[3](18,4) = -75.0/4096.0;
  prolongation[3](18,5) = -225.0/8192.0;
  prolongation[3](18,6) = 25.0/4096.0;
  prolongation[3](18,7) = 21.0/4096.0;
  prolongation[3](18,8) = -105.0/8192.0;
  prolongation[3](18,9) = 105.0/4096.0;
  prolongation[3](18,10) = 525.0/4096.0;
  prolongation[3](18,11) = 1575.0/8192.0;
  prolongation[3](18,12) = -175.0/4096.0;
  prolongation[3](18,13) = -35.0/4096.0;
  prolongation[3](18,14) = 175.0/8192.0;
  prolongation[3](18,15) = -175.0/4096.0;
  prolongation[3](18,16) = 105.0/1024.0;
  prolongation[3](18,17) = -35.0/1024.0;
  prolongation[3](18,18) = -175.0/1024.0;
  prolongation[3](18,19) = 525.0/1024.0;
  prolongation[3](18,20) = 315.0/2048.0;
  prolongation[3](18,21) = 175.0/2048.0;
  prolongation[3](18,22) = 1575.0/2048.0;
  prolongation[3](18,23) = -525.0/2048.0;
  prolongation[3](18,24) = -1575.0/4096.0;
  prolongation[3](19,0) = -175.0/16384.0;
  prolongation[3](19,1) = 25.0/16384.0;
  prolongation[3](19,2) = -175.0/16384.0;
  prolongation[3](19,3) = 1225.0/16384.0;
  prolongation[3](19,4) = -175.0/4096.0;
  prolongation[3](19,5) = 175.0/8192.0;
  prolongation[3](19,6) = -35.0/4096.0;
  prolongation[3](19,7) = -35.0/4096.0;
  prolongation[3](19,8) = 175.0/8192.0;
  prolongation[3](19,9) = -175.0/4096.0;
  prolongation[3](19,10) = 1225.0/4096.0;
  prolongation[3](19,11) = -1225.0/8192.0;
  prolongation[3](19,12) = 245.0/4096.0;
  prolongation[3](19,13) = 245.0/4096.0;
  prolongation[3](19,14) = -1225.0/8192.0;
  prolongation[3](19,15) = 1225.0/4096.0;
  prolongation[3](19,16) = 245.0/1024.0;
  prolongation[3](19,17) = 49.0/1024.0;
  prolongation[3](19,18) = 245.0/1024.0;
  prolongation[3](19,19) = 1225.0/1024.0;
  prolongation[3](19,20) = -245.0/2048.0;
  prolongation[3](19,21) = -245.0/2048.0;
  prolongation[3](19,22) = -1225.0/2048.0;
  prolongation[3](19,23) = -1225.0/2048.0;
  prolongation[3](19,24) = 1225.0/4096.0;
  prolongation[3](20,4) = 3.0/128.0;
  prolongation[3](20,10) = -5.0/128.0;
  prolongation[3](20,16) = -5.0/32.0;
  prolongation[3](20,19) = 15.0/32.0;
  prolongation[3](20,23) = 45.0/64.0;
  prolongation[3](21,9) = 3.0/128.0;
  prolongation[3](21,15) = -5.0/128.0;
  prolongation[3](21,18) = -5.0/32.0;
  prolongation[3](21,19) = 15.0/32.0;
  prolongation[3](21,22) = 45.0/64.0;
  prolongation[3](22,4) = -5.0/128.0;
  prolongation[3](22,10) = 35.0/128.0;
  prolongation[3](22,16) = 7.0/32.0;
  prolongation[3](22,19) = 35.0/32.0;
  prolongation[3](22,23) = -35.0/64.0;
  prolongation[3](23,9) = -5.0/128.0;
  prolongation[3](23,15) = 35.0/128.0;
  prolongation[3](23,18) = 7.0/32.0;
  prolongation[3](23,19) = 35.0/32.0;
  prolongation[3](23,22) = -35.0/64.0;
  prolongation[3](24,19) = 1.0;

  restriction[0](0,0) = 1.0;
  restriction[0](4,5) = 1.0;
  restriction[0](5,1) = 1.0;
  restriction[0](13,14) = 1.0;
  restriction[0](14,3) = 1.0;
  restriction[0](16,24) = 1.0;
  restriction[0](20,8) = 1.0;
  restriction[0](23,11) = 1.0;
  restriction[0](24,2) = 1.0;
  restriction[1](1,1) = 1.0;
  restriction[1](5,0) = 1.0;
  restriction[1](6,5) = 1.0;
  restriction[1](7,8) = 1.0;
  restriction[1](8,2) = 1.0;
  restriction[1](17,24) = 1.0;
  restriction[1](20,14) = 1.0;
  restriction[1](21,11) = 1.0;
  restriction[1](24,3) = 1.0;
  restriction[2](2,2) = 1.0;
  restriction[2](8,1) = 1.0;
  restriction[2](9,8) = 1.0;
  restriction[2](11,3) = 1.0;
  restriction[2](12,11) = 1.0;
  restriction[2](18,24) = 1.0;
  restriction[2](21,5) = 1.0;
  restriction[2](22,14) = 1.0;
  restriction[2](24,0) = 1.0;
  restriction[3](3,3) = 1.0;
  restriction[3](10,11) = 1.0;
  restriction[3](11,2) = 1.0;
  restriction[3](14,0) = 1.0;
  restriction[3](15,14) = 1.0;
  restriction[3](19,24) = 1.0;
  restriction[3](22,8) = 1.0;
  restriction[3](23,5) = 1.0;
  restriction[3](24,1) = 1.0;
};


template <>
double
FEQ4<2>::shape_value (const unsigned int i,
			      const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  switch (i)
    {
      case 0: return 1.0-25.0/3.0*xi+70.0/3.0*xi*xi-80.0/3.0*xi*xi*xi+32.0/
		3.0*xi*xi*xi*xi+(-25.0/3.0+625.0/9.0*xi-1750.0/9.0*xi*xi+2000.0/9.0*xi*xi*xi
				 -800.0/9.0*xi*xi*xi*xi)*eta+(70.0/3.0-1750.0/9.0*xi+4900.0/9.0*xi*xi-5600.0/9.0
							      *xi*xi*xi+2240.0/9.0*xi*xi*xi*xi)*eta*eta+(-80.0/3.0+2000.0/9.0*xi-5600.0/9.0*
													 xi*xi+6400.0/9.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(32.0/3.0-800.0/
																					9.0*xi+2240.0/9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta
		*eta;
      case 1: return -xi+22.0/3.0*xi*xi-16.0*xi*xi*xi+32.0/3.0*xi*xi*xi*xi+(
	25.0/3.0*xi-550.0/9.0*xi*xi+400.0/3.0*xi*xi*xi-800.0/9.0*xi*xi*xi*xi)*eta+(
	  -70.0/3.0*xi+1540.0/9.0*xi*xi-1120.0/3.0*xi*xi*xi+2240.0/9.0*xi*xi*xi*xi)*eta*
		eta+(80.0/3.0*xi-1760.0/9.0*xi*xi+1280.0/3.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*
		eta*eta*eta+(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*
			     xi*xi)*eta*eta*eta*eta;
      case 2: return (xi-22.0/3.0*xi*xi+16.0*xi*xi*xi-32.0/3.0*xi*xi*xi*xi)*
		eta+(-22.0/3.0*xi+484.0/9.0*xi*xi-352.0/3.0*xi*xi*xi+704.0/9.0*xi*xi*xi*xi)*eta
		*eta+(16.0*xi-352.0/3.0*xi*xi+256.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta*eta*eta
		+(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*
		eta*eta*eta;
      case 3: return (-1.0+25.0/3.0*xi-70.0/3.0*xi*xi+80.0/3.0*xi*xi*xi-32.0/
		      3.0*xi*xi*xi*xi)*eta+(22.0/3.0-550.0/9.0*xi+1540.0/9.0*xi*xi-1760.0/9.0*xi*xi*
					    xi+704.0/9.0*xi*xi*xi*xi)*eta*eta+(-16.0+400.0/3.0*xi-1120.0/3.0*xi*xi+1280.0/
									       3.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta*eta*eta+(32.0/3.0-800.0/9.0*xi+2240.0/
																9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 4: return 16.0*xi-208.0/3.0*xi*xi+96.0*xi*xi*xi-128.0/3.0*xi*xi*xi
		*xi+(-400.0/3.0*xi+5200.0/9.0*xi*xi-800.0*xi*xi*xi+3200.0/9.0*xi*xi*xi*xi)*eta+
		(1120.0/3.0*xi-14560.0/9.0*xi*xi+2240.0*xi*xi*xi-8960.0/9.0*xi*xi*xi*xi)*eta*
		eta+(-1280.0/3.0*xi+16640.0/9.0*xi*xi-2560.0*xi*xi*xi+10240.0/9.0*xi*xi*xi*xi)*
		eta*eta*eta+(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*
			     xi)*eta*eta*eta*eta;
      case 5: return -12.0*xi+76.0*xi*xi-128.0*xi*xi*xi+64.0*xi*xi*xi*xi+(
	100.0*xi-1900.0/3.0*xi*xi+3200.0/3.0*xi*xi*xi-1600.0/3.0*xi*xi*xi*xi)*eta+(
	  -280.0*xi+5320.0/3.0*xi*xi-8960.0/3.0*xi*xi*xi+4480.0/3.0*xi*xi*xi*xi)*eta*eta+
		(320.0*xi-6080.0/3.0*xi*xi+10240.0/3.0*xi*xi*xi-5120.0/3.0*xi*xi*xi*xi)*eta*eta
		*eta+(-128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*
		eta*eta*eta*eta;
      case 6: return 16.0/3.0*xi-112.0/3.0*xi*xi+224.0/3.0*xi*xi*xi-128.0/3.0
		*xi*xi*xi*xi+(-400.0/9.0*xi+2800.0/9.0*xi*xi-5600.0/9.0*xi*xi*xi+3200.0/9.0*xi*
			      xi*xi*xi)*eta+(1120.0/9.0*xi-7840.0/9.0*xi*xi+15680.0/9.0*xi*xi*xi-8960.0/9.0*
					     xi*xi*xi*xi)*eta*eta+(-1280.0/9.0*xi+8960.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+
								   10240.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*
													 xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 7: return (-16.0*xi+352.0/3.0*xi*xi-256.0*xi*xi*xi+512.0/3.0*xi*xi
		      *xi*xi)*eta+(208.0/3.0*xi-4576.0/9.0*xi*xi+3328.0/3.0*xi*xi*xi-6656.0/9.0*xi*xi
				   *xi*xi)*eta*eta+(-96.0*xi+704.0*xi*xi-1536.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta*
		eta*eta+(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*
			 xi)*eta*eta*eta*eta;
      case 8: return (12.0*xi-88.0*xi*xi+192.0*xi*xi*xi-128.0*xi*xi*xi*xi)*
		eta+(-76.0*xi+1672.0/3.0*xi*xi-1216.0*xi*xi*xi+2432.0/3.0*xi*xi*xi*xi)*eta*eta+
		(128.0*xi-2816.0/3.0*xi*xi+2048.0*xi*xi*xi-4096.0/3.0*xi*xi*xi*xi)*eta*eta*eta+
		(-64.0*xi+1408.0/3.0*xi*xi-1024.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta*
		eta;
      case 9: return (-16.0/3.0*xi+352.0/9.0*xi*xi-256.0/3.0*xi*xi*xi+512.0/
		      9.0*xi*xi*xi*xi)*eta+(112.0/3.0*xi-2464.0/9.0*xi*xi+1792.0/3.0*xi*xi*xi-3584.0/
					    9.0*xi*xi*xi*xi)*eta*eta+(-224.0/3.0*xi+4928.0/9.0*xi*xi-3584.0/3.0*xi*xi*xi+
								      7168.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*
													   xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 10: return (-16.0*xi+208.0/3.0*xi*xi-96.0*xi*xi*xi+128.0/3.0*xi*xi
		       *xi*xi)*eta+(352.0/3.0*xi-4576.0/9.0*xi*xi+704.0*xi*xi*xi-2816.0/9.0*xi*xi*xi*
				    xi)*eta*eta+(-256.0*xi+3328.0/3.0*xi*xi-1536.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)
		 *eta*eta*eta+(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi
			       *xi)*eta*eta*eta*eta;
      case 11: return (12.0*xi-76.0*xi*xi+128.0*xi*xi*xi-64.0*xi*xi*xi*xi)*
		 eta+(-88.0*xi+1672.0/3.0*xi*xi-2816.0/3.0*xi*xi*xi+1408.0/3.0*xi*xi*xi*xi)*eta*
		 eta+(192.0*xi-1216.0*xi*xi+2048.0*xi*xi*xi-1024.0*xi*xi*xi*xi)*eta*eta*eta+(
		   -128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*
		 eta*eta;
      case 12: return (-16.0/3.0*xi+112.0/3.0*xi*xi-224.0/3.0*xi*xi*xi+128.0/
		       3.0*xi*xi*xi*xi)*eta+(352.0/9.0*xi-2464.0/9.0*xi*xi+4928.0/9.0*xi*xi*xi-2816.0/
					     9.0*xi*xi*xi*xi)*eta*eta+(-256.0/3.0*xi+1792.0/3.0*xi*xi-3584.0/3.0*xi*xi*xi+
								       2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta+(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*
													    xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 13: return (16.0-400.0/3.0*xi+1120.0/3.0*xi*xi-1280.0/3.0*xi*xi*xi
		       +512.0/3.0*xi*xi*xi*xi)*eta+(-208.0/3.0+5200.0/9.0*xi-14560.0/9.0*xi*xi+16640.0
						    /9.0*xi*xi*xi-6656.0/9.0*xi*xi*xi*xi)*eta*eta+(96.0-800.0*xi+2240.0*xi*xi
												   -2560.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta*eta*eta+(-128.0/3.0+3200.0/9.0*xi
																		     -8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 14: return (-12.0+100.0*xi-280.0*xi*xi+320.0*xi*xi*xi-128.0*xi*xi*
		       xi*xi)*eta+(76.0-1900.0/3.0*xi+5320.0/3.0*xi*xi-6080.0/3.0*xi*xi*xi+2432.0/3.0*
				   xi*xi*xi*xi)*eta*eta+(-128.0+3200.0/3.0*xi-8960.0/3.0*xi*xi+10240.0/3.0*xi*xi*
							 xi-4096.0/3.0*xi*xi*xi*xi)*eta*eta*eta+(64.0-1600.0/3.0*xi+4480.0/3.0*xi*xi
												 -5120.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 15: return (16.0/3.0-400.0/9.0*xi+1120.0/9.0*xi*xi-1280.0/9.0*xi*
		       xi*xi+512.0/9.0*xi*xi*xi*xi)*eta+(-112.0/3.0+2800.0/9.0*xi-7840.0/9.0*xi*xi+
							 8960.0/9.0*xi*xi*xi-3584.0/9.0*xi*xi*xi*xi)*eta*eta+(224.0/3.0-5600.0/9.0*xi+
													      15680.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+7168.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(
														-128.0/3.0+3200.0/9.0*xi-8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi
														*xi*xi)*eta*eta*eta*eta;
      case 16: return (256.0*xi-3328.0/3.0*xi*xi+1536.0*xi*xi*xi-2048.0/3.0*
		       xi*xi*xi*xi)*eta+(-3328.0/3.0*xi+43264.0/9.0*xi*xi-6656.0*xi*xi*xi+26624.0/9.0*
					 xi*xi*xi*xi)*eta*eta+(1536.0*xi-6656.0*xi*xi+9216.0*xi*xi*xi-4096.0*xi*xi*xi*xi
					 )*eta*eta*eta+(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*xi*xi*xi+16384.0/9.0*xi*
							xi*xi*xi)*eta*eta*eta*eta;
      case 17: return (256.0/3.0*xi-1792.0/3.0*xi*xi+3584.0/3.0*xi*xi*xi
		       -2048.0/3.0*xi*xi*xi*xi)*eta+(-3328.0/9.0*xi+23296.0/9.0*xi*xi-46592.0/9.0*xi*
						     xi*xi+26624.0/9.0*xi*xi*xi*xi)*eta*eta+(512.0*xi-3584.0*xi*xi+7168.0*xi*xi*xi
											     -4096.0*xi*xi*xi*xi)*eta*eta*eta+(-2048.0/9.0*xi+14336.0/9.0*xi*xi-28672.0/9.0*
															       xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 18: return (256.0/9.0*xi-1792.0/9.0*xi*xi+3584.0/9.0*xi*xi*xi
		       -2048.0/9.0*xi*xi*xi*xi)*eta+(-1792.0/9.0*xi+12544.0/9.0*xi*xi-25088.0/9.0*xi*
						     xi*xi+14336.0/9.0*xi*xi*xi*xi)*eta*eta+(3584.0/9.0*xi-25088.0/9.0*xi*xi+50176.0
											     /9.0*xi*xi*xi-28672.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(-2048.0/9.0*xi+14336.0/9.0*
																		 xi*xi-28672.0/9.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 19: return (256.0/3.0*xi-3328.0/9.0*xi*xi+512.0*xi*xi*xi-2048.0/
		       9.0*xi*xi*xi*xi)*eta+(-1792.0/3.0*xi+23296.0/9.0*xi*xi-3584.0*xi*xi*xi+14336.0/
					     9.0*xi*xi*xi*xi)*eta*eta+(3584.0/3.0*xi-46592.0/9.0*xi*xi+7168.0*xi*xi*xi
								       -28672.0/9.0*xi*xi*xi*xi)*eta*eta*eta+(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*
													      xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 20: return (-192.0*xi+1216.0*xi*xi-2048.0*xi*xi*xi+1024.0*xi*xi*xi
		       *xi)*eta+(832.0*xi-15808.0/3.0*xi*xi+26624.0/3.0*xi*xi*xi-13312.0/3.0*xi*xi*xi*
				 xi)*eta*eta+(-1152.0*xi+7296.0*xi*xi-12288.0*xi*xi*xi+6144.0*xi*xi*xi*xi)*eta*
		 eta*eta+(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)
		 *eta*eta*eta*eta;
      case 21: return (-64.0*xi+448.0*xi*xi-896.0*xi*xi*xi+512.0*xi*xi*xi*xi)
		 *eta+(1216.0/3.0*xi-8512.0/3.0*xi*xi+17024.0/3.0*xi*xi*xi-9728.0/3.0*xi*xi*xi*
		       xi)*eta*eta+(-2048.0/3.0*xi+14336.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+16384.0/3.0*
				    xi*xi*xi*xi)*eta*eta*eta+(1024.0/3.0*xi-7168.0/3.0*xi*xi+14336.0/3.0*xi*xi*xi
							      -8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 22: return (-64.0*xi+1216.0/3.0*xi*xi-2048.0/3.0*xi*xi*xi+1024.0/
		       3.0*xi*xi*xi*xi)*eta+(448.0*xi-8512.0/3.0*xi*xi+14336.0/3.0*xi*xi*xi-7168.0/3.0
					     *xi*xi*xi*xi)*eta*eta+(-896.0*xi+17024.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+14336.0
								    /3.0*xi*xi*xi*xi)*eta*eta*eta+(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi
												   -8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta*eta;
      case 23: return (-192.0*xi+832.0*xi*xi-1152.0*xi*xi*xi+512.0*xi*xi*xi*
		       xi)*eta+(1216.0*xi-15808.0/3.0*xi*xi+7296.0*xi*xi*xi-9728.0/3.0*xi*xi*xi*xi)*
		 eta*eta+(-2048.0*xi+26624.0/3.0*xi*xi-12288.0*xi*xi*xi+16384.0/3.0*xi*xi*xi*xi)
		 *eta*eta*eta+(1024.0*xi-13312.0/3.0*xi*xi+6144.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*
			       xi)*eta*eta*eta*eta;
      case 24: return (144.0*xi-912.0*xi*xi+1536.0*xi*xi*xi-768.0*xi*xi*xi*xi
      )*eta+(-912.0*xi+5776.0*xi*xi-9728.0*xi*xi*xi+4864.0*xi*xi*xi*xi)*eta*eta+(
	1536.0*xi-9728.0*xi*xi+16384.0*xi*xi*xi-8192.0*xi*xi*xi*xi)*eta*eta*eta+(-768.0
										 *xi+4864.0*xi*xi-8192.0*xi*xi*xi+4096.0*xi*xi*xi*xi)*eta*eta*eta*eta;
    };
  return 0;
};



template <>
Tensor<1,2>
FEQ4<2>::shape_grad (const unsigned int i,
			     const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
				   // originally, the return type of the
				   // function was Point<dim>, so we
				   // still construct it as that. it should
				   // make no difference in practice,
				   // however
  switch (i)
    {
      case 0: return Point<2>(-25.0/3.0+140.0/3.0*xi-80.0*xi*xi+128.0/3.0*xi*xi*xi+(625.0/9.0-3500.0/9.0*xi+2000.0/3.0*xi*xi-3200.0/9.0*xi*xi*xi)*eta+(-1750.0/9.0+9800.0/9.0*xi-5600.0/3.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta*eta+(2000.0/9.0-11200.0/9.0*xi+6400.0/3.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta*eta+(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -25.0/3.0+625.0/9.0*xi-1750.0/9.0*xi*xi+2000.0/9.0*xi*xi*xi-800.0/9.0*xi*xi*xi*xi+2.0*(70.0/3.0-1750.0/9.0*xi+4900.0/9.0*xi*xi-5600.0/9.0*xi*xi*xi+2240.0/9.0*xi*xi*xi*xi)*eta+3.0*(-80.0/3.0+2000.0/9.0*xi-5600.0/9.0*xi*xi+6400.0/9.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(32.0/3.0-800.0/9.0*xi+2240.0/9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 1: return Point<2>(-1.0+44.0/3.0*xi-48.0*xi*xi+128.0/3.0*xi*xi*xi+(25.0/3.0-1100.0/9.0*xi+400.0*xi*xi-3200.0/9.0*xi*xi*xi)*eta+(-70.0/3.0+3080.0/9.0*xi-1120.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta*eta+(80.0/3.0-3520.0/9.0*xi+1280.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta*eta+(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      25.0/3.0*xi-550.0/9.0*xi*xi+400.0/3.0*xi*xi*xi-800.0/9.0*xi*xi*xi*xi+2.0*(-70.0/3.0*xi+1540.0/9.0*xi*xi-1120.0/3.0*xi*xi*xi+2240.0/9.0*xi*xi*xi*xi)*eta+3.0*(80.0/3.0*xi-1760.0/9.0*xi*xi+1280.0/3.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 2: return Point<2>((1.0-44.0/3.0*xi+48.0*xi*xi-128.0/3.0*xi*xi*xi)*eta+(-22.0/3.0+968.0/9.0*xi-352.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta*eta+(16.0-704.0/3.0*xi+768.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta*eta+(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      xi-22.0/3.0*xi*xi+16.0*xi*xi*xi-32.0/3.0*xi*xi*xi*xi+2.0*(-22.0/3.0*xi+484.0/9.0*xi*xi-352.0/3.0*xi*xi*xi+704.0/9.0*xi*xi*xi*xi)*eta+3.0*(16.0*xi-352.0/3.0*xi*xi+256.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 3: return Point<2>((25.0/3.0-140.0/3.0*xi+80.0*xi*xi-128.0/3.0*xi*xi*xi)*eta+(-550.0/9.0+3080.0/9.0*xi-1760.0/3.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta*eta+(400.0/3.0-2240.0/3.0*xi+1280.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta*eta+(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -1.0+25.0/3.0*xi-70.0/3.0*xi*xi+80.0/3.0*xi*xi*xi-32.0/3.0*xi*xi*xi*xi+2.0*(22.0/3.0-550.0/9.0*xi+1540.0/9.0*xi*xi-1760.0/9.0*xi*xi*xi+704.0/9.0*xi*xi*xi*xi)*eta+3.0*(-16.0+400.0/3.0*xi-1120.0/3.0*xi*xi+1280.0/3.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(32.0/3.0-800.0/9.0*xi+2240.0/9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 4: return Point<2>(16.0-416.0/3.0*xi+288.0*xi*xi-512.0/3.0*xi*xi*xi+(-400.0/3.0+10400.0/9.0*xi-2400.0*xi*xi+12800.0/9.0*xi*xi*xi)*eta+(1120.0/3.0-29120.0/9.0*xi+6720.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta*eta+(-1280.0/3.0+33280.0/9.0*xi-7680.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta*eta+(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -400.0/3.0*xi+5200.0/9.0*xi*xi-800.0*xi*xi*xi+3200.0/9.0*xi*xi*xi*xi+2.0*(1120.0/3.0*xi-14560.0/9.0*xi*xi+2240.0*xi*xi*xi-8960.0/9.0*xi*xi*xi*xi)*eta+3.0*(-1280.0/3.0*xi+16640.0/9.0*xi*xi-2560.0*xi*xi*xi+10240.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 5: return Point<2>(-12.0+152.0*xi-384.0*xi*xi+256.0*xi*xi*xi+(100.0-3800.0/3.0*xi+3200.0*xi*xi-6400.0/3.0*xi*xi*xi)*eta+(-280.0+10640.0/3.0*xi-8960.0*xi*xi+17920.0/3.0*xi*xi*xi)*eta*eta+(320.0-12160.0/3.0*xi+10240.0*xi*xi-20480.0/3.0*xi*xi*xi)*eta*eta*eta+(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			      100.0*xi-1900.0/3.0*xi*xi+3200.0/3.0*xi*xi*xi-1600.0/3.0*xi*xi*xi*xi+2.0*(-280.0*xi+5320.0/3.0*xi*xi-8960.0/3.0*xi*xi*xi+4480.0/3.0*xi*xi*xi*xi)*eta+3.0*(320.0*xi-6080.0/3.0*xi*xi+10240.0/3.0*xi*xi*xi-5120.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(-128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 6: return Point<2>(16.0/3.0-224.0/3.0*xi+224.0*xi*xi-512.0/3.0*xi*xi*xi+(-400.0/9.0+5600.0/9.0*xi-5600.0/3.0*xi*xi+12800.0/9.0*xi*xi*xi)*eta+(1120.0/9.0-15680.0/9.0*xi+15680.0/3.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta*eta+(-1280.0/9.0+17920.0/9.0*xi-17920.0/3.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta*eta+(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -400.0/9.0*xi+2800.0/9.0*xi*xi-5600.0/9.0*xi*xi*xi+3200.0/9.0*xi*xi*xi*xi+2.0*(1120.0/9.0*xi-7840.0/9.0*xi*xi+15680.0/9.0*xi*xi*xi-8960.0/9.0*xi*xi*xi*xi)*eta+3.0*(-1280.0/9.0*xi+8960.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+10240.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 7: return Point<2>((-16.0+704.0/3.0*xi-768.0*xi*xi+2048.0/3.0*xi*xi*xi)*eta+(208.0/3.0-9152.0/9.0*xi+3328.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta*eta+(-96.0+1408.0*xi-4608.0*xi*xi+4096.0*xi*xi*xi)*eta*eta*eta+(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -16.0*xi+352.0/3.0*xi*xi-256.0*xi*xi*xi+512.0/3.0*xi*xi*xi*xi+2.0*(208.0/3.0*xi-4576.0/9.0*xi*xi+3328.0/3.0*xi*xi*xi-6656.0/9.0*xi*xi*xi*xi)*eta+3.0*(-96.0*xi+704.0*xi*xi-1536.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta*eta+4.0*(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 8: return Point<2>((12.0-176.0*xi+576.0*xi*xi-512.0*xi*xi*xi)*eta+(-76.0+3344.0/3.0*xi-3648.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta*eta+(128.0-5632.0/3.0*xi+6144.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta*eta+(-64.0+2816.0/3.0*xi-3072.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			      12.0*xi-88.0*xi*xi+192.0*xi*xi*xi-128.0*xi*xi*xi*xi+2.0*(-76.0*xi+1672.0/3.0*xi*xi-1216.0*xi*xi*xi+2432.0/3.0*xi*xi*xi*xi)*eta+3.0*(128.0*xi-2816.0/3.0*xi*xi+2048.0*xi*xi*xi-4096.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(-64.0*xi+1408.0/3.0*xi*xi-1024.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 9: return Point<2>((-16.0/3.0+704.0/9.0*xi-256.0*xi*xi+2048.0/9.0*xi*xi*xi)*eta+(112.0/3.0-4928.0/9.0*xi+1792.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta*eta+(-224.0/3.0+9856.0/9.0*xi-3584.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta*eta+(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			      -16.0/3.0*xi+352.0/9.0*xi*xi-256.0/3.0*xi*xi*xi+512.0/9.0*xi*xi*xi*xi+2.0*(112.0/3.0*xi-2464.0/9.0*xi*xi+1792.0/3.0*xi*xi*xi-3584.0/9.0*xi*xi*xi*xi)*eta+3.0*(-224.0/3.0*xi+4928.0/9.0*xi*xi-3584.0/3.0*xi*xi*xi+7168.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 10: return Point<2>((-16.0+416.0/3.0*xi-288.0*xi*xi+512.0/3.0*xi*xi*xi)*eta+(352.0/3.0-9152.0/9.0*xi+2112.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta*eta+(-256.0+6656.0/3.0*xi-4608.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta+(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       -16.0*xi+208.0/3.0*xi*xi-96.0*xi*xi*xi+128.0/3.0*xi*xi*xi*xi+2.0*(352.0/3.0*xi-4576.0/9.0*xi*xi+704.0*xi*xi*xi-2816.0/9.0*xi*xi*xi*xi)*eta+3.0*(-256.0*xi+3328.0/3.0*xi*xi-1536.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 11: return Point<2>((12.0-152.0*xi+384.0*xi*xi-256.0*xi*xi*xi)*eta+(-88.0+3344.0/3.0*xi-2816.0*xi*xi+5632.0/3.0*xi*xi*xi)*eta*eta+(192.0-2432.0*xi+6144.0*xi*xi-4096.0*xi*xi*xi)*eta*eta*eta+(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       12.0*xi-76.0*xi*xi+128.0*xi*xi*xi-64.0*xi*xi*xi*xi+2.0*(-88.0*xi+1672.0/3.0*xi*xi-2816.0/3.0*xi*xi*xi+1408.0/3.0*xi*xi*xi*xi)*eta+3.0*(192.0*xi-1216.0*xi*xi+2048.0*xi*xi*xi-1024.0*xi*xi*xi*xi)*eta*eta+4.0*(-128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 12: return Point<2>((-16.0/3.0+224.0/3.0*xi-224.0*xi*xi+512.0/3.0*xi*xi*xi)*eta+(352.0/9.0-4928.0/9.0*xi+4928.0/3.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta*eta+(-256.0/3.0+3584.0/3.0*xi-3584.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta+(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       -16.0/3.0*xi+112.0/3.0*xi*xi-224.0/3.0*xi*xi*xi+128.0/3.0*xi*xi*xi*xi+2.0*(352.0/9.0*xi-2464.0/9.0*xi*xi+4928.0/9.0*xi*xi*xi-2816.0/9.0*xi*xi*xi*xi)*eta+3.0*(-256.0/3.0*xi+1792.0/3.0*xi*xi-3584.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 13: return Point<2>((-400.0/3.0+2240.0/3.0*xi-1280.0*xi*xi+2048.0/3.0*xi*xi*xi)*eta+(5200.0/9.0-29120.0/9.0*xi+16640.0/3.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta*eta+(-800.0+4480.0*xi-7680.0*xi*xi+4096.0*xi*xi*xi)*eta*eta*eta+(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       16.0-400.0/3.0*xi+1120.0/3.0*xi*xi-1280.0/3.0*xi*xi*xi+512.0/3.0*xi*xi*xi*xi+2.0*(-208.0/3.0+5200.0/9.0*xi-14560.0/9.0*xi*xi+16640.0/9.0*xi*xi*xi-6656.0/9.0*xi*xi*xi*xi)*eta+3.0*(96.0-800.0*xi+2240.0*xi*xi-2560.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta*eta+4.0*(-128.0/3.0+3200.0/9.0*xi-8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 14: return Point<2>((100.0-560.0*xi+960.0*xi*xi-512.0*xi*xi*xi)*eta+(-1900.0/3.0+10640.0/3.0*xi-6080.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta*eta+(3200.0/3.0-17920.0/3.0*xi+10240.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta*eta+(-1600.0/3.0+8960.0/3.0*xi-5120.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       -12.0+100.0*xi-280.0*xi*xi+320.0*xi*xi*xi-128.0*xi*xi*xi*xi+2.0*(76.0-1900.0/3.0*xi+5320.0/3.0*xi*xi-6080.0/3.0*xi*xi*xi+2432.0/3.0*xi*xi*xi*xi)*eta+3.0*(-128.0+3200.0/3.0*xi-8960.0/3.0*xi*xi+10240.0/3.0*xi*xi*xi-4096.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(64.0-1600.0/3.0*xi+4480.0/3.0*xi*xi-5120.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 15: return Point<2>((-400.0/9.0+2240.0/9.0*xi-1280.0/3.0*xi*xi+2048.0/9.0*xi*xi*xi)*eta+(2800.0/9.0-15680.0/9.0*xi+8960.0/3.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta*eta+(-5600.0/9.0+31360.0/9.0*xi-17920.0/3.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta*eta+(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       16.0/3.0-400.0/9.0*xi+1120.0/9.0*xi*xi-1280.0/9.0*xi*xi*xi+512.0/9.0*xi*xi*xi*xi+2.0*(-112.0/3.0+2800.0/9.0*xi-7840.0/9.0*xi*xi+8960.0/9.0*xi*xi*xi-3584.0/9.0*xi*xi*xi*xi)*eta+3.0*(224.0/3.0-5600.0/9.0*xi+15680.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+7168.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(-128.0/3.0+3200.0/9.0*xi-8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 16: return Point<2>((256.0-6656.0/3.0*xi+4608.0*xi*xi-8192.0/3.0*xi*xi*xi)*eta+(-3328.0/3.0+86528.0/9.0*xi-19968.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta*eta+(1536.0-13312.0*xi+27648.0*xi*xi-16384.0*xi*xi*xi)*eta*eta*eta+(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       256.0*xi-3328.0/3.0*xi*xi+1536.0*xi*xi*xi-2048.0/3.0*xi*xi*xi*xi+2.0*(-3328.0/3.0*xi+43264.0/9.0*xi*xi-6656.0*xi*xi*xi+26624.0/9.0*xi*xi*xi*xi)*eta+3.0*(1536.0*xi-6656.0*xi*xi+9216.0*xi*xi*xi-4096.0*xi*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 17: return Point<2>((256.0/3.0-3584.0/3.0*xi+3584.0*xi*xi-8192.0/3.0*xi*xi*xi)*eta+(-3328.0/9.0+46592.0/9.0*xi-46592.0/3.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta*eta+(512.0-7168.0*xi+21504.0*xi*xi-16384.0*xi*xi*xi)*eta*eta*eta+(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       256.0/3.0*xi-1792.0/3.0*xi*xi+3584.0/3.0*xi*xi*xi-2048.0/3.0*xi*xi*xi*xi+2.0*(-3328.0/9.0*xi+23296.0/9.0*xi*xi-46592.0/9.0*xi*xi*xi+26624.0/9.0*xi*xi*xi*xi)*eta+3.0*(512.0*xi-3584.0*xi*xi+7168.0*xi*xi*xi-4096.0*xi*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0*xi+14336.0/9.0*xi*xi-28672.0/9.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 18: return Point<2>((256.0/9.0-3584.0/9.0*xi+3584.0/3.0*xi*xi-8192.0/9.0*xi*xi*xi)*eta+(-1792.0/9.0+25088.0/9.0*xi-25088.0/3.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta*eta+(3584.0/9.0-50176.0/9.0*xi+50176.0/3.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta*eta+(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       256.0/9.0*xi-1792.0/9.0*xi*xi+3584.0/9.0*xi*xi*xi-2048.0/9.0*xi*xi*xi*xi+2.0*(-1792.0/9.0*xi+12544.0/9.0*xi*xi-25088.0/9.0*xi*xi*xi+14336.0/9.0*xi*xi*xi*xi)*eta+3.0*(3584.0/9.0*xi-25088.0/9.0*xi*xi+50176.0/9.0*xi*xi*xi-28672.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0*xi+14336.0/9.0*xi*xi-28672.0/9.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 19: return Point<2>((256.0/3.0-6656.0/9.0*xi+1536.0*xi*xi-8192.0/9.0*xi*xi*xi)*eta+(-1792.0/3.0+46592.0/9.0*xi-10752.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta*eta+(3584.0/3.0-93184.0/9.0*xi+21504.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta*eta+(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta*eta,
			       256.0/3.0*xi-3328.0/9.0*xi*xi+512.0*xi*xi*xi-2048.0/9.0*xi*xi*xi*xi+2.0*(-1792.0/3.0*xi+23296.0/9.0*xi*xi-3584.0*xi*xi*xi+14336.0/9.0*xi*xi*xi*xi)*eta+3.0*(3584.0/3.0*xi-46592.0/9.0*xi*xi+7168.0*xi*xi*xi-28672.0/9.0*xi*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta*eta);
      case 20: return Point<2>((-192.0+2432.0*xi-6144.0*xi*xi+4096.0*xi*xi*xi)*eta+(832.0-31616.0/3.0*xi+26624.0*xi*xi-53248.0/3.0*xi*xi*xi)*eta*eta+(-1152.0+14592.0*xi-36864.0*xi*xi+24576.0*xi*xi*xi)*eta*eta*eta+(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       -192.0*xi+1216.0*xi*xi-2048.0*xi*xi*xi+1024.0*xi*xi*xi*xi+2.0*(832.0*xi-15808.0/3.0*xi*xi+26624.0/3.0*xi*xi*xi-13312.0/3.0*xi*xi*xi*xi)*eta+3.0*(-1152.0*xi+7296.0*xi*xi-12288.0*xi*xi*xi+6144.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 21: return Point<2>((-64.0+896.0*xi-2688.0*xi*xi+2048.0*xi*xi*xi)*eta+(1216.0/3.0-17024.0/3.0*xi+17024.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta*eta+(-2048.0/3.0+28672.0/3.0*xi-28672.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta*eta+(1024.0/3.0-14336.0/3.0*xi+14336.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       -64.0*xi+448.0*xi*xi-896.0*xi*xi*xi+512.0*xi*xi*xi*xi+2.0*(1216.0/3.0*xi-8512.0/3.0*xi*xi+17024.0/3.0*xi*xi*xi-9728.0/3.0*xi*xi*xi*xi)*eta+3.0*(-2048.0/3.0*xi+14336.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+16384.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(1024.0/3.0*xi-7168.0/3.0*xi*xi+14336.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 22: return Point<2>((-64.0+2432.0/3.0*xi-2048.0*xi*xi+4096.0/3.0*xi*xi*xi)*eta+(448.0-17024.0/3.0*xi+14336.0*xi*xi-28672.0/3.0*xi*xi*xi)*eta*eta+(-896.0+34048.0/3.0*xi-28672.0*xi*xi+57344.0/3.0*xi*xi*xi)*eta*eta*eta+(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       -64.0*xi+1216.0/3.0*xi*xi-2048.0/3.0*xi*xi*xi+1024.0/3.0*xi*xi*xi*xi+2.0*(448.0*xi-8512.0/3.0*xi*xi+14336.0/3.0*xi*xi*xi-7168.0/3.0*xi*xi*xi*xi)*eta+3.0*(-896.0*xi+17024.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+14336.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 23: return Point<2>((-192.0+1664.0*xi-3456.0*xi*xi+2048.0*xi*xi*xi)*eta+(1216.0-31616.0/3.0*xi+21888.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta*eta+(-2048.0+53248.0/3.0*xi-36864.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta*eta+(1024.0-26624.0/3.0*xi+18432.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta*eta,
			       -192.0*xi+832.0*xi*xi-1152.0*xi*xi*xi+512.0*xi*xi*xi*xi+2.0*(1216.0*xi-15808.0/3.0*xi*xi+7296.0*xi*xi*xi-9728.0/3.0*xi*xi*xi*xi)*eta+3.0*(-2048.0*xi+26624.0/3.0*xi*xi-12288.0*xi*xi*xi+16384.0/3.0*xi*xi*xi*xi)*eta*eta+4.0*(1024.0*xi-13312.0/3.0*xi*xi+6144.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta*eta);
      case 24: return Point<2>((144.0-1824.0*xi+4608.0*xi*xi-3072.0*xi*xi*xi)*eta+(-912.0+11552.0*xi-29184.0*xi*xi+19456.0*xi*xi*xi)*eta*eta+(1536.0-19456.0*xi+49152.0*xi*xi-32768.0*xi*xi*xi)*eta*eta*eta+(-768.0+9728.0*xi-24576.0*xi*xi+16384.0*xi*xi*xi)*eta*eta*eta*eta,
			       144.0*xi-912.0*xi*xi+1536.0*xi*xi*xi-768.0*xi*xi*xi*xi+2.0*(-912.0*xi+5776.0*xi*xi-9728.0*xi*xi*xi+4864.0*xi*xi*xi*xi)*eta+3.0*(1536.0*xi-9728.0*xi*xi+16384.0*xi*xi*xi-8192.0*xi*xi*xi*xi)*eta*eta+4.0*(-768.0*xi+4864.0*xi*xi-8192.0*xi*xi*xi+4096.0*xi*xi*xi*xi)*eta*eta*eta);
    };
  return Point<2> ();
};



template <>
Tensor<2,2>
FEQ4<2>::shape_grad_grad (const unsigned int i,
				  const Point<2>    &p) const
{
  Assert (i<total_dofs, ExcInvalidIndex(i));

  const double xi = p(0),
	       eta= p(1);
  Tensor<2,2> return_value;
  
  switch (i)
    {
      case 0:
	    return_value[0][0] = 140.0/3.0-160.0*xi+128.0*xi*xi+(-3500.0/9.0+4000.0/3.0*xi-3200.0/3.0*xi*xi)*eta+(9800.0/9.0-11200.0/3.0*xi+8960.0/3.0*xi*xi)*eta*eta+(-11200.0/9.0+12800.0/3.0*xi-10240.0/3.0*xi*xi)*eta*eta*eta+(4480.0/9.0-5120.0/3.0*xi+4096.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 625.0/9.0-3500.0/9.0*xi+2000.0/3.0*xi*xi-3200.0/9.0*xi*xi*xi+2.0*(-1750.0/9.0+9800.0/9.0*xi-5600.0/3.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta+3.0*(2000.0/9.0-11200.0/9.0*xi+6400.0/3.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta+4.0*(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 625.0/9.0-3500.0/9.0*xi+2000.0/3.0*xi*xi-3200.0/9.0*xi*xi*xi+2.0*(-1750.0/9.0+9800.0/9.0*xi-5600.0/3.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta+3.0*(2000.0/9.0-11200.0/9.0*xi+6400.0/3.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta+4.0*(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 140.0/3.0-3500.0/9.0*xi+9800.0/9.0*xi*xi-11200.0/9.0*xi*xi*xi+4480.0/9.0*xi*xi*xi*xi+6.0*(-80.0/3.0+2000.0/9.0*xi-5600.0/9.0*xi*xi+6400.0/9.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*eta+12.0*(32.0/3.0-800.0/9.0*xi+2240.0/9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 1:
	    return_value[0][0] = 44.0/3.0-96.0*xi+128.0*xi*xi+(-1100.0/9.0+800.0*xi-3200.0/3.0*xi*xi)*eta+(3080.0/9.0-2240.0*xi+8960.0/3.0*xi*xi)*eta*eta+(-3520.0/9.0+2560.0*xi-10240.0/3.0*xi*xi)*eta*eta*eta+(1408.0/9.0-1024.0*xi+4096.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 25.0/3.0-1100.0/9.0*xi+400.0*xi*xi-3200.0/9.0*xi*xi*xi+2.0*(-70.0/3.0+3080.0/9.0*xi-1120.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta+3.0*(80.0/3.0-3520.0/9.0*xi+1280.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 25.0/3.0-1100.0/9.0*xi+400.0*xi*xi-3200.0/9.0*xi*xi*xi+2.0*(-70.0/3.0+3080.0/9.0*xi-1120.0*xi*xi+8960.0/9.0*xi*xi*xi)*eta+3.0*(80.0/3.0-3520.0/9.0*xi+1280.0*xi*xi-10240.0/9.0*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -140.0/3.0*xi+3080.0/9.0*xi*xi-2240.0/3.0*xi*xi*xi+4480.0/9.0*xi*xi*xi*xi+6.0*(80.0/3.0*xi-1760.0/9.0*xi*xi+1280.0/3.0*xi*xi*xi-2560.0/9.0*xi*xi*xi*xi)*eta+12.0*(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 2:
	    return_value[0][0] = (-44.0/3.0+96.0*xi-128.0*xi*xi)*eta+(968.0/9.0-704.0*xi+2816.0/3.0*xi*xi)*eta*eta+(-704.0/3.0+1536.0*xi-2048.0*xi*xi)*eta*eta*eta+(1408.0/9.0-1024.0*xi+4096.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 1.0-44.0/3.0*xi+48.0*xi*xi-128.0/3.0*xi*xi*xi+2.0*(-22.0/3.0+968.0/9.0*xi-352.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta+3.0*(16.0-704.0/3.0*xi+768.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 1.0-44.0/3.0*xi+48.0*xi*xi-128.0/3.0*xi*xi*xi+2.0*(-22.0/3.0+968.0/9.0*xi-352.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta+3.0*(16.0-704.0/3.0*xi+768.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta+4.0*(-32.0/3.0+1408.0/9.0*xi-512.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -44.0/3.0*xi+968.0/9.0*xi*xi-704.0/3.0*xi*xi*xi+1408.0/9.0*xi*xi*xi*xi+6.0*(16.0*xi-352.0/3.0*xi*xi+256.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta+12.0*(-32.0/3.0*xi+704.0/9.0*xi*xi-512.0/3.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 3:
	    return_value[0][0] = (-140.0/3.0+160.0*xi-128.0*xi*xi)*eta+(3080.0/9.0-3520.0/3.0*xi+2816.0/3.0*xi*xi)*eta*eta+(-2240.0/3.0+2560.0*xi-2048.0*xi*xi)*eta*eta*eta+(4480.0/9.0-5120.0/3.0*xi+4096.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 25.0/3.0-140.0/3.0*xi+80.0*xi*xi-128.0/3.0*xi*xi*xi+2.0*(-550.0/9.0+3080.0/9.0*xi-1760.0/3.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta+3.0*(400.0/3.0-2240.0/3.0*xi+1280.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta+4.0*(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 25.0/3.0-140.0/3.0*xi+80.0*xi*xi-128.0/3.0*xi*xi*xi+2.0*(-550.0/9.0+3080.0/9.0*xi-1760.0/3.0*xi*xi+2816.0/9.0*xi*xi*xi)*eta+3.0*(400.0/3.0-2240.0/3.0*xi+1280.0*xi*xi-2048.0/3.0*xi*xi*xi)*eta*eta+4.0*(-800.0/9.0+4480.0/9.0*xi-2560.0/3.0*xi*xi+4096.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 44.0/3.0-1100.0/9.0*xi+3080.0/9.0*xi*xi-3520.0/9.0*xi*xi*xi+1408.0/9.0*xi*xi*xi*xi+6.0*(-16.0+400.0/3.0*xi-1120.0/3.0*xi*xi+1280.0/3.0*xi*xi*xi-512.0/3.0*xi*xi*xi*xi)*eta+12.0*(32.0/3.0-800.0/9.0*xi+2240.0/9.0*xi*xi-2560.0/9.0*xi*xi*xi+1024.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 4:
	    return_value[0][0] = -416.0/3.0+576.0*xi-512.0*xi*xi+(10400.0/9.0-4800.0*xi+12800.0/3.0*xi*xi)*eta+(-29120.0/9.0+13440.0*xi-35840.0/3.0*xi*xi)*eta*eta+(33280.0/9.0-15360.0*xi+40960.0/3.0*xi*xi)*eta*eta*eta+(-13312.0/9.0+6144.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -400.0/3.0+10400.0/9.0*xi-2400.0*xi*xi+12800.0/9.0*xi*xi*xi+2.0*(1120.0/3.0-29120.0/9.0*xi+6720.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta+3.0*(-1280.0/3.0+33280.0/9.0*xi-7680.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta+4.0*(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -400.0/3.0+10400.0/9.0*xi-2400.0*xi*xi+12800.0/9.0*xi*xi*xi+2.0*(1120.0/3.0-29120.0/9.0*xi+6720.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta+3.0*(-1280.0/3.0+33280.0/9.0*xi-7680.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta+4.0*(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 2240.0/3.0*xi-29120.0/9.0*xi*xi+4480.0*xi*xi*xi-17920.0/9.0*xi*xi*xi*xi+6.0*(-1280.0/3.0*xi+16640.0/9.0*xi*xi-2560.0*xi*xi*xi+10240.0/9.0*xi*xi*xi*xi)*eta+12.0*(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 5:
	    return_value[0][0] = 152.0-768.0*xi+768.0*xi*xi+(-3800.0/3.0+6400.0*xi-6400.0*xi*xi)*eta+(10640.0/3.0-17920.0*xi+17920.0*xi*xi)*eta*eta+(-12160.0/3.0+20480.0*xi-20480.0*xi*xi)*eta*eta*eta+(4864.0/3.0-8192.0*xi+8192.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 100.0-3800.0/3.0*xi+3200.0*xi*xi-6400.0/3.0*xi*xi*xi+2.0*(-280.0+10640.0/3.0*xi-8960.0*xi*xi+17920.0/3.0*xi*xi*xi)*eta+3.0*(320.0-12160.0/3.0*xi+10240.0*xi*xi-20480.0/3.0*xi*xi*xi)*eta*eta+4.0*(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 100.0-3800.0/3.0*xi+3200.0*xi*xi-6400.0/3.0*xi*xi*xi+2.0*(-280.0+10640.0/3.0*xi-8960.0*xi*xi+17920.0/3.0*xi*xi*xi)*eta+3.0*(320.0-12160.0/3.0*xi+10240.0*xi*xi-20480.0/3.0*xi*xi*xi)*eta*eta+4.0*(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -560.0*xi+10640.0/3.0*xi*xi-17920.0/3.0*xi*xi*xi+8960.0/3.0*xi*xi*xi*xi+6.0*(320.0*xi-6080.0/3.0*xi*xi+10240.0/3.0*xi*xi*xi-5120.0/3.0*xi*xi*xi*xi)*eta+12.0*(-128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 6:
	    return_value[0][0] = -224.0/3.0+448.0*xi-512.0*xi*xi+(5600.0/9.0-11200.0/3.0*xi+12800.0/3.0*xi*xi)*eta+(-15680.0/9.0+31360.0/3.0*xi-35840.0/3.0*xi*xi)*eta*eta+(17920.0/9.0-35840.0/3.0*xi+40960.0/3.0*xi*xi)*eta*eta*eta+(-7168.0/9.0+14336.0/3.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -400.0/9.0+5600.0/9.0*xi-5600.0/3.0*xi*xi+12800.0/9.0*xi*xi*xi+2.0*(1120.0/9.0-15680.0/9.0*xi+15680.0/3.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta+3.0*(-1280.0/9.0+17920.0/9.0*xi-17920.0/3.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta+4.0*(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -400.0/9.0+5600.0/9.0*xi-5600.0/3.0*xi*xi+12800.0/9.0*xi*xi*xi+2.0*(1120.0/9.0-15680.0/9.0*xi+15680.0/3.0*xi*xi-35840.0/9.0*xi*xi*xi)*eta+3.0*(-1280.0/9.0+17920.0/9.0*xi-17920.0/3.0*xi*xi+40960.0/9.0*xi*xi*xi)*eta*eta+4.0*(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 2240.0/9.0*xi-15680.0/9.0*xi*xi+31360.0/9.0*xi*xi*xi-17920.0/9.0*xi*xi*xi*xi+6.0*(-1280.0/9.0*xi+8960.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+10240.0/9.0*xi*xi*xi*xi)*eta+12.0*(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 7:
	    return_value[0][0] = (704.0/3.0-1536.0*xi+2048.0*xi*xi)*eta+(-9152.0/9.0+6656.0*xi-26624.0/3.0*xi*xi)*eta*eta+(1408.0-9216.0*xi+12288.0*xi*xi)*eta*eta*eta+(-5632.0/9.0+4096.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -16.0+704.0/3.0*xi-768.0*xi*xi+2048.0/3.0*xi*xi*xi+2.0*(208.0/3.0-9152.0/9.0*xi+3328.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta+3.0*(-96.0+1408.0*xi-4608.0*xi*xi+4096.0*xi*xi*xi)*eta*eta+4.0*(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -16.0+704.0/3.0*xi-768.0*xi*xi+2048.0/3.0*xi*xi*xi+2.0*(208.0/3.0-9152.0/9.0*xi+3328.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta+3.0*(-96.0+1408.0*xi-4608.0*xi*xi+4096.0*xi*xi*xi)*eta*eta+4.0*(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 416.0/3.0*xi-9152.0/9.0*xi*xi+6656.0/3.0*xi*xi*xi-13312.0/9.0*xi*xi*xi*xi+6.0*(-96.0*xi+704.0*xi*xi-1536.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta+12.0*(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 8:
	    return_value[0][0] = (-176.0+1152.0*xi-1536.0*xi*xi)*eta+(3344.0/3.0-7296.0*xi+9728.0*xi*xi)*eta*eta+(-5632.0/3.0+12288.0*xi-16384.0*xi*xi)*eta*eta*eta+(2816.0/3.0-6144.0*xi+8192.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 12.0-176.0*xi+576.0*xi*xi-512.0*xi*xi*xi+2.0*(-76.0+3344.0/3.0*xi-3648.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta+3.0*(128.0-5632.0/3.0*xi+6144.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta+4.0*(-64.0+2816.0/3.0*xi-3072.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 12.0-176.0*xi+576.0*xi*xi-512.0*xi*xi*xi+2.0*(-76.0+3344.0/3.0*xi-3648.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta+3.0*(128.0-5632.0/3.0*xi+6144.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta+4.0*(-64.0+2816.0/3.0*xi-3072.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -152.0*xi+3344.0/3.0*xi*xi-2432.0*xi*xi*xi+4864.0/3.0*xi*xi*xi*xi+6.0*(128.0*xi-2816.0/3.0*xi*xi+2048.0*xi*xi*xi-4096.0/3.0*xi*xi*xi*xi)*eta+12.0*(-64.0*xi+1408.0/3.0*xi*xi-1024.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta;
	    break;
      case 9:
	    return_value[0][0] = (704.0/9.0-512.0*xi+2048.0/3.0*xi*xi)*eta+(-4928.0/9.0+3584.0*xi-14336.0/3.0*xi*xi)*eta*eta+(9856.0/9.0-7168.0*xi+28672.0/3.0*xi*xi)*eta*eta*eta+(-5632.0/9.0+4096.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -16.0/3.0+704.0/9.0*xi-256.0*xi*xi+2048.0/9.0*xi*xi*xi+2.0*(112.0/3.0-4928.0/9.0*xi+1792.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta+3.0*(-224.0/3.0+9856.0/9.0*xi-3584.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta+4.0*(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -16.0/3.0+704.0/9.0*xi-256.0*xi*xi+2048.0/9.0*xi*xi*xi+2.0*(112.0/3.0-4928.0/9.0*xi+1792.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta+3.0*(-224.0/3.0+9856.0/9.0*xi-3584.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta+4.0*(128.0/3.0-5632.0/9.0*xi+2048.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 224.0/3.0*xi-4928.0/9.0*xi*xi+3584.0/3.0*xi*xi*xi-7168.0/9.0*xi*xi*xi*xi+6.0*(-224.0/3.0*xi+4928.0/9.0*xi*xi-3584.0/3.0*xi*xi*xi+7168.0/9.0*xi*xi*xi*xi)*eta+12.0*(128.0/3.0*xi-2816.0/9.0*xi*xi+2048.0/3.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (416.0/3.0-576.0*xi+512.0*xi*xi)*eta+(-9152.0/9.0+4224.0*xi-11264.0/3.0*xi*xi)*eta*eta+(6656.0/3.0-9216.0*xi+8192.0*xi*xi)*eta*eta*eta+(-13312.0/9.0+6144.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -16.0+416.0/3.0*xi-288.0*xi*xi+512.0/3.0*xi*xi*xi+2.0*(352.0/3.0-9152.0/9.0*xi+2112.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta+3.0*(-256.0+6656.0/3.0*xi-4608.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -16.0+416.0/3.0*xi-288.0*xi*xi+512.0/3.0*xi*xi*xi+2.0*(352.0/3.0-9152.0/9.0*xi+2112.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta+3.0*(-256.0+6656.0/3.0*xi-4608.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0/3.0-13312.0/9.0*xi+3072.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 704.0/3.0*xi-9152.0/9.0*xi*xi+1408.0*xi*xi*xi-5632.0/9.0*xi*xi*xi*xi+6.0*(-256.0*xi+3328.0/3.0*xi*xi-1536.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta+12.0*(512.0/3.0*xi-6656.0/9.0*xi*xi+1024.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-152.0+768.0*xi-768.0*xi*xi)*eta+(3344.0/3.0-5632.0*xi+5632.0*xi*xi)*eta*eta+(-2432.0+12288.0*xi-12288.0*xi*xi)*eta*eta*eta+(4864.0/3.0-8192.0*xi+8192.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 12.0-152.0*xi+384.0*xi*xi-256.0*xi*xi*xi+2.0*(-88.0+3344.0/3.0*xi-2816.0*xi*xi+5632.0/3.0*xi*xi*xi)*eta+3.0*(192.0-2432.0*xi+6144.0*xi*xi-4096.0*xi*xi*xi)*eta*eta+4.0*(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 12.0-152.0*xi+384.0*xi*xi-256.0*xi*xi*xi+2.0*(-88.0+3344.0/3.0*xi-2816.0*xi*xi+5632.0/3.0*xi*xi*xi)*eta+3.0*(192.0-2432.0*xi+6144.0*xi*xi-4096.0*xi*xi*xi)*eta*eta+4.0*(-128.0+4864.0/3.0*xi-4096.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -176.0*xi+3344.0/3.0*xi*xi-5632.0/3.0*xi*xi*xi+2816.0/3.0*xi*xi*xi*xi+6.0*(192.0*xi-1216.0*xi*xi+2048.0*xi*xi*xi-1024.0*xi*xi*xi*xi)*eta+12.0*(-128.0*xi+2432.0/3.0*xi*xi-4096.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (224.0/3.0-448.0*xi+512.0*xi*xi)*eta+(-4928.0/9.0+9856.0/3.0*xi-11264.0/3.0*xi*xi)*eta*eta+(3584.0/3.0-7168.0*xi+8192.0*xi*xi)*eta*eta*eta+(-7168.0/9.0+14336.0/3.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -16.0/3.0+224.0/3.0*xi-224.0*xi*xi+512.0/3.0*xi*xi*xi+2.0*(352.0/9.0-4928.0/9.0*xi+4928.0/3.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta+3.0*(-256.0/3.0+3584.0/3.0*xi-3584.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -16.0/3.0+224.0/3.0*xi-224.0*xi*xi+512.0/3.0*xi*xi*xi+2.0*(352.0/9.0-4928.0/9.0*xi+4928.0/3.0*xi*xi-11264.0/9.0*xi*xi*xi)*eta+3.0*(-256.0/3.0+3584.0/3.0*xi-3584.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0/9.0-7168.0/9.0*xi+7168.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 704.0/9.0*xi-4928.0/9.0*xi*xi+9856.0/9.0*xi*xi*xi-5632.0/9.0*xi*xi*xi*xi+6.0*(-256.0/3.0*xi+1792.0/3.0*xi*xi-3584.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta+12.0*(512.0/9.0*xi-3584.0/9.0*xi*xi+7168.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (2240.0/3.0-2560.0*xi+2048.0*xi*xi)*eta+(-29120.0/9.0+33280.0/3.0*xi-26624.0/3.0*xi*xi)*eta*eta+(4480.0-15360.0*xi+12288.0*xi*xi)*eta*eta*eta+(-17920.0/9.0+20480.0/3.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -400.0/3.0+2240.0/3.0*xi-1280.0*xi*xi+2048.0/3.0*xi*xi*xi+2.0*(5200.0/9.0-29120.0/9.0*xi+16640.0/3.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta+3.0*(-800.0+4480.0*xi-7680.0*xi*xi+4096.0*xi*xi*xi)*eta*eta+4.0*(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -400.0/3.0+2240.0/3.0*xi-1280.0*xi*xi+2048.0/3.0*xi*xi*xi+2.0*(5200.0/9.0-29120.0/9.0*xi+16640.0/3.0*xi*xi-26624.0/9.0*xi*xi*xi)*eta+3.0*(-800.0+4480.0*xi-7680.0*xi*xi+4096.0*xi*xi*xi)*eta*eta+4.0*(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -416.0/3.0+10400.0/9.0*xi-29120.0/9.0*xi*xi+33280.0/9.0*xi*xi*xi-13312.0/9.0*xi*xi*xi*xi+6.0*(96.0-800.0*xi+2240.0*xi*xi-2560.0*xi*xi*xi+1024.0*xi*xi*xi*xi)*eta+12.0*(-128.0/3.0+3200.0/9.0*xi-8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-560.0+1920.0*xi-1536.0*xi*xi)*eta+(10640.0/3.0-12160.0*xi+9728.0*xi*xi)*eta*eta+(-17920.0/3.0+20480.0*xi-16384.0*xi*xi)*eta*eta*eta+(8960.0/3.0-10240.0*xi+8192.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 100.0-560.0*xi+960.0*xi*xi-512.0*xi*xi*xi+2.0*(-1900.0/3.0+10640.0/3.0*xi-6080.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta+3.0*(3200.0/3.0-17920.0/3.0*xi+10240.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta+4.0*(-1600.0/3.0+8960.0/3.0*xi-5120.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 100.0-560.0*xi+960.0*xi*xi-512.0*xi*xi*xi+2.0*(-1900.0/3.0+10640.0/3.0*xi-6080.0*xi*xi+9728.0/3.0*xi*xi*xi)*eta+3.0*(3200.0/3.0-17920.0/3.0*xi+10240.0*xi*xi-16384.0/3.0*xi*xi*xi)*eta*eta+4.0*(-1600.0/3.0+8960.0/3.0*xi-5120.0*xi*xi+8192.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 152.0-3800.0/3.0*xi+10640.0/3.0*xi*xi-12160.0/3.0*xi*xi*xi+4864.0/3.0*xi*xi*xi*xi+6.0*(-128.0+3200.0/3.0*xi-8960.0/3.0*xi*xi+10240.0/3.0*xi*xi*xi-4096.0/3.0*xi*xi*xi*xi)*eta+12.0*(64.0-1600.0/3.0*xi+4480.0/3.0*xi*xi-5120.0/3.0*xi*xi*xi+2048.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (2240.0/9.0-2560.0/3.0*xi+2048.0/3.0*xi*xi)*eta+(-15680.0/9.0+17920.0/3.0*xi-14336.0/3.0*xi*xi)*eta*eta+(31360.0/9.0-35840.0/3.0*xi+28672.0/3.0*xi*xi)*eta*eta*eta+(-17920.0/9.0+20480.0/3.0*xi-16384.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -400.0/9.0+2240.0/9.0*xi-1280.0/3.0*xi*xi+2048.0/9.0*xi*xi*xi+2.0*(2800.0/9.0-15680.0/9.0*xi+8960.0/3.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta+3.0*(-5600.0/9.0+31360.0/9.0*xi-17920.0/3.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta+4.0*(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -400.0/9.0+2240.0/9.0*xi-1280.0/3.0*xi*xi+2048.0/9.0*xi*xi*xi+2.0*(2800.0/9.0-15680.0/9.0*xi+8960.0/3.0*xi*xi-14336.0/9.0*xi*xi*xi)*eta+3.0*(-5600.0/9.0+31360.0/9.0*xi-17920.0/3.0*xi*xi+28672.0/9.0*xi*xi*xi)*eta*eta+4.0*(3200.0/9.0-17920.0/9.0*xi+10240.0/3.0*xi*xi-16384.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -224.0/3.0+5600.0/9.0*xi-15680.0/9.0*xi*xi+17920.0/9.0*xi*xi*xi-7168.0/9.0*xi*xi*xi*xi+6.0*(224.0/3.0-5600.0/9.0*xi+15680.0/9.0*xi*xi-17920.0/9.0*xi*xi*xi+7168.0/9.0*xi*xi*xi*xi)*eta+12.0*(-128.0/3.0+3200.0/9.0*xi-8960.0/9.0*xi*xi+10240.0/9.0*xi*xi*xi-4096.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-6656.0/3.0+9216.0*xi-8192.0*xi*xi)*eta+(86528.0/9.0-39936.0*xi+106496.0/3.0*xi*xi)*eta*eta+(-13312.0+55296.0*xi-49152.0*xi*xi)*eta*eta*eta+(53248.0/9.0-24576.0*xi+65536.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 256.0-6656.0/3.0*xi+4608.0*xi*xi-8192.0/3.0*xi*xi*xi+2.0*(-3328.0/3.0+86528.0/9.0*xi-19968.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta+3.0*(1536.0-13312.0*xi+27648.0*xi*xi-16384.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 256.0-6656.0/3.0*xi+4608.0*xi*xi-8192.0/3.0*xi*xi*xi+2.0*(-3328.0/3.0+86528.0/9.0*xi-19968.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta+3.0*(1536.0-13312.0*xi+27648.0*xi*xi-16384.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -6656.0/3.0*xi+86528.0/9.0*xi*xi-13312.0*xi*xi*xi+53248.0/9.0*xi*xi*xi*xi+6.0*(1536.0*xi-6656.0*xi*xi+9216.0*xi*xi*xi-4096.0*xi*xi*xi*xi)*eta+12.0*(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-3584.0/3.0+7168.0*xi-8192.0*xi*xi)*eta+(46592.0/9.0-93184.0/3.0*xi+106496.0/3.0*xi*xi)*eta*eta+(-7168.0+43008.0*xi-49152.0*xi*xi)*eta*eta*eta+(28672.0/9.0-57344.0/3.0*xi+65536.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 256.0/3.0-3584.0/3.0*xi+3584.0*xi*xi-8192.0/3.0*xi*xi*xi+2.0*(-3328.0/9.0+46592.0/9.0*xi-46592.0/3.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta+3.0*(512.0-7168.0*xi+21504.0*xi*xi-16384.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 256.0/3.0-3584.0/3.0*xi+3584.0*xi*xi-8192.0/3.0*xi*xi*xi+2.0*(-3328.0/9.0+46592.0/9.0*xi-46592.0/3.0*xi*xi+106496.0/9.0*xi*xi*xi)*eta+3.0*(512.0-7168.0*xi+21504.0*xi*xi-16384.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -6656.0/9.0*xi+46592.0/9.0*xi*xi-93184.0/9.0*xi*xi*xi+53248.0/9.0*xi*xi*xi*xi+6.0*(512.0*xi-3584.0*xi*xi+7168.0*xi*xi*xi-4096.0*xi*xi*xi*xi)*eta+12.0*(-2048.0/9.0*xi+14336.0/9.0*xi*xi-28672.0/9.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-3584.0/9.0+7168.0/3.0*xi-8192.0/3.0*xi*xi)*eta+(25088.0/9.0-50176.0/3.0*xi+57344.0/3.0*xi*xi)*eta*eta+(-50176.0/9.0+100352.0/3.0*xi-114688.0/3.0*xi*xi)*eta*eta*eta+(28672.0/9.0-57344.0/3.0*xi+65536.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 256.0/9.0-3584.0/9.0*xi+3584.0/3.0*xi*xi-8192.0/9.0*xi*xi*xi+2.0*(-1792.0/9.0+25088.0/9.0*xi-25088.0/3.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta+3.0*(3584.0/9.0-50176.0/9.0*xi+50176.0/3.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 256.0/9.0-3584.0/9.0*xi+3584.0/3.0*xi*xi-8192.0/9.0*xi*xi*xi+2.0*(-1792.0/9.0+25088.0/9.0*xi-25088.0/3.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta+3.0*(3584.0/9.0-50176.0/9.0*xi+50176.0/3.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/9.0+28672.0/9.0*xi-28672.0/3.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -3584.0/9.0*xi+25088.0/9.0*xi*xi-50176.0/9.0*xi*xi*xi+28672.0/9.0*xi*xi*xi*xi+6.0*(3584.0/9.0*xi-25088.0/9.0*xi*xi+50176.0/9.0*xi*xi*xi-28672.0/9.0*xi*xi*xi*xi)*eta+12.0*(-2048.0/9.0*xi+14336.0/9.0*xi*xi-28672.0/9.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-6656.0/9.0+3072.0*xi-8192.0/3.0*xi*xi)*eta+(46592.0/9.0-21504.0*xi+57344.0/3.0*xi*xi)*eta*eta+(-93184.0/9.0+43008.0*xi-114688.0/3.0*xi*xi)*eta*eta*eta+(53248.0/9.0-24576.0*xi+65536.0/3.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 256.0/3.0-6656.0/9.0*xi+1536.0*xi*xi-8192.0/9.0*xi*xi*xi+2.0*(-1792.0/3.0+46592.0/9.0*xi-10752.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta+3.0*(3584.0/3.0-93184.0/9.0*xi+21504.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 256.0/3.0-6656.0/9.0*xi+1536.0*xi*xi-8192.0/9.0*xi*xi*xi+2.0*(-1792.0/3.0+46592.0/9.0*xi-10752.0*xi*xi+57344.0/9.0*xi*xi*xi)*eta+3.0*(3584.0/3.0-93184.0/9.0*xi+21504.0*xi*xi-114688.0/9.0*xi*xi*xi)*eta*eta+4.0*(-2048.0/3.0+53248.0/9.0*xi-12288.0*xi*xi+65536.0/9.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -3584.0/3.0*xi+46592.0/9.0*xi*xi-7168.0*xi*xi*xi+28672.0/9.0*xi*xi*xi*xi+6.0*(3584.0/3.0*xi-46592.0/9.0*xi*xi+7168.0*xi*xi*xi-28672.0/9.0*xi*xi*xi*xi)*eta+12.0*(-2048.0/3.0*xi+26624.0/9.0*xi*xi-4096.0*xi*xi*xi+16384.0/9.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (2432.0-12288.0*xi+12288.0*xi*xi)*eta+(-31616.0/3.0+53248.0*xi-53248.0*xi*xi)*eta*eta+(14592.0-73728.0*xi+73728.0*xi*xi)*eta*eta*eta+(-19456.0/3.0+32768.0*xi-32768.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -192.0+2432.0*xi-6144.0*xi*xi+4096.0*xi*xi*xi+2.0*(832.0-31616.0/3.0*xi+26624.0*xi*xi-53248.0/3.0*xi*xi*xi)*eta+3.0*(-1152.0+14592.0*xi-36864.0*xi*xi+24576.0*xi*xi*xi)*eta*eta+4.0*(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -192.0+2432.0*xi-6144.0*xi*xi+4096.0*xi*xi*xi+2.0*(832.0-31616.0/3.0*xi+26624.0*xi*xi-53248.0/3.0*xi*xi*xi)*eta+3.0*(-1152.0+14592.0*xi-36864.0*xi*xi+24576.0*xi*xi*xi)*eta*eta+4.0*(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 1664.0*xi-31616.0/3.0*xi*xi+53248.0/3.0*xi*xi*xi-26624.0/3.0*xi*xi*xi*xi+6.0*(-1152.0*xi+7296.0*xi*xi-12288.0*xi*xi*xi+6144.0*xi*xi*xi*xi)*eta+12.0*(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (896.0-5376.0*xi+6144.0*xi*xi)*eta+(-17024.0/3.0+34048.0*xi-38912.0*xi*xi)*eta*eta+(28672.0/3.0-57344.0*xi+65536.0*xi*xi)*eta*eta*eta+(-14336.0/3.0+28672.0*xi-32768.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -64.0+896.0*xi-2688.0*xi*xi+2048.0*xi*xi*xi+2.0*(1216.0/3.0-17024.0/3.0*xi+17024.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta+3.0*(-2048.0/3.0+28672.0/3.0*xi-28672.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta+4.0*(1024.0/3.0-14336.0/3.0*xi+14336.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -64.0+896.0*xi-2688.0*xi*xi+2048.0*xi*xi*xi+2.0*(1216.0/3.0-17024.0/3.0*xi+17024.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta+3.0*(-2048.0/3.0+28672.0/3.0*xi-28672.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta+4.0*(1024.0/3.0-14336.0/3.0*xi+14336.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 2432.0/3.0*xi-17024.0/3.0*xi*xi+34048.0/3.0*xi*xi*xi-19456.0/3.0*xi*xi*xi*xi+6.0*(-2048.0/3.0*xi+14336.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+16384.0/3.0*xi*xi*xi*xi)*eta+12.0*(1024.0/3.0*xi-7168.0/3.0*xi*xi+14336.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (2432.0/3.0-4096.0*xi+4096.0*xi*xi)*eta+(-17024.0/3.0+28672.0*xi-28672.0*xi*xi)*eta*eta+(34048.0/3.0-57344.0*xi+57344.0*xi*xi)*eta*eta*eta+(-19456.0/3.0+32768.0*xi-32768.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -64.0+2432.0/3.0*xi-2048.0*xi*xi+4096.0/3.0*xi*xi*xi+2.0*(448.0-17024.0/3.0*xi+14336.0*xi*xi-28672.0/3.0*xi*xi*xi)*eta+3.0*(-896.0+34048.0/3.0*xi-28672.0*xi*xi+57344.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -64.0+2432.0/3.0*xi-2048.0*xi*xi+4096.0/3.0*xi*xi*xi+2.0*(448.0-17024.0/3.0*xi+14336.0*xi*xi-28672.0/3.0*xi*xi*xi)*eta+3.0*(-896.0+34048.0/3.0*xi-28672.0*xi*xi+57344.0/3.0*xi*xi*xi)*eta*eta+4.0*(512.0-19456.0/3.0*xi+16384.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 896.0*xi-17024.0/3.0*xi*xi+28672.0/3.0*xi*xi*xi-14336.0/3.0*xi*xi*xi*xi+6.0*(-896.0*xi+17024.0/3.0*xi*xi-28672.0/3.0*xi*xi*xi+14336.0/3.0*xi*xi*xi*xi)*eta+12.0*(512.0*xi-9728.0/3.0*xi*xi+16384.0/3.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (1664.0-6912.0*xi+6144.0*xi*xi)*eta+(-31616.0/3.0+43776.0*xi-38912.0*xi*xi)*eta*eta+(53248.0/3.0-73728.0*xi+65536.0*xi*xi)*eta*eta*eta+(-26624.0/3.0+36864.0*xi-32768.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = -192.0+1664.0*xi-3456.0*xi*xi+2048.0*xi*xi*xi+2.0*(1216.0-31616.0/3.0*xi+21888.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta+3.0*(-2048.0+53248.0/3.0*xi-36864.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta+4.0*(1024.0-26624.0/3.0*xi+18432.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = -192.0+1664.0*xi-3456.0*xi*xi+2048.0*xi*xi*xi+2.0*(1216.0-31616.0/3.0*xi+21888.0*xi*xi-38912.0/3.0*xi*xi*xi)*eta+3.0*(-2048.0+53248.0/3.0*xi-36864.0*xi*xi+65536.0/3.0*xi*xi*xi)*eta*eta+4.0*(1024.0-26624.0/3.0*xi+18432.0*xi*xi-32768.0/3.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = 2432.0*xi-31616.0/3.0*xi*xi+14592.0*xi*xi*xi-19456.0/3.0*xi*xi*xi*xi+6.0*(-2048.0*xi+26624.0/3.0*xi*xi-12288.0*xi*xi*xi+16384.0/3.0*xi*xi*xi*xi)*eta+12.0*(1024.0*xi-13312.0/3.0*xi*xi+6144.0*xi*xi*xi-8192.0/3.0*xi*xi*xi*xi)*eta*eta;
	    return_value[0][0] = (-1824.0+9216.0*xi-9216.0*xi*xi)*eta+(11552.0-58368.0*xi+58368.0*xi*xi)*eta*eta+(-19456.0+98304.0*xi-98304.0*xi*xi)*eta*eta*eta+(9728.0-49152.0*xi+49152.0*xi*xi)*eta*eta*eta*eta;
	    return_value[0][1] = 144.0-1824.0*xi+4608.0*xi*xi-3072.0*xi*xi*xi+2.0*(-912.0+11552.0*xi-29184.0*xi*xi+19456.0*xi*xi*xi)*eta+3.0*(1536.0-19456.0*xi+49152.0*xi*xi-32768.0*xi*xi*xi)*eta*eta+4.0*(-768.0+9728.0*xi-24576.0*xi*xi+16384.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][0] = 144.0-1824.0*xi+4608.0*xi*xi-3072.0*xi*xi*xi+2.0*(-912.0+11552.0*xi-29184.0*xi*xi+19456.0*xi*xi*xi)*eta+3.0*(1536.0-19456.0*xi+49152.0*xi*xi-32768.0*xi*xi*xi)*eta*eta+4.0*(-768.0+9728.0*xi-24576.0*xi*xi+16384.0*xi*xi*xi)*eta*eta*eta;
	    return_value[1][1] = -1824.0*xi+11552.0*xi*xi-19456.0*xi*xi*xi+9728.0*xi*xi*xi*xi+6.0*(1536.0*xi-9728.0*xi*xi+16384.0*xi*xi*xi-8192.0*xi*xi*xi*xi)*eta+12.0*(-768.0*xi+4864.0*xi*xi-8192.0*xi*xi*xi+4096.0*xi*xi*xi*xi)*eta*eta;
	    break;
    };
  return return_value;
};



template <>
void FEQ4<2>::get_local_mass_matrix (const DoFHandler<2>::cell_iterator &cell,
					     FullMatrix<double> &local_mass_matrix) const {
  Assert (local_mass_matrix.n() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.n(),total_dofs));
  Assert (local_mass_matrix.m() == total_dofs,
	  ExcWrongFieldDimension(local_mass_matrix.m(),total_dofs));

  const double x[4] = { cell->vertex(0)(0),
			cell->vertex(1)(0),
			cell->vertex(2)(0),
			cell->vertex(3)(0)  };
  const double y[4] = { cell->vertex(0)(1),
			cell->vertex(1)(1),
			cell->vertex(2)(1),
			cell->vertex(3)(1)  };
  
/* check that the Jacobi determinant

   t0 = (-x[0]*(1.0-eta)+x[1]*(1.0-eta)+x[2]*eta-x[3]*eta) *
   (-y[0]*(1.0-xi)-y[1]*xi+y[2]*xi+y[3]*(1.0-xi))        -
   (-x[0]*(1.0-xi)-x[1]*xi+x[2]*xi+x[3]*(1.0-xi)) *
   (-y[0]*(1.0-eta)+y[1]*(1.0-eta)+y[2]*eta-y[3]*eta)

   has the right sign.	 
	 
   We do not attempt to check its (hopefully) positive sign at all points
   on the unit cell, but we check that it is positive in the four corners,
   which is sufficient since $det J$ is a bilinear function.
*/
  Assert ((-x[0]+x[1])*(-y[0]+y[3])-(-x[0]+x[3])*(-y[0]+y[1]),  // xi=eta=0
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((x[2]-x[3])*(-y[0]+y[3])-(-x[0]+x[3])*(y[2]-y[3]),    // xi=0, eta=1
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((x[2]-x[3])*(-y[1]+y[2])-(-x[1]+x[2])*(y[2]-y[3]),    // xi=eta=1
	  ExcJacobiDeterminantHasWrongSign());
  Assert ((-x[0]+x[1])*(-y[1]+y[2])-(-x[1]+x[2])*(-y[0]+y[1]),  // xi=1, eta=0
	  ExcJacobiDeterminantHasWrongSign());

  const double t1 = -x[0]+x[1];
  const double t2 = y[0]-y[1]+y[2]-y[3];
  const double t3 = t1*t2;
  const double t4 = 2117.0/16074450.0*t3;
  const double t5 = x[0]-x[1]+x[2]-x[3];
  const double t6 = -y[0]+y[1];
  const double t7 = t5*t6;
  const double t8 = 2117.0/16074450.0*t7;
  const double t9 = -y[0]+y[3];
  const double t10 = t1*t9;
  const double t11 = 21316.0/8037225.0*t10;
  const double t12 = t5*t9;
  const double t13 = 2117.0/16074450.0*t12;
  const double t14 = -x[0]+x[3];
  const double t15 = t14*t6;
  const double t16 = 21316.0/8037225.0*t15;
  const double t17 = t14*t2;
  const double t18 = 2117.0/16074450.0*t17;
  const double t20 = 2117.0/8037225.0*t10;
  const double t21 = 841.0/64297800.0*t12;
  const double t22 = 2117.0/8037225.0*t15;
  const double t23 = 841.0/64297800.0*t17;
  const double t24 = -t4+t8-t20-t21+t22+t23;
  const double t25 = 841.0/64297800.0*t3;
  const double t26 = 841.0/64297800.0*t7;
  const double t29 = t25-t26+841.0/32148900.0*t10+t21-841.0/32148900.0*t15-t23;
  const double t30 = -t25+t26-t20-t13+t22+t18;
  const double t31 = 1168.0/8037225.0*t3;
  const double t32 = 1168.0/8037225.0*t7;
  const double t33 = 21608.0/8037225.0*t10;
  const double t34 = 1073.0/8037225.0*t12;
  const double t35 = 21608.0/8037225.0*t15;
  const double t36 = 1073.0/8037225.0*t17;
  const double t37 = t31-t32+t33+t34-t35-t36;
  const double t38 = 4234.0/2679075.0*t10;
  const double t39 = 841.0/10716300.0*t12;
  const double t40 = 4234.0/2679075.0*t15;
  const double t41 = 841.0/10716300.0*t17;
  const double t42 = -t38-t39+t40+t41;
  const double t43 = 584.0/1148175.0*t10;
  const double t44 = 29.0/1148175.0*t12;
  const double t45 = 584.0/1148175.0*t15;
  const double t46 = 29.0/1148175.0*t17;
  const double t47 = -t31+t32+t43+t44-t45-t46;
  const double t48 = 1073.0/8037225.0*t3;
  const double t49 = 1073.0/8037225.0*t7;
  const double t50 = 2146.0/8037225.0*t10;
  const double t51 = 116.0/8037225.0*t12;
  const double t52 = 2146.0/8037225.0*t15;
  const double t53 = 116.0/8037225.0*t17;
  const double t54 = -t48+t49-t50-t51+t52+t53;
  const double t55 = 841.0/10716300.0*t3;
  const double t56 = 841.0/10716300.0*t7;
  const double t57 = 841.0/5358150.0*t10;
  const double t58 = 841.0/5358150.0*t15;
  const double t59 = t55-t56+t57-t58;
  const double t60 = 29.0/1148175.0*t3;
  const double t61 = 29.0/1148175.0*t7;
  const double t62 = 58.0/1148175.0*t10;
  const double t63 = 58.0/1148175.0*t15;
  const double t64 = -t60+t61-t62+t51+t63-t53;
  const double t65 = 116.0/8037225.0*t3;
  const double t66 = 116.0/8037225.0*t7;
  const double t67 = -t65+t66-t50-t34+t52+t36;
  const double t68 = t57+t39-t58-t41;
  const double t69 = t65-t66-t62-t44+t63+t46;
  const double t70 = 1168.0/8037225.0*t12;
  const double t71 = 1168.0/8037225.0*t17;
  const double t72 = t48-t49+t33+t70-t35-t71;
  const double t73 = -t55+t56-t38+t40;
  const double t74 = t60-t61+t43-t70-t45+t71;
  const double t75 = 1184.0/8037225.0*t3;
  const double t76 = 1184.0/8037225.0*t7;
  const double t77 = 21904.0/8037225.0*t10;
  const double t78 = 1184.0/8037225.0*t12;
  const double t79 = 21904.0/8037225.0*t15;
  const double t80 = 1184.0/8037225.0*t17;
  const double t81 = t75-t76+t77+t78-t79-t80;
  const double t82 = 592.0/1148175.0*t10;
  const double t83 = 32.0/1148175.0*t12;
  const double t84 = 592.0/1148175.0*t15;
  const double t85 = 32.0/1148175.0*t17;
  const double t86 = -t75+t76+t82+t83-t84-t85;
  const double t87 = 32.0/1148175.0*t3;
  const double t88 = 32.0/1148175.0*t7;
  const double t89 = 16.0/164025.0*t10;
  const double t90 = 16.0/164025.0*t15;
  const double t91 = -t87+t88+t89-t83-t90+t85;
  const double t92 = t87-t88+t82-t78-t84+t80;
  const double t93 = 4292.0/2679075.0*t10;
  const double t94 = 232.0/2679075.0*t12;
  const double t95 = 4292.0/2679075.0*t15;
  const double t96 = 232.0/2679075.0*t17;
  const double t97 = -t93-t94+t95+t96;
  const double t98 = 232.0/2679075.0*t3;
  const double t99 = 232.0/2679075.0*t7;
  const double t100 = 116.0/382725.0*t10;
  const double t101 = 116.0/382725.0*t15;
  const double t102 = t98-t99-t100+t101;
  const double t103 = -t100+t94+t101-t96;
  const double t104 = -t98+t99-t93+t95;
  const double t105 = t10-t15;
  const double t106 = 2701.0/1071630.0*t3;
  const double t107 = 2701.0/1071630.0*t7;
  const double t111 = -1073.0/4286520.0*t3+1073.0/4286520.0*t7-t20-t13+t22+t18;
  const double t112 = 584.0/893025.0*t3;
  const double t113 = 584.0/893025.0*t7;
  const double t114 = t112-t113+t43+t44-t45-t46;
  const double t115 = 4234.0/2679075.0*t3;
  const double t116 = 4234.0/2679075.0*t7;
  const double t117 = -t115+t116-t38-t39+t40+t41;
  const double t118 = 584.0/229635.0*t3;
  const double t119 = 584.0/229635.0*t7;
  const double t120 = t118-t119+t33+t34-t35-t36;
  const double t121 = 1369.0/535815.0*t3;
  const double t122 = 1369.0/535815.0*t7;
  const double t123 = t121-t122+t33+t70-t35-t71;
  const double t124 = 1073.0/714420.0*t3;
  const double t125 = 1073.0/714420.0*t7;
  const double t126 = -t124+t125-t38+t40;
  const double t127 = 37.0/76545.0*t3;
  const double t128 = 37.0/76545.0*t7;
  const double t129 = t127-t128+t43-t70-t45+t71;
  const double t132 = -58.0/893025.0*t3+58.0/893025.0*t7-t62-t44+t63+t46;
  const double t135 = 841.0/5358150.0*t3-841.0/5358150.0*t7+t57+t39-t58-t41;
  const double t138 = -58.0/229635.0*t3+58.0/229635.0*t7-t50-t34+t52+t36;
  const double t139 = 592.0/893025.0*t3;
  const double t140 = 592.0/893025.0*t7;
  const double t141 = t139-t140+t82+t83-t84-t85;
  const double t142 = 592.0/229635.0*t3;
  const double t143 = 592.0/229635.0*t7;
  const double t144 = t142-t143+t77+t78-t79-t80;
  const double t145 = 16.0/32805.0*t3;
  const double t146 = 16.0/32805.0*t7;
  const double t147 = t145-t146+t82-t78-t84+t80;
  const double t148 = 16.0/127575.0*t3;
  const double t149 = 16.0/127575.0*t7;
  const double t150 = t148-t149+t89-t83-t90+t85;
  const double t151 = 4292.0/2679075.0*t3;
  const double t152 = 4292.0/2679075.0*t7;
  const double t153 = -t151+t152-t93-t94+t95+t96;
  const double t154 = 116.0/76545.0*t3;
  const double t155 = 116.0/76545.0*t7;
  const double t156 = -t154+t155-t93+t95;
  const double t157 = 116.0/382725.0*t3;
  const double t158 = 116.0/382725.0*t7;
  const double t159 = -t157+t158-t100+t94+t101-t96;
  const double t160 = 116.0/297675.0*t3;
  const double t161 = 116.0/297675.0*t7;
  const double t162 = -t160+t161-t100+t101;
  const double t163 = t3-t7+t10-t15;
  const double t164 = 2701.0/1071630.0*t12;
  const double t165 = 2701.0/1071630.0*t17;
  const double t169 = -t4+t8-t20-1073.0/4286520.0*t12+t22+1073.0/4286520.0*t17;
  const double t170 = 584.0/893025.0*t12;
  const double t171 = 584.0/893025.0*t17;
  const double t172 = t127-t128+t43+t170-t45-t171;
  const double t173 = 4234.0/2679075.0*t12;
  const double t174 = 4234.0/2679075.0*t17;
  const double t175 = -t124+t125-t38-t173+t40+t174;
  const double t176 = 584.0/229635.0*t12;
  const double t177 = 584.0/229635.0*t17;
  const double t178 = t121-t122+t33+t176-t35-t177;
  const double t179 = 37.0/76545.0*t12;
  const double t180 = 37.0/76545.0*t17;
  const double t181 = t112-t113+t43+t179-t45-t180;
  const double t182 = 1073.0/714420.0*t12;
  const double t183 = 1073.0/714420.0*t17;
  const double t184 = -t115+t116-t38-t182+t40+t183;
  const double t185 = 1369.0/535815.0*t12;
  const double t186 = 1369.0/535815.0*t17;
  const double t187 = t118-t119+t33+t185-t35-t186;
  const double t190 = -t60+t61-t62-58.0/893025.0*t12+t63+58.0/893025.0*t17;
  const double t193 = t55-t56+t57+841.0/5358150.0*t12-t58-841.0/5358150.0*t17;
  const double t196 = -t48+t49-t50-58.0/229635.0*t12+t52+58.0/229635.0*t17;
  const double t197 = 16.0/127575.0*t12;
  const double t198 = 16.0/127575.0*t17;
  const double t199 = t148-t149+t89+t197-t90-t198;
  const double t200 = 592.0/893025.0*t12;
  const double t201 = 592.0/893025.0*t17;
  const double t202 = t145-t146+t82+t200-t84-t201;
  const double t203 = 592.0/229635.0*t12;
  const double t204 = 592.0/229635.0*t17;
  const double t205 = t142-t143+t77+t203-t79-t204;
  const double t206 = 16.0/32805.0*t12;
  const double t207 = 16.0/32805.0*t17;
  const double t208 = t139-t140+t82+t206-t84-t207;
  const double t209 = 116.0/297675.0*t12;
  const double t210 = 116.0/297675.0*t17;
  const double t211 = -t157+t158-t100-t209+t101+t210;
  const double t212 = 4292.0/2679075.0*t12;
  const double t213 = 4292.0/2679075.0*t17;
  const double t214 = -t154+t155-t93-t212+t95+t213;
  const double t215 = 116.0/76545.0*t12;
  const double t216 = 116.0/76545.0*t17;
  const double t217 = -t151+t152-t93-t215+t95+t216;
  const double t218 = 116.0/382725.0*t12;
  const double t219 = 116.0/382725.0*t17;
  const double t220 = -t160+t161-t100-t218+t101+t219;
  const double t221 = t3-t7+t10+t12-t15-t17;
  const double t223 = t31-t32+t33+t185-t35-t186;
  const double t224 = -t38-t182+t40+t183;
  const double t225 = -t31+t32+t43+t179-t45-t180;
  const double t226 = t60-t61+t43+t170-t45-t171;
  const double t227 = -t55+t56-t38-t173+t40+t174;
  const double t228 = t48-t49+t33+t176-t35-t177;
  const double t229 = t87-t88+t82+t200-t84-t201;
  const double t230 = -t87+t88+t89+t197-t90-t198;
  const double t231 = -t75+t76+t82+t206-t84-t207;
  const double t232 = t75-t76+t77+t203-t79-t204;
  const double t233 = -t100-t209+t101+t210;
  const double t234 = t98-t99-t100-t218+t101+t219;
  const double t235 = -t93-t215+t95+t216;
  const double t236 = -t98+t99-t93-t212+t95+t213;
  const double t237 = t10+t12-t15-t17;
  const double t238 = 9344.0/2679075.0*t3;
  const double t239 = 9344.0/2679075.0*t7;
  const double t240 = 18688.0/1148175.0*t10;
  const double t241 = 928.0/1148175.0*t12;
  const double t242 = 18688.0/1148175.0*t15;
  const double t243 = 928.0/1148175.0*t17;
  const double t245 = 2336.0/2679075.0*t3;
  const double t246 = 2336.0/2679075.0*t7;
  const double t247 = 9344.0/2679075.0*t10;
  const double t248 = 464.0/2679075.0*t12;
  const double t249 = 9344.0/2679075.0*t15;
  const double t250 = 464.0/2679075.0*t17;
  const double t251 = -t245+t246-t247-t248+t249+t250;
  const double t252 = 9344.0/8037225.0*t3;
  const double t253 = 9344.0/8037225.0*t7;
  const double t254 = 18688.0/8037225.0*t10;
  const double t255 = 928.0/8037225.0*t12;
  const double t256 = 18688.0/8037225.0*t15;
  const double t257 = 928.0/8037225.0*t17;
  const double t258 = t252-t253+t254+t255-t256-t257;
  const double t261 = 1856.0/1148175.0*t10;
  const double t262 = 1856.0/1148175.0*t15;
  const double t263 = -928.0/2679075.0*t3+928.0/2679075.0*t7-t261-t241+t262+t243;
  const double t264 = 928.0/2679075.0*t10;
  const double t265 = 928.0/2679075.0*t15;
  const double t266 = t98-t99+t264+t248-t265-t250;
  const double t267 = 928.0/8037225.0*t3;
  const double t268 = 928.0/8037225.0*t7;
  const double t271 = -t267+t268-1856.0/8037225.0*t10-t255+1856.0/8037225.0*t15+t257;
  const double t272 = 9472.0/2679075.0*t3;
  const double t273 = 9472.0/2679075.0*t7;
  const double t274 = 18944.0/1148175.0*t10;
  const double t275 = 1024.0/1148175.0*t12;
  const double t276 = 18944.0/1148175.0*t15;
  const double t277 = 1024.0/1148175.0*t17;
  const double t278 = t272-t273+t274+t275-t276-t277;
  const double t279 = 9472.0/8037225.0*t3;
  const double t280 = 9472.0/8037225.0*t7;
  const double t281 = 18944.0/8037225.0*t10;
  const double t282 = 1024.0/8037225.0*t12;
  const double t283 = 18944.0/8037225.0*t15;
  const double t284 = 1024.0/8037225.0*t17;
  const double t285 = t279-t280+t281+t282-t283-t284;
  const double t286 = 256.0/1148175.0*t3;
  const double t287 = 256.0/1148175.0*t7;
  const double t288 = 512.0/1148175.0*t10;
  const double t289 = 512.0/1148175.0*t15;
  const double t290 = t286-t287+t288-t282-t289+t284;
  const double t291 = 256.0/382725.0*t3;
  const double t292 = 256.0/382725.0*t7;
  const double t293 = 512.0/164025.0*t10;
  const double t294 = 512.0/164025.0*t15;
  const double t295 = t291-t292+t293-t275-t294+t277;
  const double t296 = 2368.0/2679075.0*t3;
  const double t297 = 2368.0/2679075.0*t7;
  const double t298 = 9472.0/2679075.0*t10;
  const double t299 = 512.0/2679075.0*t12;
  const double t300 = 9472.0/2679075.0*t15;
  const double t301 = 512.0/2679075.0*t17;
  const double t302 = -t296+t297-t298-t299+t300+t301;
  const double t303 = 1856.0/2679075.0*t3;
  const double t304 = 1856.0/2679075.0*t7;
  const double t305 = 3712.0/2679075.0*t10;
  const double t306 = 3712.0/2679075.0*t15;
  const double t307 = -t303+t304-t305+t306;
  const double t308 = 64.0/382725.0*t3;
  const double t309 = 64.0/382725.0*t7;
  const double t310 = 256.0/382725.0*t10;
  const double t311 = 256.0/382725.0*t15;
  const double t312 = -t308+t309-t310+t299+t311-t301;
  const double t313 = 1856.0/893025.0*t3;
  const double t314 = 1856.0/893025.0*t7;
  const double t315 = 3712.0/382725.0*t10;
  const double t316 = 3712.0/382725.0*t15;
  const double t317 = -t313+t314-t315+t316;
  const double t318 = 464.0/893025.0*t3;
  const double t319 = 464.0/893025.0*t7;
  const double t320 = 1856.0/893025.0*t10;
  const double t321 = 1856.0/893025.0*t15;
  const double t322 = t318-t319+t320-t321;
  const double t323 = 7592.0/893025.0*t3;
  const double t324 = 7592.0/893025.0*t7;
  const double t325 = 15184.0/893025.0*t10;
  const double t326 = 754.0/893025.0*t12;
  const double t327 = 15184.0/893025.0*t15;
  const double t328 = 754.0/893025.0*t17;
  const double t330 = 2336.0/893025.0*t3;
  const double t331 = 2336.0/893025.0*t7;
  const double t332 = -t330+t331-t247-t248+t249+t250;
  const double t333 = 754.0/893025.0*t3;
  const double t334 = 754.0/893025.0*t7;
  const double t337 = -t333+t334-1508.0/893025.0*t10-t326+1508.0/893025.0*t15+t328;
  const double t340 = 232.0/893025.0*t3-232.0/893025.0*t7+t264+t248-t265-t250;
  const double t341 = 2368.0/893025.0*t3;
  const double t342 = 2368.0/893025.0*t7;
  const double t343 = -t341+t342-t298-t299+t300+t301;
  const double t344 = 64.0/127575.0*t3;
  const double t345 = 64.0/127575.0*t7;
  const double t346 = -t344+t345-t310+t299+t311-t301;
  const double t347 = 7696.0/893025.0*t3;
  const double t348 = 7696.0/893025.0*t7;
  const double t349 = 15392.0/893025.0*t10;
  const double t350 = 832.0/893025.0*t12;
  const double t351 = 15392.0/893025.0*t15;
  const double t352 = 832.0/893025.0*t17;
  const double t353 = t347-t348+t349+t350-t351-t352;
  const double t354 = 464.0/297675.0*t3;
  const double t355 = 464.0/297675.0*t7;
  const double t356 = t354-t355+t320-t321;
  const double t357 = 208.0/127575.0*t3;
  const double t358 = 208.0/127575.0*t7;
  const double t359 = 416.0/127575.0*t10;
  const double t360 = 416.0/127575.0*t15;
  const double t361 = t357-t358+t359-t350-t360+t352;
  const double t362 = 1508.0/297675.0*t3;
  const double t363 = 1508.0/297675.0*t7;
  const double t364 = 3016.0/297675.0*t10;
  const double t365 = 3016.0/297675.0*t15;
  const double t366 = -t362+t363-t364+t365;
  const double t367 = 102784.0/8037225.0*t3;
  const double t368 = 102784.0/8037225.0*t7;
  const double t372 = -10208.0/8037225.0*t3+10208.0/8037225.0*t7-t261-t241+t262+t243;
  const double t373 = 104192.0/8037225.0*t3;
  const double t374 = 104192.0/8037225.0*t7;
  const double t375 = t373-t374+t274+t275-t276-t277;
  const double t376 = 2816.0/1148175.0*t3;
  const double t377 = 2816.0/1148175.0*t7;
  const double t378 = t376-t377+t293-t275-t294+t277;
  const double t379 = 20416.0/2679075.0*t3;
  const double t380 = 20416.0/2679075.0*t7;
  const double t381 = -t379+t380-t315+t316;
  const double t382 = 1184.0/76545.0*t3;
  const double t383 = 1184.0/76545.0*t7;
  const double t384 = 9344.0/2679075.0*t12;
  const double t385 = 9344.0/2679075.0*t17;
  const double t387 = 592.0/178605.0*t3;
  const double t388 = 592.0/178605.0*t7;
  const double t389 = 2336.0/2679075.0*t12;
  const double t390 = 2336.0/2679075.0*t17;
  const double t391 = -t387+t388-t247-t389+t249+t390;
  const double t394 = 9344.0/8037225.0*t12;
  const double t395 = 9344.0/8037225.0*t17;
  const double t396 = 1184.0/535815.0*t3-1184.0/535815.0*t7+t254+t394-t256-t395;
  const double t397 = 928.0/1148175.0*t3;
  const double t398 = 928.0/1148175.0*t7;
  const double t401 = -t397+t398-t261-928.0/2679075.0*t12+t262+928.0/2679075.0*t17;
  const double t402 = 464.0/2679075.0*t3;
  const double t403 = 464.0/2679075.0*t7;
  const double t404 = t402-t403+t264+t94-t265-t96;
  const double t405 = 512.0/127575.0*t3;
  const double t406 = 512.0/127575.0*t7;
  const double t407 = 256.0/382725.0*t12;
  const double t408 = 256.0/382725.0*t17;
  const double t409 = t405-t406+t293+t407-t294-t408;
  const double t410 = 512.0/32805.0*t3;
  const double t411 = 512.0/32805.0*t7;
  const double t412 = 9472.0/2679075.0*t12;
  const double t413 = 9472.0/2679075.0*t17;
  const double t414 = t410-t411+t274+t412-t276-t413;
  const double t417 = 9472.0/8037225.0*t12;
  const double t418 = 9472.0/8037225.0*t17;
  const double t419 = 512.0/229635.0*t3-512.0/229635.0*t7+t281+t417-t283-t418;
  const double t422 = 256.0/1148175.0*t12;
  const double t423 = 256.0/1148175.0*t17;
  const double t424 = 512.0/893025.0*t3-512.0/893025.0*t7+t288+t422-t289-t423;
  const double t425 = 3712.0/382725.0*t3;
  const double t426 = 3712.0/382725.0*t7;
  const double t427 = 1856.0/893025.0*t12;
  const double t428 = 1856.0/893025.0*t17;
  const double t429 = -t425+t426-t315-t427+t316+t428;
  const double t430 = 256.0/76545.0*t3;
  const double t431 = 256.0/76545.0*t7;
  const double t432 = 2368.0/2679075.0*t12;
  const double t433 = 2368.0/2679075.0*t17;
  const double t434 = -t430+t431-t298-t432+t300+t433;
  const double t437 = 1856.0/2679075.0*t12;
  const double t438 = 1856.0/2679075.0*t17;
  const double t439 = -3712.0/2679075.0*t3+3712.0/2679075.0*t7-t305-t437+t306+t438;
  const double t440 = 256.0/297675.0*t3;
  const double t441 = 256.0/297675.0*t7;
  const double t442 = 64.0/382725.0*t12;
  const double t443 = 64.0/382725.0*t17;
  const double t444 = -t440+t441-t310-t442+t311+t443;
  const double t445 = 464.0/893025.0*t12;
  const double t446 = 464.0/893025.0*t17;
  const double t447 = t313-t314+t320+t445-t321-t446;
  const double t450 = 7592.0/893025.0*t12;
  const double t451 = 7592.0/893025.0*t17;
  const double t453 = 2336.0/893025.0*t12;
  const double t454 = 2336.0/893025.0*t17;
  const double t455 = -t387+t388-t247-t453+t249+t454;
  const double t458 = t402-t403+t264+232.0/893025.0*t12-t265-232.0/893025.0*t17;
  const double t459 = 2368.0/893025.0*t12;
  const double t460 = 2368.0/893025.0*t17;
  const double t461 = -t430+t431-t298-t459+t300+t460;
  const double t462 = 64.0/127575.0*t12;
  const double t463 = 64.0/127575.0*t17;
  const double t464 = -t440+t441-t310-t462+t311+t463;
  const double t467 = 7696.0/893025.0*t12;
  const double t468 = 7696.0/893025.0*t17;
  const double t469 = 416.0/25515.0*t3-416.0/25515.0*t7+t349+t467-t351-t468;
  const double t470 = 464.0/297675.0*t12;
  const double t471 = 464.0/297675.0*t17;
  const double t472 = t313-t314+t320+t470-t321-t471;
  const double t475 = 208.0/127575.0*t12;
  const double t476 = 208.0/127575.0*t17;
  const double t477 = 416.0/99225.0*t3-416.0/99225.0*t7+t359+t475-t360-t476;
  const double t480 = 1508.0/297675.0*t12;
  const double t481 = 1508.0/297675.0*t17;
  const double t482 = -3016.0/297675.0*t3+3016.0/297675.0*t7-t364-t480+t365+t481;
  const double t483 = 102784.0/8037225.0*t12;
  const double t484 = 102784.0/8037225.0*t17;
  const double t488 = -t397+t398-t261-10208.0/8037225.0*t12+t262+10208.0/8037225.0*t17;
  const double t489 = 104192.0/8037225.0*t12;
  const double t490 = 104192.0/8037225.0*t17;
  const double t491 = t410-t411+t274+t489-t276-t490;
  const double t492 = 2816.0/1148175.0*t12;
  const double t493 = 2816.0/1148175.0*t17;
  const double t494 = t405-t406+t293+t492-t294-t493;
  const double t495 = 20416.0/2679075.0*t12;
  const double t496 = 20416.0/2679075.0*t17;
  const double t497 = -t425+t426-t315-t495+t316+t496;
  const double t498 = 1184.0/76545.0*t12;
  const double t499 = 1184.0/76545.0*t17;
  const double t501 = 592.0/178605.0*t12;
  const double t502 = 592.0/178605.0*t17;
  const double t503 = -t245+t246-t247-t501+t249+t502;
  const double t506 = t252-t253+t254+1184.0/535815.0*t12-t256-1184.0/535815.0*t17;
  const double t507 = 512.0/127575.0*t12;
  const double t508 = 512.0/127575.0*t17;
  const double t509 = t291-t292+t293+t507-t294-t508;
  const double t512 = t286-t287+t288+512.0/893025.0*t12-t289-512.0/893025.0*t17;
  const double t515 = t279-t280+t281+512.0/229635.0*t12-t283-512.0/229635.0*t17;
  const double t516 = 512.0/32805.0*t12;
  const double t517 = 512.0/32805.0*t17;
  const double t518 = t272-t273+t274+t516-t276-t517;
  const double t519 = 256.0/297675.0*t12;
  const double t520 = 256.0/297675.0*t17;
  const double t521 = -t308+t309-t310-t519+t311+t520;
  const double t524 = -t303+t304-t305-3712.0/2679075.0*t12+t306+3712.0/2679075.0*t17;
  const double t525 = 256.0/76545.0*t12;
  const double t526 = 256.0/76545.0*t17;
  const double t527 = -t296+t297-t298-t525+t300+t526;
  const double t528 = 3712.0/382725.0*t12;
  const double t529 = 3712.0/382725.0*t17;
  const double t530 = -t313+t314-t315-t528+t316+t529;
  const double t531 = t318-t319+t320+t427-t321-t428;
  const double t535 = -t330+t331-t247-t501+t249+t502;
  const double t536 = -t344+t345-t310-t519+t311+t520;
  const double t537 = -t341+t342-t298-t525+t300+t526;
  const double t540 = t357-t358+t359+416.0/99225.0*t12-t360-416.0/99225.0*t17;
  const double t541 = t354-t355+t320+t427-t321-t428;
  const double t544 = t347-t348+t349+416.0/25515.0*t12-t351-416.0/25515.0*t17;
  const double t547 = -t362+t363-t364-3016.0/297675.0*t12+t365+3016.0/297675.0*t17;
  const double t549 = t376-t377+t293+t507-t294-t508;
  const double t550 = t373-t374+t274+t516-t276-t517;
  const double t551 = -t379+t380-t315-t528+t316+t529;
  const double t553 = -t402+t403-t247-t389+t249+t390;
  const double t554 = t267-t268+t254+t394-t256-t395;
  const double t555 = 1024.0/1148175.0*t3;
  const double t556 = 1024.0/1148175.0*t7;
  const double t557 = t555-t556+t274+t412-t276-t413;
  const double t558 = -t555+t556+t293+t407-t294-t408;
  const double t559 = 1024.0/8037225.0*t3;
  const double t560 = 1024.0/8037225.0*t7;
  const double t561 = -t559+t560+t288+t422-t289-t423;
  const double t562 = t559-t560+t281+t417-t283-t418;
  const double t563 = -t315-t427+t316+t428;
  const double t564 = 512.0/2679075.0*t3;
  const double t565 = 512.0/2679075.0*t7;
  const double t566 = t564-t565-t310-t442+t311+t443;
  const double t567 = -t305-t437+t306+t438;
  const double t568 = -t564+t565-t298-t432+t300+t433;
  const double t569 = t320+t445-t321-t446;
  const double t571 = -t402+t403-t247-t453+t249+t454;
  const double t572 = t564-t565-t310-t462+t311+t463;
  const double t573 = -t564+t565-t298-t459+t300+t460;
  const double t574 = 832.0/893025.0*t3;
  const double t575 = 832.0/893025.0*t7;
  const double t576 = -t574+t575+t359+t475-t360-t476;
  const double t577 = t320+t470-t321-t471;
  const double t578 = t574-t575+t349+t467-t351-t468;
  const double t579 = -t364-t480+t365+t481;
  const double t581 = -t555+t556+t293+t492-t294-t493;
  const double t582 = t555-t556+t274+t489-t276-t490;
  const double t583 = -t315-t495+t316+t496;
  const double t584 = 8192.0/382725.0*t3;
  const double t585 = 8192.0/382725.0*t7;
  const double t586 = 16384.0/164025.0*t10;
  const double t587 = 8192.0/382725.0*t12;
  const double t588 = 16384.0/164025.0*t15;
  const double t589 = 8192.0/382725.0*t17;
  const double t591 = 8192.0/1148175.0*t3;
  const double t592 = 8192.0/1148175.0*t7;
  const double t593 = 16384.0/1148175.0*t10;
  const double t595 = 16384.0/1148175.0*t15;
  const double t597 = t591-t592+t593+8192.0/2679075.0*t12-t595-8192.0/2679075.0*t17;
  const double t604 = 8192.0/8037225.0*t3-8192.0/8037225.0*t7+16384.0/8037225.0*t10+
		      8192.0/8037225.0*t12-16384.0/8037225.0*t15-8192.0/8037225.0*t17;
  const double t607 = 8192.0/1148175.0*t12;
  const double t608 = 8192.0/1148175.0*t17;
  const double t609 = 8192.0/2679075.0*t3-8192.0/2679075.0*t7+t593+t607-t595-t608;
  const double t610 = 2048.0/382725.0*t3;
  const double t611 = 2048.0/382725.0*t7;
  const double t612 = 8192.0/382725.0*t10;
  const double t613 = 4096.0/893025.0*t12;
  const double t614 = 8192.0/382725.0*t15;
  const double t615 = 4096.0/893025.0*t17;
  const double t616 = -t610+t611-t612-t613+t614+t615;
  const double t617 = 4096.0/2679075.0*t3;
  const double t618 = 4096.0/2679075.0*t7;
  const double t619 = 8192.0/2679075.0*t10;
  const double t621 = 8192.0/2679075.0*t15;
  const double t623 = -t617+t618-t619-2048.0/2679075.0*t12+t621+2048.0/2679075.0*t17;
  const double t626 = 4096.0/2679075.0*t12;
  const double t627 = 4096.0/2679075.0*t17;
  const double t628 = -2048.0/2679075.0*t3+2048.0/2679075.0*t7-t619-t626+t621+t627;
  const double t629 = 4096.0/893025.0*t3;
  const double t630 = 4096.0/893025.0*t7;
  const double t631 = 2048.0/382725.0*t12;
  const double t632 = 2048.0/382725.0*t17;
  const double t633 = -t629+t630-t612-t631+t614+t632;
  const double t634 = 1024.0/893025.0*t3;
  const double t635 = 1024.0/893025.0*t7;
  const double t636 = 4096.0/893025.0*t10;
  const double t637 = 1024.0/893025.0*t12;
  const double t638 = 4096.0/893025.0*t15;
  const double t639 = 1024.0/893025.0*t17;
  const double t640 = t634-t635+t636+t637-t638-t639;
  const double t641 = 90112.0/1148175.0*t3;
  const double t642 = 90112.0/1148175.0*t7;
  const double t646 = 90112.0/8037225.0*t3-90112.0/8037225.0*t7+t593+t607-t595-t608;
  const double t647 = 2048.0/127575.0*t3;
  const double t648 = 2048.0/127575.0*t7;
  const double t649 = -t647+t648-t612-t613+t614+t615;
  const double t650 = 45056.0/2679075.0*t3;
  const double t651 = 45056.0/2679075.0*t7;
  const double t652 = -t650+t651-t612-t631+t614+t632;
  const double t655 = -2048.0/893025.0*t3+2048.0/893025.0*t7-t619-t626+t621+t627;
  const double t656 = 1024.0/297675.0*t3;
  const double t657 = 1024.0/297675.0*t7;
  const double t658 = t656-t657+t636+t637-t638-t639;
  const double t659 = 90112.0/1148175.0*t12;
  const double t660 = 90112.0/1148175.0*t17;
  const double t664 = t591-t592+t593+90112.0/8037225.0*t12-t595-90112.0/8037225.0*t17;
  const double t665 = 2048.0/127575.0*t12;
  const double t666 = 2048.0/127575.0*t17;
  const double t667 = -t650+t651-t612-t665+t614+t666;
  const double t668 = 45056.0/2679075.0*t12;
  const double t669 = 45056.0/2679075.0*t17;
  const double t670 = -t647+t648-t612-t668+t614+t669;
  const double t673 = -t617+t618-t619-2048.0/893025.0*t12+t621+2048.0/893025.0*t17;
  const double t674 = 1024.0/297675.0*t12;
  const double t675 = 1024.0/297675.0*t17;
  const double t676 = t656-t657+t636+t674-t638-t675;
  const double t678 = -t610+t611-t612-t668+t614+t669;
  const double t679 = -t629+t630-t612-t665+t614+t666;
  const double t680 = t634-t635+t636+t674-t638-t675;
  const double t681 = 6656.0/127575.0*t3;
  const double t682 = 6656.0/127575.0*t7;
  const double t683 = 13312.0/127575.0*t10;
  const double t685 = 13312.0/127575.0*t15;
  const double t694 = 6656.0/893025.0*t3-6656.0/893025.0*t7+13312.0/893025.0*t10+6656.0/
		      893025.0*t12-13312.0/893025.0*t15-6656.0/893025.0*t17;
  const double t695 = 3328.0/297675.0*t3;
  const double t696 = 3328.0/297675.0*t7;
  const double t697 = 6656.0/297675.0*t10;
  const double t699 = 6656.0/297675.0*t15;
  const double t701 = -t695+t696-t697-1664.0/297675.0*t12+t699+1664.0/297675.0*t17;
  const double t704 = 6656.0/127575.0*t12;
  const double t705 = 6656.0/127575.0*t17;
  const double t709 = 3328.0/297675.0*t12;
  const double t710 = 3328.0/297675.0*t17;
  const double t711 = -1664.0/99225.0*t3+1664.0/99225.0*t7-t697-t709+t699+t710;
  const double t717 = -t695+t696-t697-1664.0/99225.0*t12+t699+1664.0/99225.0*t17;
  const double t723 = -1664.0/297675.0*t3+1664.0/297675.0*t7-t697-t709+t699+t710;
  local_mass_matrix(0,0) = t4-t8+t11+t13-t16-t18;
  local_mass_matrix(0,1) = t24;
  local_mass_matrix(0,2) = t29;
  local_mass_matrix(0,3) = t30;
  local_mass_matrix(0,4) = t37;
  local_mass_matrix(0,5) = t42;
  local_mass_matrix(0,6) = t47;
  local_mass_matrix(0,7) = t54;
  local_mass_matrix(0,8) = t59;
  local_mass_matrix(0,9) = t64;
  local_mass_matrix(0,10) = t67;
  local_mass_matrix(0,11) = t68;
  local_mass_matrix(0,12) = t69;
  local_mass_matrix(0,13) = t72;
  local_mass_matrix(0,14) = t73;
  local_mass_matrix(0,15) = t74;
  local_mass_matrix(0,16) = t81;
  local_mass_matrix(0,17) = t86;
  local_mass_matrix(0,18) = t91;
  local_mass_matrix(0,19) = t92;
  local_mass_matrix(0,20) = t97;
  local_mass_matrix(0,21) = t102;
  local_mass_matrix(0,22) = t103;
  local_mass_matrix(0,23) = t104;
  local_mass_matrix(0,24) = 841.0/893025.0*t105;
  local_mass_matrix(1,0) = t24;
  local_mass_matrix(1,1) = t106-t107+t11+t13-t16-t18;
  local_mass_matrix(1,2) = t111;
  local_mass_matrix(1,3) = t29;
  local_mass_matrix(1,4) = t114;
  local_mass_matrix(1,5) = t117;
  local_mass_matrix(1,6) = t120;
  local_mass_matrix(1,7) = t123;
  local_mass_matrix(1,8) = t126;
  local_mass_matrix(1,9) = t129;
  local_mass_matrix(1,10) = t132;
  local_mass_matrix(1,11) = t135;
  local_mass_matrix(1,12) = t138;
  local_mass_matrix(1,13) = t54;
  local_mass_matrix(1,14) = t59;
  local_mass_matrix(1,15) = t64;
  local_mass_matrix(1,16) = t141;
  local_mass_matrix(1,17) = t144;
  local_mass_matrix(1,18) = t147;
  local_mass_matrix(1,19) = t150;
  local_mass_matrix(1,20) = t153;
  local_mass_matrix(1,21) = t156;
  local_mass_matrix(1,22) = t159;
  local_mass_matrix(1,23) = t162;
  local_mass_matrix(1,24) = 841.0/893025.0*t163;
  local_mass_matrix(2,0) = t29;
  local_mass_matrix(2,1) = t111;
  local_mass_matrix(2,2) = t106-t107+t11+t164-t16-t165;
  local_mass_matrix(2,3) = t169;
  local_mass_matrix(2,4) = t132;
  local_mass_matrix(2,5) = t135;
  local_mass_matrix(2,6) = t138;
  local_mass_matrix(2,7) = t172;
  local_mass_matrix(2,8) = t175;
  local_mass_matrix(2,9) = t178;
  local_mass_matrix(2,10) = t181;
  local_mass_matrix(2,11) = t184;
  local_mass_matrix(2,12) = t187;
  local_mass_matrix(2,13) = t190;
  local_mass_matrix(2,14) = t193;
  local_mass_matrix(2,15) = t196;
  local_mass_matrix(2,16) = t199;
  local_mass_matrix(2,17) = t202;
  local_mass_matrix(2,18) = t205;
  local_mass_matrix(2,19) = t208;
  local_mass_matrix(2,20) = t211;
  local_mass_matrix(2,21) = t214;
  local_mass_matrix(2,22) = t217;
  local_mass_matrix(2,23) = t220;
  local_mass_matrix(2,24) = 841.0/893025.0*t221;
  local_mass_matrix(3,0) = t30;
  local_mass_matrix(3,1) = t29;
  local_mass_matrix(3,2) = t169;
  local_mass_matrix(3,3) = t4-t8+t11+t164-t16-t165;
  local_mass_matrix(3,4) = t67;
  local_mass_matrix(3,5) = t68;
  local_mass_matrix(3,6) = t69;
  local_mass_matrix(3,7) = t190;
  local_mass_matrix(3,8) = t193;
  local_mass_matrix(3,9) = t196;
  local_mass_matrix(3,10) = t223;
  local_mass_matrix(3,11) = t224;
  local_mass_matrix(3,12) = t225;
  local_mass_matrix(3,13) = t226;
  local_mass_matrix(3,14) = t227;
  local_mass_matrix(3,15) = t228;
  local_mass_matrix(3,16) = t229;
  local_mass_matrix(3,17) = t230;
  local_mass_matrix(3,18) = t231;
  local_mass_matrix(3,19) = t232;
  local_mass_matrix(3,20) = t233;
  local_mass_matrix(3,21) = t234;
  local_mass_matrix(3,22) = t235;
  local_mass_matrix(3,23) = t236;
  local_mass_matrix(3,24) = 841.0/893025.0*t237;
  local_mass_matrix(4,0) = t37;
  local_mass_matrix(4,1) = t114;
  local_mass_matrix(4,2) = t132;
  local_mass_matrix(4,3) = t67;
  local_mass_matrix(4,4) = t238-t239+t240+t241-t242-t243;
  local_mass_matrix(4,5) = t251;
  local_mass_matrix(4,6) = t258;
  local_mass_matrix(4,7) = t141;
  local_mass_matrix(4,8) = t162;
  local_mass_matrix(4,9) = t150;
  local_mass_matrix(4,10) = t263;
  local_mass_matrix(4,11) = t266;
  local_mass_matrix(4,12) = t271;
  local_mass_matrix(4,13) = t81;
  local_mass_matrix(4,14) = t104;
  local_mass_matrix(4,15) = t92;
  local_mass_matrix(4,16) = t278;
  local_mass_matrix(4,17) = t285;
  local_mass_matrix(4,18) = t290;
  local_mass_matrix(4,19) = t295;
  local_mass_matrix(4,20) = t302;
  local_mass_matrix(4,21) = t307;
  local_mass_matrix(4,22) = t312;
  local_mass_matrix(4,23) = t317;
  local_mass_matrix(4,24) = t322;
  local_mass_matrix(5,0) = t42;
  local_mass_matrix(5,1) = t117;
  local_mass_matrix(5,2) = t135;
  local_mass_matrix(5,3) = t68;
  local_mass_matrix(5,4) = t251;
  local_mass_matrix(5,5) = t323-t324+t325+t326-t327-t328;
  local_mass_matrix(5,6) = t332;
  local_mass_matrix(5,7) = t153;
  local_mass_matrix(5,8) = 841.0/893025.0*t163;
  local_mass_matrix(5,9) = t159;
  local_mass_matrix(5,10) = t266;
  local_mass_matrix(5,11) = t337;
  local_mass_matrix(5,12) = t340;
  local_mass_matrix(5,13) = t97;
  local_mass_matrix(5,14) = 841.0/893025.0*t105;
  local_mass_matrix(5,15) = t103;
  local_mass_matrix(5,16) = t302;
  local_mass_matrix(5,17) = t343;
  local_mass_matrix(5,18) = t346;
  local_mass_matrix(5,19) = t312;
  local_mass_matrix(5,20) = t353;
  local_mass_matrix(5,21) = t356;
  local_mass_matrix(5,22) = t361;
  local_mass_matrix(5,23) = t322;
  local_mass_matrix(5,24) = t366;
  local_mass_matrix(6,0) = t47;
  local_mass_matrix(6,1) = t120;
  local_mass_matrix(6,2) = t138;
  local_mass_matrix(6,3) = t69;
  local_mass_matrix(6,4) = t258;
  local_mass_matrix(6,5) = t332;
  local_mass_matrix(6,6) = t367-t368+t240+t241-t242-t243;
  local_mass_matrix(6,7) = t144;
  local_mass_matrix(6,8) = t156;
  local_mass_matrix(6,9) = t147;
  local_mass_matrix(6,10) = t271;
  local_mass_matrix(6,11) = t340;
  local_mass_matrix(6,12) = t372;
  local_mass_matrix(6,13) = t86;
  local_mass_matrix(6,14) = t102;
  local_mass_matrix(6,15) = t91;
  local_mass_matrix(6,16) = t285;
  local_mass_matrix(6,17) = t375;
  local_mass_matrix(6,18) = t378;
  local_mass_matrix(6,19) = t290;
  local_mass_matrix(6,20) = t343;
  local_mass_matrix(6,21) = t381;
  local_mass_matrix(6,22) = t346;
  local_mass_matrix(6,23) = t307;
  local_mass_matrix(6,24) = t356;
  local_mass_matrix(7,0) = t54;
  local_mass_matrix(7,1) = t123;
  local_mass_matrix(7,2) = t172;
  local_mass_matrix(7,3) = t190;
  local_mass_matrix(7,4) = t141;
  local_mass_matrix(7,5) = t153;
  local_mass_matrix(7,6) = t144;
  local_mass_matrix(7,7) = t382-t383+t240+t384-t242-t385;
  local_mass_matrix(7,8) = t391;
  local_mass_matrix(7,9) = t396;
  local_mass_matrix(7,10) = t199;
  local_mass_matrix(7,11) = t211;
  local_mass_matrix(7,12) = t202;
  local_mass_matrix(7,13) = t401;
  local_mass_matrix(7,14) = t404;
  local_mass_matrix(7,15) = t271;
  local_mass_matrix(7,16) = t409;
  local_mass_matrix(7,17) = t414;
  local_mass_matrix(7,18) = t419;
  local_mass_matrix(7,19) = t424;
  local_mass_matrix(7,20) = t429;
  local_mass_matrix(7,21) = t434;
  local_mass_matrix(7,22) = t439;
  local_mass_matrix(7,23) = t444;
  local_mass_matrix(7,24) = t447;
  local_mass_matrix(8,0) = t59;
  local_mass_matrix(8,1) = t126;
  local_mass_matrix(8,2) = t175;
  local_mass_matrix(8,3) = t193;
  local_mass_matrix(8,4) = t162;
  local_mass_matrix(8,5) = 841.0/893025.0*t163;
  local_mass_matrix(8,6) = t156;
  local_mass_matrix(8,7) = t391;
  local_mass_matrix(8,8) = 962.0/59535.0*t3-962.0/59535.0*t7+t325+t450-
			   t327-t451;
  local_mass_matrix(8,9) = t455;
  local_mass_matrix(8,10) = t220;
  local_mass_matrix(8,11) = 841.0/893025.0*t221;
  local_mass_matrix(8,12) = t214;
  local_mass_matrix(8,13) = t404;
  local_mass_matrix(8,14) = t337;
  local_mass_matrix(8,15) = t458;
  local_mass_matrix(8,16) = t444;
  local_mass_matrix(8,17) = t434;
  local_mass_matrix(8,18) = t461;
  local_mass_matrix(8,19) = t464;
  local_mass_matrix(8,20) = t447;
  local_mass_matrix(8,21) = t469;
  local_mass_matrix(8,22) = t472;
  local_mass_matrix(8,23) = t477;
  local_mass_matrix(8,24) = t482;
  local_mass_matrix(9,0) = t64;
  local_mass_matrix(9,1) = t129;
  local_mass_matrix(9,2) = t178;
  local_mass_matrix(9,3) = t196;
  local_mass_matrix(9,4) = t150;
  local_mass_matrix(9,5) = t159;
  local_mass_matrix(9,6) = t147;
  local_mass_matrix(9,7) = t396;
  local_mass_matrix(9,8) = t455;
  local_mass_matrix(9,9) = t382-t383+t240+t483-t242-t484;
  local_mass_matrix(9,10) = t208;
  local_mass_matrix(9,11) = t217;
  local_mass_matrix(9,12) = t205;
  local_mass_matrix(9,13) = t271;
  local_mass_matrix(9,14) = t458;
  local_mass_matrix(9,15) = t488;
  local_mass_matrix(9,16) = t424;
  local_mass_matrix(9,17) = t419;
  local_mass_matrix(9,18) = t491;
  local_mass_matrix(9,19) = t494;
  local_mass_matrix(9,20) = t439;
  local_mass_matrix(9,21) = t461;
  local_mass_matrix(9,22) = t497;
  local_mass_matrix(9,23) = t464;
  local_mass_matrix(9,24) = t472;
  local_mass_matrix(10,0) = t67;
  local_mass_matrix(10,1) = t132;
  local_mass_matrix(10,2) = t181;
  local_mass_matrix(10,3) = t223;
  local_mass_matrix(10,4) = t263;
  local_mass_matrix(10,5) = t266;
  local_mass_matrix(10,6) = t271;
  local_mass_matrix(10,7) = t199;
  local_mass_matrix(10,8) = t220;
  local_mass_matrix(10,9) = t208;
  local_mass_matrix(10,10) = t238-t239+t240+t498-t242-t499;
  local_mass_matrix(10,11) = t503;
  local_mass_matrix(10,12) = t506;
  local_mass_matrix(10,13) = t229;
  local_mass_matrix(10,14) = t236;
  local_mass_matrix(10,15) = t232;
  local_mass_matrix(10,16) = t509;
  local_mass_matrix(10,17) = t512;
  local_mass_matrix(10,18) = t515;
  local_mass_matrix(10,19) = t518;
  local_mass_matrix(10,20) = t521;
  local_mass_matrix(10,21) = t524;
  local_mass_matrix(10,22) = t527;
  local_mass_matrix(10,23) = t530;
  local_mass_matrix(10,24) = t531;
  local_mass_matrix(11,0) = t68;
  local_mass_matrix(11,1) = t135;
  local_mass_matrix(11,2) = t184;
  local_mass_matrix(11,3) = t224;
  local_mass_matrix(11,4) = t266;
  local_mass_matrix(11,5) = t337;
  local_mass_matrix(11,6) = t340;
  local_mass_matrix(11,7) = t211;
  local_mass_matrix(11,8) = 841.0/893025.0*t221;
  local_mass_matrix(11,9) = t217;
  local_mass_matrix(11,10) = t503;
  local_mass_matrix(11,11) = t323-t324+t325+962.0/59535.0*t12-t327-962.0/
			     59535.0*t17;
  local_mass_matrix(11,12) = t535;
  local_mass_matrix(11,13) = t233;
  local_mass_matrix(11,14) = 841.0/893025.0*t237;
  local_mass_matrix(11,15) = t235;
  local_mass_matrix(11,16) = t521;
  local_mass_matrix(11,17) = t536;
  local_mass_matrix(11,18) = t537;
  local_mass_matrix(11,19) = t527;
  local_mass_matrix(11,20) = t540;
  local_mass_matrix(11,21) = t541;
  local_mass_matrix(11,22) = t544;
  local_mass_matrix(11,23) = t531;
  local_mass_matrix(11,24) = t547;
  local_mass_matrix(12,0) = t69;
  local_mass_matrix(12,1) = t138;
  local_mass_matrix(12,2) = t187;
  local_mass_matrix(12,3) = t225;
  local_mass_matrix(12,4) = t271;
  local_mass_matrix(12,5) = t340;
  local_mass_matrix(12,6) = t372;
  local_mass_matrix(12,7) = t202;
  local_mass_matrix(12,8) = t214;
  local_mass_matrix(12,9) = t205;
  local_mass_matrix(12,10) = t506;
  local_mass_matrix(12,11) = t535;
  local_mass_matrix(12,12) = t367-t368+t240+t498-t242-t499;
  local_mass_matrix(12,13) = t230;
  local_mass_matrix(12,14) = t234;
  local_mass_matrix(12,15) = t231;
  local_mass_matrix(12,16) = t512;
  local_mass_matrix(12,17) = t549;
  local_mass_matrix(12,18) = t550;
  local_mass_matrix(12,19) = t515;
  local_mass_matrix(12,20) = t536;
  local_mass_matrix(12,21) = t551;
  local_mass_matrix(12,22) = t537;
  local_mass_matrix(12,23) = t524;
  local_mass_matrix(12,24) = t541;
  local_mass_matrix(13,0) = t72;
  local_mass_matrix(13,1) = t54;
  local_mass_matrix(13,2) = t190;
  local_mass_matrix(13,3) = t226;
  local_mass_matrix(13,4) = t81;
  local_mass_matrix(13,5) = t97;
  local_mass_matrix(13,6) = t86;
  local_mass_matrix(13,7) = t401;
  local_mass_matrix(13,8) = t404;
  local_mass_matrix(13,9) = t271;
  local_mass_matrix(13,10) = t229;
  local_mass_matrix(13,11) = t233;
  local_mass_matrix(13,12) = t230;
  local_mass_matrix(13,13) = t397-t398+t240+t384-t242-t385;
  local_mass_matrix(13,14) = t553;
  local_mass_matrix(13,15) = t554;
  local_mass_matrix(13,16) = t557;
  local_mass_matrix(13,17) = t558;
  local_mass_matrix(13,18) = t561;
  local_mass_matrix(13,19) = t562;
  local_mass_matrix(13,20) = t563;
  local_mass_matrix(13,21) = t566;
  local_mass_matrix(13,22) = t567;
  local_mass_matrix(13,23) = t568;
  local_mass_matrix(13,24) = t569;
  local_mass_matrix(14,0) = t73;
  local_mass_matrix(14,1) = t59;
  local_mass_matrix(14,2) = t193;
  local_mass_matrix(14,3) = t227;
  local_mass_matrix(14,4) = t104;
  local_mass_matrix(14,5) = 841.0/893025.0*t105;
  local_mass_matrix(14,6) = t102;
  local_mass_matrix(14,7) = t404;
  local_mass_matrix(14,8) = t337;
  local_mass_matrix(14,9) = t458;
  local_mass_matrix(14,10) = t236;
  local_mass_matrix(14,11) = 841.0/893025.0*t237;
  local_mass_matrix(14,12) = t234;
  local_mass_matrix(14,13) = t553;
  local_mass_matrix(14,14) = t333-t334+t325+t450-t327-t451;
  local_mass_matrix(14,15) = t571;
  local_mass_matrix(14,16) = t568;
  local_mass_matrix(14,17) = t566;
  local_mass_matrix(14,18) = t572;
  local_mass_matrix(14,19) = t573;
  local_mass_matrix(14,20) = t569;
  local_mass_matrix(14,21) = t576;
  local_mass_matrix(14,22) = t577;
  local_mass_matrix(14,23) = t578;
  local_mass_matrix(14,24) = t579;
  local_mass_matrix(15,0) = t74;
  local_mass_matrix(15,1) = t64;
  local_mass_matrix(15,2) = t196;
  local_mass_matrix(15,3) = t228;
  local_mass_matrix(15,4) = t92;
  local_mass_matrix(15,5) = t103;
  local_mass_matrix(15,6) = t91;
  local_mass_matrix(15,7) = t271;
  local_mass_matrix(15,8) = t458;
  local_mass_matrix(15,9) = t488;
  local_mass_matrix(15,10) = t232;
  local_mass_matrix(15,11) = t235;
  local_mass_matrix(15,12) = t231;
  local_mass_matrix(15,13) = t554;
  local_mass_matrix(15,14) = t571;
  local_mass_matrix(15,15) = t397-t398+t240+t483-t242-t484;
  local_mass_matrix(15,16) = t562;
  local_mass_matrix(15,17) = t561;
  local_mass_matrix(15,18) = t581;
  local_mass_matrix(15,19) = t582;
  local_mass_matrix(15,20) = t567;
  local_mass_matrix(15,21) = t572;
  local_mass_matrix(15,22) = t583;
  local_mass_matrix(15,23) = t573;
  local_mass_matrix(15,24) = t577;
  local_mass_matrix(16,0) = t81;
  local_mass_matrix(16,1) = t141;
  local_mass_matrix(16,2) = t199;
  local_mass_matrix(16,3) = t229;
  local_mass_matrix(16,4) = t278;
  local_mass_matrix(16,5) = t302;
  local_mass_matrix(16,6) = t285;
  local_mass_matrix(16,7) = t409;
  local_mass_matrix(16,8) = t444;
  local_mass_matrix(16,9) = t424;
  local_mass_matrix(16,10) = t509;
  local_mass_matrix(16,11) = t521;
  local_mass_matrix(16,12) = t512;
  local_mass_matrix(16,13) = t557;
  local_mass_matrix(16,14) = t568;
  local_mass_matrix(16,15) = t562;
  local_mass_matrix(16,16) = t584-t585+t586+t587-t588-t589;
  local_mass_matrix(16,17) = t597;
  local_mass_matrix(16,18) = t604;
  local_mass_matrix(16,19) = t609;
  local_mass_matrix(16,20) = t616;
  local_mass_matrix(16,21) = t623;
  local_mass_matrix(16,22) = t628;
  local_mass_matrix(16,23) = t633;
  local_mass_matrix(16,24) = t640;
  local_mass_matrix(17,0) = t86;
  local_mass_matrix(17,1) = t144;
  local_mass_matrix(17,2) = t202;
  local_mass_matrix(17,3) = t230;
  local_mass_matrix(17,4) = t285;
  local_mass_matrix(17,5) = t343;
  local_mass_matrix(17,6) = t375;
  local_mass_matrix(17,7) = t414;
  local_mass_matrix(17,8) = t434;
  local_mass_matrix(17,9) = t419;
  local_mass_matrix(17,10) = t512;
  local_mass_matrix(17,11) = t536;
  local_mass_matrix(17,12) = t549;
  local_mass_matrix(17,13) = t558;
  local_mass_matrix(17,14) = t566;
  local_mass_matrix(17,15) = t561;
  local_mass_matrix(17,16) = t597;
  local_mass_matrix(17,17) = t641-t642+t586+t587-t588-t589;
  local_mass_matrix(17,18) = t646;
  local_mass_matrix(17,19) = t604;
  local_mass_matrix(17,20) = t649;
  local_mass_matrix(17,21) = t652;
  local_mass_matrix(17,22) = t655;
  local_mass_matrix(17,23) = t623;
  local_mass_matrix(17,24) = t658;
  local_mass_matrix(18,0) = t91;
  local_mass_matrix(18,1) = t147;
  local_mass_matrix(18,2) = t205;
  local_mass_matrix(18,3) = t231;
  local_mass_matrix(18,4) = t290;
  local_mass_matrix(18,5) = t346;
  local_mass_matrix(18,6) = t378;
  local_mass_matrix(18,7) = t419;
  local_mass_matrix(18,8) = t461;
  local_mass_matrix(18,9) = t491;
  local_mass_matrix(18,10) = t515;
  local_mass_matrix(18,11) = t537;
  local_mass_matrix(18,12) = t550;
  local_mass_matrix(18,13) = t561;
  local_mass_matrix(18,14) = t572;
  local_mass_matrix(18,15) = t581;
  local_mass_matrix(18,16) = t604;
  local_mass_matrix(18,17) = t646;
  local_mass_matrix(18,18) = t641-t642+t586+t659-t588-t660;
  local_mass_matrix(18,19) = t664;
  local_mass_matrix(18,20) = t655;
  local_mass_matrix(18,21) = t667;
  local_mass_matrix(18,22) = t670;
  local_mass_matrix(18,23) = t673;
  local_mass_matrix(18,24) = t676;
  local_mass_matrix(19,0) = t92;
  local_mass_matrix(19,1) = t150;
  local_mass_matrix(19,2) = t208;
  local_mass_matrix(19,3) = t232;
  local_mass_matrix(19,4) = t295;
  local_mass_matrix(19,5) = t312;
  local_mass_matrix(19,6) = t290;
  local_mass_matrix(19,7) = t424;
  local_mass_matrix(19,8) = t464;
  local_mass_matrix(19,9) = t494;
  local_mass_matrix(19,10) = t518;
  local_mass_matrix(19,11) = t527;
  local_mass_matrix(19,12) = t515;
  local_mass_matrix(19,13) = t562;
  local_mass_matrix(19,14) = t573;
  local_mass_matrix(19,15) = t582;
  local_mass_matrix(19,16) = t609;
  local_mass_matrix(19,17) = t604;
  local_mass_matrix(19,18) = t664;
  local_mass_matrix(19,19) = t584-t585+t586+t659-t588-t660;
  local_mass_matrix(19,20) = t628;
  local_mass_matrix(19,21) = t673;
  local_mass_matrix(19,22) = t678;
  local_mass_matrix(19,23) = t679;
  local_mass_matrix(19,24) = t680;
  local_mass_matrix(20,0) = t97;
  local_mass_matrix(20,1) = t153;
  local_mass_matrix(20,2) = t211;
  local_mass_matrix(20,3) = t233;
  local_mass_matrix(20,4) = t302;
  local_mass_matrix(20,5) = t353;
  local_mass_matrix(20,6) = t343;
  local_mass_matrix(20,7) = t429;
  local_mass_matrix(20,8) = t447;
  local_mass_matrix(20,9) = t439;
  local_mass_matrix(20,10) = t521;
  local_mass_matrix(20,11) = t540;
  local_mass_matrix(20,12) = t536;
  local_mass_matrix(20,13) = t563;
  local_mass_matrix(20,14) = t569;
  local_mass_matrix(20,15) = t567;
  local_mass_matrix(20,16) = t616;
  local_mass_matrix(20,17) = t649;
  local_mass_matrix(20,18) = t655;
  local_mass_matrix(20,19) = t628;
  local_mass_matrix(20,20) = t681-t682+t683+6656.0/297675.0*t12-t685
			     -6656.0/297675.0*t17;
  local_mass_matrix(20,21) = t658;
  local_mass_matrix(20,22) = t694;
  local_mass_matrix(20,23) = t640;
  local_mass_matrix(20,24) = t701;
  local_mass_matrix(21,0) = t102;
  local_mass_matrix(21,1) = t156;
  local_mass_matrix(21,2) = t214;
  local_mass_matrix(21,3) = t234;
  local_mass_matrix(21,4) = t307;
  local_mass_matrix(21,5) = t356;
  local_mass_matrix(21,6) = t381;
  local_mass_matrix(21,7) = t434;
  local_mass_matrix(21,8) = t469;
  local_mass_matrix(21,9) = t461;
  local_mass_matrix(21,10) = t524;
  local_mass_matrix(21,11) = t541;
  local_mass_matrix(21,12) = t551;
  local_mass_matrix(21,13) = t566;
  local_mass_matrix(21,14) = t576;
  local_mass_matrix(21,15) = t572;
  local_mass_matrix(21,16) = t623;
  local_mass_matrix(21,17) = t652;
  local_mass_matrix(21,18) = t667;
  local_mass_matrix(21,19) = t673;
  local_mass_matrix(21,20) = t658;
  local_mass_matrix(21,21) = 73216.0/893025.0*t3-73216.0/893025.0*t7+t683+
			     t704-t685-t705;
  local_mass_matrix(21,22) = t676;
  local_mass_matrix(21,23) = t694;
  local_mass_matrix(21,24) = t711;
  local_mass_matrix(22,0) = t103;
  local_mass_matrix(22,1) = t159;
  local_mass_matrix(22,2) = t217;
  local_mass_matrix(22,3) = t235;
  local_mass_matrix(22,4) = t312;
  local_mass_matrix(22,5) = t361;
  local_mass_matrix(22,6) = t346;
  local_mass_matrix(22,7) = t439;
  local_mass_matrix(22,8) = t472;
  local_mass_matrix(22,9) = t497;
  local_mass_matrix(22,10) = t527;
  local_mass_matrix(22,11) = t544;
  local_mass_matrix(22,12) = t537;
  local_mass_matrix(22,13) = t567;
  local_mass_matrix(22,14) = t577;
  local_mass_matrix(22,15) = t583;
  local_mass_matrix(22,16) = t628;
  local_mass_matrix(22,17) = t655;
  local_mass_matrix(22,18) = t670;
  local_mass_matrix(22,19) = t678;
  local_mass_matrix(22,20) = t694;
  local_mass_matrix(22,21) = t676;
  local_mass_matrix(22,22) = t681-t682+t683+73216.0/893025.0*t12-t685
			     -73216.0/893025.0*t17;
  local_mass_matrix(22,23) = t680;
  local_mass_matrix(22,24) = t717;
  local_mass_matrix(23,0) = t104;
  local_mass_matrix(23,1) = t162;
  local_mass_matrix(23,2) = t220;
  local_mass_matrix(23,3) = t236;
  local_mass_matrix(23,4) = t317;
  local_mass_matrix(23,5) = t322;
  local_mass_matrix(23,6) = t307;
  local_mass_matrix(23,7) = t444;
  local_mass_matrix(23,8) = t477;
  local_mass_matrix(23,9) = t464;
  local_mass_matrix(23,10) = t530;
  local_mass_matrix(23,11) = t531;
  local_mass_matrix(23,12) = t524;
  local_mass_matrix(23,13) = t568;
  local_mass_matrix(23,14) = t578;
  local_mass_matrix(23,15) = t573;
  local_mass_matrix(23,16) = t633;
  local_mass_matrix(23,17) = t623;
  local_mass_matrix(23,18) = t673;
  local_mass_matrix(23,19) = t679;
  local_mass_matrix(23,20) = t640;
  local_mass_matrix(23,21) = t694;
  local_mass_matrix(23,22) = t680;
  local_mass_matrix(23,23) = 6656.0/297675.0*t3-6656.0/297675.0*t7+t683+
			     t704-t685-t705;
  local_mass_matrix(23,24) = t723;
  local_mass_matrix(24,0) = 841.0/893025.0*t105;
  local_mass_matrix(24,1) = 841.0/893025.0*t163;
  local_mass_matrix(24,2) = 841.0/893025.0*t221;
  local_mass_matrix(24,3) = 841.0/893025.0*t237;
  local_mass_matrix(24,4) = t322;
  local_mass_matrix(24,5) = t366;
  local_mass_matrix(24,6) = t356;
  local_mass_matrix(24,7) = t447;
  local_mass_matrix(24,8) = t482;
  local_mass_matrix(24,9) = t472;
  local_mass_matrix(24,10) = t531;
  local_mass_matrix(24,11) = t547;
  local_mass_matrix(24,12) = t541;
  local_mass_matrix(24,13) = t569;
  local_mass_matrix(24,14) = t579;
  local_mass_matrix(24,15) = t577;
  local_mass_matrix(24,16) = t640;
  local_mass_matrix(24,17) = t658;
  local_mass_matrix(24,18) = t676;
  local_mass_matrix(24,19) = t680;
  local_mass_matrix(24,20) = t701;
  local_mass_matrix(24,21) = t711;
  local_mass_matrix(24,22) = t717;
  local_mass_matrix(24,23) = t723;
  local_mass_matrix(24,24) = 5408.0/99225.0*t3-5408.0/99225.0*t7+10816.0/
			     99225.0*t10+5408.0/99225.0*t12-10816.0/99225.0*t15-5408.0/99225.0*t17;
};



template <>
void FEQ4<2>::get_unit_support_points (vector<Point<2> > &unit_points) const {
  Assert (unit_points.size() == total_dofs,
	  ExcWrongFieldDimension (unit_points.size(), total_dofs));

  unit_points[0] = Point<2>(0,0);
  unit_points[1] = Point<2>(1,0);
  unit_points[2] = Point<2>(1,1);
  unit_points[3] = Point<2>(0,1);
  unit_points[4] = Point<2>(1./4,0);
  unit_points[5] = Point<2>(2./4,0);
  unit_points[6] = Point<2>(3./4,0);
  unit_points[7] = Point<2>(1,1./4);
  unit_points[8] = Point<2>(1,2./4);
  unit_points[9] = Point<2>(1,3./4);
  unit_points[10] = Point<2>(1./4,1);
  unit_points[11] = Point<2>(2./4,1);
  unit_points[12] = Point<2>(3./4,1);
  unit_points[13] = Point<2>(0,1./4);
  unit_points[14] = Point<2>(0,2./4);
  unit_points[15] = Point<2>(0,3./4);
  unit_points[16] = Point<2>(1./4,1./4);
  unit_points[17] = Point<2>(3./4,1./4);
  unit_points[18] = Point<2>(3./4,3./4);
  unit_points[19] = Point<2>(1./4,3./4);
  unit_points[20] = Point<2>(1./2,1./4);
  unit_points[21] = Point<2>(3./4,1./2);
  unit_points[22] = Point<2>(1./2,3./4);
  unit_points[23] = Point<2>(1./4,1./2);
  unit_points[24] = Point<2>(1./2,1./2);
};



template <>
void FEQ4<2>::get_support_points (const DoFHandler<2>::cell_iterator &cell,
				  vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == total_dofs,
	  ExcWrongFieldDimension (support_points.size(), total_dofs));

  const double x[4] = { cell->vertex(0)(0),
			cell->vertex(1)(0),
			cell->vertex(2)(0),
			cell->vertex(3)(0)  };
  const double y[4] = { cell->vertex(0)(1),
			cell->vertex(1)(1),
			cell->vertex(2)(1),
			cell->vertex(3)(1)  };

  const double t1 = 3.0/4.0*x[0];
  const double t2 = x[1]/4.0;
  const double t4 = 3.0/4.0*y[0];
  const double t5 = y[1]/4.0;
  const double t9 = x[0]/4.0;
  const double t10 = 3.0/4.0*x[1];
  const double t12 = y[0]/4.0;
  const double t13 = 3.0/4.0*y[1];
  const double t15 = x[2]/4.0;
  const double t17 = y[2]/4.0;
  const double t21 = 3.0/4.0*x[2];
  const double t23 = 3.0/4.0*y[2];
  const double t25 = 3.0/4.0*x[3];
  const double t27 = 3.0/4.0*y[3];
  const double t31 = x[3]/4.0;
  const double t33 = y[3]/4.0;
  const double t42 = 3.0/16.0*x[1];
  const double t44 = 3.0/16.0*x[3];
  const double t47 = 3.0/16.0*y[1];
  const double t49 = 3.0/16.0*y[3];
  const double t51 = 3.0/16.0*x[0];
  const double t53 = 3.0/16.0*x[2];
  const double t56 = 3.0/16.0*y[0];
  const double t58 = 3.0/16.0*y[2];
  const double t73 = 3.0/8.0*x[0];
  const double t74 = 3.0/8.0*x[1];
  const double t75 = x[2]/8.0;
  const double t76 = x[3]/8.0;
  const double t78 = 3.0/8.0*y[0];
  const double t79 = 3.0/8.0*y[1];
  const double t80 = y[2]/8.0;
  const double t81 = y[3]/8.0;
  const double t83 = x[0]/8.0;
  const double t84 = 3.0/8.0*x[2];
  const double t86 = y[0]/8.0;
  const double t87 = 3.0/8.0*y[2];
  const double t89 = x[1]/8.0;
  const double t90 = 3.0/8.0*x[3];
  const double t92 = y[1]/8.0;
  const double t93 = 3.0/8.0*y[3];
  support_points[0](0) = x[0];
  support_points[0](1) = y[0];
  support_points[1](0) = x[1];
  support_points[1](1) = y[1];
  support_points[2](0) = x[2];
  support_points[2](1) = y[2];
  support_points[3](0) = x[3];
  support_points[3](1) = y[3];
  support_points[4](0) = t1+t2;
  support_points[4](1) = t4+t5;
  support_points[5](0) = x[0]/2.0+x[1]/2.0;
  support_points[5](1) = y[0]/2.0+y[1]/2.0;
  support_points[6](0) = t9+t10;
  support_points[6](1) = t12+t13;
  support_points[7](0) = t10+t15;
  support_points[7](1) = t13+t17;
  support_points[8](0) = x[1]/2.0+x[2]/2.0;
  support_points[8](1) = y[1]/2.0+y[2]/2.0;
  support_points[9](0) = t2+t21;
  support_points[9](1) = t5+t23;
  support_points[10](0) = t15+t25;
  support_points[10](1) = t17+t27;
  support_points[11](0) = x[2]/2.0+x[3]/2.0;
  support_points[11](1) = y[2]/2.0+y[3]/2.0;
  support_points[12](0) = t21+t31;
  support_points[12](1) = t23+t33;
  support_points[13](0) = t1+t31;
  support_points[13](1) = t4+t33;
  support_points[14](0) = x[0]/2.0+x[3]/2.0;
  support_points[14](1) = y[0]/2.0+y[3]/2.0;
  support_points[15](0) = t9+t25;
  support_points[15](1) = t12+t27;
  support_points[16](0) = 9.0/16.0*x[0]+t42+x[2]/16.0+t44;
  support_points[16](1) = 9.0/16.0*y[0]+t47+y[2]/16.0+t49;
  support_points[17](0) = t51+9.0/16.0*x[1]+t53+x[3]/16.0;
  support_points[17](1) = t56+9.0/16.0*y[1]+t58+y[3]/16.0;
  support_points[18](0) = x[0]/16.0+t42+9.0/16.0*x[2]+t44;
  support_points[18](1) = y[0]/16.0+t47+9.0/16.0*y[2]+t49;
  support_points[19](0) = t51+x[1]/16.0+t53+9.0/16.0*x[3];
  support_points[19](1) = t56+y[1]/16.0+t58+9.0/16.0*y[3];
  support_points[20](0) = t73+t74+t75+t76;
  support_points[20](1) = t78+t79+t80+t81;
  support_points[21](0) = t83+t74+t84+t76;
  support_points[21](1) = t86+t79+t87+t81;
  support_points[22](0) = t83+t89+t84+t90;
  support_points[22](1) = t86+t92+t87+t93;
  support_points[23](0) = t73+t89+t75+t90;
  support_points[23](1) = t78+t92+t80+t93;
  support_points[24](0) = x[0]/4.0+x[1]/4.0+x[2]/4.0+x[3]/4.0;
  support_points[24](1) = y[0]/4.0+y[1]/4.0+y[2]/4.0+y[3]/4.0;
};



template <>
void FEQ4<2>::get_face_support_points (const DoFHandler<2>::face_iterator &face,
				       vector<Point<2> >  &support_points) const {
  Assert (support_points.size() == dofs_per_face,
	  ExcWrongFieldDimension (support_points.size(), dofs_per_face));

  for (unsigned int vertex=0; vertex<2; ++vertex)
    support_points[vertex] = face->vertex(vertex);
  support_points[2] = (3*support_points[0] + support_points[1]) / 4;
  support_points[3] = (support_points[0] + support_points[1]) / 2;
  support_points[4] = (support_points[0] + 3*support_points[1]) / 4;
};



#endif







// explicit instantiations

template class FEQ4<deal_II_dimension>;

