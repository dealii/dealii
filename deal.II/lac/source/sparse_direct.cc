//----------------------------  sparse_direct.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_direct.cc  ---------------------------


#include <lac/sparse_direct.h>
#include <base/memory_consumption.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>

// if we know that at least one of the HSL functions are there,
// include the respective include file. Otherwise save some CPU cycles
// in the compiler
#if HAVE_HSL_MA27 || HAVE_HSL_MA47
#  include <hsl/hsl.h>
#endif

// if the HSL functions are not there, define them empty and throw an
// exception
#ifndef HAVE_HSL_MA27
namespace HSL
{
  namespace MA27
  {
    extern "C"
    void ma27ad_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  unsigned int       *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

    
    extern "C"
    void ma27bd_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };


    extern "C"
    void ma27cd_ (const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };


    extern "C" void ma27x1_ (unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
      

    extern "C" void ma27x2_ (unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
      
    
    extern "C" void ma27x3_ (const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
  };
};
#endif  // ifndef HAVE_HSL_MA27


#ifndef HAVE_HSL_MA47
namespace HSL
{
  namespace MA47
  {
    extern "C"
    void ma47id_ (double       *,
		  unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };
    

    extern "C"
    void ma47ad_ (const unsigned int *,
		  const unsigned int *,
		  unsigned int       *,
		  unsigned int       *,
		  unsigned int       *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  double             *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

      
    extern "C"
    void ma47bd_ (const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  const unsigned int *,
		  unsigned int       *,
		  const unsigned int *,
		  const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  unsigned int       *,
		  double             *,
		  int                *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };

    
    extern "C"
    void ma47cd_ (const unsigned int *,
		  const double       *,
		  const unsigned int *,
		  const unsigned int *,
		  const unsigned int *,
		  double             *,
		  double             *,
		  unsigned int       *,
		  const unsigned int *)
    {
      AssertThrow (false, ExcMessage("You can only use the HSL functions after putting "
				     "the respective files in the right place, "
				     "re-configuring the library and re-building it!"));
    };    
  };
};
#endif   // ifndef HAVE_HSL_MA47





/* -------------------------- MA27 ---------------------------- */


SparseDirectMA27::SparseDirectMA27 (const double LIW_factor_1,
				    const double LIW_factor_2,
				    const double LA_factor,
				    const double LIW_increase_factor_1,
				    const double LIW_increase_factor_2,
				    const double LA_increase_factor,
				    const bool   suppress_output) :
		LIW_factor_1 (LIW_factor_1),
		LIW_factor_2 (LIW_factor_2),
		LA_factor (LA_factor),
		LIW_increase_factor_1 (LIW_increase_factor_1),
		LIW_increase_factor_2 (LIW_increase_factor_2),
		LA_increase_factor (LA_increase_factor),
		initialize_called (false),
		factorize_called (false),
		sparsity_pattern (0)
{
				   // suppress error output if
				   // requested
  if (suppress_output)
    {
      const unsigned int LP = 0;
      HSL::MA27::ma27x3_ (&LP);
    };
};



void
SparseDirectMA27::initialize (const SparsityPattern &sp)
{
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());

  sparsity_pattern = &sp;
  
  const unsigned int n_rows = sparsity_pattern->n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern->get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern->get_column_numbers();

				   // first count number of nonzero
				   // elements in the upper right
				   // part. the matrix is symmetric,
				   // so this suffices
  n_nonzero_elements = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
      if (row <= *col)
	++n_nonzero_elements;
  

				   // fill the row numbers and column
				   // numbers arrays from the sparsity
				   // pattern. note that we have
				   // Fortran convention, i.e. indices
				   // need to be 1-base, as opposed to
				   // C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if (row <= *col)
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  row_numbers[global_index] = row+1;
	  column_numbers[global_index] = *col+1;
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());
  
				   // initialize scratch arrays and
				   // variables
  LIW = static_cast<unsigned int>((2*n_nonzero_elements + 3*n_rows + 1) *
				  LIW_factor_1);
  IW.resize (LIW);
  IKEEP.resize (3*n_rows);
  IW1.resize (2*n_rows);

				   // no output please
  IFLAG = 0;

				   // loop until memory requirements
				   // are satisfied or we are not
				   // allowed to allocate more memory
				   // no more
  bool call_succeeded = true;
  do 
    {
      HSL::MA27::ma27ad_(&n_rows, &n_nonzero_elements,
			 &row_numbers[0], &column_numbers[0],
			 &IW[0], &LIW, &IKEEP[0],
			 &IW1[0], &NSTEPS, &IFLAG);
      call_succeeded = (IFLAG==0);

				       // if enough memory or no
				       // increase allowed: exit loop
      if (call_succeeded || (LIW_increase_factor_1 <= 1))
	break;
      
				       // otherwise: increase LIW and retry
      LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);
    }
  while (true);

				   // if we were not allowed to
				   // allocate more memory, then throw
				   // an exception
  if (!call_succeeded)
    throw ExcMA27AFailed(IFLAG);

				   // catch returned values from the
				   // COMMON block. we need these
				   // values in order to set array
				   // sizes in the next function
  HSL::MA27::ma27x1_(&NRLNEC);
  HSL::MA27::ma27x2_(&NIRNEC);

				   // note that we have already been
				   // in this function
  initialize_called = true;
};



void
SparseDirectMA27::factorize (const SparseMatrix<double> &matrix)
{
				   // if necessary, initialize process
  if (initialize_called == false)
    initialize (matrix.get_sparsity_pattern());

				   // make sure the sparsity patterns
				   // are the same
  Assert (sparsity_pattern == &matrix.get_sparsity_pattern(),
	  ExcDifferentSparsityPatterns());
  
  
				   // set LA and fill the A array of
				   // values
  LA=static_cast<int>(NRLNEC * LA_factor);
  A.resize (LA);
  fill_A (matrix);

				   // if necessary extend IW
  if (LIW < NIRNEC * LIW_factor_2)
    {
      LIW = static_cast<unsigned int>(NIRNEC * LIW_factor_2);
      IW.resize (LIW);
    };
  
  const unsigned int n_rows = matrix.get_sparsity_pattern().n_rows();
  
				   // loop until memory requirements
				   // are satisfied or we are not
				   // allowed to allocate more memory
				   // no more
  bool call_succeeded = true;
  do 
    {
      HSL::MA27::ma27bd_(&n_rows, &n_nonzero_elements,
			 &row_numbers[0], &column_numbers[0],
			 &A[0], &LA,
			 &IW[0], &LIW, &IKEEP[0], &NSTEPS, &MAXFRT,
			 &IW1[0], &IFLAG);
      call_succeeded = (IFLAG==0);

				       // if enough memory or no
				       // increase allowed: exit loop
      if (call_succeeded)
	break;

				       // otherwise: increase LIW or
				       // LA if that is allowed and
				       // retry
      switch (IFLAG)
	{
	  case -3:
	  {
	    if (LIW_increase_factor_2 <= 1)
	      goto exit_loop;
	    
	    LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_2);
	    IW.resize (LIW);
	    break;
	  };

	  case -4:
	  {
	    if (LA_increase_factor <= 1)
	      goto exit_loop;
					     // increase A. note that
					     // since the function has
					     // already part of the
					     // array @p{A}, we have
					     // to re-fill it with the
					     // original values. minor
					     // clue: since the old
					     // entries are no more
					     // needed, we can discard
					     // them; we use this to
					     // first release all
					     // memory (through the
					     // call to @p{swap} and
					     // the subsequent call to
					     // the destructor of the
					     // @p{tmp} object) and
					     // only then re-allocate
					     // it. If we called
					     // @p{resize} directly,
					     // this would first
					     // allocate more memory,
					     // then copy the old
					     // contents, and only
					     // then release the old
					     // memory, but keeping
					     // both memory regions at
					     // the same time could
					     // sometimes be more than
					     // we can do, leading to
					     // an exception on the
					     // allocation.
	    std::cout << "<*>" << std::flush;
	    
	    LA  = static_cast<unsigned int>(LA * LA_increase_factor);
	    if (true)
	      {
		std::vector<double> tmp;
		A.swap (tmp);
	      };
	    
	    A.resize (LA);
	    fill_A (matrix);
	    
	    break;
	  };
	   
					    // ups, other return
					    // value, don't know
					    // what to do here
	  default:
		throw ExcMA27BFailed(IFLAG);
	};
      continue;

      exit_loop:
      break;
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA27BFailed(IFLAG);

				   // note that we have been here
				   // already
  factorize_called = true;
};



void
SparseDirectMA27::solve (Vector<double> &rhs_and_solution) const
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());
  
  const unsigned int n_rows = rhs_and_solution.size();
  std::vector<double> W(MAXFRT);
  HSL::MA27::ma27cd_(&n_rows, &A[0], &LA,
		     &IW[0], &LIW, &W[0], &MAXFRT,
		     &rhs_and_solution(0), &IW1[0], &NSTEPS);
};



void
SparseDirectMA27::solve (const SparseMatrix<double> &matrix,
			 Vector<double>             &rhs_and_solution)
{
  initialize (matrix.get_sparsity_pattern());
  factorize (matrix);
  solve (rhs_and_solution);
};



unsigned int
SparseDirectMA27::memory_consumption () const
{
  return (sizeof(*this) +
	  MemoryConsumption::memory_consumption (row_numbers) +
	  MemoryConsumption::memory_consumption (column_numbers) +
	  MemoryConsumption::memory_consumption (A) +
	  MemoryConsumption::memory_consumption (IW) +
	  MemoryConsumption::memory_consumption (IKEEP) +
	  MemoryConsumption::memory_consumption (IW1));
};



void
SparseDirectMA27::fill_A (const SparseMatrix<double> &matrix)
{
  
  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();
  
  const unsigned int n_rows = sparsity_pattern.n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern.get_column_numbers();

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if (row <= *col)
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  A[global_index] = matrix(row,*col);
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());  
};



/* -------------------------- MA47 ---------------------------- */




SparseDirectMA47::SparseDirectMA47 (const double LIW_factor_1,
				    const double LIW_factor_2,
				    const double LA_factor,
				    const double LIW_increase_factor_1,
				    const double LIW_increase_factor_2,
				    const double LA_increase_factor,
				    const bool   suppress_output) :
		LIW_factor_1 (LIW_factor_1),
		LIW_factor_2 (LIW_factor_2),
		LA_factor (LA_factor),
		LIW_increase_factor_1 (LIW_increase_factor_1),
		LIW_increase_factor_2 (LIW_increase_factor_2),
		LA_increase_factor (LA_increase_factor),
		initialize_called (false),
		factorize_called (false),
		matrix (0)
{
  HSL::MA47::ma47id_ (CNTL, ICNTL);

				   // suppress error output if
				   // requested
  if (suppress_output)
    ICNTL[0] = 0;
};



void
SparseDirectMA47::initialize (const SparseMatrix<double> &m)
{
  Assert (initialize_called == false, ExcInitializeAlreadyCalled());

  matrix = &m;
  const SparsityPattern &sparsity_pattern = matrix->get_sparsity_pattern();
  
  const unsigned int n_rows = sparsity_pattern.n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern.get_column_numbers();

				   // first count number of nonzero
				   // elements in the upper right
				   // part. the matrix is symmetric,
				   // so this suffices
  n_nonzero_elements = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // skip zero elements, as
				       // required by the docs of MA47
      if ((row <= *col) && (m(row,*col) != 0))
	++n_nonzero_elements;
  

				   // fill the row numbers and column
				   // numbers arrays from the sparsity
				   // pattern. note that we have
				   // Fortran convention, i.e. indices
				   // need to be 1-base, as opposed to
				   // C's 0-based convention!
  row_numbers.resize (n_nonzero_elements);
  column_numbers.resize (n_nonzero_elements);

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if ((row <= *col) && (m(row,*col) != 0))
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  row_numbers[global_index] = row+1;
	  column_numbers[global_index] = *col+1;
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());
  
				   // initialize scratch arrays and
				   // variables
  LIW = static_cast<unsigned int>((2*n_nonzero_elements + 5*n_rows + 4) *
				  LIW_factor_1);
  IW.resize (LIW);
  KEEP.resize (n_nonzero_elements + 5*n_rows + 2);

				   // declare output info fields
  double RINFO[4];
  bool call_succeeded;
  do
    {
      HSL::MA47::ma47ad_(&n_rows, &n_nonzero_elements,
			 &row_numbers[0], &column_numbers[0],
			 &IW[0], &LIW, &KEEP[0],
			 &ICNTL[0], &RINFO[0], &INFO[0]);
      call_succeeded = (INFO[0] == 0);

				       // if enough memory or no
				       // increase allowed: exit loop
      if (call_succeeded || (LIW_increase_factor_1 <= 1))
	break;
      
				       // otherwise: increase LIW and retry
      LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_1);
      IW.resize (LIW);      
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA47AFailed(INFO[0]);

				   // note that we have already been
				   // in this function
  initialize_called = true;
};



void
SparseDirectMA47::factorize (const SparseMatrix<double> &m)
{
  Assert (factorize_called == false, ExcCantFactorizeAgain());
  
				   // if necessary, initialize process
  if (initialize_called == false)
    initialize (m);

				   // make sure the matrices
				   // are the same
  Assert (matrix == &m, ExcDifferentMatrices());
  
  
				   // set LA and fill the A array of
				   // values
  LA=static_cast<int>(INFO[5] * LA_factor);
  A.resize (LA);
  fill_A (m);
  
				   // if necessary extend IW
  if (LIW < INFO[6] * LIW_factor_2)
    {
      LIW = static_cast<unsigned int>(INFO[6] * LIW_factor_2);
      IW.resize (LIW);
    };

  const unsigned int n_rows = m.get_sparsity_pattern().n_rows();
  IW1.resize (2*n_rows+2);

				   // output info flags
  double RINFO[4];

  bool call_succeeded;
  do 
    {    
      HSL::MA47::ma47bd_(&n_rows, &n_nonzero_elements, &column_numbers[0],
			 &A[0], &LA,
			 &IW[0], &LIW, &KEEP[0], &CNTL[0], &ICNTL[0],
			 &IW1[0], &RINFO[0], &INFO[0]);
      call_succeeded = (INFO[0] == 0);

				       // if enough memory or no
				       // increase allowed: exit loop
      if (call_succeeded)
	break;

				       // otherwise: increase LIW or
				       // LA if that is allowed and
				       // retry
      switch (INFO[0])
	{
	  case -3:
	  {
	    if (LIW_increase_factor_2 <= 1)
	      goto exit_loop;
	    
	    LIW = static_cast<unsigned int>(LIW * LIW_increase_factor_2);
	    IW.resize (LIW);
	    break;
	  };

	  case -4:
	  {
	    if (LA_increase_factor <= 1)
	      goto exit_loop;
					     // increase A. note that
					     // since the function has
					     // already part of the
					     // array @p{A}, we have
					     // to re-fill it with the
					     // original values. minor
					     // clue: since the old
					     // entries are no more
					     // needed, we can discard
					     // them; we use this to
					     // first release all
					     // memory (through the
					     // call to @p{swap} and
					     // the subsequent call to
					     // the destructor of the
					     // @p{tmp} object) and
					     // only then re-allocate
					     // it. If we called
					     // @p{resize} directly,
					     // this would first
					     // allocate more memory,
					     // then copy the old
					     // contents, and only
					     // then release the old
					     // memory, but keeping
					     // both memory regions at
					     // the same time could
					     // sometimes be more than
					     // we can do, leading to
					     // an exception on the
					     // allocation.
	    std::cout << "<*>" << std::flush;
	    
	    LA  = static_cast<unsigned int>(LA * LA_increase_factor);
	    if (true)
	      {
		std::vector<double> tmp;
		A.swap (tmp);
	      };
	    
	    A.resize (LA);
	    fill_A (m);
	    
	    break;
	  };
	   
					    // ups, other return
					    // value, don't know
					    // what to do here
	  default:
		throw ExcMA47BFailed(INFO[0]);
	};
      continue;

      exit_loop:
      break;
    }
  while (true);

  if (!call_succeeded)
    throw ExcMA47BFailed(INFO[0]);

				   // note that we have been here
				   // already
  factorize_called = true;
};



void
SparseDirectMA47::solve (Vector<double> &rhs_and_solution)
{
  Assert (factorize_called == true, ExcFactorizeNotCalled());
  
  const unsigned int n_rows = rhs_and_solution.size();
  std::vector<double> W(n_rows);
  HSL::MA47::ma47cd_(&n_rows, &A[0], &LA,
		     &IW[0], &LIW, &W[0],
		     &rhs_and_solution(0), &IW1[0], &ICNTL[0]);
};



void
SparseDirectMA47::solve (const SparseMatrix<double> &matrix,
			 Vector<double>             &rhs_and_solution)
{
  initialize (matrix);
  factorize (matrix);
  solve (rhs_and_solution);
};



unsigned int
SparseDirectMA47::memory_consumption () const
{
  return (sizeof(*this) +
	  MemoryConsumption::memory_consumption (row_numbers) +
	  MemoryConsumption::memory_consumption (column_numbers) +
	  MemoryConsumption::memory_consumption (A) +
	  MemoryConsumption::memory_consumption (IW) +
	  MemoryConsumption::memory_consumption (KEEP) +
	  MemoryConsumption::memory_consumption (IW1));
};



void
SparseDirectMA47::fill_A (const SparseMatrix<double> &matrix)
{
  const SparsityPattern &sparsity_pattern = matrix.get_sparsity_pattern ();
  
  const unsigned int n_rows = sparsity_pattern.n_rows();
  const unsigned int *rowstart_indices = sparsity_pattern.get_rowstart_indices();
  const unsigned int *col_nums         = sparsity_pattern.get_column_numbers();

  unsigned int global_index = 0;
  for (unsigned int row=0; row<n_rows; ++row)
    for (const unsigned int *col=&col_nums[rowstart_indices[row]];
	 col != &col_nums[rowstart_indices[row+1]];
	 ++col)
				       // note that the matrix must be
				       // symmetric, so only treat the
				       // upper right part
      if ((row <= *col) && (matrix(row,*col) != 0))
	{
	  Assert (global_index < n_nonzero_elements, ExcInternalError());
	  
	  A[global_index] = matrix(row,*col);
	  ++global_index;
	};
  Assert (global_index == n_nonzero_elements, ExcInternalError());  
};
