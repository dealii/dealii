//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_vector_base.h>

#include <cmath>

#ifdef DEAL_II_USE_TRILINOS

#include <Epetra_Import.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace internal
  {
    VectorReference::operator TrilinosScalar () const
    {
      Assert (index < vector.size(),
	      ExcIndexRange (index, 0, vector.size()));

                                        // Trilinos allows for vectors
                                        // to be referenced by the [] or
                                        // () operators but only ()
                                        // checks index bounds. We check
                                        // these bounds by ourselves, so
                                        // we can use []. Note that we
                                        // can only get local values.

      const int local_index = vector.vector->Map().LID(index);
      Assert (local_index >= 0,
	      ExcAccessToNonLocalElement (index, 
					  vector.vector->Map().MinMyGID(),
					  vector.vector->Map().MaxMyGID()));


      return (*(vector.vector))[0][local_index];
    }
  }



  VectorBase::VectorBase ()
                        :
                        last_action (Zero),
			compressed  (true),
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
			vector(std::auto_ptr<Epetra_FEVector> 
			  (new Epetra_FEVector(
				 Epetra_Map(0,0,Epetra_MpiComm(MPI_COMM_WORLD)))))
#else
			vector(std::auto_ptr<Epetra_FEVector> 
                          (new Epetra_FEVector(
				 Epetra_Map(0,0,Epetra_SerialComm()))))
#endif
  {}
  

  
  VectorBase::VectorBase (const VectorBase &v)
                        :
			Subscriptor(),
			last_action (Zero),
			compressed (true),
			vector(std::auto_ptr<Epetra_FEVector> 
			       (new Epetra_FEVector(*v.vector)))
  {}
  


  VectorBase::~VectorBase ()
  {}



  void
  VectorBase::clear ()
  {
                                     // When we clear the vector,
				     // reset the pointer and generate
				     // an empty vector.
    vector.reset();

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    Epetra_Map map (0, 0, Epetra_MpiComm(MPI_COMM_WORLD));
#else
    Epetra_Map map (0, 0, Epetra_SerialComm());
#endif

    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
    last_action = Zero;
  }



  void
  VectorBase::reinit (const VectorBase &v,
		      const bool        fast)
  {
    Assert (&*vector != 0,
	    ExcMessage("Vector has not been constructed properly."));

    if (fast == false || local_range() != v.local_range())
      vector = std::auto_ptr<Epetra_FEVector>(new Epetra_FEVector(*v.vector));
  }



  void
  VectorBase::compress ()
  {
				 // Now pass over the information about
				 // what we did last to the vector.
    const int ierr = vector->GlobalAssemble(last_action);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    last_action = Zero;
  
    compressed = true;
  }



  VectorBase &
  VectorBase::operator = (const TrilinosScalar s)
  {

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but either "
		       "infinite or Not A Number (NaN)"));

    compress();

    const int ierr = vector->PutScalar(s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator = (const VectorBase &v)
  {
    Assert (&*vector != 0,
	    ExcMessage("Vector is not constructed properly."));

    if (local_range() != v.local_range())
      {
	vector.reset();
	last_action = Zero;
	vector = std::auto_ptr<Epetra_FEVector>(new Epetra_FEVector(*v.vector));
      }
    else
      {
	Assert (vector->Map().SameAs(v.vector->Map()) == true,
		ExcMessage ("The Epetra maps in the assignment operator ="
			    " do not match, even though the local_range "
			    " seems to be the same. Check vector setup!"));
	int ierr;
	ierr = vector->GlobalAssemble(last_action);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	ierr = vector->Update(1.0, *v.vector, 0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	last_action = Zero;
      }

    return *this;
  }



  template <typename number>
  VectorBase &
  VectorBase::operator = (const ::dealii::Vector<number> &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));
    
				     // this is probably not very efficient
				     // but works. in particular, we could do
				     // better if we know that
				     // number==TrilinosScalar because then we
				     // could elide the copying of elements
				     //
				     // let's hope this isn't a
				     // particularly frequent operation
    std::pair<unsigned int, unsigned int>
      local_range = this->local_range ();
    std::vector<unsigned int>   indices (local_range.second -
					 local_range.first);
    std::vector<TrilinosScalar> values (local_range.second -
					local_range.first);

    for (unsigned int i=0; i<local_range.second-local_range.first; ++i)
      {
	indices[i] = i + local_range.first;
	values[i] = v(i + local_range.first);
      }
    
    set (indices.size(), &indices[0], &values[0]);
      
    return *this;
  }
  


  bool
  VectorBase::operator == (const VectorBase &v) const
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));
    unsigned int i;
    for (i=0; i<size(); i++) 
      if ((*(v.vector))[0][i]!=(*vector)[0][i]) return false;

    return true;
  }



  bool
  VectorBase::operator != (const VectorBase &v) const
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    return (!(*this==v));
  }



  unsigned int
  VectorBase::size () const
  {
    return (unsigned int) vector->Map().NumGlobalElements();
  }



  unsigned int
  VectorBase::local_size () const
  {
    return (unsigned int) vector->Map().NumMyElements();
  }



  TrilinosScalar
  VectorBase::el (const unsigned int index) const
  {
                                        // Extract local indices in
                                        // the vector.
    int trilinos_i = vector->Map().LID(index);
    TrilinosScalar value = 0.;

				        // If the element is not
				        // present on the current
				        // processor, we can't
				        // continue. Just print out 0.

				        // TODO: Is this reasonable?
    if (trilinos_i == -1 )
      {
	return 0.;
        //Assert (false, ExcAccessToNonlocalElement(index, local_range().first,
	//				  local_range().second));
      }
    else
      value = (*vector)[0][trilinos_i];

    return value;
  }



  TrilinosScalar
  VectorBase::operator () (const unsigned int index) const
  {
                                        // Extract local indices in
                                        // the vector.
    int trilinos_i = vector->Map().LID(index);
    TrilinosScalar value = 0.;

				        // If the element is not present
				        // on the current processor, we
				        // can't continue. This is the
				        // main difference to the el()
				        // function.
    if (trilinos_i == -1 )
      {
	Assert (false, ExcAccessToNonlocalElement(index, local_range().first,
						  local_range().second));
      }
    else
      value = (*vector)[0][trilinos_i];

    return value;
  }



  TrilinosScalar
  VectorBase::operator * (const VectorBase &vec) const
  {
    Assert (local_range() == vec.local_range(),
	    ExcDimensionMismatch(size(), vec.size()));

    TrilinosScalar result;

    const int ierr = vector->Dot(*(vec.vector), &result);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return result;
  }



  VectorBase::real_type
  VectorBase::norm_sqr () const
  {
    const TrilinosScalar d = l2_norm();
    return d*d;
  }



  TrilinosScalar
  VectorBase::mean_value () const
  {
    TrilinosScalar mean;

    const int ierr = vector->MeanValue (&mean);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return mean;
  }



  VectorBase::real_type
  VectorBase::l1_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->Norm1 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::l2_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->Norm2 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::lp_norm (const TrilinosScalar p) const
  {
                                        // get a representation of the
                                        // vector and loop over all
                                        // the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    TrilinosScalar norm = 0;
    TrilinosScalar sum=0;

    const TrilinosScalar * ptr  = start_ptr;

                                       // add up elements 
				       // TODO: This
				       // won't work in parallel like
				       // this. Find out a better way to
				       // this in that case.
    while (ptr != start_ptr+size())
      sum += std::pow(std::fabs(*ptr++), p);

    norm = std::pow(sum, static_cast<TrilinosScalar>(1./p));

    return norm;
  }



  VectorBase::real_type
  VectorBase::linfty_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->NormInf (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  bool
  VectorBase::all_zero () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

				       // TODO: This
				       // won't work in parallel like
				       // this. Find out a better way to
				       // this in that case.
    const TrilinosScalar *ptr  = start_ptr,
                         *eptr = start_ptr + size();
    bool flag = true;
    while (ptr != eptr)
      {
	if (*ptr != 0)
	  {
	    flag = false;
	    break;
	  }
	++ptr;
      }

    return flag;
  }



  bool
  VectorBase::is_non_negative () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

				       // TODO: This
				       // won't work in parallel like
				       // this. Find out a better way to
				       // this in that case.
    const TrilinosScalar *ptr  = start_ptr,
                         *eptr = start_ptr + size();
    bool flag = true;
    while (ptr != eptr)
      {
	if (*ptr < 0.0)
	  {
	    flag = false;
	    break;
	  }
	++ptr;
      }

    return flag;
  }



  VectorBase &
  VectorBase::operator *= (const TrilinosScalar a)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Scale(a);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator /= (const TrilinosScalar a)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const TrilinosScalar factor = 1./a;

    Assert (numbers::is_finite(factor),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Scale(factor);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator += (const VectorBase &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    const int ierr = vector->Update (1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator -= (const VectorBase &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    const int ierr = vector->Update (-1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  void
  VectorBase::add (const TrilinosScalar s)
  {

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    unsigned int n_local = local_size();
    int ierr;

    for (unsigned int i=0; i<n_local; i++)
      {
	ierr = vector->SumIntoMyValue(i,0,s);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
  }



  void
  VectorBase::add (const VectorBase &v,
		   const bool        allow_different_maps)
  {
    if (allow_different_maps == false)
      *this += v;
    else
      {
	AssertThrow (size() == v.size(),
		     ExcDimensionMismatch (size(), v.size()));

	Epetra_Import data_exchange (vector->Map(), v.vector->Map());

	int ierr = vector->Import(*v.vector, data_exchange, Add);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	last_action = Insert;
      }
  }



  void
  VectorBase::add (const TrilinosScalar  a,
		   const VectorBase     &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(a, *(v.vector), 1.);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::add (const TrilinosScalar  a,
		   const VectorBase     &v,
		   const TrilinosScalar  b,
		   const VectorBase     &w)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));
    Assert (size() == w.size(),
	    ExcDimensionMismatch(size(), w.size()));

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(b),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), 1.);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const VectorBase     &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(1., *(v.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const TrilinosScalar  a,
		    const VectorBase     &v)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(a, *(v.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const TrilinosScalar  a,
		    const VectorBase     &v,
		    const TrilinosScalar  b,
		    const VectorBase     &w)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));
    Assert (size() == w.size(),
	    ExcDimensionMismatch(size(), w.size()));


    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(b),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::sadd (const TrilinosScalar  s,
		    const TrilinosScalar  a,
		    const VectorBase     &v,
		    const TrilinosScalar  b,
		    const VectorBase     &w,
		    const TrilinosScalar  c,
		    const VectorBase     &x)
  {
    Assert (size() == v.size(),
	    ExcDimensionMismatch(size(), v.size()));
    Assert (size() == w.size(),
	    ExcDimensionMismatch(size(), w.size()));
    Assert (size() == x.size(),
	    ExcDimensionMismatch(size(), x.size()));

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(b),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(c),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

                                        // Update member can only
				        // input two other vectors so
				        // do it in two steps
    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    const int jerr = vector->Update(c, *(x.vector), 1.);
    Assert (jerr == 0, ExcTrilinosError(jerr));
  }



  void
  VectorBase::scale (const VectorBase &factors)
  {
    Assert (size() == factors.size(),
	    ExcDimensionMismatch(size(), factors.size()));

    const int ierr = vector->Multiply (1.0, *(factors.vector), *vector, 0.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  VectorBase::equ (const TrilinosScalar  a,
		   const VectorBase     &v)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

				   // If we don't have the same map, copy.
    if (local_range() != v.local_range())
      {
	vector.reset();
	last_action = Zero;
	*vector = *v.vector;
	*this *= a;
      }
    else
      {
				   // Otherwise, just update
	Assert (vector->Map().SameAs(v.vector->Map()) == true,
		ExcMessage ("The Epetra maps in the assignment operator ="
			    " do not match, even though the local_range "
			    " seems to be the same. Check vector setup!"));
	int ierr;
	ierr = vector->GlobalAssemble(last_action);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	ierr = vector->Update(a, *v.vector, 0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	last_action = Zero;
      }

  }



  void
  VectorBase::equ (const TrilinosScalar  a,
		   const VectorBase     &v,
		   const TrilinosScalar  b,
		   const VectorBase     &w)
  {

    Assert (v.size() == w.size(),
	    ExcDimensionMismatch (v.size(), w.size()));

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(b),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

				   // If we don't have the same map, copy.
    if (local_range() != v.local_range())
      {
	vector.reset();
	last_action = Zero;
	*vector = *v.vector;
	sadd (a, b, w);
      }
    else
      {
				   // Otherwise, just update
	Assert (vector->Map().SameAs(v.vector->Map()) == true,
		ExcMessage ("The Epetra maps in the assignment operator ="
			    " do not match, even though the local_range "
			    " seems to be the same. Check vector setup!"));
	int ierr;
	ierr = vector->GlobalAssemble(last_action);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	ierr = vector->Update(a, *v.vector, b, *w.vector, 0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	last_action = Zero;
      }
  }



  void
  VectorBase::ratio (const VectorBase &v,
		     const VectorBase &w)
  {
    Assert (v.size() == w.size(),
	    ExcDimensionMismatch (v.size(), w.size()));

    Assert (size() == w.size(),
	    ExcDimensionMismatch (size(), w.size()));

    const int ierr = vector->ReciprocalMultiply(1.0, *(w.vector), 
						*(v.vector), 0.0);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }


  
				        // TODO: up to now only local
                                        // data printed out! Find a
                                        // way to neatly output
                                        // distributed data...
  void
  VectorBase::print (const char *format) const
  {
    Assert (vector->GlobalLength()!=0, ExcEmptyObject());

    for (unsigned int j=0; j<size(); ++j)
      {
	double t = (*vector)[0][j];

	if (format != 0)
	  std::printf (format, t);
	else
	  std::printf (" %5.2f", double(t));
      }
    std::printf ("\n");
  }



  void
  VectorBase::print (std::ostream      &out,
		     const unsigned int precision,
		     const bool         scientific,
		     const bool         across) const
  {
    AssertThrow (out, ExcIO());

                                        // get a representation of the
                                        // vector and loop over all
                                        // the elements TODO: up to
                                        // now only local data printed
                                        // out! Find a way to neatly
                                        // output distributed data...
    TrilinosScalar *val;
    int leading_dimension;
    int ierr = vector->ExtractView (&val, &leading_dimension);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    out.precision (precision);
    if (scientific)
      out.setf (std::ios::scientific, std::ios::floatfield);
    else
      out.setf (std::ios::fixed, std::ios::floatfield);

    if (across)
      for (unsigned int i=0; i<size(); ++i)
	out << static_cast<double>(val[i]) << ' ';
    else
      for (unsigned int i=0; i<size(); ++i)
	out << static_cast<double>(val[i]) << std::endl;
    out << std::endl;

                                        // restore the representation
                                        // of the vector
    AssertThrow (out, ExcIO());
  }



  void
  VectorBase::swap (VectorBase &v)
  {
                                        // Just swap the pointers to
                                        // the two Epetra vectors that
                                        // hold all the data.
    VectorBase *p_v = &v, *p_this = this;
    VectorBase* tmp = p_v;
    p_v = p_this;
    p_this = tmp;
  }



  unsigned int
  VectorBase::memory_consumption () const
  {
    AssertThrow(false, ExcNotImplemented() );
    return 0;
  }

} /* end of namespace TrilinosWrappers */


namespace TrilinosWrappers
{
#include "trilinos_vector_base.inst"
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
