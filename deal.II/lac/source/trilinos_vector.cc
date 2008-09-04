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


#include <lac/trilinos_vector.h>

#include <cmath>

#ifdef DEAL_II_USE_TRILINOS

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace internal
  {
    VectorReference::operator TrilinosScalar () const // i believe useless with
                                                      // trilinos
    {
      Assert (index < vector.size(),
              ExcIndexRange (index, 0, vector.size()));

                                       // Trilinos allows for vectors to be
                                       // referenced by the [] or () operators
                                       // but only () checks index bounds
                                       // Also, can only get local values

      AssertThrow ((static_cast<signed int>(index) >= vector.map.MinMyGID()) &&
		   (static_cast<signed int>(index) <= vector.map.MaxMyGID()),
		   ExcAccessToNonLocalElement (index, vector.map.MinMyGID(),
					       vector.map.MaxMyGID()-1));

      return (*(vector.vector))[0][index];
    }
  }

  Vector::Vector ()
                  :
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
		  map (0,0,Epetra_MpiComm(MPI_COMM_WORLD)),
#else
		  map (0,0,Epetra_SerialComm()),
#endif
		  vector(std::auto_ptr<Epetra_FEVector> 
			 (new Epetra_FEVector(map))),
                  last_action (Insert)
  {}

  Vector::Vector (unsigned int GlobalSize, Epetra_Comm &Comm)
                  :
		  map (GlobalSize, 0, Comm),
		  vector (std::auto_ptr<Epetra_FEVector> 
			  (new Epetra_FEVector(map))),
                  last_action (Insert)
  {}

  
  Vector::Vector (const Epetra_Map &InputMap)
                  :
		  map (InputMap),
		  vector (std::auto_ptr<Epetra_FEVector> 
			  (new Epetra_FEVector(map))),
                  last_action (Insert)
  {}
  
  
  Vector::Vector (const Vector &v,
		  const bool    fast)
                  :
		  map (v.map),
		  vector(std::auto_ptr<Epetra_FEVector> 
			 (new Epetra_FEVector(v.map,!fast))),
                  last_action (Insert)
  {}
  


  Vector::~Vector ()
  {}



  void
  Vector::reinit (const Epetra_Map &input_map)
  {
    vector.reset();
    map = input_map;

    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(input_map));
    last_action = Insert;
  }



  void
  Vector::reinit (const Vector &v,
		  const bool    fast)
  {
    vector.reset();

    if (!map.SameAs(v.map))
      map = v.map;

    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(v.map,!fast));
    last_action = Insert;
  }



  void
  Vector::clear ()
  {
                                     // When we clear the matrix,
				     // reset the pointer and 
				     // generate an empty matrix.
    vector.reset();
#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
    map = Epetra_Map (0,0,Epetra_MpiComm(MPI_COMM_WORLD)),
#else
    map = Epetra_Map (0,0,Epetra_SerialComm()),
#endif

    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
  }



  void
  Vector::compress ()
  {
				 // Now pass over the information
				 // about what we did last to Trilinos.
    const int ierr = vector->GlobalAssemble(last_action);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  Vector &
  Vector::operator = (const TrilinosScalar s)
  {

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but either "
		       "infinite or Not A Number (NaN)"));

    const int ierr = vector->PutScalar(s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  bool
  Vector::operator == (const Vector &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));
    unsigned int i;
    for (i=0; i<size(); i++) 
      if ((*(v.vector))[0][i]!=(*vector)[0][i]) return false;

    return true;
  }



  bool
  Vector::operator != (const Vector &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));

    return (!(*this==v));
  }



  unsigned int
  Vector::size () const
  {
    return (unsigned int) vector->Map().NumGlobalElements();
  }



  unsigned int
  Vector::local_size () const
  {
    return (unsigned int) vector->Map().NumMyElements();
  }



  std::pair<unsigned int, unsigned int>
  Vector::local_range () const
  {
    int begin, end;
    begin = vector->Map().MinMyGID();
    end = vector->Map().MaxMyGID()+1;
    return std::make_pair (begin, end);
  }



  TrilinosScalar
  Vector::el (const unsigned int index) const
  {
                                      // Extract local indices in
                                      // the vector.
    int trilinos_i = map.LID(index);
    TrilinosScalar value = 0.;
    if (trilinos_i == -1 )
      {
        Assert (false, ExcAccessToNonlocalElement(index, local_range().first,
						  local_range().second));
      }
    else
      value = (*vector)[0][trilinos_i];

    return value;
  }



  void
  Vector::set (const std::vector<unsigned int>    &indices,
	       const std::vector<TrilinosScalar>  &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], &values[0]);
  }



  void
  Vector::set (const std::vector<unsigned int>        &indices,
	       const ::dealii::Vector<TrilinosScalar> &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    set (indices.size(), &indices[0], values.begin());
  }



  void
  Vector::set (const unsigned int    n_elements,
	       const unsigned int   *indices,
	       const TrilinosScalar *values)
  {
    if (last_action == Add)
      {
	vector->GlobalAssemble(Add);
	last_action = Insert;
      }

    const int ierr= vector->ReplaceGlobalValues (n_elements, 
				      (int*)(const_cast<unsigned int*>(indices)), 
				      const_cast<TrilinosScalar*>(values));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  Vector::add (const std::vector<unsigned int>    &indices,
	       const std::vector<TrilinosScalar>  &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    add (indices.size(), &indices[0], &values[0]);
  }



  void
  Vector::add (const std::vector<unsigned int>        &indices,
	       const ::dealii::Vector<TrilinosScalar> &values)
  {
    Assert (indices.size() == values.size(),
	    ExcDimensionMismatch(indices.size(),values.size()));

    add (indices.size(), &indices[0], values.begin());
  }



  void
  Vector::add (const unsigned int    n_elements,
	       const unsigned int   *indices,
	       const TrilinosScalar *values)
  {
    if (last_action == Insert)
      {
	vector->GlobalAssemble(Insert);
	last_action = Add;
      }

    const int ierr= vector->SumIntoGlobalValues (n_elements, 
				    (int*)(const_cast<unsigned int*>(indices)), 
				    const_cast<TrilinosScalar*>(values));
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  TrilinosScalar
  Vector::operator * (const Vector &vec) const
  {
    Assert (size() == vec.size(),
            ExcDimensionMismatch(size(), vec.size()));

    TrilinosScalar result;

    const int ierr = vector->Dot(*(vec.vector), &result);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return result;
  }



  Vector::real_type
  Vector::norm_sqr () const
  {
    const TrilinosScalar d = l2_norm();
    return d*d;
  }



  TrilinosScalar
  Vector::mean_value () const
  {
    TrilinosScalar mean;

    const int ierr = vector->MeanValue (&mean);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return mean;
  }



  Vector::real_type
  Vector::l1_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->Norm1 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  Vector::real_type
  Vector::l2_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->Norm2 (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  Vector::real_type
  Vector::lp_norm (const TrilinosScalar p) const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    TrilinosScalar norm = 0;
    TrilinosScalar sum=0;

    const TrilinosScalar * ptr  = start_ptr;
                                       // add up elements
    while (ptr != start_ptr+size())
      sum += std::pow(std::fabs(*ptr++), p);

    norm = std::pow(sum, static_cast<TrilinosScalar>(1./p));

    return norm;
  }



  Vector::real_type
  Vector::linfty_norm () const
  {
    TrilinosScalar d;

    const int ierr = vector->NormInf (&d);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return d;
  }



  bool
  Vector::all_zero () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

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
  Vector::is_non_negative () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    TrilinosScalar *start_ptr;
    int leading_dimension;
    int ierr = vector->ExtractView (&start_ptr, &leading_dimension);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

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



  Vector &
  Vector::operator *= (const TrilinosScalar a)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Scale(a);
    Assert (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  Vector &
  Vector::operator /= (const TrilinosScalar a)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const TrilinosScalar factor = 1./a;

    Assert (numbers::is_finite(factor),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Scale(factor);
    Assert (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  Vector &
  Vector::operator += (const Vector &v)
  {
    const int ierr = vector->Update (1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
    return *this;
  }



  Vector &
  Vector::operator -= (const Vector &v)
  {
    const int ierr = vector->Update (-1.0, *(v.vector), 1.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    return *this;
  }



  void
  Vector::add (const TrilinosScalar s)
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
  Vector::add (const Vector &v)
  {
    *this += v;
  }



  void
  Vector::add (const TrilinosScalar a,
	       const Vector        &v)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(a, *(v.vector), 1.);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  Vector::add (const TrilinosScalar a,
	       const Vector        &v,
	       const TrilinosScalar b,
	       const Vector        &w)
  {

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
  Vector::sadd (const TrilinosScalar s,
		const Vector        &v)
  {

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    const int ierr = vector->Update(1., *(v.vector), s);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  Vector::sadd (const TrilinosScalar s,
		const TrilinosScalar a,
		const Vector        &v)
  {

    Assert (numbers::is_finite(s),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

     const int ierr = vector->Update(a, *(v.vector), s);

     AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }



  void
  Vector::sadd (const TrilinosScalar s,
		const TrilinosScalar a,
		const Vector        &v,
		const TrilinosScalar b,
		const Vector        &w)
  {

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
  Vector::sadd (const TrilinosScalar s,
		const TrilinosScalar a,
		const Vector        &v,
		const TrilinosScalar b,
		const Vector        &w,
		const TrilinosScalar c,
		const Vector        &x)
  {

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

                                     // Update member can only input 
				     // two other vectors so
                                     // do it in two steps
    const int ierr = vector->Update(a, *(v.vector), b, *(w.vector), s);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

    const int jerr = vector->Update(c, *(x.vector), 1.);
    AssertThrow (jerr == 0, ExcTrilinosError(jerr));
  }



  void
  Vector::scale (const Vector &factors)
  {
    const int ierr = vector->Multiply (1.0, *(factors.vector), *vector, 0.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }



  void
  Vector::equ (const TrilinosScalar a,
	       const Vector        &v)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    *vector = *v.vector;
    map = v.map;
    
    *this *= a;
  }



  void
  Vector::equ (const TrilinosScalar a,
	       const Vector        &v,
	       const TrilinosScalar b,
	       const Vector        &w)
  {

    Assert (numbers::is_finite(a),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));
    Assert (numbers::is_finite(b),
	    ExcMessage("The given value is not finite but "
		       "either infinite or Not A Number (NaN)"));

    Assert (v.size() == w.size(),
            ExcDimensionMismatch (v.size(), w.size()));

    *vector = *v.vector;
    map = v.map;
    sadd (a, b, w);
  }



  void
  Vector::ratio (const Vector &v,
		 const Vector &w)
  {
    Assert (v.size() == w.size(),
            ExcDimensionMismatch (v.size(), w.size()));

    Assert (size() == w.size(),
            ExcDimensionMismatch (size(), w.size()));

    const int ierr = vector->ReciprocalMultiply(1.0, *(w.vector), 
						*(v.vector), 0.0);

    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
  }


  
				     // TODO: up to now only local data
                                     // printed out! Find a way to neatly
				     // output distributed data...
  void
  Vector::print (const char *format) const
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
  Vector::print (std::ostream      &out,
		 const unsigned int precision,
		 const bool         scientific,
		 const bool         across) const
  {
    AssertThrow (out, ExcIO());

                                     // get a representation of the vector and
                                     // loop over all the elements 
                                     // TODO: up to now only local data
                                     // printed out! Find a way to neatly
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

                                     // restore the representation of the
                                     // vector
    AssertThrow (out, ExcIO());
  }



  void
  Vector::swap (Vector &v)
  {
                                    // Just swap the pointers to the 
                                    // two Epetra vectors that hold all
                                    // the data.
    Vector *p_v = &v, *p_this = this;
    Vector* tmp = p_v;
    p_v = p_this;
    p_this = tmp;
  }


  unsigned int
  Vector::memory_consumption () const
  {
    AssertThrow(false, ExcNotImplemented() );
    return 0;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
