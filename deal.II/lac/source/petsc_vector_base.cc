//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/petsc_vector_base.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_parallel_vector.h>

#include <cmath>

#ifdef DEAL_II_USE_PETSC


namespace PETScWrappers
{
  namespace internal
  {
    VectorReference::operator PetscScalar () const
    {
      Assert (index < vector.size(),
              ExcIndexRange (index, 0, vector.size()));
              
                                       // this is clumsy: there is no simple
                                       // way in PETSc to read an element from
                                       // a vector, i.e. there is no function
                                       // VecGetValue or so. The only way is
                                       // to obtain a pointer to a contiguous
                                       // representation of the vector and
                                       // read from it. Subsequently, the
                                       // vector representation has to be
                                       // restored. In addition, we can only
                                       // get access to the local part of the
                                       // vector, so we have to guard against
                                       // that
      if (dynamic_cast<const PETScWrappers::Vector *>(&vector) != 0)
        {
          PetscScalar *ptr;
          int ierr
            = VecGetArray (static_cast<const Vec &>(vector), &ptr);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          
          const PetscScalar value = *(ptr+index);
          
          ierr = VecRestoreArray (static_cast<const Vec &>(vector), &ptr);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          
          return value;
        }
      else if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&vector) != 0)
        {
                                           // first verify that the requested
                                           // element is actually locally
                                           // available
          int ierr;
          int begin, end;
          ierr = VecGetOwnershipRange (static_cast<const Vec &>(vector),
                                       &begin, &end);
          AssertThrow (ierr == 0, ExcPETScError(ierr));

          AssertThrow ((index >= static_cast<unsigned int>(begin)) &&
                       (index < static_cast<unsigned int>(end)),
                       ExcAccessToNonlocalElement (index, begin, end-1));

                                           // then access it
          PetscScalar *ptr;
          ierr = VecGetArray (static_cast<const Vec &>(vector), &ptr);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          
          const PetscScalar value = *(ptr+index-begin);
          
          ierr = VecRestoreArray (static_cast<const Vec &>(vector), &ptr);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          
          return value;
        }
      else
                                         // what? what other kind of vector
                                         // exists there?
        Assert (false, ExcInternalError());
      return -1e20;
    }  
  }
  
  VectorBase::VectorBase ()
                  :
                  last_action (LastAction::none)
  {}
  

  
  VectorBase::VectorBase (const VectorBase &v)
                  :
                  last_action (LastAction::none)
  {
    int ierr = VecDuplicate (v.vector, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = VecCopy (v.vector, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }

  
  VectorBase::~VectorBase ()
  {
    const int ierr = VecDestroy (vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  VectorBase &
  VectorBase::operator = (const PetscScalar s)
  {
                                     // flush previously cached elements. this
                                     // seems to be necessary since petsc
                                     // 2.2.1, at least for parallel vectors
                                     // (see test bits/petsc_65)
    compress ();
    
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecSet (&s, vector);
#else
    const int ierr = VecSet (vector, s);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  bool
  VectorBase::operator == (const VectorBase &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));
    
    PetscTruth flag;
    
    const int ierr = VecEqual (vector, v.vector, &flag);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_TRUE);
  }


  
  bool
  VectorBase::operator != (const VectorBase &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));

    PetscTruth flag;
    
    const int ierr = VecEqual (vector, v.vector, &flag);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_FALSE);
  }


  
  unsigned int
  VectorBase::size () const
  {
    int sz;
    const int ierr = VecGetSize (vector, &sz);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  unsigned int
  VectorBase::local_size () const
  {
    int sz;
    const int ierr = VecGetLocalSize (vector, &sz);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  std::pair<unsigned int, unsigned int>
  VectorBase::local_range () const
  {
    int begin, end;
    const int ierr = VecGetOwnershipRange (static_cast<const Vec &>(vector),
					   &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return std::make_pair (begin, end);
  }



  void
  VectorBase::set (const std::vector<unsigned int> &indices,
                   const std::vector<PetscScalar>  &values)
  {
    Assert (indices.size() == values.size(),
            ExcMessage ("Function called with arguments of different sizes"));

    if (last_action != VectorBase::LastAction::insert)
      {
        int ierr;
        ierr = VecAssemblyBegin (vector);
        AssertThrow (ierr == 0, ExcPETScError(ierr));

        ierr = VecAssemblyEnd (vector);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }

#if (PETSC_VERSION_MAJOR <= 2) && \
    ((PETSC_VERSION_MINOR < 2) ||  \
     ((PETSC_VERSION_MINOR == 2) && (PETSC_VERSION_SUBMINOR == 0)))
    const std::vector<int> petsc_indices (indices.begin(),
                                          indices.end());
#else
    const std::vector<PetscInt> petsc_indices (indices.begin(),
                                               indices.end());
#endif
    
    const int ierr
      = VecSetValues (vector, indices.size(),
                      &petsc_indices[0], &values[0],
                      INSERT_VALUES);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    

    last_action = LastAction::insert;
  }
  


  PetscScalar
  VectorBase::operator * (const VectorBase &vec) const
  {
    Assert (size() == vec.size(),
            ExcDimensionMismatch(size(), vec.size()));

    PetscScalar result;

    const int ierr = VecDot (vector, vec.vector, &result);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return result;
  }
  


  void
  VectorBase::compress ()
  {
    int ierr;
    ierr = VecAssemblyBegin(vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    ierr = VecAssemblyEnd(vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  PetscScalar
  VectorBase::norm_sqr () const
  {
    const PetscScalar d = l2_norm();
    return d*d;
  }



  PetscScalar
  VectorBase::mean_value () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    PetscScalar *start_ptr;
    int ierr = VecGetArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    PetscScalar mean = 0;
    {
      PetscScalar sum0 = 0,
                  sum1 = 0,
                  sum2 = 0,
                  sum3 = 0;

                                       // use modern processors better by
                                       // allowing pipelined commands to be
                                       // executed in parallel
      const PetscScalar * ptr  = start_ptr;
      const PetscScalar * eptr = ptr + (size()/4)*4;
      while (ptr!=eptr)
        {
          sum0 += *ptr++;
          sum1 += *ptr++;
          sum2 += *ptr++;
          sum3 += *ptr++;
        };
                                       // add up remaining elements
      while (ptr != start_ptr+size())
        sum0 += *ptr++;
  
      mean = (sum0+sum1+sum2+sum3)/size();
    }
    
                                     // restore the representation of the
                                     // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
      
    return mean;
  }
  


  PetscScalar
  VectorBase::l1_norm () const
  {
    PetscScalar d;

    const int ierr = VecNorm (vector, NORM_1, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }
  


  PetscScalar
  VectorBase::l2_norm () const
  {
    PetscScalar d;

    const int ierr = VecNorm (vector, NORM_2, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }
  


  PetscScalar
  VectorBase::lp_norm (const PetscScalar p) const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    PetscScalar *start_ptr;
    int ierr = VecGetArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    PetscScalar norm = 0;
    {
      PetscScalar sum0 = 0,
                  sum1 = 0,
                  sum2 = 0,
                  sum3 = 0;

                                       // use modern processors better by
                                       // allowing pipelined commands to be
                                       // executed in parallel
      const PetscScalar * ptr  = start_ptr;
      const PetscScalar * eptr = ptr + (size()/4)*4;
      while (ptr!=eptr)
        {
          sum0 += std::pow(std::fabs(*ptr++), p);
          sum1 += std::pow(std::fabs(*ptr++), p);
          sum2 += std::pow(std::fabs(*ptr++), p);
          sum3 += std::pow(std::fabs(*ptr++), p);
        };
                                       // add up remaining elements
      while (ptr != start_ptr+size())
        sum0 += std::pow(std::fabs(*ptr++), p);
  
      norm = std::pow(sum0+sum1+sum2+sum3,
                      static_cast<PetscScalar>(1./p));
    }
    
                                     // restore the representation of the
                                     // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
      
    return norm;
  }



  PetscScalar
  VectorBase::linfty_norm () const
  {
    PetscScalar d;

    const int ierr = VecNorm (vector, NORM_INFINITY, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }



  bool
  VectorBase::all_zero () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    PetscScalar *start_ptr;
    int ierr = VecGetArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    const PetscScalar *ptr  = start_ptr,
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
    
                                     // restore the representation of the
                                     // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
      
    return flag;
  }



  bool
  VectorBase::is_non_negative () const
  {
                                     // get a representation of the vector and
                                     // loop over all the elements
    PetscScalar *start_ptr;
    int ierr = VecGetArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    const PetscScalar *ptr  = start_ptr,
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
    
                                     // restore the representation of the
                                     // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
      
    return flag;
  }
  


  VectorBase &
  VectorBase::operator *= (const PetscScalar a)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecScale (&a, vector);
#else
    const int ierr = VecScale (vector, a);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator /= (const PetscScalar a)
  {
    const PetscScalar factor = 1./a;
    
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecScale (&factor, vector);
#else
    const int ierr = VecScale (vector, factor);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator += (const VectorBase &v)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const PetscScalar one = 1.0;    
    const int ierr = VecAXPY (&one, v, vector);
#else
    const int ierr = VecAXPY (vector, 1, v);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator -= (const VectorBase &v)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const PetscScalar minus_one = -1.0;    
    const int ierr = VecAXPY (&minus_one, v, vector);
#else
    const int ierr = VecAXPY (vector, -1, v);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  void
  VectorBase::add (const PetscScalar s)
  {
 #if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecShift (&s, vector);
#else
    const int ierr = VecShift (vector, s);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  void
  VectorBase::add (const VectorBase &v)
  {
    *this += v;
  }
  


  void
  VectorBase::add (const PetscScalar a,
                   const VectorBase     &v)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecAXPY (&a, v, vector);
#else
    const int ierr = VecAXPY (vector, a, v);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  void
  VectorBase::add (const PetscScalar a,
                   const VectorBase &v,
                   const PetscScalar b,
                   const VectorBase &w)
  {
    const PetscScalar weights[2] = {a,b};
    Vec               addends[2] = {v.vector, w.vector};
    
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecMAXPY (2, weights, vector, addends);
#else
    const int ierr = VecMAXPY (vector, 2, weights, addends);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  void
  VectorBase::sadd (const PetscScalar s,
                    const VectorBase &v)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecAYPX (&s, v, vector);
#else
    const int ierr = VecAYPX (vector, s, v);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  void
  VectorBase::sadd (const PetscScalar s,
                    const PetscScalar a,
                    const VectorBase     &v)
  {
                                     // there is nothing like a AXPAY
                                     // operation in Petsc, so do it in two
                                     // steps
    *this *= s;
    add (a,v);
  }
  


  void
  VectorBase::sadd (const PetscScalar s,
                    const PetscScalar a,
                    const VectorBase     &v,
                    const PetscScalar b,
                    const VectorBase     &w)
  {
                                     // there is no operation like MAXPAY, so
                                     // do it in two steps
    *this *= s;
    
    const PetscScalar weights[2] = {a,b};
    Vec               addends[2] = {v.vector,w.vector};
    
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecMAXPY (2, weights, vector, addends);
#else
    const int ierr = VecMAXPY (vector, 2, weights, addends);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::sadd (const PetscScalar s,
                    const PetscScalar a,
                    const VectorBase     &v,
                    const PetscScalar b,
                    const VectorBase     &w,
                    const PetscScalar c,
                    const VectorBase     &x)
  {
                                     // there is no operation like MAXPAY, so
                                     // do it in two steps
    *this *= s;

    const PetscScalar weights[3] = {a,b,c};
    Vec               addends[3] = {v.vector, w.vector, x.vector};
    
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecMAXPY (3, weights, vector, addends);
#else
    const int ierr = VecMAXPY (vector, 3, weights, addends);
#endif
    
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::scale (const VectorBase &factors)
  {
    const int ierr
      = VecPointwiseMult (vector, factors, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::equ (const PetscScalar a,
                   const VectorBase &v)
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));

                                     // there is no simple operation for this
                                     // in PETSc. there are multiple ways to
                                     // emulate it, we choose this one:
    const int ierr = VecCopy (v.vector, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    *this *= a;
  }
  


  void
  VectorBase::equ (const PetscScalar a,
                   const VectorBase &v,
                   const PetscScalar b,
                   const VectorBase &w)
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch (size(), v.size()));

                                     // there is no simple operation for this
                                     // in PETSc. there are multiple ways to
                                     // emulate it, we choose this one:
    const int ierr = VecCopy (v.vector, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    sadd (a, b, w);
  }



  void
  VectorBase::ratio (const VectorBase &a,
                     const VectorBase &b)
  {
#if (PETSC_VERSION_MAJOR <= 2) && (PETSC_VERSION_MINOR < 3)
    const int ierr = VecPointwiseDivide (a, b, vector);
#else
    const int ierr = VecPointwiseDivide (vector, a, b);
#endif

    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::print (std::ostream      &out,
                     const unsigned int precision,
                     const bool         scientific,
                     const bool         across) const
  {
    AssertThrow (out, ExcIO());

                                     // get a representation of the vector and
                                     // loop over all the elements
    PetscScalar *val;
    int ierr = VecGetArray (vector, &val);

    AssertThrow (ierr == 0, ExcPETScError(ierr));
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
    ierr = VecRestoreArray (vector, &val);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    AssertThrow (out, ExcIO());
  }


  
  void
  VectorBase::swap (VectorBase &v)
  {
    const int ierr = VecSwap (vector, v.vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }
  


  VectorBase::operator const Vec & () const
  {
    return vector;
  }
  
}

#else
// On gcc2.95 on Alpha OSF1, the native assembler does not like empty
// files, so provide some dummy code
namespace { void dummy () {} }
#endif // DEAL_II_USE_PETSC
