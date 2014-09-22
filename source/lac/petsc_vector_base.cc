// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/petsc_vector_base.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/memory_consumption.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#  include <cmath>
#  include <deal.II/base/multithread_info.h>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace internal
  {
    VectorReference::operator PetscScalar () const
    {
      Assert (index < vector.size(),
              ExcIndexRange (index, 0, vector.size()));

      // if the vector is local, then
      // simply access the element we
      // are interested in
      if (dynamic_cast<const PETScWrappers::Vector *>(&vector) != 0)
        {
          PetscInt idx = index;
          PetscScalar value;
          int ierr = VecGetValues(vector.vector, 1, &idx, &value);
          AssertThrow (ierr == 0, ExcPETScError(ierr));
          return value;
        }
      // else see if we are dealing
      // with a parallel vector
      else if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&vector) != 0)
        {
          int ierr;

          // there is the possibility
          // that the vector has
          // ghost elements. in that
          // case, we first need to
          // figure out which
          // elements we own locally,
          // then get a pointer to
          // the elements that are
          // stored here (both the
          // ones we own as well as
          // the ghost elements). in
          // this array, the locally
          // owned elements come
          // first followed by the
          // ghost elements whose
          // position we can get from
          // an index set
          if (vector.ghosted)
            {
              PetscInt begin, end;
              ierr = VecGetOwnershipRange (vector.vector, &begin, &end);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              Vec locally_stored_elements = PETSC_NULL;
              ierr = VecGhostGetLocalForm(vector.vector, &locally_stored_elements);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              PetscInt lsize;
              ierr = VecGetSize(locally_stored_elements, &lsize);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              PetscScalar *ptr;
              ierr = VecGetArray(locally_stored_elements, &ptr);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              PetscScalar value;

              if ( index>=static_cast<size_type>(begin)
                   && index<static_cast<size_type>(end) )
                {
                  //local entry
                  value = *(ptr+index-begin);
                }
              else
                {
                  //ghost entry
                  const size_type ghostidx
                    = vector.ghost_indices.index_within_set(index);

                  Assert(ghostidx+end-begin<(size_type)lsize, ExcInternalError());
                  value = *(ptr+ghostidx+end-begin);
                }

              ierr = VecRestoreArray(locally_stored_elements, &ptr);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              ierr = VecGhostRestoreLocalForm(vector.vector, &locally_stored_elements);
              AssertThrow (ierr == 0, ExcPETScError(ierr));

              return value;
            }


          // first verify that the requested
          // element is actually locally
          // available
          PetscInt begin, end;

          ierr = VecGetOwnershipRange (vector.vector, &begin, &end);
          AssertThrow (ierr == 0, ExcPETScError(ierr));



          AssertThrow ((index >= static_cast<size_type>(begin)) &&
                       (index < static_cast<size_type>(end)),
                       ExcAccessToNonlocalElement (index, begin, end-1));

          // old version which only work with
          // VecGetArray()...
          PetscInt idx = index;
          PetscScalar value;
          ierr = VecGetValues(vector.vector, 1, &idx, &value);
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
    ghosted(false),
    last_action (::dealii::VectorOperation::unknown),
    attained_ownership(true)
  {
    Assert( multithread_info.is_running_single_threaded(),
            ExcMessage("PETSc does not support multi-threaded access, set "
                       "the thread limit to 1 in MPI_InitFinalize()."));
  }



  VectorBase::VectorBase (const VectorBase &v)
    :
    Subscriptor (),
    ghosted(v.ghosted),
    ghost_indices(v.ghost_indices),
    last_action (::dealii::VectorOperation::unknown),
    attained_ownership(true)
  {
    Assert( multithread_info.is_running_single_threaded(),
            ExcMessage("PETSc does not support multi-threaded access, set "
                       "the thread limit to 1 in MPI_InitFinalize()."));

    int ierr = VecDuplicate (v.vector, &vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    ierr = VecCopy (v.vector, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  VectorBase::VectorBase (const Vec &v)
    :
    Subscriptor (),
    vector(v),
    ghosted(false),
    last_action (::dealii::VectorOperation::unknown),
    attained_ownership(false)
  {
    Assert( multithread_info.is_running_single_threaded(),
            ExcMessage("PETSc does not support multi-threaded access, set "
                       "the thread limit to 1 in MPI_InitFinalize()."));
  }



  VectorBase::~VectorBase ()
  {
    if (attained_ownership)
      {
#if DEAL_II_PETSC_VERSION_LT(3,2,0)
        const int ierr = VecDestroy (vector);
#else
        const int ierr = VecDestroy (&vector);
#endif
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }
  }



  VectorBase &
  VectorBase::operator = (const PetscScalar s)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    //TODO[TH]: assert(is_compressed())

    const int ierr = VecSet (vector, s);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  bool
  VectorBase::operator == (const VectorBase &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    PetscTruth
#else
    PetscBool
#endif
    flag;

    const int ierr = VecEqual (vector, v.vector, &flag);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_TRUE);
  }



  bool
  VectorBase::operator != (const VectorBase &v) const
  {
    Assert (size() == v.size(),
            ExcDimensionMismatch(size(), v.size()));

#if DEAL_II_PETSC_VERSION_LT(3,2,0)
    PetscTruth
#else
    PetscBool
#endif
    flag;

    const int ierr = VecEqual (vector, v.vector, &flag);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_FALSE);
  }



  VectorBase::size_type
  VectorBase::size () const
  {
    PetscInt sz;
    const int ierr = VecGetSize (vector, &sz);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  VectorBase::size_type
  VectorBase::local_size () const
  {
    PetscInt sz;
    const int ierr = VecGetLocalSize (vector, &sz);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  std::pair<VectorBase::size_type, VectorBase::size_type>
  VectorBase::local_range () const
  {
    PetscInt begin, end;
    const int ierr = VecGetOwnershipRange (static_cast<const Vec &>(vector),
                                           &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return std::make_pair (begin, end);
  }



  void
  VectorBase::set (const std::vector<size_type> &indices,
                   const std::vector<PetscScalar>  &values)
  {
    Assert (indices.size() == values.size(),
            ExcMessage ("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), &indices[0], &values[0], false);
  }



  void
  VectorBase::add (const std::vector<size_type> &indices,
                   const std::vector<PetscScalar>  &values)
  {
    Assert (indices.size() == values.size(),
            ExcMessage ("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), &indices[0], &values[0], true);
  }



  void
  VectorBase::add (const std::vector<size_type>    &indices,
                   const ::dealii::Vector<PetscScalar> &values)
  {
    Assert (indices.size() == values.size(),
            ExcMessage ("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), &indices[0], values.begin(), true);
  }



  void
  VectorBase::add (const size_type    n_elements,
                   const size_type   *indices,
                   const PetscScalar *values)
  {
    do_set_add_operation(n_elements, indices, values, true);
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
  VectorBase::compress (::dealii::VectorOperation::values)
  {
    // note that one may think that
    // we only need to do something
    // if in fact the state is
    // anything but
    // last_action::unknown. but
    // that's not true: one
    // frequently gets into
    // situations where only one
    // processor (or a subset of
    // processors) actually writes
    // something into a vector, but
    // we still need to call
    // VecAssemblyBegin/End on all
    // processors.
    int ierr;
    ierr = VecAssemblyBegin(vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    ierr = VecAssemblyEnd(vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    // reset the last action field to
    // indicate that we're back to a
    // pristine state
    last_action = ::dealii::VectorOperation::unknown;
  }



  void
  VectorBase::compress ()
  {
    compress(VectorOperation::unknown);
  }



  VectorBase::real_type
  VectorBase::norm_sqr () const
  {
    const real_type d = l2_norm();
    return d*d;
  }



  PetscScalar
  VectorBase::mean_value () const
  {
    int ierr;

    // We can only use our more efficient
    // routine in the serial case.
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(this) != 0)
      {
        PetscScalar sum;
        ierr = VecSum(vector, &sum);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
        return sum/static_cast<PetscReal>(size());
      }

    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *start_ptr;
    ierr = VecGetArray (vector, &start_ptr);
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
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (size()/4)*4;
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

      mean = (sum0+sum1+sum2+sum3)/static_cast<PetscReal>(size());
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return mean;
  }


  VectorBase::real_type
  VectorBase::l1_norm () const
  {
    real_type d;

    const int ierr = VecNorm (vector, NORM_1, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::l2_norm () const
  {
    real_type d;

    const int ierr = VecNorm (vector, NORM_2, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::lp_norm (const real_type p) const
  {
    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *start_ptr;
    int ierr = VecGetArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    real_type norm = 0;
    {
      real_type sum0 = 0,
                sum1 = 0,
                sum2 = 0,
                sum3 = 0;

      // use modern processors better by
      // allowing pipelined commands to be
      // executed in parallel
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (size()/4)*4;
      while (ptr!=eptr)
        {
          sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum1 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum2 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum3 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
        }
      // add up remaining elements
      while (ptr != start_ptr+size())
        sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);

      norm = std::pow(sum0+sum1+sum2+sum3,
                      1./p);
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArray (vector, &start_ptr);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return norm;
  }



  VectorBase::real_type
  VectorBase::linfty_norm () const
  {
    real_type d;

    const int ierr = VecNorm (vector, NORM_INFINITY, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::normalize () const
  {
    real_type d;
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecNormalize (vector, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }


  VectorBase::real_type
  VectorBase::min ()  const
  {
    PetscInt  p;
    real_type d;

    const int ierr = VecMin (vector, &p, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }


  VectorBase::real_type
  VectorBase::max ()  const
  {
    PetscInt  p;
    real_type d;

    const int ierr = VecMax (vector, &p, &d);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return d;
  }


  VectorBase &
  VectorBase::abs ()
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    const int ierr = VecAbs (vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::conjugate ()
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    const int ierr = VecConjugate (vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::mult ()
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());

    const int ierr = VecPointwiseMult (vector,vector,vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }


  VectorBase &
  VectorBase::mult (const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecPointwiseMult (vector,vector,v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }


  VectorBase &
  VectorBase::mult (const VectorBase &u,
                    const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecPointwiseMult (vector,u,v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
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
                       *eptr = start_ptr + local_size();
    bool flag = true;
    while (ptr != eptr)
      {
        if (*ptr != value_type())
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


  namespace internal
  {
    template <typename T>
    bool is_non_negative (const T &t)
    {
      return t >= 0;
    }



    template <typename T>
    bool is_non_negative (const std::complex<T> &)
    {
      Assert (false,
              ExcMessage ("You can't ask a complex value "
                          "whether it is non-negative."))
      return true;
    }
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
                       *eptr = start_ptr + local_size();
    bool flag = true;
    while (ptr != eptr)
      {
        if (! internal::is_non_negative(*ptr))
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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const int ierr = VecScale (vector, a);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator /= (const PetscScalar a)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const PetscScalar factor = 1./a;
    Assert (numbers::is_finite(factor), ExcNumberNotFinite());

    const int ierr = VecScale (vector, factor);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator += (const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecAXPY (vector, 1, v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator -= (const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecAXPY (vector, -1, v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  void
  VectorBase::add (const PetscScalar s)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    const int ierr = VecShift (vector, s);
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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

    const int ierr = VecAXPY (vector, a, v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::add (const PetscScalar a,
                   const VectorBase &v,
                   const PetscScalar b,
                   const VectorBase &w)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    const PetscScalar weights[2] = {a,b};
    Vec               addends[2] = {v.vector, w.vector};

    const int ierr = VecMAXPY (vector, 2, weights, addends);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::sadd (const PetscScalar s,
                    const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());

    const int ierr = VecAYPX (vector, s, v);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::sadd (const PetscScalar s,
                    const PetscScalar a,
                    const VectorBase     &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

    // there is no operation like MAXPAY, so
    // do it in two steps
    *this *= s;

    const PetscScalar weights[2] = {a,b};
    Vec               addends[2] = {v.vector,w.vector};

    const int ierr = VecMAXPY (vector, 2, weights, addends);
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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(s), ExcNumberNotFinite());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());
    Assert (numbers::is_finite(c), ExcNumberNotFinite());

    // there is no operation like MAXPAY, so
    // do it in two steps
    *this *= s;

    const PetscScalar weights[3] = {a,b,c};
    Vec               addends[3] = {v.vector, w.vector, x.vector};

    const int ierr = VecMAXPY (vector, 3, weights, addends);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::scale (const VectorBase &factors)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr
      = VecPointwiseMult (vector, factors, vector);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::equ (const PetscScalar a,
                   const VectorBase &v)
  {
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());

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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    Assert (numbers::is_finite(a), ExcNumberNotFinite());
    Assert (numbers::is_finite(b), ExcNumberNotFinite());

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
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    const int ierr = VecPointwiseDivide (vector, a, b);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::write_ascii (const PetscViewerFormat format)
  {
    //TODO[TH]:assert(is_compressed())

    // Set options
    PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                          format);

    // Write to screen
    VecView (vector, PETSC_VIEWER_STDOUT_WORLD);
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

    // save the state of out stream
    const std::ios::fmtflags old_flags = out.flags();
    const unsigned int old_precision = out.precision (precision);

    out.precision (precision);
    if (scientific)
      out.setf (std::ios::scientific, std::ios::floatfield);
    else
      out.setf (std::ios::fixed, std::ios::floatfield);

    if (across)
      for (size_type i=0; i<local_size(); ++i)
        out << val[i] << ' ';
    else
      for (size_type i=0; i<local_size(); ++i)
        out << val[i] << std::endl;
    out << std::endl;

    // reset output format
    out.flags (old_flags);
    out.precision(old_precision);

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



  VectorBase::operator const Vec &() const
  {
    return vector;
  }


  std::size_t
  VectorBase::memory_consumption () const
  {
    std::size_t mem = sizeof(Vec)+sizeof(last_action)
                      +MemoryConsumption::memory_consumption(ghosted)
                      +MemoryConsumption::memory_consumption(ghost_indices);

    // TH: I am relatively sure that PETSc is
    // storing the local data in a contiguous
    // block without indices:
    mem += local_size()*sizeof(PetscScalar);
    // assume that PETSc is storing one index
    // and one double per ghost element
    if (ghosted)
      mem += ghost_indices.n_elements()*(sizeof(PetscScalar)+sizeof(int));

    //TODO[TH]: size of constant memory for PETSc?
    return mem;
  }



  void
  VectorBase::do_set_add_operation (const size_type    n_elements,
                                    const size_type   *indices,
                                    const PetscScalar *values,
                                    const bool         add_values)
  {
    ::dealii::VectorOperation::values action = (add_values ?
                                                ::dealii::VectorOperation::add :
                                                ::dealii::VectorOperation::insert);
    Assert ((last_action == action)
            ||
            (last_action == ::dealii::VectorOperation::unknown),
            internal::VectorReference::ExcWrongMode (action,
                                                     last_action));
    Assert (!has_ghost_elements(), ExcGhostsPresent());
    // VecSetValues complains if we
    // come with an empty
    // vector. however, it is not a
    // collective operation, so we
    // can skip the call if necessary
    // (unlike the above calls)
    if (n_elements != 0)
      {
#ifdef PETSC_USE_64BIT_INDICES
        std::vector<PetscInt> petsc_ind (n_elements);
        for (size_type i=0; i<n_elements; ++i)
          petsc_ind[i] = indices[i];
        const PetscInt *petsc_indices = &petsc_ind[0];
#else
        const int *petsc_indices = (const int *)indices;
#endif

        const InsertMode mode = (add_values ? ADD_VALUES : INSERT_VALUES);
        const int ierr
          = VecSetValues (vector, n_elements, petsc_indices, values,
                          mode);
        AssertThrow (ierr == 0, ExcPETScError(ierr));
      }

    // set the mode here, independent of whether we have actually
    // written elements or whether the list was empty
    last_action = action;
  }


  void
  VectorBase::update_ghost_values() const
  {
    // generate an error for not ghosted
    // vectors
    if (!ghosted)
      AssertThrow (false, ExcInternalError());

    int ierr;

    ierr = VecGhostUpdateBegin(vector, INSERT_VALUES, SCATTER_FORWARD);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    ierr = VecGhostUpdateEnd(vector, INSERT_VALUES, SCATTER_FORWARD);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
  }



}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
