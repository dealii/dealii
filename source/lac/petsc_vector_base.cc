// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/petsc_vector_base.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/memory_consumption.h>
#  include <deal.II/base/multithread_info.h>

#  include <deal.II/lac/exceptions.h>
#  include <deal.II/lac/petsc_compatibility.h>
#  include <deal.II/lac/petsc_vector.h>

#  include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  namespace internal
  {
#  ifndef DOXYGEN
    VectorReference::operator PetscScalar() const
    {
      AssertIndexRange(index, vector.size());

      // The vector may have ghost entries. In that case, we first need to
      // figure out which elements we own locally, then get a pointer to the
      // elements that are stored here (both the ones we own as well as the
      // ghost elements). In this array, the locally owned elements come first
      // followed by the ghost elements whose position we can get from an
      // index set.
      if (vector.ghosted)
        {
          PetscInt       begin, end;
          PetscErrorCode ierr =
            VecGetOwnershipRange(vector.vector, &begin, &end);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          Vec locally_stored_elements = PETSC_NULL;
          ierr = VecGhostGetLocalForm(vector.vector, &locally_stored_elements);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          PetscInt lsize;
          ierr = VecGetSize(locally_stored_elements, &lsize);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          PetscScalar *ptr;
          ierr = VecGetArray(locally_stored_elements, &ptr);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          PetscScalar value;

          if (index >= static_cast<size_type>(begin) &&
              index < static_cast<size_type>(end))
            {
              // local entry
              value = *(ptr + index - begin);
            }
          else
            {
              // ghost entry
              const size_type ghostidx =
                vector.ghost_indices.index_within_set(index);

              AssertIndexRange(ghostidx + end - begin, lsize);
              value = *(ptr + ghostidx + end - begin);
            }

          ierr = VecRestoreArray(locally_stored_elements, &ptr);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          ierr =
            VecGhostRestoreLocalForm(vector.vector, &locally_stored_elements);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          return value;
        }


      // first verify that the requested
      // element is actually locally
      // available
      PetscInt begin, end;

      PetscErrorCode ierr = VecGetOwnershipRange(vector.vector, &begin, &end);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      AssertThrow((index >= static_cast<size_type>(begin)) &&
                    (index < static_cast<size_type>(end)),
                  ExcAccessToNonlocalElement(index, begin, end - 1));

      PetscInt    idx = index;
      PetscScalar value;
      ierr = VecGetValues(vector.vector, 1, &idx, &value);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return value;
    }
#  endif
  } // namespace internal

  VectorBase::VectorBase()
    : vector(nullptr)
    , ghosted(false)
    , last_action(::dealii::VectorOperation::unknown)
    , obtained_ownership(true)
  {
    Assert(MultithreadInfo::is_running_single_threaded(),
           ExcMessage("PETSc does not support multi-threaded access, set "
                      "the thread limit to 1 in MPI_InitFinalize()."));
  }



  VectorBase::VectorBase(const VectorBase &v)
    : Subscriptor()
    , ghosted(v.ghosted)
    , ghost_indices(v.ghost_indices)
    , last_action(::dealii::VectorOperation::unknown)
    , obtained_ownership(true)
  {
    Assert(MultithreadInfo::is_running_single_threaded(),
           ExcMessage("PETSc does not support multi-threaded access, set "
                      "the thread limit to 1 in MPI_InitFinalize()."));

    PetscErrorCode ierr = VecDuplicate(v.vector, &vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = VecCopy(v.vector, vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  VectorBase::VectorBase(const Vec &v)
    : Subscriptor()
    , vector(v)
    , ghosted(false)
    , last_action(::dealii::VectorOperation::unknown)
    , obtained_ownership(false)
  {
    Assert(MultithreadInfo::is_running_single_threaded(),
           ExcMessage("PETSc does not support multi-threaded access, set "
                      "the thread limit to 1 in MPI_InitFinalize()."));
  }



  VectorBase::~VectorBase()
  {
    if (obtained_ownership)
      {
        const PetscErrorCode ierr = VecDestroy(&vector);
        AssertNothrow(ierr == 0, ExcPETScError(ierr));
        (void)ierr;
      }
  }



  void
  VectorBase::clear()
  {
    if (obtained_ownership)
      {
        const PetscErrorCode ierr = VecDestroy(&vector);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    ghosted = false;
    ghost_indices.clear();
    last_action        = ::dealii::VectorOperation::unknown;
    obtained_ownership = true;
  }



  VectorBase &
  VectorBase::operator=(const PetscScalar s)
  {
    AssertIsFinite(s);

    // TODO[TH]: assert(is_compressed())

    PetscErrorCode ierr = VecSet(vector, s);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    if (has_ghost_elements())
      {
        Vec ghost = PETSC_NULL;
        ierr      = VecGhostGetLocalForm(vector, &ghost);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecSet(ghost, s);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecGhostRestoreLocalForm(vector, &ghost);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    return *this;
  }



  bool
  VectorBase::operator==(const VectorBase &v) const
  {
    Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

    PetscBool            flag;
    const PetscErrorCode ierr = VecEqual(vector, v.vector, &flag);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_TRUE);
  }



  bool
  VectorBase::operator!=(const VectorBase &v) const
  {
    Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

    PetscBool            flag;
    const PetscErrorCode ierr = VecEqual(vector, v.vector, &flag);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return (flag == PETSC_FALSE);
  }



  VectorBase::size_type
  VectorBase::size() const
  {
    PetscInt             sz;
    const PetscErrorCode ierr = VecGetSize(vector, &sz);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  VectorBase::size_type
  VectorBase::local_size() const
  {
    PetscInt             sz;
    const PetscErrorCode ierr = VecGetLocalSize(vector, &sz);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return sz;
  }



  std::pair<VectorBase::size_type, VectorBase::size_type>
  VectorBase::local_range() const
  {
    PetscInt             begin, end;
    const PetscErrorCode ierr =
      VecGetOwnershipRange(static_cast<const Vec &>(vector), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return std::make_pair(begin, end);
  }



  void
  VectorBase::set(const std::vector<size_type> &  indices,
                  const std::vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.data(), false);
  }



  void
  VectorBase::add(const std::vector<size_type> &  indices,
                  const std::vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.data(), true);
  }



  void
  VectorBase::add(const std::vector<size_type> &       indices,
                  const ::dealii::Vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.begin(), true);
  }



  void
  VectorBase::add(const size_type    n_elements,
                  const size_type *  indices,
                  const PetscScalar *values)
  {
    do_set_add_operation(n_elements, indices, values, true);
  }



  PetscScalar VectorBase::operator*(const VectorBase &vec) const
  {
    Assert(size() == vec.size(), ExcDimensionMismatch(size(), vec.size()));

    PetscScalar result;

    // For complex vectors, VecDot() computes
    //    val = (x,y) = y^H x,
    // where y^H denotes the conjugate transpose of y.
    // Note that this corresponds to the usual "mathematicians'"
    // complex inner product where the SECOND argument gets the
    // complex conjugate, which is also how we document this
    // operation.
    const PetscErrorCode ierr = VecDot(vec.vector, vector, &result);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return result;
  }



  PetscScalar
  VectorBase::add_and_dot(const PetscScalar a,
                          const VectorBase &V,
                          const VectorBase &W)
  {
    this->add(a, V);
    return *this * W;
  }



  void
  VectorBase::compress(const VectorOperation::values operation)
  {
    {
#  ifdef DEBUG
#    ifdef DEAL_II_WITH_MPI
      // Check that all processors agree that last_action is the same (or none!)

      int my_int_last_action = last_action;
      int all_int_last_action;

      const int ierr = MPI_Allreduce(&my_int_last_action,
                                     &all_int_last_action,
                                     1,
                                     MPI_INT,
                                     MPI_BOR,
                                     get_mpi_communicator());
      AssertThrowMPI(ierr);

      AssertThrow(all_int_last_action != (::dealii::VectorOperation::add |
                                          ::dealii::VectorOperation::insert),
                  ExcMessage("Error: not all processors agree on the last "
                             "VectorOperation before this compress() call."));
#    endif
#  endif
    }

    AssertThrow(
      last_action == ::dealii::VectorOperation::unknown ||
        last_action == operation,
      ExcMessage(
        "Missing compress() or calling with wrong VectorOperation argument."));

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
    PetscErrorCode ierr = VecAssemblyBegin(vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = VecAssemblyEnd(vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // reset the last action field to
    // indicate that we're back to a
    // pristine state
    last_action = ::dealii::VectorOperation::unknown;
  }



  VectorBase::real_type
  VectorBase::norm_sqr() const
  {
    const real_type d = l2_norm();
    return d * d;
  }



  PetscScalar
  VectorBase::mean_value() const
  {
    // We can only use our more efficient
    // routine in the serial case.
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(this) != nullptr)
      {
        PetscScalar          sum;
        const PetscErrorCode ierr = VecSum(vector, &sum);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        return sum / static_cast<PetscReal>(size());
      }

    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *  start_ptr;
    PetscErrorCode ierr = VecGetArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    PetscScalar mean = 0;
    {
      PetscScalar sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

      // use modern processors better by
      // allowing pipelined commands to be
      // executed in parallel
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (size() / 4) * 4;
      while (ptr != eptr)
        {
          sum0 += *ptr++;
          sum1 += *ptr++;
          sum2 += *ptr++;
          sum3 += *ptr++;
        }
      // add up remaining elements
      while (ptr != start_ptr + size())
        sum0 += *ptr++;

      mean = (sum0 + sum1 + sum2 + sum3) / static_cast<PetscReal>(size());
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return mean;
  }


  VectorBase::real_type
  VectorBase::l1_norm() const
  {
    real_type d;

    const PetscErrorCode ierr = VecNorm(vector, NORM_1, &d);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::l2_norm() const
  {
    real_type d;

    const PetscErrorCode ierr = VecNorm(vector, NORM_2, &d);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::lp_norm(const real_type p) const
  {
    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *  start_ptr;
    PetscErrorCode ierr = VecGetArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    real_type norm = 0;
    {
      real_type sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

      // use modern processors better by
      // allowing pipelined commands to be
      // executed in parallel
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (size() / 4) * 4;
      while (ptr != eptr)
        {
          sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum1 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum2 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum3 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
        }
      // add up remaining elements
      while (ptr != start_ptr + size())
        sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);

      norm = std::pow(sum0 + sum1 + sum2 + sum3, 1. / p);
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return norm;
  }



  VectorBase::real_type
  VectorBase::linfty_norm() const
  {
    real_type d;

    const PetscErrorCode ierr = VecNorm(vector, NORM_INFINITY, &d);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return d;
  }



  VectorBase::real_type
  VectorBase::min() const
  {
    PetscInt  p;
    real_type d;

    const PetscErrorCode ierr = VecMin(vector, &p, &d);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return d;
  }


  VectorBase::real_type
  VectorBase::max() const
  {
    PetscInt  p;
    real_type d;

    const PetscErrorCode ierr = VecMax(vector, &p, &d);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return d;
  }



  bool
  VectorBase::all_zero() const
  {
    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *  start_ptr;
    PetscErrorCode ierr = VecGetArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const PetscScalar *ptr = start_ptr, *eptr = start_ptr + local_size();
    bool               flag = true;
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
    ierr = VecRestoreArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return flag;
  }


  namespace internal
  {
    template <typename T>
    bool
    is_non_negative(const T &t)
    {
      return t >= 0;
    }



    template <typename T>
    bool
    is_non_negative(const std::complex<T> &)
    {
      Assert(false,
             ExcMessage("You can't ask a complex value "
                        "whether it is non-negative.")) return true;
    }
  } // namespace internal



  bool
  VectorBase::is_non_negative() const
  {
    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *  start_ptr;
    PetscErrorCode ierr = VecGetArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const PetscScalar *ptr = start_ptr, *eptr = start_ptr + local_size();
    bool               flag = true;
    while (ptr != eptr)
      {
        if (!internal::is_non_negative(*ptr))
          {
            flag = false;
            break;
          }
        ++ptr;
      }

    // restore the representation of the
    // vector
    ierr = VecRestoreArray(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return flag;
  }



  VectorBase &
  VectorBase::operator*=(const PetscScalar a)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);

    const PetscErrorCode ierr = VecScale(vector, a);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator/=(const PetscScalar a)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);

    const PetscScalar factor = 1. / a;
    AssertIsFinite(factor);

    const PetscErrorCode ierr = VecScale(vector, factor);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator+=(const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    const PetscErrorCode ierr = VecAXPY(vector, 1, v);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  VectorBase &
  VectorBase::operator-=(const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    const PetscErrorCode ierr = VecAXPY(vector, -1, v);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
  }



  void
  VectorBase::add(const PetscScalar s)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(s);

    const PetscErrorCode ierr = VecShift(vector, s);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::add(const PetscScalar a, const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);

    const PetscErrorCode ierr = VecAXPY(vector, a, v);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::add(const PetscScalar a,
                  const VectorBase &v,
                  const PetscScalar b,
                  const VectorBase &w)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);
    AssertIsFinite(b);

    const PetscScalar weights[2] = {a, b};
    Vec               addends[2] = {v.vector, w.vector};

    const PetscErrorCode ierr = VecMAXPY(vector, 2, weights, addends);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::sadd(const PetscScalar s, const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(s);

    const PetscErrorCode ierr = VecAYPX(vector, s, v);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::sadd(const PetscScalar s,
                   const PetscScalar a,
                   const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(s);
    AssertIsFinite(a);

    // there is nothing like a AXPAY
    // operation in Petsc, so do it in two
    // steps
    *this *= s;
    add(a, v);
  }



  void
  VectorBase::scale(const VectorBase &factors)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    const PetscErrorCode ierr = VecPointwiseMult(vector, factors, vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::equ(const PetscScalar a, const VectorBase &v)
  {
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    AssertIsFinite(a);

    Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

    // there is no simple operation for this
    // in PETSc. there are multiple ways to
    // emulate it, we choose this one:
    const PetscErrorCode ierr = VecCopy(v.vector, vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    *this *= a;
  }



  void
  VectorBase::write_ascii(const PetscViewerFormat format)
  {
    // TODO[TH]:assert(is_compressed())

    // Set options
    PetscErrorCode ierr =
      PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, format);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // Write to screen
    ierr = VecView(vector, PETSC_VIEWER_STDOUT_WORLD);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::print(std::ostream &     out,
                    const unsigned int precision,
                    const bool         scientific,
                    const bool         across) const
  {
    AssertThrow(out, ExcIO());

    // get a representation of the vector and
    // loop over all the elements
    PetscScalar *  val;
    PetscErrorCode ierr = VecGetArray(vector, &val);

    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // save the state of out stream
    const std::ios::fmtflags old_flags     = out.flags();
    const unsigned int       old_precision = out.precision(precision);

    out.precision(precision);
    if (scientific)
      out.setf(std::ios::scientific, std::ios::floatfield);
    else
      out.setf(std::ios::fixed, std::ios::floatfield);

    if (across)
      for (size_type i = 0; i < local_size(); ++i)
        out << val[i] << ' ';
    else
      for (size_type i = 0; i < local_size(); ++i)
        out << val[i] << std::endl;
    out << std::endl;

    // reset output format
    out.flags(old_flags);
    out.precision(old_precision);

    // restore the representation of the
    // vector
    ierr = VecRestoreArray(vector, &val);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    AssertThrow(out, ExcIO());
  }



  void
  VectorBase::swap(VectorBase &v)
  {
    const PetscErrorCode ierr = VecSwap(vector, v.vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  VectorBase::operator const Vec &() const
  {
    return vector;
  }


  std::size_t
  VectorBase::memory_consumption() const
  {
    std::size_t mem = sizeof(Vec) + sizeof(last_action) +
                      MemoryConsumption::memory_consumption(ghosted) +
                      MemoryConsumption::memory_consumption(ghost_indices);

    // TH: I am relatively sure that PETSc is
    // storing the local data in a contiguous
    // block without indices:
    mem += local_size() * sizeof(PetscScalar);
    // assume that PETSc is storing one index
    // and one double per ghost element
    if (ghosted)
      mem += ghost_indices.n_elements() * (sizeof(PetscScalar) + sizeof(int));

    // TODO[TH]: size of constant memory for PETSc?
    return mem;
  }



  void
  VectorBase::do_set_add_operation(const size_type    n_elements,
                                   const size_type *  indices,
                                   const PetscScalar *values,
                                   const bool         add_values)
  {
    ::dealii::VectorOperation::values action =
      (add_values ? ::dealii::VectorOperation::add :
                    ::dealii::VectorOperation::insert);
    Assert((last_action == action) ||
             (last_action == ::dealii::VectorOperation::unknown),
           internal::VectorReference::ExcWrongMode(action, last_action));
    Assert(!has_ghost_elements(), ExcGhostsPresent());
    // VecSetValues complains if we
    // come with an empty
    // vector. however, it is not a
    // collective operation, so we
    // can skip the call if necessary
    // (unlike the above calls)
    if (n_elements != 0)
      {
        const PetscInt *petsc_indices =
          reinterpret_cast<const PetscInt *>(indices);

        const InsertMode     mode = (add_values ? ADD_VALUES : INSERT_VALUES);
        const PetscErrorCode ierr =
          VecSetValues(vector, n_elements, petsc_indices, values, mode);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }

    // set the mode here, independent of whether we have actually
    // written elements or whether the list was empty
    last_action = action;
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
