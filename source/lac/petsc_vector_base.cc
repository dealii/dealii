// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/petsc_vector_base.h>

#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/base/memory_consumption.h>

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

          Vec locally_stored_elements = nullptr;
          ierr = VecGhostGetLocalForm(vector.vector, &locally_stored_elements);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          PetscInt lsize;
          ierr = VecGetSize(locally_stored_elements, &lsize);
          AssertThrow(ierr == 0, ExcPETScError(ierr));

          const PetscScalar *ptr;
          ierr = VecGetArrayRead(locally_stored_elements, &ptr);
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
              Assert(vector.ghost_indices.is_element(index),
                     ExcMessage(
                       "You are trying to access an element of a vector "
                       "that is neither a locally owned element nor a "
                       "ghost element of the vector."));
              const size_type ghostidx =
                vector.ghost_indices.index_within_set(index);

              AssertIndexRange(ghostidx + end - begin, lsize);
              value = *(ptr + ghostidx + end - begin);
            }

          ierr = VecRestoreArrayRead(locally_stored_elements, &ptr);
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

      const PetscScalar *ptr;
      PetscScalar        value;
      ierr = VecGetArrayRead(vector.vector, &ptr);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      value = *(ptr + index - begin);
      ierr  = VecRestoreArrayRead(vector.vector, &ptr);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return value;
    }
#  endif
  } // namespace internal

  VectorBase::VectorBase()
    : vector(nullptr)
    , ghosted(false)
    , last_action(VectorOperation::unknown)
  {}



  VectorBase::VectorBase(const VectorBase &v)
    : ghosted(v.ghosted)
    , ghost_indices(v.ghost_indices)
    , last_action(VectorOperation::unknown)
  {
    PetscErrorCode ierr = VecDuplicate(v.vector, &vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = VecCopy(v.vector, vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  VectorBase::VectorBase(const Vec &v)
    : vector(v)
    , ghosted(false)
    , last_action(VectorOperation::unknown)
  {
    const PetscErrorCode ierr =
      PetscObjectReference(reinterpret_cast<PetscObject>(vector));
    AssertNothrow(ierr == 0, ExcPETScError(ierr));
    this->determine_ghost_indices();
  }



  VectorBase::~VectorBase()
  {
    const PetscErrorCode ierr = VecDestroy(&vector);
    AssertNothrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::reinit(Vec v)
  {
    AssertThrow(last_action == VectorOperation::unknown,
                ExcMessage("Cannot assign a new Vec"));
    PetscErrorCode ierr =
      PetscObjectReference(reinterpret_cast<PetscObject>(v));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = VecDestroy(&vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    vector = v;
    this->determine_ghost_indices();
  }



  namespace
  {
    template <typename Iterator, typename OutType>
    class ConvertingIterator
    {
      Iterator m_iterator;

    public:
      using difference_type =
        typename std::iterator_traits<Iterator>::difference_type;
      using value_type        = OutType;
      using pointer           = OutType *;
      using reference         = OutType &;
      using iterator_category = std::forward_iterator_tag;

      ConvertingIterator(const Iterator &iterator)
        : m_iterator(iterator)
      {}

      OutType
      operator*() const
      {
        return static_cast<OutType>(std::real(*m_iterator));
      }

      ConvertingIterator &
      operator++()
      {
        ++m_iterator;
        return *this;
      }

      ConvertingIterator
      operator++(int)
      {
        ConvertingIterator old = *this;
        ++m_iterator;
        return old;
      }

      bool
      operator==(const ConvertingIterator &other) const
      {
        return this->m_iterator == other.m_iterator;
      }

      bool
      operator!=(const ConvertingIterator &other) const
      {
        return this->m_iterator != other.m_iterator;
      }
    };
  } // namespace



  void
  VectorBase::determine_ghost_indices()
  {
    // Reset ghost data
    ghosted = false;
    ghost_indices.clear();

    // There's no API to infer ghost indices from a PETSc Vec which
    // unfortunately doesn't allow integer entries. We use the
    // "ConvertingIterator" class above to do an implicit conversion when
    // sorting and adding ghost indices below.
    PetscErrorCode ierr;
    Vec            ghosted_vec;
    ierr = VecGhostGetLocalForm(vector, &ghosted_vec);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    if (ghosted_vec && ghosted_vec != vector)
      {
        Vec          tvector;
        PetscScalar *array;
        PetscInt     ghost_start_index, end_index, n_elements_stored_locally;

        ierr = VecGhostRestoreLocalForm(vector, &ghosted_vec);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecGetOwnershipRange(vector, &ghost_start_index, &end_index);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecDuplicate(vector, &tvector);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGetArray(tvector, &array);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        // Store the indices we care about in the vector, so that we can then
        // exchange this information between processes. It is unfortunate that
        // we have to store integers in floating point numbers. Let's at least
        // make sure we do that in a way that ensures that when we get these
        // numbers back as integers later on, we get the same thing.
        for (PetscInt i = 0; i < end_index - ghost_start_index; i++)
          {
            Assert(static_cast<PetscInt>(std::real(static_cast<PetscScalar>(
                     ghost_start_index + i))) == (ghost_start_index + i),
                   ExcInternalError());
            array[i] = ghost_start_index + i;
          }

        ierr = VecRestoreArray(tvector, &array);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGhostUpdateBegin(tvector, INSERT_VALUES, SCATTER_FORWARD);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGhostUpdateEnd(tvector, INSERT_VALUES, SCATTER_FORWARD);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGhostGetLocalForm(tvector, &ghosted_vec);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGetLocalSize(ghosted_vec, &n_elements_stored_locally);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGetArrayRead(ghosted_vec, (const PetscScalar **)&array);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        // Populate the 'ghosted' flag and the ghost_indices variable. The
        // latter is an index set that is most efficiently filled by
        // sorting the indices to add. At the same time, we don't want to
        // sort the indices stored in a PETSc-owned array; so if the array
        // is already sorted, pass that to the IndexSet variable, and if
        // not then copy the indices, sort them, and then add those.
        ghosted = true;
        ghost_indices.set_size(this->size());

        ConvertingIterator<PetscScalar *, types::global_dof_index> begin_ghosts(
          &array[end_index - ghost_start_index]);
        ConvertingIterator<PetscScalar *, types::global_dof_index> end_ghosts(
          &array[n_elements_stored_locally]);
        if (std::is_sorted(&array[end_index - ghost_start_index],
                           &array[n_elements_stored_locally],
                           [](PetscScalar left, PetscScalar right) {
                             return static_cast<PetscInt>(std::real(left)) <
                                    static_cast<PetscInt>(std::real(right));
                           }))
          {
            ghost_indices.add_indices(begin_ghosts, end_ghosts);
          }
        else
          {
            std::vector<PetscInt> sorted_indices(begin_ghosts, end_ghosts);
            std::sort(sorted_indices.begin(), sorted_indices.end());
            ghost_indices.add_indices(sorted_indices.begin(),
                                      sorted_indices.end());
          }
        ghost_indices.compress();

        ierr = VecRestoreArrayRead(ghosted_vec, (const PetscScalar **)&array);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecGhostRestoreLocalForm(tvector, &ghosted_vec);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
        ierr = VecDestroy(&tvector);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    else
      {
        ierr = VecGhostRestoreLocalForm(vector, &ghosted_vec);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
  }


  void
  VectorBase::clear()
  {
    const PetscErrorCode ierr = VecDestroy(&vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ghosted = false;
    ghost_indices.clear();
    last_action = VectorOperation::unknown;
  }



  VectorBase &
  VectorBase::operator=(const VectorBase &v)
  {
    Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

    PetscErrorCode ierr = VecCopy(v, vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return *this;
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
        Vec ghost = nullptr;
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
  VectorBase::locally_owned_size() const
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
  VectorBase::set(const std::vector<size_type>   &indices,
                  const std::vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.data(), false);
  }



  void
  VectorBase::add(const std::vector<size_type>   &indices,
                  const std::vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.data(), true);
  }



  void
  VectorBase::add(const std::vector<size_type>        &indices,
                  const ::dealii::Vector<PetscScalar> &values)
  {
    Assert(indices.size() == values.size(),
           ExcMessage("Function called with arguments of different sizes"));
    do_set_add_operation(indices.size(), indices.data(), values.begin(), true);
  }



  void
  VectorBase::add(const size_type    n_elements,
                  const size_type   *indices,
                  const PetscScalar *values)
  {
    do_set_add_operation(n_elements, indices, values, true);
  }



  PetscScalar
  VectorBase::operator*(const VectorBase &vec) const
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
    Assert(has_ghost_elements() == false,
           ExcMessage("Calling compress() is only useful if a vector "
                      "has been written into, but this is a vector with ghost "
                      "elements and consequently is read-only. It does "
                      "not make sense to call compress() for such "
                      "vectors."));

    {
      if constexpr (running_in_debug_mode())
        {
          // Check that all processors agree that last_action is the same (or
          // none!)

          int my_int_last_action = last_action;
          int all_int_last_action;

          const int ierr = MPI_Allreduce(&my_int_last_action,
                                         &all_int_last_action,
                                         1,
                                         MPI_INT,
                                         MPI_BOR,
                                         get_mpi_communicator());
          AssertThrowMPI(ierr);

          AssertThrow(all_int_last_action !=
                        (VectorOperation::add | VectorOperation::insert),
                      ExcMessage(
                        "Error: not all processors agree on the last "
                        "VectorOperation before this compress() call."));
        }
    }

    AssertThrow(
      last_action == VectorOperation::unknown || last_action == operation,
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
    last_action = VectorOperation::unknown;
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
    const PetscScalar *start_ptr;
    PetscErrorCode     ierr = VecGetArrayRead(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    PetscScalar mean = 0;
    {
      PetscScalar sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

      // use modern processors better by
      // allowing pipelined commands to be
      // executed in parallel
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (locally_owned_size() / 4) * 4;
      while (ptr != eptr)
        {
          sum0 += *ptr++;
          sum1 += *ptr++;
          sum2 += *ptr++;
          sum3 += *ptr++;
        }
      // add up remaining elements
      while (ptr != start_ptr + locally_owned_size())
        sum0 += *ptr++;

      mean =
        Utilities::MPI::sum(sum0 + sum1 + sum2 + sum3, get_mpi_communicator()) /
        static_cast<PetscReal>(size());
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArrayRead(vector, &start_ptr);
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
    const PetscScalar *start_ptr;
    PetscErrorCode     ierr = VecGetArrayRead(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    real_type norm = 0;
    {
      real_type sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

      // use modern processors better by
      // allowing pipelined commands to be
      // executed in parallel
      const PetscScalar *ptr  = start_ptr;
      const PetscScalar *eptr = ptr + (locally_owned_size() / 4) * 4;
      while (ptr != eptr)
        {
          sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum1 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum2 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
          sum3 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);
        }
      // add up remaining elements
      while (ptr != start_ptr + locally_owned_size())
        sum0 += std::pow(numbers::NumberTraits<value_type>::abs(*ptr++), p);

      norm = std::pow(Utilities::MPI::sum(sum0 + sum1 + sum2 + sum3,
                                          get_mpi_communicator()),
                      1. / p);
    }

    // restore the representation of the
    // vector
    ierr = VecRestoreArrayRead(vector, &start_ptr);
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



  bool
  VectorBase::all_zero() const
  {
    // get a representation of the vector and
    // loop over all the elements
    const PetscScalar *start_ptr;
    PetscErrorCode     ierr = VecGetArrayRead(vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const PetscScalar *ptr  = start_ptr,
                      *eptr = start_ptr + locally_owned_size();
    bool flag               = true;
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
    ierr = VecRestoreArrayRead(vector, &start_ptr);
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
                        "whether it is non-negative."));
      return true;
    }
  } // namespace internal



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
    // operation in PETSc, so do it in two
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

    const PetscErrorCode ierr = VecAXPBY(vector, a, 0.0, v.vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::write_ascii(const PetscViewerFormat format)
  {
    // TODO[TH]:assert(is_compressed())
    MPI_Comm comm = PetscObjectComm((PetscObject)vector);

    // Set options
    PetscErrorCode ierr =
      PetscViewerPushFormat(PETSC_VIEWER_STDOUT_(comm), format);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    // Write to screen
    ierr = VecView(vector, PETSC_VIEWER_STDOUT_(comm));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_(comm));
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  void
  VectorBase::print(std::ostream      &out,
                    const unsigned int precision,
                    const bool         scientific,
                    const bool         across) const
  {
    AssertThrow(out.fail() == false, ExcIO());

    // get a representation of the vector and
    // loop over all the elements
    const PetscScalar *val;
    PetscErrorCode     ierr = VecGetArrayRead(vector, &val);

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
      for (size_type i = 0; i < locally_owned_size(); ++i)
        out << val[i] << ' ';
    else
      for (size_type i = 0; i < locally_owned_size(); ++i)
        out << val[i] << std::endl;
    out << std::endl;

    // reset output format
    out.flags(old_flags);
    out.precision(old_precision);

    // restore the representation of the
    // vector
    ierr = VecRestoreArrayRead(vector, &val);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    AssertThrow(out.fail() == false, ExcIO());
  }



  void
  VectorBase::swap(VectorBase &v) noexcept
  {
    std::swap(this->vector, v.vector);
    std::swap(this->ghosted, v.ghosted);
    std::swap(this->last_action, v.last_action);
    // missing swap for IndexSet
    IndexSet t(this->ghost_indices);
    this->ghost_indices = v.ghost_indices;
    v.ghost_indices     = t;
  }


  VectorBase::operator const Vec &() const
  {
    return vector;
  }


  Vec &
  VectorBase::petsc_vector()
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
    mem += locally_owned_size() * sizeof(PetscScalar);
    // assume that PETSc is storing one index
    // and one double per ghost element
    if (ghosted)
      mem += ghost_indices.n_elements() * (sizeof(PetscScalar) + sizeof(int));

    // TODO[TH]: size of constant memory for PETSc?
    return mem;
  }



  void
  VectorBase::do_set_add_operation(const size_type    n_elements,
                                   const size_type   *indices,
                                   const PetscScalar *values,
                                   const bool         add_values)
  {
    VectorOperation::values action =
      (add_values ? VectorOperation::add : VectorOperation::insert);
    Assert((last_action == action) || (last_action == VectorOperation::unknown),
           internal::VectorReference::ExcWrongMode(action, last_action));
    Assert(!has_ghost_elements(), ExcGhostsPresent());

    std::vector<PetscInt> petsc_indices(n_elements);
    for (size_type i = 0; i < n_elements; ++i)
      {
        const auto petsc_index = static_cast<PetscInt>(indices[i]);
        AssertIntegerConversion(petsc_index, indices[i]);
        petsc_indices[i] = petsc_index;
      }

    const InsertMode     mode = (add_values ? ADD_VALUES : INSERT_VALUES);
    const PetscErrorCode ierr = VecSetValues(
      vector, petsc_indices.size(), petsc_indices.data(), values, mode);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    last_action = action;
  }

} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
