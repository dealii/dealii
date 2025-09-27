// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/trilinos_epetra_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/trilinos_utilities.h>

#  include <deal.II/lac/read_write_vector.h>

#  include <boost/io/ios_state.hpp>

#  include <Epetra_Import.h>
#  include <Epetra_Map.h>
#  include <Epetra_MpiComm.h>

#  include <memory>


DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
#  ifndef DOXYGEN
    namespace internal
    {
      VectorReference::operator value_type() const
      {
        AssertIndexRange(index, vector.size());

        // Trilinos allows for vectors to be referenced by the [] or ()
        // operators but only () checks index bounds. We check these bounds by
        // ourselves, so we can use []. Note that we can only get local values.

        const TrilinosWrappers::types::int_type local_index =
          vector.vector->Map().LID(
            static_cast<TrilinosWrappers::types::int_type>(index));

#    ifndef DEAL_II_WITH_64BIT_INDICES
        Assert(local_index >= 0,
               ExcAccessToNonLocalElement(index,
                                          vector.vector->Map().NumMyElements(),
                                          vector.vector->Map().MinMyGID(),
                                          vector.vector->Map().MaxMyGID()));
#    else
        Assert(local_index >= 0,
               ExcAccessToNonLocalElement(index,
                                          vector.vector->Map().NumMyElements(),
                                          vector.vector->Map().MinMyGID64(),
                                          vector.vector->Map().MaxMyGID64()));
#    endif

        return (*(vector.vector))[0][local_index];
      }
    } // namespace internal
#  endif


    // Check that the class we declare here satisfies the
    // vector-space-vector concept. If we catch it here,
    // any mistake in the vector class declaration would
    // show up in uses of this class later on as well.
#  ifdef DEAL_II_HAVE_CXX20
    static_assert(concepts::is_vector_space_vector<Vector>);
#  endif

    Vector::Vector()
      : vector(new Epetra_FEVector(
          Epetra_Map(0, 0, 0, Utilities::Trilinos::comm_self())))
    {}



    Vector::Vector(const Vector &V)
      : vector(new Epetra_FEVector(V.trilinos_vector()))
    {}



    Vector::Vector(const IndexSet &parallel_partitioner,
                   const MPI_Comm  communicator)
      : vector(new Epetra_FEVector(
          parallel_partitioner.make_trilinos_map(communicator, false)))
    {}



    void
    Vector::reinit(const IndexSet &parallel_partitioner,
                   const MPI_Comm  communicator,
                   const bool      omit_zeroing_entries)
    {
      Epetra_Map input_map =
        parallel_partitioner.make_trilinos_map(communicator, false);
      if (vector->Map().SameAs(input_map) == false)
        vector = std::make_unique<Epetra_FEVector>(input_map);
      else if (omit_zeroing_entries == false)
        {
          const int ierr = vector->PutScalar(0.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
        }
    }



    void
    Vector::reinit(const Vector &V, const bool omit_zeroing_entries)
    {
      reinit(V.locally_owned_elements(),
             V.get_mpi_communicator(),
             omit_zeroing_entries);
    }



    void
    Vector::extract_subvector_to(
      const ArrayView<const types::global_dof_index> &indices,
      const ArrayView<double>                        &elements) const
    {
      AssertDimension(indices.size(), elements.size());
      const auto &vector = trilinos_vector();
      const auto &map    = vector.Map();

      for (unsigned int i = 0; i < indices.size(); ++i)
        {
          AssertIndexRange(indices[i], size());
          const auto trilinos_i =
            map.LID(static_cast<TrilinosWrappers::types::int_type>(indices[i]));
          elements[i] = vector[0][trilinos_i];
        }
    }



    Vector &
    Vector::operator=(const Vector &V)
    {
      // Distinguish three cases:
      //  - First case: both vectors have the same layout.
      //  - Second case: both vectors have the same size but different layout.
      //  - Third case: the vectors have different size.
      if (vector->Map().SameAs(V.trilinos_vector().Map()))
        *vector = V.trilinos_vector();
      else
        {
          if (size() == V.size())
            {
              Epetra_Import data_exchange(vector->Map(),
                                          V.trilinos_vector().Map());

              const int ierr =
                vector->Import(V.trilinos_vector(), data_exchange, Insert);
              Assert(ierr == 0, ExcTrilinosError(ierr));
            }
          else
            vector = std::make_unique<Epetra_FEVector>(V.trilinos_vector());
        }

      return *this;
    }



    Vector &
    Vector::operator=(const double s)
    {
      Assert(s == 0., ExcMessage("Only 0 can be assigned to a vector."));

      const int ierr = vector->PutScalar(s);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return *this;
    }



    void
    Vector::import_elements(
      const ReadWriteVector<double> &V,
      VectorOperation::values        operation,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern)
    {
      // If no communication pattern is given, create one. Otherwise, use the
      // one given.
      if (communication_pattern == nullptr)
        {
          // The first time import is called, a communication pattern is
          // created. Check if the communication pattern already exists and if
          // it can be reused.
          if ((source_stored_elements.size() !=
               V.get_stored_elements().size()) ||
              (source_stored_elements != V.get_stored_elements()))
            {
              create_epetra_comm_pattern(
                V.get_stored_elements(),
                dynamic_cast<const Epetra_MpiComm &>(vector->Comm()).Comm());
            }
        }
      else
        {
          epetra_comm_pattern =
            std::dynamic_pointer_cast<const CommunicationPattern>(
              communication_pattern);
          AssertThrow(
            epetra_comm_pattern != nullptr,
            ExcMessage("The communication pattern is not of type "
                       "LinearAlgebra::EpetraWrappers::CommunicationPattern."));
        }

      Epetra_Import import_map(epetra_comm_pattern->get_epetra_import());

      // The TargetMap and the SourceMap have their roles inverted.
      Epetra_FEVector source_vector(import_map.TargetMap());
      double         *values = source_vector.Values();
      std::copy(V.begin(), V.end(), values);

      if (operation == VectorOperation::insert)
        vector->Export(source_vector, import_map, Insert);
      else if (operation == VectorOperation::add)
        vector->Export(source_vector, import_map, Add);
      else if (operation == VectorOperation::max)
        vector->Export(source_vector, import_map, Epetra_Max);
      else if (operation == VectorOperation::min)
        vector->Export(source_vector, import_map, Epetra_Min);
      else
        AssertThrow(false, ExcNotImplemented());
    }



    Vector &
    Vector::operator*=(const double factor)
    {
      AssertIsFinite(factor);
      vector->Scale(factor);

      return *this;
    }



    Vector &
    Vector::operator/=(const double factor)
    {
      AssertIsFinite(factor);
      Assert(factor != 0., ExcZero());
      *this *= 1. / factor;

      return *this;
    }



    Vector &
    Vector::operator+=(const Vector &V)
    {
      // If the maps are the same we can Update right away.
      if (vector->Map().SameAs(V.trilinos_vector().Map()))
        {
          const int ierr = vector->Update(1., V.trilinos_vector(), 1.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
        }
      else
        {
          Assert(this->size() == V.size(),
                 ExcDimensionMismatch(this->size(), V.size()));

          Epetra_Import data_exchange(vector->Map(), V.trilinos_vector().Map());
          const int     ierr = vector->Import(V.trilinos_vector(),
                                          data_exchange,
                                          Epetra_AddLocalAlso);
          Assert(ierr == 0, ExcTrilinosError(ierr));
        }

      return *this;
    }



    Vector &
    Vector::operator-=(const Vector &V)
    {
      this->add(-1., V);

      return *this;
    }



    double
    Vector::operator*(const Vector &V) const
    {
      Assert(this->size() == V.size(),
             ExcDimensionMismatch(this->size(), V.size()));
      Assert(vector->Map().SameAs(V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      double    result(0.);
      const int ierr = vector->Dot(V.trilinos_vector(), &result);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return result;
    }



    void
    Vector::add(const double a)
    {
      AssertIsFinite(a);
      const unsigned local_size(vector->MyLength());
      for (unsigned int i = 0; i < local_size; ++i)
        (*vector)[0][i] += a;
    }



    void
    Vector::add(const double a, const Vector &V)
    {
      AssertIsFinite(a);
      Assert(vector->Map().SameAs(V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr = vector->Update(a, V.trilinos_vector(), 1.);
      Assert(ierr == 0, ExcTrilinosError(ierr));
    }



    void
    Vector::add(const double  a,
                const Vector &V,
                const double  b,
                const Vector &W)
    {
      Assert(vector->Map().SameAs(V.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());
      Assert(vector->Map().SameAs(W.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());
      AssertIsFinite(a);
      AssertIsFinite(b);

      const int ierr =
        vector->Update(a, V.trilinos_vector(), b, W.trilinos_vector(), 1.);
      Assert(ierr == 0, ExcTrilinosError(ierr));
    }



    void
    Vector::sadd(const double s, const double a, const Vector &V)
    {
      *this *= s;
      Vector tmp(V);
      tmp *= a;
      *this += tmp;
    }



    void
    Vector::scale(const Vector &scaling_factors)
    {
      Assert(vector->Map().SameAs(scaling_factors.trilinos_vector().Map()),
             ExcDifferentParallelPartitioning());

      const int ierr =
        vector->Multiply(1.0, scaling_factors.trilinos_vector(), *vector, 0.0);
      Assert(ierr == 0, ExcTrilinosError(ierr));
    }



    void
    Vector::equ(const double a, const Vector &V)
    {
      // If we don't have the same map, copy.
      if (vector->Map().SameAs(V.trilinos_vector().Map()) == false)
        this->sadd(0., a, V);
      else
        {
          // Otherwise, just update
          int ierr = vector->Update(a, V.trilinos_vector(), 0.);
          Assert(ierr == 0, ExcTrilinosError(ierr));
        }
    }



    bool
    Vector::all_zero() const
    {
      // get a representation of the vector and
      // loop over all the elements
      double       *start_ptr = (*vector)[0];
      const double *ptr = start_ptr, *eptr = start_ptr + vector->MyLength();
      unsigned int  flag = 0;
      while (ptr != eptr)
        {
          if (*ptr != 0)
            {
              flag = 1;
              break;
            }
          ++ptr;
        }

      // Check that the vector is zero on _all_ processors.
      const Epetra_MpiComm *mpi_comm =
        dynamic_cast<const Epetra_MpiComm *>(&vector->Map().Comm());
      Assert(mpi_comm != nullptr, ExcInternalError());
      unsigned int num_nonzero = Utilities::MPI::sum(flag, mpi_comm->Comm());

      return num_nonzero == 0;
    }



    double
    Vector::mean_value() const
    {
      double mean_value(0.);

      int ierr = vector->MeanValue(&mean_value);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return mean_value;
    }



    double
    Vector::l1_norm() const
    {
      double norm(0.);
      int    ierr = vector->Norm1(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return norm;
    }



    double
    Vector::l2_norm() const
    {
      double norm(0.);
      int    ierr = vector->Norm2(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return norm;
    }



    double
    Vector::linfty_norm() const
    {
      double norm(0.);
      int    ierr = vector->NormInf(&norm);
      Assert(ierr == 0, ExcTrilinosError(ierr));

      return norm;
    }



    double
    Vector::add_and_dot(const double a, const Vector &V, const Vector &W)
    {
      this->add(a, V);

      return *this * W;
    }



    Vector::size_type
    Vector::size() const
    {
#  ifndef DEAL_II_WITH_64BIT_INDICES
      return vector->GlobalLength();
#  else
      return vector->GlobalLength64();
#  endif
    }



    Vector::size_type
    Vector::locally_owned_size() const
    {
      return vector->MyLength();
    }



    MPI_Comm
    Vector::get_mpi_communicator() const
    {
      const Epetra_MpiComm *epetra_comm =
        dynamic_cast<const Epetra_MpiComm *>(&(vector->Comm()));
      Assert(epetra_comm != nullptr, ExcInternalError());
      return epetra_comm->GetMpiComm();
    }



    ::dealii::IndexSet
    Vector::locally_owned_elements() const
    {
      IndexSet is(size());

      // easy case: local range is contiguous
      if (vector->Map().LinearMap())
        {
#  ifndef DEAL_II_WITH_64BIT_INDICES
          is.add_range(vector->Map().MinMyGID(), vector->Map().MaxMyGID() + 1);
#  else
          is.add_range(vector->Map().MinMyGID64(),
                       vector->Map().MaxMyGID64() + 1);
#  endif
        }
      else if (vector->Map().NumMyElements() > 0)
        {
          const size_type n_indices = vector->Map().NumMyElements();
#  ifndef DEAL_II_WITH_64BIT_INDICES
          unsigned int *vector_indices =
            reinterpret_cast<unsigned int *>(vector->Map().MyGlobalElements());
#  else
          size_type *vector_indices =
            reinterpret_cast<size_type *>(vector->Map().MyGlobalElements64());
#  endif
          is.add_indices(vector_indices, vector_indices + n_indices);
        }
      is.compress();

      return is;
    }


    void
    Vector::compress(const VectorOperation::values /*operation*/)
    {}



    const Epetra_FEVector &
    Vector::trilinos_vector() const
    {
      return *vector;
    }



    Epetra_FEVector &
    Vector::trilinos_vector()
    {
      return *vector;
    }



    void
    Vector::print(std::ostream      &out,
                  const unsigned int precision,
                  const bool         scientific,
                  const bool         across) const
    {
      AssertThrow(out.fail() == false, ExcIO());
      boost::io::ios_flags_saver restore_flags(out);

      // Get a representation of the vector and loop over all
      // the elements
      double *val;
      int     leading_dimension;
      int     ierr = vector->ExtractView(&val, &leading_dimension);

      Assert(ierr == 0, ExcTrilinosError(ierr));
      out.precision(precision);
      if (scientific)
        out.setf(std::ios::scientific, std::ios::floatfield);
      else
        out.setf(std::ios::fixed, std::ios::floatfield);

      if (across)
        for (int i = 0; i < vector->MyLength(); ++i)
          out << val[i] << ' ';
      else
        for (int i = 0; i < vector->MyLength(); ++i)
          out << val[i] << std::endl;
      out << std::endl;

      // restore the representation
      // of the vector
      AssertThrow(out.fail() == false, ExcIO());
    }



    std::size_t
    Vector::memory_consumption() const
    {
      return sizeof(*this) +
             vector->MyLength() *
               (sizeof(double) + sizeof(TrilinosWrappers::types::int_type));
    }



    void
    Vector::create_epetra_comm_pattern(const IndexSet &source_index_set,
                                       const MPI_Comm  mpi_comm)
    {
      source_stored_elements = source_index_set;
      epetra_comm_pattern =
        std::make_shared<CommunicationPattern>(locally_owned_elements(),
                                               source_index_set,
                                               mpi_comm);
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
