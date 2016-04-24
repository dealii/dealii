// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2015 by the deal.II authors
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

#include <deal.II/lac/trilinos_vector.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_block_vector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <Epetra_Import.h>
#  include <Epetra_Vector.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <cmath>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace
  {
#ifndef DEAL_II_WITH_64BIT_INDICES
    // define a helper function that queries the size of an Epetra_BlockMap object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    int n_global_elements (const Epetra_BlockMap &map)
    {
      return map.NumGlobalElements();
    }
    // define a helper function that queries the pointer to internal array
    // containing list of global IDs assigned to the calling processor
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    int *my_global_elements(const Epetra_BlockMap &map)
    {
      return map.MyGlobalElements();
    }
    // define a helper function that queries the global vector length of an
    // Epetra_FEVector object  by calling either the 32- or 64-bit
    // function necessary.
    int global_length(const Epetra_FEVector &vector)
    {
      return vector.GlobalLength();
    }
#else
    // define a helper function that queries the size of an Epetra_BlockMap object
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    long long int n_global_elements (const Epetra_BlockMap &map)
    {
      return map.NumGlobalElements64();
    }
    // define a helper function that queries the pointer to internal array
    // containing list of global IDs assigned to the calling processor
    // by calling either the 32- or 64-bit function necessary, and returns the
    // result in the correct data type so that we can use it in calling other
    // Epetra member functions that are overloaded by index type
    long long int *my_global_elements(const Epetra_BlockMap &map)
    {
      return map.MyGlobalElements64();
    }
    // define a helper function that queries the global vector length of an
    // Epetra_FEVector object  by calling either the 32- or 64-bit
    // function necessary.
    long long int global_length(const Epetra_FEVector &vector)
    {
      return vector.GlobalLength64();
    }
#endif
  }

  namespace MPI
  {


    Vector::Vector ()
    {
      last_action = Zero;
      vector.reset(new Epetra_FEVector(Epetra_Map(0,0,0,Utilities::Trilinos::comm_self())));
    }



    Vector::Vector (const Epetra_Map &parallel_partitioning)
    {
      reinit (parallel_partitioning);
    }



    Vector::Vector (const IndexSet &parallel_partitioning,
                    const MPI_Comm &communicator)
    {
      reinit (parallel_partitioning, communicator);
    }



    Vector::Vector (const Vector &v)
      :
      VectorBase()
    {
      last_action = Zero;
      vector.reset (new Epetra_FEVector(*v.vector));
      has_ghosts = v.has_ghosts;
    }



#ifdef DEAL_II_WITH_CXX11
    Vector::Vector (Vector &&v)
    {
      // initialize a minimal, valid object and swap
      last_action = Zero;
      vector.reset(new Epetra_FEVector(Epetra_Map(0,0,0,Utilities::Trilinos::comm_self())));

      swap(v);
    }
#endif



    Vector::Vector (const Epetra_Map &input_map,
                    const VectorBase &v)
      :
      VectorBase()
    {
      AssertThrow (n_global_elements(input_map) == n_global_elements(v.vector->Map()),
                   ExcDimensionMismatch (n_global_elements(input_map),
                                         n_global_elements(v.vector->Map())));

      last_action = Zero;

      if (input_map.SameAs(v.vector->Map()) == true)
        vector.reset (new Epetra_FEVector(*v.vector));
      else
        {
          vector.reset (new Epetra_FEVector(input_map));
          reinit (v, false, true);
        }
    }



    Vector::Vector (const IndexSet   &parallel_partitioner,
                    const VectorBase &v,
                    const MPI_Comm   &communicator)
      :
      VectorBase()
    {
      AssertThrow (parallel_partitioner.size() ==
                   static_cast<size_type>(n_global_elements(v.vector->Map())),
                   ExcDimensionMismatch (parallel_partitioner.size(),
                                         n_global_elements(v.vector->Map())));

      last_action = Zero;

      vector.reset (new Epetra_FEVector
                    (parallel_partitioner.make_trilinos_map(communicator,
                                                            true)));
      reinit (v, false, true);
    }

    Vector::Vector (const IndexSet &local,
                    const IndexSet &ghost,
                    const MPI_Comm &communicator)
      :
      VectorBase()
    {
      IndexSet parallel_partitioning = local;
      parallel_partitioning.add_indices(ghost);
      reinit(parallel_partitioning, communicator);
    }



    Vector::~Vector ()
    {}



    void
    Vector::reinit (const Epetra_Map &input_map,
                    const bool      /*omit_zeroing_entries*/)
    {
      nonlocal_vector.reset();

      vector.reset (new Epetra_FEVector(input_map));

      has_ghosts = vector->Map().UniqueGIDs()==false;
      last_action = Zero;
    }



    void
    Vector::reinit (const IndexSet &parallel_partitioner,
                    const MPI_Comm &communicator,
                    const bool      omit_zeroing_entries)
    {
      nonlocal_vector.reset();

      Epetra_Map map = parallel_partitioner.make_trilinos_map (communicator,
                                                               true);
      reinit (map, omit_zeroing_entries);
    }



    void
    Vector::reinit (const VectorBase &v,
                    const bool        omit_zeroing_entries,
                    const bool        allow_different_maps)
    {
      nonlocal_vector.reset();

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (allow_different_maps == false)
        {
          // check equality for MPI communicators: We can only choose the fast
          // version in case the underlying Epetra_MpiComm object is the same,
          // otherwise we might access an MPI_Comm object that has been
          // deleted
#ifdef DEAL_II_WITH_MPI
          const Epetra_MpiComm *my_comm = dynamic_cast<const Epetra_MpiComm *>(&vector->Comm());
          const Epetra_MpiComm *v_comm = dynamic_cast<const Epetra_MpiComm *>(&v.vector->Comm());
          const bool same_communicators = my_comm != NULL && v_comm != NULL &&
                                          my_comm->DataPtr() == v_comm->DataPtr();
#else
          const bool same_communicators = true;
#endif
          if (!same_communicators || vector->Map().SameAs(v.vector->Map()) == false)
            {
              vector.reset (new Epetra_FEVector(v.vector->Map()));
              has_ghosts = v.has_ghosts;
              last_action = Zero;
            }
          else if (omit_zeroing_entries == false)
            {
              // old and new vectors have exactly the same map, i.e. size and
              // parallel distribution
              int ierr;
              ierr = vector->GlobalAssemble (last_action);
              (void)ierr;
              Assert (ierr == 0, ExcTrilinosError(ierr));

              ierr = vector->PutScalar(0.0);
              Assert (ierr == 0, ExcTrilinosError(ierr));

              last_action = Zero;
            }
        }

      // Otherwise, we have to check that the two vectors are already of the
      // same size, create an object for the data exchange and then insert all
      // the data. The first assertion is only a check whether the user knows
      // what she is doing.
      else
        {
          Assert (omit_zeroing_entries == false,
                  ExcMessage ("It is not possible to exchange data with the "
                              "option 'omit_zeroing_entries' set, which would not write "
                              "elements."));

          AssertThrow (size() == v.size(),
                       ExcDimensionMismatch (size(), v.size()));

          Epetra_Import data_exchange (vector->Map(), v.vector->Map());

          const int ierr = vector->Import(*v.vector, data_exchange, Insert);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));

          last_action = Insert;
        }

    }



    void
    Vector::reinit (const BlockVector &v,
                    const bool         import_data)
    {
      nonlocal_vector.reset();

      // In case we do not allow to have different maps, this call means that
      // we have to reset the vector. So clear the vector, initialize our map
      // with the map in v, and generate the vector.
      if (v.n_blocks() == 0)
        return;

      // create a vector that holds all the elements contained in the block
      // vector. need to manually create an Epetra_Map.
      size_type n_elements = 0, added_elements = 0, block_offset = 0;
      for (size_type block=0; block<v.n_blocks(); ++block)
        n_elements += v.block(block).local_size();
      std::vector<TrilinosWrappers::types::int_type> global_ids (n_elements, -1);
      for (size_type block=0; block<v.n_blocks(); ++block)
        {
          TrilinosWrappers::types::int_type *glob_elements =
            my_global_elements(v.block(block).vector_partitioner());
          for (size_type i=0; i<v.block(block).local_size(); ++i)
            global_ids[added_elements++] = glob_elements[i] + block_offset;
          block_offset += v.block(block).size();
        }

      Assert (n_elements == added_elements, ExcInternalError());
      Epetra_Map new_map (v.size(), n_elements, &global_ids[0], 0,
                          v.block(0).vector_partitioner().Comm());

      std_cxx11::shared_ptr<Epetra_FEVector> actual_vec;
      if ( import_data == true )
        actual_vec.reset (new Epetra_FEVector (new_map));
      else
        {
          vector.reset (new Epetra_FEVector (new_map));
          actual_vec = vector;
        }

      TrilinosScalar *entries = (*actual_vec)[0];
      block_offset = 0;
      for (size_type block=0; block<v.n_blocks(); ++block)
        {
          v.block(block).trilinos_vector().ExtractCopy (entries, 0);
          entries += v.block(block).local_size();
        }

      if (import_data == true)
        {
          AssertThrow (static_cast<size_type>(global_length(*actual_vec))
                       == v.size(),
                       ExcDimensionMismatch (global_length(*actual_vec),
                                             v.size()));

          Epetra_Import data_exchange (vector->Map(), actual_vec->Map());

          const int ierr = vector->Import(*actual_vec, data_exchange, Insert);
          AssertThrow (ierr == 0, ExcTrilinosError(ierr));

          last_action = Insert;
        }

    }


    void Vector::reinit(const IndexSet &locally_owned_entries,
                        const IndexSet &ghost_entries,
                        const MPI_Comm &communicator,
                        const bool      vector_writable)
    {
      nonlocal_vector.reset();
      if (vector_writable == false)
        {
          IndexSet parallel_partitioning = locally_owned_entries;
          parallel_partitioning.add_indices(ghost_entries);
          reinit(parallel_partitioning, communicator);
        }
      else
        {
          Epetra_Map map = locally_owned_entries.make_trilinos_map (communicator,
                                                                    true);
          Assert (map.IsOneToOne(),
                  ExcMessage("A writable vector must not have ghost entries in "
                             "its parallel partitioning"));
          reinit (map);

          IndexSet nonlocal_entries(ghost_entries);
          nonlocal_entries.subtract_set(locally_owned_entries);
          if (Utilities::MPI::n_mpi_processes(communicator) > 1)
            {
              Epetra_Map nonlocal_map =
                nonlocal_entries.make_trilinos_map(communicator, true);
              nonlocal_vector.reset(new Epetra_MultiVector(nonlocal_map, 1));
            }
        }
    }


    Vector &
    Vector::operator = (const Vector &v)
    {
      // check equality for MPI communicators to avoid accessing a possibly
      // invalid MPI_Comm object
#ifdef DEAL_II_WITH_MPI
      const Epetra_MpiComm *my_comm = dynamic_cast<const Epetra_MpiComm *>(&vector->Comm());
      const Epetra_MpiComm *v_comm = dynamic_cast<const Epetra_MpiComm *>(&v.vector->Comm());
      const bool same_communicators = my_comm != NULL && v_comm != NULL &&
                                      my_comm->DataPtr() == v_comm->DataPtr();
#else
      const bool same_communicators = true;
#endif

      // distinguish three cases. First case: both vectors have the same
      // layout (just need to copy the local data, not reset the memory and
      // the underlying Epetra_Map). The third case means that we have to
      // rebuild the calling vector.
      if (same_communicators && vector->Map().SameAs(v.vector->Map()))
        {
          *vector = *v.vector;
          if (v.nonlocal_vector.get() != 0)
            nonlocal_vector.reset(new Epetra_MultiVector(v.nonlocal_vector->Map(), 1));
          last_action = Zero;
        }
      // Second case: vectors have the same global
      // size, but different parallel layouts (and
      // one of them a one-to-one mapping). Then we
      // can call the import/export functionality.
      else if (size() == v.size() &&
               (v.vector->Map().UniqueGIDs() || vector->Map().UniqueGIDs()))
        {
          reinit (v, false, true);
        }
      // Third case: Vectors do not have the same
      // size.
      else
        {
          vector.reset (new Epetra_FEVector(*v.vector));
          last_action = Zero;
          has_ghosts = v.has_ghosts;
        }

      if (v.nonlocal_vector.get() != 0)
        nonlocal_vector.reset(new Epetra_MultiVector(v.nonlocal_vector->Map(), 1));

      return *this;
    }



#ifdef DEAL_II_WITH_CXX11
    Vector &Vector::operator= (Vector &&v)
    {
      swap(v);
      return *this;
    }
#endif



    Vector &
    Vector::operator = (const TrilinosWrappers::Vector &v)
    {
      nonlocal_vector.reset();

      Assert (size() == v.size(), ExcDimensionMismatch(size(), v.size()));

      Epetra_Import data_exchange (vector->Map(), v.vector->Map());
      const int ierr = vector->Import(*v.vector, data_exchange, Insert);

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      last_action = Insert;

      return *this;
    }



    void
    Vector::import_nonlocal_data_for_fe (const TrilinosWrappers::SparseMatrix &m,
                                         const Vector                         &v)
    {
      Assert (m.trilinos_matrix().Filled() == true,
              ExcMessage ("Matrix is not compressed. "
                          "Cannot find exchange information!"));
      Assert (v.vector->Map().UniqueGIDs() == true,
              ExcMessage ("The input vector has overlapping data, "
                          "which is not allowed."));

      if (vector->Map().SameAs(m.trilinos_matrix().ColMap()) == false)
        {
          vector.reset (new Epetra_FEVector(
                          m.trilinos_matrix().ColMap()
                        ));
        }

      Epetra_Import data_exchange (vector->Map(), v.vector->Map());
      const int ierr = vector->Import(*v.vector, data_exchange, Insert);

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      last_action = Insert;
    }

  } /* end of namespace MPI */




  Vector::Vector ()
  {
    last_action = Zero;
    Epetra_LocalMap map (0, 0, Utilities::Trilinos::comm_self());
    vector.reset (new Epetra_FEVector(map));
  }



  Vector::Vector (const size_type n)
  {
    last_action = Zero;
    Epetra_LocalMap map ((TrilinosWrappers::types::int_type)n, 0, Utilities::Trilinos::comm_self());
    vector.reset (new Epetra_FEVector (map));
  }



  Vector::Vector (const Epetra_Map &input_map)
  {
    last_action = Zero;
    Epetra_LocalMap map (n_global_elements(input_map),
                         input_map.IndexBase(),
                         input_map.Comm());
    vector.reset (new Epetra_FEVector(map));
  }



  Vector::Vector (const IndexSet &partitioning,
                  const MPI_Comm &communicator)
  {
    last_action = Zero;
    Epetra_LocalMap map (static_cast<TrilinosWrappers::types::int_type>(partitioning.size()),
                         0,
#ifdef DEAL_II_WITH_MPI
                         Epetra_MpiComm(communicator));
#else
                         Epetra_SerialComm());
    (void)communicator;
#endif
    vector.reset (new Epetra_FEVector(map));
  }



  Vector::Vector (const VectorBase &v)
  {
    last_action = Zero;
    Epetra_LocalMap map (n_global_elements(v.vector->Map()),
                         v.vector->Map().IndexBase(),
                         v.vector->Map().Comm());
    vector.reset (new Epetra_FEVector(map));

    if (vector->Map().SameAs(v.vector->Map()) == true)
      {
        const int ierr = vector->Update(1.0, *v.vector, 0.0);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      reinit (v, false, true);

  }



  void
  Vector::reinit (const size_type n,
                  const bool    /*omit_zeroing_entries*/)
  {
    Epetra_LocalMap map ((TrilinosWrappers::types::int_type)n, 0,
                         Utilities::Trilinos::comm_self());
    vector.reset (new Epetra_FEVector (map));

    last_action = Zero;
  }



  void
  Vector::reinit (const Epetra_Map &input_map,
                  const bool     /*omit_zeroing_entries*/)
  {
    Epetra_LocalMap map (n_global_elements(input_map),
                         input_map.IndexBase(),
                         input_map.Comm());
    vector.reset (new Epetra_FEVector (map));

    last_action = Zero;
  }



  void
  Vector::reinit (const IndexSet &partitioning,
                  const MPI_Comm &communicator,
                  const bool    /*omit_zeroing_entries*/)
  {
    Epetra_LocalMap map (static_cast<TrilinosWrappers::types::int_type>(partitioning.size()),
                         0,
#ifdef DEAL_II_WITH_MPI
                         Epetra_MpiComm(communicator));
#else
                         Epetra_SerialComm());
    (void)communicator;
#endif
    vector.reset (new Epetra_FEVector(map));

    last_action = Zero;
  }



  void
  Vector::reinit (const VectorBase &v,
                  const bool        omit_zeroing_entries,
                  const bool        allow_different_maps)
  {
    // In case we do not allow to
    // have different maps, this
    // call means that we have to
    // reset the vector. So clear
    // the vector, initialize our
    // map with the map in v, and
    // generate the vector.
    if (allow_different_maps == false)
      {
        // check equality for MPI communicators
#ifdef DEAL_II_WITH_MPI
        const Epetra_MpiComm *my_comm = dynamic_cast<const Epetra_MpiComm *>(&vector->Comm());
        const Epetra_MpiComm *v_comm = dynamic_cast<const Epetra_MpiComm *>(&v.vector->Comm());
        const bool same_communicators = my_comm != NULL && v_comm != NULL &&
                                        my_comm->DataPtr() == v_comm->DataPtr();
#else
        const bool same_communicators = true;
#endif
        if (!same_communicators || local_range() != v.local_range())
          {
            Epetra_LocalMap map (global_length(*(v.vector)),
                                 v.vector->Map().IndexBase(),
                                 v.vector->Comm());
            vector.reset (new Epetra_FEVector(map));
          }
        else if (omit_zeroing_entries)
          {
            int ierr;
            Assert (vector->Map().SameAs(v.vector->Map()) == true,
                    ExcMessage ("The Epetra maps in the assignment operator ="
                                " do not match, even though the local_range "
                                " seems to be the same. Check vector setup!"));

            ierr = vector->GlobalAssemble(last_action);
            (void)ierr;
            Assert (ierr == 0, ExcTrilinosError(ierr));

            ierr = vector->PutScalar(0.0);
            Assert (ierr == 0, ExcTrilinosError(ierr));
          }
        last_action = Zero;
      }

    // Otherwise, we have to check
    // that the two vectors are
    // already of the same size,
    // create an object for the data
    // exchange and then insert all
    // the data.
    else
      {
        Assert (omit_zeroing_entries == false,
                ExcMessage ("It is not possible to exchange data with the "
                            "option 'omit_zeroing_entries' set, which would not write "
                            "elements."));

        AssertThrow (size() == v.size(),
                     ExcDimensionMismatch (size(), v.size()));

        Epetra_Import data_exchange (vector->Map(), v.vector->Map());

        const int ierr = vector->Import(*v.vector, data_exchange, Insert);
        AssertThrow (ierr == 0, ExcTrilinosError(ierr));

        last_action = Insert;
      }

  }



  Vector &
  Vector::operator = (const MPI::Vector &v)
  {
    if (size() != v.size())
      {
        Epetra_LocalMap map (n_global_elements(v.vector->Map()),
                             v.vector->Map().IndexBase(),
                             v.vector->Comm());
        vector.reset (new Epetra_FEVector(map));
      }

    reinit (v, false, true);
    return *this;
  }



  Vector &
  Vector::operator = (const Vector &v)
  {
    if (size() != v.size())
      {
        Epetra_LocalMap map (n_global_elements(v.vector->Map()),
                             v.vector->Map().IndexBase(),
                             v.vector->Comm());
        vector.reset (new Epetra_FEVector(map));
      }

    const int ierr = vector->Update(1.0, *v.vector, 0.0);
    Assert (ierr == 0, ExcTrilinosError(ierr));
    (void)ierr;

    return *this;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
