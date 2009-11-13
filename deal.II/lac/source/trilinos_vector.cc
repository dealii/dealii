//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/trilinos_vector.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/trilinos_sparse_matrix.h>
#  include <cmath>
#  include <Epetra_Import.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace MPI
  {


    Vector::Vector ()
    {
      last_action = Zero;
      vector = std::auto_ptr<Epetra_FEVector> 
	(new Epetra_FEVector(Epetra_Map(0,0,0,Utilities::Trilinos::comm_self())));
    }


  
    Vector::Vector (const Epetra_Map &input_map)
    {
      reinit (input_map);
    }


  
    Vector::Vector (const IndexSet &parallel_partitioner,
		    const MPI_Comm &communicator)
    {
      reinit (parallel_partitioner, communicator);
    }
  

  
    Vector::Vector (const Vector &v)
                    :
                    VectorBase()
    {
      last_action = Zero;
      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(*v.vector));
    }



    Vector::Vector (const Epetra_Map &input_map,
		    const VectorBase &v)
                    :
                    VectorBase()
    {
      AssertThrow (input_map.NumGlobalElements() == v.vector->Map().NumGlobalElements(),
		   ExcDimensionMismatch (input_map.NumGlobalElements(),
					 v.vector->Map().NumGlobalElements()));

      last_action = Zero;
      
      if (input_map.SameAs(v.vector->Map()) == true)
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(*v.vector));
      else
	{
	  vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(input_map));
	  reinit (v, false, true);
	}
    }



    Vector::Vector (const IndexSet   &parallel_partitioner,
		    const VectorBase &v,
		    const MPI_Comm   &communicator)
                    :
                    VectorBase()
    {
      AssertThrow ((int)parallel_partitioner.size() == v.vector->Map().NumGlobalElements(),
		   ExcDimensionMismatch (parallel_partitioner.size(),
					 v.vector->Map().NumGlobalElements()));

      last_action = Zero;
      
      vector = std::auto_ptr<Epetra_FEVector>
	(new Epetra_FEVector(parallel_partitioner.make_trilinos_map(communicator,
								    true)));
      reinit (v, false, true);
    }



    Vector::~Vector ()
    {}
    


    void
    Vector::reinit (const Epetra_Map &input_map,
		    const bool        fast)
    {
      if (!vector->Map().SameAs(input_map))
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(input_map));
      else if (fast == false)
	{
	  const int ierr = vector->PutScalar(0.);
	  Assert (ierr == 0, ExcTrilinosError(ierr));
	}
  
      last_action = Zero;
    }
    


    void
    Vector::reinit (const IndexSet &parallel_partitioner,
		    const MPI_Comm &communicator,
		    const bool      fast)
    {
      Epetra_Map map = parallel_partitioner.make_trilinos_map (communicator,
							       true);
      reinit (map, fast);
    }



    void
    Vector::reinit (const VectorBase &v,
		    const bool        fast,
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
	  if (vector->Map().SameAs(v.vector->Map()) == false)
	    {
	      vector.reset();

	      vector = std::auto_ptr<Epetra_FEVector> 
		(new Epetra_FEVector(v.vector->Map()));
	      last_action = Zero;
	    }
	  else if (fast == false)
	    {
					       // old and new vectors
					       // have exactly the
					       // same map, i.e. size
					       // and parallel
					       // distribution
	      int ierr;
	      ierr = vector->GlobalAssemble (last_action);	
	      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	      ierr = vector->PutScalar(0.0);
	      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	      last_action = Zero;
	    }
	}

					// Otherwise, we have to check
					// that the two vectors are
					// already of the same size,
					// create an object for the data
					// exchange and then insert all
					// the data. The first assertion
					// is only a check whether the
					// user knows what she is doing.
      else
        {
	  Assert (fast == false,
		  ExcMessage ("It is not possible to exchange data with the "
			      "option fast set, which would not write "
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
    Vector::operator = (const Vector &v)
    {
				// distinguish three cases. First case: both
				// vectors have the same layout (just need to
				// copy the local data). Second case: vectors
				// have the same size, but different
				// layout. The calling vector has a wider
				// local range than the input vector, and the
				// input vector has a 1-to-1 map (need to
				// import data). The third case means that we
				// have to rebuild the calling vector.
      if (size() == v.size() && 
	  local_range() == v.local_range())
	{
	  Assert (vector->Map().SameAs(v.vector->Map()) == true,
		  ExcMessage ("The Epetra maps in the assignment operator ="
			      " do not match, even though the local_range "
			      " seems to be the same. Check vector setup!"));

	  const int ierr = vector->Update(1.0, *v.vector, 0.0);
	  AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	  last_action = Zero;
	}
      else if (size() == v.size() && local_range().first<=v.local_range().first &&
	       local_range().second>=v.local_range().second &&
	       v.vector->Map().UniqueGIDs())
	{
	  Epetra_Import data_exchange (vector->Map(), v.vector->Map());
	  const int ierr = vector->Import(*v.vector, data_exchange, Insert);
	  AssertThrow (ierr == 0, ExcTrilinosError(ierr));
	}
      else
	{
	  vector.reset();
	  vector = std::auto_ptr<Epetra_FEVector> 
	                    (new Epetra_FEVector(*v.vector));
	  last_action = Zero;
	}

      return *this;
    }



    Vector &
    Vector::operator = (const TrilinosWrappers::Vector &v)
    {
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

      if (vector->Map().SameAs(m.col_partitioner()) == false)
	{
	  Epetra_Map map = m.col_partitioner();
	  vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
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
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
  }



  Vector::Vector (const unsigned int n)
  {
    last_action = Zero;
    Epetra_LocalMap map ((int)n, 0, Utilities::Trilinos::comm_self());
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector (map));
  }



  Vector::Vector (const Epetra_Map &input_map)
  {
    last_action = Zero;
    Epetra_LocalMap map (input_map.NumGlobalElements(), 
			 input_map.IndexBase(), 
			 input_map.Comm());
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
  }



  Vector::Vector (const VectorBase &v)
  {
    last_action = Zero;
    Epetra_LocalMap map (v.vector->Map().NumGlobalElements(), 
			 v.vector->Map().IndexBase(), 
			 v.vector->Map().Comm());
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));

    if (vector->Map().SameAs(v.vector->Map()) == true)
      {
	const int ierr = vector->Update(1.0, *v.vector, 0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    else
      reinit (v, false, true);

  }



  void
  Vector::reinit (const unsigned int n,
		  const bool         fast)
  {
    if (size() != n)
      {
	vector.reset();

	Epetra_LocalMap map ((int)n, 0,
			     Utilities::Trilinos::comm_self());
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector (map));
      }
    else if (fast == false)
      {
	int ierr;
	ierr = vector->GlobalAssemble(last_action);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	ierr = vector->PutScalar(0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }
    
    last_action = Zero;
  }



  void
  Vector::reinit (const Epetra_Map &input_map,
                  const bool        fast)
  {
    if (vector->Map().NumGlobalElements() != input_map.NumGlobalElements())
      {
	vector.reset();
	Epetra_LocalMap map (input_map.NumGlobalElements(),
			     input_map.IndexBase(),
			     input_map.Comm());
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector (map));
      }
    else if (fast == false)
      {
	int ierr;
	ierr = vector->GlobalAssemble(last_action);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	ierr = vector->PutScalar(0.0);
	AssertThrow (ierr == 0, ExcTrilinosError(ierr));
      }

    last_action = Zero;
  }


  void
  Vector::reinit (const VectorBase &v,
		  const bool        fast,
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
	if (local_range() != v.local_range())
	  {
	    vector.reset();
	    Epetra_LocalMap map (v.vector->GlobalLength(),
				 v.vector->Map().IndexBase(),
				 v.vector->Comm());
	    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
	  }
	else
	  {
	    int ierr;
	    Assert (vector->Map().SameAs(v.vector->Map()) == true,
		    ExcMessage ("The Epetra maps in the assignment operator ="
				" do not match, even though the local_range "
				" seems to be the same. Check vector setup!"));

	    ierr = vector->GlobalAssemble(last_action);
	    AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	    ierr = vector->PutScalar(0.0);
	    AssertThrow (ierr == 0, ExcTrilinosError(ierr));
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
	Assert (fast == false,
		ExcMessage ("It is not possible to exchange data with the "
			    "option fast set, which would not write "
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
	vector.reset();
	Epetra_LocalMap map (v.vector->Map().NumGlobalElements(), 
			     v.vector->Map().IndexBase(),
			     v.vector->Comm());
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
      }

    reinit (v, false, true);
    return *this;
  }



  Vector &
  Vector::operator = (const Vector &v)
  {
    if (size() != v.size())
      {
	Epetra_LocalMap map (v.vector->Map().NumGlobalElements(), 
			     v.vector->Map().IndexBase(),
			     v.vector->Comm());
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
      }

    const int ierr = vector->Update(1.0, *v.vector, 0.0);
    AssertThrow (ierr == 0, ExcTrilinosError(ierr));	

    return *this;
  }
  
}

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_USE_TRILINOS
