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


#include <base/utilities.h>
#include <lac/trilinos_vector.h>

#include <cmath>

#ifdef DEAL_II_USE_TRILINOS

#include <Epetra_Import.h>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  namespace MPI
  {


    Vector::Vector ()
		    :
		    communicator (Utilities::Trilinos::
				  duplicate_communicator(Utilities::Trilinos::comm_self())),
		    map (0, 0, *communicator)
    {
      last_action = Zero;
      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
    }


  
    Vector::Vector (const Epetra_Map &InputMap)
                    :
		    communicator (Utilities::Trilinos::
				  duplicate_communicator (InputMap.Comm())),
		    map (Utilities::Trilinos::
			 duplicate_map (InputMap, *communicator))
    {
      last_action = Zero;
      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
    }
  

  
    Vector::Vector (const Vector &v)
                    :
                    VectorBase(),
		    communicator (Utilities::Trilinos::
				  duplicate_communicator (*v.communicator)),
		    map (Utilities::Trilinos::duplicate_map (v.map,
							     *communicator))
    {
      last_action = Zero;
      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(*v.vector));
    }



    Vector::Vector (const Epetra_Map &InputMap,
		    const VectorBase &v)
                    :
                    VectorBase(),
		    communicator (Utilities::Trilinos::
				  duplicate_communicator (InputMap.Comm())),
		    map (Utilities::Trilinos::
			 duplicate_map (InputMap, *communicator))
    {
      AssertThrow (map.NumGlobalElements() == v.vector->Map().NumGlobalElements(),
		   ExcDimensionMismatch (map.NumGlobalElements(),
					 v.vector->Map().NumGlobalElements()));

      last_action = Zero;
      
      if (map.SameAs(v.vector->Map()) == true)
	vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(*v.vector));
      else
	{
	  vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
	  reinit (v, false, true);
	}
    }



    void
    Vector::reinit (const Epetra_Map &input_map,
		    const bool        fast)
    {
      vector.reset();

      communicator.reset (Utilities::Trilinos::
			  duplicate_communicator(input_map.Comm()));
      map = Utilities::Trilinos::duplicate_map (input_map, *communicator);
      
      
      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));

      if (fast == false)
	{
	  const int ierr = vector->PutScalar(0.);
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
	  if (vector->Map().SameAs(v.vector->Map()) == false)
	    {
	      vector.reset();
	      communicator.reset (Utilities::Trilinos::
				  duplicate_communicator(v.trilinos_vector()
							 .Map().Comm()));
	      map = Utilities::Trilinos::
		    duplicate_map (v.trilinos_vector().Map(), *communicator);

	      vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
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
      if (local_range() == v.local_range())
	{
	  Assert (vector->Map().SameAs(v.vector->Map()) == true,
		  ExcMessage ("The Epetra maps in the assignment operator ="
			      " do not match, even though the local_range "
			      " seems to be the same. Check vector setup!"));

	  const int ierr = vector->Update(1.0, *v.vector, 0.0);
	  AssertThrow (ierr == 0, ExcTrilinosError(ierr));

	  last_action = Zero;
	}
      else
	{
	  vector.reset();

	  communicator.reset (Utilities::Trilinos::
			      duplicate_communicator(v.trilinos_vector().
						     Map().Comm()));
	  map = Utilities::Trilinos::
		duplicate_map (v.trilinos_vector().Map(),
			       *communicator);

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
	  map = m.col_partitioner();
	  vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
	}

      Epetra_Import data_exchange (vector->Map(), v.vector->Map());
      const int ierr = vector->Import(*v.vector, data_exchange, Insert);

      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      last_action = Insert;
    }

  } /* end of namespace MPI */




  Vector::Vector ()
		  :
		  map (0, 0,
		       Utilities::Trilinos::comm_self())
  {
    last_action = Zero;
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
  }



  Vector::Vector (const unsigned int n)
                 :
		  map ((int)n, 0,
		       Utilities::Trilinos::comm_self())
  {
    last_action = Zero;
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector (map));
  }



  Vector::Vector (const Epetra_Map &InputMap)
                 :
                 map (InputMap.NumGlobalElements(), InputMap.IndexBase(), 
		      InputMap.Comm())
  {
    last_action = Zero;
    vector = std::auto_ptr<Epetra_FEVector> (new Epetra_FEVector(map));
  }



  Vector::Vector (const VectorBase &v)
                 :
                 map (v.vector->Map().NumGlobalElements(), 0, v.vector->Comm())
  {
    last_action = Zero;
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

	map = Epetra_LocalMap ((int)n, 0,
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
    if (map.NumGlobalElements() != input_map.NumGlobalElements())
      {
	vector.reset();
	map = Epetra_LocalMap (input_map.NumGlobalElements(),
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
	    map = Epetra_LocalMap (v.vector->GlobalLength(),
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
	map = Epetra_LocalMap (v.vector->Map().NumGlobalElements(), 
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
	map = Epetra_LocalMap (v.vector->Map().NumGlobalElements(), 
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
