//----------------------------  swappable_vector.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  swappable_vector.templates.h  ---------------------------
#ifndef __deal2__swappable_vector_templates_h
#define __deal2__swappable_vector_templates_h


#include <lac/swappable_vector.h>
#include <fstream>
#include <iostream>



#ifdef DEAL_II_USE_MT
				 // allocate static variable
template <typename number>
Threads::ThreadManager SwappableVector<number>::thread_manager;
#endif




template <typename number>
SwappableVector<number>::SwappableVector ()
		:
		data_is_preloaded (false)
{};



template <typename number>
SwappableVector<number>::SwappableVector (const SwappableVector<number> &v) :
		Vector<number>(v),
                filename (),
                data_is_preloaded (false)
{
  Assert (v.filename == "", ExcInvalidCopyOperation());
};



template <typename number>
SwappableVector<number>::~SwappableVector ()
{
				   // if the vector was stored in a
				   // file previously, and that has
				   // not been deleted in the
				   // meantime, then we kill that file
				   // first, before killing the vector
				   // itself

  if (filename != "")
    kill_file ();
};



template <typename number>
SwappableVector<number> &
SwappableVector<number>::operator= (const SwappableVector<number> &v)
{
				   // if necessary, first delete data
  if (filename != "")
    kill_file ();

				   // if in MT mode, block all other
				   // operations
#ifdef DEAL_II_USE_MT
  lock.acquire ();
#endif
  
  Vector<number>::operator = (v);
  data_is_preloaded = false;
  
#ifdef DEAL_II_USE_MT
  lock.release ();
#endif
  
  
  return *this;
};



template <typename number>
void SwappableVector<number>::swap_out (const string &name)
{
				   // if the vector was stored in
				   // another file previously, and
				   // that has not been deleted in the
				   // meantime, then we kill that file
				   // first
  if (filename != "")
    kill_file ();
  
  filename = name;
  
  Assert (size() != 0, ExcSizeZero());

				   // if in MT mode, block all other
				   // operations
#ifdef DEAL_II_USE_MT
  lock.acquire ();
#endif
  
				   //  check that we have not called
				   //  #alert# without the respective
				   //  #reload# function
  Assert (data_is_preloaded == false, ExcInternalError());
  
  ofstream tmp_out(filename.c_str());
  block_write (tmp_out);
  tmp_out.close ();
	
  reinit (0);

#ifdef DEAL_II_USE_MT
  lock.release ();
#endif
};



template <typename number>
void SwappableVector<number>::reload () 
{
				   // if in MT mode: synchronise with
				   // possibly existing #alert# calls
#ifdef DEAL_II_USE_MT
  lock.acquire ();
#endif
				   // if data was already preloaded,
				   // then there is no more need to
				   // load it
  if (data_is_preloaded == false)
				     // reload data. note that this
				     // function also releases the
				     // lock
    reload_vector (false);
  else
    {
				       // clear flag since no more
				       // needed
      data_is_preloaded = false;

				       // release lock. the lock is
				       // also released in the other
				       // branch of the if-clause
#ifdef DEAL_II_USE_MT
      lock.release ();
#endif
    };
};



template <typename number>
void SwappableVector<number>::alert ()
{
				   // note: this function does nothing
				   // in non-MT mode
#ifdef DEAL_II_USE_MT
  
				   // synchronise with possible other
				   // invokations of this function and
				   // other functions in this class
  lock.acquire ();
  
				   // calling this function multiple
				   // times does no harm:
  if ( (data_is_preloaded == true) ||
				   // calling this function while the
				   // vector is active does no harm
				   // either
       (size() != 0))
    lock.release ();
  else
				     // data has not been preloaded so
				     // far, so go on!
    Threads::spawn (thread_manager,
		    Threads::encapsulate (&SwappableVector<number>::reload_vector)
		    .collect_args(this, true));
				   // note that reload_vector also
				   // releases the lock
#endif
};



template <typename number>
void SwappableVector<number>::reload_vector (const bool set_flag)
{
  Assert (filename != "", ExcInvalidFilename (filename));
  Assert (size() == 0, ExcSizeNonzero());
  
  ifstream tmp_in(filename.c_str());
  block_read (tmp_in);
  tmp_in.close ();

				   // release the lock that was
				   // acquired by the calling
				   // functions
#ifdef DEAL_II_USE_MT
				   // set the flag if so required
  if (set_flag)
    data_is_preloaded = true;
  lock.release ();
#else
  (void)set_flag;
#endif
};



template <typename number>
void SwappableVector<number>::kill_file () 
{
				   // if in MT mode, wait for other
				   // operations to finish first
				   // (there should be none, but who
				   // knows)
#ifdef DEAL_II_USE_MT
  lock.acquire();
#endif

  Assert (data_is_preloaded == false, ExcInternalError());
  
  if (filename != "")
    {
      int status = remove (filename.c_str());
      AssertThrow (status == 0, ExcInternalError());
  
      filename = "";
    };

#ifdef DEAL_II_USE_MT
  lock.release ();
#endif
};



template <typename number>
const string &
SwappableVector<number>::get_filename () const 
{
  return filename;
};




#endif // __deal2__swappable_vector_templates_h
