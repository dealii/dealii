//----------------------------  swappable_vector.h  ---------------------------
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
//----------------------------  swappable_vector.h  ---------------------------
#ifndef __deal2__swappable_vector_h
#define __deal2__swappable_vector_h


#include <lac/vector.h>
#include <base/thread_management.h>
#include <string>



/**
 * This class is a wrapper to the @p{Vector} class which allows to swap
 * out the data to a file and reload it later on. It handles the
 * management of the name of the file where the data is to be stored
 * temporarily and removes the file is necessary at the end of the
 * lifetime of the vector.
 *
 * There are functions to swap the data to a file, and to reload
 * it. There is also a function @p{alert} that can be used to signal to
 * an object of this class that the data will be needed shortly, so
 * the object can initiate that the data be loaded already. While in
 * non-multithreading mode, this function has no effect since @p{reload}
 * has to be called afterwards anyway. On the other hand, in
 * multithreading mode, the data is preloaded in the background using
 * a thread on its own, and may already be available at the time when
 * @p{reload} is called. If it is not available, @p{reload} waits until
 * the detached thread has loaded the data.
 *
 *
 * @sect2{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Wolfgang Bangerth, 1999, 2000
 */
template <typename number>
class SwappableVector : public Vector<number>
{
  public:
				     /**
				      * Constructor. Does nothing
				      * apart from calling the
				      * constructor of the underlying
				      * class.
				      */
    SwappableVector ();

				     /**
				      * Copy constructor. Copies the
				      * data if @p{v} contains some, but
				      * does not do so if @p{v} is empty
				      * or has its data swapped
				      * out. In the latter case, warn
				      * about that. In particular do
				      * not take over ownership of any
				      * files that @p{v} might own (as,
				      * for example, @p{C++}'s
				      * @p{auto_ptr} objects would do).
				      */
    SwappableVector (const SwappableVector &v);

				     /**
				      * Destructor. If this class
				      * still owns a file to which
				      * temporary data was stored,
				      * then it is deleted.
				      */
    virtual ~SwappableVector ();
    
				     /**
				      * Copy operator. Do mostly the
				      * same as the copy constructor
				      * does; if necessary, delete
				      * temporary files owned by this
				      * object at first.
				      */
    SwappableVector & operator = (const SwappableVector &);

				     /**
				      * Swap out the data of this
				      * vector to the file of which
				      * the name is given. It is
				      * assumed that the file can be
				      * overwritten if it exists, and
				      * ownership of this file is
				      * assumed by this object. The
				      * file is deleted either upon
				      * calling @p{kill_file}, or on
				      * destruction of this object.
				      *
				      * The content of the vector is
				      * cleared and it size is reset
				      * to zero.
				      *
				      * If this object owns another
				      * file, for example when
				      * @p{swap_out} but no @p{kill_file}
				      * has previously been called,
				      * then that is deleted first.
				      */
    void swap_out (const string &filename);

				     /**
				      * Reload the data of this vector
				      * from the file to which it has
				      * been stored previously using
				      * @p{swap_out}. Since the file is
				      * not deleted after reloading,
				      * you can call @p{reload} multiple
				      * times, in between you may do
				      * everything with the vector,
				      * including changing its size.
				      *
				      * This function resets the size
				      * of the vector to the number of
				      * elements there were upon
				      * calling @p{swap_out} before.
				      */
    void reload ();

				     /**
				      * Calling this function can be
				      * used to alert this vector that
				      * it will need to reload its
				      * data soon. It has no effect if
				      * the data of the vector is
				      * already present, and it has no
				      * effect in single-thread mode
				      * as well, but in multithread
				      * mode, it spawns another thread
				      * that reads the data in
				      * parallel to the usual
				      * execution of the program, such
				      * that when @p{reload} is called,
				      * the data may eventually be
				      * available already. It might
				      * therefore be wirthwhile to
				      * call this function some time
				      * in advance, if you know that
				      * the data will be needed, and
				      * loading takes some time, for
				      * instance if the file to which
				      * the data was written is not in
				      * a local tmp directory.
				      *
				      * Calling this function multiple
				      * times before calling @p{reload}
				      * is allowed and has no effect
				      * for subsequent calls. Calling
				      * this function while the data
				      * is still or already in memory
				      * is allowed and has no effect.
				      *
				      * It is noted that versions of
				      * gcc up to at least version
				      * 2.95.2 are not thread-safe in
				      * the C++ standard
				      * library. Thus, parallel in-
				      * and/or output may sometimes
				      * lead to unreproducable
				      * crashes, segmentation faults,
				      * or bus errors of programs. You
				      * should therefore only call
				      * this function if either you
				      * are sure that the I/O library
				      * distributed with your compiler
				      * supports parallel
				      * writing/reading, or if you are
				      * sure that during the run-time
				      * of this function no other I/O
				      * operations are active.
				      */
    void alert ();
    

				     /**
				      * Remove the file to which the
				      * data has been stored the last
				      * time. After this, the object
				      * does not own any file any
				      * more, so of course you can't
				      * call @p{reload} no more.
				      *
				      * If this object does not own a
				      * file, for example since
				      * @p{swap_out} was not called, or
				      * because @p{kill_file} has been
				      * called previously, then this
				      * function does nothing.
				      */
    void kill_file ();

				     /**
				      * Return the name of the file to
				      * which the data was stored the
				      * last time you called
				      * @p{swap_out}. If @p{swap_out} was
				      * not called, or if in between
				      * @p{kill_file} was called, then
				      * the filename is an empty
				      * string.
				      */
    const string & get_filename () const;
   
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception.
				      */ 
    DeclException0 (ExcSizeZero);
				     /**
				      * Exception.
				      */ 
    DeclException0 (ExcSizeNonzero);
				     /**
				      * Exception.
				      */ 
    DeclException1 (ExcInvalidFilename,
		    string,
		    << "The filename <" << arg1
		    << "> is not a valid one here.");
				     /**
				      * Exception.
				      */ 
    DeclException0 (ExcInvalidCopyOperation);
    
  private:
				     /**
				      * Name of the file to which data
				      * was swapped out. If no data is
				      * presently swapped out
				      * (i.e. before calling
				      * @p{swap_out} and after
				      * @p{kill_file}), the string is
				      * empty, indicating no ownership
				      * of files.
				      */
    string filename;

				     /**
				      * If in multithread mode, then
				      * the @p{alert} function has
				      * functionality, but needs to
				      * coordinate with the @p{reload}
				      * function. This is done through
				      * the following lock.
				      *
				      * If not in MT mode, then the
				      * class used here is empty, and
				      * we can as well get away
				      * with it.
				      */
    Threads::ThreadMutex lock;

				     /**
				      * Declare a thread manager that
				      * is used to spawn threads in
				      * @p{alert} detached.
				      *
				      * If not in MT mode, then the
				      * class used here is empty, and
				      * we can as well get away
				      * without it.
				      */
    static Threads::ThreadManager thread_manager;

				     /**
				      * Flag by which the @p{alert}
				      * function signifies that the
				      * data has been preloaded
				      * already. This flag is always
				      * @p{false} in non-MT mode.
				      */
    bool data_is_preloaded;

				     /**
				      * Internal function that
				      * actually reloads the
				      * vector. Called from @p{reload}
				      * and @p{alert}.
				      *
				      * The parameter specifies
				      * whether the function shall set
				      * @p{data_is_preloaded} or
				      * not. The calling functions
				      * can't sometimes do this
				      * themselves, since they call
				      * this function detached, so
				      * this function has to signal
				      * success itself if this is
				      * required.
				      */
    void reload_vector (const bool set_flag);
};



/*----------------------------   swappable_vector.h     ---------------------------*/
/* end of #ifndef __deal2__swappable_vector_h */
#endif
/*----------------------------   swappable_vector.h     ---------------------------*/
