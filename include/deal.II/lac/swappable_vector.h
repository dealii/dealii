// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#ifndef __deal2__swappable_vector_h
#define __deal2__swappable_vector_h


#include <deal.II/base/config.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <string>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Vectors
 *@{
 */

/**
 * This class is a wrapper to the @p Vector class which allows to swap
 * out the data to a file and reload it later on. It handles the
 * management of the name of the file where the data is to be stored
 * temporarily and removes the file is necessary at the end of the
 * lifetime of the vector.
 *
 * There are functions to swap the data to a file, and to reload
 * it. There is also a function @p alert that can be used to signal to
 * an object of this class that the data will be needed shortly, so
 * the object can initiate that the data be loaded already. While in
 * non-multithreading mode, this function has no effect since @p reload
 * has to be called afterwards anyway. On the other hand, in
 * multithreading mode, the data is preloaded in the background using
 * a thread on its own, and may already be available at the time when
 * @p reload is called. If it is not available, @p reload waits until
 * the detached thread has loaded the data.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
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
   * data if @p v contains some, but
   * does not do so if @p v is empty
   * or has its data swapped
   * out. In the latter case, warn
   * about that. In particular do
   * not take over ownership of any
   * files that @p v might own.
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
  SwappableVector &operator = (const SwappableVector &);

  /**
   * Swap out the data of this
   * vector to the file of which
   * the name is given. It is
   * assumed that the file can be
   * overwritten if it exists, and
   * ownership of this file is
   * assumed by this object. The
   * file is deleted either upon
   * calling @p kill_file, or on
   * destruction of this object.
   *
   * The content of the vector is
   * cleared and it size is reset
   * to zero.
   *
   * If this object owns another
   * file, for example when
   * @p swap_out but no @p kill_file
   * has previously been called,
   * then that is deleted first.
   */
  void swap_out (const std::string &filename);

  /**
   * Reload the data of this vector
   * from the file to which it has
   * been stored previously using
   * @p swap_out. Since the file is
   * not deleted after reloading,
   * you can call @p reload multiple
   * times, in between you may do
   * everything with the vector,
   * including changing its size.
   *
   * This function resets the size
   * of the vector to the number of
   * elements there were upon
   * calling @p swap_out before.
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
   * that when @p reload is called,
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
   * times before calling @p reload
   * is allowed and has no effect
   * for subsequent calls. Calling
   * this function while the data
   * is still or already in memory
   * is allowed and has no effect.
   */
  void alert ();


  /**
   * Remove the file to which the
   * data has been stored the last
   * time. After this, the object
   * does not own any file any
   * more, so of course you can't
   * call @p reload no more.
   *
   * If this object does not own a
   * file, for example since
   * @p swap_out was not called, or
   * because @p kill_file has been
   * called previously, then this
   * function does nothing.
   */
  void kill_file ();

  /**
   * Return the name of the file to
   * which the data was stored the
   * last time you called
   * @p swap_out. If @p swap_out was
   * not called, or if in between
   * @p kill_file was called, then
   * the filename is an empty
   * string.
   */
  const std::string &get_filename () const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /** @addtogroup Exceptions
   * @{ */
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
                  std::string,
                  << "The filename <" << arg1
                  << "> is not a valid one here.");
  /**
   * Exception.
   */
  DeclException0 (ExcInvalidCopyOperation);
  //@}
private:
  /**
   * Name of the file to which data
   * was swapped out. If no data is
   * presently swapped out
   * (i.e. before calling
   * @p swap_out and after
   * @p kill_file), the string is
   * empty, indicating no ownership
   * of files.
   */
  std::string filename;

  /**
   * If in multithread mode, then
   * the @p alert function has
   * functionality, but needs to
   * coordinate with the @p reload
   * function. This is done through
   * the following lock.
   *
   * If not in MT mode, then the
   * class used here is empty, and
   * we can as well get away
   * with it.
   */
  Threads::Mutex lock;

  /**
   * Flag by which the @p alert
   * function signifies that the
   * data has been preloaded
   * already. This flag is always
   * @p false in non-MT mode.
   */
  bool data_is_preloaded;

  /**
   * Internal function that
   * actually reloads the
   * vector. Called from @p reload
   * and @p alert.
   *
   * The parameter specifies
   * whether the function shall set
   * @p data_is_preloaded or
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

/*@}*/
/*----------------------------   swappable_vector.h     ---------------------------*/
/* end of #ifndef __deal2__swappable_vector_h */
DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------   swappable_vector.h     ---------------------------*/
