// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2019 by the deal.II authors
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

#ifndef dealii_swappable_vector_templates_h
#define dealii_swappable_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/swappable_vector.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN


template <typename number>
SwappableVector<number>::SwappableVector()
  : data_is_preloaded(false)
{}



template <typename number>
SwappableVector<number>::SwappableVector(const SwappableVector<number> &v)
  : Vector<number>(v)
  , filename()
  , data_is_preloaded(false)
{
  Assert(v.filename.empty(), ExcInvalidCopyOperation());
}



template <typename number>
SwappableVector<number>::~SwappableVector()
{
  // if the vector was stored in a
  // file previously, and that has
  // not been deleted in the
  // meantime, then we kill that file
  // first, before killing the vector
  // itself

  if (!filename.empty())
    {
      try
        {
          kill_file();
        }
      catch (...)
        {}
    }
}



template <typename number>
SwappableVector<number> &
SwappableVector<number>::operator=(const SwappableVector<number> &v)
{
  // if necessary, first delete data
  if (!filename.empty())
    kill_file();

  // if in MT mode, block all other
  // operations. if not in MT mode,
  // this is a no-op
  std::lock_guard<std::mutex> lock(this->lock);

  Vector<number>::operator=(v);
  data_is_preloaded       = false;

  return *this;
}



template <typename number>
void
SwappableVector<number>::swap_out(const std::string &name)
{
  // if the vector was stored in
  // another file previously, and
  // that has not been deleted in the
  // meantime, then we kill that file
  // first
  if (!filename.empty())
    kill_file();

  filename = name;

  Assert(this->size() != 0, ExcSizeZero());

  // if in MT mode, block all other
  // operations. if not in MT mode,
  // this is a no-op
  std::lock_guard<std::mutex> lock(this->lock);

  //  check that we have not called
  //  @p alert without the respective
  //  @p reload function
  Assert(data_is_preloaded == false, ExcInternalError());

  std::ofstream tmp_out(filename.c_str());
  this->block_write(tmp_out);
  tmp_out.close();

  this->reinit(0);
}



template <typename number>
void
SwappableVector<number>::reload()
{
  // if in MT mode: synchronize with
  // possibly existing @p alert
  // calls. if not in MT mode, this
  // is a no-op
  lock.lock();

  // if data was already preloaded,
  // then there is no more need to
  // load it
  if (data_is_preloaded == false)
    // reload data. note that this
    // function also releases the
    // lock
    reload_vector(false);
  else
    {
      // clear flag since no more
      // needed
      data_is_preloaded = false;

      // release lock. the lock is
      // also released in the other
      // branch of the if-clause
      lock.unlock();
    }
}



template <typename number>
void
SwappableVector<number>::alert()
{
#ifndef DEAL_II_WITH_THREADS
  // note: this function does nothing
  // in non-MT mode
  return;
#else

  // synchronize with possible other
  // invocations of this function and
  // other functions in this class
  lock.lock();

  // calling this function multiple
  // times does no harm:
  if ((data_is_preloaded == true) ||
      // calling this function while the
      // vector is active does no harm
      // either
      (this->size() != 0))
    lock.unlock();
  else
    // data has not been preloaded so
    // far, so go on! For this, start
    // a detached thread
    Threads::new_thread(&SwappableVector<number>::reload_vector, *this, true);
    // note that reload_vector also
    // releases the lock
#endif
}



template <typename number>
void
SwappableVector<number>::reload_vector(const bool set_flag)
{
  (void)set_flag;

  Assert(!filename.empty(), ExcInvalidFilename(filename));
  Assert(this->size() == 0, ExcSizeNonzero());

  std::ifstream tmp_in(filename.c_str());
  this->block_read(tmp_in);
  tmp_in.close();

#ifdef DEAL_II_WITH_THREADS
  // release the lock that was acquired by the calling functions

  // set the flag if so required
  if (set_flag)
    data_is_preloaded = true;
  lock.unlock();
#endif
}



template <typename number>
void
SwappableVector<number>::kill_file()
{
  // if in MT mode, wait for other
  // operations to finish first
  // (there should be none, but who
  // knows). if not in MT mode,
  // this is a no-op
  std::lock_guard<std::mutex> lock(this->lock);

  // this is too bad: someone
  // requested the vector in advance,
  // but never got to fetch it. this
  // is most probably an error, not?
  Assert(data_is_preloaded == false, ExcInternalError());

  if (!filename.empty())
    {
      int status = std::remove(filename.c_str());
      (void)status;
      AssertThrow(status == 0, ExcIO());

      filename = "";
    };
}



template <typename number>
const std::string &
SwappableVector<number>::get_filename() const
{
  return filename;
}



template <typename number>
std::size_t
SwappableVector<number>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(filename) + sizeof(lock) +
          MemoryConsumption::memory_consumption(data_is_preloaded) +
          Vector<number>::memory_consumption());
}



DEAL_II_NAMESPACE_CLOSE

#endif // dealii_swappable_vector_templates_h
