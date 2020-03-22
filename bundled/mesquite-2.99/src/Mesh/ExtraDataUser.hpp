/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file ExtraDataUser.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_EXTRA_DATA_USER_HPP
#define MSQ_EXTRA_DATA_USER_HPP

#include "Mesquite.hpp"
#include "ExtraData.hpp"
#include <assert.h>

namespace MESQUITE_NS {

template <typename T> class ExtraUserData;

/**\brief Manage extra data attached to PatchData instances
 *
 * This class manages the details of using the ExtraData mechanism
 * for attaching data to and/or observing a PatchData.  The template
 * parameter is the type definint the data to be stored on the 
 * PatchData.  
 * 
 * To use this class, define a type (struct or whatever) to contain
 * the data that will be stored on PatchData instances, and create
 * a subclass of this class, where the template parameter for this
 * class is the data type.  Provide implementations of the pure
 * virtual (abstract) methods for notification of changes to the PatchData.
 */
template <typename T> class ExtraDataUser
{
  public:
  
    ExtraDataUser();
    
    virtual ~ExtraDataUser();
  
    typedef T DataType;
    
    /** returns null if data hasn't been set via set_data */
    T* get_data_ptr( PatchData& patch );
    
    /** creates data if doesn't exist */
    T& get_data( PatchData& patch );
    
    void set_data( PatchData& patch, const T& data );
    
    void notify_patch_destroyed( ExtraUserData<T>* data );
    
    virtual void notify_patch_destroyed( T& data ) = 0;
    
    virtual void notify_new_patch( PatchData& patch, T& data ) = 0;
    
    virtual void notify_sub_patch( PatchData& patch, 
                                   T& data, 
                                   PatchData& sub_patch, 
                                   const size_t* vertex_index_map,
                                   const size_t* element_index_map,
                                   MsqError& err ) = 0;
  
  private:
  
    ExtraUserData<T>* listHead;
};

template <typename T> class ExtraUserData : public ExtraData
{
  public:
    
    ExtraDataUser<T>* dataOwner;
    ExtraUserData<T>* userNext;
    T userData;
    
    ExtraUserData( PatchData& patch, 
                   ExtraDataUser<T>* owner, 
                   ExtraUserData<T>* next,
                   const T& data )
      : ExtraData(patch), dataOwner(owner), userNext(next), userData(data)
      {}
    
    ExtraUserData( PatchData& patch, 
                   ExtraDataUser<T>* owner, 
                   ExtraUserData<T>* next )
      : ExtraData(patch), dataOwner(owner), userNext(next)
      {}
    
    virtual void notify_patch_destroyed();
    
    virtual void notify_new_patch();
    
    virtual void notify_sub_patch( PatchData& sub_patch, 
                                   const size_t* vtx_index_map,
                                   const size_t* elm_index_map,
                                   MsqError& err );
};
    
template <typename T> 
void ExtraUserData<T>::notify_patch_destroyed()
  { dataOwner->notify_patch_destroyed( this ); }

template <typename T> 
void ExtraUserData<T>::notify_new_patch()
  { dataOwner->notify_new_patch( *get_patch_data(), userData ); }

template <typename T> 
void ExtraUserData<T>::notify_sub_patch( PatchData& sub, 
                                         const size_t* vertex_map,
                                         const size_t* element_map,
                                         MsqError& err )
  { dataOwner->notify_sub_patch( *get_patch_data(), userData, sub, vertex_map, element_map, err ); }

template <typename T> 
ExtraDataUser<T>::ExtraDataUser() : listHead(0) {}

template <typename T> 
ExtraDataUser<T>::~ExtraDataUser()
{
  while (ExtraUserData<T>* dead_ptr = listHead) {
    listHead = dead_ptr->userNext;
    dead_ptr->userNext = 0;
    delete dead_ptr;
  }
}

template <typename T> 
T* ExtraDataUser<T>::get_data_ptr( PatchData& patch )
{
  for (ExtraUserData<T>* ptr = listHead; ptr; ptr = ptr->userNext)
    if (ptr->get_patch_data() == &patch)
      return &(ptr->userData);
  return 0;
}

template <typename T> 
T& ExtraDataUser<T>::get_data( PatchData& patch )
{
  T* ptr = get_data_ptr( patch );
  if (!ptr) {
    listHead = new ExtraUserData<T>( patch, this, listHead );
    ptr = &(listHead->userData);
  }
  return *ptr;
}

template <typename T> 
void ExtraDataUser<T>::set_data( PatchData& patch, const T& data )
{
  T* ptr = get_data_ptr( patch );
  if (ptr) 
    *ptr = data;
  else
    listHead = new ExtraUserData<T>( patch, this, listHead, data );
}

template <typename T> 
void ExtraDataUser<T>::notify_patch_destroyed( ExtraUserData<T>* data )
{
    // remove from list
  assert(listHead != 0);
  if (listHead == data)
    listHead = data->userNext;
  else {
    ExtraUserData<T>* prev;
    for (prev = listHead; prev->userNext != data; prev = prev->userNext)
      assert( prev->userNext != 0 );
    prev->userNext = data->userNext;
  }
  data->userNext = 0;
  
    // notify concrete class of event
  notify_patch_destroyed( data->userData );
  
  delete data;
}

} // namespace Mesquite

#endif
