//----------------------------  thread_management.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  thread_management.h  ---------------------------
#ifndef __deal2__thread_management_h
#define __deal2__thread_management_h


#ifdef DEAL_II_USE_MT

#include <base/exceptions.h>
#include <ace/Thread_Manager.h>
#include <ace/Synch.h>

#include <string>
#include <utility>
#include <vector>
#include <iterator>




namespace Threads 
{
				   // forward declarations
  class FunDataBase;


/**
 * Class used to store a pointer temporarily and delete the object
 * pointed to upon destruction of this object. For more information on
 * use and internals of this class see the report on this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  class FunEncapsulation
  {
    private:
				       /**
					* Default constructor. Construct
					* the object empty, i.e. set
					* #data==0#. Since this is not
					* very useful, disallow it by
					* declaring this constructor
					* #private#.
					*/
      FunEncapsulation ();

    public:
				       /**
					* Copy constructor. Clone the
					* object pointed to by
					* #fun_data.fun_data_base#.
					*/
      FunEncapsulation (const FunEncapsulation &fun_encapsulation);

				       /**
					* This is the usual
					* constructor. Set #fun_data_base# to
					* #fun_data_base#. This is what
					* the #fun_data_*# functions
					* use.
					*/
      FunEncapsulation (FunDataBase *fun_data_base);

				       /**
					* Destructor. Delete the object
					* pointed to by #fun_data_base#.
					*/
      ~FunEncapsulation ();

				       /**
					* Copy another object of this
					* type by cloning its #fun_data_base#
					* object.
					*/
      const FunEncapsulation & operator = (const FunEncapsulation &fun_encapsulation);
    
				       /**
					* Pointer to the object which
					* contains all the parameters.
					*/
      const FunDataBase * fun_data_base;
  };



/**
 * Abstract base class for those classes that actually store
 * parameters of functions. For more information on use and internals
 * of this class see the report on this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  class FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a
					* function that satifies the
					* requirements of thread entry
					* points.
					*/
      typedef void * (*ThreadEntryPoint) (void *);

				       /**
					* Default constructor. Store the
					* pointer to the function which
					* we will use as thread entry
					* point for the derived class.
					*/
      FunDataBase (const ThreadEntryPoint thread_entry_point);

				       /**
					* Copy constructor.
					*/
      FunDataBase (const FunDataBase &);

				       /**
					* Destructor. Needs to be
					* virtual to make destruction of
					* derived classes through base
					* class pointers.
					*/
      virtual ~FunDataBase ();

				       /**
					* Virtual constructor. Needed to
					* copy an object of which we
					* only have a pointer to the
					* base class. Copying such
					* objects is necessary to
					* guarantee memory consistency.
					*/
      virtual FunDataBase * clone () const = 0;

				       /**
					* Lock to be used when starting
					* a thread and which is released
					* after the data of this object
					* is copied and therefore no
					* more needed. This ensures that
					* no data is deleted when it is
					* still in use.
					*/
      mutable ACE_Thread_Mutex lock;
    
    private:
				       /**
					* Pointer to the thread entry
					* point function. The address of
					* that function is passed from
					* the derived classes to the
					* constructor of this class.
					*/
      ThreadEntryPoint thread_entry_point;

				       /**
					* Make the thread starter
					* function a friend, since it
					* needs to have access to the
					* #thread_entry_point# variable.
					*/
      friend void spawn (ACE_Thread_Manager       &thread_manager,
			 const FunEncapsulation   &fun_data);
  };



/**
 * Class to store the parameters of a void function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename RetType>
  class FunData0 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) ();

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData0 (FunPtr fun_ptr);

				       /**
					* Copy constructor.
					*/
      FunData0 (const FunData0 &fun_data0);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData0<RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    *
					    * Since the function in
					    * question here does not
					    * take parameters, this
					    * function also does
					    * nothing. It is only
					    * present for
					    * orthogonality of thread
					    * creation.
					    */
	  FunEncapsulation collect_args ();
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a unary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename Arg1, typename RetType>
  class FunData1 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) (Arg1);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData1 (FunPtr fun_ptr,
		Arg1   arg1);

				       /**
					* Copy constructor.
					*/
      FunData1 (const FunData1 &fun_data1);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;

				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData1<Arg1,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Arg1 arg1);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a binary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename Arg1, typename Arg2, typename RetType>
  class FunData2 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) (Arg1, Arg2);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData2 (FunPtr fun_ptr,
		Arg1   arg1,
		Arg2   arg2);

				       /**
					* Copy constructor.
					*/
      FunData2 (const FunData2 &fun_data2);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData2<Arg1,Arg2,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Arg1 arg1,
					 Arg2 arg2);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };

  
/**
 * Class to store the parameters of a ternary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  class FunData3 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) (Arg1, Arg2, Arg3);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData3 (FunPtr fun_ptr,
		Arg1   arg1,
		Arg2   arg2,
		Arg3   arg3);

				       /**
					* Copy constructor.
					*/
      FunData3 (const FunData3 &fun_data3);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData3<Arg1,Arg2,Arg3,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };


/**
 * Class to store the parameters of a quaternary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  class FunData4 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) (Arg1, Arg2, Arg3, Arg4);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData4 (FunPtr fun_ptr,
		Arg1   arg1,
		Arg2   arg2,
		Arg3   arg3,
		Arg4   arg4);

				       /**
					* Copy constructor.
					*/
      FunData4 (const FunData4 &fun_data4);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      Arg4   arg4;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3,
					 Arg4 arg4);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a quintary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  class FunData5 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (*FunPtr) (Arg1, Arg2, Arg3, Arg4, Arg5);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      FunData5 (FunPtr fun_ptr,
		Arg1   arg1,
		Arg2   arg2,
		Arg3   arg3,
		Arg4   arg4,
		Arg5   arg5);

				       /**
					* Copy constructor.
					*/
      FunData5 (const FunData5 &fun_data5);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      Arg4   arg4;
      Arg5   arg5;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3,
					 Arg4 arg4,
					 Arg5 arg5);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };

  
  
/**
 * Class to store the parameters of a void function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename RetType>
  class MemFunData0 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) ();

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData0 (FunPtr fun_ptr,
		   Class *object);

				       /**
					* Copy constructor.
					*/
      MemFunData0 (const MemFunData0 &fun_data0);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
 
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData0<Class,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a unary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename Arg1, typename RetType>
  class MemFunData1 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) (Arg1);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData1 (FunPtr fun_ptr,
		   Class *object,
		   Arg1   arg1);

				       /**
					* Copy constructor.
					*/
      MemFunData1 (const MemFunData1 &fun_data1);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
      
				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;

				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData1<Class,Arg1,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object,
					 Arg1 arg1);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a binary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename Arg1, typename Arg2, typename RetType>
  class MemFunData2 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) (Arg1, Arg2);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData2 (FunPtr fun_ptr,
		   Class *object,
		   Arg1   arg1,
		   Arg2   arg2);

				       /**
					* Copy constructor.
					*/
      MemFunData2 (const MemFunData2 &fun_data2);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
 
				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData2<Class,Arg1,Arg2,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object,
					 Arg1 arg1,
					 Arg2 arg2);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };

  
/**
 * Class to store the parameters of a ternary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  class MemFunData3 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) (Arg1, Arg2, Arg3);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData3 (FunPtr fun_ptr,
		   Class *object,
		   Arg1   arg1,
		   Arg2   arg2,
		   Arg3   arg3);

				       /**
					* Copy constructor.
					*/
      MemFunData3 (const MemFunData3 &fun_data3);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
 
				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object,
					 Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a quaternary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  class MemFunData4 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) (Arg1, Arg2, Arg3, Arg4);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData4 (FunPtr fun_ptr,
		   Class *object,
		   Arg1   arg1,
		   Arg2   arg2,
		   Arg3   arg3,
		   Arg4   arg4);

				       /**
					* Copy constructor.
					*/
      MemFunData4 (const MemFunData4 &fun_data4);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
 
				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      Arg4   arg4;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object,
					 Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3,
					 Arg4 arg4);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };



/**
 * Class to store the parameters of a quintary function. For more
 * information on use and internals of this class see the report on
 * this subject.
 *
 * @author Wolfgang Bangerth, 2000
 */
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  class MemFunData5 : public FunDataBase
  {
    public:
				       /**
					* Typedef a pointer to a global
					* function which we will call
					* from this class.
					*/
      typedef RetType (Class::*FunPtr) (Arg1, Arg2, Arg3, Arg4, Arg5);

				       /**
					* Constructor. Store pointer to
					* the function and the values of
					* the arguments.
					*/
      MemFunData5 (FunPtr fun_ptr,
		   Class *object,
		   Arg1   arg1,
		   Arg2   arg2,
		   Arg3   arg3,
		   Arg4   arg4,
		   Arg5   arg5);

				       /**
					* Copy constructor.
					*/
      MemFunData5 (const MemFunData5 &fun_data5);

				       /**
					* Virtual constructor.
					*/
      virtual FunDataBase * clone () const;

    private:

				       /**
					* Pointer to the function to be
					* called and values of the
					* arguments to be passed.
					*/
      FunPtr fun_ptr;

				       /**
					* Pointer to the object which
					* we shall work on.
					*/
      Class *object;
 
				       /**
					* Values of the arguments of the
					* function to be called.
					*/
      Arg1   arg1;
      Arg2   arg2;
      Arg3   arg3;
      Arg4   arg4;
      Arg5   arg5;
      
				       /**
					* Static function used as entry
					* point for the new thread.
					*/
      static void * thread_entry_point (void *arg);

				       /**
					* Helper class, used to collect
					* the values of the parameters
					* which we will pass to the
					* function, once we know its
					* type.
					*/
      class ArgCollector
      {
	public:
					   /**
					    * Typedef the function
					    * pointer type of the
					    * underlying class to a
					    * local type.
					    */
	  typedef MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::FunPtr FunPtr;
	
					   /**
					    * Constructor. Take and store a
					    * pointer to the function which
					    * is to be called.
					    */
	  ArgCollector (FunPtr fun_ptr);
    
					   /**
					    * Take the arguments with
					    * which we want to call the
					    * function and produce and
					    * object of the desired type
					    * from that.
					    */
	  FunEncapsulation collect_args (Class *object,
					 Arg1 arg1,
					 Arg2 arg2,
					 Arg3 arg3,
					 Arg4 arg4,
					 Arg5 arg5);
    
	private:
					   /**
					    * Space to temporarily store
					    * the function pointer.
					    */
	  FunPtr fun_ptr;
      };
  };
  
  
				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  FunData0<void>::ArgCollector
  encapsulate (void (*fun_ptr)());

				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <typename Arg1>
  typename FunData1<Arg1,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1));


				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <typename Arg1, typename Arg2>
  typename FunData2<Arg1,Arg2,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2));
  

				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <typename Arg1, typename Arg2, typename Arg3>
  typename FunData3<Arg1,Arg2,Arg3,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3));
  

				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  typename FunData4<Arg1,Arg2,Arg3,Arg4,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3, Arg4));



				   /**
				    * Encapsulate a function pointer
				    * into an object with which a new
				    * thread can later be spawned.
				    * For more information on use and
				    * internals of this class see the
				    * report on this subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  typename FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3, Arg4,Arg5));


				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class>
  typename MemFunData0<Class,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)());

				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class, typename Arg1>
  typename MemFunData1<Class,Arg1,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1));


				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class, typename Arg1, typename Arg2>
  typename MemFunData2<Class,Arg1,Arg2,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2));
  

				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class, typename Arg1, typename Arg2, typename Arg3>
  typename MemFunData3<Class,Arg1,Arg2,Arg3,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3));
  

				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  typename MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3, Arg4));
  

				   /**
				    * Encapsulate a member function
				    * pointer into an object with
				    * which a new thread can later be
				    * spawned.  For more information
				    * on use and internals of this
				    * class see the report on this
				    * subject.
				    *
				    * This function exists once for
				    * each number of parameters.
				    */
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  typename MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3, Arg4, Arg5));



				   /**
				    * Spawn a new thread using the
				    * function and parameters
				    * encapsulated in #fun_data#, and
				    * using the given thread manager
				    * object.
				    */
  void spawn (ACE_Thread_Manager     &thread_manager,
	      const FunEncapsulation &fun_data);


				   /**
				    * Spawn several threads at once,
				    * using the same parameters and
				    * calling the same function.
				    */
  void spawn_n (ACE_Thread_Manager     &thread_manager,
		const FunEncapsulation &fun_encapsulation,
		const unsigned int      n_threads);


				   /**
				    * Split the range #[begin,end)#
				    * into #n_intervals# subintervals
				    * of equal size. The last interval
				    * will be a little bit larger, if
				    * the number of elements in the
				    * whole range is not exactly
				    * divisible by #n_intervals#. The
				    * type of the iterators has to
				    * fulfill the requirements of a
				    * forward iterator,
				    * i.e. #operator++# must be
				    * available, and of course it must
				    * be assignable.
				    *
				    * A list of subintervals is
				    * returned as a vector of pairs of
				    * iterators, where each pair
				    * denotes the range
				    * #[begin[i],end[i])#.
				    */
  template <typename ForwardIterator>
  vector<pair<ForwardIterator,ForwardIterator> >
  split_range (const ForwardIterator &begin,
	       const ForwardIterator &end,
	       const unsigned int n_intervals);

				   /**
				    * Split the interval #[begin,end)#
				    * into subintervals of (almost)
				    * equal size. This function works
				    * mostly as the one before, with
				    * the difference that instead of
				    * iterators, now values are taken
				    * that define the whole interval.
				    */
  vector<pair<unsigned int,unsigned int> >
  split_interval (const unsigned int begin,
		  const unsigned int end,
		  const unsigned int n_intervals);
  
  

/**
 * This class is used to make some sanity checks on the numbers of
 * objects of some types related with thread spawning, which are
 * created and deleted. This is a helpful thing when trying to
 * implement the data copying using #clone# functions etc, in order to
 * avoid that there are some objects which are copied but not deleted.
 *
 * It basically only monitors the number of objects which is alive at
 * each time, and complains if the number is nonzero when the counting
 * object is deleted. Since one will probably want to use one global
 * counter, the complaint is raised at the end of the program, and
 * then means that somewhen within the lifetime of your program there
 * has occured a memory leak.
 *
 * This class is not meant for public use.
 *
 * @author Wolfgang Bangerth, 2000
 */
  struct FunDataCounter
  {
				       /**
					* Constructor. Sets all
					* counters to zero.
					*/
      FunDataCounter ();
      
				       /**
					* Destructor. Check whether
					* the total number of objects
					* is zero, otherwise throw an
					* exception.
					*/
      ~FunDataCounter ();

				       /**
					* Counters for the two types
					* of objects which we
					* presently monitor.
					*/
      unsigned int n_fun_encapsulation_objects;
      unsigned int n_fun_data_base_objects;

				       /**
					* Exception
					*/
      DeclException2 (ExcObjectsExist,
		      string, int,
		      << "There are still " << arg2 << " objects of type "
		      << arg1 << " alive. You probably have a memory "
		      << "leak somewhere.");
  };

};   // end declarations of namespace Threads







/* ----------- implementation of functions in namespace Threads ---------- */
namespace Threads 
{

/* ---------------------- FunData0 implementation ------------------------ */

  template <typename RetType>
  FunData0<RetType>::FunData0 (FunPtr fun_ptr) :
		  FunDataBase (&FunData0<RetType>::thread_entry_point),
		  fun_ptr (fun_ptr)
  {};



  template <typename RetType>
  FunData0<RetType>::FunData0 (const FunData0 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr)
  {};



  template <typename RetType>
  FunDataBase *
  FunData0<RetType>::clone () const 
  {
    return new FunData0 (*this);
  };



  template <typename RetType>
  void *
  FunData0<RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData0<RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)();
  
    return 0;
  };



  template <typename RetType>
  FunData0<RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <typename RetType>
  FunEncapsulation
  FunData0<RetType>::ArgCollector::collect_args ()
  {
    return new FunData0<void>(fun_ptr);
  };
 


/* ---------------------- FunData1 implementation ------------------------ */

  template <typename Arg1, typename RetType>
  FunData1<Arg1,RetType>::FunData1 (FunPtr fun_ptr,
				    Arg1   arg1) :
		  FunDataBase (&FunData1<Arg1,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  arg1 (arg1)
  {};



  template <typename Arg1, typename RetType>
  FunData1<Arg1,RetType>::FunData1 (const FunData1 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  arg1 (fun_data.arg1)
  {};



  template <typename Arg1, typename RetType>
  FunDataBase *
  FunData1<Arg1,RetType>::clone () const 
  {
    return new FunData1 (*this);
  };



  template <typename Arg1, typename RetType>
  void *
  FunData1<Arg1,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData1<Arg1,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Arg1              arg1    = fun_data->arg1;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)(arg1);
  
    return 0;
  };



  template <typename Arg1, typename RetType>
  FunData1<Arg1,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  template <typename Arg1, typename RetType>
  FunEncapsulation
  FunData1<Arg1,RetType>::ArgCollector::collect_args (Arg1 arg1)
  {
    return new FunData1<Arg1,void>(fun_ptr, arg1);
  };
 
  


/* ---------------------- FunData2 implementation ------------------------ */

  template <typename Arg1, typename Arg2, typename RetType>
  FunData2<Arg1,Arg2,RetType>::FunData2 (FunPtr fun_ptr,
					 Arg1   arg1,
					 Arg2   arg2) :
		  FunDataBase (&FunData2<Arg1,Arg2,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  arg1 (arg1),
		  arg2 (arg2)
  {};



  template <typename Arg1, typename Arg2, typename RetType>
  FunData2<Arg1,Arg2,RetType>::FunData2 (const FunData2 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2)
  {};



  template <typename Arg1, typename Arg2, typename RetType>
  FunDataBase *
  FunData2<Arg1,Arg2,RetType>::clone () const 
  {
    return new FunData2 (*this);
  };



  template <typename Arg1, typename Arg2, typename RetType>
  void *
  FunData2<Arg1,Arg2,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData2<Arg1,Arg2,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)(arg1, arg2);
  
    return 0;
  };



  template <typename Arg1, typename Arg2, typename RetType>
  FunData2<Arg1,Arg2,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <typename Arg1, typename Arg2, typename RetType>
  FunEncapsulation
  FunData2<Arg1,Arg2,RetType>::ArgCollector::collect_args (Arg1 arg1,
							   Arg2 arg2)
  {
    return new FunData2<Arg1,Arg2,void>(fun_ptr, arg1, arg2);
  };
 
  


/* ---------------------- FunData3 implementation ------------------------ */

  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunData3<Arg1,Arg2,Arg3,RetType>::FunData3 (FunPtr fun_ptr,
					      Arg1   arg1,
					      Arg2   arg2,
					      Arg3   arg3) :
		  FunDataBase (&FunData3<Arg1,Arg2,Arg3,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunData3<Arg1,Arg2,Arg3,RetType>::FunData3 (const FunData3 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunDataBase *
  FunData3<Arg1,Arg2,Arg3,RetType>::clone () const 
  {
    return new FunData3 (*this);
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  void *
  FunData3<Arg1,Arg2,Arg3,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData3<Arg1,Arg2,Arg3,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)(arg1, arg2, arg3);
  
    return 0;
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunData3<Arg1,Arg2,Arg3,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunEncapsulation
  FunData3<Arg1,Arg2,Arg3,RetType>::ArgCollector::collect_args (Arg1 arg1,
								Arg2 arg2,
								Arg3 arg3)
  {
    return new FunData3<Arg1,Arg2,Arg3,void>(fun_ptr, arg1, arg2, arg3);
  };
 
  


/* ---------------------- FunData4 implementation ------------------------ */

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::FunData4 (FunPtr fun_ptr,
						   Arg1   arg1,
						   Arg2   arg2,
						   Arg3   arg3,
						   Arg4   arg4) :
		  FunDataBase (&FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3),
		  arg4 (arg4)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::FunData4 (const FunData4 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3),
		  arg4 (fun_data.arg4)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunDataBase *
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::clone () const 
  {
    return new FunData4 (*this);
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  void *
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData4<Arg1,Arg2,Arg3,Arg4,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;
    Arg4              arg4    = fun_data->arg4;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)(arg1, arg2, arg3, arg4);
  
    return 0;
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunEncapsulation
  FunData4<Arg1,Arg2,Arg3,Arg4,RetType>::ArgCollector::collect_args (Arg1 arg1,
								     Arg2 arg2,
								     Arg3 arg3,
								     Arg4 arg4)
  {
    return new FunData4<Arg1,Arg2,Arg3,Arg4,void>(fun_ptr, arg1, arg2, arg3, arg4);
  };
 
  


/* ---------------------- FunData5 implementation ------------------------ */

  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::FunData5 (FunPtr fun_ptr,
							Arg1   arg1,
							Arg2   arg2,
							Arg3   arg3,
							Arg4   arg4,
							Arg5   arg5) :
		  FunDataBase (&FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3),
		  arg4 (arg4),
		  arg5 (arg5)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::FunData5 (const FunData5 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3),
		  arg4 (fun_data.arg4),
		  arg5 (fun_data.arg5)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunDataBase *
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::clone () const 
  {
    return new FunData5 (*this);
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  void *
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;
    Arg4              arg4    = fun_data->arg4;
    Arg5              arg5    = fun_data->arg5;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (*fun_ptr)(arg1, arg2, arg3, arg4, arg5);
  
    return 0;
  };



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};



  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunEncapsulation
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::ArgCollector::collect_args (Arg1 arg1,
									  Arg2 arg2,
									  Arg3 arg3,
									  Arg4 arg4,
									  Arg5 arg5)
  {
    return new FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,void>(fun_ptr, arg1, arg2,
						       arg3, arg4, arg5);
  };
 
  



/* ---------------------- MemFunData0 implementation ------------------------ */

  template <class Class, typename RetType>
  MemFunData0<Class,RetType>::MemFunData0 (FunPtr  fun_ptr,
					   Class  *object) :
		  FunDataBase (&MemFunData0<Class,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object)
  {};



  template <class Class, typename RetType>
  MemFunData0<Class,RetType>::MemFunData0 (const MemFunData0 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object)
  {};



  template <class Class, typename RetType>
  FunDataBase *
  MemFunData0<Class,RetType>::clone () const 
  {
    return new MemFunData0 (*this);
  };



  template <class Class, typename RetType>
  void *
  MemFunData0<Class,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData0<Class,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;

				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)();
  
    return 0;
  };



  template <class Class, typename RetType>
  MemFunData0<Class,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <class Class, typename RetType>
  FunEncapsulation
  MemFunData0<Class,RetType>::ArgCollector::collect_args (Class *object)
  {
    return new MemFunData0<Class,void>(fun_ptr, object);
  };
 


/* ---------------------- MemFunData1 implementation ------------------------ */

  template <class Class, typename Arg1, typename RetType>
  MemFunData1<Class,Arg1,RetType>::MemFunData1 (FunPtr  fun_ptr,
						Class  *object,
						Arg1    arg1) :
		  FunDataBase (&MemFunData1<Class,Arg1,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object),
		  arg1 (arg1)
  {};



  template <class Class, typename Arg1, typename RetType>
  MemFunData1<Class,Arg1,RetType>::MemFunData1 (const MemFunData1 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object),
		  arg1 (fun_data.arg1)
  {};



  template <class Class, typename Arg1, typename RetType>
  FunDataBase *
  MemFunData1<Class,Arg1,RetType>::clone () const 
  {
    return new MemFunData1 (*this);
  };



  template <class Class, typename Arg1, typename RetType>
  void *
  MemFunData1<Class,Arg1,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData1<Class,Arg1,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;
    Arg1              arg1    = fun_data->arg1;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)(arg1);
  
    return 0;
  };



  template <class Class, typename Arg1, typename RetType>
  MemFunData1<Class,Arg1,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  template <class Class, typename Arg1, typename RetType>
  FunEncapsulation
  MemFunData1<Class,Arg1,RetType>::ArgCollector::collect_args (Class *object,
							       Arg1   arg1)
  {
    return new MemFunData1<Class,Arg1,void>(fun_ptr, object, arg1);
  };
 
  


/* ---------------------- MemFunData2 implementation ------------------------ */

  template <class Class, typename Arg1, typename Arg2, typename RetType>
  MemFunData2<Class,Arg1,Arg2,RetType>::MemFunData2 (FunPtr  fun_ptr,
						     Class  *object,
						     Arg1    arg1,
						     Arg2    arg2) :
		  FunDataBase (&MemFunData2<Class,Arg1,Arg2,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object),
		  arg1 (arg1),
		  arg2 (arg2)
  {};



  template <class Class, typename Arg1, typename Arg2, typename RetType>
  MemFunData2<Class,Arg1,Arg2,RetType>::MemFunData2 (const MemFunData2 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2)
  {};



  template <class Class, typename Arg1, typename Arg2, typename RetType>
  FunDataBase *
  MemFunData2<Class,Arg1,Arg2,RetType>::clone () const 
  {
    return new MemFunData2 (*this);
  };



  template <class Class, typename Arg1, typename Arg2, typename RetType>
  void *
  MemFunData2<Class,Arg1,Arg2,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData2<Class,Arg1,Arg2,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)(arg1, arg2);
  
    return 0;
  };



  template <class Class, typename Arg1, typename Arg2, typename RetType>
  MemFunData2<Class,Arg1,Arg2,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <class Class, typename Arg1, typename Arg2, typename RetType>
  FunEncapsulation
  MemFunData2<Class,Arg1,Arg2,RetType>::ArgCollector::collect_args (Class *object,
								    Arg1   arg1,
								    Arg2   arg2)
  {
    return new MemFunData2<Class,Arg1,Arg2,void>(fun_ptr, object, arg1, arg2);
  };
 
  


/* ---------------------- MemFunData3 implementation ------------------------ */

  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::MemFunData3 (FunPtr  fun_ptr,
							  Class  *object,
							  Arg1    arg1,
							  Arg2    arg2,
							  Arg3    arg3) :
		  FunDataBase (&MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::MemFunData3 (const MemFunData3 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunDataBase *
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::clone () const 
  {
    return new MemFunData3 (*this);
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  void *
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData3<Class,Arg1,Arg2,Arg3,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)(arg1, arg2, arg3);
  
    return 0;
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};
 	

  
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename RetType>
  FunEncapsulation
  MemFunData3<Class,Arg1,Arg2,Arg3,RetType>::ArgCollector::collect_args (Class *object,
									 Arg1   arg1,
									 Arg2   arg2,
									 Arg3   arg3)
  {
    return new MemFunData3<Class,Arg1,Arg2,Arg3,void>(fun_ptr, object,
						      arg1, arg2, arg3);
  };
 
  


/* ---------------------- MemFunData4 implementation ------------------------ */

  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::MemFunData4 (FunPtr  fun_ptr,
							       Class  *object,
							       Arg1    arg1,
							       Arg2    arg2,
							       Arg3    arg3,
							       Arg4    arg4) :
		  FunDataBase (&MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3),
		  arg4 (arg4)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::MemFunData4 (const MemFunData4 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3),
		  arg4 (fun_data.arg4)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunDataBase *
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::clone () const 
  {
    return new MemFunData4 (*this);
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  void *
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;
    Arg4              arg4    = fun_data->arg4;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)(arg1, arg2, arg3, arg4);
  
    return 0;
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename RetType>
  FunEncapsulation
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,RetType>::ArgCollector::collect_args (Class *object,
									      Arg1   arg1,
									      Arg2   arg2,
									      Arg3   arg3,
									      Arg4   arg4)
  {
    return new MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,void>(fun_ptr, object,
							   arg1, arg2, arg3, arg4);
  };
 
  


/* ---------------------- MemFunData5 implementation ------------------------ */

  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::MemFunData5 (FunPtr  fun_ptr,
								    Class  *object,
								    Arg1    arg1,
								    Arg2    arg2,
								    Arg3    arg3,
								    Arg4    arg4,
								    Arg5    arg5) :
		  FunDataBase (&MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::thread_entry_point),
		  fun_ptr (fun_ptr),
		  object (object),
		  arg1 (arg1),
		  arg2 (arg2),
		  arg3 (arg3),
		  arg4 (arg4),
		  arg5 (arg5)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::MemFunData5 (const MemFunData5 &fun_data) :
		  FunDataBase (fun_data),
		  fun_ptr (fun_data.fun_ptr),
		  object (fun_data.object),
		  arg1 (fun_data.arg1),
		  arg2 (fun_data.arg2),
		  arg3 (fun_data.arg3),
		  arg4 (fun_data.arg4),
		  arg5 (fun_data.arg5)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunDataBase *
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::clone () const 
  {
    return new MemFunData5 (*this);
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  void *
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::thread_entry_point (void *arg) 
  {
				     // convenience typedef, since we
				     // will need that class name
				     // several times below
    typedef MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType> ThisClass;
  
    FunEncapsulation *fun_encapsulation
      = reinterpret_cast<FunEncapsulation*>(arg);
    const ThisClass *fun_data
      = dynamic_cast<const ThisClass*> (fun_encapsulation->fun_data_base);

				   // copy the parameters
    ThisClass::FunPtr fun_ptr = fun_data->fun_ptr;
    Class            *object  = fun_data->object;
    Arg1              arg1    = fun_data->arg1;
    Arg2              arg2    = fun_data->arg2;
    Arg3              arg3    = fun_data->arg3;
    Arg4              arg4    = fun_data->arg4;
    Arg5              arg5    = fun_data->arg5;


				   // copying of parameters is done,
				   // now we can release the lock on
				   // #fun_data#
    fun_data->lock.release ();

				     // call the function
    (object->*fun_ptr)(arg1, arg2, arg3, arg4, arg5);
  
    return 0;
  };



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::ArgCollector::ArgCollector (FunPtr fun_ptr) :
		  fun_ptr (fun_ptr)
  {};



  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5, typename RetType>
  FunEncapsulation
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,RetType>::ArgCollector::collect_args (Class *object,
										   Arg1   arg1,
										   Arg2   arg2,
										   Arg3   arg3,
										   Arg4   arg4,
										   Arg5   arg5)
  {
    return new MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,void>(fun_ptr, object,
								arg1, arg2, arg3, arg4, arg5);
  };
 
  


/* ---------------------------------------------------------------- */

  inline
  FunData0<void>::ArgCollector
  encapsulate (void (*fun_ptr)())
  {
    return fun_ptr;
  };



  template <typename Arg1>
  FunData1<Arg1,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1))
  {
    return fun_ptr;
  };
  

  
  template <typename Arg1, typename Arg2>
  FunData2<Arg1,Arg2,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2))
  {
    return fun_ptr;
  };
  

  
  template <typename Arg1, typename Arg2, typename Arg3>
  FunData3<Arg1,Arg2,Arg3,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3))
  {
    return fun_ptr;
  };
  

  
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  FunData4<Arg1,Arg2,Arg3,Arg4,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3, Arg4))
  {
    return fun_ptr;
  };
  
  
  template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  FunData5<Arg1,Arg2,Arg3,Arg4,Arg5,void>::ArgCollector
  encapsulate (void (*fun_ptr)(Arg1, Arg2, Arg3, Arg4, Arg5))
  {
    return fun_ptr;
  };
  
  
  template <class Class>
  MemFunData0<Class,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)())
  {
    return fun_ptr;
  };



  template <class Class, typename Arg1>
  MemFunData1<Class,Arg1,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1))
  {
    return fun_ptr;
  };
  

  
  template <class Class, typename Arg1, typename Arg2>
  MemFunData2<Class,Arg1,Arg2,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2))
  {
    return fun_ptr;
  };
  

  
  template <class Class, typename Arg1, typename Arg2, typename Arg3>
  MemFunData3<Class,Arg1,Arg2,Arg3,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3))
  {
    return fun_ptr;
  };
  

  
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
  MemFunData4<Class,Arg1,Arg2,Arg3,Arg4,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3, Arg4))
  {
    return fun_ptr;
  };


  
  template <class Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
  MemFunData5<Class,Arg1,Arg2,Arg3,Arg4,Arg5,void>::ArgCollector
  encapsulate (void (Class::*fun_ptr)(Arg1, Arg2, Arg3, Arg4, Arg5))
  {
    return fun_ptr;
  };


  
  template <typename ForwardIterator>
  vector<pair<ForwardIterator,ForwardIterator> >
  split_range (const ForwardIterator &begin, const ForwardIterator &end,
	       const unsigned int n_intervals)
  {
    const unsigned int n_elements              = distance (begin, end);
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;
    
    vector<pair<ForwardIterator,ForwardIterator> >
      return_values (n_intervals);

    return_values[0].first = begin;
    for (unsigned int i=0; i<n_intervals; ++i)
      {
	if (i != n_intervals-1) 
	  {
	    return_values[i].second = return_values[i].first;
					     // note: the cast is
					     // performed to avoid a
					     // warning of gcc that in
					     // the library `dist>=0'
					     // is checked (dist has a
					     // template type, which
					     // here is unsigned if no
					     // cast is performed)
	    advance (return_values[i].second,
		     static_cast<signed int>(n_elements_per_interval));
					     // distribute residual in
					     // division equally among
					     // the first few
					     // subintervals
	    if (i < residual)
	      ++return_values[i].second;
	    
	    return_values[i+1].first = return_values[i].second;
	  }
	else
	  return_values[i].second = end;
      };
    return return_values;
  };  


	    
};   // end of implementation of namespace Threads


#endif // DEAL_II_USE_MT


//----------------------------   thread_management.h     ---------------------------
// end of #ifndef __deal2__thread_management_h
#endif
//----------------------------   thread_management.h     ---------------------------
