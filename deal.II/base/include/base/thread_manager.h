/*----------------------------   thread_manager.h     ---------------------------*/
/*      $Id$                 */
#ifndef __thread_manager_H
#define __thread_manager_H
/*----------------------------   thread_manager.h     ---------------------------*/


#ifdef DEAL_II_USE_MT

#include <ace/Thread_Manager.h>


/**
 * This class only wraps up some functionality not included in the
 * #ACE_Thread_Manager# class, namely the possibility to call member
 * functions when spawning threads.
 *
 * @author Wolfgang Bangerth, 1999
 */
class ThreadManager : public ACE_Thread_Manager 
{
  public:
				     /**
				      * This class is used to package all
				      * data needed to call a specific unary
				      * member function of an object. It is
				      * used to pass it over to a static
				      * function which serves as the entry
				      * points for a new thread; that function
				      * in turn calls the member function
				      * with its object and argument.
				      */
    template <typename Class, typename Arg>
    struct Mem_Fun_Data1
    {
					 /**
					  * Convenience #typedef# for the
					  * member functions data type.
					  */
	typedef void * (Class::*MemFun) (Arg);

					 /**
					  * Pointer to the object for which
					  * the member function is to be
					  * called.
					  */
	Class *object;

					 /**
					  * Argument for the function call.
					  */
	Arg    arg;
	
					 /**
					  * Pointer to the member function.
					  */
	MemFun mem_fun;

	Mem_Fun_Data1 (Class *object,
		       Arg    arg,
		       MemFun mem_fun) :
			object (object),
			arg (arg),
			mem_fun (mem_fun) {};
    };

    template <typename Class, typename Arg1, typename Arg2>
    struct Mem_Fun_Data2
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	MemFun mem_fun;

	Mem_Fun_Data2 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			mem_fun (mem_fun) {};
    };


    template <typename Class, typename Arg1, typename Arg2, typename Arg3>
    struct Mem_Fun_Data3
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	MemFun mem_fun;

	Mem_Fun_Data3 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			mem_fun (mem_fun) {};
    };
    

    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4>
    struct Mem_Fun_Data4
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	MemFun mem_fun;

	Mem_Fun_Data4 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       Arg4   arg4,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			mem_fun (mem_fun) {};
    };

    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5>
    struct Mem_Fun_Data5
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	MemFun mem_fun;

	Mem_Fun_Data5 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       Arg4   arg4,
		       Arg5   arg5,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			arg5 (arg5),
			mem_fun (mem_fun) {};
    };

    
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6>
    struct Mem_Fun_Data6
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5, Arg6);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	MemFun mem_fun;

	Mem_Fun_Data6 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       Arg4   arg4,
		       Arg5   arg5,
		       Arg6   arg6,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			arg5 (arg5),
			arg6 (arg6),
			mem_fun (mem_fun) {};
    };
    


				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking a single argument.
				      */
    template <typename ObjectClass, typename Arg>
    int spawn (Mem_Fun_Data1<ObjectClass,Arg> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking two arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2>
    int spawn (Mem_Fun_Data2<ObjectClass,Arg1,Arg2> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking three arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3>
    int spawn (Mem_Fun_Data3<ObjectClass,Arg1,Arg2,Arg3> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking four arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    int spawn (Mem_Fun_Data4<ObjectClass,Arg1,Arg2,Arg3,Arg4> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking five arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5>
    int spawn (Mem_Fun_Data5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking six arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6>
    int spawn (Mem_Fun_Data6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);
    
				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking a single argument.
				      */
    template <typename ObjectClass, typename Arg>
    int spawn_n (size_t n,
		 Mem_Fun_Data1<ObjectClass,Arg> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);

				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking two arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2>
    int spawn_n (size_t n,
		 Mem_Fun_Data2<ObjectClass,Arg1,Arg2> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);

				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking three arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3>
    int spawn_n (size_t n,
		 Mem_Fun_Data3<ObjectClass,Arg1,Arg2,Arg3> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);

				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking four arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    int spawn_n (size_t n,
		 Mem_Fun_Data4<ObjectClass,Arg1,Arg2,Arg3,Arg4> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);

				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking five arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5>
    int spawn_n (size_t n,
		 Mem_Fun_Data5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);

				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking six arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6>
    int spawn_n (size_t n,
		 Mem_Fun_Data6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *mem_fun_data,
		 long flags = THR_NEW_LWP | THR_JOINABLE,
		 long priority = ACE_DEFAULT_THREAD_PRIORITY,
		 int grp_id = -1,
		 ACE_Task_Base *task = 0,
		 ACE_hthread_t thread_handles[] = 0,
		 void *stack[] = 0,
		 size_t stack_size[] = 0);
    
  private:


				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * unary member function.
				      */
    template <typename Class, typename Arg>
    static void * thread_entry_point1 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * binary member function.
				      */
    template <typename Class, typename Arg1, typename Arg2>
    static void * thread_entry_point2 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * ternary member function.
				      */
    template <typename Class, typename Arg1, typename Arg2, typename Arg3>
    static void * thread_entry_point3 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with four parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    static void * thread_entry_point4 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with five parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5>
    static void * thread_entry_point5 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with six parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5, typename Arg6>
    static void * thread_entry_point6 (void *_arg);
};



/* ------------------------------ Template functions -------------------------------- */


template <typename ObjectClass, typename Arg>
int ThreadManager::spawn (Mem_Fun_Data1<ObjectClass,Arg> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point1<ObjectClass,Arg>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2>
int ThreadManager::spawn (Mem_Fun_Data2<ObjectClass,Arg1,Arg2> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point2<ObjectClass,Arg1,Arg2>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3>
int ThreadManager::spawn (Mem_Fun_Data3<ObjectClass,Arg1,Arg2,Arg3> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point3<ObjectClass,Arg1,Arg2,Arg3>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
int ThreadManager::spawn (Mem_Fun_Data4<ObjectClass,Arg1,Arg2,Arg3,Arg4> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point4<ObjectClass,Arg1,Arg2,Arg3,Arg4>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5>
int ThreadManager::spawn (Mem_Fun_Data5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6>
int ThreadManager::spawn (Mem_Fun_Data6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&thread_entry_point6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};



template <typename ObjectClass, typename Arg>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data1<ObjectClass,Arg> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point1<ObjectClass,Arg>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data2<ObjectClass,Arg1,Arg2> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point2<ObjectClass,Arg1,Arg2>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data3<ObjectClass,Arg1,Arg2,Arg3> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point3<ObjectClass,Arg1,Arg2,Arg3>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data4<ObjectClass,Arg1,Arg2,Arg3,Arg4> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point4<ObjectClass,Arg1,Arg2,Arg3,Arg4>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename ObjectClass, typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &thread_entry_point6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};



template <typename Class, typename Arg>
void * ThreadManager::thread_entry_point1 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data1<Class,Arg> *arg = reinterpret_cast<Mem_Fun_Data1<Class,Arg> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg);
};



template <typename Class, typename Arg1, typename Arg2>
void * ThreadManager::thread_entry_point2 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data2<Class,Arg1,Arg2> *arg
    = reinterpret_cast<Mem_Fun_Data2<Class,Arg1,Arg2> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1, arg->arg2);
};



template <typename Class, typename Arg1, typename Arg2, typename Arg3>
void * ThreadManager::thread_entry_point3 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data3<Class,Arg1,Arg2,Arg3> *arg
    = reinterpret_cast<Mem_Fun_Data3<Class,Arg1,Arg2,Arg3> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3);
};



template <typename Class, typename Arg1, typename Arg2, typename Arg3, typename Arg4>
void * ThreadManager::thread_entry_point4 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data4<Class,Arg1,Arg2,Arg3,Arg4> *arg
    = reinterpret_cast<Mem_Fun_Data4<Class,Arg1,Arg2,Arg3,Arg4> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4);
};



template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5>
void * ThreadManager::thread_entry_point5 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data5<Class,Arg1,Arg2,Arg3,Arg4,Arg5> *arg
    = reinterpret_cast<Mem_Fun_Data5<Class,Arg1,Arg2,Arg3,Arg4,Arg5> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5);
};



template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6>
void * ThreadManager::thread_entry_point6 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data6<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *arg
    = reinterpret_cast<Mem_Fun_Data6<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5,
					arg->arg6);
};


#endif

/*----------------------------   thread_manager.h     ---------------------------*/
/* end of #ifndef __thread_manager_H */
#endif
/*----------------------------   thread_manager.h     ---------------------------*/
