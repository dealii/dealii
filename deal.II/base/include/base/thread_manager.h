//----------------------------  thread_manager.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  thread_manager.h  ---------------------------
#ifndef __deal2__thread_manager_h
#define __deal2__thread_manager_h


#ifdef DEAL_II_USE_MT

#include <ace/Thread_Manager.h>


/**
 * This class only wraps up some functionality not included in the
 * #ACE_Thread_Manager# class, namely the possibility to call member
 * functions when spawning threads. Fuethermore, it provides ways to
 * pass arguments to functions taking more than a single parameter.
 *
 * The basic problem with member functions is that they implicitely
 * want to get a pointer to the object they are to work on. Therefore,
 * a member function taking one parameter is in fact a binary
 * function, since the #this# pointer is passed as well. This is
 * reflected by the fact that member function pointers need to be
 * called like this: #object_pointer->(*mem_fun_ptr) (args)#.
 *
 * On the other hand, the thread creation routines (in both the
 * operating system and ACE) only accept unary, global functions as
 * thread entry points. Their argument is a #void*#, which might be
 * used to pass a pointer to a structure to the function; that may
 * then include values of parameters, the address of the member
 * function to be called, and finally the address upon which the
 * function is to operate (i.e. what will be the function's #this#
 * pointer).
 *
 * In practice, this usually leads to code like the following:
 * \begin{verbatim}
 *   class TestClass {
 *            void   spawn_thread ();
 *     static void * thread_entry (void *);
 *            void * threaded_function (int arg1, double arg2);
 *  };
 *
 *  struct ParametersForThreadEntry {
 *    int    arg1;
 *    double arg2;
 *
 *    TestClass * object;
 *  };
 *
 *  void TestClass:spawn_thread () {
 *    ParametersForThreadEntry params;
 *    params.object = this;
 *    params.arg1   = 1;
 *    params.arg2   = 3.1415926;
 *
 *    spawn (..., &thread_entry, (void*)&params, ...);
 *  };
 *
 *  void * TestClass::thread_entry (void * params_tmp) {
 *                        // cast parameters object to the
 *                        // type it really has
 *    ParametersForThreadEntry * params =
 *               (ParametersForThreadEntry*)params_tmp;
 *                        // unpack entries of params
 *                        // and call member function
 *    void * return_value =
 *         (params->object)->threaded_function (params.arg1, params.arg2);
 *    return return_value;
 *  };
 *
 *  void * TestClass::threaded_function (int arg1, double arg2) {
 *    ... // do something in parallel to the main thread
 *  };
 * \end{verbatim}
 *
 * Note that the #static# in the declaration of #thread_entry# means
 * that is not a true member function in that it does not take a
 * #this# pointer and therefore is truely a unary function. However,
 * due to this, it can't access member variables.
 *
 * This program above suffers from several problems:
 * \begin{enumerate}
 *   \item One has to have a different structure for passing arguments
 *     to functions (if there would be another function taking two
 *     integers instead of the argument list above, one would have to have
 *     another parameter structure including two integers).
 *
 *   \item One would also need a different #thread_entry# function for
 *     all functions to be called, since each of them can call only
 *     one member function (one could include the address of the member
 *     function into the structure wrapped around the parameter values;
 *     then one would only need one #thread_entry# function for each
 *     different parameter list of functions).
 *
 * \item The program has a bug: #params# in #spawn_thread# is a local
 *     variable, which is destroyed at the end of the function. If now
 *     the operating system returns after executing the #spawn# call,
 *     but takes some time for internal actions before passing control
 *     to the #thread_entry# function, then the pointer given to
 *     #thread_entry# might not point to anything reasonable any more,
 *     thus producing random results.
 *
 *     The problem gets even worse, if the parameters in
 *     #ParametersForThreadEntry# contain not only values, but also
 *     references of pointers to other data elements; one will then
 *     have to guarantee that not only the #params# variable is
 *     preserved until its values are copied by #thread_entry# to the
 *     places where function parameters of #threaded_function# are
 *     stored (i.e. usually on the stack), but also the lifetime of
 *     the variables to which the elements of #params# point, have to
 *     be long enough to guarantee that they still exist whenever
 *     #threaded_function# uses them.
 * \end{enumerate}
 *
 * The present class at leasts solves the first two problems:
 * \begin{enumerate}
 * \item By providing standardized templated class for parameter
 *     transfer, it is no more necessary to declare parameter
 *     structures one-by-one. For example, #ParametersForThreadEntry#
 *     above could be replaced by
 *     #ThreadManager::Mem_Fun_Data2<TestClass,int,double>#, with
 *     template parameters representing the type of the class for
 *     which a member function shall be called, and the parameters to
 *     this function.
 *
 * \item The #thread_entry# function above was declared as a static
 *     member function in order to not clutter up the global namespace
 *     with thread entry functions for each and every purpose. It
 *     could, however, also be a global function, or a static member
 *     function of another class. This class provides a thread entry
 *     function for each possible list of parameters for a member
 *     function; these thread entry functions obviously are made
 *     #static#, i.e. you need not create an object of class
 *     #ThreadManager# to call them, and they satisfy the requirement
 *     of #spawn# to be of data type (function pointer type) #void *
 *     (*) (void *)#.
 * \end{enumerate}
 *
 * The third problem mentioned above cannot be solved so easily,
 * however. It is still up to the calling function to guarantee that
 * the #params# structure exists long enough and that objects to which
 * elements of #params# point exist long enough. There are at least
 * two strategies for this problem:
 * \begin{enumerate}
 * \item Allocate the parameters object on the stack: one could modify
 *     the example as follows:
 *     \begin{verbatim}
 *       void TestClass:spawn_thread () {
 *         ParametersForThreadEntry *params = new ParametersForThreadEntry;
 *         params->object = this;
 *         params->arg1   = 1;
 *         params->arg2   = 3.1415926;
 *
 *         spawn (..., &thread_entry, (void*)params, ...);
 *       };
 *     \end{verbatim}
 *
 *     Thus, the parameters object is on the heap instead of on the
 *     stack, and its lifetime is until it is #delete#d some
 *     time. Again, the #spawn_thread# function can't do that since it
 *     does not exactly know at which point the second thread does not
 *     need the data any more. However, the #thread_entry# function
 *     could do that:
 *     \begin{verbatim}
 *       void * TestClass::thread_entry (void * params_tmp) {
 *                             // cast parameters object to the
 *                             // type it really has
 *         ParametersForThreadEntry * params =
 *                    (ParametersForThreadEntry*)params_tmp;
 *                             // unpack entries of params
 *                             // and call member function
 *         void * return_value =
 *              (params->object)->threaded_function (params.arg1, params.arg2);
 *
 *                             // delete parameters object
 *         delete params;
 * 
 *         return return_value;
 *       };
 *     \end{verbatim}
 *     This is safe, since the parameters object is deleted only after
 *     the member function #threaded_function# returns; the parameters
 *     are therefore no more needed.
 *
 *     The downside here is that there is another system function
 *     which is commonly called: #spawn_n#, which creates #n# threads
 *     at the same time, i.e. it jumps into #thread_entry# #n# times
 *     at once. However, the #delete# operation must only be performed
 *     once, namely by the thread which exits last; the code in
 *     #thread_entry# would therefore have to synchronize which thread
 *     calls the #delete# and when. This is feasible, but difficult.
 *
 * \item Blocking the first thread after spawning other threads: this would
 *     yield an implementation of #spawn_thread# like this:
 *     \begin{verbatim}
 *       void TestClass:spawn_thread () {
 *         ParametersForThreadEntry params;
 *         params.object = this;
 *         params.arg1   = 1;
 *         params.arg2   = 3.1415926;
 *
 *         spawn (..., &thread_entry, (void*)&params, ...);
 *
 *         ... // some code which waits until the spawned thread returns
 *       };
 *     \end{verbatim}
 *     Since execution of #spawn_treads# is suspended until the spawned
 *     thread exits, so is the destruction of the #params# object. It is
 *     therefore guaranteed that it exists longer than the lifetime of the
 *     thread which might use it.
 *
 *     Obviously, the above function is useless as stated here, since
 *     if we start a new threads and then stop the old one until the
 *     new one returns, we could as well have called the member
 *     function directly, without need to create a new thread. This
 *     approach therefore is only useful, if the function creates more
 *     than one thread and waits for them all to return. Thread
 *     creation may happen using several #spawn# calls (and maybe
 *     different parameter objects), as well as using the #spawn_n#
 *     function. Destruction of the parameter object remains with the
 *     calling function, as in the original example.
 * \end{enumerate}
 *
 *
 * \subsection{Example of use of this class}
 *
 * The following example shows how to use the elements of this class.
 * \begin{verbatim}
 *  void TestClass:spawn_thread () {
 *                     // create ThreadManager object
 *    ThreadManager thread_manager;
 *                     // generate an object to pass
 *                     // the two parameters
 *    const ThreadManager::Mem_Fun_Data2<const SparseMatrix<number>,
 *	                                 unsigned int,
 *	                                 unsigned int> 
 *	    mem_fun_data (this,     // object to operate on
 *                        1,        // first parameter
 *                        2,        // second parameter
 *                                  // address of member function
 *	 	          &TestClass::threaded_function 
 *                       );
 *
 *                     // spawn a thread
 *     thread_manager.spawn (&mem_fun_data);
 *
 *     ...             // do something more, start more threads, etc
 *
 *
 *		       // ... and wait until they're finished:
 *     thread_manager.wait ();
 *  };
 * \end{verbatim}
 *
 * Note that in this example, there is no need for the #thread_entry#
 * function and the structure encapsulating parameters. Often,
 * functions will not have a return value, i.e. their return type is
 * #void#. In this case, the interface above would force us to declare
 * them as returning a #void*# anyway. Therefore, this class also
 * allows for functions without a return type; a Null pointer will
 * then be returned to the caller.
 *
 * Usage is exactly as above, but that the #threaded_function# may
 * also be a #void# function. This is done by a second constructor of
 * the #Mem_Fun_Data# structures, where either a pointer to a #void#
 * function or a pointer to a #void*# function is set. The other
 * pointer is always set to zero. If the function returns no value,
 * then the #thread_entry# function of this class returns a Null
 * pointer instead.
 *
 *
 * \subsection{Using this class on global functions}
 *
 * The structures provided by this class are easier to use than the
 * methods of the basic ThreadManager. It is also possible to call global
 * functions (without a this pointer) as a thread. The struct used to
 * package all data needed to call a global functions is called #Fun_Data#. 
 *
 * The following example shows how to call a global function
 * \begin{verbatim}
 *  void TestClass::spawn_thread () {
 *                     // create ThreadManager object
 *    ThreadManager thread_manager;
 *                     // generate an object to pass
 *                     // the two parameters
 *    const ThreadManager::Fun_Data2<unsigned int,
 *	                             unsigned int> 
 *	    fun_data ( 1,        // first parameter
 *                     2,        // second parameter
 *                               // address of global function
 *	               &threaded_function 
 *                   );
 *
 *                     // spawn a thread
 *     thread_manager.spawn (&fun_data);
 *
 *     ...             // do something more, start more threads, etc
 *
 *
 *		       // ... and wait until they're finished:
 *     thread_manager.wait ();
 *  };
 * \end{verbatim}
 * @author Wolfgang Bangerth, 1999, 2000. Extensions by Thomas Richter, 2000
 */
class ThreadManager : public ACE_Thread_Manager 
{
  public:
				     /**
				      * This class is used to package
				      * all data needed to call a
				      * specific void member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class>
    struct Mem_Fun_Data0
    {
	typedef void * (Class::*MemFun) ();

	typedef void (Class::*VoidMemFun) ();

	Class     *object;
	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

	Mem_Fun_Data0 (Class *object,
		       MemFun mem_fun) :
			object (object),
			mem_fun (mem_fun),
			void_mem_fun (0) {};

	Mem_Fun_Data0 (Class     *object,
		       VoidMemFun void_mem_fun) :
			object (object),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };


/**
				      * This class is used to package all
				      * data needed to call a specific unary
				      * member function of an object. It is
				      * used to pass it over to a static
				      * function which serves as the entry
				      * points for a new thread; that function
				      * in turn calls the member function
				      * with its object and argument.
				      *
				      * By the way, the strange naming
				      * of this class with capital
				      * latters and underscores tries
				      * to mimic the corresponding
				      * classes of the STL for
				      * non-member functions.
				      */
    template <typename Class, typename Arg>
    struct Mem_Fun_Data1
    {
					 /**
					  * Convenience #typedef# for
					  * the member functions data
					  * type. This is for member
					  * functions that return a
					  * #void *#, which is
					  * returned to the caller.
					  */
	typedef void * (Class::*MemFun) (Arg);

					 /**
					  * Convenience #typedef# for
					  * the member functions data
					  * type. This is for member
					  * functions without return
					  * value. A #(void*)0# is
					  * then returned to the
					  * caller.
					  */
	typedef void (Class::*VoidMemFun) (Arg);
	
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
					  * Pointer to the member
					  * function if that returns a
					  * #void*#. If the function
					  * has no return value, then
					  * this pointer is zero.
					  */
	MemFun     mem_fun;

					 /**
					  * If the function has no
					  * return value, then this
					  * pointer is used. If the
					  * function has a return
					  * value, then this pointer
					  * is zero.
					  */
	VoidMemFun void_mem_fun;

					 /**
					  * Constructor for functions
					  * with a return value.
					  */
	Mem_Fun_Data1 (Class *object,
		       Arg    arg,
		       MemFun mem_fun) :
			object (object),
			arg (arg),
			mem_fun (mem_fun),
			void_mem_fun (0) {};

					 /**
					  * Constructor for functions
					  * without a return value.
					  */
	Mem_Fun_Data1 (Class *object,
		       Arg    arg,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg (arg),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };

				     /**
				      * This class is used to package
				      * all data needed to call a
				      * specific binary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class, typename Arg1, typename Arg2>
    struct Mem_Fun_Data2
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2);

	typedef void (Class::*VoidMemFun) (Arg1, Arg2);

	Class     *object;
	Arg1       arg1;
	Arg2       arg2;

	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

	Mem_Fun_Data2 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			mem_fun (mem_fun),
			void_mem_fun (0) {};

	Mem_Fun_Data2 (Class     *object,
		       Arg1       arg1,
		       Arg2       arg2,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * specific ternary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class, typename Arg1, typename Arg2, typename Arg3>
    struct Mem_Fun_Data3
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3);
	typedef void   (Class::*VoidMemFun) (Arg1, Arg2, Arg3);

	Class     *object;
	Arg1       arg1;
	Arg2       arg2;
	Arg3       arg3;
	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

	Mem_Fun_Data3 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       MemFun mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			mem_fun (mem_fun),
			void_mem_fun (0) {};

	Mem_Fun_Data3 (Class     *object,
		       Arg1       arg1,
		       Arg2       arg2,
		       Arg3       arg3,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * specific quaternary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4>
    struct Mem_Fun_Data4
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4);
	typedef void   (Class::*VoidMemFun) (Arg1, Arg2, Arg3, Arg4);
	
	Class     *object;
	Arg1       arg1;
	Arg2       arg2;
	Arg3       arg3;
	Arg4       arg4;
	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

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
			mem_fun (mem_fun),
			void_mem_fun (0) {};
	
	Mem_Fun_Data4 (Class     *object,
		       Arg1       arg1,
		       Arg2       arg2,
		       Arg3       arg3,
		       Arg4       arg4,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };

				     /**
				      * This class is used to package
				      * all data needed to call a
				      * specific quintary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5>
    struct Mem_Fun_Data5
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5);
	typedef void   (Class::*VoidMemFun) (Arg1, Arg2, Arg3, Arg4, Arg5);
	Class     *object;
	Arg1       arg1;
	Arg2       arg2;
	Arg3       arg3;
	Arg4       arg4;
	Arg5       arg5;
	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

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
			mem_fun (mem_fun),
			void_mem_fun (0) {};

	Mem_Fun_Data5 (Class     *object,
		       Arg1       arg1,
		       Arg2       arg2,
		       Arg3       arg3,
		       Arg4       arg4,
		       Arg5       arg5,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			arg5 (arg5),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * specific hexary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6>
    struct Mem_Fun_Data6
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5, Arg6);
	typedef void   (Class::*VoidMemFun) (Arg1, Arg2, Arg3, Arg4, Arg5, Arg6);

	Class     *object;
	Arg1       arg1;
	Arg2       arg2;
	Arg3       arg3;
	Arg4       arg4;
	Arg5       arg5;
	Arg6       arg6;
	MemFun     mem_fun;
	VoidMemFun void_mem_fun;

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
			mem_fun (mem_fun),
			void_mem_fun (0) {};

	Mem_Fun_Data6 (Class     *object,
		       Arg1       arg1,
		       Arg2       arg2,
		       Arg3       arg3,
		       Arg4       arg4,
		       Arg5       arg5,
		       Arg6       arg6,
		       VoidMemFun void_mem_fun) :
			object (object),
			arg1 (arg1),
			arg2 (arg2),
			arg3 (arg3),
			arg4 (arg4),
			arg5 (arg5),
			arg6 (arg6),
			mem_fun (0),
			void_mem_fun (void_mem_fun) {};
    };
    
				     /**
				      * This class is used to package
				      * all data needed to call a
				      * specific heptary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7>
      struct Mem_Fun_Data7
      {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5,
					 Arg6, Arg7);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	Arg7   arg7;
	MemFun mem_fun;

	Mem_Fun_Data7 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       Arg4   arg4,
		       Arg5   arg5,
		       Arg6   arg6,
		       Arg7   arg7,
		       MemFun mem_fun) :
			object (object),
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			arg6  (arg6),
                        arg7  (arg7),
			mem_fun (mem_fun) {};
    };
    
				     /**
				      * This class is used to package
				      * all data needed to call a
				      * specific octary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8>
    struct Mem_Fun_Data8
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5,
					 Arg6, Arg7, Arg8);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	Arg7   arg7;
	Arg8   arg8;
	MemFun mem_fun;

	Mem_Fun_Data8 (Class *object,
		       Arg1   arg1,
		       Arg2   arg2,
		       Arg3   arg3,
		       Arg4   arg4,
		       Arg5   arg5,
		       Arg6   arg6,
		       Arg7   arg7,
		       Arg8   arg8,
		       MemFun mem_fun) :
			object (object),
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			arg6  (arg6),
                        arg7  (arg7),
			arg8  (arg8),
                        mem_fun (mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * specific nonary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8,
              typename Arg9>
    struct Mem_Fun_Data9
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5,
					 Arg6, Arg7, Arg8, Arg9);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	Arg7   arg7;
	Arg8   arg8;
	Arg9   arg9;
	MemFun mem_fun;

	Mem_Fun_Data9 (Class *object,
			Arg1   arg1,
			Arg2   arg2,
			Arg3   arg3,
			Arg4   arg4,
			Arg5   arg5,
			Arg6   arg6,
			Arg7   arg7,
			Arg8   arg8,
			Arg9   arg9,
		       MemFun mem_fun) :
			object (object),
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			arg6  (arg6),
                        arg7  (arg7),
			arg8  (arg8),
			arg9  (arg9),
			mem_fun (mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * specific decary member
				      * function of an object. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Class,
              typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8,
              typename Arg9, typename Arg10>
    struct Mem_Fun_Data10
    {
	typedef void * (Class::*MemFun) (Arg1, Arg2, Arg3, Arg4, Arg5,
					 Arg6, Arg7, Arg8, Arg9, Arg10);
	Class *object;
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	Arg7   arg7;
	Arg8   arg8;
	Arg9   arg9;
	Arg10  arg10;
	MemFun mem_fun;

	Mem_Fun_Data10 (Class *object,
			Arg1   arg1,
			Arg2   arg2,
			Arg3   arg3,
			Arg4   arg4,
			Arg5   arg5,
			Arg6   arg6,
			Arg7   arg7,
			Arg8   arg8,
			Arg9   arg9,
			Arg10  arg10,
			MemFun mem_fun) :
			object (object),
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			arg6  (arg6),
                        arg7  (arg7),
			arg8  (arg8),
			arg9  (arg9),
			arg10 (arg10),
			mem_fun (mem_fun) {};
    };


/**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1>
      struct Fun_Data1
      {
 	  typedef void * (*FunPtr) (Arg1);
	  typedef void   (*VoidFunPtr) (Arg1);
	  Arg1   arg1;
	  FunPtr fun_ptr;
	  VoidFunPtr void_fun_ptr;
	  
	  Fun_Data1 (Arg1   arg1,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  fun_ptr (fun_ptr),
			  void_fun_ptr (0) {};
	  
	  Fun_Data1 (Arg1   arg1,
		     VoidFunPtr void_fun_ptr) :
			  arg1  (arg1),
			  fun_ptr (0),
			  void_fun_ptr (void_fun_ptr) {};
      };


/**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2>
      struct Fun_Data2
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2);
 	  typedef void   (*VoidFunPtr) (Arg1 , Arg2);
	  Arg1       arg1;
	  Arg2       arg2;
	  FunPtr     fun_ptr;
	  VoidFunPtr void_fun_ptr;
	  
	  Fun_Data2 (Arg1   arg1,
		     Arg2   arg2,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  fun_ptr (fun_ptr),
			  void_fun_ptr (0) {};

	  Fun_Data2 (Arg1       arg1,
		     Arg2       arg2,
		     VoidFunPtr void_fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  fun_ptr (0),
			  void_fun_ptr (void_fun_ptr) {};
      };

    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3>
      struct Fun_Data3
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3);
 	  typedef void   (*VoidFunPtr) (Arg1 , Arg2, Arg3);
	  Arg1   arg1;
	  Arg2   arg2;
	  Arg3   arg3;
	  FunPtr fun_ptr;
	  VoidFunPtr void_fun_ptr;
	  
	  Fun_Data3 (Arg1   arg1,
		     Arg2   arg2,
		     Arg3   arg3,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  fun_ptr (fun_ptr),
			  void_fun_ptr (0) {};

	  Fun_Data3 (Arg1   arg1,
		     Arg2   arg2,
		     Arg3   arg3,
		     VoidFunPtr void_fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  fun_ptr (0),
			  void_fun_ptr (void_fun_ptr) {};
      };
    
    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4>
      struct Fun_Data4
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4);
 	  typedef void   (*VoidFunPtr) (Arg1 , Arg2, Arg3, Arg4);

	  Arg1       arg1;
	  Arg2       arg2;
	  Arg3       arg3;
	  Arg4       arg4;
	  FunPtr     fun_ptr;
	  VoidFunPtr void_fun_ptr;
	  
	  Fun_Data4 (Arg1   arg1,
		     Arg2   arg2,
		     Arg3   arg3,
		     Arg4   arg4,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  arg4  (arg4),
			  fun_ptr (fun_ptr),
			  void_fun_ptr (0) {};
	  Fun_Data4 (Arg1       arg1,
		     Arg2       arg2,
		     Arg3       arg3,
		     Arg4       arg4,
		     VoidFunPtr void_fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  arg4  (arg4),
			  fun_ptr (0),
			  void_fun_ptr (void_fun_ptr) {};
      };

    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5>
    struct Fun_Data5
    {
	typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5);
	typedef void   (*VoidFunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5);

	Arg1       arg1;
	Arg2       arg2;
	Arg3       arg3;
	Arg4       arg4;
	Arg5       arg5;
	FunPtr     fun_ptr;
	VoidFunPtr void_fun_ptr;
	  
	Fun_Data5 (Arg1   arg1,
		   Arg2   arg2,
		   Arg3   arg3,
		   Arg4   arg4,
		   Arg5   arg5,
		   FunPtr fun_ptr) :
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			fun_ptr (fun_ptr),
			void_fun_ptr (0) {};
	Fun_Data5 (Arg1       arg1,
		   Arg2       arg2,
		   Arg3       arg3,
		   Arg4       arg4,
		   Arg5       arg5,
		   VoidFunPtr void_fun_ptr) :
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			fun_ptr (0),
			void_fun_ptr (void_fun_ptr) {};
    };

    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6>
      struct Fun_Data6
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5,
				    Arg6);
	  Arg1   arg1;
	  Arg2   arg2;
	  Arg3   arg3;
	  Arg4   arg4;
	  Arg5   arg5;
	  Arg6   arg6;
	  FunPtr fun_ptr;
	  
	  Fun_Data6 (Arg1   arg1,
		     Arg2   arg2,
		     Arg3   arg3,
		     Arg4   arg4,
		     Arg5   arg5,
		     Arg6   arg6,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  arg4  (arg4),
			  arg5  (arg5),
			  arg6  (arg6),
			  fun_ptr (fun_ptr) {};
      };

    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7>
    struct Fun_Data7
    {
	typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5,
				  Arg6, Arg7);
	Arg1   arg1;
	Arg2   arg2;
	Arg3   arg3;
	Arg4   arg4;
	Arg5   arg5;
	Arg6   arg6;
	Arg7   arg7;
	FunPtr fun_ptr;
	  
	Fun_Data7 (Arg1   arg1,
		   Arg2   arg2,
		   Arg3   arg3,
		   Arg4   arg4,
		   Arg5   arg5,
		   Arg6   arg6,
		   Arg7   arg7,
		   FunPtr fun_ptr) :
			arg1  (arg1),
			arg2  (arg2),
			arg3  (arg3),
			arg4  (arg4),
			arg5  (arg5),
			arg6  (arg6),
			arg7  (arg7),
			fun_ptr (fun_ptr) {};
    };
    
    				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8>
      struct Fun_Data8
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5,
				    Arg6, Arg7, Arg8);
	  Arg1   arg1;
	  Arg2   arg2;
	  Arg3   arg3;
	  Arg4   arg4;
	  Arg5   arg5;
	  Arg6   arg6;
	  Arg7   arg7;
	  Arg8   arg8;
	  FunPtr fun_ptr;
	  
	  Fun_Data8 (Arg1   arg1,
		     Arg2   arg2,
		     Arg3   arg3,
		     Arg4   arg4,
		     Arg5   arg5,
		     Arg6   arg6,
		     Arg7   arg7,
		     Arg8   arg8,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
			  arg2  (arg2),
			  arg3  (arg3),
			  arg4  (arg4),
			  arg5  (arg5),
			  arg6  (arg6),
			  arg7  (arg7),
			  arg8  (arg8),
			  fun_ptr (fun_ptr) {};
      };
    
				     /**
				      * This class is used to package
				      * all data needed to call a
				      * global function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8,
      typename Arg9>
      struct Fun_Data9
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5,
				    Arg6, Arg7, Arg8, Arg9);
	  Arg1   arg1;
	  Arg2   arg2;
	  Arg3   arg3;
	  Arg4   arg4;
	  Arg5   arg5;
	  Arg6   arg6;
	  Arg7   arg7;
	  Arg8   arg8;
	  Arg9   arg9;
	  FunPtr fun_ptr;
	  
	  Fun_Data9 (Arg1   arg1,
		      Arg2   arg2,
		      Arg3   arg3,
		      Arg4   arg4,
		      Arg5   arg5,
		      Arg6   arg6,
		      Arg7   arg7,
		      Arg8   arg8,
		      Arg9   arg9,
		     FunPtr fun_ptr) :
			  arg1  (arg1),
	    arg2  (arg2),
	    arg3  (arg3),
	    arg4  (arg4),
	    arg5  (arg5),
	    arg6  (arg6),
	    arg7  (arg7),
	    arg8  (arg8),
	    arg9  (arg9),
	    fun_ptr (fun_ptr) {};
      };
    				     
				     /**
				      * This class is used to package
				      * all data needed to call a
				      * function. See the
				      * general documentation of the
				      * #ThreadManager# class or of
				      * the class
				      * #ThreadManager::Mem_Fun_Data1#
				      * for more information.
				      */
    template <typename Arg1, typename Arg2,
              typename Arg3, typename Arg4,
              typename Arg5, typename Arg6,
              typename Arg7, typename Arg8,
      typename Arg9, typename Arg10>
      struct Fun_Data10
      {
 	  typedef void * (*FunPtr) (Arg1 , Arg2, Arg3, Arg4, Arg5,
				    Arg6, Arg7, Arg8, Arg9, Arg10);
	  Arg1   arg1;
	  Arg2   arg2;
	  Arg3   arg3;
	  Arg4   arg4;
	  Arg5   arg5;
	  Arg6   arg6;
	  Arg7   arg7;
	  Arg8   arg8;
	  Arg9   arg9;
	  Arg10  arg10;
	  FunPtr fun_ptr;
	  
	  Fun_Data10 (Arg1   arg1,
		      Arg2   arg2,
		      Arg3   arg3,
		      Arg4   arg4,
		      Arg5   arg5,
		      Arg6   arg6,
		      Arg7   arg7,
		      Arg8   arg8,
		      Arg9   arg9,
		      Arg10  arg10,
		      FunPtr fun_ptr) :
			  arg1  (arg1),
	    arg2  (arg2),
	    arg3  (arg3),
	    arg4  (arg4),
	    arg5  (arg5),
	    arg6  (arg6),
	    arg7  (arg7),
	    arg8  (arg8),
	    arg9  (arg9),
	    arg10 (arg10),
	    fun_ptr (fun_ptr) {};
      };


/**
				      * Wrapper function to allow spawning
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking no arguments.
				      */
    template <typename ObjectClass>
    int spawn (Mem_Fun_Data0<ObjectClass> *mem_fun_data,
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
				      * threads for member functions as well,
				      * rather than for global functions only.
				      *
				      * This version is for member functions
				      * taking seven arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7>
    int spawn (Mem_Fun_Data7<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
	                                 Arg6,Arg7> *mem_fun_data,
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
				      * taking eight arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8>
      int spawn (Mem_Fun_Data8<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
	                                   Arg6,Arg7,Arg8> *mem_fun_data,
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
				      * taking nine arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9>
    int spawn (Mem_Fun_Data9<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
	                                  Arg6,Arg7,Arg8,Arg9> *mem_fun_data,
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
				      * taking ten arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9, typename Arg10>
    int spawn (Mem_Fun_Data10<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
	                                  Arg6,Arg7,Arg8,Arg9,Arg10> *mem_fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);


/**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking one arguments
				      */
    template <typename Arg1>
    int spawn (Fun_Data1<Arg1> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking two arguments
				      */
    template <typename Arg1, typename Arg2>
    int spawn (Fun_Data2<Arg1,Arg2> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking three arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3>
    int spawn (Fun_Data3<Arg1,Arg2,Arg3> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking four arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    int spawn (Fun_Data4<Arg1,Arg2,Arg3,Arg4> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking five arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5>
    int spawn (Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking six arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5, typename Arg6>
    int spawn (Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,
	       Arg6> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking seven arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7>
    int spawn (Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,
	       Arg6,Arg7> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking eigth arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8>
    int spawn (Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,
	       Arg6,Arg7,Arg8> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

    				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking nine arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8,
      typename Arg9>
    int spawn (Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,
	       Arg6,Arg7,Arg8,Arg9> *fun_data,
	       long flags = THR_NEW_LWP | THR_JOINABLE,
	       ACE_thread_t * = 0,
	       ACE_hthread_t *t_handle = 0,
	       long priority = ACE_DEFAULT_THREAD_PRIORITY,
	       int grp_id = -1,
	       void *stack = 0,
	       size_t stack_size = 0);

				     /**
				      * Wrapper function to allow spawning
				      * threads for funtions.
				      *
				      * This version is for member functions
				      * taking ten arguments
				      */
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8,
      typename Arg9, typename Arg10>
    int spawn (Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,
	       Arg6,Arg7,Arg8,Arg9,Arg10> *fun_data,
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
    template <typename ObjectClass>
    int spawn_n (size_t n,
		 Mem_Fun_Data0<ObjectClass> *mem_fun_data,
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
    
				     /**
				      * Wrapper function to allow spawning
				      * multiple threads for member functions
				      * as well, rather than for global
				      * functions only.
				      *
				      * This version is for member functions
				      * taking seven arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7>
    int spawn_n (size_t n,
		 Mem_Fun_Data7<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
		                            Arg6,Arg7> *mem_fun_data,
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
				      * taking eight arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8>
    int spawn_n (size_t n,
		 Mem_Fun_Data8<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
		                           Arg6,Arg7,Arg8> *mem_fun_data,
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
				      * taking nine arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9>
    int spawn_n (size_t n,
		 Mem_Fun_Data9<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
		                            Arg6,Arg7,Arg8,Arg9> *mem_fun_data,
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
				      * taking ten arguments
				      */
    template <typename ObjectClass, typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9, typename Arg10>
    int spawn_n (size_t n,
		 Mem_Fun_Data10<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
		                            Arg6,Arg7,Arg8,Arg9,Arg10> *mem_fun_data,
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
				      * taking one arguments
				      */
    template < typename Arg1>
    int spawn_n (size_t n,
		 Fun_Data1<Arg1> *fun_data,
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
    template < typename Arg1, typename Arg2>
    int spawn_n (size_t n,
		 Fun_Data2<Arg1,Arg2> *fun_data,
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
    template < typename Arg1, typename Arg2,
                                    typename Arg3>
    int spawn_n (size_t n,
		 Fun_Data3<Arg1,Arg2,Arg3> *fun_data,
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
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4>
    int spawn_n (size_t n,
		 Fun_Data4<Arg1,Arg2,Arg3,Arg4> *fun_data,
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
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5>
    int spawn_n (size_t n,
		 Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *fun_data,
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
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6>
    int spawn_n (size_t n,
		 Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,
		 Arg6> *fun_data,
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
				      * taking seven arguments
				      */
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7>
    int spawn_n (size_t n,
		 Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,
		 Arg6,Arg7> *fun_data,
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
				      * taking eight arguments
				      */
    template < typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8>
    int spawn_n (size_t n,
		 Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,
		 Arg6,Arg7,Arg8> *fun_data,
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
				      * taking nine arguments
				      */
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9>
    int spawn_n (size_t n,
		 Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,
		 Arg6,Arg7,Arg8,Arg9> *fun_data,
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
				      * taking ten arguments
				      */
    template < typename Arg1, typename Arg2,
                                    typename Arg3, typename Arg4,
                                    typename Arg5, typename Arg6,
                                    typename Arg7, typename Arg8,
                                    typename Arg9, typename Arg10>
    int spawn_n (size_t n,
		 Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,
		 Arg6,Arg7,Arg8,Arg9,Arg10> *fun_data,
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
				      * void member function.
				      */
    template <typename Class>
    static void * thread_entry_point0 (void *_arg);

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

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with seven parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5, typename Arg6,
                              typename Arg7>
    static void * thread_entry_point7 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with eight parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5, typename Arg6,
                              typename Arg7, typename Arg8>
    static void * thread_entry_point8 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with nine parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5, typename Arg6,
                              typename Arg7, typename Arg8,
                              typename Arg9>
    static void * thread_entry_point9 (void *_arg);

				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with ten parameters.
				      */
    template <typename Class, typename Arg1, typename Arg2,
                              typename Arg3, typename Arg4,
                              typename Arg5, typename Arg6,
                              typename Arg7, typename Arg8,
                              typename Arg9, typename Arg10>
    static void * thread_entry_point10 (void *_arg);


/**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with one parameters.
				      */
    template <typename Arg1>
    static void * thread_entry_point_1 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with two parameters.
				      */
    template <typename Arg1, typename Arg2>
    static void * thread_entry_point_2 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with three parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3>
    static void * thread_entry_point_3 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with four parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4>
    static void * thread_entry_point_4 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with five parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5>
    static void * thread_entry_point_5 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with six parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6>
    static void * thread_entry_point_6 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with seven parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7>
    static void * thread_entry_point_7 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with eight parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8>
    static void * thread_entry_point_8 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with nine parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8,
      typename Arg9>
    static void * thread_entry_point_9 (void *_arg);

    				     /**
				      * This is a function satisfying the
				      * requirements for thread entry points.
				      * It takes as argument all the
				      * information necessary to call a
				      * member function with ten parameters.
				      */
    template <typename Arg1, typename Arg2,
      typename Arg3, typename Arg4,
      typename Arg5, typename Arg6,
      typename Arg7, typename Arg8,
      typename Arg9, typename Arg10>
    static void * thread_entry_point_10 (void *_arg);

};


/* ------------------------------ Template functions -------------------------------- */


template <typename ObjectClass>
int ThreadManager::spawn (Mem_Fun_Data0<ObjectClass> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point0<ObjectClass>,
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
int ThreadManager::spawn (Mem_Fun_Data1<ObjectClass,Arg> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point1<ObjectClass,Arg>,
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
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point2<ObjectClass,Arg1,Arg2>,
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
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point3<ObjectClass,Arg1,Arg2,Arg3>,
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
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point4<ObjectClass,Arg1,Arg2,Arg3,Arg4>,
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
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5>,
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
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7>
int ThreadManager::spawn (Mem_Fun_Data7<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                            Arg6,Arg7> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point7<ObjectClass,
				    Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8>
int ThreadManager::spawn (Mem_Fun_Data8<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                            Arg6,Arg7,Arg8> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point8<ObjectClass,
				    Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8>,
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
  typename Arg5, typename Arg6,
  typename Arg7, typename Arg8,
  typename Arg9>
int ThreadManager::spawn (Mem_Fun_Data9<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6,Arg7,Arg8,Arg9> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point9<ObjectClass,
				    Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9, typename Arg10>
int ThreadManager::spawn (Mem_Fun_Data10<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                            Arg6,Arg7,Arg8,Arg9,Arg10> *mem_fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point10<ObjectClass,
				    Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10>,
				    (void*)mem_fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1>
int ThreadManager::spawn (Fun_Data1<Arg1> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_1<
				    Arg1>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2>
int ThreadManager::spawn (Fun_Data2<Arg1,Arg2> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_2
				    <Arg1,Arg2>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3>
int ThreadManager::spawn (Fun_Data3<Arg1,Arg2,Arg3> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_3
				    <Arg1,Arg2,Arg3>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4>
int ThreadManager::spawn (Fun_Data4<Arg1,Arg2,Arg3,Arg4> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_4<
				    Arg1,Arg2,Arg3,Arg4>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5>
int ThreadManager::spawn (Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_5
				    <Arg1,Arg2,Arg3,Arg4,Arg5>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6>
int ThreadManager::spawn (Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_6
				    <Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7>
int ThreadManager::spawn (Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6,Arg7> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_7
				    <Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8>
int ThreadManager::spawn (Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6,Arg7,Arg8> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_8
				    <Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9>
int ThreadManager::spawn (Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6,Arg7,Arg8,Arg9> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_9
				    <Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9, typename Arg10>
int ThreadManager::spawn (Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,
			  Arg6,Arg7,Arg8,Arg9,Arg10> *fun_data,
			  long flags,
			  ACE_thread_t *t,
			  ACE_hthread_t *t_handle,
			  long priority,
			  int grp_id,
			  void *stack,
			  size_t stack_size)
{
  return ACE_Thread_Manager::spawn (&ThreadManager::template thread_entry_point_10
				    <Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10>,
				    (void*)fun_data,
				    flags,
				    t,
				    t_handle,
				    priority,
				    grp_id,
				    stack,
				    stack_size);
};


template <typename ObjectClass>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data0<ObjectClass> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point0<ObjectClass>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
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
				      &ThreadManager::template thread_entry_point1<ObjectClass,Arg>,
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
				      &ThreadManager::template thread_entry_point2<ObjectClass,Arg1,Arg2>,
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
				      &ThreadManager::template thread_entry_point3<ObjectClass,Arg1,Arg2,Arg3>,
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
				      &ThreadManager::template thread_entry_point4<ObjectClass,Arg1,Arg2,Arg3,Arg4>,
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
				      &ThreadManager::template thread_entry_point5<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5>,
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
				      &ThreadManager::template thread_entry_point6<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data7<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                               Arg6,Arg7> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point7<ObjectClass,
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data8<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                               Arg6,Arg7,Arg8> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point8<ObjectClass,
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data9<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                              Arg6,Arg7,Arg8,Arg9> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point9<ObjectClass,
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9>,
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
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9, typename Arg10>
int ThreadManager::spawn_n (size_t n,
			    Mem_Fun_Data10<ObjectClass,Arg1,Arg2,Arg3,Arg4,Arg5,
			                               Arg6,Arg7,Arg8,Arg9,Arg10> *mem_fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point10<ObjectClass,
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10>,
				      (void*)mem_fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};


template <typename Arg1>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data1<Arg1> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_1<
				      Arg1>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data2<Arg1,Arg2> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_2<
				      Arg1,Arg2>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data3<Arg1,Arg2,Arg3> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_3<
				      Arg1,Arg2,Arg3>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data4<Arg1,Arg2,Arg3,Arg4> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_4<
				      Arg1,Arg2,Arg3,Arg4>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_5<
				      Arg1,Arg2,Arg3,Arg4,Arg5>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,
			    Arg6> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_6<
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,
			    Arg6,Arg7> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_7<
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,
			    Arg6,Arg7,Arg8> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_8<
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};

template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,
			    Arg6,Arg7,Arg8,Arg9> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_9<
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};


template <typename Arg1, typename Arg2,
                                typename Arg3, typename Arg4,
                                typename Arg5, typename Arg6,
                                typename Arg7, typename Arg8,
                                typename Arg9, typename Arg10>
int ThreadManager::spawn_n (size_t n,
			    Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,
			    Arg6,Arg7,Arg8,Arg9,Arg10> *fun_data,
			    long flags,
			    long priority,
			    int grp_id,
			    ACE_Task_Base *task,
			    ACE_hthread_t thread_handles[],
			    void *stack[],
			    size_t stack_size[]) 
{
  return ACE_Thread_Manager::spawn_n (n,
				      &ThreadManager::template thread_entry_point_10<
				      Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10>,
				      (void*)fun_data,
				      flags,
				      priority,
				      grp_id,
				      task,
				      thread_handles,
				      stack,
				      stack_size);
};


template <typename Class>
void * ThreadManager::thread_entry_point0 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data0<Class> *arg = reinterpret_cast<Mem_Fun_Data0<Class> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))();
  else
    {
      (arg->object->*(arg->void_mem_fun))();
      return 0;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg);
      return 0;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg1, arg->arg2);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg1, arg->arg2);
      return 0;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3);
      return 0;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4);
      return 0;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4,
					  arg->arg5);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4,
					  arg->arg5);
      return;
    };
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
  if (arg->mem_fun != 0)
    return (arg->object->*(arg->mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4,
					  arg->arg5,
					  arg->arg6);
  else
    {
      (arg->object->*(arg->void_mem_fun))(arg->arg1,
					  arg->arg2,
					  arg->arg3,
					  arg->arg4,
					  arg->arg5,
					  arg->arg6);
      return 0;
    };
};


template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7>
void * ThreadManager::thread_entry_point7 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data7<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7> *arg
    = reinterpret_cast<Mem_Fun_Data7<Class,Arg1,Arg2,Arg3,Arg4,Arg5,
                                            Arg6,Arg7> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5,
					arg->arg6,
					arg->arg7);
};


template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8>
void * ThreadManager::thread_entry_point8 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data8<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8> *arg
    = reinterpret_cast<Mem_Fun_Data8<Class,Arg1,Arg2,Arg3,Arg4,Arg5,
                                            Arg6,Arg7,Arg8> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5,
					arg->arg6,
					arg->arg7,
					arg->arg8);
};


template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8,
                          typename Arg9>
void * ThreadManager::thread_entry_point9 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data9<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9> *arg
    = reinterpret_cast<Mem_Fun_Data9<Class,Arg1,Arg2,Arg3,Arg4,Arg5,
                                            Arg6,Arg7,Arg8,Arg9> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5,
					arg->arg6,
					arg->arg7,
					arg->arg8,
					arg->arg9);
};


template <typename Class, typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8,
                          typename Arg9, typename Arg10>
void * ThreadManager::thread_entry_point10 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Mem_Fun_Data10<Class,Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10> *arg
    = reinterpret_cast<Mem_Fun_Data10<Class,Arg1,Arg2,Arg3,Arg4,Arg5,
                                            Arg6,Arg7,Arg8,Arg9,Arg10> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->object->*(arg->mem_fun))(arg->arg1,
					arg->arg2,
					arg->arg3,
					arg->arg4,
					arg->arg5,
					arg->arg6,
					arg->arg7,
					arg->arg8,
					arg->arg9,
					arg->arg10);
};


template <typename Arg1>
void * ThreadManager::thread_entry_point_1 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data1<Arg1> *arg
    = reinterpret_cast<Fun_Data1<Arg1> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->fun_ptr != 0)
    return (arg->fun_ptr)(arg->arg1);
   else
    {
      (arg->void_fun_ptr)(arg->arg1);
      return 0;
    }; 
};


template <typename Arg1, typename Arg2>
void * ThreadManager::thread_entry_point_2 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data2<Arg1,Arg2> *arg
    = reinterpret_cast<Fun_Data2<Arg1,Arg2> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->fun_ptr != 0)
    return (arg->fun_ptr)(arg->arg1,
			  arg->arg2);
  else
    {
      (arg->void_fun_ptr)(arg->arg1,
			  arg->arg2);
      return 0;
    };
};


template <typename Arg1, typename Arg2,
  typename Arg3>
void * ThreadManager::thread_entry_point_3 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data3<Arg1,Arg2,Arg3> *arg
    = reinterpret_cast<Fun_Data3<Arg1,Arg2,Arg3> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->fun_ptr != 0)
    return (arg->fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3);
  else
    {
      (arg->void_fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3);
      return 0;
    };
};


template <typename Arg1, typename Arg2,
  typename Arg3, typename Arg4>
void * ThreadManager::thread_entry_point_4 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data4<Arg1,Arg2,Arg3,Arg4> *arg
    = reinterpret_cast<Fun_Data4<Arg1,Arg2,Arg3,Arg4> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->fun_ptr != 0)
    return (arg->fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3,
			  arg->arg4);
  else
    {
      (arg->void_fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3,
			  arg->arg4);
      return 0;
    };
};


template <typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5>
void * ThreadManager::thread_entry_point_5 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *arg
    = reinterpret_cast<Fun_Data5<Arg1,Arg2,Arg3,Arg4,Arg5> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  if (arg->fun_ptr != 0)
    return (arg->fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3,
			  arg->arg4,
			  arg->arg5);
  else
    {
      (arg->void_fun_ptr)(arg->arg1,
			  arg->arg2,
			  arg->arg3,
			  arg->arg4,
			  arg->arg5);
      return 0;
    };
};


template <typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6>
void * ThreadManager::thread_entry_point_6 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,Arg6> *arg
    = reinterpret_cast<Fun_Data6<Arg1,Arg2,Arg3,Arg4,Arg5,
    Arg6> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->fun_ptr)(arg->arg1,
			arg->arg2,
			arg->arg3,
			arg->arg4,
			arg->arg5,
			arg->arg6);
};


template <typename Arg1, typename Arg2,
  typename Arg3, typename Arg4,
  typename Arg5, typename Arg6,
  typename Arg7>
void * ThreadManager::thread_entry_point_7 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7> *arg
    = reinterpret_cast<Fun_Data7<Arg1,Arg2,Arg3,Arg4,Arg5,
    Arg6,Arg7> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->fun_ptr)(arg->arg1,
			 arg->arg2,
			 arg->arg3,
			 arg->arg4,
			 arg->arg5,
			 arg->arg6,
			arg->arg7);
};


template <typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8>
void * ThreadManager::thread_entry_point_8 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8> *arg
    = reinterpret_cast<Fun_Data8<Arg1,Arg2,Arg3,Arg4,Arg5,
    Arg6,Arg7,Arg8> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->fun_ptr)(arg->arg1,
			 arg->arg2,
			 arg->arg3,
			 arg->arg4,
			 arg->arg5,
			 arg->arg6,
			 arg->arg7,
			 arg->arg8);
};


template <typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8,
                          typename Arg9>
void * ThreadManager::thread_entry_point_9 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9> *arg
    = reinterpret_cast<Fun_Data9<Arg1,Arg2,Arg3,Arg4,Arg5,
    Arg6,Arg7,Arg8,Arg9> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->fun_ptr)(arg->arg1,
			 arg->arg2,
			 arg->arg3,
			 arg->arg4,
			 arg->arg5,
			 arg->arg6,
			 arg->arg7,
			 arg->arg8,
			 arg->arg9);
};

template <typename Arg1, typename Arg2,
                          typename Arg3, typename Arg4,
                          typename Arg5, typename Arg6,
                          typename Arg7, typename Arg8,
                          typename Arg9, typename Arg10>
void * ThreadManager::thread_entry_point_10 (void *_arg)
{
				   // reinterpret the given pointer as
				   // a pointer to the structure
				   // containing all the necessary
				   // information
  Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,Arg6,Arg7,Arg8,Arg9,Arg10> *arg
    = reinterpret_cast<Fun_Data10<Arg1,Arg2,Arg3,Arg4,Arg5,
    Arg6,Arg7,Arg8,Arg9,Arg10> *>(_arg);

				   // extract function pointer, object
				   // and argument and dispatch the
				   // call
  return (arg->fun_ptr)(arg->arg1,
			 arg->arg2,
			 arg->arg3,
			 arg->arg4,
			 arg->arg5,
			 arg->arg6,
			 arg->arg7,
			 arg->arg8,
			 arg->arg9,
			 arg->arg10);
};


#endif


#endif
