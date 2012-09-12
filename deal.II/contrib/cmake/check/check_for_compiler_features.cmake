INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckIncludeFiles)


#
# Check for various compiler features.
#



#
# Check whether the std::vector::iterator is just a plain pointer
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <vector>
  template <typename T> void f(T) {}
  template void f(int *);
  template void f(std::vector<int>::iterator);
  int main(){return 0;}
  "
  DEAL_II_VECTOR_ITERATOR_IS_NOT_POINTER)

IF(NOT DEAL_II_VECTOR_ITERATOR_IS_NOT_POINTER)
  SET(DEAL_II_VECTOR_ITERATOR_IS_POINTER TRUE)
ENDIF()



#
# Check for existence of the __builtin_expect facility of newer
# gcc compilers. This can be used to hint the compiler's branch
# prediction unit in some cases. We use it in the AssertThrow
# macros.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  bool f() {}
  int main(){ if (__builtin_expect(f(),false)) ; }
  "
  HAVE_BUILTIN_EXPECT)



#
# Newer versions of gcc have a very nice feature: you can set
# a verbose terminate handler, that not only aborts a program
# when an exception is thrown and not caught somewhere, but
# before aborting it prints that an exception has been thrown,
# and possibly what the std::exception::what() function has to
# say. Since many people run into the trap of not having a
# catch clause in main(), they wonder where that abort may be
# coming from.  The terminate handler then at least says what is
# missing in their program.
#
# This test checks whether this feature is available.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <exception>
  namespace __gnu_cxx
  {
    extern void __verbose_terminate_handler ();
  }
  struct preload_terminate_dummy
  {
    preload_terminate_dummy()
    {
      std::set_terminate (__gnu_cxx::__verbose_terminate_handler);
    }
  };
  static preload_terminate_dummy dummy;
  int main() { throw 1; return 0; }
  "
  HAVE_VERBOSE_TERMINATE)



#
# Check whether glibc-like stacktrace information is available
# for the Exception class. If it is, then try to also determine
# whether the compiler accepts the -rdynamic flag, since that is
# recommended for linking if one wants to have meaningful
# backtraces.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <execinfo.h>
  #include <stdlib.h>
  void * array[25];
  int nSize = backtrace(array, 25);
  char ** symbols = backtrace_symbols(array, nSize);
  int main(){ free(symbols); return 0; }
  "
  HAVE_GLIBC_STACKTRACE)

# On Mac OS X, -rdynamic is accepted by the compiler (i.e.
# it doesn't produce an error) but we always get a warning
# that it isn't supported. That's pretty stupid because
# we can't test for it. Consequently, only run the test
# if not on OS X.

IF(HAVE_GLIBC_STACKTRACE AND NOT CMAKE_SYSTEM_NAME MATCHES "Darwin")

  LIST(APPEND CMAKE_REQUIRED_FLAGS "-rdynamic")

  CHECK_CXX_SOURCE_COMPILES(
    "
    int main() { return 0; }
    "
    DEAL_II_COMPILER_HAS_RDYNAMIC)

  LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "-rdynamic")

  IF(DEAL_II_COMPILER_HAS_RDYNAMIC)
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -rdynamic")
  ENDIF()

ENDIF()



#
# Check whether the compiler offers a way to demangle symbols
# from within the program. Used inside the exception stacktrace
# mechanism.
#
# The example code is taken from
#   http://gcc.gnu.org/onlinedocs/libstdc++/18_support/howto.html#6
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <exception>
  #include <iostream>
  #include <cxxabi.h>
  #include <cstdlib>

  struct empty { };

  template <typename T, int N>
  struct bar { };

  int     status;
  char   *realname;

  int main()
  {
    // exception classes not in <stdexcept>, thrown by the implementation
    // instead of the user
    std::bad_exception  e;
    realname = abi::__cxa_demangle(e.what(), 0, 0, &status);
    free(realname);


    // typeid
    bar<empty,17>          u;
    const std::type_info  &ti = typeid(u);

    realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
    free(realname);

      return 0;
  }
  "
  HAVE_LIBSTDCXX_DEMANGLER)
