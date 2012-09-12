INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckIncludeFiles)

#
# Check for various C++ language features
#


CHECK_CXX_SOURCE_COMPILES( # TODO: Phase out
  "
  #include <ostream>
  int main(){return 0;}
  "
  HAVE_STD_OSTREAM_HEADER)

CHECK_CXX_SOURCE_COMPILES(
  "
  #include <iosfwd>
  int main(){return 0;}
  "
  HAVE_STD_IOSFWD_HEADER)




#
# C++0x Support:
#


# See if there is a compiler flag to enable C++0x features
# (Only test for -std=c++0x for the moment.)

# Set CMAKE_REQUIRED_FLAGS for the unit tests
LIST(APPEND CMAKE_REQUIRED_FLAGS "-std=c++0x")

CHECK_CXX_SOURCE_COMPILES(
  "
  int main(){ return 0; }
  "
  DEAL_II_HAVE_CXX0X_FLAG)



IF(DEAL_II_HAVE_CXX0X_FLAG)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <array>
    std::array<int,3> p;
    int main(){  p[0]; return 0; }
    "
    DEAL_II_HAVE_CXX0X_ARRAY)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <condition_variable>
    std::condition_variable c;
    int main(){ c.notify_all(); return 0; }
    "
    DEAL_II_HAVE_CXX0X_CONDITION_VARIABLE)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <functional>
    void f(int, double){}
    std::function<void (int)> g = std::bind (f, std::placeholders::_1,1.1);
    int main(){ return 0; }
    "
    DEAL_II_HAVE_CXX0X_FUNCTIONAL)

  # Make sure we don't run into GCC bug 35569
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <functional>
    void f(int){}
    using namespace std;
    using namespace std::placeholders;
    int main(){ bind(multiplies<int>(),4,_1)(5); return 0; }
    "
    DEAL_II_HAVE_CXX0X_FUNCTIONAL_GCCBUG35569_OK)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <memory>
    std::shared_ptr<int> p(new int(3));
    int main(){ return 0; }
    "
    DEAL_II_HAVE_CXX0X_SHARED_PTR)

  LIST(APPEND CMAKE_REQUIRED_FLAGS "-lpthreads")

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <thread>
    void f(int){}
    int main(){ std::thread t(f,1); t.join(); return 0; }
    "
    DEAL_II_HAVE_CXX0X_THREAD)

  #On some systems with gcc 4.5.0, we can compile the code
  #above but it will throw an exception when run. So test
  #that as well. The test will only be successful if we have
  #libpthread available, so link it in for this test. If
  #multithreading is requested, it will be added to CXXFLAGS
  #later on so there is no need to do this here.
  CHECK_CXX_SOURCE_RUNS(
    "
    #include <thread>
    void f(int){}
    int main(){ std::thread t(f,1); t.join(); return 0; }
    "
    DEAL_II_HAVE_CXX0X_THREAD_RUN_OK)

  LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "-lpthreads")

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <mutex>
    std::mutex m;
    int main(){ m.lock(); return 0; }
    "
    DEAL_II_HAVE_CXX0X_MUTEX)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <tuple>],
    std::tuple<int,double,char> p(1,1.1,'a');
    int main(){ return 0; }
    "
    DEAL_II_HAVE_CXX0X_TUPLE)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <type_traits>
    const bool m0 = std::is_trivial<double>::value;
    const bool m1 = std::is_standard_layout<double>::value;
    const bool m2 = std::is_pod<double>::value;
    int main(){ return 0; }
    "
    DEAL_II_HAVE_CXX0X_TYPE_TRAITS)

  IF( DEAL_II_HAVE_CXX0X_ARRAY AND
      DEAL_II_HAVE_CXX0X_CONDITION_VARIABLE AND
      DEAL_II_HAVE_CXX0X_FUNCTIONAL AND
      DEAL_II_HAVE_CXX0X_FUNCTIONAL_GCCBUG35569_OK AND
      DEAL_II_HAVE_CXX0X_SHARED_PTR AND
      DEAL_II_HAVE_CXX0X_THREAD AND
      DEAL_II_HAVE_CXX0X_THREAD_RUN_OK AND
      DEAL_II_HAVE_CXX0X_MUTEX AND
      DEAL_II_HAVE_CXX0X_TUPLE AND
      DEAL_II_HAVE_CXX0X_TYPE_TRAITS )

    MESSAGE(STATUS "Sufficient C++0x support. Enabling -std=c++0x.")

    SET(DEAL_II_CAN_USE_CXX1X TRUE)

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") # TODO

  ELSE()

    MESSAGE(STATUS "Insufficient C++0x support. Disabling -std=c++0x.")
  ENDIF()

  IF(DEAL_II_CAN_USE_CXX1X)
    #
    # Also test for a couple C++0x things that we don't use in the library
    # but that users may want to use in their applications and that we
    # might want to test in the testsuite
    #

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <vector>
      std::vector<int> v;
      int main(){ auto i = v.begin(); *i; return 0;}
      "
      DEAL_II_HAVE_CXX0X_AUTO_TYPE)

    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <vector>],
      std::vector<int> v;
      int main(){ for (std::vector<int>::iterator i : v) *i; return 0;}
      "
      DEAL_II_HAVE_CXX0X_RANGE_BASED_FOR)

    IF( DEAL_II_HAVE_CXX0X_AUTO_TYPE AND
        DEAL_II_HAVE_CXX0X_RANGE_BASED_FOR )

      MESSAGE(STATUS "Additional C++0x support available.")

      SET(DEAL_II_CAN_USE_ADDITIONAL_CXX1X_FEATURES)

    ENDIF()
  ENDIF()

ENDIF()

LIST(REMOVE_ITEM CMAKE_REQUIRED_FLAGS "-std=c++0x")



#
# Checks for various header files:
#

CHECK_INCLUDE_FILES("stdint.h" HAVE_STDINT_H)
CHECK_INCLUDE_FILES("stdlib.h" HAVE_STDLIB_H)
CHECK_INCLUDE_FILES("strings.h" HAVE_STRINGS_H)
CHECK_INCLUDE_FILES("string.h" HAVE_STRING_H)
CHECK_INCLUDE_FILES("sys/stat.h" HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES("sys/syscall.h" HAVE_SYS_SYSCALL_H)
CHECK_INCLUDE_FILES("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES("sys/types.h" HAVE_SYS_TYPES_H)

