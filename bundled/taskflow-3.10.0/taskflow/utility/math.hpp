#pragma once

#include <atomic>
#include <chrono>

namespace tf {

/**
 * @brief rounds the given 64-bit unsigned integer to the nearest power of 2
 */
template <typename T, std::enable_if_t<
  (std::is_unsigned_v<std::decay_t<T>> && sizeof(T) == 8), void
>* = nullptr>
constexpr T next_pow2(T x) {
  if(x == 0) return 1;
  x--;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  x++;
  return x;
}

/**
 * @brief rounds the given 32-bit unsigned integer to the nearest power of 2
 */
template <typename T, std::enable_if_t<
  (std::is_unsigned_v<std::decay_t<T>> && sizeof(T) == 4), void
>* = nullptr>
constexpr T next_pow2(T y) {
  if(y == 0) return 1;
  y--;
  y |= y >> 1;
  y |= y >> 2;
  y |= y >> 4;
  y |= y >> 8;
  y |= y >> 16;
  y++;
  return y;
}

/**
 * @brief checks if the given number is a power of 2
 *
 * This function determines if the given integer is a power of 2.
 *
 * @tparam T The type of the input. Must be an integral type.
 * @param x The integer to check.
 * @return `true` if `x` is a power of 2, otherwise `false`.
 *
 * @attention This function is constexpr and can be evaluated at compile time.
 *
 */
template <typename T, std::enable_if_t<
  std::is_integral_v<std::decay_t<T>>, void>* = nullptr
>
constexpr bool is_pow2(const T& x) {
  return x && (!(x&(x-1)));
}

/**
 * @brief computes the floor of the base-2 logarithm of a number using count-leading-zeros (CTL).
 *
 * This function efficiently calculates the floor of `log2(n)` for both 32-bit and 64-bit integers.
 *
 * @tparam T integer type (uint32_t or uint64_t).
 * @param n input number.
 * @return floor of `log2(n)`
 */
template <typename T>
constexpr size_t floor_log2(T n) {

   static_assert(std::is_unsigned_v<T>, "log2 only supports unsigned integer types");

#if defined(_MSC_VER)
  unsigned long index;
  if constexpr (sizeof(T) == 8) {
    _BitScanReverse64(&index, n);
  } else {
    _BitScanReverse(&index, static_cast<unsigned long>(n));
  }
  return static_cast<size_t>(index);
#elif defined(__GNUC__) || defined(__clang__)
  if constexpr (sizeof(T) == 8) {
    return 63 - __builtin_clzll(n);
  } else {
    return 31 - __builtin_clz(n);
  }
#else
  // Portable fallback: Uses bit shifts to count leading zeros manually
  size_t log = 0;
  while (n >>= 1) {
    ++log;
  }
  return log;
#endif
}

/**
@brief returns the floor of `log2(N)` at compile time 
*/
template<size_t N>
constexpr size_t static_floor_log2() {
  return (N < 2) ? 0 : 1 + static_floor_log2<N / 2>();
  //auto log = 0;
  //while (N >>= 1) {
  //  ++log;
  //}
  //return log;
}

/**
 * @brief finds the median of three numbers pointed to by iterators using the given comparator
 *
 * This function determines the median value of the elements pointed to by
 * three random-access iterators using the provided comparator.
 *
 * @tparam RandItr The type of the random-access iterator.
 * @tparam C The type of the comparator.
 * @param l Iterator to the first element.
 * @param m Iterator to the second element.
 * @param r Iterator to the third element.
 * @param cmp The comparator used to compare the dereferenced iterator values.
 * @return The iterator pointing to the median value among the three elements.
 *
 */
template <typename RandItr, typename C>
RandItr median_of_three(RandItr l, RandItr m, RandItr r, C cmp) {
  return cmp(*l, *m) ? (cmp(*m, *r) ? m : (cmp(*l, *r) ? r : l ))
                     : (cmp(*r, *m) ? m : (cmp(*r, *l) ? r : l ));
}

/**
 * @brief finds the pseudo median of a range of items using a spread of nine numbers
 *
 * This function computes an approximate median of a range of items by sampling
 * nine values spread across the range and finding their median. It uses a
 * combination of the `median_of_three` function to determine the pseudo median.
 *
 * @tparam RandItr The type of the random-access iterator.
 * @tparam C The type of the comparator.
 * @param beg Iterator to the beginning of the range.
 * @param end Iterator to the end of the range.
 * @param cmp The comparator used to compare the dereferenced iterator values.
 * @return The iterator pointing to the pseudo median of the range.
 *
 * @attention The pseudo median is an approximation of the true median and may not
 *       be the exact middle value of the range.
 *
 */
template <typename RandItr, typename C>
RandItr pseudo_median_of_nine(RandItr beg, RandItr end, C cmp) {
  size_t N = std::distance(beg, end);
  size_t offset = N >> 3;
  return median_of_three(
    median_of_three(beg, beg+offset, beg+(offset*2), cmp),
    median_of_three(beg+(offset*3), beg+(offset*4), beg+(offset*5), cmp),
    median_of_three(beg+(offset*6), beg+(offset*7), end-1, cmp),
    cmp
  );
}

/**
 * @brief sorts two elements of dereferenced iterators using the given comparison function
 *
 * This function compares two elements pointed to by iterators and swaps them
 * if they are out of order according to the provided comparator.
 *
 * @tparam Iter The type of the iterator.
 * @tparam Compare The type of the comparator.
 * @param a Iterator to the first element.
 * @param b Iterator to the second element.
 * @param comp The comparator used to compare the dereferenced iterator values.
 *
 */
template<typename Iter, typename Compare>
void sort2(Iter a, Iter b, Compare comp) {
  if (comp(*b, *a)) std::iter_swap(a, b);
}

/**
 * @brief Sorts three elements of dereferenced iterators using the given comparison function.
 *
 * This function sorts three elements pointed to by iterators in ascending order
 * according to the provided comparator. The sorting is performed using a sequence
 * of calls to the `sort2` function to ensure the correct order of elements.
 *
 * @tparam Iter The type of the iterator.
 * @tparam Compare The type of the comparator.
 * @param a Iterator to the first element.
 * @param b Iterator to the second element.
 * @param c Iterator to the third element.
 * @param comp The comparator used to compare the dereferenced iterator values.
 *
 */
template<typename Iter, typename Compare>
void sort3(Iter a, Iter b, Iter c, Compare comp) {
  sort2(a, b, comp);
  sort2(b, c, comp);
  sort2(a, b, comp);
}

/**
 * @brief generates a program-wide unique ID of the given type in a thread-safe manner
 *
 * This function provides a globally unique identifier of the specified integral type.
 * It uses a static `std::atomic` counter to ensure thread safety and increments the
 * counter in a relaxed memory ordering for efficiency.
 *
 * @tparam T The type of the ID to generate. Must be an integral type.
 * @return A unique ID of type `T`.
 *
 * @attention The uniqueness of the ID is guaranteed only within the program's lifetime.
 * @attention The function does not throw exceptions.
 *
 */
template <typename T, std::enable_if_t<std::is_integral_v<T>, void>* = nullptr>
T unique_id() {
  static std::atomic<T> counter{0};
  return counter.fetch_add(1, std::memory_order_relaxed);
}

/**
 * @brief updates an atomic variable with the maximum value
 *
 * This function atomically updates the provided atomic variable `v` to hold
 * the maximum of its current value and `max_v`. The update is performed using
 * a relaxed memory ordering for efficiency in non-synchronizing contexts.
 *
 * @tparam T The type of the atomic variable. Must be trivially copyable and comparable.
 * @param v The atomic variable to update.
 * @param max_v The value to compare with the current value of `v`.
 *
 * @attention If multiple threads call this function concurrently, the value of `v`
 *       will be the maximum value seen across all threads.
 *
 */
template <typename T>
inline void atomic_max(std::atomic<T>& v, const T& max_v) noexcept {
  T prev = v.load(std::memory_order_relaxed);
  while(prev < max_v && 
        !v.compare_exchange_weak(prev, max_v, std::memory_order_relaxed,
                                              std::memory_order_relaxed)) {
  }
}

/**
 * @brief updates an atomic variable with the minimum value
 *
 * This function atomically updates the provided atomic variable `v` to hold
 * the minimum of its current value and `min_v`. The update is performed using 
 * a relaxed memory ordering for efficiency in non-synchronizing contexts.
 *
 * @tparam T The type of the atomic variable. Must be trivially copyable and comparable.
 * @param v The atomic variable to update.
 * @param min_v The value to compare with the current value of `v`.
 *
 * @attention If multiple threads call this function concurrently, the value of `v` 
 *       will be the minimum value seen across all threads.
 *
 */
template <typename T>
inline void atomic_min(std::atomic<T>& v, const T& min_v) noexcept {
  T prev = v.load(std::memory_order_relaxed);
  while(prev > min_v && 
        !v.compare_exchange_weak(prev, min_v, std::memory_order_relaxed,
                                              std::memory_order_relaxed)) {
  }
}

/**
 * @brief generates a random seed based on the current system clock
 *
 * This function returns a seed value derived from the number of clock ticks
 * since the epoch as measured by the system clock. The seed can be used
 * to initialize random number generators.
 *
 * @tparam T The type of the returned seed. Must be an integral type.
 * @return A seed value based on the system clock.
 *
 */
template <typename T>
inline T seed() noexcept {
  return std::chrono::system_clock::now().time_since_epoch().count();
}

/**
 * @brief counts the number of trailing zeros in an integer.
 *
 * This function provides a portable implementation for counting the number of 
 * trailing zeros across different platforms and integer sizes (32-bit and 64-bit).
 *
 * @tparam T integer type (32-bit or 64-bit).
 * @param x non-zero integer to count trailing zeros from
 * @return the number of trailing zeros in @c x
 *
 * @attention
 * The behavior is undefined when @c x is 0.
 */
template <typename T, typename = std::enable_if_t<std::is_unsigned_v<T>>>
auto ctz(T x) {

  #if defined(_MSC_VER)
    unsigned long index;
    if constexpr (sizeof(T) == 8) {
      _BitScanForward64(&index, x);
    } else {
      _BitScanForward(&index, (unsigned long)x);
    }
    return index;
  #elif defined(__GNUC__) || defined(__clang__)
    if constexpr (sizeof(T) == 8) {
      return __builtin_ctzll(x);
    } else {
      return __builtin_ctz(x);
    }
  #else 
    size_t r = 0;
    while ((x & 1) == 0) {
      x >>= 1;
      r++;
    }
    return r;
  #endif
}

// ------------------------------------------------------------------------------------------------
// coprime
// ------------------------------------------------------------------------------------------------

/**
 * @brief computes a coprime of a given number
 *
 * This function finds the largest number less than N that is coprime (i.e., has a greatest common divisor of 1) with @c N.
 * If @c N is less than 3, it returns 1 as a default coprime.
 *
 * @param N input number for which a coprime is to be found.
 * @return the largest number < @c N that is coprime to N
 */
constexpr size_t coprime(size_t N) {
  if(N < 3) {
    return 1;
  }
  for (size_t x = N; --x > 0;) {
    if (std::gcd(x, N) == 1) {
      return x;
    }
  }
  return 1;
}

/**
 * @brief generates a compile-time array of coprimes for numbers from 0 to N-1
 *
 * This function constructs a constexpr array where each element at index `i` contains a coprime of `i`
 * (the largest number less than `i` that is coprime to it).
 *
 * @tparam N the size of the array to generate (should be greater than 0).
 * @return a constexpr array of size @c N where each index holds a coprime of its value.
 */
template <size_t N>
constexpr std::array<size_t, N> make_coprime_lut() {
  static_assert(N>0, "N must be greater than 0");
  std::array<size_t, N> coprimes{};
  for (size_t n = 0; n < N; ++n) {
    coprimes[n] = coprime(n);
  }
  return coprimes;
}


//class XorShift64 {
//
//  public:
//  
//  explicit XorShift64(uint64_t seed) : _state(seed) {}
//
//  uint64_t next() {
//    _state ^= _state >> 12;
//    _state ^= _state << 25;
//    _state ^= _state >> 27;
//    return _state * 0x2545F4914F6CDD1DULL; // Scramble for better randomness
//  }
//
//  size_t random_range(size_t min, size_t max) {
//    return min + (next() % (max - min + 1));
//  }
//
//  private:
//
//  uint64_t _state;
//};

//inline int generate_random_excluding(int worker_id, int W, XorShift64& rng) {
//    int random_number = rng.random_range(0, 2 * W - 2); // Range: [0, 2W-2]
//    return random_number + (random_number >= worker_id); // Skip worker_id
//}
//
//
//class Xoroshiro128Plus {
//
//  public:
//
//    explicit Xoroshiro128Plus(uint64_t seed1, uint64_t seed2) : _state{seed1, seed2} {}
//
//    uint64_t next() {
//      uint64_t s0 = _state[0];
//      uint64_t s1 = _state[1];
//      uint64_t result = s0 + s1;
//
//      s1 ^= s0;
//      _state[0] = _rotl(s0, 55) ^ s1 ^ (s1 << 14); // Scramble _state
//      _state[1] = _rotl(s1, 36);
//
//      return result;
//    }
//
//    int random_range(int min, int max) {
//      return min + (next() % (max - min + 1));
//    }
//
//  private:
//
//    std::array<uint64_t, 2> _state;
//
//    static uint64_t _rotl(uint64_t x, int k) {
//      return (x << k) | (x >> (64 - k));
//    }
//};


}  // end of namespace tf -----------------------------------------------------



