#pragma once

#include "../taskflow.hpp"

namespace tf {

// Function: make_reduce_task
template <typename B, typename E, typename T, typename O, typename P = DefaultPartitioner>
auto make_reduce_task(B b, E e, T& init, O bop, P part = P()) {
  
  using namespace std::string_literals;

  using B_t = std::decay_t<unwrap_ref_decay_t<B>>;
  using E_t = std::decay_t<unwrap_ref_decay_t<E>>;

  return [=, &init] (Runtime& rt) mutable {

    // fetch the iterator values
    B_t beg = b;
    E_t end = e;

    size_t W = rt.executor().num_workers();
    size_t N = std::distance(beg, end);

    // only myself - no need to spawn another graph
    if(W <= 1 || N <= part.chunk_size()) {
      part([=, &init] () mutable { for(; beg!=end; init = bop(init, *beg++)); })();
      return;
    }
    
    PreemptionGuard preemption_guard(rt);

    if(N < W) {
      W = N;
    }

    auto mutex = std::make_shared<std::mutex>();

    // static partitioner
    if constexpr(part.type() == PartitionerType::STATIC) {
      
      for(size_t w=0, curr_b=0; w<W && curr_b < N;) {
        
        // we force chunk size to be at least two because the temporary
        // variable sum need to avoid copy at the first step
        auto chunk_size = std::max(size_t{2}, part.adjusted_chunk_size(N, W, w));
        
        auto task = part([=, &init] () mutable {

          std::advance(beg, curr_b);

          if(N - curr_b == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            init = bop(init, *beg);
            return;
          }

          auto beg1 = beg++;
          auto beg2 = beg++;
          T sum = bop(*beg1, *beg2);
        
          // loop reduce
          part.loop(N, W, curr_b, chunk_size,
            [=, &sum, prev_e=curr_b+2](size_t part_b, size_t part_e) mutable {

              if(part_b > prev_e) {
                std::advance(beg, part_b - prev_e);
              }
              else {
                part_b = prev_e;
              }

              for(size_t x=part_b; x<part_e; x++, beg++) {
                sum = bop(sum, *beg);
              }
              prev_e = part_e;
            }
          ); 
          
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          init = bop(init, sum);
        });

        (++w == W || (curr_b += chunk_size) >= N) ? task() : rt.silent_async(task);
      }
    }
    // dynamic partitioner
    else {
      auto next = std::make_shared<std::atomic<size_t>>(0);
      
      for(size_t w=0; w<W;) {

        auto task = part([=, &init] () mutable {
          // pre-reduce
          size_t s0 = next->fetch_add(2, std::memory_order_relaxed);

          if(s0 >= N) {
            return;
          }

          std::advance(beg, s0);

          if(N - s0 == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            init = bop(init, *beg);
            return;
          }

          auto beg1 = beg++;
          auto beg2 = beg++;

          T sum = bop(*beg1, *beg2);
          
          // loop reduce
          part.loop(N, W, *next, 
            [=, &sum, prev_e=s0+2](size_t curr_b, size_t curr_e) mutable {
              std::advance(beg, curr_b - prev_e);
              for(size_t x=curr_b; x<curr_e; x++, beg++) {
                sum = bop(sum, *beg);
              }
              prev_e = curr_e;
            }
          ); 
          
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          init = bop(init, sum);
        });
        (++w == W) ? task() : rt.silent_async(task);
      }
    }
  };
}

// Function: make_transform_reduce_task
template <
  typename B, typename E, typename T, typename BOP, typename UOP, 
  typename P = DefaultPartitioner
>
auto make_transform_reduce_task(B b, E e, T& init, BOP bop, UOP uop, P part = P()) {

  using namespace std::string_literals;
  using B_t = std::decay_t<unwrap_ref_decay_t<B>>;
  using E_t = std::decay_t<unwrap_ref_decay_t<E>>;

  return [=, &init] (Runtime& rt) mutable {

    // fetch the iterator values
    B_t beg = b;
    E_t end = e;

    size_t W = rt.executor().num_workers();
    size_t N = std::distance(beg, end);

    // only myself - no need to spawn another graph
    if(W <= 1 || N <= part.chunk_size()) {
      part([=, &init] () mutable { for(; beg!=end; init = bop(std::move(init), uop(*beg++))); })();
      return;
    }
    
    PreemptionGuard preemption_guard(rt);

    if(N < W) {
      W = N;
    }

    auto mutex = std::make_shared<std::mutex>();
    
    // static partitioner
    if constexpr(part.type() == PartitionerType::STATIC) {

      for(size_t w=0, curr_b=0; w<W && curr_b < N;) {
      
        auto chunk_size = part.adjusted_chunk_size(N, W, w);

        auto task = part([=, &init] () mutable {

          std::advance(beg, curr_b);

          if(N - curr_b == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            init = bop(std::move(init), uop(*beg));
            return;
          }

          //auto beg1 = beg++;
          //auto beg2 = beg++;
          //T sum = bop(uop(*beg1), uop(*beg2));

          T sum = (chunk_size == 1) ? uop(*beg++) : bop(uop(*beg++), uop(*beg++));
        
          // loop reduce
          part.loop(N, W, curr_b, chunk_size,
            [=, &sum, prev_e=curr_b+(chunk_size == 1 ? 1 : 2)]
            (size_t part_b, size_t part_e) mutable {
              if(part_b > prev_e) {
                std::advance(beg, part_b - prev_e);
              }
              else {
                part_b = prev_e;
              }
              for(size_t x=part_b; x<part_e; x++, beg++) {
                sum = bop(std::move(sum), uop(*beg));
              }
              prev_e = part_e;
            }
          ); 
          
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          init = bop(std::move(init), std::move(sum));
        });

        (++w == W || (curr_b += chunk_size) >= N) ? task() : rt.silent_async(task);
      }
    }
    // dynamic partitioner
    else {
      auto next = std::make_shared<std::atomic<size_t>>(0);
      for(size_t w=0; w<W;) {
        auto task = part([=, &init] () mutable {

          // pre-reduce
          size_t s0 = next->fetch_add(2, std::memory_order_relaxed);

          if(s0 >= N) {
            return;
          }

          std::advance(beg, s0);

          if(N - s0 == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            init = bop(std::move(init), uop(*beg));
            return;
          }

          auto beg1 = beg++;
          auto beg2 = beg++;

          T sum = bop(uop(*beg1), uop(*beg2));
          
          // loop reduce
          part.loop(N, W, *next, 
            [=, &sum, prev_e=s0+2](size_t curr_b, size_t curr_e) mutable {
              std::advance(beg, curr_b - prev_e);
              for(size_t x=curr_b; x<curr_e; x++, beg++) {
                sum = bop(std::move(sum), uop(*beg));
              }
              prev_e = curr_e;
            }
          ); 
          
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          init = bop(std::move(init), std::move(sum));
        });
        (++w == W) ? task() : rt.silent_async(task);
      }
    }
  };
}

// Function: make_transform_reduce_task with two binary operation
template <
  typename B1, typename E1, typename B2, typename T, typename BOP_R, typename BOP_T, 
  typename P = DefaultPartitioner,
  std::enable_if_t<!is_partitioner_v<std::decay_t<BOP_T>>, void>* = nullptr
>
auto make_transform_reduce_task(
  B1 b1, E1 e1, B2 b2, T& init, BOP_R bop_r, BOP_T bop_t, P part = P()
) {
  
  using namespace std::string_literals;

  using B1_t = std::decay_t<unwrap_ref_decay_t<B1>>;
  using E1_t = std::decay_t<unwrap_ref_decay_t<E1>>;
  using B2_t = std::decay_t<unwrap_ref_decay_t<B2>>;

  return [=, &r=init] (Runtime& rt) mutable {

    // fetch the iterator values
    B1_t beg1 = b1;
    E1_t end1 = e1;
    B2_t beg2 = b2; 

    size_t W = rt.executor().num_workers();
    size_t N = std::distance(beg1, end1);

    // only myself - no need to spawn another graph
    if(W <= 1 || N <= part.chunk_size()) {
      part([=, &r] () mutable { for(; beg1!=end1; r = bop_r(std::move(r), bop_t(*beg1++, *beg2++))); })();
      return;
    }   
    
    PreemptionGuard preemption_guard(rt);

    if(N < W) {
      W = N;
    }   

    auto mutex = std::make_shared<std::mutex>();
    
    // static partitioner
    if constexpr(part.type() == PartitionerType::STATIC) {

      for(size_t w=0, curr_b=0; w<W && curr_b < N;) {
    
        auto chunk_size = part.adjusted_chunk_size(N, W, w); 

        auto task = part([=, &r] () mutable {
          std::advance(beg1, curr_b);
          std::advance(beg2, curr_b);

          if(N - curr_b == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            r = bop_r(std::move(r), bop_t(*beg1, *beg2));
            return;
          }   

          T sum = (chunk_size == 1) ? bop_t(*beg1++, *beg2++) : 
            bop_r(bop_t(*beg1++, *beg2++), bop_t(*beg1++, *beg2++));
    
          // loop reduce
          part.loop(N, W, curr_b, chunk_size,
            [=, &sum, prev_e=curr_b+(chunk_size == 1 ? 1 : 2)] 
            (size_t part_b, size_t part_e) mutable {
              if(part_b > prev_e) {
                std::advance(beg1, part_b - prev_e);
                std::advance(beg2, part_b - prev_e);
              }   
              else {
                part_b = prev_e;
              }   
              for(size_t x=part_b; x<part_e; x++, beg1++, beg2++) { 
                sum = bop_r(std::move(sum), bop_t(*beg1, *beg2));
              }   
              prev_e = part_e;
            }   
          );  
    
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          r = bop_r(std::move(r), std::move(sum));
        });

        (++w == W || (curr_b += chunk_size) >= N) ? task() : rt.silent_async(task);
      }   
    }   
    // dynamic partitioner
    else {
      auto next = std::make_shared<std::atomic<size_t>>(0);
    
      for(size_t w=0; w<W;) {

        auto task = part([=, &r] () mutable {

          // pre-reduce
          size_t s0 = next->fetch_add(2, std::memory_order_relaxed);

          if(s0 >= N) {
            return;
          }   

          std::advance(beg1, s0);
          std::advance(beg2, s0);

          if(N - s0 == 1) {
            std::lock_guard<std::mutex> lock(*mutex);
            r = bop_r(std::move(r), bop_t(*beg1, *beg2));
            return;
          }   

          auto beg11 = beg1++;
          auto beg12 = beg1++;
          auto beg21 = beg2++;
          auto beg22 = beg2++;

          T sum = bop_r(bop_t(*beg11, *beg21), bop_t(*beg12, *beg22));

          // loop reduce
          part.loop(N, W, *next, 
            [=, &sum, prev_e=s0+2](size_t curr_b, size_t curr_e) mutable {
              std::advance(beg1, curr_b - prev_e);
              std::advance(beg2, curr_b - prev_e);
              for(size_t x=curr_b; x<curr_e; x++, beg1++, beg2++) {
                sum = bop_r(std::move(sum), bop_t(*beg1, *beg2));
              }   
              prev_e = curr_e;
            }   
          );  
    
          // final reduce
          std::lock_guard<std::mutex> lock(*mutex);
          r = bop_r(std::move(r), std::move(sum));
        });
        
        (++w == W) ? task() : rt.silent_async(task);
      }
    }   
  };  
}


// Function: make_reduce_by_index_task
template <typename R, typename T, typename L, typename G, typename P = DefaultPartitioner>
auto make_reduce_by_index_task(R range, T& init, L lop, G gop, P part = P()) {
  
  using range_type = std::decay_t<unwrap_ref_decay_t<R>>;

  return [=, &init] (Runtime& rt) mutable {

    // fetch the iterator values
    range_type r = range;
    
    // nothing to be done if the range is invalid
    if(is_index_range_invalid(r.begin(), r.end(), r.step_size())) {
      return;
    }

    size_t W = rt.executor().num_workers();
    size_t N = r.size();

    // only myself - no need to spawn another graph
    if(W <= 1 || N <= part.chunk_size()) {
      part([=, &init] () mutable { init = lop(r, std::move(init)); })();
      return;
    }
    
    PreemptionGuard preemption_guard(rt);

    if(N < W) {
      W = N;
    }

    auto mutex = std::make_shared<std::mutex>();

    // static partitioner
    if constexpr(part.type() == PartitionerType::STATIC) {
      
      for(size_t w=0, curr_b=0; w<W && curr_b < N;) {
        
        // we force chunk size to be at least two because the temporary
        // variable sum need to avoid copy at the first step
        auto chunk_size = part.adjusted_chunk_size(N, W, w);
        
        auto task = part([=, &init] () mutable {

          // temporary result so far
          std::optional<T> tmp;

          // loop reduce
          part.loop(N, W, curr_b, chunk_size, [=, &tmp](size_t part_b, size_t part_e) mutable {
            tmp = lop(r.discrete_domain(part_b, part_e), std::move(tmp));
          }); 
          
          // final reduce - tmp is guaranteed to have value
          // assert(tmp.has_value());
          std::lock_guard<std::mutex> lock(*mutex);
          init = gop(std::move(init), std::move(*tmp));
        });

        (++w == W || (curr_b += chunk_size) >= N) ? task() : rt.silent_async(task);
      }
    }
    // dynamic partitioner
    else {
      auto next = std::make_shared<std::atomic<size_t>>(0);
      
      for(size_t w=0; w<W;) {

        auto task = part([=, &init] () mutable {
          
          // temporary result so far
          std::optional<T> tmp;
          
          // loop reduce
          part.loop(N, W, *next, [=, &tmp](size_t part_b, size_t part_e) mutable {
            tmp = lop(r.discrete_domain(part_b, part_e), std::move(tmp));
          }); 
          
          // final reduce - need to check if the running total has value since
          // this is a dynamic scheduler; the worker may not actually acquire any work
          if(tmp) {
            std::lock_guard<std::mutex> lock(*mutex);
            init = gop(std::move(init), std::move(*tmp));
          }
        });
        (++w == W) ? task() : rt.silent_async(task);
      }
    }
  };
}

// ------------------------------------------------------------------------------------------------
// default reduction
// ------------------------------------------------------------------------------------------------

// Function: reduce
template <typename B, typename E, typename T, typename O, typename P>
Task FlowBuilder::reduce(B beg, E end, T& init, O bop, P part) {
  return emplace(make_reduce_task(beg, end, init, bop, part));
}

// ------------------------------------------------------------------------------------------------
// default transform and reduction
// ------------------------------------------------------------------------------------------------

// Function: transform_reduce
template <typename B, typename E, typename T, typename BOP, typename UOP, typename P,
  std::enable_if_t<is_partitioner_v<std::decay_t<P>>, void>*
>
Task FlowBuilder::transform_reduce(
  B beg, E end, T& init, BOP bop, UOP uop, P part
) {
  return emplace(make_transform_reduce_task(beg, end, init, bop, uop, part));
}

// Function: transform_reduce
template <
  typename B1, typename E1, typename B2, typename T, typename BOP_R, typename BOP_T, 
  typename P,
  std::enable_if_t<!is_partitioner_v<std::decay_t<BOP_T>>, void>*
>
Task FlowBuilder::transform_reduce(
  B1 beg1, E1 end1, B2 beg2, T& init, BOP_R bop_r, BOP_T bop_t, P part
) {
  return emplace(make_transform_reduce_task(beg1, end1, beg2, init, bop_r, bop_t, part));
}

// ------------------------------------------------------------------------------------------------
// default reduce_by_key
// ------------------------------------------------------------------------------------------------

// Function: make_index_reduce_task
template <typename R, typename T, typename L, typename G, typename P>
Task FlowBuilder::reduce_by_index(R range, T& init, L lop, G gop, P part) {
  return emplace(make_reduce_by_index_task(range, init, lop, gop, part));
}

}  // end of namespace tf -------------------------------------------------------------------------




