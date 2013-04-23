/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

#ifndef __TBB__flow_graph_or_impl_H
#define __TBB__flow_graph_or_impl_H

#ifndef __TBB_flow_graph_H
#error Do not #include this internal file directly; use public TBB headers instead.
#endif

#if TBB_PREVIEW_GRAPH_NODES
#include "tbb/internal/_flow_graph_types_impl.h"

namespace internal {

    // Output of the or_node is a struct containing a tbb::flow::tuple, and will be of
    // the form
    //
    //  struct {
    //     size_t indx;
    //     tuple_types result;
    //  };
    //
    //  where the value of indx will indicate which result was put to the
    //  successor. So if oval is the output to the successor, indx == 0 
    //  means tbb::flow::get<0>(oval.result) is the output, and so on.
    //
    //  tuple_types is the tuple that specified the possible outputs (and
    //  the corresponding inputs to the or_node.)
    //
    //  the types of each element are represented by tuple_types, a typedef
    //  in the or_node.  So the 2nd type in the union that is the
    //  output type for an or_node OrType is
    //
    //      tbb::flow::tuple_element<1,OrType::tuple_types>::type

    // the struct has an OutputTuple default constructed, with element index assigned
    // the actual output value.
    template<typename OutputTuple>
    struct or_output_type {
        typedef OutputTuple tuple_types;
        struct type {
            size_t indx;
            OutputTuple result;

// The LLVM libc++ that ships with OS X* 10.7 has a bug in tuple that disables
// the copy assignment operator (LLVM bug #11921).
//TODO: introduce according broken macro.
//it can not be done right now, as tbb_config.h does not allowed to include other headers,
//and without this it is not possible to detect libc++ version, as compiler version for clang
//is vendor specific
#ifdef _LIBCPP_TUPLE
            type &operator=(type const &x) {
                indx = x.indx;
                result = const_cast<OutputTuple&>(x.result);
                return *this;
            }
#endif
        };
    };

    template<typename TupleTypes,int N>
    struct or_item_helper {
        template<typename OutputType>
        static inline void create_output_value(OutputType &o, void *v) {
            o.indx = N;
            tbb::flow::get<N>(o.result) = *(reinterpret_cast<typename tbb::flow::tuple_element<N, TupleTypes>::type *>(v));
        }
    };

    template<typename TupleTypes,int N>
    struct or_helper {
        template<typename OutputType>
        static inline void create_output(OutputType &o, size_t i, void* v) {
            if(i == N-1) {
                or_item_helper<TupleTypes,N-1>::create_output_value(o,v);
            }
            else
                or_helper<TupleTypes,N-1>::create_output(o,i,v);
        }
        template<typename PortTuple, typename PutBase>
        static inline void set_or_node_pointer(PortTuple &my_input, PutBase *p) {
            tbb::flow::get<N-1>(my_input).set_up(p, N-1);
            or_helper<TupleTypes,N-1>::set_or_node_pointer(my_input, p);
        }
    };

    template<typename TupleTypes>
    struct or_helper<TupleTypes,1> {
        template<typename OutputType>
        static inline void create_output(OutputType &o, size_t i, void* v) {
            if(i == 0) {
                or_item_helper<TupleTypes,0>::create_output_value(o,v);
            }
        }
        template<typename PortTuple, typename PutBase>
        static inline void set_or_node_pointer(PortTuple &my_input, PutBase *p) {
            tbb::flow::get<0>(my_input).set_up(p, 0);
        }
    };

    struct put_base {
        // virtual bool try_put_with_index(size_t index, void *v) = 0;
        virtual task * try_put_task_with_index(size_t index, void *v) = 0;
        virtual ~put_base() { }
    };

    template<typename T>
    class or_input_port : public receiver<T> {
    private:
        size_t my_index;
        put_base *my_or_node;
    public:
        void set_up(put_base *p, size_t i) { my_index = i; my_or_node = p; }
    protected:
        template< typename R, typename B > friend class run_and_put_task;
        template<typename X, typename Y> friend class internal::broadcast_cache;
        template<typename X, typename Y> friend class internal::round_robin_cache;
        task *try_put_task(const T &v) {
            return my_or_node->try_put_task_with_index(my_index, reinterpret_cast<void *>(const_cast<T*>(&v)));
        }
        /*override*/void reset_receiver() {}
    };

    template<typename InputTuple, typename OutputType, typename StructTypes>
    class or_node_FE : public put_base {
    public:
        static const int N = tbb::flow::tuple_size<InputTuple>::value;
        typedef OutputType output_type;
        typedef InputTuple input_type;

        or_node_FE( ) {
            or_helper<StructTypes,N>::set_or_node_pointer(my_inputs, this);
        }

        input_type &input_ports() { return my_inputs; }
    protected:
        input_type my_inputs;
    };

    //! or_node_base
    template<typename InputTuple, typename OutputType, typename StructTypes>
    class or_node_base : public graph_node, public or_node_FE<InputTuple, OutputType,StructTypes>,
                           public sender<OutputType> {
    protected:
       using graph_node::my_graph;
    public:
        static const size_t N = tbb::flow::tuple_size<InputTuple>::value;
        typedef OutputType output_type;
        typedef StructTypes tuple_types;
        typedef receiver<output_type> successor_type;
        typedef or_node_FE<InputTuple, output_type,StructTypes> input_ports_type;

    private:
        // ----------- Aggregator ------------
        enum op_type { reg_succ, rem_succ, try__put_task };
        enum op_stat {WAIT=0, SUCCEEDED, FAILED};
        typedef or_node_base<InputTuple,output_type,StructTypes> my_class;

        class or_node_base_operation : public aggregated_operation<or_node_base_operation> {
        public:
            char type;
            size_t indx;
            union {
                void *my_arg;
                successor_type *my_succ;
                task *bypass_t;
            };
            or_node_base_operation(size_t i, const void* e, op_type t) :
                type(char(t)), indx(i), my_arg(const_cast<void *>(e)) {}
            or_node_base_operation(const successor_type &s, op_type t) : type(char(t)), 
                my_succ(const_cast<successor_type *>(&s)) {}
            or_node_base_operation(op_type t) : type(char(t)) {}
        };

        typedef internal::aggregating_functor<my_class, or_node_base_operation> my_handler;
        friend class internal::aggregating_functor<my_class, or_node_base_operation>;
        aggregator<my_handler, or_node_base_operation> my_aggregator;

        void handle_operations(or_node_base_operation* op_list) {
            or_node_base_operation *current;
            while(op_list) {
                current = op_list;
                op_list = op_list->next;
                switch(current->type) {

                case reg_succ:
                    my_successors.register_successor(*(current->my_succ));
                    __TBB_store_with_release(current->status, SUCCEEDED);
                    break;

                case rem_succ:
                    my_successors.remove_successor(*(current->my_succ));
                    __TBB_store_with_release(current->status, SUCCEEDED);
                    break;
                case try__put_task: {
                        output_type oo;
                        or_helper<tuple_types,N>::create_output(oo, current->indx, current->my_arg);
                        current->bypass_t = my_successors.try_put_task(oo);
                        __TBB_store_with_release(current->status, SUCCEEDED);  // return of try_put_task actual return value
                    }
                    break;
                }
            }
        }
        // ---------- end aggregator -----------
    public:
        or_node_base(graph& g) : graph_node(g), input_ports_type() {
            my_successors.set_owner(this);
            my_aggregator.initialize_handler(my_handler(this));
        }

        or_node_base(const or_node_base& other) : graph_node(other.my_graph), input_ports_type(), sender<output_type>() {
            my_successors.set_owner(this);
            my_aggregator.initialize_handler(my_handler(this));
        }

        bool register_successor(successor_type &r) {
            or_node_base_operation op_data(r, reg_succ);
            my_aggregator.execute(&op_data);
            return op_data.status == SUCCEEDED;
        }

        bool remove_successor( successor_type &r) {
            or_node_base_operation op_data(r, rem_succ);
            my_aggregator.execute(&op_data);
            return op_data.status == SUCCEEDED;
        }

        task * try_put_task_with_index(size_t indx, void *v) {
            or_node_base_operation op_data(indx, v, try__put_task);
            my_aggregator.execute(&op_data);
            return op_data.bypass_t;
        }

    protected:
        /*override*/void reset() {}

    private:
        broadcast_cache<output_type, null_rw_mutex> my_successors;
    };

    // type generators
    template<typename OutputTuple>
    struct or_types {
        static const int N = tbb::flow::tuple_size<OutputTuple>::value;
        typedef typename wrap_tuple_elements<N,or_input_port,OutputTuple>::type input_ports_type;
        typedef typename or_output_type<OutputTuple>::type output_type;
        typedef internal::or_node_FE<input_ports_type,output_type,OutputTuple> or_FE_type;
        typedef internal::or_node_base<input_ports_type, output_type, OutputTuple> or_base_type;
    };

    template<class OutputTuple>
    class unfolded_or_node : public or_types<OutputTuple>::or_base_type {
    public:
        typedef typename or_types<OutputTuple>::input_ports_type input_ports_type;
        typedef OutputTuple tuple_types;
        typedef typename or_types<OutputTuple>::output_type output_type;
    private:
        typedef typename or_types<OutputTuple>::or_base_type base_type;
    public:
        unfolded_or_node(graph& g) : base_type(g) {}
        unfolded_or_node(const unfolded_or_node &other) : base_type(other) {}
    };


} /* namespace internal */
#endif  // TBB_PREVIEW_GRAPH_NODES

#endif  /* __TBB__flow_graph_or_impl_H */
