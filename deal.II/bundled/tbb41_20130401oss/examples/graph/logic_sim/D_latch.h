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

#ifndef __TBBexample_graph_logicsim_dlatch_H
#define __TBBexample_graph_logicsim_dlatch_H 1

#include "basics.h"

class D_latch {
    broadcast_node<signal_t> D_port;
    broadcast_node<signal_t> E_port;
    not_gate a_not;
    and_gate<two_input> first_and;
    and_gate<two_input> second_and;
    nor_gate<two_input> first_nor;
    nor_gate<two_input> second_nor;
    graph& my_graph;
 public:
    D_latch(graph& g) : my_graph(g), D_port(g), E_port(g), a_not(g), first_and(g), second_and(g), 
                        first_nor(g), second_nor(g) 
    {
        make_edge(D_port, a_not.get_in(0));
        make_edge(D_port, second_and.get_in(1));
        make_edge(E_port, first_and.get_in(1));
        make_edge(E_port, second_and.get_in(0));
        make_edge(a_not.get_out(), first_and.get_in(0));
        make_edge(first_and.get_out(), first_nor.get_in(0));
        make_edge(second_and.get_out(), second_nor.get_in(1));
        make_edge(first_nor.get_out(), second_nor.get_in(0));
        make_edge(second_nor.get_out(), first_nor.get_in(1));
    }
    ~D_latch() {}
    receiver<signal_t>& get_D() { return D_port; }
    receiver<signal_t>& get_E() { return E_port; }
    sender<signal_t>& get_Q() { return first_nor.get_out(); }
    sender<signal_t>& get_notQ() { return second_nor.get_out(); }
};

#endif /* __TBBexample_graph_logicsim_dlatch_H */
