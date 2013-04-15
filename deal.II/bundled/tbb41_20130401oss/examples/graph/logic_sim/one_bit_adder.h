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

#ifndef __TBBexample_graph_logicsim_oba_H
#define __TBBexample_graph_logicsim_oba_H 1

#include "basics.h"

class one_bit_adder {
    broadcast_node<signal_t> A_port;
    broadcast_node<signal_t> B_port;
    broadcast_node<signal_t> CI_port;
    xor_gate<two_input> FirstXOR;
    xor_gate<two_input> SecondXOR;
    and_gate<two_input> FirstAND;
    and_gate<two_input> SecondAND;
    or_gate<two_input> FirstOR;
    graph& my_graph;
public:
    one_bit_adder(graph& g) : my_graph(g), A_port(g), B_port(g), CI_port(g), FirstXOR(g), 
                              SecondXOR(g), FirstAND(g), SecondAND(g), FirstOR(g) {
        make_connections();
    }
    one_bit_adder(const one_bit_adder& src) : 
        my_graph(src.my_graph), A_port(src.my_graph), B_port(src.my_graph), 
        CI_port(src.my_graph), FirstXOR(src.my_graph), SecondXOR(src.my_graph), 
        FirstAND(src.my_graph), SecondAND(src.my_graph), FirstOR(src.my_graph) 
    {
        make_connections();
    }
        
    ~one_bit_adder() {}
    receiver<signal_t>& get_A() { return A_port; }
    receiver<signal_t>& get_B() { return B_port; }
    receiver<signal_t>& get_CI() { return CI_port; }
    sender<signal_t>& get_out() {
        return SecondXOR.get_out();
    }
    sender<signal_t>& get_CO() {
        return FirstOR.get_out();
    }
private:
    void make_connections() {
        make_edge(A_port, FirstXOR.get_in(0));
        make_edge(A_port, FirstAND.get_in(0));
        make_edge(B_port, FirstXOR.get_in(1));
        make_edge(B_port, FirstAND.get_in(1));
        make_edge(CI_port, SecondXOR.get_in(1));
        make_edge(CI_port, SecondAND.get_in(1));
        make_edge(FirstXOR.get_out(), SecondXOR.get_in(0));
        make_edge(FirstXOR.get_out(), SecondAND.get_in(0));
        make_edge(SecondAND.get_out(), FirstOR.get_in(0));
        make_edge(FirstAND.get_out(), FirstOR.get_in(1));
    }
};

#endif /* __TBBexample_graph_logicsim_oba_H */
