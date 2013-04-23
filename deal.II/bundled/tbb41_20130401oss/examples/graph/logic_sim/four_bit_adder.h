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

#ifndef __TBBexample_graph_logicsim_fba_H
#define __TBBexample_graph_logicsim_fba_H 1

#include "one_bit_adder.h"

class four_bit_adder {
    graph& my_graph;
    std::vector<one_bit_adder> four_adders; 
 public:
    four_bit_adder(graph& g) : my_graph(g), four_adders(4, one_bit_adder(g)) {
        make_connections();
    }
    four_bit_adder(const four_bit_adder& src) : 
        my_graph(src.my_graph), four_adders(4, one_bit_adder(src.my_graph)) 
    {
        make_connections();
    }
    ~four_bit_adder() {}
    receiver<signal_t>& get_A(size_t bit) {
        return four_adders[bit].get_A();
    }
    receiver<signal_t>& get_B(size_t bit) {
        return four_adders[bit].get_B();
    }
    receiver<signal_t>& get_CI() {
        return four_adders[0].get_CI();
    }
    sender<signal_t>& get_out(size_t bit) {
        return four_adders[bit].get_out();
    }
    sender<signal_t>& get_CO() {
        return four_adders[3].get_CO();
    }
private:
    void make_connections() {
        make_edge(four_adders[0].get_CO(), four_adders[1].get_CI());
        make_edge(four_adders[1].get_CO(), four_adders[2].get_CI());
        make_edge(four_adders[2].get_CO(), four_adders[3].get_CI());
    }
};

#endif /* __TBBexample_graph_logicsim_fba_H */
