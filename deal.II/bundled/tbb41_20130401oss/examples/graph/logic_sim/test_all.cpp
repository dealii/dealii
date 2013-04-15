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

#if _MSC_VER
#pragma warning (disable: 4503) // Suppress "decorated name length exceeded, name was truncated" warning
#endif

#include "basics.h"
#include "one_bit_adder.h"
#include "four_bit_adder.h"
#include "D_latch.h"
#include <cassert>

// User-specified globals with default values
bool verbose = false;            // prints bin details and other diagnostics to screen
bool silent = false;             // suppress all output except for time

int get_default_num_threads() {
    static int threads = 0;
    if (threads == 0)
        threads = tbb::task_scheduler_init::default_num_threads();
    return threads;
}

int main(int argc, char *argv[]) {
    try {
        utility::thread_number_range threads(get_default_num_threads);
        utility::parse_cli_arguments(argc, argv,
                                     utility::cli_argument_pack()
                                     //"-h" option for displaying help is present implicitly
                                     .positional_arg(threads,"#threads",utility::thread_number_range_desc)
                                     .arg(verbose,"verbose","   print diagnostic output to screen")
                                     .arg(silent,"silent","    limits output to timing info; overrides verbose")
        );

        if (silent) verbose = false;  // make silent override verbose

        tick_count start = tick_count::now();
        for(int p = threads.first; p <= threads.last; p = threads.step(p)) {
            task_scheduler_init init(p);
            if (!silent)  cout << "graph test running on " << p << " threads.\n";
            
            graph g;

            { // test buffer: 0, 1
                buffer b(g);
                toggle input(g);
                led output(g, "OUTPUT", false); // false means we will explicitly call display to see LED
                
                make_edge(input.get_out(), b.get_in(0));
                make_edge(b.get_out(), output.get_in());
                
                if (!silent) printf("Testing buffer...\n");
                input.activate(); // 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input.flip(); // 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
            }

            { // test not_gate: 0, 1
                not_gate n(g);
                toggle input(g);
                led output(g, "OUTPUT", false);
                
                make_edge(input.get_out(), n.get_in(0));
                make_edge(n.get_out(), output.get_in());
                
                if (!silent) printf("Testing not_gate...\n");
                input.activate(); // 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input.flip(); // 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
            }

            { // test two-input and_gate: 00, 01, 10, 11
                and_gate<two_input> a(g);
                toggle input0(g);
                toggle input1(g);
                led output(g, "OUTPUT", false);
                
                make_edge(input0.get_out(), a.get_in(0));
                make_edge(input1.get_out(), a.get_in(1));
                make_edge(a.get_out(), output.get_in());
                
                if (!silent) printf("Testing and_gate...\n");
                input1.activate();  input0.activate();  // 0 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input0.flip();  // 0 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input1.flip(); input0.flip();  // 1 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input0.flip();  // 1 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
            }

            { // test three-input or_gate: 000, 001, 010, 100, 011, 101, 110, 111
                or_gate<three_input> o(g);
                toggle input0(g);
                toggle input1(g);
                toggle input2(g);
                led output(g, "OUTPUT", false);
                
                make_edge(input0.get_out(), o.get_in(0));
                make_edge(input1.get_out(), o.get_in(1));
                make_edge(input2.get_out(), o.get_in(2));
                make_edge(o.get_out(), output.get_in());
                
                if (!silent) printf("Testing or_gate...\n");
                input2.activate();  input1.activate();  input0.activate();  // 0 0 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input0.flip();  // 0 0 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input1.flip(); input0.flip();  // 0 1 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input2.flip();  input1.flip();  // 1 0 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input2.flip();  input1.flip();  input0.flip();  // 0 1 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input2.flip();  input1.flip();  // 1 0 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input1.flip();  input0.flip();  // 1 1 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input0.flip();  // 1 1 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
            }

            { // test two-input xor_gate: 00, 01, 10, 11
                xor_gate<two_input> x(g);
                toggle input0(g);
                toggle input1(g);
                led output(g, "OUTPUT", false);
                
                make_edge(input0.get_out(), x.get_in(0));
                make_edge(input1.get_out(), x.get_in(1));
                make_edge(x.get_out(), output.get_in());
                
                if (!silent) printf("Testing xor_gate...\n");
                input1.activate();  input0.activate();  // 0 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input0.flip();  // 0 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input1.flip();  input0.flip();  // 1 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input0.flip();  // 1 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
            }


            { // test two-input nor_gate: 00, 01, 10, 11
                nor_gate<two_input> n(g);
                toggle input0(g);
                toggle input1(g);
                led output(g, "OUTPUT", false);
                
                make_edge(input0.get_out(), n.get_in(0));
                make_edge(input1.get_out(), n.get_in(1));
                make_edge(n.get_out(), output.get_in());
                
                if (!silent) printf("Testing nor_gate...\n");
                input1.activate();  input0.activate();  // 0 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == high);
                input0.flip();  // 0 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input1.flip();  input0.flip();  // 1 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
                input0.flip();  // 1 1
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == low);
            }

            { // test steady_signal and digit
                steady_signal input0(g, high);
                steady_signal input1(g, low);
                and_gate<two_input> a(g);
                or_gate<two_input> o(g);
                xor_gate<two_input> x(g);
                nor_gate<two_input> n(g);
                digit output(g, "OUTPUT", false);
                
                make_edge(input0.get_out(), a.get_in(0));
                make_edge(input1.get_out(), a.get_in(1));
                make_edge(a.get_out(), output.get_in(0));

                make_edge(input0.get_out(), o.get_in(0));
                make_edge(input1.get_out(), o.get_in(1));
                make_edge(o.get_out(), output.get_in(1));

                make_edge(input0.get_out(), x.get_in(0));
                make_edge(input1.get_out(), x.get_in(1));
                make_edge(x.get_out(), output.get_in(2));

                make_edge(input0.get_out(), n.get_in(0));
                make_edge(input1.get_out(), n.get_in(1));
                make_edge(n.get_out(), output.get_in(3));
                
                if (!silent) printf("Testing steady_signal...\n");
                input0.activate();  // 1
                input1.activate();  // 0
                g.wait_for_all();
                if (!silent) output.display();
                assert(output.get_value() == 6);
            }

            { // test push_button
                push_button p(g);
                buffer b(g);
                led output(g, "OUTPUT", !silent); // true means print all LED state changes

                make_edge(p.get_out(), b.get_in(0));
                make_edge(b.get_out(), output.get_in());

                if (!silent) printf("Testing push_button...\n");
                p.press();
                p.release();
                p.press();
                p.release();
                g.wait_for_all();
            }

            { // test one_bit_adder
                one_bit_adder my_adder(g);
                toggle A(g);
                toggle B(g);
                toggle CarryIN(g);
                led Sum(g, "SUM");
                led CarryOUT(g, "CarryOUT");
                
                make_edge(A.get_out(), my_adder.get_A());
                make_edge(B.get_out(), my_adder.get_B());
                make_edge(CarryIN.get_out(), my_adder.get_CI());
                make_edge(my_adder.get_out(), Sum.get_in());
                make_edge(my_adder.get_CO(), CarryOUT.get_in());
                
                A.activate();
                B.activate();
                CarryIN.activate();
                
                if (!silent) printf("A on\n");
                A.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == high) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("A off\n");
                A.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("B on\n");
                B.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == high) && (CarryOUT.get_value() == low));
                if (!silent) printf("B off\n");
                B.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("CarryIN on\n");
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == high) && (CarryOUT.get_value() == low));
                if (!silent) printf("CarryIN off\n");
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("A&B on\n");
                A.flip();
                B.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == high));
                if (!silent) printf("A&B off\n");
                A.flip();
                B.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("A&CarryIN on\n");
                A.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == high));
                if (!silent) printf("A&CarryIN off\n");
                A.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("B&CarryIN on\n");
                B.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == high));
                if (!silent) printf("B&CarryIN off\n");
                B.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("A&B&CarryIN on\n");
                A.flip();
                B.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == high) && (CarryOUT.get_value() == high));
                if (!silent) printf("A&B&CarryIN off\n");
                A.flip();
                B.flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == low) && (CarryOUT.get_value() == low));
            }

            { // test four_bit_adder
                four_bit_adder four_adder(g);
                std::vector<toggle> A(4, toggle(g));
                std::vector<toggle> B(4, toggle(g));
                toggle CarryIN(g);
                digit Sum(g, "SUM");
                led CarryOUT(g, "CarryOUT");
                
                for (int i=0; i<4; ++i) {
                    make_edge(A[i].get_out(), four_adder.get_A(i));
                    make_edge(B[i].get_out(), four_adder.get_B(i));
                    make_edge(four_adder.get_out(i), Sum.get_in(i));
                }
                make_edge(CarryIN.get_out(), four_adder.get_CI());
                make_edge(four_adder.get_CO(), CarryOUT.get_in());
                
                // Activate all switches at low state
                for (int i=0; i<4; ++i) {
                    A[i].activate();
                    B[i].activate();
                }
                CarryIN.activate();
                
                if (!silent) printf("1+0\n");
                A[0].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 1) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("0+1\n");
                A[0].flip();
                B[0].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 1) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("3+4\n");
                A[0].flip();
                A[1].flip();
                B[0].flip();
                B[2].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 7) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("6+1\n");
                A[0].flip();
                A[2].flip();
                B[0].flip();
                B[2].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 7) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("0+0+carry\n");
                A[1].flip();
                A[2].flip();
                B[0].flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 1) && (CarryOUT.get_value() == low));
                
                if (!silent) printf("15+15+carry\n");
                A[0].flip();
                A[1].flip();
                A[2].flip();
                A[3].flip();
                B[0].flip();
                B[1].flip();
                B[2].flip();
                B[3].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 0xf) && (CarryOUT.get_value() == high));
                
                if (!silent) printf("8+8\n");
                A[0].flip();
                A[1].flip();
                A[2].flip();
                B[0].flip();
                B[1].flip();
                B[2].flip();
                CarryIN.flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 0) && (CarryOUT.get_value() == high));
                
                if (!silent) printf("0+0\n");
                A[3].flip();
                B[3].flip();
                g.wait_for_all();
                if (!silent) Sum.display();
                if (!silent) CarryOUT.display();
                assert((Sum.get_value() == 0) && (CarryOUT.get_value() == low));
            }

            { // test D_latch
                D_latch my_d_latch(g);
                toggle D(g);
                pulse E(g, 500, 4); // clock changes every 500ms; stops after 4 changes
                led Q(g, " Q", verbose); // if true, LEDs print at every state change
                led notQ(g, "~Q", verbose);
                
                make_edge(D.get_out(), my_d_latch.get_D());
                make_edge(E.get_out(), my_d_latch.get_E());
                make_edge(my_d_latch.get_Q(), Q.get_in());
                make_edge(my_d_latch.get_notQ(), notQ.get_in());
                
                D.activate();
                
                if (!silent) printf("Toggling D\n");
                E.activate();
                D.flip();
                g.wait_for_all();
                if (!silent && !verbose) { Q.display(); notQ.display(); }
                assert((Q.get_value() == high) && (notQ.get_value() == low));
                E.reset();
                
                if (!silent) printf("Toggling D\n");
                E.activate();
                D.flip();
                g.wait_for_all();
                if (!silent && !verbose) { Q.display(); notQ.display(); }
                assert((Q.get_value() == low) && (notQ.get_value() == high));
                E.reset();
                
                if (!silent) printf("Toggling D\n");
                E.activate();
                D.flip();
                g.wait_for_all();
                if (!silent && !verbose) { Q.display(); notQ.display(); }
                assert((Q.get_value() == high) && (notQ.get_value() == low));
                E.reset();
                
                if (!silent) printf("Toggling D\n");
                E.activate();
                D.flip();
                g.wait_for_all();
                if (!silent && !verbose) { Q.display(); notQ.display(); }
                assert((Q.get_value() == low) && (notQ.get_value() == high));
                E.reset();
                
                if (!silent) printf("Toggling D\n");
                E.activate();
                D.flip();
                g.wait_for_all();
                if (!silent && !verbose) { Q.display(); notQ.display(); }
                assert((Q.get_value() == high) && (notQ.get_value() == low));
            }
        }
        utility::report_elapsed_time((tbb::tick_count::now() - start).seconds());
        return 0;
    } catch(std::exception& e) {
        cerr<<"error occurred. error text is :\"" <<e.what()<<"\"\n";
        return 1;
    }
}
