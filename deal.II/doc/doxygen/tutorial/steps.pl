######################################################################
# $Id$
######################################################################
#
# Copyright (c) the deal.II authors 2009
#
######################################################################

use strict;

my @steps = (1,2,3,4,5,6,7,8,9,
	     10,11,12,13,14,15,16,17,18,19,
	     20,21,22,23,24,25,   27,28,29,
	     30,31,   33,34,35,36,      39);
;

# List of additional node attributes to highlight purpose and state of the example
my %style = (
 "basic" => ',height=.7,width=.7,shape="octagon",fillcolor="green"',
 "techniques" => ',fillcolor="orange"',
 "fluids" => ',fillcolor="yellow"',
 "solids" => ',fillcolor="lightblue"',
 "time dependent" => ',fillcolor="blue"',
 "unfinished" => ',style="dashed"'
    );


# Print a preamble setting common attributes

print << 'EOT'
digraph StepsMap
{
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="black",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
EOT
    ;

# print all nodes

foreach (@steps)
{
    # read first line of tooltip file
    open TF, "../../../examples/step-$_/doc/tooltip"
	or die "Can't open tooltip file step-$_/doc/tooltip";
    my $tooltip = <TF>;
    close TF;
    chop $tooltip;

    printf "Step%02d [label=\"$_\", URL=\"../deal.II/step_$_.html\", tooltip=\"$tooltip\"", $_;

    open KF, "../../../examples/step-$_/doc/kind"
	or die "Can't open kind file step-$_/doc/kind";
    my $kind = <KF>;
    close KF;
    chop $kind;

    die "Unknown kind '$kind' in file step-$_/doc/kind" if (! defined $style{$kind});
    print "$style{$kind}";

    print "];\n";
}

# Print all edges
# Keep sorted by second node on edge!

print << 'EOT'

Step01 -> Step02;
Step02 -> Step03;
Step03 -> Step04;
Step04 -> Step05;
Step05 -> Step06;

Step06 -> Step07;
Step06 -> Step08;
Step06 -> Step09;

Step04 -> Step10;
Step10 -> Step11;

Step07 -> Step12;

Step06 -> Step16;

Step12 -> Step33;

Step06 -> Step13;
Step13 -> Step14;

Step04 -> Step15;

Step08 -> Step17;
Step17 -> Step18;

Step04 -> Step20;
Step20 -> Step21;

Step06 -> Step22;
Step22 -> Step31;

Step04 -> Step23;
Step23 -> Step24;
Step24 -> Step25;

Step06 -> Step27;

Step06 -> Step28;

Step06 -> Step39;

Step04 -> Step29;

Step12 -> Step30;
Step12 -> Step39;

Step04 -> Step34;

Step21 -> Step22;
Step22 -> Step35;

Step04 -> Step36;

Step39 -> Step12;
}

EOT
    ;

