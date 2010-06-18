######################################################################
# $Id$
######################################################################
#
# Copyright (c) the deal.II authors 2009
#
######################################################################

use strict;

my $laststep = 45;

my $essential = ',fillcolor="red"';
my $technique = ',fillcolor="orange"';
my $fluidapplication = ',fillcolor="yellow"';
my $solidsapplication = ',fillcolor="lightblue"';
my $timeapplication = ',fillcolor="blue"';
my $unfinished = ',style="dashed"';

# List of additional node attributes to highlight purpose and state of the example

my %attribute = (
    1 => $essential,
    2 => $essential,
    3 => $essential,
    4 => $essential,
    5 => $essential,
    6 => $essential,

    7 => $technique,
    8 => $technique,
    9 => $technique,
    10 => $technique,
    11 => $technique,
    12 => $technique,
    13 => $technique,
    14 => $technique,
    15 => $technique,
    16 => $technique,
    19 => $technique,
    27 => $technique,
    29 => $technique,
    30 => $technique,
    36 => $technique,
    37 => $technique,
    39 => $technique,
    45 => $technique,

    17 => $solidsapplication,
    18 => $solidsapplication,

    20 => $fluidapplication,
    21 => $fluidapplication,
    22 => $fluidapplication,
    31 => $fluidapplication,
    32 => $fluidapplication,
    33 => $fluidapplication,
    34 => $fluidapplication,
    35 => $fluidapplication,

    23 => $timeapplication,
    24 => $timeapplication,
    25 => $timeapplication,
    );

# Print a preamble setting common attributes

print << 'EOT'
digraph G
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
        shape=record,
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
EOT
    ;

# print all nodes

for (my $i=1; $i<=$laststep;++$i)
{
    printf 'Step%02d [label="%d", URL="step_%d.html", tooltip="@step%d@"', $i, $i, $i, $i;
    print $attribute{$i};
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
Step16 -> Step37;
Step05 -> Step37;

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
Step31 -> Step32;

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
Step03 -> Step45;

Step39 -> Step12;
}

EOT
    ;

