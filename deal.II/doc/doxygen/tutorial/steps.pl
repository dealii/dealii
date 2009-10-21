######################################################################
# $Id$
######################################################################
#
# Copyright (c) the deal.II authors 2009
#
######################################################################

use strict;

my $laststep = 38;

my $essential = ',fillcolor="red"';
my $technique = ',fillcolor="orange"';
my $application = ',fillcolor="yellow"';
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
    17 => $technique,
    35 => $application,
    38 => $unfinished
    );

# Print a preamble setting common attributes

print << 'EOT'
digraph G
{
  edge [fontname="FreeSans",fontsize="10",labelfontname="FreeSans",labelfontsize="10",color="black",style="solid"];
  node [fontname="FreeSans",fontsize="10",shape=record,height=0.2,width=0.4,color="black",fillcolor="white",style="filled"];
  
EOT
    ;

# print all nodes

for (my $i=1; $i<=$laststep;++$i)
{
    printf 'Step%02d [label="Step %d", URL="step_%d.html"', $i, $i, $i;
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
Step06 -> Step10;
Step10 -> Step11;
Step06 -> Step12;
Step05 -> Step16;
Step12 -> Step38;
Step38 -> Step12;
Step06 -> Step13;
Step06 -> Step14;
Step06 -> Step15;
Step06 -> Step16;
Step06 -> Step17;
Step08 -> Step18;
Step06 -> Step19;
}

EOT
    ;

