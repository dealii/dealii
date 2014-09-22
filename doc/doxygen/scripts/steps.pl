## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

use strict;

my $tutorial_file = shift;
open TUTORIAL, "<$tutorial_file";

# Print the first part of tutorial.h.in
while (my $line = <TUTORIAL>)
{
  last if($line =~ m/\@\@MAP\@\@/);
  print $line;
}

# List of additional node attributes to highlight purpose and state of the example
my %style = (
 "basic"          => ',height=.8,width=.8,shape="octagon",fillcolor="green"',
 "techniques"     => ',height=.35,width=.35,fillcolor="orange"',
 "fluids"         => ',height=.25,width=.25,fillcolor="yellow"',
 "solids"         => ',height=.25,width=.25,fillcolor="lightblue"',
 "time dependent" => ',height=.25,width=.25,fillcolor="blue"',
 "unfinished"     => ',height=.25,width=.25,style="dashed"'
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

my $step;
foreach $step (@ARGV)
{
    my $number = $step;
    $number =~ s/^.*-//;

    # read first line of tooltip file
    open TF, "$step/doc/tooltip"
	or die "Can't open tooltip file $step/doc/tooltip";
    my $tooltip = <TF>;
    close TF;
    chop $tooltip;

    printf "Step$number [label=\"$number\", URL=\"\\ref step_$number\", tooltip=\"$tooltip\"";


    # read first line of 'kind' file
    open KF, "$step/doc/kind"
	or die "Can't open kind file $step/doc/kind";
    my $kind = <KF>;
    close KF;
    chop $kind;

    die "Unknown kind '$kind' in file $step/doc/kind" if (! defined $style{$kind});
    print "$style{$kind}";

    print "];\n";
}

# Print all edges
# Keep sorted by second node on edge!

my $step;
foreach $step (@ARGV)
{
    my $number = $step;
    $number =~ s/^.*-//;

    # read first line of dependency file
    open BF, "$step/doc/builds-on"
	or die "Can't open builds-on file $step/doc/builds-on";
    my $buildson = <BF>;
    close BF;
    chop $buildson;

    my $source;
    foreach $source (split ' ', $buildson) {
	$source =~ s/step-/Step/g;
	print "$source -> Step$number\n";
    }
}

print "}\n";

# Print the rest of tutorial.h.in
while (my $line = <TUTORIAL>)
{
  print $line;
}
close TUTORIAL;
