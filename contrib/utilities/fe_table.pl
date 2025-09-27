## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2005 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Author: Guido Kanschat

# Create a table of finite element capabilities out of
# tests/fe/fe_data_test.output

use strict;

print <<'EOF'
<html>
<head>
<title>Finite element capabilities</title>
</head>
<body>
<table border="1">
<tr>
<th rowspan="2">Element</th>
<th rowspan="2">degree</th>
<th rowspan="2">DoFs</th>
<th colspan="4">DoFs on</th>
<th rowspan="2">conforms</th>
<th rowspan="2">components</th>
<th colspan="4">Support points</th>
</tr>
<tr>
<th>V</th><th>L</th><th>Q</th><th>H</th>
<th>uc</th><th>uf</th><th>gc</th><th>gf</th>
</tr>
EOF
    ;

my @field = ('degree', 'dofs_per_cell',
	     'dofs_per_vertex', 'dofs_per_line', 'dofs_per_quad', 'dofs_per_hex',
	     'conformity', 'components',
	     'unit_support_points', 'unit_face_support_points',
	     'generalized_support_points', 'generalized_face_support_points');

my %fe;
my $hashref;

while(<>)
{
    if (/DEAL::fe_data.*:(.*)/)
    {
	$fe{$1} =  { 'set' => 't' } unless ($1 =~ m/FESystem/);
	$hashref = $fe{$1};
    }
    foreach my $entry (@field)
    {
	$hashref->{$entry} = $1 if /DEAL::$entry=(.*)/;
    }
}

foreach (sort keys %fe)
{
    print '<tr><td>',$_,"</td>\n";
    $hashref = $fe{$_};
    foreach (@field)
    {
	print "<td>", $hashref->{$_}, "</td>\n";
    }
    print "</tr>\n";
}

print <<'EOF'
</table>
</body>
</html>
EOF
    ;
