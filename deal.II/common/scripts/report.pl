print <<'EOT'
<HTML>
<head>
    <link href="screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="author" content="The deal.II authors">
    <meta name="keywords" content="deal.II"></head>
<title>Regression tests</title>
</head>
<body>
<h1>Regression tests</h1>
<h2>Results</h2>
EOT
    ;


# read in list of test results
while (<>)
{
    $dir = $1 if (m/=====Report: ([^ ]+)/);
    if (/([\d-]+) ([\d:]*) ([+-]) (.*)/)
    {
	$date   = $1;
	$time   = $2;
	$result = $3;
	$name   = $4;

	$results{$date}{$dir.':'.$name} = '<img src="pictures/ok.gif" size="1">'
	    if ($result eq '+');
	$results{$date}{$dir.':'.$name} = '<img src="pictures/fail.gif" size="1">'
	    if ($result eq '-');
    }
}


# generate a list of test case names and assign a number
# in alphabetical order to them
foreach $date (keys %results) {
    foreach $name (keys %{ $results{$date} }) {
	$testcase{$name} = 0;
    }
}
$next_index = 1;
foreach $name (sort keys %testcase) {
    $testcase{$name} = $next_index++;
}


print <<'EOT'
<table>
<tr><th>Date
EOT
    ;

for ($i=1;$i<$next_index;$i++)
{
    printf "<th><small>%02d</small>", $i;
}
print "\n";


# finally output a table of results
foreach $date (sort {$b cmp $a} keys %results)
{
    # if this is not the first iteration, and if the month has changed,
    # then put in a break into the table to avoid overly long tables
    # which browsers take infinitely long to render
    $date =~ /(\d+)-(\d+)-\d+/;
    $this_year  = $1;
    $this_month = $2;
    if ($this_month != $old_month) {
	if (defined $old_month) {
	    print "</table>\n";

	    print <<'EOT'
		<table>
   	        <tr><th>Date
EOT
			    ;

	    for ($i=1;$i<$next_index;$i++)
	    {
		printf "<th><small>%02d</small>", $i;
	    }
	    print "\n";
	}
    }

    print "<tr><td>$date  ";
    foreach $name (sort keys %testcase)
    {
	$_ = $results{$date}{$name};
	print '<td> ', $_, '</td>';
    }
    print "</tr></td>\n";

    # store old month name for the next iteration of the loop
    $old_month = $this_month;
}

print << 'EOT'
</table>

<h2>Names of test programs</h2>
<table>
EOT
    ;


# output the list of test cases. always put four in a row
$col = 0;
foreach $name (sort keys %testcase) {
    if ($col == 0) {
	print "<tr>\n";
    }
    print "   <td>$testcase{$name} <td>$name\n";

    # next column. if at end, wrap around
    $col = ($col+1)%4;
}

print <<'EOT'
</table>
</body>
</html>
EOT
    ;


