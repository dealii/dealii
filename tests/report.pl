print <<'EOT'
<HTML>
<head>
<title>Regression tests</title>
</head>
<body>
<h1>Regression tests</h1>
<h2>Results</h2>
<table>
<tr><th>Date
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

	$results{$date}{$dir.':'.$name} = $result;
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

for ($i=1;$i<$next_index;$i++)
{
    printf "<th>%02d", $i;
}

# finally output a table of results
foreach $date (sort keys %results) {
    print "<tr><td>$date  ";
    foreach $name (sort keys %testcase) {
	$_ = $results{$date}{$name};
	if (/-/)
	{
	    print '<th> <font color="red">X</font>';
	} elsif (/\+/)
	{
	    print '<th> <font color="green">O</font>';
	} else {
	    print '<th>?';
	}
    }
    print "\n";
}

print << 'EOT'
</table>

<h2>Names of test programs</h2>
<table>
EOT
    ;


# output the list of test cases
foreach $name (sort keys %testcase) {
    print "<tr><td>$testcase{$name} <td>$name\n";
}

print <<'EOT'
</table>
</body>
</html>
EOT
    ;
