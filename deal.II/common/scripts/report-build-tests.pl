############################################################
# first read in list of test results
while (<>)
{
    if (/([\d-]+) ([\d:]*) ([^\s]+) (.*)/)
    {
	$date   = $1;
	$time   = $2;
	$name   = $3;
	$result = $4;

	$total_testcases{$date}++;

	if ($result eq '+') {
	    $results{$date}{$name} 
	    = '<img src="pictures/ok.gif" size="1">';
	}
	else
	{
	    # TODO: the result is actually the name of a file with logs
	    # in it. Use that and cross-link to that file so that
	    # one can see what has gone wrong
	    $results{$date}{$name}
	    = '<img src="pictures/fail.gif" size="1">';
	    $failed_testcases{$date}++;
	};
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




###########################################################
# then generate output for the three frames and the panels for 
# the different months

# first get last active month for the default panel
foreach $date (sort {$b cmp $a} keys %results)
{
    $date =~ /(\d+)-(\d+)-\d+/;
    $default_year  = $1;
    $default_month = $2;
    last;
}



open REPORT_FILE, ">builds_report.html";

print REPORT_FILE <<"EOT"
<HTML>
<head>
    <link href="screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="author" content="The deal.II authors">
    <meta name="keywords" content="deal.II"></head>
<title>Build tests</title>
</head>

<frameset rows="120,*,120" border=2>
<frame name="report_head"    src="builds_report_head.html"    frameborder="yes">
<frame name="report_results" src="builds_report_${default_year}_${default_month}.html" frameborder="yes">
<frame name="report_names"   src="builds_report_names.html"   frameborder="yes">
</frameset>
</html>
EOT
    ;


open HEAD_FILE, ">builds_report_head.html";
print HEAD_FILE << 'EOT'
<HTML>
<head>
    <link href="screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="author" content="The deal.II authors">
    <meta name="keywords" content="deal.II"></head>
<title>Build tests head</title>
</head>
<body>
<h1 align="center">Build tests</h1>
<p align="center">
Select results for one of the following months:<br>
EOT
    ;

open NAMES_FILE, ">builds_report_names.html";
print NAMES_FILE << 'EOT'
<HTML>
<head>
    <link href="screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="author" content="The deal.II authors">
    <meta name="keywords" content="deal.II"></head>
<title>Build tests head</title>
</head>
<body>
EOT
    ;




# finally output a table of results
$number_of_present_month = 0;
foreach $date (sort {$b cmp $a} keys %results)
{
    # if the month has changed (or if this is the first month we deal
    # with), open a new file and write the file and table heads. also write
    # a note into the head panel
    #
    # if this is not the first month, then put in a break into the table
    # to avoid overly long tables which browsers take infinitely
    # long to render
    $date =~ /(\d+)-(\d+)-\d+/;
    $this_year  = $1;
    $this_month = $2;
    if ($this_month != $old_month) {
	$number_of_present_month++;

	if (defined $old_month) {
	    print TABLE_FILE "</table>\n";
	    print TABLE_FILE "</body>\n</html>\n";
	    close TABLE_FILE;
	}
	$file = "builds_report_${this_year}_${this_month}.html";
	use Cwd;
	$dir = cwd();
	open TABLE_FILE, ">$file" or die "Can't open output file $file in $dir\n";

print TABLE_FILE <<"EOT"
<HTML>
<head>
    <link href="screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="author" content="The deal.II authors">
    <meta name="keywords" content="deal.II"></head>
<title>Build tests for $year/$date</title>
</head>
<body>

<h3 align="center">Results for $this_year/$this_month</h3>
<table align="center">
<tr><th>Date <th> Fail
EOT
    ;

	for ($i=1;$i<$next_index;$i++)
	{
	    printf TABLE_FILE "<th><small>%02d</small>", $i;
	}
	print TABLE_FILE "\n";

	
	# only write up to 8 months in a row into the month
	# listing. if more, insert a linebreak
	if (($number_of_present_month % 8 == 0) &&
	    ($number_of_present_month != 0)) 
	{
	    print HEAD_FILE "<br>\n";
	}
	# now write link to month file
	print HEAD_FILE "<a href=\"$file\" target=\"report_results\">";
	print HEAD_FILE "${this_year}/${this_month}</a>&nbsp;&nbsp;\n";
    }

    print TABLE_FILE "<tr><td>$date  ";

    $failed_testcases{$date} = 0 if (!defined $failed_testcases{$date});

    print TABLE_FILE "<td><b style=\"color:blue;\">" if ($failed_testcases{$date} == 0);
    print TABLE_FILE "<td><b style=\"color:red;\">" if ($failed_testcases{$date} != 0);

    print TABLE_FILE "$failed_testcases{$date}/$total_testcases{$date}</b></td>";

    foreach $name (sort keys %testcase)
    {
	$_ = $results{$date}{$name};
	print TABLE_FILE '<td> ', $_, '</td>';
    }
    print TABLE_FILE "</tr></td>\n";

    # store old month name for the next iteration of the loop
    $old_month = $this_month;
}

print TABLE_FILE << 'EOT'
</table>
</body>
</html>
EOT
    ;


# now print the names of the tests, in two blocks of 3 columns each. have an
# empty column of specified width in between, to separate the blocks
print NAMES_FILE << 'EOT'
<h3 align="center">Names of build tests</h3>
<table align="center">
<tr>
  <td align="center"><strong> Test number </strong></td>
  <td align="center"><strong> System name </strong></td>
  <td align="center"><strong> Compiler    </strong></td>
  <td style="center">------------------</td>
  <td align="center"><strong> Test number </strong></td>
  <td align="center"><strong> System name </strong></td>
  <td align="center"><strong> Compiler    </strong></td>
</tr>
EOT
    ;


# output the list of test cases. always put two in a row
$col = 0;
foreach $name (sort keys %testcase) {
    if ($col == 0) {
	print NAMES_FILE "<tr>\n";
    }

    $test_number   = $testcase{$name};
    $test_system   = $name;
    $test_system   =~ s/\+.*//;
    $test_compiler = $name;
    $test_compiler =~ s/.*\+//g;
    $compiler_table{"gcc295"} = "gcc 2.95";
    $compiler_table{"gcc30"} = "gcc 3.0";
    $compiler_table{"icc50"} = "Intel ICC 5.0";
    if (defined $compiler_table{$test_compiler} ) {
	$test_compiler = $compiler_table{$test_compiler};
    }

    print NAMES_FILE "   <td align=\"center\"><strong>${test_number}</strong></td>",
                     "   <td align=\"center\">${test_system}</td>",
                     "   <td align=\"center\">${test_compiler}</td>\n";

    # next column. if at end, wrap around
    $col = ($col+1)%2;

    # if at end of row, then end it, otherwise insert the empty column
    if ($col == 0) {
	print NAMES_FILE "</tr>\n";
    } else { 
	print NAMES_FILE "<td></td>"; 
    }
}

print NAMES_FILE <<'EOT'
</table>
</body>
</html>
EOT
    ;



print HEAD_FILE << 'EOT'
</p>
</body>
</html>
EOT
    ;
