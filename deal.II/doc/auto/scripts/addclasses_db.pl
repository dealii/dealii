#!/usr/local/bin/perl

# Update Postgres database entries from kdoc sources

use DBI;

$db = DBI->connect("DBI:Pg:dbname=deal",undef,undef);

die $DBI::errstr if (!$db);

print "Connected\n";

$db->do("UPDATE classes SET persistent=FALSE;")
    || die $db->errstr;

while(<>)
{
    next unless(/=/);
    next if (/\#/);
    s/=.*//g;
    chop;
    print $_, " ";

    # Try to find class name

    $sth = $db->prepare("UPDATE classes set persistent=TRUE where name = '$_';")
	|| die $db->errstr;
    $sth->execute || die $sth->errstr;

    $rows = $sth->rows;
    $sth->finish;

    # Insert name if it did not exist.

    if ($rows == 0)
    {
	print "insert";
	
	$db->do("INSERT INTO classes ( name, display, persistent, manually ) values ( '$_', TRUE, TRUE, FALSE);")|| die $db->errstr;;
    }
    print "\n";
}
