#!/usr/bin/perl

while(<>)
{
    1 while s/<[^<>]*>//g;
    print;
}

