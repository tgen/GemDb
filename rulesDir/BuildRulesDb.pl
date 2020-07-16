#!/usr/bin/perl -w
use POSIX;
use FindBin;                 # locate this script
require "$FindBin::Bin/insertRules.pm";  # use the parent directory
require "$FindBin::Bin/initMongo.pm";  # use the parent directory

$hostName=#your hostname;

foreach my $arg (@ARGV) {
    chomp($arg);
    if ($arg=~/(.*)=(.*)/) {
       $key=$1;$file=$2;
       if ($key eq "rules") {
            print "inserting rules";
	   insertRules->loadTo ($file,'rules',$hostName,'AnnotateBy');
       }elsif ($key eq "bbb") {
          insertRules->joinTo ($file,'rules',$hostName,'AnnotateBy',"drug");
       }
    }
}
