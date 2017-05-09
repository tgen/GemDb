#!/usr/bin/perl
%stuff=();
while (<>) {
if (/gene_name "(.*?)";/) {
   $stuff{$1}=1;
}
}
foreach $st (keys %stuff) {
  print "$st\n";
}
