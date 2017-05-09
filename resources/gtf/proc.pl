#!/usr/bin/perl
while (<>) {
    chomp;
    @fields=split(/\t/);
    unless ($fields[2] eq "gene") {next}
    unless (/gene_biotype "protein_coding"/) {next}
    unless (/gene_source "ensembl_havana";/) {next}
    if (/; gene_name \"(.*?)\";/) {
        $gene{$1}=1;
    }
}
foreach $g (sort keys %gene) {print "$g\n";}
