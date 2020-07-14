#!/usr/bin/perl -w
package insertRefseq;
use initMongo;

use Storable qw(dclone);
use Data::Dumper;

sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"db/refseq/refLink.txt");
      %lookup=();
      while (<FILE>) {
          @temp=split;
          $lookup{$temp[4]}=$temp[0];
      }
      close (FILE);
      
      open (FILE,"$filePath") or die "Can't open file: $filePath\n";
     LOOP1: while (<FILE>) {
        $line=$_;
        chomp($line);
        @temp=split(/\t/,$line);
        if (exists($lookup{$temp[0]}) && defined($temp[2])) {
           $gene=$lookup{$temp[0]};
           if (length($temp[2])>1) {
           $text=$temp[2];
           $text=~s/Publication Note.*//g;
           $text=~s/Sequence Note.*//g;
           $text=~s/##Evidence.*//g;
               #print "$gene:  $text\n";
               $GENE->update({"gene"=>$gene},{'$set'=>{'refseq'=>$text}},{safe=>1});
           }
        }
      }
     close (FILE);
}
1;
