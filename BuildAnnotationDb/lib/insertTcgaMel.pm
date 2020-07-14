#!/usr/bin/perl 
package insertTcgaMel;
use initMongo;
sub loadTo {
      print "-Module: @_\n";
      $function=shift;
      $filePath=shift;
      $collection=shift;
      $hostName=shift;
      ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"$filePath") or die "Can't open file\n";
     $line=(<FILE>);
     $line=lc($line);
     @head=split(/,/,$line);
     $head[0]="gene";
     while (<FILE>) {
        chomp;
        @temp=split(/,/,$_);
        %hins=();
        $hins{'gene'}=$temp[0];
        $hins{'name'}=$temp[0];
        $hins{'coord'}="$temp[4]:$temp[5]:$temp[10]:$temp[11]";
        $GENE->insert(\%hins);
     }
     close (FILE);
}
1;
