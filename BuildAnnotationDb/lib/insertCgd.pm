#!/usr/bin/perl 
package insertCgd;
use initMongo;
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      print "insertCgd: filePath: $filePath collection: $collection hostname: $hostName\n";
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"$filePath") or die "Can't open file\n";
      my $line=(<FILE>);
      while (<FILE>) {
        chomp;
        @temp=split(/\t/);
        $gene=$temp[0];
        $condition=$temp[3];
        $inheritance=$temp[4];
        $age_group=$temp[5];
        $manifestation="$temp[7]";
        $intervention="$temp[8]|$temp[9]|$temp[10]";
        $references=$temp[11];
        $set= $GENE->find({'gene'=>$gene});
        if ($doc=$set->next) {
          $id=$doc->{'_id'};
          $id=$GENE->update({"_id" => $id},{'$set'=>{'cgdCondition'=>$condition,'cgdInheritance'=>$inheritance,'cgdAgeGroup'=>$age_group,'cgdManifestation'=>$manifestation,'cgdIntervention'=>$intervention,'name'=>$gene,'cgdRef'=>$references}});
        }
      }
     close (FILE);
}
1;
