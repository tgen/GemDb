#!/usr/bin/perl 
package insertInteractions;
use initMongo;
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      print "insertInteractions: filepath: $filePath collection: $collection hostname: $hostName\n";
      my ($conn,$db,$genes)=initMongo->new('AnnotateBy',$collection,$hostName);
      my ($connUpdate,$dbUpdate,$genesUpdate)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"$filePath") or die "Can't open file\n";
      my $file=$filePath;
      while (<FILE>) {
        chomp;
        my $geneName=$_;
        my $cur=$genes->find({"interactWith.$geneName"=>{'$exists' => 1}});
        $c=0;
        while (my $col=$cur->next) {
          ++$c;
          $id=$col->{'_id'};
          if (exists ($col->{'listInteraction'})) {
            $prev=$col->{'listInteraction'};
            $genesUpdate->update({'_id'=>$id},{'$set'=>{'listInteraction'=>"$prev,$geneName in $file"}},{'upsert' => 1});
          } else {
            $genesUpdate->update({'_id'=>$id},{'$set'=>{'listInteraction'=>"$geneName in $file"}},{'upsert' => 1});
          }
         }
      }
      close (FILE);
}

1;
