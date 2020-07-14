#!/usr/bin/perl 
package insertClinvarNames;
use initMongo;
sub loadTo {
      print "-Module: @_\n";
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my ($conn,$db,$genes)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"$filePath") or die "Can't open file\n";
      my $line=(<FILE>);
      while (<FILE>) {
        chomp;
        my @temp=split(/\t/);
        $gene=$temp[1];
        $clin=$temp[3];
        $set= $genes->find({'gene'=>$gene});
        if ($doc=$set->next) {
          $id=$doc->{'_id'};
          if (exists($doc->{'clinvar'})) {
            $id=$genes->update({"_id" => $id},{'$set'=>{'clinvar' => join("| ",$clin,$doc->{'clinvar'})}});
          } else {
            $id=$genes->update({"_id" => $id},{'$set'=>{'clinvar' => $clin}});
          }
        } else {$id=$genes->insert({'clinvar'=>$clin,'gene'=>$gene});}
     }
}

1;
