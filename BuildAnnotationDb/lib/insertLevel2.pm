#!/usr/bin/perl -w
   package insertToLevel2;
   use initMongo;
   use POSIX;
sub removeEntries {
     $function=shift;
      $collection=shift;
      $hostName=shift;
      print "--insertToPanels::reset $function $collection $hostName\n";
      ($conn,$db,$dbsnp)=initMongo->new('AnnotateBy',$collection,$hostName);
      $dbsnp -> remove();
}
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $label=shift;
$label='level2';
      print "InsertToPanels $filePath $collection $hostName $label\n";
      my ($conn,$db,$genes)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (DBSNP,"$filePath") or die "DIE:  Can't open file:  $filePath\n";
      my %head=();
      my %done=();
      my $idt="";
      my $id="";
      LOOP: while (<DBSNP>) {
        chomp;
        $id="";
        @temp=split(/\t/);
        $gene=$temp[0];
        if (exists($done{$gene})) {next LOOP}
        $set=$genes->find({'gene'=>$gene});
        if ($doc=$set->next) {
          $id=$doc->{'_id'};
          $idt=$genes->update({"_id" => $id},{'$set'=>{'level2'=>$label}});
        }
        $done{$gene}=1;
      }
      close (DBSNP);
}
1;
