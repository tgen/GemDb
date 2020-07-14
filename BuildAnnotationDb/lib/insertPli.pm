#!/usr/bin/perl 
package insertPli;
use initMongo;
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $key=shift;
      print "insertTwoCol: key: $key filepath:$filePath collection:$collection hostName: $hostName\n";
      my %done=();
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (GENE,"db/genes/germline.txt") or die "Can't open germline.txt\n";
      while (<GENE>) {
        chomp;
        $isGene{$_}=1;
      }
      open (TWOCOL,"$filePath") or die "Can't open file\n";
      LOOP: while (<TWOCOL>) {
        chomp;
        $line=$_;
        if ($line=~/(.*?)\t(.*)/) {
          if (!defined $1 && !defined $2) {next LOOP}
          my $gene=$1;
          my $label=$2;
          if (exists($done{$gene})) {next LOOP}
          unless (exists($isGene{$gene})) {next LOOP}
          $set= $GENE->find({'gene'=>$gene});
          if ($doc=$set->next) {
            $id=$doc->{'_id'};
#            if (exists($doc->{$key})) {
#              $idt=$GENE->update({"_id" => $id},{'$set'=>{$key=>join("|",$doc->{$key},$label)}});
#            } else {
              $idt=$GENE->update({"_id" => $id},{'$set'=>{$key=>$label}});
#            }
          } else {
#            $idt=$GENE->insert({$key=>$label,'gene'=>$gene,'name'=>$gene});
          }
          $done{$gene}=1;
        }
      }
      close (TWOCOL);
    }
1;
