#!/usr/bin/perl 
package insertTsv;
use initMongo;
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $col=shift;
      print "insertTSV: filePath: $filePath collection: $collection hostname: $hostName col: $col\n";
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (DBSNP,"$filePath") or die "Can't open file\n";
      $line=(<FILE>);
      $line=lc($line);
      my @head=split(/,/,$line);
      $head[$col]="gene";
      LOOP: while (<DBSNP>) {
        chomp;
        $line=$_;
        if ($line=~/(.*)\t(.*)/) {
          $gene=$1;
          $label=$2;
          if (exists($done{$gene})) {next LOOP}
          $set= $GENE->find({'gene'=>$gene});
          if ($doc=$set->next) {
            $id=$doc->{'_id'};
            if (exists($doc->{$key})) {
              $idt=$GENE->update({"_id" => $id},{'$set'=>{$key=>join("|",$doc->{$key},$label)}});
            } else {
              $idt=$GENE->update({"_id" => $id},{'$set'=>{$key=>$label}});
            }
          } else {
            $idt=$GENE->insert({$key=>$label,'gene'=>$gene,'name'=>$gene});
          }
          $done{$gene}=1;
        }
      }
      close (DBSNP);
      $GENE->ensure_index({'gene'=>1});
    }
}

1;
