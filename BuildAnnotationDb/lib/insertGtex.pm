#!/usr/bin/perl 
package insertGtex;
use initMongo;
sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filename=shift;
      my $collection=shift;
      my $hostName=shift;
      print "insertGtex: collection: $collection hostname: $hostName\n";
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      $GENE->ensure_index({"RS"=>1});
      @ls=`ls db/gtex/Assoc-total*.txt`;
      foreach $file (@ls) {
          if ($file=~/Assoc-total\.(.*?)\.txt/) {
              $tissue=$1;
              open (FILE,"$file") or die "Can't open file\n";
              my $line=(<FILE>);
              print "Reading GTEX: $file\n";
              while (<FILE>) {
                  chomp;
                  @temp=split(/\t/);
                  $dbsnp=$temp[0];
                  if ($dbsnp=~/^rs(.*)/) {
                     $snp=$1;
                     $snp=~s/rs//g;
                     $ttest=$temp[9]*1.00;
                     $pvalue=$temp[10];
                     MongoDB::force_double($ttest);
#                     MongoDB::force_double($pvalue);
                     $name1="gtex_" . $tissue . "_test";
#$cur=$GENE->find({"RS" => $snp});
#print "snp:$snp\n";
#while ($doc=$cur->next){print "here:$doc->{'_id'} \n"}
                     $id=$GENE->update({"RS" => int($snp)},{'$set'=>{$name1=>$ttest}},{safe=>1});
                  }
              }
              close (FILE);
          }
      }
}
1;
