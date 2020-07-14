#!/usr/bin/perl -w
$|=1;
package insertGerp;
use initMongo;
use POSIX;
    sub loadTo {
      my $filesInsert=addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $dbName=shift;
      print "-Module insertVcf $function $filePath $collection $hostName\n";
      my ($CONN,$DB,$VCF)=initMongo->new('AnnotateBy',$collection,$hostName);
      my ($fileconn,$filedb,$fileCollections)=initMongo->new('AnnotateBy','collections',$hostName);
#      $VCF -> remove({'filePath'=>$filePath});
#      $VCF -> update({cadd=>{'$exists'=>1}},{{'$unset'=>{'cadd'=>""}},{safe => 1});
      open (CADD,"$filePath") or die "Can't open file\n";
      while (<CADD>) {
          $line=$_;
          if ($line=~/^#/) {next}
          chomp($line);
          ($coord,$cadd)=split(/\t/,$line);
          #($chr,$pos,$ref,$alt,$score,$cadd)=split(/\t/,$line);
          MongoDB::force_double($cadd);
          #$coord="chr$chr:$pos:$ref:$alt";
          $varInfo={"coord"=>$coord,"GERP_RS"=>$cadd};
          $VCF->update({'coord'=>$coord},{'$set'=>$varInfo},{'upsert'=>1,'safe'=>1});
      }
      $fileCollections->update({'md5sum'=>$fileInsert->{'md5sum'}},$fileInsert,{'upsert'=>1,'safe'=>1});
      $idx = Tie::IxHash->new("coord" => 1); 
      $VCF->ensure_index($idx);
      return 1;
    }
sub addToCollections {
      my $function=shift;
      my $filePath=shift;
     my  $collection=shift;
     my  $hostName=shift;
      my $type='coord';
      my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
      my $md5sum=`md5sum $filePath`;
      chomp($md5sum);
      $md5sum=~s/ .*//g;
      my $wc=`wc -l $filePath`;
      $wc=~s/ .*//g; 
      my $filesInsert={'date'=>$date,'type'=>$type,'filename'=>$filePath,'md5sum'=>$md5sum,'collection'=>$collection,'lines'=>$wc,'md5sum'=>$md5sum};
      return $filesInsert;
}
1;
