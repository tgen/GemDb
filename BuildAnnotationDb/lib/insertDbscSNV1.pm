#!/usr/bin/perl -w
#dbnsfp2.4
package insertDbscSNV1;
use initMongo;
use POSIX;
use Time::HiRes qw( usleep );
sub loadTo {
      print "-Module: @_\n";
#      my $md5=addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;    
      my ($conn,$db,$dbnsfp)=initMongo->new('AnnotateBy',$collection,$hostName);
      my @FILES=`ls -1 db/dbnsfp3.4/dbscSNV1.1.chr*`;
      foreach $file (@FILES) {
              chomp($file);
              print "-dbnsfp opening $file\n";
              open(FILE,"<$file") or die "Can't open $file\n";
              $header=(<FILE>);
              chomp($header);
              #$head[0]='chr'; hg19
              #$head[1]='pos'; hg19
              $head[2]='ref';
              $head[3]='alt';
              $head[4]='chr';
	      $head[5]='pos';
	      $head[16]='adaSnv';
              $head[17]='rfSnv';
$line=<FILE>;
LOOP: while(<FILE>) {
          $line=$_;
          chomp($line);
          @temp = split("\t",$line);                   
          %hash=();
          $hash{'coord'}=join(':',"chr" . $temp[4],$temp[5],$temp[2],$temp[3]);                  
          for $i (16,17){
              if ($temp[$i] ne ".") {
                 $temp[$i]=sprintf("%0.3f",$temp[$i])*1.00;
                 MongoDB::force_double($temp[$i]);
                 $hash{$head[$i]}=$temp[$i];
              }
          }
#          $dbnsfp->insert(\%hash,{'safe'=>1});
          $dbnsfp->update({'coord'=>$hash{'coord'}},{'$set'=>\%hash},{'upsert'=>1,'safe'=>1});
#print $hash{'coord'};die;
      }
      close (FILE);
  }
  $idx = Tie::IxHash->new("coord"=> 1); 
  $dbnsfp->ensure_index($idx);
}

sub addToCollections {
      my $function=shift;
      my $filePath=shift;
     my  $collection=shift;
     my  $hostName=shift;
      my $type='coord';
     # $filePath='db/dbnsfp/dbNSFP2.8_variant.chr22';
      my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
#      my $md5sum=`md5sum $filePath`;
#      $md5sum=~s/ .*//g;
#      chomp($md5sum);
      my $wc=`wc -l $filePath`;
      $wc=~s/ .*//g;
      my ($conn,$db,$collections)=initMongo->new('denote','collections',$hostName);
      $collections->update({'md5sum'=>$md5sum},{'date'=>$date,'type'=>$type,'collection'=>$collection},{'upsert'=>1,'safe'=>1});
      return;
}
1;

