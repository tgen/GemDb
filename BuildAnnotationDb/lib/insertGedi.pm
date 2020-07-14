#!/usr/bin/perl -w
$|=1;
package insertGedi;
use initMongo;
use POSIX;

sub loadTo {
     my $filesInsert=addToCollections(@_);
     my $function=shift;
     my $filePath=shift;
     my $collection=shift;
     my $hostName=shift;
     print "-Module insertGedi.pm $function $filePath $collection $hostName\n";
     my ($CONN,$DB,$RULES)=initMongo->new('AnnotateBy',$collection,$hostName);
     my ($fileconn,$filedb,$fileCollections)=initMongo->new('AnnotateBy','collections',$hostName);
     $RULES -> remove({'filePath'=>$filePath});
     $RULES->ensure_index({'biomarker'=>1});
     open (GEDI,"$filePath") or die "Can't open file\n";
     my $head=(<GEDI>);
     my @headFields=split(/\t/,$head);        
     LOOP: while (<GEDI>) {
          $line=$_;
          chomp($line);
          my @fields=split(/\t/,$line);
          my $hashRef={};
          LOOP2: for ($i=0;$i<=$#fields;++$i) {
               if ($fields[$i] eq "EMPTY") {next LOOP2}
               $hashRef->{"$headFields[$i]"}=$fields[$i];
          }
          if (exists ($hashRef->{'biomarker_symbol'})) {
               if ($hashRef->{'biomarker_symbol'} ne "EMPTY") {
                  $biomarker=$hashRef->{'biomarker_symbol'};
                  $RULES->update({'biomarker'=>$biomarker},{'$set'=>$hashRef},{'safe'=>1});
               } elsif ($biomarker=$hashRef->{'modifier_symbol'} ne "EMPTY") {
                  $biomarker=$hashRef->{'modifier_symbol'};
                  $RULES->update({'biomarker'=>$biomarker},{'$set'=>$hashRef},{'safe'=>1});
               }

          }
     }
     close (GEDI);
     $fileCollections->update({'md5sum'=>$fileInsert->{'md5sum'}},$fileInsert,{'upsert'=>1,'safe'=>1});
     return 1;
}

sub addToCollections {
     my $function=shift;
     my $filePath=shift;
     my  $collection=shift;
     my  $hostName=shift;
     my $type='coord';
     print "filePath:$filePath\n";
     my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
     my $md5sum=`md5sum $filePath`;
     chomp($md5sum);
     $md5sum=~s/ .*//g;
     my $wc=`wc -l $filePath`;
     $wc=~s/ .*//g; 
     print "mdsub: $md5sum wc: $wc\n";
     my $filesInsert={'date'=>$date,'type'=>$type,'filename'=>$filePath,'md5sum'=>$md5sum,'collection'=>$collection,'lines'=>$wc,'md5sum'=>$md5sum};
     return $filesInsert;
}
1;
