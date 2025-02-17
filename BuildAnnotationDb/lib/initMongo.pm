#!/usr/bin/perl 
package initMongo;
use MongoDB;
use MongoDB::OID;
#use insertVcf;
use POSIX;
$MongoDB::BSON::use_binary = 4;

sub new {
		my $function=shift;
		my $database=shift;
		my $collection=shift;
		my $hostName=shift;
		#my $conn = MongoDB::Connection->new(host => "$hostName");
		my $timeOut=1000000;
		my $socketTimeOut=10000001;
		my $conn = MongoDB::MongoClient->new( host => "$hostName", connectTimeoutMS=> "$timeOut", socketTimeoutMS => "$socketTimeOut", socketCheckIntervalMS => "$timeOut");
		my $db = $conn->get_database($database);
		my $genes = $db->get_collection($collection);
		return ($conn,$db,$genes);
}

sub removeEntries {
     my $function=shift;
     my  $database=shift;
     my  $collection=shift;
     my  $hostName=shift;
      print "Removing Entries: function: $function collection: $collection hostName: $hostName\n";
      ($conn,$db,$dbsnp)=initMongo->new($database,$collection,$hostName);
      $dbsnp -> remove();
      $collection='collections';
      ($conn,$db,$dbsnp)=initMongo->new($database,$collection,$hostName);
      $dbsnp ->{'collection'=>$collection};
}

sub addToCollections {
      my $function=shift;
      $function=shift;
      my $filePath=shift;
      my  $collection=shift;
      my  $hostName=shift;
      my  $label=shift;
      my $type='gene';
      $type=$collection;
      my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
      my $md5sum='5';#`md5sum $filePath`;
      chomp($md5sum);
      $md5sum=~s/ .*//g;
      my $wc=15;#`wc -l $filePath`;
      $wc=~s/ .*//g;
      my $fileInsert={'date'=>$date,'type'=>$type,'filename'=>$filePath,'md5sum'=>$md5sum,'collection'=>$collection,'lines'=>$wc,'md5sum'=>$md5sum};
      my ($fileconn,$filedb,$fileCollections)=initMongo->new('AnnotateBy','collections',$hostName);
      $fileCollections->update({'md5sum'=>$fileInsert->{'md5sum'}},$fileInsert,{'upsert'=>1,'safe'=>1});
      return $fileInsert;
    }

1;
