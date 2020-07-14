#!/usr/bin/perl -w
package insertGnomad_b38;
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
      my ($conn,$db,$gnomad)=initMongo->new('AnnotateBy',$collection,$hostName);
      my @FILES=`ls -1 db/gnomad/b38/gnomad.genomes.r2.1.1.sites*bgz`;
	  foreach $file (@FILES) {
		  my $done={};
		  chomp($file);
		  print "-gnomad38 opening $file\n";
		  open(FILE,"gunzip -c $file |") or die "gunzip $file: $!\n";
		  $header=(<FILE>);
		  chomp($header);
		  #&addHead;
		  LOOP: while(<FILE>) {
			  next if $_ =~ m/^#/g;
			  $line = $_;
			  @temp = split("\t",$line);                   
			  %hash=();
			  %hashAnnotate=();
			  $hash{'coord'}=join(':',$temp[0],$temp[1],$temp[3],$temp[4]);                  
			  if (exists($done->{$hash{'coord'}})){next LOOP}               
			  $infoLine=$temp[7]; 
			  @infoFields = split(";",$infoLine);
			  foreach $infoField (@infoFields) {
			  	@pair=split("=",$infoField);
				if ($pair[0] eq "AC" || $pair[0] eq "AN" || $pair[0] eq "AF" || $pair[0] eq "nhomalt") {
					if ($pair[0] eq "AF"){
						$af_num = $pair[1];
						if ($af_num =~ /e/){
							$af_num = sn_to_dec($af_num);
						}
						$pair[1] = $af_num * 1;
					}
					else{
						$pair[1] = int($pair[1]);
					}
					$hash{"gnomad_$pair[0]"}=$pair[1];
				}
			  }
			  #print map { "$_ => $hash{$_}\n" } keys %hash; 
			 # $done->{$hash{'coord'}}=1;
			 $gnomad->update_one({'coord'=>$hash{'coord'}},{'$set'=>\%hash},{'upsert'=>1,'safe'=>1});
     		}
		close (FILE);
	 }
  $idx = Tie::IxHash->new("coord"=> 1); 
  $gnomad->ensure_index($idx);
}

sub addToCollections {
      my $function=shift;
      my $filePath=shift;
     my  $collection=shift;
     my  $hostName=shift;
      my $type='coord';
     # $filePath='db/dbnsfp/dbNSFP2.8_variant.chr22';
     $filepath=
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


sub sn_to_dec {
    my $num = shift;

    if ($num =~ /^([+-]?)(\d*)(\.?)(\d*)[Ee]([-+]?\d+)$/) {
        my ($sign, $int, $period, $dec, $exp) = ($1, $2, $3, $4, $5);

        if ($exp < 0) {
            my $len = 1 - $exp;
            $int = ('0' x ($len - length $int)) . $int if $len > length $int;
            substr $int, $exp, 0, '.';
            return $sign.$int.$dec;

        } elsif ($exp > 0) {
            $dec .= '0' x ($exp - length $dec) if $exp > length $dec;
            substr $dec, $exp, 0, '.' if $exp < length $dec;
            return $sign.$int.$dec;

        } else {
            return $sign.$int.$period.$dec;
        }
    }

    return $num;
}


1;

