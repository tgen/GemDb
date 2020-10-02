#!/usr/bin/perl 
#/*************************************************************************
# *
# * TGen CONFIDENTIAL
# * __________________
# *
# *  [2010] - [2014] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * NOTICE:  All information contained herein is, and remains
# * the property of Translational Genomics Research Institute (TGen).
# * The intellectual and technical concepts contained herein are proprietary
# * to  TGen and may be covered by U.S. and Foreign Patents,
# * patents in process, and are protected by trade secret or copyright law.
# * Dissemination of this information, application or reproduction
# * is strictly forbidden unless prior written permission is obtained
# * from TGen.
# *
# * Major Contributor(s):
#    David Craig
#/
######################################################################################

$| = 1;

package insertVcf;
use initMongo;
#use MongoDB;
#use MongoDB::Connection;
use POSIX;
use Scalar::Util qw(looks_like_number);

sub loadTo {
    my $filesInsert = addToCollections(@_);
    my $function    = shift;
    my $filePath    = shift;
    my $collection  = shift;
    my $hostName    = shift;
    my $dbName      = shift;
    print "-Module insertVcf $function $filePath $collection $hostName\n";
    my ( $CONN,     $DB,     $VCF )             = initMongo->new( 'AnnotateBy', $collection,   $hostName );
    my ( $fileconn, $filedb, $fileCollections ) = initMongo->new( 'AnnotateBy', 'collections', $hostName );

    #      $VCF -> remove({'filePath'=>$filePath});
          $VCF->ensure_index({'coord'=>1});
    open( DBSNP, "$filePath" ) or die "gunzip $filePath $!";
    my %head        = ();
    my $c           = 0;
    my %types       = ();
    my %typesNumber = ();
    $types{"AF"}=2;
    $types{"1000Gp3_AF"}=2;
    $typesNumber{"AF"}="A";
    $typesNumber{"1000Gp3_AF"}="A";
  LOOP3: while (<DBSNP>) {
        $line = $_;
        if ( $line =~ /#CHR/ ) { last LOOP3 }
        if ( $line =~ /##INFO=<ID=(.*?),Number=(.*?),Type=(.*?),/ ) {
            $field  = $1;
            $number = $2;
            $type   = $3;
            if ( $type eq "Integer" ) {
                $types{$field} = 1;
            }  elsif ( $type eq "Float" ) {
                $types{$field} = 2;
            } else {
                $types{$field} = 3;
            }
            $typesNumber{$field} = $number;
        }

    }
  LOOP: while (<DBSNP>) { 
        ++$c;
        if ( $c % 10000 == 0 ) { print "Counter is at : $c\n"; }
        $line = $_;
        chomp($line);
        @temp = split( /\t/, $line );
        @alts = ();
        ( $chr, $hg19pos, $dbsnp, $ref, $alt, $qual, $filter, $info ) =
          ( $temp[0], int( $temp[1] ), $temp[2], $temp[3], $temp[4], $temp[5], $temp[6], $temp[7] );
        if ( $filter =~ /VQSR/ ) { next LOOP }
        $chr =~ s/chr//g;
        if ( $alt =~ /\,/ ) { @alts = split( /\,/, $alt ); }
        else                { $alts[0] = $alt; }
        if ( $ref =~ /\,/ ) { @refs = split( /\,/, $ref ); }
        else                { $refs[0] = $ref; }
        my $altCount = -1;
      LOOP1: foreach $alt (@alts) {
            ++$altCount;
            foreach $ref (@refs) {
              my $coord = "chr$chr:$hg19pos:$ref:$alt";
              my $varInfo = { "coord" => $coord };
              @words = split( /\;/, $info );
              my %pop=();                
              LOOP2: foreach $word (@words) {
                    if ( $word =~ /(.*?)=(.*)/ ) {
                        $key = $1;
                        $val = $2;
                        $key =~ s/\./_/g;
                        $key =~ s/\ /_/g;
                        if ( $key eq "CLNSIG" ) {
                            $insert = 1;
                            $varInfo->{'clinvarSignificance'} = $val;
                        }
                        elsif ( $key eq "CLNDBN" ) {
                            $varInfo->{'clinvarDiseaseKnownVariant'} = $val;
                            $insert = 1;
                        }
                        elsif ( $key eq "CLNREVSTAT" ) {
                            $varInfo->{'clinvarReviewStatus'} = $val;
                            $insert = 1;
                        }
                        elsif ( $key eq "CLNSIGINCL" ) {
                            $varInfo->{'clinvarClinicalSignificance'} = $val;
                            $insert = 1;
                        }
                        elsif ( $key eq "ORIGIN" ) {
                            @vals = split( /\||\,/, $val );
                            foreach $val (@vals) {
                                if ( $val == 512 ) {
                                    $varInfo->{'clinvarOrigin'} .= "tested-inconclusive|";
                                }
                                elsif ( $val == 1073741824 ) {
                                    $varInfo->{'clinvarOrigin'} .= "other|";
                                }
                                elsif ( $val == 128 ) {
                                    $varInfo->{'clinvarOrigin'} .= "uniparental|";
                                }
                                elsif ( $val == 64 ) {
                                    $varInfo->{'clinvarOrigin'} .= "biparental|";
                                }
                                elsif ( $val == 32 ) {
                                    $varInfo->{'clinvarOrigin'} .= "de-novo|";
                                }
                                elsif ( $val == 16 ) {
                                    $varInfo->{'clinvarOrigin'} .= "maternal|";
                                }
                                elsif ( $val == 8 ) {
                                    $varInfo->{'clinvarOrigin'} .= "paternal|";
                                }
                                elsif ( $val == 4 ) {
                                    $varInfo->{'clinvarOrigin'} .= "somatic|";
                                }
                                elsif ( $val == 2 ) {
                                    $varInfo->{'clinvarOrigin'} .= "germline|";
                                }
                                elsif ( $val == 1 ) {
                                    $varInfo->{'clinvarOrigin'} .= "unknown|";
                                }
                            }
                        }
			elsif ( index($dbName, "cosmic") != -1 && $key eq "ID") {
				$varInfo->{'cosmicID'} = $val;
				$insert		       = 1;
			}
                        elsif ( index($dbName, "cosmic") != -1 && $key eq "CNT" ) {
                            $varInfo->{'cosmicCount'} = int($val);
                            $insert                   = 1;
                            $varInfo->{$key}          = $val;
                        } elsif ($key eq "GENEINFO" || $key eq "CLNALLE" || $key eq "CLNRCID" || $key eq "CLNACC" || $key eq "RSPOS" || $key eq "VC" || $key eq "WGT" || $key eq "SAO" || $key eq "VP"  || $key eq "SSR" || $key eq "CLNDSDBID" || $key eq "CLNSRCID" ) {
                            next LOOP2;
                        } elsif ($key eq "CSQ" || $key eq "GQ_HIST" || $key eq "DP_HIST" || $key eq "ReadPosRankSum" || $key eq "VQSLOD" || $key eq "culprit" ) { next LOOP2; 
                        } elsif ($key =~/AN_(.*)/) {
                            $t=$1;
						   if ( $typesNumber{$key} eq "A" ) {
							   my @ind = split( ",", $val );
							   $val = @ind[$altCount];
						   }                        
                            $pop{$t}{'AN'}=$val;                            
                        } elsif ($key =~/Het_(.*)/) {
                            $t=$1;                        
						   if ( $typesNumber{$key} eq "A" ) {
							   my @ind = split( ",", $val );
							   $val = @ind[$altCount];
						   }                                                
                            $pop{$t}{'HET'}=$val;                            
                        }elsif ($key =~/Hom_(.*)/) {
                            $t=$1;                        
						   if ( $typesNumber{$key} eq "A" ) {
							   my @ind = split( ",", $val );
							   $val = @ind[$altCount];
						   }                                                
                            $pop{$t}{'HOM'}=$val;                            
                        }elsif ($key =~/Hemi_(.*)/) {
                            $t=$1;                        
						   if ( $typesNumber{$key} eq "A" ) {
							   my @ind = split( ",", $val );
							   $val = @ind[$altCount];
						   }                        
                            $pop{$t}{'HEM'}=$val;                            
                        }else {
                            if ( exists( $types{$key} ) ) {
                                my $index = 0;
                                if ( $typesNumber{$key} eq "A" ) {
                                    my @ind = split( ",", $val );
                                    $val = @ind[$altCount];
                                }
                                if ( $types{$key} == 1 ) {
                                    if ($dbName eq 'ExAc') { $key="ExAC_$key";}
				    if ($dbName eq 'gnomad') { $key="gnomad_$key";}
                                    if ($dbName eq 'phase3KG') { $key="1000Gp3_$key";}
                                    $varInfo->{$key} = int($val);
                                }
                                elsif ( $types{$key} == 2 ) {
                                    if ($dbName eq 'ExAc') { $key="ExAC_$key";}
				    if ($dbName eq 'gnomad') { $key="gnomad_$key";}
                                    if ($dbName eq 'phase3KG') { $key="1000Gp3_$key";}
                                    $varInfo->{$key} = $val;
				    if (looks_like_number($varInfo->{$key})) {	
					MongoDB::force_double( $varInfo->{$key} );
				    }
                                }
                                else {
                                    if ($dbName eq 'phase3KG') { $key="1000Gp3_$key";}
				    if ($dbName eq 'gnomad') { $key="gnomad_$key";}
                                    if ($dbName eq 'ExAc') { $key="ExAC_$key";}
                                    $varInfo->{$key} = $val;
                                }
                            }
                            else {
                                $varInfo->{$key} = $val;
                            }
                        }
                        if ($key="RS") {$variants{'RS'}=$val}
			if ($key="dbSNPBuildID") {$variants{'dbSNPBuildID'}=$val}
                    }
                    else {
                        $varInfo->{'clinvarNote'} = "";
                        if ( $word eq "PM" ) {
                            $varInfo->{'clinvarNote'} .= "Variant is pubmed cited";
                        }
                        if ( $word eq "TPA" ) {
                            $varInfo->{'clinvarNote'} .= "Provisional Third Party Annotation";
                        }
                        if ( $word eq "PMC" ) {
                            $varInfo->{'clinvarNote'} .= "Links exist to PubMed Central article";
                        }
                        if ( $word eq "MUT" ) {
                            $varInfo->{'clinvarNote'} .= "Cited in journal and other reputable sources";
                        }
                        if ( $word eq "CDA" ) {
                            $varInfo->{'clinvarNote'} .= "Variation is interrogated in a clinical diagnostic assay";
                        }
                        if ( $word eq "OM" ) {
                            $varInfo->{'clinvarNote'} .= "Variation has Omim entry";
                        }
                        if ( $word eq "G5" ) {
                            $varInfo->{'MinPopFreq'} = 0.05;
			    	MongoDB::force_double( $varInfo->{'MinPopFreq'} );
                        }
                        if ( $word eq "G5A" ) {
                            $varInfo->{'MinPopFreq'} = 0.05;
                                MongoDB::force_double( $varInfo->{'MinPopFreq'} );
                        }
                        if ( $word eq "KGPROD" ) {
                            $varInfo->{'submitter'} .= "KGPROD|";
                        }
                        if ( $word eq "OTHERKG" ) {
                            $varInfo->{'submitter'} .= "OTHERKG|";
                        }
                        if ( $word eq "PH3" ) {
                            $varInfo->{'submitter'} .= "PH3";
                        }
                    }
                }
                ++$c;
                foreach $popCount (keys %pop) {
#print "pop: $popCount $pop{$popCount}{'AN'} $coord\n";
                    if (exists($pop{$popCount}{'AN'})&& exists($pop{$popCount}{'HOM'}) &&  exists($pop{$popCount}{'HEM'})&&  exists($pop{$popCount}{'HET'})) {
                        $varInfo->{$popCount . "_cnts"} = join("/",$pop{$popCount}{'AN'},$pop{$popCount}{'HOM'},$pop{$popCount}{'HET'},$pop{$popCount}{'HEM'});
                        $freq=($pop{$popCount}{'HOM'}*2+$pop{$popCount}{'HEM'}+$pop{$popCount}{'HET'})/($pop{$popCount}{'AN'}+1);
                        $freq=sprintf("%2.6f",$freq);
			MongoDB::force_double($varInfo->{$popCount . "_maf"});
                    } elsif (exists($pop{$popCount}{'AN'}) && exists($pop{$popCount}{'HET'}) && exists($pop{$popCount}{'HOM'})) {
                        $varInfo->{$popCount . "_cnts"} = join("/",$pop{$popCount}{'AN'},$pop{$popCount}{'HOM'},$pop{$popCount}{'HET'});
                        $freq=($pop{$popCount}{'HOM'}*2+$pop{$popCount}{'HET'})/($pop{$popCount}{'AN'}+1);
                        $freq=sprintf("%2.6f",$freq);
                        $varInfo->{$popCount . "_maf"}=$freq;
			MongoDB::force_double($varInfo->{$popCount . "_maf"});
                    }
                }
                

                undef( $varInfo->{""} );
                delete( $varInfo->{""} );
                if ( length($coord) > 200 ) { next LOOP; }
                $VCF->update( { 'coord' => $coord }, { '$set' => $varInfo }, { 'upsert' => 1, 'safe' => 1 } );
                if ( $c % 100000 == 0 ) { print "Count: $c\n"; }
            }
        }
    }
    $fileCollections->update( { 'lines' => $fileInsert->{'lines'} }, $fileInsert, { 'upsert' => 1, 'safe' => 1 } );
    $idx = Tie::IxHash->new( "coord" => 1 );
    $VCF->ensure_index($idx);
    return 1;
}

sub addToCollections {
    my $function   = shift;
    my $filePath   = shift;
    my $collection = shift;
    my $hostName   = shift;
    my $type       = 'coord';
    print "filePath:$filePath\n";
   # my $date = "";#POSIX::strftime( "%d%m%y", localtime( ( stat $filePath )[9] ) );

    my $wc = 15;#`wc -l $filePath`;
   my $wc="";
    $wc =~ s/ .*//g;
    my $filesInsert = {
        'type'       => $type,
        'filename'   => $filePath,
        'collection' => $collection,
        'lines'      => $wc,
        'md5sum'     => $md5sum
    };
    return $filesInsert;
}
1;
