#!/usr/bin/perl -w
package insertUniprot;
use initMongo;

use Storable qw(dclone);
use Data::Dumper;
use XML::XML2JSON;

sub loadTo {
      my $fileInsert=initMongo->addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      $Data::Dumper::Indent = 0;     
      $Data::Dumper::Pair = " : ";    
      my $XML2JSON = XML::XML2JSON->new();      
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (FILE,"$filePath") or die "Can't open file: $filePath\n";

     LOOP1: while (<FILE>) {
        $line=$_;
        $buf="";
        if ($line=~/\<entry/) {
               $buf=$buf . $_;
               LOOP: while (<FILE>) {
                    $buf=$buf . $_;
                   if (/\<\/entry\>/) {
                        my $Obj = $XML2JSON->xml2obj($buf);            
                       # $XML2JSON->sanitize($Obj);
                        my $original=dclone($Obj);
                        my $new={}; 
                        #my $JSON = $XML2JSON->convert($buf);       
                        #my $str = Dumper $JSON;
                        #$str=~s/\$VAR1 = //;
                        #print "$str\n\n";      
                                                                    
                        unless (  ref($original->{'entry'}->{'gene'}) eq 'HASH'){next LOOP1}
                        $new->{'gene'}="unknown";
                        if (ref( $original->{'entry'}->{'gene'}->{'name'}) eq 'ARRAY') {
                            $new->{'gene'}=$original->{'entry'}->{'gene'}->{'name'}->[0]->{'$t'} ;                   
                        } else {next LOOP1} #elsif( ref($original->{'entry'}->{'gene'}->{'name'}) eq 'HASH') {
                          #  $new->{'gene'}=$original->{'entry'}->{'gene'}->{'name'}->[0]->{'$t'} ;   
                        #}
                        unless ( $original->{'entry'}->{'name'}->{'$t'}=~/_HUMAN/) {next LOOP1}
              
                        if (ref($original->{'entry'}->{'feature'}) eq 'ARRAY' ) {
                            foreach my $feature (@{$original->{'entry'}->{'feature'}}) {     
                                if (exists ($feature->{'location'}) && exists($feature->{'original'}) && exists($feature->{'variation'}) && exists($feature->{'original'}) && exists ($feature->{'location'}->{'position'})) {                                
                                    if (ref($feature->{'variation'}) eq "ARRAY") {
                                         foreach my $val (@{$feature->{'variation'}}) {
                                             my $ref="";
                                             if (exists($feature->{'@evidence'})) {
                                                 foreach my $eviNum (split(/ /,$feature->{'@evidence'})) {
                                                     my $eviNumA=$eviNum-1;
                                                     if (ref($original->{'entry'}->{'evidence'}) eq "ARRAY") {
                                                         if (exists($original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'})) {
                                                            $ref=$ref . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@type'} . 
                                                            ":" . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@id'};
                                                         }   
                                                      }
                                                }                                                  
                                             }
                                             if (length($ref) <2 ) {$ref="Uniprot"}
                                             my $effectVal=$new->{'gene'}  . "(" .  $feature->{'original'}->{'$t'} . $feature->{'location'}->{'position'}->{'@position'} . $val->{'$t'}. ")";
                                             if (length($effectVal)<100) {
                                                 push(@{$new->{'effect'}}, {'codon'=>$feature->{'location'}->{'position'}->{'@position'},'variant'=>$effectVal,'description'=>$feature->{'@type'} . ":" . $feature->{'@description'} . "(". $ref .")" });
                                            }
                                         }
                                    } else {
                                        my $ref="";
                                        if (exists($feature->{'@evidence'})) {
                                             foreach my $eviNum (split(/ /,$feature->{'@evidence'})) {
                                                 my $eviNumA=$eviNum-1;
                                                 if (ref($original->{'entry'}->{'evidence'}) eq "ARRAY") {
                                                     if (exists($original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'})) {
                                                        $ref=$ref . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@type'} . 
                                                        ":" . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@id'};
                                                     }   
                                                  } 
                                            }                                                  
                                         }
                                         if (length($ref) <2 ) {$ref="Uniprot"}  
                                          my $effectVal=$new->{'gene'}  . "(" .  $feature->{'original'}->{'$t'} . $feature->{'location'}->{'position'}->{'@position'} . $feature->{'variation'}->{'$t'}. ")";                                
                                          if (length($effectVal)<100) {                                           
                                             push(@{$new->{'effect'}}, {'codon'=>$feature->{'location'}->{'position'}->{'@position'},'variant'=>$effectVal,'description'=>$feature->{'@type'} . ":" . $feature->{'@description'}. "(". $ref .")"});
                                          }
                                         
                                    }
                                 }
                             }                                                  
                         }
                         $new->{'UniprotDisease'}=""; 
                         $new->{'UniprotFunction'}="";                                                      
                        if (ref($original->{'entry'}->{'comment'}) eq 'ARRAY' ) {
                           foreach my $comment (@{$original->{'entry'}->{'comment'}}) {
                               if (exists($comment->{'text'}->{'$t'})) {
                                   my $ref="";
                                   if (exists($comment->{'@evidence'})) {
                                      foreach my $eviNum (split(/ /,$comment->{'@evidence'})) {
                                         my $eviNumA=$eviNum-1; 
                                         if (ref($original->{'entry'}->{'evidence'}) eq "ARRAY") {                                        
                                             if (exists($original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'})) {
                                                $ref=$ref . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@type'} . 
                                                ":" . $original->{'entry'}->{'evidence'}->[$eviNumA]->{'source'}->{'dbReference'}->{'@id'} . " ";
                                             }
                                         }
                                      }
                                   }
                                   if (length($ref) <2 ) {$ref="Uniprot"}
                                   if ($comment->{'@type'} eq "disease") {
                                       if (exists($comment->{'disease'}->{'name'})) {
                                           $new->{'UniprotDisease'}=$new->{'UniprotDisease'} . "  " . $comment->{'disease'}->{'name'}->{'$t'};
                                       }
                                       if (exists($comment->{'disease'}->{'description'}->{'$t'} )) {
                                           $new->{'UniprotDisease'}=$new->{'UniprotDisease'} . "  " .$comment->{'disease'}->{'description'}->{'$t'};
                                       }
                                       if (exists($comment->{'text'})) {
                                           $new->{'UniprotDisease'}=$new->{'UniprotDisease'} .  $comment->{'text'}->{'$t'};
                                       }
                                       $new->{'UniprotDisease'}=$new->{'UniprotDisease'} . " (" . $ref . ");"
                                   } 
                                   if ($comment->{'@type'} eq "function") {
                                       $new->{'UniprotFunction'}=$new->{'UniprotFunction'} . "  " . $comment->{'text'}->{'$t'} . " (" . $ref . ")"
                                   }     
                                   push(@{$new->{'comment'}},{'text'=>$comment->{'text'}->{'$t'} . "(" . $ref . ")","type"=>$comment->{'@type'}});
                               }
                           }
                         }                           
                         $GENE->update({"gene" => $new->{'gene'}},{'$set'=>$new},{safe=>1});
                         print "Gene Added: $new->{'gene'}\n";                                             
                         last LOOP;
                     }
               }
         }        
    }
     close (FILE);
     $GENE->ensure_index({'effect.variant'});
     $GENE->ensure_index({'effect.codon'});
     
    
}
1;
