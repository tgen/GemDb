#!/usr/bin/perl -w
package insertUniprotAnnotate;
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
      my $XML2JSON = XML::XML2JSON->new();      print "$hostName\n";
      my ($conn,$db,$GENE)=initMongo->new('AnnotateBy','gene',$hostName);
      my ($conn,$db,$GENECODON)=initMongo->new('AnnotateBy','geneCodon',$hostName);
      
      open (FILE,"$filePath") or die "Can't open file: $filePath\n";
     $GENE->ensure_index({'gene'=>1});     
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
                        my $codons={};
                                                                    
                        unless (  ref($original->{'entry'}->{'gene'}) eq 'HASH'){next LOOP1}
                        $new->{'gene'}="unknown";
                        if (ref( $original->{'entry'}->{'gene'}->{'name'}) eq 'ARRAY') {
                            $new->{'gene'}=$original->{'entry'}->{'gene'}->{'name'}->[0]->{'$t'} ;                   
                        } else {next LOOP1}
                        unless ( $original->{'entry'}->{'name'}->{'$t'}=~/_HUMAN/) {next LOOP1}
              
                        if (ref($original->{'entry'}->{'feature'}) eq 'ARRAY' ) {
                            LOOP: foreach my $feature (@{$original->{'entry'}->{'feature'}}) {   
                                if ($feature->{'@type'} eq 'chain') {
                                     my $b=$feature->{'location'}->{'begin'}->{'@position'};
                                     my $e=$feature->{'location'}->{'end'}->{'@position'};                                                                       
                                    push (@{$new->{'features'}},{'type'=>'chain','codons'=>"$b - $e"});
                                    next LOOP;
                                }
                                if (exists ($feature->{'location'}) && exists($feature->{'original'}) && exists($feature->{'variation'}) && exists ($feature->{'location'}->{'position'})) {                     
                                    if (ref($feature->{'variation'}) eq "ARRAY") {
                                         foreach my $val (@{$feature->{'variation'}}) {
                                             my $effectVal=$new->{'gene'}  . " (" .  $feature->{'original'}->{'$t'} . $feature->{'location'}->{'position'}->{'@position'} . $val->{'$t'}. ")";
                                             if (length($effectVal)<100) {
                                                 my $gc=$new->{'gene'} . '_' . $feature->{'location'}->{'position'}->{'@position'};
                                                # push (@{$new->{'features'}},{'type'=>'variation','codons'=>$feature->{'location'}->{'position'}->{'@position'}});                                                                                                                                                           
                                                 if (!exists($codons->{ $gc })) {
                                                     $codons->{ $gc }= {'geneCodon'=> $gc,'codon'=>$feature->{'location'}->{'position'}->{'@position'},'description'=>$feature->{'@type'} . " (" . $feature->{'@description'} . " $effectVal)" };
                                                 } else {
                                                     if (exists($codons->{ $gc }->{'description'})) {
                                                        $codons->{ $gc }->{'description'} .= ", " . $feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" ;
                                                     } else {
                                                        $codons->{ $gc }->{'description'} = $feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" ;                                                     
                                                     }
                                                 }
                                                 $codons->{ $gc }->{'gene'}=$new->{'gene'};
                                            }
                                         }
                                    } else { 
                                          my $effectVal=$new->{'gene'}  . "(" .  $feature->{'original'}->{'$t'} . $feature->{'location'}->{'position'}->{'@position'} . $feature->{'variation'}->{'$t'}. ")";   
                                          if (length($effectVal)<100) { 
                                                 my $gc=$new->{'gene'} . '_' . $feature->{'location'}->{'position'}->{'@position'};
                                                 if (!exists($codons->{ $gc })) {
                                                     $codons->{ $gc }= {'geneCodon'=> $gc,'codon'=>$feature->{'location'}->{'position'}->{'@position'},'description'=>$feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" };
                                                 } else {
                                                     if (exists($codons->{ $gc }->{'description'})) {
                                                        $codons->{ $gc }->{'description'}  .= ", " . $feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" ;
                                                     } else {
                                                        $codons->{ $gc }->{'description'} = $feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" ;                                                     
                                                     }
                                                 }
                                                 $codons->{ $gc }->{'gene'}=$new->{'gene'};                                                 
                                                #push (@{$new->{'features'}},{'type'=>'variation','codons'=>$feature->{'location'}->{'position'}->{'@position'}});
                                          }
                                         
                                    }
                                 }  elsif ( exists ($feature->{'location'}->{'begin'}) && exists ($feature->{'location'}->{'end'})) {
                                     my $b=$feature->{'location'}->{'begin'}->{'@position'};
                                     my $e=$feature->{'location'}->{'end'}->{'@position'};  
                                     my $val=$feature->{'@type'} . " (" . $feature->{'@description'} . ")" ;
                                     $val=~s/\(\)//;                                     
                                     if ($feature->{'@type'} eq 'domain') {
                                         push (@{$new->{'features'}},{'type'=>$val ,'codons'=>"$b - $e"});                                                                                                     
                                     }
                                     for (my $i=$b;$i<=$e;++$i) {
											 my $gc=$new->{'gene'} . '_' . $i;  
											 if (!exists($codons->{ $gc })) {
												 $codons->{ $gc }= {'geneCodon'=> $gc,'codon'=>$i,'description'=>$feature->{'@type'} . " (" . $feature->{'@description'} .  " $effectVal)" };
											 } else {
												 if (exists($codons->{ $gc }->{'description'})) {
													$codons->{ $gc }->{'description'} .= ", " . $feature->{'@type'} . " (" . $feature->{'@description'} . " $effectVal)" ;
												 } else {
													$codons->{ $gc }->{'description'} = $feature->{'@type'} . " (" . $feature->{'@description'} . " $effectVal)" ;                                                     
												 }
											 }    
                                             $codons->{ $gc }->{'gene'}=$new->{'gene'};											                                  
                                    }
                                 } elsif ( exists ($feature->{'location'}->{'position'})) {
                                     my $val=$feature->{'@type'} . " (" . $feature->{'@description'} . ")" ;
                                     $val=~s/\(\)//;
                                    # push (@{$new->{'features'}},{'type'=> $val ,'codons'=>$feature->{'location'}->{'position'}->{'@position'}});
										 my $gc=$new->{'gene'} . '_' . $feature->{'location'}->{'position'}->{'@position'};
										 if (!exists($codons->{ $gc })) {
											 $codons->{ $gc }= {'geneCodon'=> $gc,'codon'=>$feature->{'location'}->{'position'}->{'@position'},'description'=>$feature->{'@type'} . " (" . $feature->{'@description'} . " $effectVal)"  };
										 } else {
											 if (exists($codons->{ $gc }->{'description'})) {
												$codons->{ $gc }->{'description'} .= ", " .  $feature->{'@type'} . " (" . $feature->{'@description'} . ")";
											 } else {
												$codons->{ $gc }->{'description'} =  $feature->{'@type'} . " (" . $feature->{'@description'} . ")";                                                     
											 }
										 }                                 
                                         $codons->{ $gc }->{'gene'}=$new->{'gene'};										 
                                  } 
                             }                                                  
                         }
                         $new->{'uniprotDisease'}=$new->{'gene'} . ". "; 
                         $new->{'uniprotFunction'}=$new->{'gene'} . ". ";
                         my $check="DrugBank";
                         if (exists(  $original->{'entry'}->{'dbReference'})) {
                               if (ref($original->{'entry'}->{'dbReference'}) eq 'ARRAY' ) {
								   foreach my $dbref (@{$original->{'entry'}->{'dbReference'}}) {
									   if ($dbref->{'@type'} eq $check) {
											if (exists($new->{'uniprotNotes'})) {
												$new->{'uniprotNotes'}.= ", $check=" .  $dbref->{'property'}->{'@value'};
											} else {
												 $new->{'uniprotNotes'} = "$check=" .  $dbref->{'property'}->{'@value'};
											}							       
									   }
								   }
							   } else {
							          $dbref=$original->{'entry'}->{'dbReference'}; 
									   if ($dbref->{'@type'}  eq $check) {
											if (exists($new->{'uniprotNotes'})) {
												$new->{'uniprotNotes'}.= ", $check=" .  $dbref->{'property'}->{'@value'};
											} else {
												 $new->{'uniprotNotes'} = "$check=" . $dbref->{'property'}->{'@value'};
											}							       
									   }							   
							   }
                         
                         
                         }
					     if (ref($original->{'entry'}->{'comment'}) eq 'ARRAY' ) {
							   foreach my $comment (@{$original->{'entry'}->{'comment'}}) {
								   if (exists($comment->{'text'}->{'$t'})) {
									   if ($comment->{'@type'} eq "disease") {
										   if (exists($comment->{'disease'}->{'name'})) {
											   $new->{'uniprotDisease'}=$new->{'uniprotDisease'} . "  " . $comment->{'disease'}->{'name'}->{'$t'};
										   }
										   if (exists($comment->{'disease'}->{'description'}->{'$t'} )) {
											   $new->{'uniprotDisease'}=$new->{'uniprotDisease'} . "  " .$comment->{'disease'}->{'description'}->{'$t'};
										   }
										   if (exists($comment->{'text'})) {
											   $new->{'uniprotDisease'}=$new->{'uniprotDisease'} .  $comment->{'text'}->{'$t'};
										   }
									   } 
									   if ($comment->{'@type'} eq "function") {
										   $new->{'uniprotFunction'}=$new->{'uniprotFunction'} . "  " . $comment->{'text'}->{'$t'} . ".";
									   }     
								   }
							   }
							   $new->{'uniprotFunction'}=~s/\.\./\./g;
							   $new->{'uniprotDisease'}=~s/\.\./\./g;
							   
                         }                           
                         foreach my $codon (keys %$codons) {
                               $codons->{$codon}->{'description'}=~s/\( \)//g;
							#	print "here $codon $codons->{$codon}->{'description'}\n";  
							if (defined(                  $codons->{$codon}->{'geneCodon'})) {           
                               #$GENECODON->update({"geneCodon" => $codons->{$codon->{'geneCodon'}}},{'$set'=>$codons->{$codon}},{safe=>1,upsert=>1});
                               $GENECODON->insert($codons->{$codon},{safe=>1,upsert=>1});
                            }
                         }
                         if (defined($new->{'gene'})) {
#                            $GENE->update({"gene" => $new->{'gene'}},{'$set'=>$new},{safe=>1,upsert=>1});                         
                         }
                         print "Gene Added: $new->{'gene'}\n";                                             
                         last LOOP;
                     }
               }
         }        
    }
     $GENECODON->ensure_index({'geneCodon'=>1});
    
     close (FILE);
     
    
}
1;
