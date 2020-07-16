#!/usr/bin/perl 
use lib "$FindBin::Bin/lib/";  # use the parent directory
use POSIX;
use updateGeneLocat;
use FindBin;                 # locate this script
use MongoDB;
use initMongo;
use insertClinvarNames;
use insertRefseq;
use insertInteractions;
use insertVcf;
use insertTcgaMel;
use insertDbnsfpGene;
use insertCgd;
use insertToPanels;
use insertTwoCol;
use insertDbnsfp;
use insertGedi;
use insertLevel2;
use insertPli;
use insertDbscSNV1;
use insertUniprotAnnotate;
use insertDbnsfpGnomad;
use insertDiseaseDescription;
use insertGnomad_b38;
$dir=""#[directory where your files are;
$hostName=#[your host name];  Ex. "localhost:25555;"

#### Examples #####


### VCF source file
##
#insertVcf->loadTo ('$dir/dbsnp/20180418/All_20180418.vcf','coord',$hostName,'dbsnp');
#insertVcf->loadTo ('$dir/clinvar/clinvar_20191021.vcf','coord',$hostName,'clinvar');


### dbNSFP 
##   	This module inserts a folder with dbNSFP files.  
##	There is a configuration file that must be used to specify which fields you want inserted.
##
#insertDbnsfp40_b38->loadTo ('$dir/dbNSFP*.chr*','coord',$hostName);


#### Remove collections.  
##    	*Make sure you have a backup of the db before doing this, 
##       as it permanently removes the collection
##
#initMongo->removeEntries('AnnotateBy','gene',$hostName);;
#initMongo->removeEntries('AnnotateBy','geneCodon',$hostName);;


#### Gene Annotations
##
#insertUniprotAnnotate->loadTo('$dir/uniprot/uniprot_sprot.xml','gene',$hostName);
#insertDiseaseDescription->loadTo('$dir/genes/descriptions.csv','gene',$hostName,'descriptions');
#insertTwoCol->loadTo('$dir/genes/orphanet.txt','gene',$hostName,'orphanet');
#insertDbnsfpGene->loadTo ('$dir/genes/dbNSFP4.0_gene.complete','gene',$hostName);
#insertToPanels->loadTo('$dir/genes/acmg.txt','gene',$hostName,'acmg');
#insertPli->loadTo('$dir/ExAc/pLI.tsv','gene',$hostName,'pLI');
#insertClinvarNames->loadTo('$dir/genes/clinvar.txt','gene',$hostName);
