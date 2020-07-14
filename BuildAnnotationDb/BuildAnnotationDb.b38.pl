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
use insertDbnsfp35_b38;
use insertGnomad_b38;
use insertDbnsfp40_b38;
$dir="/labs/ngd-data/prodCentralDB/BuildAnnotationDb";
##$hostName='pbc-dcraig-login.tgen.org:27040';
##$hostName='pbc-dcraig-db1:26565';
##$hostName='localhost:27030';
##$hostName='labdb01.tgen.org:27227'; #new build 38 labdb01
$hostName='labdb01.tgen.org:26321'; # build38 development

#insertVcf->loadTo ('/data/db/annotations/temp/dbsnp/20180418/All_20180418.vcf','coord',$hostName,'dbsnp');
#insertVcf->loadTo ('/data/db/annotations/temp/clinvar/clinvar_20191021.vcf','coord',$hostName,'clinvar');
#insertVcf->loadTo ( "$dir/db/clinvar/clinvar_20200622.vcf.gz",'coord',$hostName,'clinvar');
#print "loaded clinvar 2020622\n";
#insertVcf->loadTo('db/1kg/KG3.head.mPASS.vcf','coord',$hostName,'phase3KG');
#insertVcf->loadTo('/data/db/annotations/temp/gnomad/2.0.2/gnomad.genomes.r2.0.2.sites.liftover.b38.autosomes_and_X.vcf','coord',$hostName,'gnomad');
##insertVcf->loadTo('db/ExAc/ExAC.r0.3.sites.vep.vcf','coord',$hostName,'ExAc');
#insertVcf->loadTo('db/gnomad/gnomad_data/vcf/genomes/gnomad.hom.header.vcf','coord',$hostName,'gnomad');
#insertGnomad_b38->loadTo('db/gnomad/b38/gnomad*bgz','coord',$hostName);
# load cosmic
#insertVcf->loadTo('/data/db/annotations/temp/cosmic/90/CosmicCodingMuts.vcf','coord',$hostName,'cosmicCoding');
#insertVcf->loadTo('/data/db/annotations/temp/cosmic/90/CosmicNonCodingVariants.vcf','coord',$hostName,'cosmicNonCoding');
#print "finished loading cosmic\n";
insertDbnsfp40_b38->loadTo ('db/dbnsfp/dbnsfp4.1links1/dbNSFP*.chr*','coord',$hostName);
print "finished loading dbnsfp4.1";
#print "finished loading cosmic v90\n";
#insertDbscSNV1->loadTo ('db/dbnsfp3.4/dbscSNV1.1*','coord',$hostName);
#print "finished loading dbscSNV1\n";
#initMongo->removeEntries('AnnotateBy','gene',$hostName);;
#initMongo->removeEntries('AnnotateBy','geneCodon',$hostName);;
#insertUniprotAnnotate->loadTo('db/uniprot/uniprot_sprot.xml','gene',$hostName);
#print "finished loading uniprot\n";
#insertDiseaseDescription->loadTo('db/genes/descriptions.csv','gene',$hostName,'descriptions');
#print "finished loading gene descriptions\n";
#insertTwoCol->loadTo('db/genes/orphanet.txt','gene',$hostName,'orphanet');
#print "finished loading orphanet\n";
#insertTwoCol->loadTo('db/genes/genecards.txt','gene',$hostName,'genecards');
#print "finished loading genecards\n";
#insertTwoCol->loadTo('db/genes/CTD.parse.txt','gene',$hostName,'diseaseKeywords');
#print "finished loading CTD.parse\n";
#insertDbnsfpGene->loadTo ('db/genes/dbNSFP4.0_gene.complete','gene',$hostName);
#print "finished loading dbNSFP3.5_gene\n";


#insertTwoCol->loadTo('db/genes/go_ontology.txt','gene',$hostName,'goOntology');
#insertTwoCol->loadTo('db/genes/isca.txt','gene',$hostName,'isca');
#insertToPanels->loadTo('db/genes/mito100.txt','gene',$hostName,'mito100');
#insertToPanels->loadTo('db/genes/acmg.txt','gene',$hostName,'acmg');
#insertToPanels->loadTo('db/genes/epilepsy.txt','gene',$hostName,'epilepsy');
#insertToPanels->loadTo('db/genes/gedi.txt','gene',$hostName,'gedi');
#insertToPanels->loadTo('db/genes/mitocarta.txt','gene',$hostName,'mitocarta');
#insertToPanels->loadTo('db/genes/foundation.txt','gene',$hostName,'foundationOne');
#insertToPanels->loadTo('db/genes/hypoplasia.txt','gene',$hostName,'hypoplasia');
#insertToPanels->loadTo('db/genes/illumina.txt','gene',$hostName,'illumina');
#insertToPanels->loadTo('db/genes/ion.txt','gene',$hostName,'ion');
#insertPli->loadTo('db/ExAc/pLI.tsv','gene',$hostName,'pLI');
#insertToPanels->loadTo('db/genes/tcga-gbm.txt','gene',$hostName,'tcga-gbm');
#insertCgd->loadTo ('db/genes/cgd.txt','gene',$hostName);
#insertTcgaMel->loadTo ('db/tcga/PR_TCGA_SKCM_PAIR_Capture_All_Pairs_QCPASS.aggregated.capture.tcga.uuid.somatic.maf.csv','tcga',$hostName);
#insertTwoCol->loadTo('db/genes/level2.txt','gene',$hostName,'level2');
#insertTwoCol->loadTo('db/genes/pharmkgb.txt','gene',$hostName,'pharmkgb');
#insertClinvarNames->loadTo('db/genes/clinvar.txt','gene',$hostName);
print "finished loading genes\n";
#insertToPanels->loadTo('db/genes/level2.txt','gene',$hostName,'level2');
#insertRefseq->loadTo('db/refseq/refSeqSummary.txt','gene',$hostName);
print "inserted everything\n";
