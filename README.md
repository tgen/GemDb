# GemDb

A mongo database for storing RNA/DNA variant data from the sequencing pipeline.

contact tizatt@tgen.org for running/installing.


## What is GemDb?
	A mongo database containing variant data (from VCF files) from genomic (RNA/DNA) samples.  Each row in the database represents a variant from one of the various vcf files for a sample.

## How to Insert Data into GemDB
	As a user there is a dropbox on isilon - /ngd-data/VCF/ - where you can place your vcf files.  The database supports various vcf files from the sequencing pipeline, both RNA and DNA.
	The pipeline has an option to automatically copy vcfs over to this folder - contact jetstream@tgen.org for how to enable this. 
	Otherwise follow these steps to manually create the folders and copy over your files.
		
	1. Create a 'project' level directory of your choosing (i.e. on unix -> mkdir /ngd-data/VCF/myProject).  
	   This is usually the name of the study you're working on (i.e. NMTRC, SU2C, etc..)
		
	2. Create 'sample' directories under the project you just made for each individual sample you have.  
	   i.e. /ngd-data/VCF/myProject/mySample1, /ngd-data/VCF/myProject/mySample2, etc...
		
	3. VCFs ->
		- Germline 	
				Create a directory named 'germline' under your sample(s) directory.  Then place the germline vcfs for the respective sample in this directory.  
				i.e. /ngd-data/VCF/myProject/mySample1/germline/example.HC_All.snpEff.vcf
		- Tumor
				Create a directory named for the tumor type or whatever you desire.  i.e. 'breast', 'lung', 'prostate'.  Then place the somatic vcfs in this directory (merged vcf, thf, cna.seg.vcf, etc...)
				i.e. /ngd-dat/VCF/myProject/mySample1/germline/example.cna.seg.vcf
		
	4. You're done!  From here the database scripts will pickup these files and place them thru a set of filters (see 'Filtering' section below) before inserting the PASS/LowQC variants into the database.   


## Filtering

	By default all variants are marked as "PASS" from a VCF file.  The other two options are "LowQC" and "FAIL".  Below is an example of the set of filters (look at conf.pm for actual values) that are applied to each row in a VCF file.  Fields from the vcf file specified are used as filters.  PASS and LowQC variants are inserted into the database, FAIL variants are excluded.

#### Example Filters:

   	**Structural Variants** 
		LowQC: 
			10 < Qual < 25
		FAIL: 
			Qual <= 10 OR
			Gap <= 1500 OR
			NormalAlleleRatio >= 0.01
	
	**CNV**
		LowQC: 
			-1 < Log2FC < 1
	
	**LumosVar CNV**
		LowQC:
			-1 < Log2FC < 1 OR
			SVLEN > 2,5000,000
		
	**Merged VCF**
		LowQC: 
			CALLERS_COUNT <= 1
    
	**GATK Haplo Caller**
		LowQC:
			Qual <= 500
		FAIL:
			Qual <= 300

	**Seurat Point Mutation Caller**
		LowQC:
			Qual <= 20 OR
                        AR1 >= 0.02 OR
                        AR2 <= 0.05 OR
                        ARDIFF <= 0.1
		FAIL:
			Qual < 20 OR
			AR1 >= 0.03 OR
			AR2 <= 0.03

	**LumosVar Point Mutation Caller**
		LowQC:
			filter != SomaticPASS
		FAIL:
			filter == REJECT OR
			filter == GermlineHetPASS OR
			filter == GermlineHetLowQC OR
			filter == GermlineHomPASS OR
			filter == GermlineHomLowQC
					 
	**Cuffdiff**
		LowQC:
			QVALUE > 0.05 OR
			SIGNIFICANT != 'yes' OR
			-2 < LOG2FC < 2
		FAIL:
			-1.9 < LOG2FC < 1.9 OR
			PVALUE > 0.01 OR
			LOG2FC == (-)inf
	
	**Deseq**
		LowQC:
			PADJ > 0.05 OR
			PADJ == 'NA' OR
			-2 < LOG2FC < 2
		FAIL:
			PVALUE > 0.01 OR
			PVALUE == NA OR
			-1.9 < LOG2FC < 1.9 OR
			LOG2FC == (-)inf

	**Tophat Fusion**
		LowQC:
			Qual < 100
		FAIL:
			Qual <= 10
     

