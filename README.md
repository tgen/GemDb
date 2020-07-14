# GemDb

A mongo database for storing and annotating VCF files.  VCF's from the Phoenix pipeline, or any other pipeline can be inserted and annotated with latest Dbnfsp, Cosmic, Clinvar, etc...

## What is GemDb?
A mongo database containing DNA/RNA events and public annotations. 

## How do I use GemDb?
*To use the existing TGen GemDb instance (containing samples from various TGen projects (germline and somatic), please follow the wiki(https://github.com/tgen/GemDb/wiki).  This will guide you through uploading your files, filtering, annotations and how to query them or extract them into excel files or reports.

## Setting up your own GemDb instance
If you would like to setup your own instance of GemDb - follow these instructions below:

1. Create a new mongo connection on your local server. Note the server and port that you created (localhost:28579) - this is your *$hostname*.

2. Load Annotations:
      
      Use the existing Annotate db in the main TGen GemDb.  
          
          - Go to isilon directory (/labs/ngd-data/prodCentralDB/backups)
            From there 'cd' into 'GemDb37' or 'GemDb38', depending on the genome build version you plan on working with.
            Then 'mongorestore' the 'AnnotateBy' folder into your mongo connection from step 1.  Be sure you are restoring the 'AnnotateBy' folder ONLY.
            
      To add additional annotations to your copy of the main TGen GemDb AnnotateBy:
          
          - Clone this git repo 
          
          - Go to the "BuildAnnotationDb" sub-folder and open the "BuildAnnotationDb*.pl" file of your choice. 
            Within this file, change the 'hostname' variable (server:port) to your server and port from step 1.  (Ex. localhost:28972)
            Following the code within "BuildAnnotationDb*.pl" - add your own annotation files by downloading them locally and adding them to this script.
            Then run BuildAnnotationDb*.pl, noting you must run this with isilon mounted (and check that the $dir variable points to the correct location where your files are). 
            
          - To add new drug rules - Go to the "rulesDir" sub-folder.  This contains drug-gene rules that can be added to the db.  Use the existing drug rule files as examples for adding your own drug rules.
            Follow the README in that folder for adding drug rules to your Annotation db.
            Make sure to change the '$hostname' within 'BuildRulesDb.pl' to the '$hostname' from Step 1. 
      
3. Load variants (VCF files) into your database
         
         - Clone the appropriate genome build version below into your local server and follow the instructions from there.  
           Remember your '$hostname' from Step 1 as it will be where variants are loaded as well.

Build 37 - https://github.com/tgen/GemDb37

Build 38 - https://github.com/tgen/GemDb38

## Wiki
https://github.com/tgen/GemDb/wiki

*contact tizatt@tgen.org for running/installing.
