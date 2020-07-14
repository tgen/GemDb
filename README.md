# What is GemDb

A mongo database for storing and annotating VCF files.  VCF's from the Phoenix pipeline, or any other pipeline can be inserted and annotated with latest Dbnfsp, Cosmic, Clinvar, etc...

## How do I use GemDb?
*To use the existing TGen GemDb instance (containing samples from various TGen projects (germline and somatic), please follow the wiki(https://github.com/tgen/GemDb/wiki).  This will guide you through uploading your files, filtering, annotations and how to query them or extract them into excel files or reports.

## Setting up your own GemDb instance
If you would like to setup your own instance of GemDb - follow these instructions below:

1. Install MongoDb
           
           - Download here: https://www.mongodb.com/try/download/community
           - Follow instructions for installing it on your server, and then start up a mongo instance - specifying a specific port (ex. localhost:28579)
             This is your '$hostname'.

2. Load Annotations:
      
      Use the existing Annotate db in the main TGen GemDb.  
          
          - cd /labs/ngd-data/prodCentralDB/backups
          
            [Using the $hostname from Step 1]
            
            mongorestore --host $hostname --db AnnotateBy GemDb37/AnnotateBy (for Build37)
            
            OR 
            
            mongorestore --host $hostname --db AnnotateBy GemDb38/AnnotateBy (for Build38)
            
            
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
