import fnmatch
import os
from pymongo import MongoClient
from lib.pharmgkb import file_management

pharmgkb_yaml = "pharmgkb.yaml"
annotate_client = MongoClient("dback-optane1.tgen.org", 28333)
annoDb = annotate_client["AnnotateBy"]
coordCollection = annoDb['coord']
testCollection = annoDb['test']
geneCollection = annoDb['gene']
pharmgkb_dir = "/labs/ngd-data/prodCentralDB/BuildAnnotationDb/db/pharmvar.pharmgkb/"
drug_label_dir = "/drugLabels"
pharmvar_file = pharmgkb_dir + "pharmvar_GRCh38_haplotypes_combined.tsv"

search_pattern = "*.tsv"
coord_files_keys = {"clinical_ann_metadata.tsv": "Location", "variant_drug_ann.tsv": "Variant",
                    "var_fa_ann.tsv": "Variant",
                    "var_pheno_ann.tsv": "Variant", "automated_annotations.tsv": "Variation Name",
                    "clinicalVariants.tsv": "variant"}

coord_short_names = {"clinical_ann_metadata.tsv": "cl_ann", "study_parameters": "", "variant_drug_ann.tsv": "var_drug",
                     "var_fa_ann.tsv": "var_fa", "var_pheno_ann.tsv": "var_pheno", "automated_annotations.tsv": "auto",
                     "clinicalVariants.tsv": "clin_var"}

gene_drug_files_keys = {"drugLabels.byGene.tsv": "", "drugLabels.tsv": ""}
cpic_pairs_keys = {"cpicPairs.tsv": "Gene"}
gene_files_keys = {"genes.tsv": "Symbol"}

file_types = ["coord", "gene"]
parse = file_management.parser()


def load_to_coord(coord_file, name):
    counter = 0
    bulk = coordCollection.initialize_unordered_bulk_op()
    search_id = coord_files_keys[name]
    short_name = coord_short_names[name]
    for record in coord_file:
        if search_id in record:
            rs_id = -1
            if record[search_id].startswith("rs"):

                rs_id = int(record[search_id].replace("rs", ""))
            else:  ### This will take some effort to do the star alleles ***
                haplotypes = record[search_id].split(",")
                for ht in haplotypes:
                    if "xN" in ht:
                        print("do something with the multiples")
                    else:
                        # Exact match
                        for pvr in pharmvar_file:
                            if record[search_id] == pvr["Haplotype Name"]:
                                rs_id = int(pvr["rsID"].replace("rs", ""))

            file_management.delete_fields(record, search_id)

            # coordCollection.update({"RS":rs_id}, {'$push': {"pharmgkb." + short_name: record}})
            bulk.find({'RS': rs_id}).update(
                {'$push': {"pharmgkb." + short_name: record}})
            counter += 1

            if counter % 500 == 0:
                print("Counter at " + str(counter))
                bulk.execute()
                bulk = coordCollection.initialize_ordered_bulk_op()

    if counter % 500 != 0:
        bulk.execute()


def load_drug_labels(drug_labels):
    label_by_gene_file = ""
    label_file = ""
    counter = 0
    bulk = geneCollection.initialize_unordered_bulk_op()
    for filepath in drug_labels:
        print(os.path.basename(filepath))
        if "drugLabels" in os.path.basename(filepath):
            if os.path.basename(filepath) == "drugLabels.byGene.tsv":
                label_by_gene_file = file_management.load_file(filepath)
            else:
                label_file = file_management.load_file(filepath)

    for drug_gene_label in label_by_gene_file:
        gene_name = drug_gene_label["Gene Symbol"]
        label_ids = drug_gene_label["Label IDs"]
        drugs = []
        for label_id in label_ids.split(";"):
            for label in label_file:
                if label_id in label["PharmGKB ID"]:
                    drugs.append(label)

        bulk.find({"gene": gene_name}).update(
            {'$set': {"pharmgkb.drug": drugs}})
        counter += 1

        if counter % 500 == 0:
            bulk.execute()
            bulk = geneCollection.initialize_ordered_bulk_op()

    if counter % 500 != 0:
        bulk.execute()


def load_cpic_labels(cpic_labels, name):
    counter = 0
    bulk = geneCollection.initialize_unordered_bulk_op()
    for cpicPair in cpic_labels:
        search_id = cpic_pairs_keys[name]
        gene_name = cpicPair[search_id]

        bulk.find({"gene": gene_name}).update(
            {'$push': {"cpic": cpicPair}}
        )
        counter += 1
        print(counter)
        if counter % 500 == 0:
            bulk.execute()
            bulk = geneCollection.initialize_ordered_bulk_op()

    if counter % 500 != 0:
        bulk.execute()


def load_to_gene(gene_file, name):
    counter = 0
    bulk = geneCollection.initialize_unordered_bulk_op()

    for gene_record in gene_file:
        search_id = gene_files_keys[name]
        gene_name = gene_record[search_id]
        gene_fields = {"pharmgkb.IsVIP": gene_record["Is VIP"],
                       "pharmgkb.Has_Variant_Annotation": gene_record["Has Variant Annotation"]}

        bulk.find({"gene": gene_name}).update(
            {'$set': gene_fields})
        counter += 1

        if counter % 500 == 0:
            bulk.execute()
            bulk = geneCollection.initialize_ordered_bulk_op()

    if counter % 500 != 0:
        bulk.execute()


def get_files(type_db, d):
    files = []
    for d_name, sd_name, f_list in os.walk(d):
        for file_name in f_list:
            if fnmatch.fnmatch(file_name, search_pattern):
                if type_db == "coord" and file_name in coord_files_keys:
                    files.append(os.path.join(d_name, file_name))
                elif type_db == "gene" and file_name in gene_files_keys:
                    files.append(os.path.join(d_name, file_name))
    return files


coord_file_names = get_files("coord", pharmgkb_dir)
gene_file_names = get_files("gene", pharmgkb_dir)

# for file in coord_file_names:
#    filename = os.path.basename(file)
#    file_dict = file_management.load_file(file)
#    print (filename)
#    load_to_coord(file_dict, filename)

#for file in gene_files_keys:
#    filename = os.path.basename(file)
#    file_dict = file_management.load_file(file)
#    load_to_gene(file_dict, filename)

gene_file = pharmgkb_dir + "genes/genes.tsv"
gene_dict = file_management.load_file(gene_file)
load_to_gene(gene_dict, "genes.tsv")

#cpic_file = pharmgkb_dir + "genes/cpicPairs.tsv"
#cpic_dict = file_management.load_file(cpic_file)
#load_cpic_labels(cpic_dict, "cpicPairs.tsv")

# labels = get_files("gene", pharmgkb_dir+drug_label_dir)
# load_by_gene(labels)

#
# Load pharmgkb files
# Send out to modules to process them
# filetypes =
#     pharmvar vcf files
#          CYP*, DPYD, NUDT15
#
#     Annotations:
#            clinical_ann_metadata.tsv (contains RS ids for matching)
#            clinical_ann.tsv (contains actual annotations)
#
#     Variant Annotations:
#            var_drug_ann.tsv
#            var_fa_ann.tsv
#            var_pheno_ann.tsv
#
#     Automated Annotations:
#            automated_annotations.tsv  (Chemical ID, Variation ID, Variation Name)
#
#     Clinical Variants:
#            clinicalVariants.tsv (variant)
#
#     Dosing Guidelines:
#             json files.
#
#     Drug Labels:
#             Should be inserted in the *gene collection*
#             Where the key is "Gene Symbol" from the drugLabels.byGene, and drugLabels give the drug name
#             drugLabels.tsv (pharmagkb id)
#             drugLabels.byGene.tsv (Gene Symbol)
#
#     Ocurrences:
#             occurrences.tsv (Object ID)
#             
#     Relationships:
#             relationships.tsv ( Entity1_id -> Entity2_id mapping)
#             This file reveals the pharmagkb ids for each gene, rs id, drug, CY ids
#
