from pymongo import MongoClient
import file_management
import os, fnmatch

annotate_client = MongoClient("labdb01.tgen.org", 26321)
annoDb = annotate_client["AnnotateBy"]
coordCollection = annoDb['coord']
testCollection = annoDb['test']
geneCollection = annoDb['gene']
pharmgkb_dir = "/Users/tizatt/Desktop/pharmvar.pharmgkb"
drug_label_dir = "/drugLabels"
search_pattern = "*.tsv"

coord_files = {"clinical_ann_metadata.tsv": "Location", "study_parameters": "", "variant_drug_ann.tsv": "Variant",
              "var_fa_ann.tsv": "Variant", "var_pheno_ann.tsv": "Variant", "automated_annotations.tsv": "Variation Name",
              "clinicalVariants.tsv": "variant"}

coord_short_names = {"clinical_ann_metadata.tsv": "cl_ann", "study_parameters": "", "variant_drug_ann.tsv": "var_drug",
              "var_fa_ann.tsv": "var_fa", "var_pheno_ann.tsv": "var_pheno", "automated_annotations.tsv": "auto",
              "clinicalVariants.tsv": "clin_var"}
#"occurrences.tsv": "", "relationships.tsv": ""

gene_files = {"drugLabels.byGene.tsv":"", "drugLabels.tsv":""}
file_types = ["coord", "gene"]
parse = file_management.parser()


def load_to_coord(coord_file, name):
    counter = 0
    bulk = coordCollection.initialize_unordered_bulk_op()
    search_id = coord_files[name]
    short_name = coord_short_names[name]
    for record in coord_file:
        if search_id in record:

            if record[search_id].startswith("rs"):
                rs_id = record[search_id].replace("rs", "")

                del record[search_id]
                if "Gene" in record:
                    del record["Gene"]
                if "Gene Symbols" in record:
                    del record["Gene Symbols"]
                if "gene" in record:
                    del record["gene"]
                if "Chromosome" in record:
                    del record["Chromosome"]

                #testCollection.update_one({"RS":rs_id}, {'$push': {"pharmgkb." + short_name: record}})
                bulk.find({'RS': rs_id}).update(
                    {'$push': {"pharmgkb." + short_name: record}})
                counter += 1

                if counter % 500 == 0:
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


def get_files(type_db, d):
    files = []
    for d_name, sd_name, f_list in os.walk(d):
        for file_name in f_list:
            if fnmatch.fnmatch(file_name, search_pattern):
                if type_db == "coord" and file_name in coord_files:
                    files.append(os.path.join(d_name, file_name))
                elif type_db == "gene" and file_name in gene_files:
                    files.append(os.path.join(d_name, file_name))
    return files


coord_file_names = get_files("coord", pharmgkb_dir)

#for file in coord_file_names:
#    filename = os.path.basename(file)
#    file_dict = file_management.load_file(file)
#    print (filename)
#    load_to_coord(file_dict, filename)

drug_labels = get_files("gene", pharmgkb_dir+drug_label_dir)
load_drug_labels(drug_labels)

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
