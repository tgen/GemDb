import fnmatch
import os
import csv
import sys
import optparse
import gzip
from pymongo import MongoClient

pharmgkb_yaml = "pharmgkb.yaml"
annotate_client = MongoClient("labdb01.tgen.org", 26321)
annoDb = annotate_client["AnnotateBy"]
coordCollection = annoDb['coord']
testCollection = annoDb['test']
geneCollection = annoDb['gene']
cosmic_dir = "/data/db/annotations/cosmic/"
cmc_file = cosmic_dir + "cmc.v92.tsv"


def load_to_coord(coord_file):
    counter = 0

    for record in coord_file:
        if record.startswith("chr"):
            coord = record.rstrip()

        coordCollection.find_one_and_update(
            {'coord': coord},
            {"$set":
                 {"Cancer Mutation Census": "yes"}
             }, upsert=True)

        counter += 1
        if counter % 500 == 0:
            print("Counter at " + str(counter)+ "; coord = "+coord)


def load_tsv_file(filename):
    if filename.endswith("tsv"):
        input = open(filename, "rt")
        tsv_file = input.readlines()
        return tsv_file


cmc_list = load_tsv_file(cmc_file)
load_to_coord(cmc_list)


