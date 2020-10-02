import csv
import sys
import optparse
import gzip
maxInt = sys.maxsize


def parser():
    parser = optparse.OptionParser()

    parser.add_option('--file',
                      action="store", dest="file",
                      help="tsv file", default="test.tsv")
    return parser


def load_tsv_file(filename):
    if filename.endswith("tsv") or filename.endswith("maf.gz"):
        tsv_file = []

        if filename.endswith("gz"):
            input = gzip.GzipFile(filename, "rb")
        else:
            input = open(filename, "rb")
        csv.field_size_limit(sys.maxsize)
        reader = csv.DictReader(input, dialect='excel-tab')
        for row in reader:
            row_dict = dict(row)
            filtered = dict((k, v) for k, v in row_dict.items() if k is not None and k != '')
            filtered2 = dict((k, v) for k, v in filtered.items() if v is not None and v != '')
            tsv_file.append(filtered2)
    return tsv_file


def delete_fields(record, search_id):
    del record[search_id]
    if "Gene" in record:
        del record["Gene"]
    if "Gene Symbols" in record:
        del record["Gene Symbols"]
    if "gene" in record:
        del record["gene"]
    if "Chromosome" in record:
        del record["Chromosome"]

    return record