import csv
import optparse


def parser():
    parser = optparse.OptionParser()

    parser.add_option('--file',
                      action="store", dest="file",
                      help="pharmagkb file", default="test.tsv")
    return parser


def load_file(filename):
    if filename.endswith("tsv"):
        tsv_file = []
        with open(filename) as file:
            reader = csv.DictReader(file, dialect='excel-tab')
            for row in reader:
                rowDict = dict(row)
                filtered = dict((k, v) for k, v in rowDict.items() if k is not None and k != '')
                filtered2 = dict((k, v) for k, v in filtered.items() if v is not None and v != '')
                tsv_file.append(filtered2)
        return tsv_file
