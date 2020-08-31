import pymongo
import pysam
import pandas as pd
import itertools
import time

hostname='labdb01.tgen.org'
port=26321
snv_file='/labs/schork/reference/CADD/whole_genome_SNVs_inclAnno.tsv.gz'
indel_file='/labs/schork/reference/CADD/gnomad.genomes.r3.0.indel_inclAnno.tsv.gz'
dtypes_file='/labs/schork/reference/CADD/CADD_dtypes.csv'
bed_file='/labs/schork/projects/USCORIEN/capture_kit/SeqCap_EZ_Exome_v3_hg38_primary_targets_padded100.bed.gz'
start_reg={'Chr':'chr14','Start':106874451,'End':106875013}

def line2update(line,names,int_idx,float_idx):
    coord=':'.join( ['chr' + line[0],line[1],line[2],line[3]])
    #print(coord)
    include=list(itertools.repeat(False,4)) + list(map(lambda x: (x!='NA') & (x!='.') & (x!='-'),line[4:]))
    include_int=list(map(lambda x,y: x and y,int_idx,include))
    int_data=list(map(int,itertools.compress(line,include_int)))
    int_keys=list(itertools.compress(names,include_int))
    include_float=list(map(lambda x,y: x and y,float_idx,include))
    float_data=list(map(float,itertools.compress(line,include_float)))
    float_keys=list(itertools.compress(names,include_float))
    include_str=list(map(lambda x,y,z: x and not y and not z,include,int_idx,float_idx))
    str_data=list(itertools.compress(line,include_str))
    str_keys=list(itertools.compress(names,include_str))
    fields=int_keys+float_keys+str_keys
    keys=list(map(lambda x: 'CADD_' + x,fields))
    records=dict(zip(keys,int_data+float_data+str_data))
    return pymongo.UpdateOne({'coord':coord},{"$set":records},upsert=True)

client=pymongo.MongoClient(hostname,port)
ab_coord=client['AnnotateBy']['coord']

fields_table=pd.read_csv(dtypes_file,index_col=0,names=['dtype'])
names=list(fields_table.index.values)
int_idx=list(map(lambda x: x=='int',fields_table.dtype.values))
float_idx=list(map(lambda x: x=='float',fields_table.dtype.values))

snv_tbx=pysam.TabixFile(snv_file)
indel_tbx=pysam.TabixFile(indel_file)
contigs=list(map(lambda x: 'chr' + x, snv_tbx.contigs))

regions=pd.read_table(bed_file,names=['Chr','Start','End'])
regions=regions.loc[regions.Chr.isin(contigs)]
regions['length']=regions.End-regions.Start
regions['cumLength']=regions['length'].cumsum()
start_idx=regions.loc[(regions['Chr']==start_reg['Chr']) & (regions['Start']==start_reg['Start']) & (regions['End']==start_reg['End'])].index
regions=regions.iloc[start_idx[0]:]

for idx,reg in regions.iterrows():
    print('region ' + reg.Chr + ':' + str(reg.Start) + '-' + str(reg.End))
    updates=[]
    for line in snv_tbx.fetch(reg.Chr.replace('chr',''),reg.Start,reg.End,parser=pysam.asTuple()):
        updates.append(line2update(line,names,int_idx,float_idx))
    for line in indel_tbx.fetch(reg.Chr.replace('chr',''),reg.Start,reg.End,parser=pysam.asTuple()):
        updates.append(line2update(line,names,int_idx,float_idx))
    result=ab_coord.bulk_write(updates)
    print('modified: ' + str(result.modified_count) + ', upserted: ' + str(result.upserted_count))
    print(time.asctime() + ': ' + str(100*reg.cumLength/regions.cumLength.iloc[-1]) + '% complete',flush=True)
