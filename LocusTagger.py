import argparse
import sys
import Bio.Seq
import pandas as pd
import numpy as np
import re
import os
from Bio import SeqIO
import gffutils
from Bio.Seq import Seq
import pyfaidx


parser = argparse.ArgumentParser(
    description='Use contig name and site to locate locus tag to which it belongs, and the sequence if provided.',
    add_help=False,
    usage='\nLocusTagger --gff [input .gff file] --fna [input .fna file] --contig [str] --site [integer] --seq [optional, if you want genomic sequence of your locus.]')

parser.add_argument(
    '--gff',
    metavar='[input.gff]',
    required=True,
    type=str)
parser.add_argument(
    '--fna',
    metavar='[input.fna]',
    required=True,
    type=str)
parser.add_argument(
    '--site',
    metavar='[Asked site]',
    required=True,
    type=int)
parser.add_argument(
    '--contig',
    metavar='[Asked contig.]',
    required=True,
    type=str)
parser.add_argument(
    '--seq',
    action = 'store_true',
    help='Provide genomic sequence of your locus.',
    required=False)

args = parser.parse_args()

if(os.path.exists('cds.fna')):
    os.remove('cds.fna')
if(os.path.exists('CDS.info')):
    os.remove('CDS.info')

gff = args.gff
fna = pyfaidx.Fasta(args.fna)

gffdb = gffutils.create_db(
    gff,dbfn='gff.db',
    force=True,
    keep_order=True,
    merge_strategy="create_unique",
    sort_attribute_values=True,
    )

for i in gffdb.features_of_type('CDS', order_by='seqid'):
    g_seq = (i.sequence(fna))
    g_seq = Seq(g_seq)
    if i.strand == '-':
        g_seq = g_seq.reverse_complement()
    f = open('cds.fna','a')
    f.write(str('>'+i.seqid+"_"+i.id.replace('cds-',''))+'\n'+str(g_seq)+'\n')
    f.close()

for l in gffdb.features_of_type('CDS', order_by='seqid'):
    cdsif = open('CDS.info','a')
    cdsif.write(str(l.id.replace('cds-', ''))+"\t"+str(l[0]) + "\t" + str(l[3]) + "\t" + str(l[4]) +'\n')
    cdsif.close()

cdst = pd.read_table(r'CDS.info',header=None,names=['Name','Contig' ,'start' ,'end'])
cdst.to_csv('CDS.info',sep='\t',index=False)

c=pd.read_csv('CDS.info',sep='\t')

cds_dict = c.set_index('Name').to_dict('dict')

def gene_select(inputcontig,inputsite):
    keycl = []
    for keyc in (cds_dict['Contig']):
        if (cds_dict['Contig'])[keyc] == inputcontig:
            keycl.append(keyc)

    keysl = []
    for keys in (cds_dict['start']):
        if int((cds_dict['start'])[keys]) <= int(inputsite):
            keysl.append(keys)

    keyel = []
    for keye in (cds_dict['end']):
        if int((cds_dict['end'])[keye]) >= int(inputsite):
            keyel.append(keye)

    set_keycl = set(keycl)
    set_keyel = set(keyel)
    set_keysl = set(keysl)
    gene = str(inputcontig) + '_' + (''.join(set_keycl & set_keyel & set_keysl))
    print('>'+gene)
    for contig in SeqIO.parse('cds.fna', "fasta"):
        if contig.id == gene:
            print(Bio.Seq.translate(contig.seq, table=11))

gene_select(args.contig,args.site)






