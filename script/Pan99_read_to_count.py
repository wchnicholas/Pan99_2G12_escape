#!/usr/bin/python
import os
import sys
import glob
import string
from Bio import SeqIO
from collections import Counter

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def readingFile2SampleID(file2SampleID_file):
  file2SampleID_dict = {}
  infile = open(file2SampleID_file,'r')
  for line in infile.xreadlines():
    if 'R1File' in line: continue
    line = line.rstrip().rsplit("\t")
    file2SampleID_dict[line[0]] = line[1]
  infile.close()
  return file2SampleID_dict

def readingBarcode2Resi(Barcode2Resi_file):
  Barcode2Resi_dict = {}
  infile = open(Barcode2Resi_file,'r')
  for line in infile.xreadlines():
    if 'Barcode' in line: continue
    line = line.rstrip().rsplit("\t")
    Barcode2Resi_dict[line[0]] = line[1]
  infile.close()
  return Barcode2Resi_dict

def readingWTcodon(WTcodon_file):
  WTcodon_dict = {}
  infile = open(WTcodon_file,'r')
  for line in infile.xreadlines():
    if 'Codon' in line: continue
    line = line.rstrip().rsplit("\t")
    WTcodon_dict[line[0]] = line[1]
  infile.close()
  return WTcodon_dict

def Processlib(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file, "fastq")
  R2records = SeqIO.parse(R2file,"fastq") 
  muts = []
  Primerlength = 23 
  R1frameshift = 1
  R2frameshift = 2
  count_record = 0
  for R1record in R1records:
    count_record += 1
    #if count_record == 1000: break
    R2record = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi  = R1seq[Primerlength+R1frameshift::]
    R2roi  = R2seq[Primerlength+R2frameshift::]
    R1overlap = R1seq[72:72+201]
    R2overlap = R2seq[85:85+201]
    if R1overlap != rc(R2overlap): continue
    if 'N' in R1roi or 'N' in R2roi: continue
    resi144 = translation(R1roi[0:3])
    resi145 = translation(R1roi[3:6])
    resi160 = translation(R1roi[48:51])
    resi172 = translation(R1roi[84:87])
    resi192 = translation(R1roi[144:147])
    resi196 = translation(R1roi[156:159])
    resi226 = translation(R1roi[246:249])
    resi246 = translation(rc(R2roi[0:3]))
    mut = resi144+resi145+resi160+resi172+resi192+resi196+resi226+resi246
    muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def Output(mutcount_dict, outfile):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['sample','mut','count'])+"\n")
  muts = [mut for sample in mutcount_dict.keys() for mut in mutcount_dict[sample].keys()]
  for mut in muts:
    for sample in mutcount_dict.keys():
      mut_count = mutcount_dict[sample][mut] if mut in mutcount_dict[sample].keys() else 0
      outfile.write("\t".join(map(str,[sample,mut,mut_count]))+"\n")
  outfile.close()

def main():
  R1filenames = glob.glob('fastq/*_R1_*.fastq')
  outfile  = 'result/Pan99_mut_count.tsv'
  print "Analyzing %i samples" % len(R1filenames)
  mutcount_dict = {}
  for R1filename in R1filenames:
    SampleID = R1filename.rsplit('Pan99lib_')[1].rsplit('_R1_')[0]
    print "Analyzing %s (%s)" % (R1filename, SampleID)
    mutcount_dict[SampleID] = Processlib(R1filename)
  Output(mutcount_dict, outfile)

if __name__ == "__main__":
  main()
