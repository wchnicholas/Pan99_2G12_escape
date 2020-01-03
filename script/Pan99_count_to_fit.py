#!/usr/bin/python
import os
import sys
import glob
import itertools
from collections import Counter, defaultdict

def reading_mut_dict(infile):
  infile = open(infile, 'r')
  mut_dict = {}
  mut_list = []
  for line in infile.xreadlines():
    if 'resi' in line: continue
    line = line.rstrip().rsplit("\t")
    mut_dict[int(line[0])] = line[1].rsplit(',')
    mut_list.append(line[1].rsplit(','))
  infile.close()
  mut_list  = [''.join(mut) for mut in list(itertools.product(*mut_list))]
  resi_list = sorted(mut_dict.keys(), key=lambda x:int(x))
  return mut_dict, resi_list, mut_list

def reading_count_dict(file_count):
  infile = open(file_count, 'r')
  count_dict = defaultdict(dict)
  for line in infile.xreadlines():
    if 'sample' in line: continue
    line = line.rstrip().rsplit("\t")
    sample = line[0]
    mut    = line[1]
    count  = line[2]
    count_dict[sample][mut] = count
  infile.close()
  return count_dict
  
def calculate_fit(mut_sel_count, mut_ipt_count, WT_sel_count, WT_ipt_count):
  enrich_mut = ((float(mut_sel_count)+1)/(float(mut_ipt_count)+1))
  enrich_WT  = ((float(WT_sel_count)+1)/(float(WT_ipt_count)+1))
  return (enrich_mut/enrich_WT)

def call_mut(mut, WT, resi_list):
  name = []
  for mutaa, WTaa, resi in zip(mut, WT, resi_list):
    if mutaa!=WTaa:
      name.append(WTaa+str(resi)+mutaa) 
  mut_num = len(name)
  name = 'WT' if len(name) == 0 else '-'.join(name)
  return name, mut_num

def count_to_fit(outfile, count_dict, mut_dict, mut_list, WT, resi_list):
  print "writing: %s" % outfile
  outfile = open(outfile,'w')
  header = "\t".join(map(str,['mut','name','num','input_count',
                              'count_noAb_rep1','count_noAb_rep2','count_2G12_rep1','count_2G12_rep2',
                              'fit_noAb_rep1','fit_noAb_rep2','fit_2G12_rep1','fit_2G12_rep2','fit_noAb','fit_2G12']))
  outfile.write(header+"\n")
  for mut in mut_list:
    name, mut_num = call_mut(mut, WT, resi_list)
    mut_input_count = count_dict['input'][mut]
    mut_noAb_rep1_count = count_dict['noAb_rep1'][mut]
    mut_noAb_rep2_count =  count_dict['noAb_rep2'][mut]
    mut_2G12_rep1_count = count_dict['2G12_rep1'][mut] 
    mut_2G12_rep2_count = count_dict['2G12_rep2'][mut]
    WT_input_count = count_dict['input'][WT]
    WT_noAb_rep1_count = count_dict['noAb_rep1'][WT]
    WT_noAb_rep2_count =  count_dict['noAb_rep2'][WT]
    WT_2G12_rep1_count = count_dict['2G12_rep1'][WT] 
    WT_2G12_rep2_count = count_dict['2G12_rep2'][WT]
    fit_noAb_rep1   = calculate_fit(mut_noAb_rep1_count, mut_input_count, WT_noAb_rep1_count, WT_input_count)
    fit_noAb_rep2   = calculate_fit(mut_noAb_rep2_count, mut_input_count, WT_noAb_rep2_count, WT_input_count)
    fit_2G12_rep1   = calculate_fit(mut_2G12_rep1_count, mut_input_count, WT_2G12_rep1_count, WT_input_count)
    fit_2G12_rep2   = calculate_fit(mut_2G12_rep2_count, mut_input_count, WT_2G12_rep2_count, WT_input_count)
    fit_noAb = (fit_noAb_rep1 + fit_noAb_rep2)/2
    fit_2G12 = (fit_2G12_rep1 + fit_2G12_rep2)/2
    out = "\t".join(map(str,[mut, name, mut_num, mut_input_count,
                             mut_noAb_rep1_count, mut_noAb_rep2_count, mut_2G12_rep1_count, mut_2G12_rep2_count,
                             fit_noAb_rep1, fit_noAb_rep2, fit_2G12_rep1, fit_2G12_rep2, fit_noAb, fit_2G12]))
    outfile.write(out+"\n")
  outfile.close()
      
def main():
  outfile        = "result/Pan99_mut_fit.tsv"
  file_count     = "result/Pan99_mut_count.tsv"
  file_mut_dict  = 'doc/mut_list.tsv'
  WT             = 'NKKEIAVN'
  mut_dict, resi_list, mut_list = reading_mut_dict(file_mut_dict)
  count_dict  = reading_count_dict(file_count)
  count_to_fit(outfile, count_dict, mut_dict, mut_list, WT, resi_list)

if __name__ == "__main__":
  main()
