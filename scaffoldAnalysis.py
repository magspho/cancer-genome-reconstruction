import os
import re
import sys
import csv
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict
from glob import glob
import matplotlib as mpl
params = {'axes.spines.top': False, 'axes.spines.right': False}
mpl.rcParams.update(params)

col2idx = {'ctgName':0,'ctgLen':1,'ctgStart':2,'ctgEnd':3,'ctgStrand':4,
           'chrName':5,'chrLen':6,'chrStart':7,'chrEnd':8,
           'ctgbLen':9,'chrbLen':10,'ctgbRatio':11,'chrbRatio':12}

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def concatenate_files_helper(filenames, out_f):
    for filename in filenames:
        chrName = filename.split('/')[-1].split('.')[0]
        chrInfo = filename.split('/')[-1].split('.')[1]
        out_f.write(f">{chrName}\t{chrInfo}\n")
        with open(filename,'r') as in_f:
            for line in in_f.readlines():
                if line.startswith('CC'):
                    continue
                out_f.write(line)

def concatenate_files(files_dir, primary=True, sample_name = ''):
    hap1_filenames = []
    hap2_filenames = []
    file_tag = "primary." if primary else ""
    hap1_concat_name = f"{sample_name}hap1.{file_tag}seqid.tsv"
    hap2_concat_name = f"{sample_name}hap2.{file_tag}seqid.tsv"
    for name in glob(f'{files_dir}/chr*.hap1.{file_tag}seqid.tsv'):
        hap1_filenames.append(name)
    for name in glob(f'{files_dir}/chr*.hap2.{file_tag}seqid.tsv'):
        hap2_filenames.append(name)
    hap1_filenames.sort(key=natural_keys)
    hap2_filenames.sort(key=natural_keys)
    hap1_concat_file = open(hap1_concat_name,'w')
    hap2_concat_file = open(hap2_concat_name,'w')
    concatenate_files_helper(hap1_filenames,hap1_concat_file)
    hap1_concat_file.close()
    concatenate_files_helper(hap2_filenames,hap2_concat_file)
    hap2_concat_file.close()
    
def hap_whole_analysis(hap_concat_file):
    hap = hap_concat_file.split('/')[-1].split('.')[0]
    curr_chr = None
    out_df   = []
    ctgbRatio_dict = defaultdict(list)
    ctgcnt_dict    = defaultdict(int)
    ctgCOcnt_dict  = defaultdict(int)
    with open(hap_concat_file,'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                curr_chr = line.lstrip('>').split('\t')[0]
                continue
            tmp = line.rstrip().split('\t')
            ctgcnt_dict[curr_chr] += 1
            lst = [atoi(x) for x in tmp[1:]]
            lst = [float(x) if isinstance(x,str) and '.' in x else x for x in lst]
            if tmp[0] == 'PR':
                out_df.append(lst)
                ctgbRatio_dict[curr_chr].append(float(lst[col2idx['ctgbRatio']]))
            if tmp[0] == 'CO':
                ctgCOcnt_dict[curr_chr] += 1
#     print(out_df)
    out_df = pd.DataFrame(out_df, columns=list(col2idx.keys()))
    return out_df,ctgbRatio_dict,ctgcnt_dict,ctgCOcnt_dict

def plot_violin(out_df):
    # violin plot all
    fig, ax = plt.subplots(figsize=(8, 15))
    sns.violinplot(out_df,x="ctgLen",y="chrName",hue="hap", 
                   cut=0, split=True, gap=.1, inner="quart", 
                   density_norm='width', linewidth=0.5, saturation=0.7)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels,loc='center right')
    # plt.yscale('log')

def plot_violin_byChr(out_df):
    # violin plot by chr
    fig, axes = plt.subplots(nrows=4,sharey=True,figsize=(10, 12))
    tmp1 = out_df.loc[out_df["chrName"].isin(['chr'+str(i) for i in range(1,7)])]
    g1 = sns.violinplot(tmp1,x="chrName",y="ctgLen",hue="hap", 
                        cut=0, split=True, gap=.1, inner="quart", 
                        linewidth=0.5, legend=False, ax=axes[0])
    g1.set(xlabel=None)
    tmp2 = out_df.loc[out_df["chrName"].isin(['chr'+str(i) for i in range(7,13)])]
    g2 = sns.violinplot(tmp2,x="chrName",y="ctgLen",hue="hap", 
                        cut=0, split=True, gap=.1, inner="quart", 
                        linewidth=0.5, legend=False, ax=axes[1])
    g2.set(xlabel=None)
    tmp3 = out_df.loc[out_df["chrName"].isin(['chr'+str(i) for i in range(13,19)])]
    g3 = sns.violinplot(tmp3,x="chrName",y="ctgLen",hue="hap", 
                        cut=0, split=True, gap=.1, inner="quart", 
                        linewidth=0.5, legend=False, ax=axes[2])
    g3.set(xlabel=None)
    tmp4 = out_df.loc[out_df["chrName"].isin(['chr'+str(i) for i in [19,20,21,22,'X']])]
    g4 = sns.violinplot(tmp4,x="chrName",y="ctgLen",hue="hap", 
                        cut=0, split=True, gap=.1, inner="quart", 
                        linewidth=0.5, legend=False, ax=axes[3])
    g4.set(xlabel=None)
    # plt.yscale('log')
    axes[2].legend([Rectangle((0,0),3,1.5, facecolor='#1f77b4',edgecolor='black',linewidth=0.5),
                     Rectangle((0,0),3,1.5, facecolor='#ff7f0e',edgecolor='black',linewidth=0.5)],
                     ['hap1', 'hap2'])

def plot_box(out_df, log=False):
    # boxplot not log-scale
    fig, ax = plt.subplots(figsize=(8, 15))
    sns.boxplot(out_df,x="ctgLen",y="chrName",hue="hap",linewidth=0.5,
                gap=.2,width=.6,saturation=0.7,flierprops={"marker": "."},)
    [plt.axhline(y, color = 'grey', linestyle='-',linewidth=0.5) for y in np.array(range(len(np.unique(out_df['chrName']))))+0.5]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels,loc='center right')
    if log:
        plt.xscale('log')
        
def plot_ctgBRatio():
    pass