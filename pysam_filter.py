from pysam import VariantFile
import pysam.bcftools
from collections import Counter
from statistics import median
import vcf, sys, pysam
import argparse
import matplotlib.pyplot as plt 
import random
import numpy as np
import pandas as pd
import os.path

af_filter_level = .1

#TODO: Can I make this more pythony?

def hard_filter(ifile, ofile, filter_str):
    pysam.bcftools.filter("-O", "v", "-o", ofile, "-i", filter_str, ifile, catch_stdout=False)

def double_genotype(f_record):
    #TODO: Check this one liner
    #Probably not worth it, this is pretty fast anyways.
    #return Counter(for sample in f_record.samples).values().count(2) == 0
    gt = []
    for sample in f_record.samples:
        if sample['GT'] in gt:
            return 0
        if sample['GT'] == './.':
            return 0
        gt.append(sample['GT'])
    return 1

def af_filter(f_record):
    for sample in f_record.samples:
        if sample['DP'] and sample['AD'][1]/sample['DP'] < af_filter_level:
            return 0
    return 1

def double_filter(f_record):
    return double_genotype(f_record) and af_filter(f_record)

def filter_str_builder(gq, dp):
    #TODO: Make this more expandable
    filter_str = "MIN(FMT/GQ) >= {0} && MIN(FMT/DP) > {1}".format(gq, dp)
    print(filter_str)
    return filter_str

def custom_filters(ifile, ofile):
    #We use pyvcf to get an iterable for the file
    #    and then use itertools filter to filter them
    #Albeit, perhaps this should've just been using the pyvcf filter tool anyways.
    #TODO: Test performance of this
    #       Addendum: Performance is very good, it's just the io that sucks
    #       But I'm bad at IO unfortunately
    vcf_reader = vcf.Reader(filename=ifile)
    #Is an iterator
    filtered_vcf = filter(double_filter, vcf_reader)
    reader_template = vcf.Reader(filename='out.vcf')
    vcf_writer = vcf.Writer(open(ofile, 'w'), reader_template)
    i = 0
    for record in filtered_vcf:
        #This is the bottleneck
        vcf_writer.write_record(record)
        if i == 500:
            vcf_writer.flush()
            i=0
        i += 1
    vcf_writer.close()

def match_positions(input_vcf1, input_vcf2, outfile, write_out):
    #TODO: Find better way to do this so that we can retreive record information
    #   Probably relatively difficult to do, depending on complexity of search vcf file
    #   This current method is O(n) based on the size of input_vcf_1
    #   Best i can think of for alternate is like O(n^2) but with io and memory issues

    vcf_reader_1=vcf.Reader(filename=input_vcf1)
    vcf_reader_2=vcf.Reader(filename=input_vcf2)


    print("Begin matching.")
    print("Reading File 2")

    pos_mem_2 = []

    for record in vcf_reader_2:
        if (record.FILTER == [] or record.FILTER == ["HighDepth"]) and record.POS not in pos_mem_2:
            pos_mem_2.append([record.POS,record.CHROM])
    print(str(len(pos_mem_2)))

    reader_template = vcf.Reader(filename=input_vcf1)
    matching_file = 'merged_files_strelka/matching/' + outfile
    unmatching_file = 'merged_files_strelka/unmatching/' + outfile
    if write_out:
        vcf_writer = vcf.Writer(open(matching_file, 'w'), reader_template)
        vcf_writer2 = vcf.Writer(open(unmatching_file, 'w'), reader_template)
    pos_mem_2 = set(pos_mem_2)
    matching_num = 0
    #TODO: Better solution for these flushes?
    i = 0
    j = 0

    print("Reading File 1")

    pos_mem_1 = 0

    #TODO: Clean this up so we're not leaving a bunch of unused writers around

    for record in vcf_reader_1:
        pos_mem_1 += 1
        if (record.POS, record.CHROM) in pos_mem_2 and (record.FILTER == [] or record.FILTER ==["HighDepth"]):
            if write_out:
                vcf_writer.write_record(record)
                if i == 500:
                    vcf_writer.flush()
                    i = 0
                i += 1
            matching_num+=1
        else:
            if write_out:
                vcf_writer2.write_record(record)
                if j == 500:
                    vcf_writer2.flush()
                    j = 0
                j += 1
    if write_out:
        vcf_writer.close()
        vcf_writer2.close()
    #TODO This isn't super helpful, at least for the second file, probably need to be checking
    #        stuff based on relative sizes of the files to contextualize the matching %.
    #        AKA normalize       
    matching_data = {
        "percent_1":  matching_num/pos_mem_1,
        "percent_2":  matching_num/len(pos_mem_2)
    }

    return matching_data
def default_mode(input_file, output_file, matching_file):
    matching_data = match_positions(input_file, matching_file, output_file, True)
    print("Match %1: {0}\nMatch %2: {1}".format(matching_data['percent_1'], matching_data['percent_2']))

def compute_bins(data, bin_size):
    min_val = np.min(data)
    max_val = np.max(data)
    min_bound = -1.0 * (min_val % bin_size - min_val)
    max_bound = max_val - max_val % bin_size + bin_size
    n_bins = int((max_bound - min_bound) / bin_size) + 1
    bins = np.linspace(min_bound, max_bound, n_bins)
    return bins

def calculate_af(t_alt, t_depth, reflect_graph):
    if not reflect_graph:
        return t_alt/t_depth
    AF = t_alt/t_depth
    if AF == 0 or AF == 1:
        return 1
    if AF > .5:
        return 1 - AF
    return AF

def create_list(input_file, reflect_graph, chr_n=False, chr_numbered=False):
    af = []

    df = pd.read_csv(input_file, header=1, sep='\t')
    if chr_numbered:
        mask = df['Chromosome'].str.contains('^chr\d+$', regex=True)
        df = df.loc[mask]
        
    if chr_n: 
        chrx = df.loc[df['Chromosome']==chr_n] #modify this line to change which chromosome is plotted
    else:
        chrx = df
    try:
        af = chrx.apply(lambda x: calculate_af(x['t_alt_count'], x['t_depth'], reflect_graph), axis=1).tolist()
        return af
    except:
        print("No values found at" + chr_n + ", returning 0")
        return [0]


def graph_mode(input_file, output_file, chr_n):
    AF = create_list(input_file)# + create_list(input_file2)

    bins=compute_bins(AF, 0.01)

    plt.xlabel("AF")
    plt.figure(figsize=(20,6))
    plt.title(output_file)
    plt.hist(AF, bins=bins, edgecolor="white", zorder=2)
    plt.xticks(bins, rotation=70)
    plt.savefig("graphs/"+output_file+".pdf")


def bulk_graph_mode(input_files, output_file, reflect_graph, chr_n):

    bins = [x * .01 for x in range(0,101)]
    #Next two lines for non dragen
    format2 = False
    if input_files[0].find("KUNUSCCLH") != -1:
        format2 = True
    if format2: 
        dict_input = {x:x[x.find("_T")+2] for x in input_files}

        input_files = sorted(dict_input.keys(), key=dict_input.get)
    
    fig, axs = plt.subplots(len(input_files),figsize=(20,18), sharex=True, sharey=True)
    fig.suptitle(output_file, size='xx-large', fontweight='heavy')
    number = 0
    for i in input_files:
        AF = create_list(i, reflect_graph, chr_n)
        axs[number].hist(AF, bins=bins, edgecolor="white", zorder=2)
        j = i.find("_T")
        #axs[number].set_title(i[2:4])
        if format2:
            j = i.find("_T")
            axs[number].set_title(i[j+1:j+3] + " " + chr_n)
        else:
            j = i.find("PASS")
            axs[number].set_title(i[j-14:j-12] + " " + chr_n)

        axs[number].set_xlabel("AF")
        axs[number].set_xticks(bins)
        axs[number].tick_params(labelrotation=70)
        
        number += 1

    plt.subplots_adjust(hspace=.2)
    plt.savefig("graphs/"+output_file+'.'+chr_n+".pdf")

def one_graph_mode(input_files, output_file, reflect_graph):

    bins = [x * .01 for x in range(0,101)]
    #Next two lines for non dragen
    format2 = False
    if input_files[0].find("KUNUSCCLH") != -1:
        format2 = True
    if format2: 
        dict_input = {x:x[x.find("_T")+2] for x in input_files}

        input_files = sorted(dict_input.keys(), key=dict_input.get)
    
    fig, axs = plt.subplots(len(input_files),figsize=(20,18), sharex=True, sharey=True)
    fig.suptitle(output_file, size='xx-large', fontweight='heavy')
    number = 0
    for i in input_files:
        AF = create_list(i, reflect_graph, 'chr'+str(1))
        for chr_n in range(2,23):
            AF.extend(create_list(i, reflect_graph, 'chr'+str(chr_n)))

            
        axs[number].hist(AF, bins=bins, edgecolor="white", zorder=2)
        j = i.find("_T")
        #axs[number].set_title(i[2:4])
        if format2:
            j = i.find("_T")
            axs[number].set_title(i[j+1:j+3])
            t = i[j+1:j+3]
            s_i = i.find("_00")
            s = i[s_i+1:s_i+5]
        else:
            j = i.find("PASS")
            axs[number].set_title(i[j-14:j-12])
            t = i[j-14:j-12]
            s = i[j-16:j-14]

        print(s+t + "\t" + str(median(create_list(i, reflect_graph, chr_numbered=True)))) 
        axs[number].set_xlabel("AF")
        axs[number].set_xticks(bins)
        axs[number].tick_params(labelrotation=70)
        
        number += 1

    plt.subplots_adjust(hspace=.2)
    plt.savefig("graphs/"+output_file+".pdf")

def create_by_pos(input_file, chr_n):
    data = {}
    df = pd.read_csv(input_file, header=1, sep='\t')
    chrx = df.loc[df['Chromosome']==chr_n]
    try:
        data['AF'] = chrx.apply(lambda x: calculate_af(x['t_alt_count'], x['t_depth'], False), axis=1).tolist()
        data['pos'] = chrx['Start_Position']
        return data
    except:
        return {}
        

def AF_by_pos(input_file, output_file):
    chrs = list(range(1,23)) + ['X','Y']
    data = {}
    fig, axs =plt.subplots(24, figsize=(20,200), sharey=True)
    fig.suptitle(output_file, size='xx-large', fontweight='heavy')
    number = 0
    for i in chrs:
        AF = create_list(input_file, False, 'chr'+str(i))
        data = create_by_pos(input_file, 'chr'+str(i))
        if data:
            axs[number].scatter(data['pos'], data['AF'])
            axs[number].set_title('chr'+str(i))
            axs[number].set_xlabel("Start Position")
            axs[number].set_ylabel("AF")
            number += 1
    plt.savefig("graphs/"+output_file+".pdf")


    

def merge_csvs(input_file1, input_file2, output_file):
    on_list = ['Chromosome', 'Start_Position']
    df1 = pd.read_csv(input_file1, header=1, sep='\t')
    df2 = pd.read_csv(input_file2, header=1, sep='\t')

    df1_len = len(df1)
    df2_len = len(df2)

    df1['Start_Position'] = df1.Start_Position.astype(str, errors='raise')
    df2['Start_Position'] = df2.Start_Position.astype(str, errors='raise')

    merged_df = pd.merge(df1, df2, how='inner', on=on_list, suffixes = ["_Dragen","_MuTect"])
    merged_df.drop_duplicates(subset=on_list, keep='first')

    combined_len = len(merged_df) 

    merged_df.to_csv(output_file, sep='\t', index=False)

    print(f'{output_file},{df1_len},{df2_len},{combined_len}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_file", help="Input File")
    #parser.add_argument("-j", "--input_file_graph", help="Temporary Input File for graphing")
    parser.add_argument("-o", "--output_file", help="Intermediary output file") #Set this flag only if you need to do the same type haplo filter again
    parser.add_argument("-m", "--matching_file", help="File to compare to")
    parser.add_argument("-g", "--make_graph", action="store_true", help="Run with a variety of gq and dp values.")
    parser.add_argument("-c", "--merge_csv", action="store_true", help="Run with csvs instead of vcf")
    parser.add_argument("--bulk_graph", nargs="*",type=str)
    parser.add_argument("--reflect_graph", action="store_true")
    parser.add_argument("--one_graph", action="store_true")
    parser.add_argument("--by_pos", action="store_true")
    parser.add_argument("--output_medians", action="store_true")

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    matching_file = args.matching_file
    merge_csv = args.merge_csv
    bulk_graph = args.bulk_graph
    reflect_graph = args.reflect_graph
    one_graph = args.one_graph
    by_pos = args.by_pos
    output_medians = args.output_medians
    #input_file_graph = args.input_file_graph

    #if True:
    #    graph_mode(input_file, output_file)
    #    exit()
    if merge_csv:
        merge_csvs(input_file, matching_file, output_file)
        exit()
    if one_graph:
        one_graph_mode(bulk_graph, output_file, reflect_graph) 
        exit()
    if by_pos:
        AF_by_pos(input_file, output_file)
        exit()
    if bulk_graph:
        for i in range(1, 27):
            bulk_graph_mode(bulk_graph, output_file, reflect_graph, 'chr'+str(i))
        exit()

            
    
    default_mode(input_file, output_file, matching_file)
