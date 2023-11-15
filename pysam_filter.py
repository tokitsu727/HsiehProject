from pysam import VariantFile
import pysam.bcftools
from collections import Counter
import vcf, sys, pysam
import argparse
import matplotlib.pyplot as plt 

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
            pos_mem_2.append(record.POS)
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
        if record.POS in pos_mem_2 and (record.FILTER == [] or record.FILTER ==["HighDepth"]):
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

def graph_mode(input_file, output_file, matching_file):
    match1 = []
    match2 = []
    for gq in range(0, 99):
        hard_filter(input_file, "out1.vcf", filter_str_builder(gq, 15))
        matching_data = match_positions("out1.vcf", matching_file, False)
        match1.append(matching_data['percent_1'])
        match2.append(matching_data['percent_2'])
    plt.figure(0)
    plt.xlabel("GQ Threshold")
    plt.title("Matching % of two files by GQ")
    plt.plot(match1, 'g', label = 'Haplotype')
    plt.plot(match2, 'r', label = 'MuTect')
    plt.savefig("gqaf.pdf")

    match1 = []
    match2 = []

    for dp in range(0, 99):
        hard_filter(input_file, "out1.vcf", filter_str_builder(99, dp))
        matching_data = match_positions("out1.vcf", matching_file, False)
        match1.append(matching_data['percent_1'])
        match2.append(matching_data['percent_2'])
    plt.figure(1)
    plt.xlabel("DP Threshold")
    plt.title("Matching % of two files by DP")
    plt.plot(match1, 'g', label = 'Haplotype')
    plt.plot(match2, 'r', label = 'MuTect')
    plt.savefig("dpaf.pdf")
       

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_file", help="Input File")
    parser.add_argument("-o", "--output_file", help="Intermediary output file") #Set this flag only if you need to do the same type haplo filter again
    parser.add_argument("-m", "--matching_file", help="File to compare to")
    parser.add_argument("-g", "--make_graph", action="store_true", help="Run with a variety of gq and dp values.")

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    matching_file = args.matching_file
    print(output_file)
    print(matching_file)
    print(input_file)
    
    default_mode(input_file, output_file, matching_file)
