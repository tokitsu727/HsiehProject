from pysam import VariantFile
import pysam.bcftools
from collections import Counter
import vcf, sys, pysam

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
    return double_genotype(f_record)# and af_filter(f_record)

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

def match_positions(input_vcf1, input_vcf2):
    #TODO: Find better way to do this so that we can retreive record information
    #   Probably relatively difficult to do, depending on complexity of search vcf file
    #   This current method is O(n) based on the size of input_vcf_1
    #   Best i can think of for alternate is like O(n^2) but with io and memory issues

    vcf_reader_1=vcf.Reader(filename=input_vcf1)
    vcf_reader_2=vcf.Reader(filename=input_vcf2)

    pos_mem_1 = []
    pos_mem_2 = []

    for record in vcf_reader_1:
        pos_mem_1.append(record.POS)
    for record in vcf_reader_2:
        pos_mem_2.append(record.POS)

    pos_mem_2 = set(pos_mem_2)
    f = open('matching_file', 'w')
    matching_num = 0

    for item in pos_mem_1:
        if item in pos_mem_2:
            f.write(str(item)+'\n')
            matching_num+=1

    matching_data = {
        "percent_1":  matching_num/len(pos_mem_1),
        "percent_2":  matching_num/len(pos_mem_2)
    }

    f.close()
    return matching_data

if __name__ == '__main__':
    #This is full pipeline. Takes about 5-10 minutes
    ifile = sys.argv[1] # Input file in vcf zipped
    ofile = sys.argv[2] # Output file in vcf unzipped format
    matchfile = sys.argv[3] # File of called data to match 

    hard_filter(ifile, ofile, filter_str_builder(90, 10))
    custom_filters(ofile, 'out1.vcf')
    matching_data = match_positions('out1.vcf', matchfile)
    print("Match %1: {0}\nMatch %2: {1}".format(matching_data['percent_1'], matching_data['percent_2']))
