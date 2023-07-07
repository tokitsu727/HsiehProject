from pysam import VariantFile
from collection import Counter
import vcf, sys

af_filter_level = .1

def hard_filter(ifile, ofile):
    filter_str = "MIN(FMT/GQ) >= 99 && MIN(FMT/DP) > 10"
    pysam.bcftools.filter("-O", "v", "-o", ofile, "-i", filter_str, ifile)

def double_genotype(f_record):
    #TODO: Check this one liner
    return Counter(getattr(sample, 'GT') for sample in f_record.samples).values().count(2) == 0

def af_filter(f_record):
    for sample in f_record.samples:
        if sample['DP'] and sample['AD'][1]/sample['DP'] < af_filter_level:
            return 0

def double_filter(f_record):
    return double_genotype and af_filter

def filter_str_builder(gq, dp):
    #TODO: Make this more expandable
    filter_str = "MIN(FMT/GQ) >= {0} && MIN(FMT/DP) > {1}".format(gq, dp)
    return filter_str

def custom_filters(ifile, ofile):
    #We use pyvcf to get an iterable for the file
    #    and then use itertools filter to filter them
    #Albeit, perhaps this should've just been using the pyvcf filter tool anyways.
    #TODO: Test performance of this
    vcf_reader = vcf.Reader(filename=ifile)
    #Is an iterator
    filtered_vcf = filter(double_filter, vcf_reader)
    reader_template = vcf.Reader(filename='format.vcf')
    vcf_writer = vcf.Writer(open(ofile, 'w'), reader_template)
    i = 0
    for record in filtered_vcf:
        vcf_writer.write_record(record)
        if i == 50:
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
    vcf_reader_2=vcf.Reader(filename=input_vc2)

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
        percent_1 = matching_num/len(pos_mem_1)
        percent_2 = matching_num/len(pos_mem_2)
    }

    f.close()
    return matching_data

if __name__ == '__main__':
    ifile = sys.argv[1] # Input file in vcf zipped
    ofile = sys.argv[2] # Output file in vcf unzipped format
    matchfile = sys.argv[3] # File of called data to match 

    hard_filter(ifile, ofile, filter_str_builder(99, 10))
    custom_filters(ofile, ofile)
    matching_data = match_positions(ofile, matchfile)










