import vcf, sys

#vcf_reader = vcf.Reader(open('KUNUSCCLH_0001_01_LibTube81_R00250S8A2M0000P0000_C1_KHWGSH_A00006-KUNUSCCLH_0001_01_LibTube11_R00250S8A2M0000P0000_T2_KHWGSH_A00002.HC_All.vcf', 'r'))

input_vcf = sys.argv[1]
vcf_reader = vcf.Reader(filename=input_vcf)

reader_template = vcf.Reader(filename=input_vcf)
vcf_writer = vcf.Writer(open(sys.argv[2], 'w'), reader_template)

i = 0

af_filter_level = .1

for record in vcf_reader:
    writev = 1 
    for sample in record.samples:
        #Filters if any sample is below filter level
        #Potentially incorrect, not last number. Need to ask.

        #ADLen = len(sample['AD'])
        if sample['DP'] and sample['AD'][1]/sample['DP'] < af_filter_level:
            writev=0
                

    if writev:
        vcf_writer.write_record(record)
        if i == 10:
            vcf_writer.flush()
            i = 0
        i += 1

vcf_writer.close()
