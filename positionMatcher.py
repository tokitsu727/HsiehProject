import vcf, sys

input_vcf1 = sys.argv[1]
vcf_reader_1=vcf.Reader(filename=input_vcf1)

input_vcf2 = sys.argv[2]
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
        f.write(str(item) + '\n')
        matching_num += 1

f.close()

print("Match \% in 1: " + str(matching_num/len(pos_mem_1)))
print("Match \% in 2: " + str(matching_num/len(pos_mem_2)))

