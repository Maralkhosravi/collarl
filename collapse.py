import numpy as np

def next_haplo_name(hh):
    if hh < 10:
        name = 'H0' + str(hh)
    else: name = 'H' + str(hh)
    return name

fname = 'sequence_data_cropped.fas'
handle = open(fname)
haplos = {} #haplos have names A, B, C etc. as keys and sequence as value
haplo_names = [] #list of haplotype names e.g. H01
sample_names = [] #list of sample names e.g. DH

#fill the fasta outfile with haplotypes
#and fill the lists haplo_names and sample_names
out1 = open('sequence_data_haplos.fas','w')
h = 1 #number for the first haplotype
while True:
    line1 = handle.readline() #the name, e.g. DH01-0001
    if not line1: break
    line1 = line1.strip('\n')
    sample = line1[1:3]
    if sample not in sample_names: #check if sample seen before, if not, add it
        sample_names.append(sample)
    line2 = handle.readline() #the DNA sequence
    line2 = line2.upper()#make all uppercase
    if line2 not in haplos.values(): #fill the haplos dictionary
        haplo_name = next_haplo_name(h)
        haplo_names.append(haplo_name)
        h += 1
        haplos[haplo_name] = line2
        line_h = '>' + haplo_name +'\n' + line2
        out1.write(line_h)
handle.close()
out1.close()

#fill the first csv outfile with the individuals' haplotypes
#will use infile fname = 'sequence_data_cropped.fas' for this
out2 = open('sequence_data_ind_haplos.csv','w')
handle = open(fname)
while True:
    line1 = handle.readline() #the name, e.g. DH01-0001
    if not line1: break
    line1 = line1.strip('\n')
    line2 = handle.readline() #the DNA sequence
    for k,v in haplos.items(): #look for the haplotype code
        if v == line2:
            print(k,v)
            line1 = line1.strip('>').split('_')
            line1 = ','.join(line1)
            line_i = line1 +','+ k + '\n'
            out2.write(line_i)
handle.close()
out2.close()

#fill csv file (out3) with frequencies per sample
#will use lists haplo_names and sample_names for this
#will also use infile fname = 'sequence_data_ind_haplos.csv' for this
#produce Arlequin infile (out4)
handle2 = open('sequence_data_ind_haplos.csv')
out3 = open('sequence_data_freqs.csv','w')
out4 = open('sequence_data_samples.arp','w')
print('haplo_names: ',haplo_names)
print('sample_names: ',sample_names)
#write the head to the Arlequin file
line_arl = '[Profile]\nTitle="A sample file designed to compute amova"\nNbSamples=12\n'\
    'GenotypicData=0\nLocusSeparator=NONE\nDataType=DNA\nCompDistMatrix=0\n[Data]\n'\
    '[[Samples]]\n'
out4.write(line_arl)
freq_table = np.zeros(shape=(len(sample_names),len(haplo_names))) #array for frequency data
freq_table.astype(int)
#fill the first row of csv file out3
#haplo_names
line = ','
for haplo_name in haplo_names:
    line = line + haplo_name + ','
line = line + '\n'
out3.write(line)

#fill the freq_table array
while True:
    line1 = handle2.readline()
    if not line1: break
    line1 = line1.strip('\n')
    ind_sample = line1[0:2]
    ind_haplo = line1[-3:]
    ind_sample_index = sample_names.index(ind_sample)
    ind_haplo_index = haplo_names.index(ind_haplo)
    freq_table[ind_sample_index,ind_haplo_index] += 1
freq_table = freq_table.astype(int)
print(freq_table)
sample_sizes = np.sum(freq_table,axis=1).tolist() #list of sample sizes
#write freq_table array to csv file out3
#write the sample data to the Arlequin file
for i in range(len(sample_names)):
    line_i = sample_names[i] + ','
    line_arl = 'SampleName="' + sample_names[i] + '"\nSampleSize=' +\
        str(sample_sizes[i]) + '\nSampleData={\n'
    for j in range(len(haplo_names)):
        line_i = line_i + str(freq_table[i,j]) +','
        line_arl = line_arl + str(haplo_names[j]) + ' ' + str(freq_table[i,j]) + ' ' +\
            haplos[haplo_names[j]]
    line_i = line_i + '\n'
    line_arl = line_arl + '}\n'
    out3.write(line_i)
    out4.write(line_arl)
#write the footer (structure) to the Arlequin file
line_arl = '[[Structure]]\nStructureName="none"\nNbGroups=1\nIndividualLevel=0\n\
    Group={\n'
for k in sample_names:
    line_arl = line_arl + '"' + k + '"\n'
line_arl = line_arl + '}\n'
out4.write(line_arl)

handle2.close()
out3.close()
out4.close()