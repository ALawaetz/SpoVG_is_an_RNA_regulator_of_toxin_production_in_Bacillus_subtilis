import os
import sys
import glob

input_folder = sys.argv[1]
output_folder = sys.argv[2]

novoindex = '/Users/andreas/Bacillus/Bioinformatics/SpoVG_CLASH_LB_Andreas/B_subtilis_NC_000964.3.novoindex'

### script to align batch of CRAC fastq files 



file1 = []
file2 = []
for file in sorted(glob.glob(input_folder + '/*.fastq')):
    if os.path.basename(file).find('trimmed_1') == -1:
        file2.append(file)
    else:
        file1.append(file)

for f1, f2 in zip(file1, file2):
    print(f'aligning {os.path.basename(f1)} and {os.path.basename(f2)}')
    os.system(f"novoalign -d {novoindex} -f {f1} {f2} -r Random > {output_folder}/{os.path.basename(f1)[:-6]}.novo") 
