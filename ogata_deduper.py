import re

test_mode=False
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description = "still need description")
	#parser.add_argument('- command line variable', '-- python variable', description)
    parser.add_argument('-f', '--file', help='sam filename to be deduped')
    parser.add_argument('-u', '--umi', help='name of umi file')
    return parser.parse_args()
args=get_args()

def UMI_builder(umi_file,umi_dict): #input UMI file and empty UMI dictionary: no output, appends to empty dictionary
    '''takes umi file and empty umi dictionary as input, note that the umi file MUST contain only UMI sequences separated by new lines'''
    u=open(umi_file, 'r')
    while True:
        umi=u.readline().strip()
        if umi == '':
            break
        umi_dict[umi]=''
    u.close()

def UMI(line): #input sam line: output UMI of read (assumes UMI is last value in QNAME)
    '''finds the UMI in the zeroeth column of SAM line'''
    UMI=line[0].split(':')[-1]
    return(UMI)

def chrom(line): #input sam line: output chromosome of read
    '''finds the chromosome number in second column of SAM line'''
    chrom=line[2]
    return(chrom)

def ref_position(line): #input sam line: output left-most start position
    '''finds the leftmost reference-based position of read in SAM line'''
    ref_position=line[3]
    return(ref_position)

def is_rev_strand(line): #input sam line: output True if rev strand else False
    flag=int(line[1])
    return True if ((flag & 16) == 16) else False

def adjusted_position(line): #input sam line: output adjusted start position -- strand specific
    '''uses cigar string to adjust the starting position
    if the strand is reverse then we use cigar string to find right-most adjusted position
    if the strand is forward then we use cigar string to find left-most adjusted position'''

    cigar=line[5]
    cigar_num=re.split("[A-Z]",cigar)
    cigar_letters=(re.split("[0-9]+",cigar))
    del cigar_letters[0] #empty value, remove
    del cigar_num[-1] #empty value, remove

    if is_rev_strand(line):
        '''if reverse strand then adjusted_position=ref_position + Ms + right Ss + Ds + Ns'''
        if cigar_letters[-1]=='S':
            S_sum=int(cigar_num[-1])
        else:
            S_sum=0

        M_indices=[i for i, x in enumerate(cigar_letters) if x == 'M']
        N_indices=[i for i, x in enumerate(cigar_letters) if x == 'N']
        D_indices=[i for i, x in enumerate(cigar_letters) if x == 'D']
       
        M_sum=0
        for i in M_indices:
            M_sum+=int(cigar_num[i])
        if test_mode:
            print('Msum is {}'.format(M_sum))
       
        N_sum=0
        for i in N_indices:
            N_sum+=int(cigar_num[i])
        if test_mode:
            print('Nsum is {}'.format(N_sum))

        if test_mode:
            print('Ssum is {}'.format(S_sum))

        D_sum=0
        for i in D_indices:
            D_sum+=int(cigar_num[i])
        if test_mode: 
            print('Dsum is {}'.format(D_sum))
        
        adjusted_position=int(ref_position(line)) + M_sum + S_sum + D_sum + N_sum
    
    else:
        '''if forward strand then adjusted_position = ref_position - left Ss'''
        if cigar_letters[0]=='S':
            adjusted_position=int(ref_position(line))-int(cigar_num[0])
        else:
            adjusted_position=int(ref_position(line))
    
    if test_mode:
        print('adjusted position is {}'.format(adjusted_position))
    
    return adjusted_position





out=(args.file + '_deduped')
print(out)
f=open(args.file, 'r')
w=open(out, 'w')

umi_dict={}
if args.umi != 'NULL':
    UMI_builder(args.umi, umi_dict)

bad_UMI_count=0
duplicate_count=0
total_count=0
unique=0

current_dict={}

prev_chrom=''

while True:
    '''reads through samfile and removes PCR duplicates'''
    line=f.readline().strip().split('\t')

    if line == ['']:
        print("The number of unique reads in chromosome {} is {}.".format(prev_chrom,unique))
        break
    
    total_count+=1
    
    if args.umi != 'NULL':
        if UMI(line) not in umi_dict: #make sure that the UMI is valid
            bad_UMI_count+=1
            prev_chrom=current_chrom
            continue

    current_chrom=chrom(line)
    
    umi=UMI(line)
    position=adjusted_position(line)
    strand=('rev' if is_rev_strand(line) else 'fwd')
        
    if current_chrom != prev_chrom: #since samfile is sorted we can clear the dictionary once we move onto the next chromosome, for memories sake.
        print("The number of unique reads in chromosome {} is {}.".format(prev_chrom,unique))
        unique=0
        current_dict.clear()

    if position in current_dict:

        if (strand, umi) in current_dict[position]:
            # print('\t'.join(line))
            prev_chrom=current_chrom
            duplicate_count+=1
            continue
        else:
            # print('\t'.join(line))
            # current_dict[position]={}
            current_dict[position][strand,umi]=''
            unique+=1
            print('\t'.join(line),file=w)
        
    else:
        # print('\t'.join(line))
        current_dict[position]={}
        current_dict[position][strand,umi]=''
        unique+=1
        print('\t'.join(line),file=w)

    prev_chrom=current_chrom
    
w.close()
f.close()

percent_bad=round((100*(bad_UMI_count+duplicate_count)/total_count), 2)
print("\n\n\n")
print("Removed all reads with UMIs not in provided UMI file, and removed all PCR duplicate reads")
print("The total number of reads in the original sam file is {}.".format(total_count))
print("The total number of reads in the deduped sam file is {}.".format(total_count-bad_UMI_count-duplicate_count))
print("The total number of bad UMIs is {}.".format(bad_UMI_count))
print("The total number of PCR duplicates is {}.".format(duplicate_count))
print("The percentage of removed reads is {}%.".format(percent_bad))
