import re

test_mode=True
def get_args(): #defines all the independent variables
	import argparse
	parser = argparse.ArgumentParser(description = "still need description")
	#parser.add_argument('- command line variable', '-- python variable', description)
	parser.add_argument('-f', '--file', help='sam filename to be deduped')
	
	return parser.parse_args()

args=get_args()

def UMI(line):
    '''finds the UMI in the zeroeth column of SAM line'''
    UMI=line[0].split(':')[7]
    return(UMI)

def chrom(line):
    '''finds the chromosome number in second column of SAM line'''
    chrom=line[2]
    return(chrom)

def ref_position(line):
    '''finds the leftmost reference-based position of read in SAM line'''
    ref_position=line[3]
    return(ref_position)

def is_rev_strand(line):
    flag=int(line[1])
    if ((flag & 16) == 16):
        return True

def adjusted_position(line):
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


f=open(args.file, 'r')
line=f.readline().strip().split('\t')

adjusted_position(line)
# UMI=UMI(line)
# chrom=chrom(line)
# ref=ref_position(line)

# cigar=line[5]
# cigar_num=re.split("[A-Z]",cigar)
# cigar_letters=(re.split("[0-9]+",cigar))
# del cigar_letters[0]

# if (is_rev_strand(line)):
#     print('reverse strand')
f.close()

# print('UMI is {}'.format(UMI))
# print('chrom is {}'.format(chrom))
# print('ref position is {}'.format(ref))
# print('cigar numbers are {}'.format(cigar_num))
# print('cigar letters are {}'.format(cigar_letters))