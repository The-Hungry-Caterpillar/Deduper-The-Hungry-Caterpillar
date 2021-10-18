## Python

**def soft_clip_adjust (samfile_line):**
> *'''read samfile line and adjusts start position if CIGAR string indicates softclipping'''* \
  (the changes are in bold) \
  *note that this function does not change the CIGAR string \
  **input:**  NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	**100**	36	2S71M	*	0	0	<sequence>	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU \
  **output:** NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	**98**	36	2S71M	*	0	0	<sequence>	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU \
  **OUT: soft_corrected_samfile_line**

  
  
**def modify_sam_file (original_samfile):**
> *'''modifies the start position of everyline in a samfile'''* \
  *read through original samfile and for each line, do **soft_clip_adjust** (ignores lines starting with "@") \
  then determine if the line is reverse compliment \
  **input:** original samfile \
  **output:** two output files: \
  (1) samfile with adjusted start positions all of which are forward strand \
  (2) samfile with adjusted start positions all of which are reverse compliment
  
modify_sam_file(original_samfile)
  
## Bash
  
Then, exit python and go to bash: \
  ```sort in-place modified_samfile by chromosome then position then UMI``` 
  
  ```sort in-place modified_samfile_rev_comp by chromosome then position then UMI```

## Python
then back into python:
  
**def deduper (modified_samfile, modified_samfile_revcomp):**
> *''' removes PCR duplicates '''* \
  read first line of modified_samfile and modified_samfile_rev_comp \
  **line** = samfile line \
  **line_revcomp** = revcomp samfile line \
  \
  **chromosome** = column 3 of samfile line \
  **chromosome_revcomp** = column 3 of revcomp samfile \
  \
  **position** = column 4 of samfile line \
  **position_revcomp** = column 4 of revcomp line \
  \
  **umi** = last part of column 1 \
  **umi_revcomp** = last part of revcomp column 1 \
  \
  then look at next line in samfile and revcomp samfile: \
  **if** chromosomes, position, and umi match, then read the next line. (do this for revcomp as well) \
  **if** chromosomes, position, and umi **do not** match, then print **line** to output file. (do this for revcomp as well) \
  \
  **input:** modified_samfile, modifed_samfile_revcomp \
  **output:** modified_deduped_samfile

deduper (modified_samfile, modified_samfile_revcomp) 

**def fix_soft_adjust (modified_deduped_samfile):** 
> *''' readjusts samfile lines to match *original unchange cigar strings* '''* \
  for each line in samfile: read CIGAR string and adjust position back to *original* position. \
  **input:** modified_deduped_samfile \
  **output:** deduped_samfile
  
fix_soft_adjust(modified_deduped_samfile)\
\
\
\
\
\
\
\
## Edit: I'm leaning more toward the method of:
  
### for each line:
  a. read and hold the line as a variable. check if feasible UMI, if not then discard line.
  
  b. fix the soft clipping 
  
  c. Look up corrected position in position_dictionary (position will be keys of a dictionary)\
     - if not in position dictionary: add it to position dictionary, write line to out file, and restart loop at (1)\
     - if in position dictionary: continue to next step
  
  d. Look up strand (strand will be a nested dictionary as the value of each position key from the position dictionary)\
    - if not in strand dictionary: add it to strand dictionary, write line to out file, and restart loop at (1)\
    - if in strand dictionary: continue to next step
 
  e. Look up chromosome (chromosome will be a nested dictionary as the value of each strand key from the strand dictionary)\
    - if not in chromosome dictionary: add it to chromosome dictionary, write line to out file, and restart loop at (1)\
    - if in chromosome dictionary: continue to next step
  
  f. Look up UMI (UMI will be a nested dictionary as the value of each chromosom key from the chromosom dictionary)\
    - if not in UMI dictionary: add it to UMI dictionary, write line to out file, and restart loop at (1)\
    - if in UMI dictionary: dicard read -- it is a PCR duplicate

  
The nested dictionary scheme is like a tree, so you only need to look as deep as to remove abiguity, and since dictionaries are ***extremely*** fast I think this will work well. 
