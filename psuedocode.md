Skip through lines in Samfile starting with "@"

create fwd_strand dictionary
create rev_strand dictionary

**def soft_clip_adjust (samfile_line):**
> *'''read samfile line and adjusts start position if CIGAR string indicates softclipping'''* \
  (the changes are in bold) \
  **input:**  NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	**76814284**	36	**2S71M**	*	0	0	<sequence>	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU \
  **output:** NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	2	**7681428288**	36	**71M**	*	0	0	<sequence>	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU 

**def
  
  
Main_loop:
  `Reads samfile line by line, and for each line it does the following:`
