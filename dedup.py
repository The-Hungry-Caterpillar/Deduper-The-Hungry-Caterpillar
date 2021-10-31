def get_args(): #defines all the independent variables
	import argparse
	parser = argparse.ArgumentParser(description = "still need description")
	#parser.add_argument('- command line variable', '-- python variable', description)
	parser.add_argument('-f', '--file', help='sam filename to be deduped')
	
	return parser.parse_args()

args=get_args()

f=open(args.file, 'r')

line=f.readline()
f.close()
print(line)
