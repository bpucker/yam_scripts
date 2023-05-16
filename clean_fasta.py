import glob, os

# --- end of imports --- #


def load_seqs( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	seqs = {}
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				seqs.update( { header: "".join( seq ) } )
				header = line[1:].strip()
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		seqs.update( { header: "".join( seq ) } )
	return seqs



input_folder = "./cleaned_peps/"
output_folder = "./orthofinder_input/"

if not os.path.exists( output_folder ):
	os.makedirs( output_folder )

filenames = glob.glob( input_folder + "*.fasta" )
for filename in filenames:
	print filename
	output_file = output_folder + filename.split('/')[-1]
	seqs = load_seqs( filename )
	with open( output_file, "w" ) as out:
		for key in seqs.keys():
			out.write( '>' + key + '\n' + seqs[ key ] + '\n' )
