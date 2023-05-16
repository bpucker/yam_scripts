
def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:].split( " " )[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_contigs_of_interest( blast_hit_file ):
	""""! @brief load all contigs of interest """
	
	contigs = []
	with open( blast_hit_file, "r" ) as f:
		line =f.readline()
		while line:
			parts = line.strip().split('\t')
			if len( parts ) > 3:
				contigs.append( parts[1] )
			line = f.readline()
	return list( set( contigs ) )


assembly_file = "./20200223_plastome_chondrome/consensusPilonx3.fasta"
blast_hit_file = "./20200223_plastome_chondrome/plastome_vs_yam.txt"

output_file = "./20200223_plastome_chondrome/putative_plastome_contigs_final.fasta"

seqs = load_sequences( assembly_file )

conitgs_of_interest = load_contigs_of_interest( blast_hit_file )

with open( output_file, "w" ) as out:
	for contig in conitgs_of_interest:
		print contig + "\t" + str( len( seqs[ contig ] ) )
		if contig == "tig00001791_166.0_160131.0:1.0-159713.0_pilon_pilon_pilon":
		#if 100000 < len( seqs[ contig ] ) < 200000:
			out.write( '>' + contig + '\n' + seqs[ contig ] + '\n' )
	
