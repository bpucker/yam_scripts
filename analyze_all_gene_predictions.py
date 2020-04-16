### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###


__usage__ = """
	python analyze_all_gene_predictions.py\n
	--in <FULL_PATH_TO_INPUT_DIRECTORY>
	--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
	
	bug and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import glob, re, sys, os
import numpy as np

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = re.findall( "tig\d+\.g\d+\.t\d+", f.readline()[1:].strip() )[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = re.findall( "tig\d+\.g\d+\.t\d+", line[1:].strip() )[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def identify_representative_peptides( aa_file ):
	"""! @brief identify representative peptides """
	
	covered_genes = {}
	sequences = {}
	
	with open( aa_file ) as f:
		header = re.findall( "tig\d+\.g\d+\.t\d+", f.readline()[1:].strip() )[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
				parent = re.findall( "tig\d+\.g\d+", header )[0]
				if not 'X' in seq:
					if not '*' in seq:
						try:
							if len( seq ) > len( sequences[ covered_genes[ parent ] ] ):
								del sequences[ covered_genes[ parent ] ]
								del covered_genes[ parent ]
								sequences.update( { header: seq } )
								covered_genes.update( { parent: header } )
						except KeyError:
							sequences.update( { header: seq } )
							covered_genes.update( { parent: header } )
				header = re.findall( "tig\d+\.g\d+\.t\d+", line[1:].strip() )[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		parent = re.findall( "tig\d+\.g\d+", header )[0]
		if not 'X' in seq:
			if not '*' in seq:
				try:
					if len( seq ) > len( sequences[ covered_genes[ parent ] ] ):
						del sequences[ covered_genes[ parent ] ]
						del covered_genes[ parent ]
						sequences.update( { header: seq } )
						covered_genes.update( { parent: header } )
				except KeyError:
					sequences.update( { header: seq } )
					covered_genes.update( { parent: header } )
	return sequences


def get_gene_length_and_exon_numbers( gff_file, repr_peps ):
	"""! @brief get lengths and number of exons of all genes in given gff3 file """
	
	valid_genes = {}
	for key in repr_peps.keys():
		valid_genes.update( { ".".join( key.split('.')[:-1] ): None } )
	
	gene_lengths = []
	exon_number = []
	exon_counter = 0
	prev_gene = ""
	
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					try:
						valid_genes[ prev_gene ]
						if exon_counter != 0:
							exon_number.append( exon_counter )
					except KeyError:
						pass
					
					try:
						valid_genes[ re.findall( "tig\d+\.g\d+", line.strip().split('ID=')[-1]  )[0] ]
						gene_lengths.append( int( parts[4] )-int( parts[3] ) )
					except KeyError:
						pass
					exon_counter = 0
					prev_gene = re.findall( "tig\d+\.g\d+", line.strip().split('ID=')[-1]  )[0]
				elif parts[2] == "exon":
					exon_counter += 1
			line = f.readline()
		try:
			if valid_genes[ prev_gene ]:
				exon_number.append( exon_counter )
		except:
			pass
	return gene_lengths, exon_number, valid_genes


def construct_clean_gff( filename, clean_gff, valid_genes ):
	"""! @brief generate new GFF3 file with all valid genes """
	
	with open( clean_gff, "w" ) as out:
		with open( filename, "r" ) as f:
			line = f.readline()
			while line:
				try:
					gene = re.findall( "tig\d+\.g\d+", line )[0]
					try:
						valid_genes[ gene ]
						out.write( line )
					except KeyError:
						pass
				except IndexError:
					pass
				line = f.readline()


def main( parameters ):
	"""! @brief calculate stats for gene predictions and generate files with representative sequences """
	
	input_dir = parameters[ parameters.index('--in')+1 ]
	output_dir = parameters[ parameters.index('--out')+1 ]
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	documentation_file = output_dir + "gene_prediction_stats.txt"
	
	gff3_files = glob.glob( input_dir + "*.gff" ) + glob.glob( input_dir + "*/*.gff" ) + glob.glob( input_dir + "*/*/*.gff" )
	print "number of gff3 files to process: " + str( len( gff3_files ) )
	with open( documentation_file, "w", 0 ) as doc:
		doc.write( "ID\tgene_number\tgene_length\texons_per_gene\taverage_mRNA_length\taverage_peptide_length\n" )
		for filename in sorted( gff3_files ):
			ID = filename.split('/')[-2]
			try:
				aa_file = filename.replace( ".gff", ".aa" )
				mrna_file = filename.replace( ".gff", ".mrna" )
				cds_file =  filename.replace( ".gff", ".codingseq" )
				
				repr_peps = identify_representative_peptides( aa_file )
				mrnas = load_sequences( mrna_file )
				cds = load_sequences( cds_file )
				repr_mrnas = {}
				repr_cds = {}
				
				for key in repr_peps.keys():
					repr_mrnas.update( { key: mrnas[ key ] } )
					repr_cds.update( { key: cds[ key ] } )
				
				# --- construct repr peptide and mRNA files --- #
				repr_mrna_file = output_dir + ID + ".repr.mrna.fasta"
				with open( repr_mrna_file, "w" ) as out:
					for key in sorted( repr_mrnas.keys() ):
						out.write( '>' + key + '\n' + repr_mrnas[ key ] + '\n' )
				repr_pep_file = output_dir + ID + ".repr.pep.fasta"
				with open( repr_pep_file, "w" ) as out:
					for key in sorted( repr_peps.keys() ):
						out.write( '>' + key + '\n' + repr_peps[ key ] + '\n' )
				repr_cds_file = output_dir + ID + ".repr.cds.fasta"
				with open( repr_cds_file, "w" ) as out:
					for key in sorted( repr_cds.keys() ):
						out.write( '>' + key + '\n' + repr_cds[ key ] + '\n' )
				
				mrna_lengths = []
				for seq in repr_mrnas.values():
					mrna_lengths.append( len( seq ) )
				pep_lengths = []
				for seq in repr_peps.values():
					pep_lengths.append( len( seq ) )
				
				gene_lengths, exon_number, valid_genes = get_gene_length_and_exon_numbers( filename,  repr_peps)
				
				clean_gff = output_dir + ID + ".repr.gff"
				construct_clean_gff( filename, clean_gff, valid_genes )
				
				doc.write( "\t".join( map( str, [ 	ID,
																	len( gene_lengths ),
																	np.mean( gene_lengths ),
																	np.mean( exon_number ),
																	np.mean( mrna_lengths ),
																	np.mean( pep_lengths ) 
																] ) ) + '\n' )
			except:
				print ID


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
