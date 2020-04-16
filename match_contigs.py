


from operator import itemgetter
import re, sys, os

# --- end of imports --- #


def load_genes2contigs( gff_file ):
	"""! @brief load contig per gene """
	
	genes2contigs = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					genes2contigs.update( { parts[-1].split('=')[1]: parts[0] } )
			line = f.readline()
	return genes2contigs


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				content.update( { header: "".join( seq ) } )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		content.update( { header: "".join( seq ) } )
	return content


def load_RBHs( input_file ):
	"""! @brief load RBHs from given input file """
	
	RBHs = {}
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			RBHs.update( { parts[0]: parts[1] } )
			RBHs.update( { parts[1]: parts[0] } )
			line = f.readline()
	return RBHs


def get_best_match( all_contig_matches ):
	"""! @brief get best match based on list of matched contigs """
	
	candidates = list( set( all_contig_matches ) )
	freq_for_sorting = []
	for candidate in candidates:
		freq_for_sorting.append( { 'contig': candidate, 'counts': all_contig_matches.count( candidate ) } )
	return sorted( freq_for_sorting, key=itemgetter('counts') )[-1]


def match_allelic_contigs( contigs, genes2contigs, RBHs ):
	""""! @brief assign contigs to each other based on number of RBHs """
	
	contig_matches = {}
	for gene in RBHs.keys():
		contig = genes2contigs[ gene ]
		try:
			contig_matches[ contig ].append( genes2contigs[ RBHs[ gene ] ] )
		except KeyError:
			contig_matches.update( { contig: [ genes2contigs[ RBHs[ gene ] ] ] } )
	
	best_match = {}
	for key in contig_matches.keys():
		best_match.update( { key: get_best_match( contig_matches[ key ] ) } )
	return best_match


def get_genes_per_contig( genes2contigs ):
	"""! @brief count genes per contig """
	
	genes_per_contig = {}
	for gene in genes2contigs.keys():
		try:
			genes_per_contig[ genes2contigs[ gene ] ] += 1
		except KeyError:
			genes_per_contig.update( { genes2contigs[ gene ]: 1 } )
	return genes_per_contig


def integrate_data( contig_matches,  genes_per_contig, contigs, doc_file, final_result_file ):
	"""! @brief integrate all information """
	
	percent_cutoff = 30
	white_list = []
	
	# --- writing data into doc file --- #
	with open( doc_file, "w" ) as out:
		for contig in contig_matches.keys():
			new_line = [ 	contig,
									contig_matches[ contig ]['counts'],
									genes_per_contig[ contig ],
									100.0*contig_matches[ contig ]['counts'] / genes_per_contig[ contig ],
									len( contigs[ contig ] ) 
								]
			out.write( "\t".join( map( str, new_line ) ) + '\n' )
			
			# --- selecting surviving contigs --- #
			if 100.0*contig_matches[ contig ]['counts'] / genes_per_contig[ contig ] < percent_cutoff:
				white_list.append( contig )
	
	# --- generate final assembly file --- #
	total_size = 0
	with open( final_result_file, "w" ) as out:
		for contig in white_list:
			total_size += len( contigs[ contig ] )
			out.write( '>' + contig + '\n' + contigs[ contig ] + '\n' )
	
	print "Summary: percent_cutoff=" + str( percent_cutoff ) + "; final assembly size=" + str( total_size )
	return white_list


def get_haplo_peps( peps, haplo_contigs, pep_output_file ):
	"""! @brief get all peptides encoded on haploid contigs """
	
	counter = 0
	with open( pep_output_file, "w" ) as out:
		for key in peps.keys():
			contig = re.findall( "contig\d+", key )[0]
			if contig in haplo_contigs:
				out.write( '>' + key + '\n' + peps[ key ] + '\n' )
				counter += 1
	print "number of surviving peptides: " + str( counter )


input_file = "/vol/gf-yam/members/bpucker/20200327_haplophase_purging/RBH_file.txt"
gff_file = "/vol/gf-yam/members/bpucker/20200327_haplophase_purging/Dioscorea_dumetorum_v1.0.gff"
fasta_file = "/vol/gf-yam/members/bpucker/20200327_haplophase_purging/Dioscorea_dumetorum_v1.0.fasta"
pep_file = "/vol/gf-yam/results/manuscript/annotation_NEW/20191207_final.aa"

output_folder = "/vol/gf-yam/members/bpucker/20200327_haplophase_purging/"
doc_file = output_folder + "20200328_docs.txt"
final_result_file = output_folder + "20200328_yam_haploid.fasta"
pep_output_file = output_folder + "20200328_haplo_peps.fasta"


RBHs = load_RBHs( input_file )
genes2contigs = load_genes2contigs( gff_file )
contigs = load_multiple_fasta_file( fasta_file )

contig_matches = match_allelic_contigs( contigs, genes2contigs, RBHs )
genes_per_contig = get_genes_per_contig( genes2contigs )

haplo_contigs = integrate_data( contig_matches,  genes_per_contig, contigs, doc_file, final_result_file )
peps = load_multiple_fasta_file( pep_file )
get_haplo_peps( peps, haplo_contigs, pep_output_file )
