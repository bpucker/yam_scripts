import os, re

# --- end of imports --- #


def load_phased_genes( phasing_status_file, min_cov, max_cov ):
	"""! @brief identify phased genes based on their average coverage """
	
	phased_genes = {}
	with open( phasing_status_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if min_cov < float( parts[1] ) < max_cov:
				phased_genes.update( { parts[0]: None } )
			line = f.readline()
	return phased_genes


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences



min_cov = 50
max_cov = 150

phasing_status_file = "/vol/gf-yam/members/bpucker/20200130_match_alleles/phasing_status.txt"
input_cds_file = "/vol/gf-yam/results/manuscript/annotation_NEW/published_annotation/Dioscorea_dumetorum_v1.0.codingseq"
outpt_cds_file = "/vol/gf-yam/members/bpucker/20200130_match_alleles/phased_cds.fasta"


phased_genes = load_phased_genes( phasing_status_file, min_cov, max_cov )
seqs = load_sequences( input_cds_file )

with open( outpt_cds_file, "w" ) as out:
	for key in seqs.keys():
		ID = re.findall( "contig\d+\.g\d+", key )[0]
		try:
			phased_genes[ ID ]
			if key[-2:] == "t1":
				out.write( '>' + ID + '\n' + seqs[ key ] + '\n' )
		except KeyError:
			pass
