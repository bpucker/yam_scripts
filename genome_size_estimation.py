### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python genome_size_estimation.py\n
					--output_dir <FULL_PATH_TO_CLUSTER_DIRECTORY>
					--fw <FULL_PATH_TO_FW_READS_FILE>
					--rv <FULL_PATH_TO_RV_READS_FILE>\n
					--kmer <INT, length of kmer>
					--read_len <INT, read length [nt]>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #
def count_kmers( jellyfish, kmer, read1_file, read2_file, jf_file ):
	"""! @brief construct count table for all seen k-mers """
	
	uncompressed_read1_file = read1_file#[:-3]
	uncompressed_read2_file = read2_file#[:-3]
	
	#os.popen( "gunzip " + read1_file )
	#os.popen( "gunzip " + read2_file )
	
	cmd = [ jellyfish, "count -C -m", str( kmer ), "-s 1000000000 -t 10", uncompressed_read1_file, uncompressed_read2_file, "-o", jf_file ]
	
	os.popen( " ".join( cmd ) )
	
	#os.popen( "gzip " + uncompressed_read1_file )
	#os.popen( "gzip " + uncompressed_read2_file )


def construct_histo( jellyfish, jf_file, histo_file ):
	"""! @brief construct histogram based on k-mer count table """
	
	cmd = [ jellyfish, "histo -t 10", jf_file, ">", histo_file ]
	
	os.popen( " ".join( cmd ) )


def run_final_analysis( Rscript, genome_scope, cluster_dir, histo_file, kmer, read_length, output_file ):
	"""! @brief run final evaluation to estimate the genome size """
	
	cmd = [ Rscript, genome_scope, histo_file, str( kmer ), str( read_length ), cluster_dir, ">", output_file ]
	
	os.popen( " ".join( cmd ) )	


def main( arguments ):
	"""! @brief controls pathway of genome size calculation """
	
	cluster_dir = arguments[ arguments.index('--output_dir')+1 ]
	read1_file = arguments[ arguments.index('--fw')+1 ]
	read2_file = arguments[ arguments.index('--rv')+1 ]
	
	if cluster_dir[-1] != '/':
		cluster_dir += "/"
	
	if not os.path.exists( cluster_dir ):
		os.makedirs( cluster_dir )
	
	if '--kmer' in arguments:
		kmer = int( arguments[ arguments.index('--kmer')+1 ] )
	else:
		kmer=21
	
	if '--read_len' in arguments:
		read_length = int( arguments[ arguments.index('--read_len')+1 ] )
	else:
		read_length = 250
	
	jellyfish = "/vol/biotools/bin/jellyfish"
	genome_scope = "/vol/cluster-data/bpucker/bin/genomescope-master/genomescope.R"
	Rscript = "/usr/bin/Rscript"
	
	jf_file = cluster_dir + "reads.jf"
	count_kmers( jellyfish, kmer, read1_file, read2_file, jf_file )
	
	histo_file = cluster_dir + "reads.histo"
	construct_histo( jellyfish, jf_file, histo_file )
	
	output_file = cluster_dir + "log.txt"
	run_final_analysis( Rscript, genome_scope, cluster_dir, histo_file, kmer, read_length, output_file )
	
	#os.popen( "rm " + histo_file )
	#os.popen( "rm " + jf_file )


if __name__ == '__main__':
	
	if '--output_dir' in sys.argv and '--fw' in sys.argv and '--rv' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
