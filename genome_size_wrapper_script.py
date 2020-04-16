import sys, glob, re, os, time, datetime, shutil

# --- end of imports --- #

def submit_jobs_to_cluster( cluster_dir_names, read_files, kmers, para_jobs, script_name, read_lengths ):
	"""! @brief submit genome size estimation jobs to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, cluster_dir_name in enumerate( cluster_dir_names ):
		ID = "G_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = cluster_dir_names[ idx ] + ID + '.sh'
		out_file = cluster_dir_names[ idx ] + ID + '.out'
		err_file = cluster_dir_names[ idx ] + ID + '.err'
		
		cmd = "python " + script_name + " --fw " + read_files[ idx ][0]  + " --rv " + read_files[ idx ][1] + " --output_dir " +  cluster_dir_name + " --kmer " + str( kmers[ idx ] ) + " --read_len " + str( read_lengths[ idx ] )
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-pe multislot 10",
																"-N",
																ID,
																"-l vf=5G",
																"-l arch=lx-amd64",
																"-P denbi",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "G_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "G_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 120 )


def load_results_from_summary_file( sum_file ):
	"""! @brief load all results from summary file """
	
	stats = {}
	
	with open( sum_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == 'k':
				stats.update( { 'k': re.findall( "\d+", line )[0] } )
			elif line[:len( "Hetero" )] == "Hetero":
				try:
					min_hetero, max_hetero = re.findall( "\d+\.\d+%", line )
				except ValueError:
					min_hetero, max_hetero = "0%", "0%"
				stats.update( { 	'hetero_min': float( min_hetero[:-1] ),
											'hetero_max': float( max_hetero[:-1] ) } )
			elif line[:len( "Genome Haploid" )] == "Genome Haploid":
				try:
					min_hap_size, max_hap_size = re.findall( "[\d,]+", line )
				except ValueError:
					min_hap_size, max_hap_size = "0", "0"
				stats.update( { 	'hap_size_min': int( min_hap_size.replace( ',', '' ) ),
											'hap_size_max': int( max_hap_size.replace( ',', '' ) ) } )
			elif line[:len( "Genome Repeat" )] == "Genome Repeat":
				try:
					min_repeat, max_repeat = re.findall( "[\d,]+", line )
				except ValueError:
					min_repeat, max_repeat = "0", "0"
				stats.update( { 	'repeat_min': int( min_repeat.replace( ',', '' ) ),
											'repeat_max': int( max_repeat.replace( ',', '' ) ) } )
			elif line[:len( "Model Fit" )] == "Model Fit":
				try:
					min_model, max_model = re.findall( "\d+\.\d+%", line )
				except ValueError:
					min_model, max_model = "0%", "0%"
				stats.update( { 	'model_min': float( min_model[:-1] ),
											'model_max': float( max_model[:-1] ) } )
			elif line[:len( "Read" )] == "Read":
				try:
					min_read_err, max_read_err = re.findall( "\d+\.\d+%", line )
				except ValueError:
					min_read_err, max_read_err = "0%", "0%"
				stats.update( { 	'read_err_min': float( min_read_err[:-1] ),
											'read_err_max': float( max_read_err[:-1] ) } )
			line = f.readline()
	return stats


def load_avg_read_lens( info_table ):
	"""! @brief load average read length per accession """
	
	read_lengths = {}
	with open( info_table, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			read_lengths.update( { parts[0]: int( parts[1] ) } )
			line = f.readline()
	return read_lengths


def get_read_pairs( input_dir ):
	"""! @brief get all read file pairs """
	
	read_pairs = []
	subdirs = [ subdir[ 0 ] for subdir in os.walk( input_dir ) ]
	for subdir in subdirs:
		filenames = sorted( glob.glob( subdir + "/*.fastq" ) )
		if len( filenames ) == 2:
			read_pairs.append( filenames )
	return read_pairs


if __name__ == '__main__':
	
	info_table = "/vol/gf-yam/members/bpucker/20200214_genome_size_estimation/read_infos.txt"
	input_dir = "/vol/cluster-data/bpucker/yam/yam_reads/"
	
	prefix = "/vol/cluster-data/bpucker/yam/yam_genome_size3/"
	
	final_result_file = "/vol/gf-yam/members/bpucker/20200214_genome_size_estimation/20200222_yam_GenomeScope.txt"
	
	read_lengths = load_avg_read_lens( info_table )
	read_file_pairs = get_read_pairs( input_dir )
	
	# --- prepare everything and run jobs --- #
	cluster_dir_names = []
	kmers = []
	read_files = []
	read_lens = []
	potential_kmers = [ 19, 21, 23, 25 ]
	for kmer in potential_kmers:
		for read_file_pair in read_file_pairs:
			ID = read_file_pair[0].split('/')[-1].split('.')[0]
			dir_name = prefix + ID + "_" + str( kmer ) + "/"
			if not os.path.exists( dir_name ):
				os.makedirs( dir_name )
				cluster_dir_names.append( dir_name )
				kmers.append( kmer )
				read_files.append( read_file_pair )
				try:
					read_lens.append( read_lengths[ ID ] )
				except KeyError:
					read_lens.append( 250 )
	
	script_name = "/vol/gf-yam/members/bpucker/20200214_genome_size_estimation/genome_size_estimation.py"
	para_jobs = 50
	print "number of input file pairs: " + str( len( read_file_pairs ) )
	submit_jobs_to_cluster( cluster_dir_names, read_files, kmers, para_jobs, script_name, read_lens )
	
	# --- collect results --- #
	cluster_dir_names = [ x[0] for x in os.walk( prefix ) ]
	print "number of different samples: " + str( len( cluster_dir_names ) )
	collected_results = {}
	for ID in cluster_dir_names:
		sum_file = ID + "/summary.txt"
		try:
			stats = load_results_from_summary_file( sum_file )
			collected_results.update( { ID: stats } )
		except IOError:
			pass
	
	# --- construct output files --- #
	properties_of_interest = sorted( collected_results[ collected_results.keys()[0] ].keys() )
	with open( final_result_file, "w" ) as out:
		out.write( "SampleName\t" + "\t".join( properties_of_interest ) + '\n' )
		for key in sorted( collected_results.keys() ):
			new_line = [ key ]
			for prop in properties_of_interest:
				try:
					new_line.append( str( collected_results[ key ][ prop ] ) )
				except KeyError:
					new_line.append( "n/a" )
			out.write( "\t".join( new_line ) + '\n' )
		
