### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

import sys, glob, re, os, time, datetime, shutil

# --- end of imports --- #

__usage__ = """ python run_interpro5_on_cluster.py\n
				--output_dir <FULL_PATH_TO_OUTPUT_DIRECTORY>\n
				--input_file <INPUT_FILE (FASTA)>\n
				--output_file <OUTPUT_FILE_NAME (TSV)>\n
				--splitting_cutoff <INT, LENGTH_OF_SEQS_IN_EACH_PART_FILE>
				
				feature requests and bug reports:
				bpucker@cebitec.uni-bielefeld.de
"""


def submit_jobs_to_cluster( prefix, query_file_names, para_jobs ):
	"""! @brief submit InterProScan5 jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "I_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		
		cmd = "/vol/biotools/bin/interproscan --output-dir " + prefix + " --goterms --input " + file_name
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=10G",
																"-l arch=lx-amd64",
																"-P rapresmabs",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep( 2 )
		os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "I_" + batch_ID + "_\d{4}", content )
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
		qstat_IDs = re.findall( "I_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				for each in content.split('\n')[2:-1]:
					if ID in each.split()[2] and not 'd' in each.split()[4]:
						waiting_status = True
		time.sleep( 10 )


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


def produce_multiple_query_files( prefix, query_file, cutoff ):
	"""! @brief produce multiple query files """
	
	sequences = load_sequences( query_file )
	
	query_file_names = []
	
	len_counter = 0
	name_counter = 1
	query_file = prefix + "0".zfill(4) + ".fasta"
	query_file_names.append( query_file )
	out = open( query_file, "w" )
	for idx, seq_id in enumerate( sorted( sequences.keys() ) ):
		if len_counter >= cutoff:
			len_counter = 0
			out.close()
			query_file = prefix + str( name_counter ).zfill(4) + ".fasta"
			query_file_names.append( query_file )
			out = open( query_file, "w" )
			name_counter += 1
		out.write( '>' + seq_id + '\n' + sequences[ seq_id ] + '\n' )
		len_counter += len( sequences[ seq_id ] )
	out.close()
	return query_file_names


def combine_files( input_dir, output_file ):
	"""! @brief combine all TSV result files """
	
	filenames = glob.glob(  input_dir + "*.tsv" )
	with open( output_file, "w", 0 ) as out:
		for filename in filenames:
			with open( filename, "r" ) as f:
				out.write( f.read() )


def main( arguments ):
	"""! @brief check inputs and call functions """
	
	t1 = datetime.datetime.now()
	if '--para_jobs' in arguments:
		para_jobs = int( arguments[ arguments.index( '--para_jobs' ) + 1 ] )
	else:
		para_jobs = 300
		
	if '--splitting_cutoff' in arguments:
		cutoff = int( arguments[ arguments.index( '--splitting_cutoff' ) + 1 ] )
	else:
		cutoff=5000
	
	query_file = arguments[ arguments.index( '--input_file' ) + 1 ]
	prefix = arguments[ arguments.index( '--output_dir' ) + 1 ]
	if not prefix[-1] == "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	output_file = arguments[ arguments.index( '--output_file' ) + 1 ]
	
	
	query_file_names = produce_multiple_query_files( prefix, query_file, cutoff )
	
	submit_jobs_to_cluster( prefix, query_file_names, para_jobs )
	
	combine_files( prefix, output_file )
	
	t2 = datetime.datetime.now()


if __name__ == '__main__':
	
	if '--input_file' in sys.argv and '--output_dir' in sys.argv and '--output_file' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
