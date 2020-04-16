



genes = []
input_file = "/vol/gf-yam/members/bpucker/20200119_BUSCO_of_final_annotaiton/busco/run_busco_run/full_table_busco_run.tsv"
with open( input_file, "r" ) as f:
	line = f.readline()
	while line:
		if line[0] != '#':
			parts = line.strip().split('\t')
			if parts[1] == "Duplicated":
				genes.append( parts[0] )
		line = f.readline()


for each in [ 2, 3, 4, 5, 6, 7, 8 ]:
	counter = []
	for gene in genes:
		if genes.count( gene ) == each:
			counter.append( gene )
	print str( each ) + " copies: " + str( len( list( set( counter ) ) ) )
