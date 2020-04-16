import re

# --- end of imports --- #

fasta_size_file = "/vol/gf-yam/members/bpucker/20200214_in_silico_map/Dioscorea_dumetorum_v1.0.fasta.sizes"
agp_file = "/vol/gf-yam/members/bpucker/20200214_in_silico_map/genetic_map.agp"
unanchored_agp_file = "/vol/gf-yam/members/bpucker/20200214_in_silico_map/genetic_map.unplaced.agp"

contig_sizes = {}
with open( fasta_size_file, "r" ) as f:
	line = f.readline()
	while line:
		parts = line.strip().split('\t')
		contig_sizes.update( { parts[0]: int( parts[1] ) } )
		line = f.readline()

with open( agp_file, "r" ) as f:
	content = f.read()
contigs = list( set( re.findall( "contig\d+", content ) ) )

with open( unanchored_agp_file, "r" ) as f:
	content = f.read()
unplaced_contigs = list( set( re.findall( "contig\d+", content ) ) )
placed_contigs = []
for contig in contigs:
	if contig not in unplaced_contigs:
		placed_contigs.append( contig )

print "number of anchored contigs: " + str( len( placed_contigs ) )
print "number of unplaced contigs: " + str( len( unplaced_contigs ) )

placed_contig_length = 0
for contig in placed_contigs:
	placed_contig_length += contig_sizes[ contig ]

unplaced_contig_length = 0
for contig in unplaced_contigs:
	unplaced_contig_length += contig_sizes[ contig ]


print "combined length of placed contigs: " + str( placed_contig_length )
print "combined length of unplaced contigs: " + str( unplaced_contig_length )

