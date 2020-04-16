

##### --- FINAL PAPER FIGURES --- #####



### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###
### based on a script included in: https://www.mdpi.com/2073-4425/10/9/671 ###

__usage__ = """
					python cov_plot_per_contig.py
					--in <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION>
					--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import sys, os, re
import matplotlib.pyplot as plt
import numpy as np

# --- end of imports --- #

def generate_plot( coverage, out_file, resolution, saturation ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots( figsize=( 10, 5 ) )
	
	max_value = 0
	collected_values = {}
	
	# --- generate list for plotting --- #
	x = []
	blocks = [ coverage[ i : i + resolution ] for i in xrange( 0, len( coverage ), resolution ) ]
	for block in blocks:
		x.append( min( [ np.mean( block ), saturation ] ) )
	max_value = max( [ max_value, max( x ) ] )
	
	# --- plot values --- #
	max_value = float( min( [ saturation, max_value ] ) )
	
	x_norm = []
	for each in x:
		x_norm.append( min( [ 1, ( each / max_value ) ] ) )
	
	ax.scatter( np.arange( 0, len( x ), 1 ), x_norm, s=1, color="lime" )
	
	ax.plot( [ 0, len( x ) ], [ 0 / max_value, 0 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 50 / max_value, 50 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 100 / max_value, 100 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 150 / max_value, 150 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 200 / max_value, 200 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 250 / max_value, 250 / max_value ], color="black" , linewidth=0.1)
	ax.plot( [ 0, len( x ) ], [ 300 / max_value, 300 / max_value ], color="black" , linewidth=0.1)
	
	ax.plot( [ 0, 0 ], [ 0, 1 ], color="black", linewidth=1, markersize=1 )
	ax.text( 0, 1, str( int( max_value ) ), ha="right", fontsize=5 )
	ax.text( 0, 0.5, str( int( max_value / 2 ) ), ha="right", fontsize=5 )
	ax.text( 0, 0, "0", ha="right", fontsize=5 )
	
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	ax.set_ylabel( "coverage" )
	
	ax.set_xlim( 0, len( x ) )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 10
	
	ax.xaxis.set_ticks( np.arange( 0, len( x ), 100 ) )
	labels = map( str, np.arange( 0, 10, 1 ) )
	ax.set_xticklabels( labels )	#[ "0", "5", "10", "15", "20", "25", "30" ]
	
	plt.subplots_adjust( left=0.03, right=0.999, top=0.99, bottom=0.1 )
	
	fig.savefig( out_file, dpi=300 )
	plt.close("all")


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file = arguments[ arguments.index( '--in' ) + 1 ]
	outdir = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--res' in arguments:
		resolution = int( arguments[ arguments.index( '--res' ) + 1 ] )
	else:
		resolution = 10000
	
	if '--sat' in arguments:
		saturation = int( arguments[ arguments.index( '--sat' ) + 1 ] )
	else:
		saturation = 300
	
	size_cutoff = 1000000
	
	if not os.path.exists( outdir ):
		os.makedirs( outdir )
	
	# --- generate per chromosome position coveage plot --- #
	with open( cov_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				if int( re.findall( "\d+", header )[-1] ) > size_cutoff:
					outputfile = outdir + header + ".pdf"
					generate_plot( tmp, outputfile, resolution, saturation )
				header = parts[0]
				tmp = []
			tmp.append( float( parts[-1] ) )
			line = f.readline()


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
