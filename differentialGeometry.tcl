set SCRIPTS "/home/marcario/scripts/"

source ${SCRIPTS}/nanodisc/nanodiscAnalysis.tcl

proc printArrayValuesScreen { arrayName } { 
	# Takes an array as input and prints
	# the keys and values of the array
	# to the screen.
	# "arrayName" A Tcl array structure

	upvar $arrayName arrName
	foreach { key value } [ array get arrName ] {
		puts "Key is $key and associated value is $value"
		unset key value
	}
}

proc printArrayValuesOutfile { arrayName outFile } {
	# Takes an array as input and prints 
	# the keys and values of the array to 
	# a txt file for later analysis/viz.
	# "arrayName" : A Tcl array structure 
	# "outFile" : The desired name of the 
	# output file. It will have a 'txt'
	# extension by default (for now).

	upvar $arrayName arrName
	set outfile [ open "$outFile.txt" w+ ]	
	foreach { key value } [ array get arrName ] {
		puts $outfile "$key,$value"
		unset key value
	}
	
	close $outfile
}

proc makeMongeSurface { radialResolution radialMax numThetaBins atomSelection referenceSelection outFilePrefix numFrames } { 
	# Takes the input nanodisc system and 
	# creates a monge patch in polar coordinates

	global M_PI
	
	# Creates array of appropriate size for measuring membrane heights
	# and creates a second array to store center of each patch.
	set numRadialBins [ expr { int( ceil( $radialMax / $radialResolution ) ) } ]
	set thetaResolution [ expr { 2 * ( $M_PI / $numThetaBins ) } ]
	array set upperBins [ initialize2DArray $numRadialBins $numThetaBins ]
	array set upperBinCounter [ initialize2DArray $numRadialBins $numThetaBins ]
	array set lowerBins [ initialize2DArray $numRadialBins $numThetaBins ]
	array set lowerBinCounter [ initialize2DArray $numRadialBins $numThetaBins ]
	
	# Fit nanodisc to xy-plane. Probably needs improvement, because right
	# now I just fit to first frame.
	alignToFirstFrame $referenceSelection $numFrames
	
	set lipids [ atomselect top $atomSelection ]
	set lipidIndices [ $lipids get index ]
	for { set i 0 } { $i < $numFrames } { incr i 1 } {
		$lipids frame $i 
		set membraneMidplane [ getMembraneMidplane $lipids ]
		set leaflets [ leafletSorter $lipids $membraneMidplane ]
		puts "Membrane midplane is $membraneMidplane"
		puts "Creating parameterized surface for frame $i"
		foreach lipidIndex $lipidIndices {
			set lipidFrameI [ atomselect top "index $lipidIndex" frame $i ]
			set cartesianCoordinates [ lindex [ $lipidFrameI get { x y } ] 0 ]
			set polarCoordinates [ cartToPolar2D [ lindex $cartesianCoordinates 0 ] [ lindex $cartesianCoordinates 1 ] ]
			set radialBin [ expr { int( [ lindex $polarCoordinates 0 ] / $radialResolution ) } ]
			set thetaBin [ expr { int( [ lindex $polarCoordinates 1 ] / $thetaResolution ) } ]
			#puts "Lipid atom is at [ lindex $polarCoordinates 0 ], [lindex $polarCoordinates 1] in the r,theta-plane \
			which is bin $radialBin, $thetaBin."
			#puts "Bin value is $bins($radialBin,$thetaBin)"
			if { [ lsearch [ lindex $leaflets 0 ] $lipidIndex ] != -1 } {
				set upperBins($radialBin,$thetaBin) [ expr { $upperBins($radialBin,$thetaBin) + [ getMembraneHeight $lipidIndex $membraneMidplane $i ] } ]
				set upperBinCounter($radialBin,$thetaBin) [ expr { $upperBinCounter($radialBin,$thetaBin) + 1  } ]
			} else {
				set lowerBins($radialBin,$thetaBin) [ expr { $lowerBins($radialBin,$thetaBin) + [ getMembraneHeight $lipidIndex $membraneMidplane $i ] } ]
				set lowerBinCounter($radialBin,$thetaBin) [ expr { $lowerBinCounter($radialBin,$thetaBin) + 1  } ]
			}
			$lipidFrameI delete
			unset cartesianCoordinates polarCoordinates
			unset radialBin thetaBin
		}
		unset membraneMidplane
		unset lipidIndex
		
	}
	$lipids delete
	unset lipidIndices i

	# Averages the heights over the length of the simulation
	for { set i 0 } { $i < $numRadialBins } { incr i 1 } { 
		for { set j 0 } { $j < $numThetaBins } { incr j 1 } {
			if { !$upperBinCounter($i,$j) } {
				set upperBins($i,$j) 0.0
			} else {
				set avgBinHeight [ expr { $upperBins($i,$j) / $upperBinCounter($i,$j) } ]
				set upperBins($i,$j) $avgBinHeight
			}
			if { !$lowerBinCounter($i,$j) } {
				set lowerBins($i,$j) 0.0
			} else {
				set avgBinHeight [ expr { $lowerBins($i,$j) / $lowerBinCounter($i,$j)  } ]
				set lowerBins($i,$j) $avgBinHeight
			}
		}
		unset -nocomplain j avgBinHeight
	}
	unset i

	printArrayValuesOutfile upperBins "$outFilePrefix.upper"
	printArrayValuesOutfile lowerBins "$outFilePrefix.lower"
	#return [ array get bins ]
}
