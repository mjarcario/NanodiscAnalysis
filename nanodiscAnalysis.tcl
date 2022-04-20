
set SCRIPTS "/home/marcario/scripts"

source $SCRIPTS/procedures.tcl

proc getAverageStandardDeviation { listOfNumbers } { 
	# Takes as input a list of numbers and 
	# returns the average and standard deviation 
	# of the data.
	# "list of Numbers" : A list of numbers which
	# will be manipulated.

	if { ![ llength $listofNumbers ] } {
		error "The list provided is empty and cannot calcualte the statistical measures."
	}

	set total 0.0; set total2 0.0
	set sum 0.0; set sum2 0.0
	set numDataPoints [ llength $listOfNumbers ]
	foreach dataPoint $listOfNumbers {
		set total [ expr { $sum + $dataPoint } ]
		set total2 [ expr { $sum2 + $datapoint**2 } ]
		unset dataPoint
	}
	set average [ expr { $sum / $numDataPoints } ]
	set average2 [ expr { $sum2 / $numDataPoints } ]
	# Handles sitation if variance is negative by setting standard deviation
	# to 0.0. Probably need a better fix later.
	if { $average2 >= [ expr { $average**2  } ] } {
		set standardDeviation [ expr { sqrt( $average2 - $average**2 )  } ]
	} else {
		puts "Standard deviation could not be calculated because "
		puts "E(X)^2 was greater than E(X^2). Do not trust the "
		puts "standard deviation of '0.0'."
		set standardDeviation 0.0
	}
	unset total total2 sum sum2
	unset numDataPoints average2

	return [ list $average $standardDeviation ]
}

proc getNanodiscDiameter { atomSelection numFrames outFileName} {
	# Obtains the instantaneous diameter of the nanodisc
	# and measures the trajectory average + std err

	# This assumes the nanodisc is circular and will output
	# the average radius. This may or may not be true. Plan 
	# on coding an elliptical fitting procedure to eliminate
	# this assumption.

	# This procedure will center the nanodisc assembly to
	# the origin for analysis, which may or may not impact 
	# later analyses. Use caution. Need a better way to 
	# orient so that the nanodisc lies in the xy-plane.

	# cartToPolar2D comes from ~/scripts/procedures.tcl and 
	# converts x,y coordinates to r,theta values in a two-
	# dimensional plane

	# "atomSelection" : The atoms to be included in calculating 
	# the radius of the nanodisc (i.e., "segname MSP1 MSP2 and 
	# name CA").
	# "numFrames" : The number of frames in the tracjetory to 
	# loop over.
	# "outFileName" : Name of the output file to write the inst-
	# antaneous diameter to as well as the average diameter at 
	# the end of the script.

	if { [ file exists ${outFileName}.txt ] } {
		puts "${outFileName}.txt exists. Will delete this file and write new data."
		file delete ${outFileName}.txt
	}
	set outFile [ open "${outFileName}.txt" w+ ] 

	puts -nonewline "Measuring the diameter of the nanodisc..."
	set all [ atomselect top "all" ]
	set nanodisc [ atomselect top "$atomSelection" ]
	set nanodiscAtomNum [ $nanodisc num ]
	set nanodiscAtoms [ $nanodisc get index ]
	set accumulatedDiameter 0.0; set accumulatedDiameter2 0.0

	# Loop over frames with each iteration calculating the estimated radius of the 
	# nanodisc based on the average distance of "atomSelection" from the origin
	for { set i 0 } { $i < $numFrames } { incr i 1 } {
		puts -nonewline "Analyzing frame $i..."
		$nanodisc frame $i
		$all frame $i
		$all moveby [ vecinvert [ measure center $nanodisc ] ]
		#set radiusValues 0.0
		foreach atom $nanodiscAtoms {
			set nanodiscAtom [ atomselect top "index $atom" frame $i ]
			set coordinates [ lindex [ $nanodiscAtom get {x y} ] 0 ]
			set polarCoordinates [ cartToPolar2D [ lindex $coordinates 0 ] [ lindex $coordinates 1 ] ]
			#set radiusValues [ expr { $radiusValues + [ lindex $polarCoordinates 0 ] } ]
			lappend radiusValues [ lindex $polarCoordinates 0 ]
			$nanodiscAtom delete
			unset coordinates polarCoordinates
		
		}
		#set frameRadius [ expr { $radiusValues / $nanodiscAtomNum } ]
		set frameRadius [ lindex [ getAverageStandardDeviation $radiusValues ] 0 ]
		set frameDiameter [ expr { 2 * $frameRadius } ]
		puts "Frame Diameter is $frameDiameter"
		puts $outFile "$i $frameDiameter"
		#set accumulatedDiameter [ expr { $accumulatedDiameter + $frameDiameter } ]
		#set accumulatedDiameter2 [ expr { $accumulatedDiameter2 + $frameDiameter**2 } ]
		lappend totalDiameter $frameDiameter
		unset radiusValues frameRadius frameDiameter
	}	
	unset i nanodiscAtomNum nanodiscAtoms
	$nanodisc delete
	$all delete

	# Calculation of average +/- standard deviation of diameter
	#set averageDiameter [ expr { $accumulatedDiameter / $numFrames } ]
	#set averageDiameter2 [ expr { $accumulatedDiameter2 / $numFrames } ]
	#puts "[expr {$averageDiameter**2}] $averageDiameter2"
	#set stdDevDiameter [ expr { sqrt( $averageDiameter2 - $averageDiameter**2 ) } ]
	set nanodiscStatistics [ getAverageStandardDeviation $totalDiameter ]
	set averageDiameter [ lindex $nanodiscStatistics 0 ]
	set stdDevDiameter [ lindex $nanodiscStatistics 1 ]
	puts "Nanodisc diameter is $averageDiameter +/- $stdDevDiameter Angstroms"
	puts $outFile "# $averageDiameter $stdDevDiameter"
	close $outFile
	
	unset accumulatedDiameter accumulatedDiameter2
	unset averageDiameter averageDiameter2 stdDevDiameter
}

proc leafletSorter { atomSelection } {
	# Sorts lipids into upper and lower leaflets
	# and returns a list of lists containing 1. a list of 
	# upper leaflet lipids and 2. a list lower leaflet lipids.

	# "atomSelection" : A (currently) one atom selection to use 
	# as the determinant for being either above or below membrane 
	# midplane.

	puts "Sorting lipids into upper and lower leaflets"
	set upperLipids ""
	set lowerLipids ""
	set lipids [ atomselect top $atomSelection ]
	set lipidNum [ $lipids num ]
	set lipidIndices [ $lipids get index ]
	set membraneCenter [ lindex [ measure center $lipids ] 2 ]
	foreach index $lipidIndices {
		set sortingAtom [ atomselect top "index $index" ]
		set atomPosition [ $sortingAtom get {z} ]
		set atomResid [ $sortingAtom get resid ]
		if { $atomPosition > $membraneCenter } {
			lappend upperLipids $atomResid
		} else {
			lappend lowerLipids $atomResid
		}
		$sortingAtom delete
		unset atomPosition atomResid
	}
	$lipids delete
	unset lipidNum lipidIndices membraneCenter

	puts "Done sorting lipids. Enjoy Coke!"	

	return [ list $upperLipids $lowerLipids ]
}


proc initialize2DArray { numI numJ } {
	# Initialzes a two-dimensional array
	# for use in sotring data

	for { set i 0 } { $i < $numI } { incr i 1 } {
		for { set j 0 } { $j < $numJ } { incr j 1 } {
			set matrix($i,$j) 0.0
		}
	}
	unset i j
	
	return [ array get matrix ]
}

proc initialize1DArray { numI } { 
	# Intializes a one-dimensional array 
	# for use in sorting data
	
	for { set i 0 } { $i < $numI } { incr i 1 } { 
		set matrix($i) 0.0
	}
	unset i

	return [ array get matrix ]
}

proc getAtomsDistanceFromProtein { atomSelection proteinSelection innerRadius outerRadius frameNum } { 
	# Obtains a list of atoms between the inner and outer radii
	# defined relative to the protein surface. Returns a list 
	# of indices containing the atoms between inner and outer 
	# radii.

	# "atomSelection" : Atoms that will be found relative 
	# to protein surface.
	# "proteinSelection" : Which atoms to use as the protein
	# surface (i.e., "backbone", "name CA", "sidechains", "all").
	# "innerRadius" : The lower bound of the radial surface 
	# in which you are looking for atoms.
	# "outerRadius" : The upper bound of the radial surface
	# in which you are looking for atoms

	if { [ string match all $proteinSelection ] } { 
		set innerSelection [ atomselect top "$atomSelection and ( within $innerRadius of protein )" frame $frameNum ]
		set outerSelection [ atomselect top "$atomSelection and ( within $outerRadius of protein )" frame $frameNum ]
	} else { 
		set innerSelection [ atomselect top "$atomSelection and ( within $innerRadius of ( protein and $proteinSelection ) )" frame $frameNum ]
		set outerSelection [ atomselect top "$atomSelection and ( within $outerRadius of ( protein and $proteinSelection ) )" frame $frameNum ]
	}

	set innerSelectionIndices [ $innerSelection get index ]
	set outerSelectionIndices [ $outerSelection get index ]
	$innerSelection delete
	$outerSelection delete

	foreach value $outerSelectionIndices {
		set matchValue [ lsearch $innerSelectionIndices $value ]
		if { $matchValue != -1 } { 
			set listIndex [ lsearch -exact $outerSelectionIndices $value ]
			set outerSelectionIndices [ lreplace $outerSelectionIndices $listIndex $listIndex ]
			unset listIndex
			}
		
	}
	unset innerSelectionIndices

	return $outerSelectionIndices
}

proc getMembraneMidplane { atomSelection } {
	# Obtains the z-coordinate for the membrane midplane
	# based on the geometric center of the atoms.
	# Returns a number containing the z-coordinate of 
	# membrane midplane for the given atomselection.
	# "atomSelection" : The atoms to use to determine 
	# the gemometric center of the membrane
	# "frame" : The frame of the trajectory to be analyzed

	# LOOK! Comments are longer than the code.

	set center [ measure center $membraneSelection ]
	return [ lindex $center 2 ]

}


proc membraneHeightFromProteinAveragedOverTheta { atomSelection nanodiscSize binResolution numFrames } {
	# Calculates the average height of an atomselection 
	# as a function of distance from the protein
	# surface averaged over all theta.

	set numBins [ expr { $distanceFromProtein / $binResolution  } ]
	set heights [ initialize1DArray $numBins ]
 
	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		
	}
}
