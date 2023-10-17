
set SCRIPTS "/home/marcario/scripts"

source $SCRIPTS/procedures.tcl

proc getAverageStandardDeviation { listOfNumbers } { 
	# Takes as input a list of numbers and 
	# returns the average and standard deviation 
	# of the data.
	# "list of Numbers" : A list of numbers which
	# will be manipulated.

	if { ![ llength $listOfNumbers ] } {
		error "The list provided is empty and cannot calcualte the statistical measures."
	}

	set sum 0.0; set sum2 0.0
	set numDataPoints [ llength $listOfNumbers ]
	foreach dataPoint $listOfNumbers {
		set sum [ expr { $sum + $dataPoint } ]
		set sum2 [ expr { $sum2 + $dataPoint**2 } ]
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
	unset sum sum2
	unset numDataPoints average2

	return [ list $average $standardDeviation ]
}

proc getNanodiscDiameter { atomSelection numFrames outFileName} {
	# CANNOT RECOMMEND USE OF THIS PROCEDURE AS NANODISC
	# IS NOT CIRCULAR, BUT ELLIPSOID	

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

	alignToFirstFrame $atomSelection $numFrames

	puts "Measuring the diameter of the nanodisc..."
	set all [ atomselect top "all" ]
	set nanodisc [ atomselect top "$atomSelection" ]
	set nanodiscAtomNum [ $nanodisc num ]
	set nanodiscAtoms [ $nanodisc get index ]
	#set accumulatedDiameter 0.0; set accumulatedDiameter2 0.0

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
		#puts $radiusValues
		set frameRadius [ lindex [ getAverageStandardDeviation $radiusValues ] 0 ]
		set frameDiameter [ expr { 2 * $frameRadius } ]
		if { 1 } {
			set radiusValues [ lsort -increasing $radiusValues ]
			set radiusValuesLength [ llength $radiusValues ]
			set startingIndex [ expr { int( floor( $radiusValuesLength * 0.5  ) ) } ]  ; # Gets the highest 1/4 of radii for calculation
			for { set j $startingIndex } { $j < $radiusValuesLength } { incr j 1 } {
				lappend newRadiusValues [ lindex $radiusValues $j ]
			}
			unset j radiusValuesLength startingIndex
			set newFrameRadius [ lindex [ getAverageStandardDeviation $newRadiusValues ] 0 ]
			set newFrameDiameter [ expr { 2 * $newFrameRadius  } ]
		}
		puts "Frame Diameter is $frameDiameter $newFrameDiameter"
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
	
	#unset accumulatedDiameter accumulatedDiameter2
	unset averageDiameter stdDevDiameter
}

proc getNanodiscEllipse { atomSelection numFrames outputFile } {
	# Obtains the xy-coordinates of the atomSelection and 
	# transfers this data to a python script which 
	# calculates the ellipse of best fit and returns 
	# the semi-major and semi-minor axes. The semi-major 
	# axis should be a good representation of the "diameter".

	# This is similar to the procedure "getNanodiscDiameter", 
	# except this does not assume the nanodisc is a perfect 
	# circle.
	# "atomSelection" : string input for use in a VMD atomselect
	# command which will be used to define the nanodisc diameter.
	# Example "segname MSP1 MSP2 and backbone"
	# "numFrames" : Integer input which tells the number of 
	# frames to be analyzed. Currently, this starts at frame 0, 
	# but could modify later to start at an arbitrary frame.
	# "outputFile" : Name desired for the output file. A '.txt'
	# will be appended to the end of the 'outputFile'
	
	puts "Generating an ellipse of best fit for the nanodisc"
	
	set outFile [ open "${outputFile}.txt" w+ ]

	alignToFirstFrame $atomSelection $numFrames
	set all [ atomselect top "all" ]
	set nanodisc [ atomselect top $atomSelection ]
	set nanodiscAtoms [ $nanodisc get index ]
	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		puts "Analyzing frame $i"
		$all frame $i
		$nanodisc frame $i
		#set tempFile [ open "nanodiscCoordinates.txt" w+ ]
		foreach index $nanodiscAtoms {
			set nanodiscAtom [ atomselect top "index $index" frame $i ]
			set coordinates [ lindex [ $nanodiscAtom get { x y } ] 0 ]
			lappend xNanodiscCoordinates [ lindex $coordinates 0 ]
			lappend yNanodiscCoordinates [ lindex $coordinates 1 ]
			#puts $tempFile "[ lindex $coordinates 0 ] [ lindex $coordinates 1 ]"
			$nanodiscAtom delete
		}
		unset index coordinates
		#close $tempFile

		puts "Fitting ellipse"
		set axes [ split [ exec python3 /home/marcario/scripts/nanodisc/ellipse.py ${xNanodiscCoordinates} ${yNanodiscCoordinates} ] ]
		#set axes [ split [ exec python3 /home/marcario/scripts/ellipse.py ] ]
		#file delete nanodiscCoordinates.txt
		#unset tempFile

		set axes [ lsort -increasing $axes ]
		set minorAxis [ expr { 2 * [ lindex $axes 0 ] } ]
		set majorAxis [ expr { 2 * [ lindex $axes 1 ] } ]
		#set eccentricity [ expr { sqrt( 1 - ( $minorAxis**2 / $majorAxis**2  )  ) } ]
		set flat [ expr { 1 - ( $minorAxis / $majorAxis ) } ]
		puts $outFile "$i $majorAxis $flat"
		puts "Frame major axis is $majorAxis. Flattening is $flat."
		lappend totalMajorAxis $majorAxis
		#lappend totalEccentricity $eccentricity
		lappend totalFlattening $flat

		unset axes minorAxis majorAxis flat
		unset -nocomplain xNanodiscCoordinates yNanodiscCoordinates
	}
	$all delete
	$nanodisc delete
	unset nanodiscAtoms

	set diameterStatistics [ getAverageStandardDeviation $totalMajorAxis ]
	#set eccentricityStatistics [ getAverageStandardDeviation $totalEccentricity ]
	set flatteningStatistics [ getAverageStandardDeviation $totalFlattening ]
	puts "Average Nanodisc Diameter is [ lindex $diameterStatistics 0 ] +/- [ lindex $diameterStatistics 1 ]."
	puts "Average Nanodisc Flattening is [ lindex $flatteningStatistics 0 ] +/- [ lindex $flatteningStatistics 1 ]."
	puts $outFile "# [ lindex $diameterStatistics 0 ] [ lindex $diameterStatistics 1 ] [lindex $flatteningStatistics 0 ] [ lindex $flatteningStatistics 1 ]"
	#close $outFile
	unset totalMajorAxis totalFlattening
	unset diameterStatistics flatteningStatistics

	close $outFile
		
}		

proc leafletSorter { atomSelection membraneMidplane } {
	# Sorts lipids into upper and lower leaflets
	# and returns a list of lists containing 1. a list of 
	# upper leaflet lipids and 2. a list lower leaflet lipids.

	# "atomSelection" : A (currently) one atom selection to use 
	# as the determinant for being either above or below membrane 
	# midplane. Example "name P" for the lipid's phosphorous atom.
	# "membraneMidplane" : This value will be derived from the function
	# 'getMembraneMidplane' below

	if { ![ $atomSelection num ] } {
		error "Your atom selection contains no atoms with which to sort. Try again."
	}

	#puts "Sorting lipids into upper and lower leaflets"
	set lipidNum [ $atomSelection num ]
	set lipidIndices [ $atomSelection get index ]
	foreach index $lipidIndices {
		set sortingAtom [ atomselect top "index $index" ]
		set atomPosition [ $sortingAtom get {z} ]
		if { $atomPosition > $membraneMidplane } {
			lappend upperLipids $index
			puts "$atomPosition upper"
		} else {
			lappend lowerLipids $index
			puts "$atomPosition lower"
		}
		$sortingAtom delete
		unset atomPosition 
	}
	unset lipidNum lipidIndices membraneMidplane

	#puts "Done sorting lipids. Enjoy Coke!"
	# Handles the cases where there may be upper lipids but no lower lipids
	# and vice versa.
	if { ![ info exists upperLipids ] } { 
		lappend upperLipids -1
	} elseif { ![ info exists lowerLipids ] } {
		lappend lowerLipids -1 
	}
	#puts "Upper lipids are $upperLipids. Lower lipids are $lowerLipids."

	return [ list $upperLipids $lowerLipids ]
}


proc initialize2DArray { numI numJ } {
	# Initialzes a two-dimensional array
	# for use in sotring data with each entry
	# as (0.0, 0.0).

	for { set i 0 } { $i < $numI } { incr i 1 } {
		for { set j 0 } { $j < $numJ } { incr j 1 } {
			set matrix($i,$j) 0.0
		}
		unset j
	}
	unset i
	
	return [ array get matrix ]
}

proc initialize1DArray { numI } { 
	# Intializes a one-dimensional array 
	# for use in sorting data with each
	# entry as 0.0.
	
	for { set i 0 } { $i < $numI } { incr i 1 } { 
		set matrix($i) 0.0
	}
	unset i

	return [ array get matrix ]
}

proc initialize1DArraySpecial { numList } { 
	# Initializes a 1D array with 
	# a specific list of numbers as
	# the keys 
	# "numList" : The list of numbers used 
	# to initialize this array.

	foreach number $numList {
		set matrix($number) 0.0
		unset number
	}
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
	# "frameNum" : The frame number for which you are binning the atoms

	#puts "Obtaining '$atomSelection' atoms $innerRadius to $outerRadius away from protein surface."
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
	#puts "Atoms within distance are $outerSelectionIndices"

	return $outerSelectionIndices
}

proc getMembraneMidplane { atomSelection } {
	# Obtains the z-coordinate for the membrane midplane
	# based on the geometric center of the atoms.
	# Returns a number containing the z-coordinate of 
	# membrane midplane for the given atomselection.
	# "atomSelection" : The atoms to use to determine 
	# the gemometric center of the membrane

	# LOOK! Comments are longer than the code.

	set center [ measure center $atomSelection ]
	return [ lindex $center 2 ]
}

proc getMembraneHeight { atomIndices membraneMidplane leaflet frameNum } { 
	# Calculates the geometric center of a group of atoms
	# and returns the z-coordinate of that center as
	# the height above the midplane.
	
	# "atomIndices" : The indices of atoms for which to measure the height
	# of the membrane above the midplane. These will be provided as input
	# from the main body of "membraneHeightFromProtein".
	# "membraneMidplane" : The midpoint of the membrane, which will be provided
	#  from the result of "getMembraneMidplane" above.
	#  "leaflet" : This requests the height of a specific leaflet above the 
	#  membrane midplane to be measured in case of asymmetries. This will be
	#  provided from the main body of "membraneHeightFromProtein"

	set membraneAtoms [ atomselect top "index $atomIndices" frame $frameNum ] 
	set membraneZ [ lindex [ measure center $membraneAtoms ] 2 ]
	if { [ string match $leaflet upper ] } {
		set membraneHeight [ expr { $membraneZ - $membraneMidplane  } ]
	} elseif { [ string match $leaflet lower ] } { 
		set membraneHeight [ expr { $membraneMidplane - $membraneZ } ]
	}
	$membraneAtoms delete
	unset membraneZ
	return $membraneHeight
}

proc alignToFirstFrame { atomSelection numFrames } { 
	# Aligns each frame to the position 
	# of "atomSelection" in the first frame 
	# of the trajectory.
	# "atomSelection" : String used for VMD atomselection.
	# Example of "protein and backbone".

	set referenceSelection [ atomselect top $atomSelection  frame 0 ]
	set allAtoms [ atomselect top "all" ]
	set trajectorySelection [ atomselect top $atomSelection ]
	for { set i 0 } { $i < $numFrames } { incr i 1 } {
		$allAtoms frame $i
		$trajectorySelection frame $i
		$allAtoms move [ measure fit $trajectorySelection $referenceSelection ]
		$allAtoms moveby [ vecinvert [ measure center $referenceSelection ] ]
	}
	$referenceSelection delete
	$allAtoms delete
	$trajectorySelection delete
	unset i
	puts "Trajectory aligned"
}

proc findZeroesInArray { arrayName  } {
	# Takes an array of "arrayName" as input and 
	# returns a list containing all keys
	# in the array that are zero.
	
	upvar $arrayName arrName
	foreach { key value } [ array get arrName ] {
		if { ![ expr { int($value)  } ] } {
			lappend zeroList $key
			puts "$key has a zero value"
		}
	}

	return $zeroList
}

proc membraneHeightFromProtein { atomSelection proteinSelection nanodiscSelection distanceFromProtein binResolution numFrames } {
	# Calculates the average height of an atomselection 
	# as a function of distance from the protein
	# surface averaged over all theta.
	# "atomSelection" : A string input for which atoms will 
	# define the height of the membrane.
	# "proteinSelection" : A string input for which atoms 
	# will define the surface of the protein.
	# "nanodiscSelection" : "A string input which represents atoms
	# of the nanodisc used for alignment in this calculation.
	# "distanceFromProtein" : A numeric input which represents 
	# the total distance from the protein surface for which 
	# membrane height will be calculated.
	# "binResolution" : A numeric input which sets how wide the bins will be.
	# "numFrames" A numeric which tells how many frames will be analyzed. 
	# For now, the analysis starts at frame 0 and marches forward
	# until "numFrames" is analyzed.
	# Outsput is a text file containing upper height +/- SD, lower height +/- SD,
	# and membrane thickness +/- SD relative to membrane surface

	puts "Calculating membrane thickness. Hold tight."
	set outFile [ open "membraneThickness.txt" w+ ]

	# Converts distanceFromProtein and binResolution to float to make math work out
	if { [ string is int $distanceFromProtein ] } { 
		set distanceFromProtein [ expr { double( $distanceFromProtein )  } ] 
	}
	if { [ string is int $binResolution ] } {
		set binResolution [ expr { double( $binResolution ) } ]
	}

	set numBins [ expr { int( ceil( $distanceFromProtein / $binResolution ) ) } ]
	#puts "Number of bins for data collection is $numBins"
	array set upperHeights [ initialize1DArray $numBins ]
	array set lowerHeights [ initialize1DArray $numBins ]
	array set upperHeights2 [ initialize1DArray $numBins ]
	array set lowerHeights2 [ initialize1DArray $numBins ]
	array set upperFrames [ initialize1DArray $numBins ]
	array set lowerFrames [ initialize1DArray $numBins ]
	
	# Creates a nested list containing the bounding distances for each bin
	for { set i 0 } { $i < $numBins } { incr i 1 } {
		if { !${i} } { 
			set lowerBound 0.0
		}
		set upperBound [ expr { $lowerBound + $binResolution } ]
		if { $upperBound > $distanceFromProtein } {
			set upperBound $distanceFromProtein
		}
		lappend bins "$lowerBound $upperBound"
		set lowerBound $upperBound
		unset upperBound
	} 
	unset lowerBound
	#puts "$bins"

	# Does not do alignment unless nanodiscSelection is set
	if { [ string match -nocase "none" $nanodiscSelection ] == 0 } {
		alignToFirstFrame $nanodiscSelection $numFrames
	}

	# Heavy lifting. This is the section that measures leaflet heights for each frame 
	# as a function of distance from the protein. 
	set membraneSelection [ atomselect top "$atomSelection" ]
	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		puts -nonewline "Analyzing frame $i..."
		$membraneSelection frame $i
		set membraneMidplane [ getMembraneMidplane $membraneSelection ]
		puts -nonewline [ format "Membrane midplane is %.2f " $membraneMidplane ]
		set j 0;  # This is a counter to place heights in right location of data array.
		foreach bin $bins {
			set indicesWithinDistance [ getAtomsDistanceFromProtein $atomSelection $proteinSelection [ lindex $bin 0 ] [ lindex $bin 1 ] $i ]
			if { ![ llength $indicesWithinDistance ] } { 
				puts "No lipids within this distance. Heading to next bin."
				incr j 1
				continue
			}
			set lipidsWithinDistance [ atomselect top "index $indicesWithinDistance" frame $i ]
			set sortedLipids [ leafletSorter $lipidsWithinDistance $membraneMidplane ]
			set upperLipidIndices [ lindex $sortedLipids 0 ]
			set lowerLipidIndices [ lindex $sortedLipids 1 ]
			# The following "if" loops only evaluate leaflet height if there are atoms 
			# returned from the leaflet sorter. If not, if will not evaluate leaflet 
			# height. The "-1" comes from the leafletSorter procedure. If there are no
			# atoms in the leaflet, the leafletSorter will return a list containing only
			# one element, -1. If no lipids are in a leaflet, a counter is increased 
			# for calculation of averages later.
			puts -nonewline "For bin [ lindex $bin 0 ] - [ lindex $bin 1 ], " 
			if { [ lsearch $upperLipidIndices -1 ] == -1 } {
				set upperHeight [ getMembraneHeight $upperLipidIndices $membraneMidplane "upper" $i ]
				set upperHeights($j) [ expr { $upperHeights($j) + $upperHeight  } ]
				set upperHeights2($j) [ expr { $upperHeights2($j) + $upperHeight**2  } ]
				set upperFrames($j) [ expr { $upperFrames($j) + 1 } ]
				puts -nonewline "there are [ llength $upperLipidIndices ] upper lipids and "
			} else {
				puts -nonewline "there are 0 upper lipids and "
			}
				
			if { [ lsearch $lowerLipidIndices -1 ] == -1 } { 
				set lowerHeight [ getMembraneHeight $lowerLipidIndices $membraneMidplane "lower" $i ]
				set lowerHeights($j) [ expr { $lowerHeights($j) + $lowerHeight  } ]
				set lowerHeights2($j) [ expr { $lowerHeights2($j) + $lowerHeight**2 } ]
				set lowerFrames($j) [ expr { $lowerFrames($j) + 1 } ]
				puts " [ llength $lowerLipidIndices ] lower lipids."
			} else {
				puts "0 lower lipids."
			}
			$lipidsWithinDistance delete
			unset indicesWithinDistance sortedLipids
			unset upperLipidIndices lowerLipidIndices 
			unset -nocomplain upperHeight lowerHeight
			incr j 1
		}
		unset membraneMidplane
		unset j bin
		puts "Finished analyzing frame $i"
	}
	unset i

	#set upperHeightZeroes [ findZeroesInArray upperHeights ]
	#set lowerHeightZeroes [ findZeroesInArray lowerHeights ]
	#set upperFirstZero [ lindex [ lsort -increasing -integer $upperHeightZeroes ] 0 ]
	#if { $upperFirstZero == 0 } { 
	#	set upperFirstZero [ lindex [ lsort -increasing -integer $upperHeightZeroes ] 1 ]
	#}
	#set lowerFirstZero [ lindex [ lsort -increasing -integer $lowerHeightZeroes ] 0 ]
	#if { $lowerFirstZero == 0 } {
	#	set lowerFirstZero [ lindex [ lsort -increasing -integer $lowerHeightZeroes ] 1 ]
	#}
	#set firstZeroDifference [ expr { $lowerFirstZero - $upperFirstZero } ]
	#puts "First time to zero is $upperFirstZero, $lowerFirstZero"
	#unset upperHeightZeroes lowerHeightZeroes
	#unset upperFirstZero lowerFirstZero

	# Calculate the average height and standard deviation for each bin.
	# This is messy. May make a procedure later to do this and clean up.
	for { set i 0 } { $i < [ array size upperHeights ] } { incr i 1 } {
		set upperDataFraction [ expr { $upperFrames($i) / $numFrames } ]
		set lowerDataFraction [ expr { $lowerFrames($i) / $numFrames } ]
		puts "Percentage of frames with data for bin $i is $upperDataFraction, $lowerDataFraction"
		#if { $upperFrames($i) } {}
		if { $upperDataFraction > 0.1 } {
			set upperBinAverage [ expr { $upperHeights($i) / $upperFrames($i) } ]
			set upperBinAverage2 [ expr { $upperHeights2($i) / $upperFrames($i) } ]
			#puts "Average of squares $upperBinAverage2, square of average [ expr { $upperBinAverage**2  } ]"		
			set upperStandardDeviation [ expr { sqrt( $upperBinAverage2 - $upperBinAverage**2  )  } ]
		} else {
			set upperBinAverage ""
			set upperStandardDeviation ""
		}

		#if { $lowerFrames($i) } {}
		if { $lowerDataFraction > 0.1 } {
			set lowerBinAverage [ expr { $lowerHeights($i) / $lowerFrames($i) } ]
			set lowerBinAverage2 [ expr { $lowerHeights2($i) / $lowerFrames($i) } ]
			#puts "Average of squares $lowerBinAverage2, square of average [ expr { $lowerBinAverage**2  } ]"
			set lowerStandardDeviation [ expr { sqrt( $lowerBinAverage2 - $lowerBinAverage**2  )  } ]
		} else {
			set lowerBinAverage ""
			set lowerStandardDeviation ""
		}
		if { [ llength $upperBinAverage ] == 0 || [ llength $lowerBinAverage ] == 0 } {
			set thicknessBinAverage ""
			set thicknessBinStandardDeviation ""
		} else {
			#set j [ expr { $i + $firstZeroDifference  } ]
			#set lowerBinForHeight [ expr { $lowerHeights($j) / $lowerFrames($j) } ]
			#set lowerBinForHeight2 [ expr { $lowerHeights2($j) / $lowerFrames($j)  } ]
			#set lowerStandardDeviationForHeight [ expr { sqrt( $lowerBinAverage2 - $lowerBinAverage**2  )  } ]
			set thicknessBinAverage [ expr { $upperBinAverage + $lowerBinAverage } ]
			set thicknessBinStandardDeviation [ expr { sqrt( $upperStandardDeviation**2 + $lowerStandardDeviation**2  )  } ]
		}
		puts $outFile " [ lindex [ lindex $bins $i ] 0 ] $upperBinAverage $upperStandardDeviation $lowerBinAverage $lowerStandardDeviation $thicknessBinAverage $thicknessBinStandardDeviation"
		puts "For bin $i, upper height is $upperBinAverage +/- $upperStandardDeviation, lower height \
		is $lowerBinAverage +/- $lowerStandardDeviation, and average thickness is $thicknessBinAverage \
		+/- $thicknessBinStandardDeviation"
	}
	unset upperHeights lowerHeights upperHeights2 lowerHeights2
	unset upperFrames lowerFrames 
	unset upperBinAverage lowerBinAverage thicknessBinAverage thicknessBinStandardDeviation
	unset -nocomplain upperBinAverage2 lowerBinAverage2 lowerBinForHeight lowerBinForHeight2 lowerStandardDeviationForHeight
	close $outFile

	puts "Calculation of membrane THICC-ness complete."
	puts "Hold on to those curves!"
}

proc proteinNanodiscInteractions { proteinSelection nanodiscSelection interactionDistance inputPDB } { 
	# This function calculates the percentage of frames that
	# the embedded protein is in contact with the nanodisc. It 
	# outputs both a plain text file containing a per residue
	# contact percentage as well as PDB demonstrating contact 
	# percentage on a molecular model.
	
	# "proteinSelection" : Which parts of the protein to consider
	# when evaluating for contact. This should be in the format of
	# a string to be used in a VMD atomselection procedure. Example
	# of "(resid 200 to 322)".
	# "nanodiscSelection" : Which parts of the nanodisc to consider
	# when evaluating for contact. This should be in the format of 
	# a string to be used in a VMD atomselection procedure. Example
	# of "segname MSP1 MSP2"
	# "interactionDistance" : The distance below which you would like
	# to consider a residue on the embedded protein to be in contact 
	# with the nanodisc. 
	# "inputPDB" : The path to the PDB file you would like to map the 
	# contact percentages onto. This is usually the strating structure.

	set numFrames [ molinfo top get numframes ]
	
	# Gets the total number of residues being analyzed in order
	# to create the array to store data
	set protein [ atomselect top $proteinSelection ]
	set proteinResidues [ $protein get resid ]
	for { set i 0 } { $i < [ llength $proteinResidues ] } { incr i 1 } { 
		if { !$i } {
			lappend proteinResidList [ lindex $proteinResidues $i ]
			continue
		}
		set resid [ lindex $proteinResidues $i ]
		if { [ lsearch $proteinResidList $resid ] == -1 } {
			lappend proteinResidList $resid
			unset resid
		}
	}
	unset i proteinResidues
	$protein delete
	
	# Sets the array for monitoring interactions and modifies it so the key
	# is the resid number and not a serial numbering
	array set interactionArray [ initialize1DArraySpecial $proteinResidList ]
	#for { set i 0 } { $i < [ llength $proteinResidList ] } { incr i 1 } { 
	#	set resid [ lindex $proteinResidList $i ]
	#	set interactionArray($resid) $interactionArray($i)
	#	unset interactionArray($i)
	#	unset i resid
	#}

	set overlappingResidues [ atomselect top "$proteinSelection and same residue as (within $interactionDistance of ($nanodiscSelection))" ]
	
	for { set i 0 } { $i < $numFrames } { incr i 1 } {
		puts "Analyzing frame $i"
		$overlappingResidues frame $i
		$overlappingResidues update
		set overlappingResids [ $overlappingResidues get resid ]
		if { ![ llength $overlappingResids ] } {
			puts "No overlapping residues for frame $i"
			continue
		}
		set overlappingResidList [ lsort -unique $overlappingResids ]
		puts "Overlapping residues are $overlappingResidList"
		foreach resid $overlappingResidList {
			set interactionArray($resid) [ expr { $interactionArray($resid) + 1.0  } ]
			unset resid
		}
		unset -nocomplain overlappingResids overlappingResidList
	}
	$overlappingResidues delete

	puts "Creating PDB of results"
	set inputMolNum [ mol load pdb $inputPDB ]
	set all [ atomselect $inputMolNum "all" ]
	$all set occupancy "0.0"
	foreach resid $proteinResidList {
		set interactionArray($resid) [ expr { $interactionArray($resid) / $numFrames  } ]
		puts "Residue $resid is in direct contact with the nanodisc scaffold $interactionArray($resid) of the time"
		set residue [ atomselect $inputMolNum "protein and resid $resid" ]
		$residue set occupancy "$interactionArray($resid)"
		$residue delete
	}
	$all writepdb proteinNanodiscContacts.pdb
	$all delete
	mol delete $inputMolNum

}

proc measureDistance { atomGroup1 atomGroup2 proteinSegments numFrames outputFile  } {
	# Measures the distance between two atom groups
	# and outputs to a txt file.

	set numberSubunits [ llength $proteinSegments ]
	array set distanceArray [ initialize1DArray $numberSubunits ]

	set outfile [ open "$outputFile.txt" w+ ]

	for { set i 0 } { $i < $numFrames } { incr i 1 } {
		set j 0	
		foreach segment $proteinSegments {
			set selection1 [ atomselect top "segname $segment and $atomGroup1" frame $i ]
			set selection2 [ atomselect top "segname $segment and $atomGroup2" frame $i ]
			set distance [ measure bond "[ $selection1 get index ] [ $selection2 get index ]" frame $i ]
			set distanceArray($j) $distance
			incr j 1
			$selection1 delete
			$selection2 delete
			unset distance segment
		}
		set outputList "$i"
		set totalDistancePerFrame "0.0"
		for { set k 0 } { $k < $numberSubunits } { incr k 1 } { 
			lappend outputList $distanceArray($k)
			set totalDistancePerFrame [ expr { $totalDistancePerFrame + $distanceArray($k)  } ]
		}
		set averageDistancePerFrame [ expr { $totalDistancePerFrame / $numberSubunits } ] 
		lappend averageDistances $averageDistancePerFrame
		puts "Average distance for frame $i is $averageDistancePerFrame"
		unset k
		puts $outfile "$outputList $averageDistancePerFrame"
		unset outputList averageDistancePerFrame totalDistancePerFrame 
	}
	set statistics [ getAverageStandardDeviation $averageDistances ]
	puts $outfile "# [ lindex $statistics 0 ] [ lindex $statistics 1 ]"
	puts "Average is [ lindex $statistics 0 ] +/- [ lindex $statistics 1 ]"
	unset i
	close $outfile
		
}

proc measureRMSF { atomSelection rmsfAtom numFrames outputFile inputPDB  } {
	# Measures the RMSF of the atomselection
	# and outputs the result to a txt file and 
	# onto a pdb backbone.
	#
	# "atomSelection" : Which parts of the protein to measure the RMSF
	#  for. Example "protein" or "(resid 10 to 35)".
	#  "rmsfAtom" : The atom you will be using to measure the RMSF. 
	#  This is usually something like the alpha carbon or some atom
	#  on the backbone.
	#  "numFrames" : The total number of frames to be analyzed. By default, 
	#  this starts at the initial frame.
	#  "outputFile" : The desired name for the output file. An extension of
	#  '.txt' will be appended to "outputName"
	# "inputPDB" : The path to the PDB file you would like to map the 
	# contact percentages onto. This is usually the strating structure.
	
	set outfile [ open "${outputFile}.txt" w+ ]

	alignToFirstFrame $atomSelection $numFrames
	set selection [ atomselect top $atomSelection ]
	set selectionResidues [ $selection get resid ]
	set selectionResidues [ lsort -increasing -unique -integer $selectionResidues ]
	array set rmsfValues [ initialize1DArraySpecial $selectionResidues ]
	$selection delete

	foreach residue $selectionResidues {
		puts -nonewline "Analyzing RMSF of residue $residue."
		set rmsfSelection [ atomselect top "$atomSelection and resid $residue and $rmsfAtom" ]
		set rmsf [ measure rmsf $rmsfSelection ]
		puts " RMSF is $rmsf."
		puts $outfile "$residue $rmsf"
		set rmsfValues($residue) $rmsf
		$rmsfSelection delete
		unset rmsf residue
	}
	close $outfile

	puts "Creating PDB of RMSF results"
	set molNum [ mol load pdb ${inputPDB} ]
	set all [ atomselect $molNum "all" ]
	$all set occupancy 0.0
	foreach residue $selectionResidues {
		set selection [ atomselect $molNum "resid $residue" ]
		$selection set occupancy $rmsfValues($residue)
		$selection delete
	}
	$all writepdb ${outputFile}.pdb
	$all delete
	mol delete $molNum

}

proc getAverageRMSF { rmsfFileList outputFile inputPDB } {
	# Calculates the average RMSF for each residue
	# given the file prefixes. To be used if the 
	# RMSF was calculated for each subunit individually.
	# "rmsfFileList" : The list of files for which the 
	# average RMSF will be calculated.
	# "inputPDB" : The path to the PDB file you would like to map the 
	# contact percentages onto. This is usually the starting structure.

	set outfile [ open "${outputFile}.txt" w+ ]

	set numEntries [ llength $rmsfFileList ]
	puts "Number of files to read is $numEntries"

	set fileData [ open [ lindex $rmsfFileList 0 ] ]
	while { [ gets ${fileData} data ] >= 0 } {
		set residue [ lindex $data 0 ]
		lappend residueList $residue
	} 
	close $fileData
	array set totalRMSF [ initialize1DArraySpecial $residueList ]

	for { set i 0 } { $i < $numEntries } { incr i 1 } {
		set fileData [ open [ lindex $rmsfFileList $i ] r ]
		while { [ gets ${fileData} data ] >= 0 } {
			set residue [ lindex $data 0 ]
			set rmsf [ lindex $data 1 ]
			set totalRMSF($residue) [ expr { $rmsf + $totalRMSF($residue) } ]
			puts "Total RMSF is $i $totalRMSF($residue)"
			unset data residue rmsf
		}
		close $fileData
	}
	unset i

	# Sorts the array for printing to file
	foreach { residue rmsf } [ array get totalRMSF ] {
		lappend totalRMSFList [ list $residue $rmsf ]
	}
	puts [ lsort -increasing -index 0 -integer $totalRMSFList ]
	foreach element [ lsort -increasing -index 0 -integer $totalRMSFList ] {
		set avgRMSF [ expr { [ lindex $element 1 ] / $numEntries  } ]
		puts "Average RMSF for residue [ lindex $element 0 ] is [ lindex $element 1 ] $avgRMSF"
		puts $outfile "[ lindex $element 0 ] $avgRMSF"
		unset avgRMSF element
	}

	close $outfile

}	

proc writeAlignedTrajectory { psfFile pdbFile listOfDCDFiles alignmentAtomSelection stride outputAtomSelection outputDCD } {
	# Takes dcd inputs and align the nanodisc
	# to a conformation where nanodisc is in the xy-plane 
	# so that it can be used to calculate lateral pressure
	# profiles post-hoc.
	#
	# THIS PROCEDURE HAS NOT BEEN TESTED AND CANNOT 
	# RECOMMEND ITS USE.

	# Load the files to be manipulated
	mol load psf $psfFile pdb $pdbFile
	foreach dcdFile $listOfDCDFiles {
		mol addfile $dcdFile type dcd step $stride waitfor all
	}
	
	set numFrames [ molinfo top get numframes ]
	alignToFirstFrame $alignmentAtomSelection $numFrames

	puts "Writing aligned trajectory to file"
	set atomsToBeWritten [ atomselect top $outputAtomSelection ]
	animate write dcd $outputDCD beg 0 end [expr { $numFrames - 1} ] skip 0 waitfor all sel $atomsToBeWritten
	$atomsToBeWritten delete
}

proc helixMotion { proteinSelection segments headgroupSelection nanodiscSelection numFrames outputFile } {
	# Calculates the spherical coordinates of the helix provided 
	# to later be plotted with the plotting software of choice.
	# This procedure calculates the inertial tensor for the helix
	# and uses the principal components to measure its spherical 
	# position. This assumes a rigid helix.
	# "proteinSelection" : A string to be used as a VMD atomselection
	# procedure to describe the helix to be measured. 
	# "segments" : The segnames of the protein which should be measured
	# in this procedure. 
	# "headgroupSelection" : Selection of the lipid headgroup which will help
	# determine the position of the helix relative to the membrane.
	# "nanodiscSelection" : A string to be used as a VMD atomselection
	# procedure to select the relevant parts of the nanodisc. 
	# "numFrames" : Total number of frames to be analyzed. By default, this
	# analysis starts at the initial frame.
	# "outputFile" : The name desired for the data output. This will be in 
	# plain text format. A '.txt' extension will be added to "outputFile".
	
	global M_PI

	#set all [ atomselect top "all" frame 0 ]
	#set pore [ atomselect top "(resid 227 to 248) and backbone" frame 0 ]
	#$all moveby [ vecinvert [ measure center $pore ] ]
	#$all delete
	#$pore delete

	alignToFirstFrame $nanodiscSelection $numFrames

	set outFile [ open "${outputFile}.txt" w+ ]

	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		puts "Measuring data for frame $i"
		set headgroup [ atomselect top $headgroupSelection frame $i ]
		set headgroupHeight [ lindex [ measure center $headgroup ] 2 ]
		# Calculates the spherical coordinates for each segment to be averaged later
		foreach segment $segments {
			#puts "Finding best atoms for segment $segment"
			set helix [ atomselect top "segname $segment and $proteinSelection" frame $i ]
			# Deriving the radial and theta coordinates
			set coordinates [ $helix get { x y z } ]
			# This function finds the atoms that are closest to the headgroup 
			# as defined in "headgroupSelection"
			for { set j 0 } { $j < [ llength $coordinates ] } { incr j 1 } { 
				set distanceFromHeadgroup [ expr { abs( $headgroupHeight - [ lindex [ lindex $coordinates $j ] 2 ]  ) } ]
				if { $j < 5 } {
					lappend zCoor [ list $j $distanceFromHeadgroup ]
					continue
				}
				lappend zCoor [ list $j $distanceFromHeadgroup ]
				#puts $zCoor
				set zCoor [ lsort -real -increasing -index 1 $zCoor ]
				#puts $zCoor
				set zCoor [ lreplace $zCoor end end ]
				#puts $zCoor
			}
			#puts $zCoor
			for { set k 0 } { $k < [ llength $zCoor ] } { incr k 1 } { ;# Extracts indices of atoms closest to headgroup plane
				lappend planeIndices [ lindex [ lindex $zCoor $k ] 0 ]
			}
			unset j k
			set totalX 0.0
			set totalY 0.0
			foreach index $planeIndices {
				set planeXYZ [ lindex $coordinates $index ]
				set totalX [ expr { $totalX + [ lindex $planeXYZ 0 ] } ]
				set totalY [ expr { $totalY + [ lindex $planeXYZ 1 ] } ]
			}
			set avgX [ expr { $totalX / 4 } ]
			set avgY [ expr { $totalY / 4 } ]
			set polarCoordinates [ cartToPolar2D $avgX $avgY ]
			unset planeIndices planeXYZ avgX avgY totalX totalY
			unset zCoor distanceFromHeadgroup coordinates
	
			# Sets the initial polar coordinate for each subunit to make relative polar changes evident
			if { $i == 0 } {
				lappend polarInitial [ lindex $polarCoordinates 1 ]
			}
			#puts "Starting polar angles are $polarInitial"
			# Obtaining the azimuthal angle
			set inertialCoordinates [ measure inertia $helix ]
			set principalAxes [ lindex $inertialCoordinates 1 ]
			set principalAxes [ lsort -increasing -index 2 $principalAxes ] ;# added after initial data analysis
			set helicalAxis [ lindex $principalAxes end ]
			set helicalAngle [ expr { acos( [ vecdot { 0 0 1 } $helicalAxis ] ) }  ]
			puts "$helicalAngle"
			lappend frameData "[ lindex $polarCoordinates 0 ] [lindex $polarCoordinates 1 ] $helicalAngle"
			unset inertialCoordinates principalAxes helicalAxis helicalAngle
			$helix delete
		}
		# Calculates the average for the frame
		set radialTotal 0.0
		set polarTotal 0.0
		set azimuthalTotal 0.0
		for { set l 0 } { $l < [ llength $frameData ] } { incr l 1 } {
			#puts "[lindex $frameData $l]"
			set radialTotal [ expr { $radialTotal + [ lindex [ lindex $frameData $l ] 0 ] } ]
			set polarTotal [ expr { $polarTotal + ( [ lindex [ lindex $frameData $l ] 1 ] - [ lindex $polarInitial $l ] ) } ]; #measures polar angle relative to initial
			set azimuthalTotal [ expr { $azimuthalTotal + [ lindex [ lindex $frameData $l ] 2 ] } ]
		}
		unset l 
		set radialAverage [ expr { $radialTotal / 5 } ]
		set polarAverage [ expr { $polarTotal / 5 } ]
		set azimuthalAverage [ expr { ( 180 / $M_PI ) * ( $azimuthalTotal / 5 ) } ]; # [ expr { ( 180 / $M_PI ) * ( $azimuthalTotal / 5 ) } ] 
		#puts "$radialAverage $polarAverage $azimuthalAverage"

		puts $outFile "$i [ lindex $frameData 0 ] [ lindex $frameData 1 ] [ lindex $frameData 2 ] [ lindex $frameData 3 ] [ lindex $frameData 4 ] $radialAverage $polarAverage $azimuthalAverage"
		unset radialAverage polarAverage azimuthalAverage
		unset radialTotal polarTotal azimuthalTotal
		unset frameData
		$headgroup delete
	}

	close $outFile
}

proc distanceFromGeometricCenterMeasurement { atomSelection numFrames } { 
	# Measures the distance of a group of atoms from their geometric center
	# atomSelection: VMD-style atom selection for which atoms should be evaluated
	# numFrames: the number of frames to be analyzed (should include start and end later)

	set totalDistance 0.0
	set totalDistance2 0.0	
	set initialAtoms [ atomselect top $atomSelection frame 0 ]
	set atomsNum [ $initialAtoms num ]
	set geometricCenter [ measure center $initialAtoms ]
	set geometricCenterXY [ lrange $geometricCenter 0 1 ]
	$initialAtoms delete
	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		set atoms [ atomselect top $atomSelection  frame $i ]
		set atomPositions [ $atoms get { x y z } ]
		set frmTotalDistance 0.0
		$atoms delete
		foreach atomPosition $atomPositions {
			set distance [ vecdist $geometricCenterXY [ lrange $atomPosition 0 1 ] ]
			set frmTotalDistance [ expr { $frmTotalDistance + $distance } ]
			unset distance
		}
		unset atomPositions atoms
		set frmAvg [ expr { $frmTotalDistance / $atomsNum } ]
		set totalDistance [ expr { $totalDistance + $frmAvg } ]
		set totalDistance2 [ expr { $totalDistance2 + $frmAvg**2 } ]
		unset frmTotalDistance frmAvg
	}
	unset geometricCenter geometricCenterXY atomsNum

	set averageDistance [ expr { $totalDistance / $numFrames } ]
	set averageDistance2 [ expr { $totalDistance2 / $numFrames } ]
	set stdDevDistance [ expr { sqrt( $averageDistance2 - $averageDistance**2 )  } ]
	puts "Distance is $averageDistance +/- $stdDevDistance"
	
	return [ list $averageDistance $stdDevDistance ]		
}

proc atomToPoreMeasurement { fragment atomName resids fittingGroup numFrames outputFile  } {
	# Measures the average distance of atom types
	# in list to the geometric center of the pore.
	# "fragment": The protein selection for which the measurement
	# applies, such as "segname PROA".
	# "atomName": the atom type name for the atoms being measured (i.e., name CA)
	# "resids": the residue numbers to be measured given as a range (i.e.,
	# "1 20"
	# "fittingGroup": the group of atoms used to fit the structure to its initial
	# configuration. This should be given as a string to be used in a VMD
	# atomselection procedure.
	# "numFrames": The total number of frames to analyze. By deafult, this analysis
	# starts at the initial frame.
	# "outputFile": The desired name of the outputfile. This will be a plain text file
	# with the format "Residue Distance". A '.txt' extnesion will be appended to "outputFile".

	alignToFirstFrame $fittingGroup $numFrames
	#puts "Trajectory aligned in procedure"
	set outFile [ open "${outputFile}.txt" w+ ]
	# Makes list of residues to be measured
	if { [ llength $resids ] != 2 } {
		error "Length of residue list is too small or too long. Try again please."
	}
	for { set i [ lindex $resids 0 ] } { $i <= [ lindex $resids end ] } { incr i 1 } {
		lappend residList $i
	}
	puts "Residues being considered are: $residList"
	unset i

	foreach resid $residList {
		puts -nonewline "Residue $resid: "
		set distanceData [ distanceFromGeometricCenterMeasurement "$fragment and $atomName and resid $resid" $numFrames ]
		puts $outFile "$resid [ lindex $distanceData 0 ] [ lindex $distanceData 1 ]"
		unset resid distanceData
	}
	puts "Finished analyzing distance from center!"
	close $outFile
}

proc extractHighContactResiduesFromPDB { pdbFile keyAtom threshold  } {
	# Takes the pdb output from the nanodisc interaction
	# procedure and makes a list of all residues with at least
	# THRESHOLD amount of contact between protein and 
	# nanodisc. Returns that list for further analysis.
	# "pdbFile": The .pdb file from "proteinNanodiscInteractions" 
	# which contains the average contact information from the simulation
	# "keyAtom": Atom type used to unqiuely define each residue such that the
	# list does not contain multiple entries from the same residue (i.e., CA
	# or N or C)
	# "threshold": The threshold percentage for which you want to measure
	# the contact types

	set molNum [ mol load pdb $pdbFile ]
	foreach letter "A B C D E" {
		[ atomselect top "segname ELC${letter}" ] set segname PRO${letter}
	}
	set highContactResidues [ atomselect top "$keyAtom and occupancy > $threshold" ]
	set highContactSeg [ $highContactResidues get segid ]
	set highContactRes [ $highContactResidues get resid ]
	$highContactResidues delete
	puts $highContactSeg
	puts $highContactRes
	mol delete $molNum

	return [ list $highContactSeg $highContactRes ] 
}

proc getUniqueContacts { segidList residList } {
	# Takes as input a list of segids and resids
	# and puts out a list of unique residues. This 
	# handles issues where the resid on MSP1 and MSP2
	# is the same.
	
	set uniqueList ""
	foreach segid $segidList resid $residList {
		puts "$segid $resid"
		if { [ lsearch $uniqueList [list $segid $resid] ] == -1 } {
			lappend uniqueList [ list $segid $resid ]
		}
	}
		puts $uniqueList
	return $uniqueList
}

proc assignContactResidueType { contactList } {
	# Takes in a contact list and outputs
	# the number of each type of interacting
	# residue in the list (i.e., apolar, polar, 
	# basic, acidic) and returns the number of 
	# in a list.
	
	foreach item $contactList {
		set contactSegid [ lindex $item 0 ]
		set contactResid [ lindex $item 1 ]
		set contactResidue [ atomselect top "segname $contactSegid and resid $contactResid" ]
		set contactResidueName [ $contactResidue get resname ]
		set residueType [ color restype [ lindex $contactResidueName 0 ] ]
		lappend residueTypeList $residueType
		unset residueType contactResidueName contactResid contactSegid
		$contactResidue delete
	}
	return $residueTypeList
}

proc measureContactTypes { proteinSegid proteinResid interactionDistance nanodiscSegments numFrames outFile } { 
	# This procedure takes as input the protein segid and resid as well as the interaction
	# distance and the nanodisc segments and outputs the relative abundance 
	# of the four main contact types "Polar, Non-polar, acidic, basic".
	#
	# "proteinSegid" : The segment ID for the residue being measured. This will be
	# provided by the main function "calculateInteractionTypes".
	# "proteinResid" : The residue ID for the residue being measured. This will be
	# provided by the main function "calculateInteractionTypes".
	# "interactionDistance" : The distance, below which, a residue from the protein
	# is said to be in contact with the nanodisc. This will be
	# provided by the main function "calculateInteractionTypes".
	# "nanodiscSegments" : The segments of the nanodisc for which to analyze potential
	# contacts. This will be provided by the main function "calculateInteractionTypes".
	# "numFrames" : The total number of frames to be analyzed by this procedure. By
	# default, this starts at the initial frame. This will be 
	# provided by the main function "calculateInteractionTypes".
	# "outFile" : The desired name for the output file. This will be
	# provided by the main function "calculateInteractionTypes".

	set numNanodiscSeg [ llength $nanodiscSegments ] ;# think about incorporating a per nanodisc interaction
	set totalNumContacts 0
	lappend outData $proteinSegid $proteinResid
	array set contactArray [ initialize1DArraySpecial [ list Nonpolar Polar Basic Acidic ] ]
	set contacts [ atomselect top "segname $nanodiscSegments and same residue as (within $interactionDistance of (segname $proteinSegid and resid $proteinResid))" ]
	for { set i 0 } { $i < $numFrames } { incr i 1 } { 
		$contacts frame $i
		$contacts update
		if { ![ $contacts num ] } { 
			#puts "No contacts for frame $i"
			continue
		}
		set contactSeg [ $contacts get segid ]
		set contactRes [ $contacts get resid ]
		set uniqueContacts [ getUniqueContacts $contactSeg $contactRes ]
		set totalNumContacts [ expr { $totalNumContacts + [ llength $uniqueContacts ] } ]
		set contactResidueList [ assignContactResidueType $uniqueContacts ]
		foreach residueType $contactResidueList {
			set contactArray($residueType) [ expr { $contactArray($residueType) + 1 } ]	
		}
		unset contactSeg contactRes uniqueContacts
		unset contactResidueList
	}
	$contacts delete
	if { $totalNumContacts != 0 } { 
		foreach type "Nonpolar Polar Basic Acidic" {
			#set interactionPercentage [ expr { double( $contactArray($type) / $totalNumContacts ) } ]
			#lappend outData $interactionPercentage
			#unset interactionPercentage
			lappend outData $contactArray($type)
		}
		lappend outData $totalNumContacts
		puts $outFile $outData
		unset outData type totalNumContacts contactArray
	}
		
}

proc calculateInteractionTypes { interactionDistance numFrames contactPercentageThreshold keyAtom pdbFile outFile } {
	# Calculates the interaction types normalized to 
	# the total number of contacts for each residue
	#
	# "interactionDistance" : The distance, below which, a residue from the protein
	# is said to be in contact with the nanodisc.
        # "numFrames" : The total number of frames to be analyzed by this procedure. By
        # default, this starts at the initial frame.
        # "contactPercentageThreshold" : The threshold, above which, the contacting residue
        # types will be analyzed. This was set to 25% before, but is now customizable.
        # "keyAtom" : Arbitrary atom used to search the contact percentage PDB to obtain the 
        # contact percentage for each residue. "name CA" works well as each residue has an alpha-
        # carbon.
        # "pdbFile" : The PDB file containing the contact percentages generated from "proteinNanodiscInteraction"
        # function.
        # "outFile" : Name of the desired output file which will be in plain text format. A '.txt' extension will
        # be appended to "outFile".
	set outFile [ open "${outFile}.txt" w+ ]
	puts $outFile "# Segment Residue Nonpolar Polar Basic Acidic Total Number"
	
	set highContactInfo [ extractHighContactResiduesFromPDB $pdbFile $keyAtom $contactPercentageThreshold ]
	set highContactSeg [ lindex $highContactInfo 0 ]
	set highContactRes [ lindex $highContactInfo 1 ]

	foreach segid $highContactSeg resid $highContactRes {
		puts "Analyzing segname $segid and resid $resid"
		measureContactTypes $segid $resid $interactionDistance "MSP1 MSP2" $numFrames $outFile
	}

	close $outFile
}
