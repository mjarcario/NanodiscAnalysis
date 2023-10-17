#####################################
# Converts 2D Cartesian Coordinates #
# to 2D Polar Coordinates           #
#####################################

proc cartToPolar2D { xCoor yCoor } {
        global M_PI
        set rCoor [ expr { sqrt( ($xCoor * $xCoor) + ($yCoor * $yCoor) ) } ]
        set thetaCoor [ expr { atan($yCoor / $xCoor) } ]
	# Handles the exceptions where atan does not exist and 
	# allows angles up to 2pi 
       if { $xCoor < 0 && $yCoor > 0 } {
                set thetaCoor [ expr { $thetaCoor + $M_PI } ]
        } elseif { $xCoor < 0 && $yCoor < 0 } {
                set thetaCoor [ expr { $thetaCoor + $M_PI } ]
        } elseif { $xCoor > 0 && $yCoor < 0 } {
                set thetaCoor [ expr { $thetaCoor + ( 2 * $M_PI ) } ]
        }
        return [ list $rCoor $thetaCoor ]

}

#################################
# Obtains the mol fraction of a #
# target lipid in a membrane    #
#################################

proc molFraction  { targetLipid bulkLipids } {
	set targetLipids [ atomselect top "resname $targetLipid" ]
	set targetLipidAtomNames [ $targetLipids get name ]
	set targetLipidAtomName [ lindex $targetLipidAtomNames 0 ]
	$targetLipids delete; set targetLipids [ atomselect top "resname $targetLipid and name $targetLipidAtomName" ]
	set numTargetLipids [ $targetLipids num ]
	$targetLipids delete
	unset targetLipidAtomNames; unset targetLipidAtomName

	# Calculates the number of bulk lipids for arbitrary membrane
	# constructs as long as "bulkLipids" is provided as a list
	set numBulk 0.0
	set numBulkLipidSpecies [ llength $bulkLipids ]
	puts -nonewline "Bulk lipid species are "
	for { set i 0 } { $i < $numBulkLipidSpecies } { incr i 1 } {
		puts -nonewline "[ lindex $bulkLipids $i ] "
		set bulkLipidSpecies [ atomselect top "resname [ lindex $bulkLipids $i ] " ] 	
		set bulkLipidAtomNames [ $bulkLipidSpecies get name ]
		set bulkLipidAtomName [ lindex $bulkLipidAtomNames 0 ]
		$bulkLipidSpecies delete; set bulkLipidSpecies [ atomselect top "resname [ lindex $bulkLipids $i ] and name $bulkLipidAtomName" ]
		set numBulkLipids [ $bulkLipidSpecies num ]
		set numBulk [ expr { $numBulk + $numBulkLipids } ]
		$bulkLipidSpecies delete
		unset numBulkLipids
		unset bulkLipidAtomNames; unset bulkLipidAtomName
	}
	set molFraction [ expr { $numTargetLipids / ( $numTargetLipids + $numBulk ) } ]
	puts "\n Number of target lipids is: $numTargetLipids \n Number of bulk lipids is: $numBulk"
	unset numTargetLipids; unset numBulk
	return $molFraction
	
}


################################
# Returns the number of lipids #
# in the species queried       #
###############################

proc numLipids { targetLipid } {
        set targetLipids [ atomselect top "resname $targetLipid" ]
        set targetLipidAtomNames [ $targetLipids get name ]
        set targetLipidAtomName [ lindex $targetLipidAtomNames 0 ]
	set numLipid [ [ atomselect top "resname $targetLipid and name $targetLipidAtomName" ] num ]
	return $numLipid
}


################################
# Returns the number of atoms/ #
# beads in a lipid type        #
################################

proc lipidAtoms { targetLipid } {
        set targetLipids [ atomselect top "resname $targetLipid" ]
	set targetLipidResid [ $targetLipids get resid ]
	set firstTargetLipid [ atomselect top "resname $targetLipid and resid [ lindex $targetLipidResid 0 ]" ]
	set numAtoms [ $firstTargetLipid num ]
	return $numAtoms
}

#########################################
# Takes a list of atomi indices from    #
# a "within" selection (i.e., lipids    #
# within 5 of protein and produces      #
# a list of unique lipid resid          #
# numbers for the input index list      #
#########################################

proc indexToResid { indexList } {
	set residList ""
	foreach index $indexList {
		set indexAtom [ atomselect top "index $index" ]
		set indexResid [ $indexAtom get resid ] 
		if { [ lsearch -exact $residList $indexResid ] == -1 } {
			lappend residList $indexResid
			#puts "Resid $indexResid added to list"
		}
		$indexAtom delete
		unset indexResid
	}
	return $residList
} 


##########################################
# Takes as input a selection and aligns  #
# the whole simulation to the frame      #
# and selection.                         #
##########################################

proc alignment { selection args } {
	puts "Aligning simulation system to $selection..."
	if { [ llength $args ] > 1 } {
		error "Too many values specified for alignment frame"
	} elseif { [ llength $args ] == 0 } {
		set frameNum 0
		puts "Using default setting of first frame for alignment"
	} elseif { [ llength $args ] == 1 } { 
		set frameNum $args
	}
	set all [ atomselect top "all" ] ;# this aligns the whole system based on the selection criteria
	set refMol [ atomselect top $selection frame $frameNum ]
	
	set nFrm [ molinfo top get numframes ]
	for { set i 1 } { $i < $nFrm } { incr i 1 } {
		#puts "in for loop"
		$all frame $i
		set frmMol [ atomselect top $selection frame $i ]	
		set moveMat [ measure fit $frmMol $refMol ]
		$all move [ measure fit $frmMol $refMol ] 
		$frmMol delete
	}
	unset nFrm
}

#########################################
# Takes as input a text file and splits #
# it and reads it into a list. This     #
# script allows you to select the       #
# column you would like to read and     #
# place in the list                     #
#########################################

proc textToList { fileName column } {
	puts -nonewline "Reading data from $fileName and sorting into list..."
	set f [ open $fileName r ]
	foreach line [ split [ read $f ] \n ] {
		# If statement takes care of extra whitespace at end of file
		if { [ string trim $line ] != "" } { 
			lappend textList [ lindex $line $column ]
		}
	}
	return $textList
}

#########################################
# Takes a PDB and outputs a string of   #
# single letter amino acid code with    #
# appropriate case for input into SCWRL #
#########################################


proc oneLetterCode { pdbFile outputFile residueSelection residuesToScan } {
	set threeLetters "ALA CYS ASP GLU PHE GLY HIS HSD HSE ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"
	set oneLetter "a c d e f g h h h i k l m n p q r s t v w y"
	
	mol load pdb $pdbFile
	set oFile [ open $outputFile w+ ]
	set CA [ atomselect top "name CA and $residueSelection" ]
	set indices [ $CA get index ]
	$CA delete
	foreach ind $indices {
		set residue [ atomselect top "index $ind" ]
		set aAcid [ $residue get resname ]
		set resNum [ $residue get resid ]
		set listPosition [ lsearch -nocase $threeLetters $aAcid ]
		set letterCode [ lindex $oneLetter $listPosition ]
		if { $resNum == $residuesToScan } {
			set letterCode [ string toupper $letterCode ]
		}
		puts "$resNum $letterCode"
		puts $oFile $letterCode
		$residue delete
		unset aAcid resNum listPosition letterCode
	}
	close $oFile
}
