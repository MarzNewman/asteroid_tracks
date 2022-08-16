#!/bin/bash

#PBS -q single
#PBS -A hpc_mishaps04
#PBS -l walltime=06:00:00
#PBS -l nodes=1:ppn=1

#regfile=$1
#headfile=$2

source $pipedir/pip2options.sh

headfile=/project/marz746/astdata/regfiles/MN_asteroid_tracking.txt

#module load python/anaconda-3.6.4

cd $workdir
mkdir astrometry_tmp
cd astrometry_tmp

#${pipedir}astreg2phot.sh $regfile $field $filter
cp ${pipedir}reg_predict.py /var/scratch/astrometry_tmp
#cp ${pipedir}astreg2phot.sh /var/scratch/astrometry_tmp

#python reg_predict.py $astname $numNewRegions $print2reg

#find number of current regions
numCurRegs=$(grep circle /project/marz746/astdata/regfiles/${astname}.reg | wc -l)
#numCurRegs=$((numCurRegs_tmp-1))

#exits if number of regions is greater than 3
if [[ $numCurRegs -gt 3 ]]
then
    echo "Regions file is too long; expecting 3 regions"
    echo "Remove additional regions and try again"
    echo "Exiting"
    exit
fi

#find number of total images
read nname chip night cfirst linked < <(grep $astname $headfile)
echo night: $night
numImages=$(grep ${night} /project/marz746/astdata/MISHAPS_F1_r.good.holdings | wc -l)
echo numImages: $numImages

#loop to make predictions over all images
numNewRegions=1
regfile=${astname}.reg

#read phot file to see where the algorithm fails
#photfile=/project/marz746/astdata/output/${astname}_2phot.txt
photfile=/project/marz746/astdata/output/${astname}_1phot.txt
readphot () {
    photdata=()
    n=0
    while read line
    do
	echo $n $line
	photdata+=( $line )
	n=$(($n+1))
    done < $photfile
    #return $photdata
}

#calculates rms of photdata
rms () {
    j=0
    sumsq=0
    nans=0
    while [[ $j -lt $numCurRegs ]]
    do
	if [[ "${photdata[$j]}" == "nan" ]]
	then
	    nans=$(($nans+1))
	    echo nan encountered during rms func
	    j=$(($j+1))
	    continue
#	else
#	    echo no nan encountered during rms func
	fi
	square=$(echo "${photdata[$j]}^2" | bc)
#	echo j and square: $j and $square
	sumsq=$(echo "$square+$sumsq" | bc)
#	echo sumsq: $sumsq
	j=$(($j+1))
    done
    #rmsPhot=$(echo "${photdata[0]}^2+${photdata[1]}^2+${photdata[2]}^2" | bc)
    numdata=$(($numCurRegs-$nans))
    #numdata=$numCurRegs
    rmsPhot=$(echo "scale=8; ${sumsq}/${numdata}" | bc)
    rmsPhot=$(echo "scale=8; sqrt($rmsPhot)" | bc)
    #echo $rmsPhot
#    return $rmsPhot
}

#function for standard devia843027.smic3tion
std () {
    k=0
    sumsq=0
    nans=0
    while [[ $k -lt $numCurRegs ]]
    do
        if [[ "${photdata[$k]}" == "nan" ]]
        then
            nans=$(($nans+1))
            echo nan encountered during std func
	    k=$(($k+1))
            continue
#	else
#	    echo no nan encountered during std func
        fi
	sub=$(echo "scale=8; ${photdata[$k]}-${rmsPhot}" | bc)
#	echo k and sub: $k and $sub
	subsq=$(echo "scale=8; $sub^2" | bc)
	sumsq=$(echo "scale=8; $subsq+$sumsq" | bc)
#	echo sumsq: $sumsq
	k=$(($k+1))
    done
    numdata=$(($numCurRegs-$nans))
#    numdata=$numCurRegs
    stdev=$(echo "scale=8; sqrt(${sumsq}/${numdata})" | bc)
    #echo $stdev
}

rm /project/marz746/astdata/output/${astname}_predInput.txt
${pipedir}astreg2phot.sh $regfile $field $filter

#calculates rms of photdata after first astrometry run
readphot
rms

#calculates standard deviation of photometry
std
echo stdev: $stdev

echo python /home/marz746/pip2/reg_predict.py $astname $numNewRegions True #debugging
python /home/marz746/pip2/reg_predict.py $astname $numNewRegions True
#exit
${pipedir}astreg2phot.sh $regfile $field $filter

numCurRegs=4

fails=0

echo 1 phot
cat /project/marz746/astdata/output/${astname}_1phot.txt
echo 2 phot
cat /project/marz746/astdata/output/${astname}_2phot.txt
echo numCur and numImages: $numCurRegs and $numImages

while [[ $numCurRegs -lt $numImages ]];
do
#    echo $numCurRegs
#    if [[ $numCurRegs -gt 71 ]]
#    then
#	echo debug exiting
#	exit
#    fi
#    echo algorithm loop
    numNewRegions=$(($numNewRegions+2))
#    echo numCur and numNewRegions: $numCurRegs and $numNewRegions
#    if [[ $numNewRegions == 13 ]];
#    then
#	numNewRegions=$(($numNewRegions-6))
#    fi
    python reg_predict.py $astname $numNewRegions True
    
    #exits if there's a negative values in the regfile
    linecount=0
    while read line
    do 
	if grep -q "-" <<< "$line"
	then
	    echo line: $line
	    echo negative value encountered within regfile
	    cp MN*.txt /project/marz746/astdata/output
	    cp MN*.photast /project/marz746/astdata/output
	    cp *.pdf /project/marz746/astdata/output
	    echo removing lines with negative values from regfile
	    #echo $linecount
	    sed -i "${linecount},$(($numCurRegs+$numNewRegions+3))d" /project/marz746/astdata/regfiles/${astname}.reg
	    echo exiting
	    exit
	fi
	linecount=$(($linecount+1))
    done < /project/marz746/astdata/regfiles/${astname}.reg

    ${pipedir}astreg2phot.sh $regfile $field $filter
    echo astrometry complete
    #if the end of night has been reached, exits script
    if [[ -f "/project/marz746/astdata/output/endOfNight.txt" ]]
    then
	echo End of night reached
	rm /project/marz746/astdata/output/endOfNight.txt
	cp MN*.txt /project/marz746/astdata/output
	cp MN*.photast /project/marz746/astdata/output
	cp *.pdf /project/marz746/astdata/output
	echo exiting
	exit
    fi
    readphot
    #cat $photfile
#    q=0
#    echo photdata[0]: ${photdata[0]}
#    while [[ $q -lt 10 ]]
#    do 
#	echo ${photdata[$q]}
#	q=$(($q+1))
#    done
    i=$(($numCurRegs))
    sigthresh=$(echo "scale=8; 3*${stdev}+${rmsPhot}" | bc)
    updateSTD=1   #updates rms, std, and numCurRegs if equal to 1
    echo numCur and numNewRegs: $numCurRegs and $numNewRegions
    echo sigthresh and stdev and rms: $sigthresh and $stdev and $rmsPhot
    while [[ $i -lt $(($numCurRegs+$numNewRegions)) ]]
    do
	echo i and photdata: $i and ${photdata[$i]}
	#enters loop if the photometry is too much dimmer than the sigma threshold and numNewRegions is greater than 1
	if (( $(echo "${photdata[$i]} > $sigthresh" | bc -l) )) && [ $numNewRegions -gt 1 ]
	then
	    #echo exiting at step $i and removing lines starting with $(($i+3)) from regfile
	    fails=$(($fails+1))
	    echo $fails failed attempts
	    if [[ $fails -gt 70 ]]
	    then
		echo Too many failed predictions
		cp MN*.txt /project/marz746/astdata/output
		cp MN*.photast /project/marz746/astdata/output
		cp *.pdf /project/marz746/astdata/output
		echo exiting
		exit
	    fi
	    #remove bad lines from previous attempt
	    #echo Removing lines starting with $(($NumCurRegs+3)) from regfile
	    echo Removing lines starting with $(($i+4)) from regfile
	    sed -i "$(($i+4)),$(($numCurRegs+$numNewRegions+3))d" /project/marz746/astdata/regfiles/${astname}.reg
#	    cat /project/marz746/astdata/regfiles/${astname}.reg
	    echo Removing bad lines from predInput file
	    sed -i "$(($i+1)),$(($numCurRegs+$numNewRegions))d" /project/marz746/astdata/output/${astname}_predInput.txt

	    #if the prediction fails, backtrack to create less new regions
	    numNewRegions=$(($numNewRegions-4))
	    #numNewRegions=3
	    
	    #sets update std to false so that the standard deviation doesn't include bad data
	    updateSTD=0
	    
	    #resets numCurRegions directly from the regions file so that I don't have to guess the number
	    numCurRegs=$(grep circle /project/marz746/astdata/regfiles/${astname}.reg | wc -l)
	    echo bugfix numCurRegs: $numCurRegs
	    break
	    #echo ${photdata[$i]}
	    #exit
	fi
	i=$(($i+1))
    done
    cat /project/marz746/astdata/regfiles/${astname}.reg
    if [[ $updateSTD == 1 ]]
    then
	rms
	std
	echo stdev: $stdev
	if (( $(echo "$stdev > 1" | bc ) ))
	then 
	    echo stdev is too high
	    cp MN*.txt /project/marz746/astdata/output
	    cp MN*.photast /project/marz746/astdata/output
	    cp *.pdf /project/marz746/astdata/output
	    echo exiting
	    exit
	fi
	numCurRegs=$(($numCurRegs+$numNewRegions))
    fi
    #echo readphot done
    #echo ${photdata[0]}
    ##rms
    #echo rms done
    #echo $rmsPhot
    ##std
    ##echo stdev: $stdev
    #echo ${pipedir}astreg2phot.sh $regfile $field $filter
    ##numCurRegs=$(($numCurRegs+$numNewRegions))
    if [[ $numNewRegions -le 0 ]]
    then
	numNewRegions=-1
    fi
    echo debug 2 numCurRegs: $numCurRegs
done


cp MN*.txt /project/marz746/astdata/output
cp MN*.photast /project/marz746/astdata/output 
cp *.pdf /project/marz746/astdata/output
cd ../
rm -r astrometry_tmp

exit
