#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <region_file> <header_file>"
    echo "Run from directory containing data"
    exit
fi

pipedir=${0%/*}

source $pipedir/pip2options.sh

#############Options###############
field=MISHAPS_F1
#############Options###############


regfile=$1
headfile=$2

astname=${regfile%.reg}


read nname chip night cfirst linked < <(grep $astname $headfile)
first=${cfirst#conv_}
filter=$(echo $first | perl -pne s'/_/ /g; s/\./ /g' | awk '{print $(NF-3)}')
#read conv c4d dday ttime ooi filter v1 cchip rest < <(echo $first | perl -pne s'/_/ /g; s/\./ /g')

echo $chip $first $filter $night

curdir=$(pwd)
isisdir=$HOME/cfht-pip2/isis/decam
refimage=${curdir##*/}.ref.fits
perl -i -pne 's#^MRJ_DIR.+$#MRJ_DIR\t'${isisdir%/*}'\tInstallation directory#; \
    s#^CONFIG_DIR.+$#CONFIG_DIR\t'$isisdir'\tWhere to find the "config" files#; \
    s#^IM_DIR.+$#IM_DIR\t'$curdir'\tDirectory with the images#; \
    s#^INFILE.+$#INFILE\t'$curdir'/dates\tDates of the frames#; \
    s#^REFERENCE.+$#REFERENCE\t'$curdir/$refimage'\tReference image for astrometry#; \
    s#^REF_STACK.+$#REF_STACK\t'$curdir/$refimage'#; \
    s#^REF_SUB.+$#REF_SUB\t'$curdir/$refimage'\tReference image for subtraction#; \
    s#^VARIABLES.+$#VARIABLES\t'$curdir'/phot.data\tCoordinates of objects#' $isisdir/process_config

awk -v first=${first%_[NS]*.fits} '{if($1 ~ first){start=1}; if(start==1) print $1,$2}' ~/Research/${field}_${filter}.good.holdings | perl -pne 's/\.fits\.fz/_'$chip'.fits/g' > dates.tmp

grep '^circle' $regfile | perl -pne 's/circle\(//g; s/,/ /g; s/\)/ /g' | awk '$3>1{print $0}' > $regfile.proc 

paste dates.tmp $regfile.proc | awk 'NF>=5{printf("%8.3f %8.3f %5.0f %5.0f %s 1 0\n",$3,$4,$3,$4,$2)}' > phot.data.tmp

head -n$(wc -l phot.data.tmp | awk '{print $1}') dates.tmp > dates

#paste $regfile.proc dates | awk '{printf("circle(%f,%f,%f) # text={%s}\n",$1,$2,$3,$4)}'
#head -n2 $regfile > 1region.reg
#cat tmpfile >> 1region.reg

rm photast*

cp phot.data.tmp ~/cfht-pip2/isis/decam/phot.data

curdir=$(pwd)
cd ~/cfht-pip2/isis/decam/
./astrom.csh
cd $curdir
cat photast*|awk 'sqrt(($1-$2)**2)<.0001{print $0}' > $astname.photast.tmp
paste $astname.photast.tmp $regfile.proc > $astname.photast
#python ~/cfht-pip2/fitastrom.py $astname
#python ~/cfht-pip2/timefitx.py $astname
#python ~/cfht-pip2/timefity.py $astname
~/cfht-pip2/convert.sh $astname 
python3 ~/cfht-pip2/all_fits.py $astname 1
~/cfht-pip2/convert.sh $astname 
~/cfht-pip2/astreg2phot2.sh $astname.txt asteroid_tracking.txt
