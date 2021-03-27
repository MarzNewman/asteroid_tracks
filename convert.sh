#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <asteroid_name>"
    echo "Run from directory containing data"
    exit
fi

astname=$1
curdir=$(pwd)
refimage=${curdir##*/}.ref.fits

curdir=${refimage%.ref.fits}

#creates new photast file with RA and Dec in degrees appended to the last two columns (used in the all_fits.py)
~/cfht-pip2/xy2sky.sh $refimage $astname.photast 7 ${curdir}.ref.obj.gaia.transform > ${astname}_1RAdec.photast

#creates photast file with RA and Dec in degrees appended to the last two columns for the optimized positions (used in the all_fits.py)
~/cfht-pip2/xy2sky.sh $refimage ${astname}_2.photast 7 ${curdir}.ref.obj.gaia.transform > ${astname}_2RAdec.photast

#extract columns with pixel data
awk '{print $1, $7, $8}' $astname.photast > ${astname}_1xy.txt

#run xy2sky, creates text file with RA and Dec data
xy2sky -aijd -n 7 $refimage @${astname}_1xy.txt > ${astname}_delete.txt
grep -v '^ *$' ${astname}_delete.txt > ${astname}_1RAdec.txt

#run skycoor for RA and Dec in degrees (needed for to_mpc.py)
#awk '{print $1, $2}' ${astname}RAdec.txt > ${astname}RAdec_deg.txt
skycoor -j @${astname}_1RAdec.txt > ${astname}_1RAdec_deg.txt


#extract columns with pixel data for optimized positions
awk '{print $1, $7, $8}' ${astname}_2.photast > ${astname}_2xy.txt

#run xy2sky, creates text file with RA and Dec data for optimized positions
xy2sky -aijd -n 7 $refimage @${astname}_2xy.txt > ${astname}_delete.txt
grep -v '^ *$' ${astname}_delete.txt > ${astname}_2RAdec.txt

#run skycoor for RA and Dec in degrees (needed for to_mpc.py) for optimized positions
skycoor -j @${astname}_2RAdec.txt > ${astname}_2RAdec_deg.txt

#text file for second run of astrometry
#awk '{print $7, $8, $17}' $astname.photast > ${astname}.txt

#concatenates optimized positions for x and y into the same file
#cat ${astname}x.txt
#cat ${astname}y.txt
awk '{print $17}' $astname.photast > regionsize_delete.txt
###paste ${astname}x.txt ${astname}y.txt regionsize.txt > ${astname}.txt 

paste ${astname}_deleteme.txt regionsize_delete.txt > ${astname}.txt 				#deleteme file comes from all_fits.py

#concatenates optimized positions for x and y into the same file (_optx and _opty come from the all_fits.py)
#paste ${astname}_optx.txt ${astname}_opty.txt regionsize.txt > ${astname}.txt  		#this is done within this all_fits.py
