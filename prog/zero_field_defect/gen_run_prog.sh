#!/bin/bash -e
STRSW=#STR#;
LSTSW=#LST#;
iter_CL=#ITCLRANGE#;
LARRAY=(#LRANGE#);
iter_DSARRAY=(#ITDSRANGE#);
PAR_XARRAY=(#PARXRANGE#);
ANGDWIDTH=#THETADWIDTH#;
INUM=0;
BOGTOLARRAY=(#BTL#);
GTOL=#GTL#;
VCUTARRAY=(#VCT#);
HCUTARRAY=(#HCT#);
EIGNUM=#EIGNUM#;
QMX=#QMX#;
QMY=#QMY#;
WM=#WM#;
XMMINARRAY=(#XMOMIN#);
XMMAXARRAY=(#XMOMAX#);
XMOMNUM=#XKGRID#;
YMMINARRAY=(#YMOMIN#);
YMMAXARRAY=(#YMOMAX#);
YMOMNUM=#YKGRID#;
MINF=#FREQMIN#;
MAXF=#FREQMAX#;
FREQNUM=#WGRID#;
IPRMIN=#IPRMIN#;
IPRMAX=#IPRMAX#;
IPRGRID=#IPRGRID#;
GWDTH=#GRNWDTH#;
CORES=#PLRUNS#;
EXEC=#BIN#;

ALENG=${#LARRAY[@]};
DLENG=${#iter_DSARRAY[@]};
PLENG=${#PAR_XARRAY[@]};
BLENG=${#BOGTOLARRAY[@]};

XMINLENG=${#XMMINARRAY[@]};
XMAXLENG=${#XMMAXARRAY[@]};
YMINLENG=${#YMMINARRAY[@]};
YMAXLENG=${#YMMAXARRAY[@]};

VCUTLENG=${#VCUTARRAY[@]};
HCUTLENG=${#HCUTARRAY[@]};

SCUTP=2;

if [[ XMINLENG -eq XMAXLENG  ]] && [[ XMINLENG -eq YMINLENG  ]] \
	&& [[ XMINLENG -eq YMAXLENG  ]] && [[ XMINLENG -eq VCUTLENG ]] \
	&& [[ XMINLENG -eq HCUTLENG ]]
then
	if [[ STRSW -eq SCUTP ]]
	then
		XLENG=${XMINLENG};
	else
		XLENG=1;
	fi
else 
	echo "Fix the brillouin zone cut inputs";
	exit 1;
fi

if [[ ALENG -ne DLENG ]]
then
	echo "Fix the L and/or iter_DS inputs";
	exit 1;
fi

if [[ ALENG -ne PLENG ]]
then
	echo "Fix the L and/or PAR_X inputs";
	exit 1;
fi

if [[ ALENG -ne BLENG ]]
then
	echo "Fix the L and/or BOGTOL inputs";
	exit 1;
fi

sed -e "s|#CITER#|${iter_CL}|g" INC/gen_trilattice_sw.hx > INC/trilattice_sw.h;
sed -i "s|#IMIN#|${IPRMIN}|g" INC/trilattice_sw.h;
sed -i "s|#IMAX#|${IPRMAX}|g" INC/trilattice_sw.h;
sed -i "s|#IGRID#|${IPRGRID}|g" INC/trilattice_sw.h;
sed -i "s|#MINF#|${MINF}|g" INC/trilattice_sw.h;
sed -i "s|#MAXF#|${MAXF}|g" INC/trilattice_sw.h;
sed -i "s|#FREQNUM#|${FREQNUM}|g" INC/trilattice_sw.h;
sed -i "s|#XMOMNUM#|${XMOMNUM}|g" INC/trilattice_sw.h;
sed -i "s|#YMOMNUM#|${YMOMNUM}|g" INC/trilattice_sw.h;
sed -i "s|#LISTCALC#|${LSTSW}|g" INC/trilattice_sw.h;
sed -i "s|#SFACTCALC#|${STRSW}|g" INC/trilattice_sw.h;
sed -i "s|#GTOL#|${GTOL}|g" INC/trilattice_sw.h;
sed -i "s|#ANGDWIDTH#|${ANGDWIDTH}|g" INC/trilattice_sw.h;
sed -i "s|#ECOUNT#|${EIGNUM}|g" INC/trilattice_sw.h;
sed -i "s|#QQX#|${QMX}|g" INC/trilattice_sw.h;
sed -i "s|#QQY#|${QMY}|g" INC/trilattice_sw.h;
sed -i "s|#WW#|${WM}|g" INC/trilattice_sw.h;
sed -i "s|#GW#|${GWDTH}|g" INC/trilattice_sw.h;

sed -e "s|#EXEC#|${EXEC}|g" gen_Makefile > Makefile;

source /opt/intel/bin/compilervars.sh intel64;

#export MKL_NUM_THREADS="40";
#export MKL_DOMAIN_NUM_THREADS="MKL_BLAS=40";
#export OMP_NUM_THREADS="40";
#export MKL_DYNAMIC="FALSE";
#export OMP_DYNAMIC="FALSE";
#export MKL_CBWR="SSE2";

make -s;

cd BIN;

for (( iL=0; iL <${ALENG}; iL++ ));
do
	for (( xL=0; xL <${XLENG}; xL++ ));
	do
		for DIS_WDTH in #DSWDRANGE#;
		do
			for JP_WDTH in #JPRANGE#;
			do

				L=${LARRAY[iL]};
				iter_DS=${iter_DSARRAY[iL]};
				PAR_X=${PAR_XARRAY[iL]};
				BOGTOL=${BOGTOLARRAY[iL]};
				XMMIN=`echo ${XMMINARRAY[xL]} | bc -l`;
				XMMAX=`echo ${XMMAXARRAY[xL]} | bc -l`;
				YMMIN=`echo ${YMMINARRAY[xL]} | bc -l`;
				YMMAX=`echo ${YMMAXARRAY[xL]} | bc -l`;
				VCUT=`echo ${VCUTARRAY[xL]} | bc -l`;
				HCUT=`echo ${HCUTARRAY[xL]} | bc -l`;

				./trilat.x ${L} ${iter_DS} ${DIS_WDTH} \
					${JP_WDTH} ${BOGTOL} ${INUM} \
					${XMMIN} ${XMMAX} ${YMMIN} ${YMMAX} \
					${VCUT} ${HCUT} ${PAR_X} &

				# checking the number of background processes
				# and halting the execution accordingly

				INUM=$((INUM+1));

				bground=( $(jobs -p) );

				if (( ${#bground[@]} >= CORES )); then
					wait -n
				fi

			done || exit 1

		done || exit 1

	done || exit 1

done || exit 1

wait

#rm *.x;
cd ..;
