#!/bin/bash -e
STRSW=0;
LSTSW=1;
iter_CL=1e+5;
LARRAY=(12 18 24 30 36 42 48 54 60);
iter_DSARRAY=(100 100 100 50 30 20 15 10 5);
PAR_XARRAY=(8 8 8 8 8 8 8 8 8);
ANGDWIDTH=0.5;
INUM=0;
BOGTOLARRAY=(0.05 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025);
GTOL=1.0e-8;
VCUTARRAY=(0 1);
HCUTARRAY=(1 2);
EIGNUM=6;
QMX=4.0*pi/3.0;
QMY=0.0;
WM=0.0;
XMMINARRAY=(0.0 0.0);
XMMAXARRAY=(1.0/3.0 2.0/3.0);
XMOMNUM=100;
YMMINARRAY=(0.0 0.0);
YMMAXARRAY=(1.0 1.0);
YMOMNUM=100;
MINF=0.0;
MAXF=0.0;
FREQNUM=1;
IPRMIN=0.0;
IPRMAX=4.00;
IPRGRID=200;
GWDTH=0.025;
CORES=10;
EXEC=trilat.x;

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
		for DIS_WDTH in 0.0 0.2 0.4 0.6 0.8 0.99;
		do
			for JP_WDTH in 0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175;
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
