#!/bin/bash
set -e
mkdir -p OBJ BIN DATA/RAW DATA/PLOT
rm -f DATA/RAW/*
source constants.txt

# virtual environment
python3 -m venv ENV
source ENV/bin/activate

pip install --upgrade pip
pip install -r requirements.txt

make -s

if (( STRX == 2 )); then
	XLENG=${#XMMINARRAY[@]}
else
	XLENG=1
fi

cd BIN
NMX=0
for (( iL=0; iL<${#LARRAY[@]}; iL++ )); do
	for (( xL=0; xL<XLENG; xL++ )); do
		for DIS in $DIS_WDTH; do
			for JP in $JP_WDTH; do
				./trilat.x ${LARRAY[iL]} ${ITER_DSARRAY[iL]} $DIS $JP \
					${BOGTOLARRAY[iL]} $NMX \
					$(echo ${XMMINARRAY[xL]} | bc -l) $(echo ${XMMAXARRAY[xL]} | bc -l) \
					$(echo ${YMMINARRAY[xL]} | bc -l) $(echo ${YMMAXARRAY[xL]} | bc -l) \
					${VCUTARRAY[xL]} ${HCUTARRAY[xL]} ${PAR_XARRAY[iL]} \
					$(echo $QQX | bc -l) $(echo $QQY | bc -l)
				NMX=$((NMX+1))
			done
		done
	done
done

cd ../DATA
if (( LSTX != 0 )); then
	cat RAW/DATA_L-* > LISTDATA.txt
	cat RAW/EIGEN_L-* > EIGENDATA.txt
	cat RAW/IPR_L-* > IPRDATA.txt
	cat RAW/ANGLE_L-* > ANGLEDATA.txt
	python ../uniform_moment.py $DELTA_FIT
fi
if (( STRX != 0 )); then
	cat RAW/GRID_L-* > GRIDDATA.txt
fi
