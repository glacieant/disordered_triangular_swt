#!/bin/bash -e
cd DATA;
LCALC=#LIST#;
SCALC=#STRFC#;
if (( ${LCALC} != 0 )); then
	> LISTDATA.dat;
	> EIGENDATA.dat;
	> IPRDATA.dat;
	> ANGLEDATA.dat;
else
	if [ -e LISTDATA.dat ]; then
		rm LISTDATA.dat;
	fi
	if [ -e EIGENDATA.dat ]; then
		rm EIGENDATA.dat;
	fi
	if [ -e IPRDATA.dat ]; then
		rm IPRDATA.dat;
	fi
	if [ -e ANGLEDATA.dat ]; then
		rm ANGLEDATA.dat;
	fi

fi
if (( ${SCALC} == 1 || ${SCALC} == 2 )); then
	> GRIDDATA.dat;
elif [ -e GRIDDATA.dat ]; then
	rm GRIDDATA.dat;
fi
if (( ${LCALC} != 0 )); then
	cat *DATA_L-* >> LISTDATA.dat;
	cat *EIGEN_L-* >> EIGENDATA.dat;
	cat *IPR_L-* >> IPRDATA.dat;
	cat *ANGLE_L-* >> ANGLEDATA.dat;
fi
if (( ${SCALC} == 1 || ${SCALC} == 2 )); then
	cat *GRID_L-* >> GRIDDATA.dat;
fi
if (( ${LCALC} != 0 )); then
	sort -k1,1n -k2,2n -k27,27n LISTDATA.dat -o LISTDATA.dat;
	sort -k1,1n -k2,2n -k3,3n EIGENDATA.dat -o EIGENDATA.dat;
	sort -k2,2n -k6,6n -k1,1n -k3,3n IPRDATA.dat -o IPRDATA.dat;
	sort -k1,1n -k2,2n -k3,3n -k4,4n ANGLEDATA.dat -o ANGLEDATA.dat;
fi
if (( ${SCALC} == 1 || ${SCALC} == 2 )); then
	sort -u -k1,1n -k2,2n -k8,8n -k3,3n -k5,5n -k4,4n GRIDDATA.dat -o GRIDDATA.dat;
fi
if (( ${LCALC} != 0 )); then
	rm *DATA_L-*;
	rm *EIGEN_L-*;
	rm *IPR_L-* 
	rm *ANGLE_L-* 
fi
if (( ${SCALC} == 1 || ${SCALC} == 2 )); then
	rm *GRID_L-*;
fi
cd ..;
