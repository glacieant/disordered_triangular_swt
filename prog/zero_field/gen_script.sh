#!/bin/bash -e
export LC_ALL=en_US.UTF8

##################################################
#  Input script file - All formats of inputs     # 
#  are supported. However, the code might        #
#  break down if integer inputs are not provided #
#  for integer parameters.                       #
#  These scripts are bash 4.3+ compliant         #
##################################################

# Specifying the lattice lengths for the calculation
L="60 66";
# Specifying the maximum number of iterations for the classical
# algorithm
iter_CL="1e+5";
# The number of disorder iterations for the specific lattice length, L
# following the same syntax of L
iter_DS="20 15";
# Number of parallel runs for the specific lattice length, L
# with the same syntax as L
PAR_X="8 8";
# The width of the randomness, DELTA of the couplings J(1+DELTA)
DIS_WDTH="0.0 0.2 0.4 0.6 0.8 0.99";
# The ratio of the couplings J_prime/J for the nnn
JP_WDTH="0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.175";
# The width of the angular fulctuation from 120 deg state
# Number of parallel execution of processes requested
PLRUNS="10";
# in the unit of pi
ANGDWIDTH="0.5";
#The minimum of frequency for IPR calculation
IMIN="0.0";
#The minimum of frequency for IPR calculation
IMAX="4.00";
#The grid of frequencies for IPR calculation
IGRID="200";
# The minimum of the frequency for BZ calculation
MINF="0.0";
# The maximum of the frequency for BZ calculation
MAXF="0.0";
# Number of frequency grid
FREQNUM="1";
# Switch for calculating observables list
LIST="1";
# Switch for structure factor calculation. STRFC=1 will output
# the structure factor inside a patch, STRFC=2 will output the
# structure factor along a cut and STRFC=0 would save time by not
# computing the structure factor. Please do not use any other
# integer in this switch as that might break the program!
STRFC="0";
# In this calculation with STRFC=2 the output would be
# produced along a cut through the brillouin zone given by
#
#	n_1(p \vec{b}_1 + q \vec_{b}_2)
#
# with \vec{b}_1=(2*pi/a,-2*pi/sqrt(3)/a) and \vec{b}_2=(0,4*Pi/sqrt(3)/a) because
# of oversampling the lattice in a rectangular grid
# In the following the p=hcut and q=vcut is specified.
hcut="1 2";
vcut="0 1";
# The minimum extent of the stretch along \hat{b}_1 within 
# the brillouin zone (in the unit of 2*pi)
XMINMO="0.0 0.0";
# The maximum extent of the stretch along \hat{b}_1 within 
# the brillouin zone (in the unit of 2*pi)
XMAXMO="1.0/3.0 2.0/3.0";
#Grid frequency along \vec{b}_1 (used only to subsample when <L)
XMOMNUM="100";
#The minimum extent of the stretch along \hat{b}_2 within
# the brillouin zone (in the unit of 2*pi) 
YMINMO="0.0 0.0";
#The maximum extent of the stretch along \hat{b}_2 within
# the brillouin zone (in the unit of 2*pi)
YMAXMO="1.0 1.0";
#Grid frequency along \vec{b}_2 (used only to subsample when <L)
YMOMNUM="100";
#The width of the pole of the greens function
GWDTH="0.025";
#Tolerance for bogoliubov subroutine for respective L
BOGTOL="0.025 0.025";
#The global tolerance
GTOL="1.0e-8";
#The wavefunction to visualise
EIGNUM="6";
#The momentum and frequency point to measure structure factor scaling
QMX="4.0*pi/3.0";
QMY="0.0";
WM="0.0";
#The name of the executable program
EXEC="trilat.x";

#### Creatingthe directory to write data into ########

mkdir -p DATA;
mkdir -p DATA/RAW;
mkdir -p SRC;
mkdir -p INC;
mkdir -p BIN;
mkdir -p OBJ;

#### Editing secondoary scripts to perform calculation #####

sed -e "s|#LRANGE#|${L}|g" gen_run_prog.sh > run_prog_temp;
sed -i "s|#ITCLRANGE#|${iter_CL}|g" run_prog_temp;
sed -i "s|#PLRUNS#|${PLRUNS}|g" run_prog_temp;
sed -i "s|#ITDSRANGE#|${iter_DS}|g" run_prog_temp;
sed -i "s|#PARXRANGE#|${PAR_X}|g" run_prog_temp;
sed -i "s|#DSWDRANGE#|${DIS_WDTH}|g" run_prog_temp;
sed -i "s|#JPRANGE#|${JP_WDTH}|g" run_prog_temp;
sed -i "s|#IPRMIN#|${IMIN}|g" run_prog_temp;
sed -i "s|#IPRMAX#|${IMAX}|g" run_prog_temp;
sed -i "s|#IPRGRID#|${IGRID}|g" run_prog_temp;
sed -i "s|#FREQMIN#|${MINF}|g" run_prog_temp;
sed -i "s|#FREQMAX#|${MAXF}|g" run_prog_temp;
sed -i "s|#WGRID#|${FREQNUM}|g" run_prog_temp;
sed -i "s|#XMOMIN#|${XMINMO}|g" run_prog_temp;
sed -i "s|#XMOMAX#|${XMAXMO}|g" run_prog_temp;
sed -i "s|#XKGRID#|${XMOMNUM}|g" run_prog_temp;
sed -i "s|#YMOMIN#|${YMINMO}|g" run_prog_temp;
sed -i "s|#YMOMAX#|${YMAXMO}|g" run_prog_temp;
sed -i "s|#YKGRID#|${YMOMNUM}|g" run_prog_temp;
sed -i "s|#LST#|${LIST}|g" run_prog_temp;
sed -i "s|#STR#|${STRFC}|g" run_prog_temp;
sed -i "s|#VCT#|${vcut}|g" run_prog_temp;
sed -i "s|#HCT#|${hcut}|g" run_prog_temp;
sed -i "s|#BTL#|${BOGTOL}|g" run_prog_temp;
sed -i "s|#GTL#|${GTOL}|g" run_prog_temp;
sed -i "s|#THETADWIDTH#|${ANGDWIDTH}|g" run_prog_temp;
sed -i "s|#EIGNUM#|${EIGNUM}|g" run_prog_temp;
sed -i "s|#QMX#|${QMX}|g" run_prog_temp;
sed -i "s|#QMY#|${QMY}|g" run_prog_temp;
sed -i "s|#WM#|${WM}|g" run_prog_temp;
sed -i "s|#BIN#|${EXEC}|g" run_prog_temp;
sed -e "s|#GRNWDTH#|${GWDTH}|g" run_prog_temp > run_prog.sh;

rm run_prog_temp;

sed -e "s|#LIST#|${LIST}|g" gen_merge_file.sh > merge_temp;
sed -e "s|#STRFC#|${STRFC}|g" merge_temp > merge_file.sh;

rm merge_temp;

chmod 777 run_prog.sh;

./run_prog.sh;

chmod 777 merge_file.sh;

./merge_file.sh;
