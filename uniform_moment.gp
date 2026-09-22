set terminal pdfcairo enhanced size 5in,5in font "Latin Modern Roman,22"
set fit quiet nolog limit 1e-14
set size square
set border lw 1.5
set tics in scale 1.5 nomirror
set key top left reverse Left box opaque font ",18" width -1

ROT = "#DE4A39"
BLAU = "#0067A5"

DELTAS = system("cut -f2 LISTDATA.txt | sort -u -g")
ALPHAS = system("cut -f27 LISTDATA.txt | sort -u -g")
ND = words(DELTAS)
NA = words(ALPHAS)

power(x) = abs(c)*x
power_3(x) = a0 + a1*x + a2*x**2
C(d,a) = ($2==d && $27==a) ? sqrt($35**2+$36**2+$37**2) : 1/0
Q(d,a) = ($2==d && $27==a) ? sqrt($38**2+$39**2+$40**2) : 1/0
S = real(system("awk '/#define S /{print $3}' ../INC/trilattice_sw.h"))
TOT(d,a) =($2==d && $27==a) ? sqrt((S*$35+$38)**2+(S*$36+$39)**2+(S*$37+$40)**2) : 1/0

array A0C[NA*ND]
array A1C[NA*ND]
array A2C[NA*ND]
array A0Q[NA*ND]
array A1Q[NA*ND]
array A2Q[NA*ND]
array A0T[NA*ND]
array A1T[NA*ND]
array A2T[NA*ND]

stats "LISTDATA.txt" using (1.0/$1) nooutput
set xrange [0:*]
set xtics 0, 2*10**floor(log10(STATS_max))
set xlabel "1/{/:Italic L}"

do for [an=1:NA] {
	A = word(ALPHAS,an)+0
	do for [dn=1:ND] {
		D = word(DELTAS,dn)+0
		k = (an-1)*ND + dn

		a0 = 1.0; a1 = 1.0; a2 = 1.0
		fit power_3(x) "LISTDATA.txt" using (1.0/$1):(C(D,A)) via a0,a1,a2
		A0C[k] = a0; A1C[k] = a1; A2C[k] = a2
		set output sprintf("PLOT/UMOM_C-DELTA_%.4f_ALPHA_%.4f.pdf",D,A)
		set ylabel sprintf("{/:Italic \U+03B4}{/:Italic m}^{(0)}({/:Italic \U+03B4}{/:Italic J}/{/:Italic J} = -%g)",D)
		plot "LISTDATA.txt" using (1.0/$1):(C(D,A)) with points pt 4 ps 1.2 lw 2 lc rgb ROT notitle, \
			power_3(x) with lines dt 2 lw 3 lc rgb BLAU title "{/:Italic a}^{(0)}+{/:Italic a}^{(1)}/{/:Italic L}+{/:Italic a}^{(2)}/{/:Italic L}^2"

		a0 = 1.0; a1 = 1.0; a2 = 1.0
		fit power_3(x) "LISTDATA.txt" using (1.0/$1):(Q(D,A)) via a0,a1,a2
		A0Q[k] = a0; A1Q[k] = a1; A2Q[k] = a2
		set output sprintf("PLOT/UMOM_QDELTA_%.4f_ALPHA_%.4f.pdf",D,A)
		set ylabel sprintf("{/:Italic \U+03B4}{/:Italic m}^{(1)}({/:Italic \U+03B4}{/:Italic J}/{/:Italic J} = -%g)",D)
		plot "LISTDATA.txt" using (1.0/$1):(Q(D,A)) with points pt 4 ps 1.2 lw 2 lc rgb ROT notitle, \
			power_3(x) with lines dt 2 lw 3 lc rgb BLAU title "{/:Italic a}^{(0)}+{/:Italic a}^{(1)}/{/:Italic L}+{/:Italic a}^{(2)}/{/:Italic L}^2"

		a0 = 1.0; a1 = 1.0; a2 = 1.0
		fit power_3(x) "LISTDATA.txt" using (1.0/$1):(TOT(D,A)) via a0,a1,a2
		A0T[k] = a0; A1T[k] = a1; A2T[k] = a2
		set output sprintf("PLOT/UMOM_TOT-DELTA_%.4f_ALPHA_%.4f.pdf",D,A)
		set ylabel sprintf("{/:Italic m}_{tot}({/:Italic \U+03B4}{/:Italic J}/{/:Italic J} = -%g)",D)
		plot "LISTDATA.txt" using (1.0/$1):(TOT(D,A)) with points pt 4 ps 1.2 lw 2 lc rgb ROT notitle, \
			power_3(x) with lines dt 2 lw 3 lc rgb BLAU title "{/:Italic a}^{(0)}+{/:Italic a}^{(1)}/{/:Italic L}+{/:Italic a}^{(2)}/{/:Italic L}^2"
	}
}

set print "UMOM.txt"
print "ALPHA DELTA CL QQ"
do for [an=1:NA] { do for [dn=1:ND] { k = (an-1)*ND + dn
	print sprintf("%g %g %.4f %.4f", word(ALPHAS,an)+0, word(DELTAS,dn)+0, A0C[k], A0Q[k]) } }
set print "L_SCALE_C.txt"
print "ALPHA DELTA a0 a1 a2"
do for [an=1:NA] { do for [dn=1:ND] { k = (an-1)*ND + dn
	print sprintf("%g %g %.4f %.4f %.4f", word(ALPHAS,an)+0, word(DELTAS,dn)+0, A0C[k], A1C[k], A2C[k]) } }
set print "L_SCALE_Q.txt"
print "ALPHA DELTA a0 a1 a2"
do for [an=1:NA] { do for [dn=1:ND] { k = (an-1)*ND + dn
	print sprintf("%g %g %.4f %.4f %.4f", word(ALPHAS,an)+0, word(DELTAS,dn)+0, A0Q[k], A1Q[k], A2Q[k]) } }
set print "L_SCALE_TOT.txt"
print "ALPHA DELTA a0 a1 a2"
do for [an=1:NA] { do for [dn=1:ND] { k = (an-1)*ND + dn
	print sprintf("%g %g %.4f %.4f %.4f", word(ALPHAS,an)+0, word(DELTAS,dn)+0, A0T[k], A1T[k], A2T[k]) } }
unset print

set xrange [*:*]
set xtics autofreq
set autoscale fix
set offsets graph 0.05, graph 0.05, graph 0.05, graph 0.05
set xlabel "-{/:Italic \U+03B4}{/:Italic J}/{/:Italic J}"
unset key

set print "DELTA_SCALE_C.txt"
print "ALPHA m0/Delta"
set print "DELTA_SCALE_Q.txt"
print "ALPHA m1/Delta"
set print "DELTA_SCALE_TOT.txt"
print "ALPHA mtot/Delta"
unset print

do for [an=1:NA] {
	A = word(ALPHAS,an)+0
	set print $SUM
	do for [dn=1:ND] { k = (an-1)*ND + dn
		print word(DELTAS,dn)+0, A0C[k], A0Q[k], dn, A0T[k] }
	unset print

	c = 1.0
	fit power(x) $SUM using 1:($4<=DELTA_FIT ? $2 : 1/0) via c
	set print "DELTA_SCALE_C.txt" append
	print sprintf("%g %.4f", A, c)
	unset print
	set output "PLOT/UMOM_C-DELTA_VS_ALPHA.pdf"
	set ylabel "{/:Italic m}_{tot}/{/:Italic S}"
	plot $SUM using 1:2 with points pt 4 ps 1.2 lw 2 lc rgb ROT, \
		power(x) with lines dt 2 lw 3 lc rgb BLAU

	c = 1.0
	fit power(x) $SUM using 1:($4<=DELTA_FIT ? $3 : 1/0) via c
	set print "DELTA_SCALE_Q.txt" append
	print sprintf("%g %.4f", A, c)
	unset print
	set output "PLOT/UMOM_Q-DELTA_VS_ALPHA.pdf"
	set ylabel "{/:Italic \U+03B4}{/:Italic m}^{(1)}"
	plot $SUM using 1:3 with points pt 4 ps 1.2 lw 2 lc rgb ROT, \
		power(x) with lines dt 2 lw 3 lc rgb BLAU

	c = 1.0
	fit power(x) $SUM using 1:($4<=DELTA_FIT ? $2-$3 : 1/0) via c
	set output "PLOT/UMOM_DELTA_VS_ALPHA.pdf"
	set ylabel "{/:Italic m}_{tot}/{/:Italic S}"
	plot $SUM using 1:($2-$3) with points pt 4 ps 1.2 lw 2 lc rgb ROT, \
		power(x) with lines dt 2 lw 3 lc rgb BLAU

	c = 1.0
	fit power(x) $SUM using 1:($4<=DELTA_FIT ? $5 : 1/0) via c
	set print "DELTA_SCALE_TOT.txt" append
	print sprintf("%g %.4f", A, c)
	unset print
	set output "PLOT/UMOM_TOT-DELTA_VS_ALPHA.pdf"
	set ylabel "{/:Italic m}_{tot}"
	plot $SUM using 1:5 with points pt 4 ps 1.2 lw 2 lc rgb ROT, \
		power(x) with lines dt 2 lw 3 lc rgb BLAU
}
