# each run with 50_000 events
for i in 1 2 4 5 6 7 8 9 10 11 12
do
        python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin1_VBF/bin/madevent launch_dmspin1_mx10_my100.sh
done

for i in 1 2 3 4 5 6
do
	python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin1_VBF/bin/madevent launch_dmspin1_mx100_my1000.sh
done

for i in 1 2 3
do
	python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin1_VBF/bin/madevent launch_dmspin1_mx1000_my5000.sh
done
