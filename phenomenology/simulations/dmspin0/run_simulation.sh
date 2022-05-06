for i in $(seq 1 $1) 
do
        python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/bin/madevent launch_dmspin0_mx10_my100.sh
done

for i in $(seq 1 $2) 
do
        python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/bin/madevent launch_dmspin0_mx100_my1000.sh
done

for i in $(seq 1 $3) 
do
	python2 /cms/mc/MG5_aMC_v3_1_1/DMsimp_spin0_VBF/bin/madevent launch_dmspin0_mx1000_my5000.sh
done