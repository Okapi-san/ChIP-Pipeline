
# -*- coding: None -*-

parallel "bedtools jaccard -a {1} -b {2} | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" ::: `ls *.bed` ::: `ls *.bed`
find . | grep jaccard | xargs grep "" | sed -e s"/\.\///" | perl -pi -e "s/.bed./.bed\t/" | perl -pi -e "s/.jaccard:/\t/" > temp
grep SRX temp > pairwise.jaccard.matrix
python make_matrix.py -i pairwise.jaccard.matrix
