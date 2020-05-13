INPUT_FILE="G4_intersect_file.G4"
START=$(cat $1.bed | cut -f 2-2 | paste -sd+ - | bc)
END=$(cat $1.bed | cut -f 3-3 | paste -sd+ - | bc)
PEAKS=$(wc -l $1.bed | awk '{ print $1 }')
((count = $END - $START))
MEAN=$(echo "scale=2 ; $count / $PEAKS" | bc)
bedtools random -n $(wc -l $1.bed | awk '{ print $1 }') -l $MEAN -g hg19.chrom.sizes > $1.random
bedtools intersect -a $1.random -b $INPUT_FILE -wa > $1.random.int
G4=$(wc -l $1.random.int | awk '{ print $1 }')
echo "scale=2 ; $G4 / $PEAKS" | bc
G4_score=$(echo "scale=2 ; $G4 / $PEAKS" | bc)
echo -e "$1 \t $G4_score" >> Output_random.txt 
echo "Created $1.random"

