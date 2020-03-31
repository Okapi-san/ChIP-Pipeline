
# -*- coding: None -*-
echo "Intersecting $1..."
INPUT_FILE="G4_intersect_file.G4"

bedtools intersect -a $1 -b $INPUT_FILE -wa > intersected.$1
BASE=$(wc $1)
G4=$(wc intersected.$1)
echo -e "$1 \t $BASE \t $G4" >> Output.txt 
#rm $1
#rm intersected.$1 
