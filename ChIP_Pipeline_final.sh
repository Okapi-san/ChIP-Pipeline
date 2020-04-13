
# -*- coding: None -*-

main(){
set_variables
split_file
extracting_IDs
formatting_file
creating_experiment_files
G4_intersecting
#Jaccard
#multiIntersect
#cleaning


}

set_variables(){
echo "Enter File name:"
read FILE_NAME
echo "Enter number of CPU cores:"
read CORE_COUNT
str=$(tail -n +2 "$FILE_NAME")
if [[ $str == track* ]]; 
then 
	echo "$(tail -n +2 "$FILE_NAME")" > "$FILE_NAME" 
else 
	echo "Header looks good."
fi

}
split_file(){
line_number=$(wc -l < $FILE_NAME)
split_number=10000000
if (( $line_number < $split_number ));
then
	echo "File is too small to split."
else
echo "File size is $line_number lines. File is split."
echo "Splitting File..."
split -l 10000000 -d $FILE_NAME
echo "File split."
fi

}

extracting_IDs(){
count=1
echo "Extracting IDs..."
line_number=$(wc -l < $FILE_NAME)
split_number=10000000
if (( $line_number < $split_number ));
then
	python3 formatting_pipeline_new.py -i $FILE_NAME >> $FILE_NAME.IDs
else
for a in x*;
	do
	python3 formatting_pipeline_new.py -i $a >> $a.IDs
	count=`expr $count + 1`
	echo "File $count processed."
	done
fi

echo "Finished processing a total of $count ID files." 
cat *.IDs >> "$FILE_NAME.IDs.final"
sort "$FILE_NAME.IDs.final" | uniq > $FILE_NAME.IDs
sed -i '/^$/d' $FILE_NAME.IDs
echo "Finished extracting IDs!"

}
formatting_file(){
line_number=$(wc -l < $FILE_NAME)
split_number=10000000
if(( $line_number < $split_number ));
then
	echo "Formatting file..."
	end_of_file=0
	while [[ $end_of_file == 0 ]]; do
  		read -r line
  		end_of_file=$?
		grep SRX | cut -f 1-4 | sed 's/;/\t/g' | cut -f 1-5 | sed 's/ID=//g' |  sed 's/Name=//g' | sed 's/%20(/\t/g' | cut -f 1-5 >> temp
		
	done < "$FILE_NAME"
	echo "File formatted."
else
	echo "Formatting file..."
	for b in x*;
	do
		end_of_file=0
		while [[ $end_of_file == 0 ]]; do
  			read -r line
  			end_of_file=$?
			grep SRX | cut -f 1-4 | sed 's/;/\t/g' | cut -f 1-5 | sed 's/ID=//g' |  sed 's/Name=//g' | sed 's/%20(/\t/g' | cut -f 1-5 >> temp
			
		done < "$b"
	echo "File formatted."
	done
	
fi

}

creating_experiment_files(){
echo "Creating experiment files..."


parallel -j $CORE_COUNT "grep {1} temp | cut -f 1-5  >> {1}.bed" ::: `grep SRX $FILE_NAME.IDs` 
FILE_COUNT=$(ls SRX* | wc -l)
echo "Created a total of $FILE_COUNT experiment files."


}

G4_intersecting(){
#Requires bedtools an GNU parallel. G4 file must be named G4_intersect.sh
echo "Intersecting Experiment files with G4s..."
ls SRX* | parallel -j $CORE_COUNT "bash G4_intersect.sh $(echo {})"

python << END
import pandas as pd
data = pd.read_csv("Output.txt", delim_whitespace = True, header = None)
data[9] = data[5] / data[1]
data = data[[0, 1, 5, 9]]
data = data.sort_values(by = 9, ascending = False)
data.to_csv("Output_sorted.txt", sep = "\t", index = None, header = None)

END

cut -f -1 Output_sorted.txt | awk -F "." '{print $1}' > ID_file.txt
end_of_file=0
while [[ $end_of_file == 0 ]]; do
  		read -r line
  		end_of_file=$?
		wget https://www.ncbi.nlm.nih.gov/sra/$line/
		sed -nr '/Sample:/ s/.*Sample:([^"]+).*/\1/p' index.html |  awk -F ">" '{print $2}' | awk -F "<" '{print $1}' >> temp_title
		rm index.html
	done < "ID_file.txt"
	
paste Output_sorted.txt temp_title > Output_final.txt

rm ID_file.txt
rm temp_title
rm Output.txt
rm Output._sorted.txt
echo "Intersecting finished."
}

Jaccard(){

ls SRX* > filenames
end_of_file=0
while [[ $end_of_file == 0 ]]; do
  	read -r line
  	end_of_file=$?
	sortBed -i $line > $line.sorted
			
done < "filenames"

file_count=$(ls *.sorted | wc -l)
echo "Calculating jaccard indices for $file_count files..."
parallel -j $CORE_COUNT "bedtools jaccard -a {1} -b {2} | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" ::: `ls *.sorted` ::: `ls *.sorted`
echo "Jaccard indices calculated."
find . | grep jaccard | xargs grep "" | sed -e s"/\.\///" | perl -pi -e "s/.bed./.bed\t/" | perl -pi -e "s/.jaccard:/\t/" > temp
grep SRX temp > pairwise.jaccard.matrix
echo "Constructing matrix..."
python make_matrix.py -i pairwise.jaccard.matrix
echo "Matrix constructed."
column_count=$(awk '{ FS = "\t" } ; { print NF}' matrix_final.tsv | head -1)
matrix_file="matrix_final.tsv"
echo "Rendering matrix plot..."
Rscript --vanilla Jaccard_matrix.R $matrix_file $column_count 
echo "Matrix plot rendered."
echo "Cleaning up..."
rm *.jaccard
rm matrix_final.tsv
rm temp 
rm pairwise.jaccard.matrix
echo "Cleaning finished."
#echo "Exiting."


}

multiIntersect(){

multiIntersectBed -i SRX* | cut -f 6- > "$FILE_NAME.matrix"
column_count=$(awk '{ FS = "\t" } ; { print NF}' $FILE_NAME.matrix | head -1)
Rscript --vanilla Jaccard_matrix.R $FILE_NAME.matrix $column_count

}
cleaning(){
echo "Removing temporary files..."
line_number=$(wc -l < $FILE_NAME)
split_number=10000000
if (( $line_number > $split_number ));
then
	rm x*
fi
rm temp
rm $FILE_NAME.IDs.final
rm $FILE_NAME.IDs
rm *.sorted
rm .sorted
rm intersected.*
rm filenames
#rm SRX*
rm Output.txt
#rm matrix

echo "Termporary files removed. Exiting."
}

main
