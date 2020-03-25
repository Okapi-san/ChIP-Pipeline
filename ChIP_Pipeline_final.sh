
# -*- coding: None -*-

main(){
set_variables
split_file
extracting_IDs
formatting_file
creating_experiment_files
G4_intersecting
cleaning


}

set_variables(){
echo "Enter File name:"
read FILE_NAME
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
INPUT_FILE="temp"
count1=0
while IFS= read -r line
	do	
		OUT_FILE="$line.bed"
		grep "$line" $INPUT_FILE | cut -f 1-5  >> "$OUT_FILE"
		count1=`expr $count1 + 1`
	done < "$FILE_NAME.IDs"
echo "Done. A total of $count1 experiment files was created."

}

G4_intersecting(){
#Requires bedtools an GNU parallel. G4 file must be named G4_intersect.sh
echo "Intersecting Experiment files with G4s..."
ls SRX* | parallel -j 3 "bash G4_intersect.sh $(echo {})"

}
cleaning(){
echo "Removing temporary files..."
line_number=$(wc -l < $FILE_NAME)
split_number=10000000
if (( $line_number > $split_number ));
then
	rm x*
fi
#rm temp
#rm $FILE_NAME.IDs.final
echo "Termporary files removed. Exiting."
}

main
