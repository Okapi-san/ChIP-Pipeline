
# -*- coding: None -*-

main(){
set_variables
formatting_file
extracting_IDs
#extracting_Targets
creating_experiment_files
#creating_target_files
G4_intersecting
Get_GSM_IDs
#Get_Sample_Names
#Jaccard
#multiIntersect
cleaning


}

set_variables(){
echo "Enter File name:"
read FILE_NAME
echo "Enter number of CPU cores:"
read CORE_COUNT
#str=$(tail -n +2 "$FILE_NAME")
#if [[ $str == track* ]]; 
#then 
#	echo "$(tail -n +2 "$FILE_NAME")" > "$FILE_NAME" 
#else 
	echo "Header looks good."
#fi
}

formatting_file(){
	echo "Formatting file..."
	end_of_file=0
	while [[ $end_of_file == 0 ]]; do
  		read -r line
  		end_of_file=$?
		grep SRX | cut -f 1-4 | sed 's/;/\t/g' | cut -f 1-5 | sed 's/ID=//g' |  sed 's/Name=//g' | sed 's/%20(/\t/g' | cut -f 1-5 >> temp
		
	done < "$FILE_NAME"
	echo "File formatted."
}

extracting_IDs(){
echo "Extracting IDs..."
grep SRX temp | cut -f 4- | cut -f -1 > $FILE_NAME.IDs.temp
sort "$FILE_NAME.IDs.temp" | uniq > $FILE_NAME.IDs
echo "Finished extracting IDs!"
}

extracting_Targets(){
echo "Extracting Targets..."
grep SRX temp | cut -f 5- > $FILE_NAME.Targets.temp
sort "$FILE_NAME.Targets.temp" | uniq > $FILE_NAME.Targets
echo "Finished extracting Targets!"
}

creating_experiment_files(){
echo "Creating experiment files..."
parallel -j $CORE_COUNT "grep {1} temp | cut -f 1-5  >> {1}.bed" ::: `grep SRX $FILE_NAME.IDs` 
FILE_COUNT=$(ls SRX* | wc -l)
echo "Created a total of $FILE_COUNT experiment files."

}

creating_target_files(){
echo "Creating Target files..."
parallel -j $CORE_COUNT "grep {1} temp | cut -f 1-5  >> {1}.bed" ::: `grep "" $FILE_NAME.Targets` 
echo " Finished creating Target files."

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



cut -f -1 Output_sorted.txt | awk -F "." '{print $1}' > ID_file
rm Output.txt
rm Output._sorted.txt
echo "Intersecting finished."
}
Get_GSM_IDs(){
echo "Creating metadata file..."
cat $FILE_NAME | cut -f 4-4 | cut -d';' -f 1-3 > temp2
sort --parallel=$CORE_COUNT temp2 | uniq > temp3
cat temp3 | sed 's/ID=//g' | sed 's/;/\t/g' | sed 's/Name=//g' | sed 's/Title=//g' | sed 's/%20/ /g' | sed 's/@//g' | sed 's/:/\t/g' | sed 's/(//g' | sed 's/)//g' > temp_metadata
echo "metadata file created."

echo "Pulling GSM IDs..."
end_of_file=0
while [[ $end_of_file == 0 ]]; do
  		read -r line
  		end_of_file=$?
		#Needs more work, but works now.
		#grep $line temp_metadata | head -1 | sed 's/;/\t/g' | sed 's/:/\t/g' | sed 's/Title=//g' | cut -f 6-6 >> GSM_file  
		echo "$line"
	done < "ID_file"

paste Output_sorted.txt GSM_file > Output_sorted_GSM.txt
#rm GSM_file
echo "GSM pulling finished."
}

Get_Sample_Names(){

end_of_file=0
while [[ $end_of_file == 0 ]]; do
  		read -r line
  		end_of_file=$?
		wget https://www.ncbi.nlm.nih.gov/sra/$line/
		sed -nr '/Sample:/ s/.*Sample:([^"]+).*/\1/p' index.html |  awk -F ">" '{print $2}' | awk -F "<" '{print $1}' >> temp_title
		rm index.html
	done < "ID_file"
	
paste Output_sorted_GSM.txt temp_title > Output_final.txt
rm ID_file.txt
rm temp_title
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

rm temp
rm $FILE_NAME.IDs.final
rm $FILE_NAME.IDs
rm $FILE_NAME.IDs.temp
rm *.sorted
rm .sorted
rm intersected.*
rm filenames
#rm SRX*
#rm Output.txt
#rm matrix

echo "Termporary files removed. Exiting."
}

main
