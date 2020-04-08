# ChIP-Pipeline
While working on my Master thesis, I was searching for ChIP-seq data to compare to my own data - I quickly discovered databases such as the GEO- and ChIP-Atlas-repository, which have close to half a million ChIP-seq experiments stored in total, quality-controlled and publicly accessible. As most journals now require submission of raw data upon publication, these repositorys will continue to grow in both size and range. With todays computational capacity, it is possible to analzye datasets of increasing size and complexity, enabling novel approaches such as untargeted clustering and principal component analysis (PCA). However, big data requires efficient automatization of data procession and analysis, which can be challenging to achieve. Here, I aim to construct a Pipeline, which automizes the ChIP-Seq analysis and provides a platform to easily implement analysis functions.

## Main Objective

The main objective of the Pipeline is to automate the process of data-gathering and data-analysis: The user should be able to insert data into the Pipeline while specify the analysis methods, and then receive a detailed evaluation on the input files' metadata (Number of experiments, quality control, spread of reads/coverage..) and an analysis based on the methods specified during the setup. For now, analysis functions of the Pipeline will focus on identifying factors possibly interacting with G4-quadruplexes and clustering these factors by similiar region targeting.

## Prerequisites & Install

**Requires:**
- R > 4.5
- Python 3
- GNU parallel
- bash or similiar

**Install:**
Clone the repository to your local workspace
```sh
mkdir PIPELINE
cd PIPELINE
git clone https://github.com/Okapi-san/ChIP-Pipeline.git
chmod u+x ChIP_Pipeline_final.sh G4_intersect.sh Jaccard_matrix.R  formatting_pipeline_new.py jaccard.sh make_matrix.py
```
Run the Pipeline with ```./ChIP_Pipeline_final.sh```.

# Workflow

## Input Data
The pipeline requires Input files in BED- or BED6-format or RAW-files from ChIP-Atlas. Input as .tsv from ChIP-Atlas looks like this (showing one line):
```
track name="EPAS1 (@ All cell types) 50" url="http://chip-atlas.org/view?id=$$" gffTags="on"
chr1	9869	10464	ID=SRX968419;Name=EPAS1%20(@%20786-O);Title=GSM1642766:%20ChIP-Seq%20of%20HIF-2a%20in%20786-O%20with%20HIF-1a%20re-expression%3B%20Homo%20sapiens%3B%20ChIP-Seq;Cell%20group=Kidney;<br>source_name=Renal%20Cancer%20Cell%20Line;cell%20line=786-O;transfection=HIF-1a%20(pRRL-HIF-1a);chip%20antibody=HIF-2alpha%20Antibody%20(PM9);	700	.	9869	10464	204,255,0
```
The files generally follow this scheme: 	
| Column        | Description           | Example  |
| ------------- |:-------------:|:-----:|
|   1     | Experimental ID (SRX, ERX, DRX) | SRX097088 |
| 2     | Genome assembly      |   hg19 |
| 3 | Antigen class     |    TFs |
|   4     | Antigen | GATA2 |
| 5     | Cell type class        |   Blood |
| 6 | Cell type    |    K-562 |
|   7     | Cell type description | Primary Tissue=Blood|Tissue Diagnosis=Leukemia Chronic Myelogenous |
| 8     | Processing logs (# of reads, % mapped, % duplicates, # of peaks [Q < 1E-05])      |   30180878,82.3,42.1,6691 |
| 9 | Title     |    GSM722415: GATA2 K562bmp r1 110325 3|
|   10     |Meta data submitted by authors| source_name=GATA2 ChIP-seq K562 BMP cell line=K562 chip antibody=GATA2 antibody catalog number=Santa Cruz SC-9008 |

This line contains all information on one read from the GSM-repository GSM722415, more precise from the experiment SRX097088, such as the chromosome, start- and stop-position, cell line and available metadata. For now, the pipeline extracts the chromosome, start, end and SRX-ID for each read.<br/>
ChIP-Atlas provides data in a concatenated format, with all reads from all experiments in one file. This usually results in huge textfiles ( >100 Gb), that cannot be handled by Python very well and/or exceed the available RAM of one workstation. Thus, the Pipeline splits the file in chunks of 10.000.000 lines into smaller files. This also allows the splitfiles to be processed in parallel (not yet implemented). 
### split(), extracting_IDs() and formatting_pipeline_new.py
```sh

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

```
The splitted files are then fed into a Python script that extracts all SRX-IDs from the input-file and writes them to a new file. Lastly, it joins all IDs from the split-files and removes duplicates: 

```python
import pandas as pd
import numpy as np
import regex as re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, dest="i", help="Specify name of input file.")
args = parser.parse_args()

ID_list=[]

def ID_Check(ID):
    if ID in ID_list:
        return False
    else:
        return True
    
with open(args.i, "r") as file:
    entries = file.readlines()
    for lines in entries:
        gene_ID = "".join(re.findall('ID=(.*?);', lines))
        if ID_Check(gene_ID) == True:
            ID_list.append(gene_ID)            
        else:
            continue
for element in ID_list:
    print(element)
```

## formatting_file()
Next up, the Pipeline removes all excess information from the mainfile to create a BED-file, with four columns for chromosome, start, end and SRX-ID specifying each read. This decreases the file size by a factor of ~8 and thus increases file handling speeds.
```sh
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
```
## creating_experiment_files()
Next, the Pipeline creates a separate file for each SRX-ID present in the mainfile and proceeds to extract each read for a given SRX-ID into the respective file. This allows the user to work with a subset of experiments, which is more ressource effective than extracting the relevant reads from the mainfile. 

```java
               +----------->  SRX19382
+--------------+
|              +----------->  SRX19383
|  GSM722415   |
|              +----------->  SRX19384
+--------------+
               +----------->  SRX19385

```
_Fig. 1: Extracting SRX files from the mainfile._

The ```creating_experiment_files()``` function reads the mainfile line by line with ```while IFS= read -r line```. While this is slower than reading the full file into the memory, it is the only option to process files exceeding the available memory of the workstation.


```sh
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
```
## UPDATE: 
While the method above is very fast for smaller files, it proved too slow for big files, as the complete file needs to be searched for every single SRX-ID. I restructured this method for files >1.000.000 lines, so that ```parallel -a``` reads a single line, and writes it to the corresponding ID file. This way, the file is only processed once, with all cores in parallel. However, for smaller files the old method still should be faster due to the fixed cost degression of opening and closing a file each line.
```sh
+----+                 +----+
|SRX1| +---+---------> |SRX1|
|    |     ^           +----+
|SRX2| +-------+-----> |SRX2|
|    |     |   ^       +----+
|SRX1| +---+   |
|    |         |       +----+
|SRX4| +-------------> |SRX4|
|    |         |       +----+
|SRX2+---------+
+----+

```
```sh
# -*- coding: None -*-	
parallel -j 60 -a temp ./make_file.sh {}
```
This is ```make_file.sh```:
```sh
# -*- coding: None -*-			
a=$(echo "$1" | cut -f 4- | cut -f -1)
b=$(echo "$1" | cut -f -3)
echo $b >> "$a.bed"
```
## Challenges
- **Big data:** 
- **Filtering for unbiased Experiments:** The most challenging problem lies in the data itself: How can the pipeline identify and remove biased experiments (stress conditions, ChIP-seqs of KO mutants ...)
- **Experiment Bias:** 
## Analysis Functions

## G4 Intersect

My main reason constructing this Pipeline was to establish a high-throughput methot identifying factors with a tendency to colocate with G4 quadruplexes *in vivo*. Technically, this question can be adressed quite easily: An intersect function such as ```bedtools intersect``` can rapidly identfy intersecting reads in an experiment file and a file with G4 reads. This is implemented in the ```G4_intersecting()``` function of the Pipeline - using ```GNU parallel``` to enhance the performance. The result of each intersect is written to a summary file, expanded with additional metadata, in the following scheme:<br/>
| SRX-ID        | Total reads       | Total G4 intersects  |Intersects / Reads | Tissue | Sample Title|
| ------------- |:-------------:| :-----:|-----:| :------:|:------------|
|SRX2346892.bed|90|190|2.1111111111111|Kidney|RCC4_Normoxia_HIF-2a (PM9)_Rep 1|

However, this approach does not account for normalization and - more importantly - experimental bias and thus works best in a controlled environment of a small number of experiments. To illustrate this problem, we will use the ChIP-Pipeline to analyze 35 ChIP-Seqs of **EPAS1** (Endothelial PAS domain-containing protein 1), as a sample analysis. After running, the Pipeline, the summary file with the intersect data looks like this (showing first ten lines):

```sh
(base) Oth.ALL.05.EPAS1.AllCell.bed>head -10 Output_sorted.txt 
SRX2346892.bed	90	190	2.111111111111111	Cardiovascular	EPAS1_ChIPSeq
SRX968419.bed	202	314	1.5544554455445545	Kidney		ChIP-Seq of HIF-2a in 786-O with HIF-1a re-expression
SRX968415.bed	277	297	1.0722021660649819	Kidney		ChIP-Seq of HIF-2a in 786-O Vector Alone
SRX4802363.bed	2208	1548	0.7010869565217391	liver		HepG2_16hrs 0.5% O2_HIF-2a (PM9)_Rep 1
SRX3346354.bed	1432	950	0.6634078212290503	liver		ChIP-Seq_786-O_HIF2A
SRX4802309.bed	279	185	0.6630824372759857	Kidney		HKC8_0.5_16hr_HIF2a_rep2
SRX4802349.bed	1663	928	0.5580276608538786	Kidney		RCC4_Normoxia_HIF-2a (PM9)_Rep 1
SRX4802308.bed	204	112	0.5490196078431373	Kidney		HKC8_0.5_16hr_HIF2a_rep1
SRX4802350.bed	1510	827	0.547682119205298	Kidney		RCC4_Normoxia_HIF-2a (PM9)_Rep 2
SRX4802364.bed	1411	744	0.5272856130403969	Kidney		HepG2_16hrs 0.5% O2_HIF-2a (PM9)_Rep 2
```
<br/>
As we can see, we have a huge spread of intersect, ranging from 2.1 - 0.5 (the experiments not shown in this list go as low as 0.05). There is a huge spread in the absolute number of peaks - obviously, SRX2346892 with 90 reads in total is not as robust as SRX4802363 with 2208 reads in total - which might require some sort of graded/adjusted intersect value. However, ```RCC4_Normoxia_HIF-2a (PM9)_Rep 1``` specifies  this would create a bias towards abundant factors. <br/>
The other problem is a general bias of experiments: The 6th corner consists of the sample titles corresponding to each experiment - following the GEO naming convention, ```ChIP-Seq_786-O_HIF2A``` indicates a plain ChIP-Seq of EPA1, whereas ```HepG2_16hrs 0.5% O2_HIF-2a (PM9)_Rep 1``` indicates some kind of treatment.

```sh
G4_intersecting(){
#Requires bedtools an GNU parallel. G4 file must be named G4_intersect.sh
echo "Intersecting Experiment files with G4s..."
ls SRX* | parallel -j 3 "bash G4_intersect.sh $(echo {})"

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


```
## Jaccard-Matrix
When analyzing multiple ChIP-seq experiments, an overall plot of similiarities between single experiments can indicate factors targeting similiar genomic regions or sites. Rather than comparing intersects, one would like to compare the *grade of similiarity*, which can be calculated from the area of overlap divided by the area of union - this is known as the Jaccard Index: <br/><br/> 
![alt text](./jaccard.jpg)
<br/><br/>
The Jaccard statistic is used in set theory to represent the ratio of the intersection of two sets to the union of the two sets. Similarly, Favorov et al. (2012) reported the use of the Jaccard statistic for genome intervals: specifically, it measures the ratio of the number of intersecting base pairs between two sets to the number of base pairs in the union of the two sets. Here, the ```bedtools jaccard``` tool is implemented to calculate this statistic, which modifies the statistic such that the length of the intersection is subtracted from the length of the union. **As a result, the final statistic ranges from 0.0 to 1.0, where 0.0 represents no overlap and 1.0 represent complete overlap.** To plot the the Jaccard indices of **n** ChIP-seq experiments, the Pipeline creates a two-dimensional matrix with an length **n** and width **n** - As this requires the calculation of **n** <sup> **2**</sup> Jaccard indices, the function rapidly outscales the processing power of the workstation and needs to be parallelized with ```GNU parallel``` for improved scalability.  

```python

     d    c    b    a
   +-------------------+        
   |0.1 |0.7 |0.8 |1.0 | a
   +-------------------+
   |0.1 |0.5 |1.0 |0.8 | b           image()
   +-------------------+        +---------------+
   |0.1 |1.0 |0.5 |0.7 | c                      |
   +-------------------+                        |
   |1.0 |0.1 |0.1 |0.1 | d                      |
   +-------------------+                        v
```
![alt text](./matrix_final.jpg)

This is the ```Jaccard()``` function, which utilizes ```GNU parallel``` to parallelize calculation of Jaccard indices. A ```Python``` script constructs an **n x n** matrix from all Jaccard indices, then feeds the matrix into an ```R``` script using  ```scan()```. Finally, ```image()``` interpretes the Jaccard indices in the matrix as colours, ranging from **0.00 (red)** to **1.00 (white)** - hence, yellow to white pixel indicate a high similiarity, whereas red indicates small to no overlap. As the matrix has all experiments on both the x- and y-axis, at each position **(n/n)** there is a white spot, where the Jaccard index was calculated for identical experiments, leading to an overlap of **1.00**. This line thus is technically a symmetry axis.<br/>
Limiting factor for the ```Jaccard()``` function is for now the ```scan()``` function: As it reads in the matrix as a whole into the memory, **n** <sup> **2**</sup> matrix elements quickly outscale available RAM, which kills the ```R``` script. Here, I might need to transform the Jaccard matrix into a sparse matrix, which would greatly reduce RAM load.  
```sh
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
parallel "bedtools jaccard -a {1} -b {2} | awk 'NR>1' | cut -f 3 > {1}.{2}.jaccard" ::: `ls *.sorted` ::: `ls *.sorted`
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
```
The ```Python``` script which constructs the Jaccard matrix using ```sort_index()``` and ```unstack()```.
```python
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, dest="i", help="Specify name of input file.")
args = parser.parse_args()

file = pd.read_csv(args.i, sep = "\t", header = None)
matrix = file.set_index([0, 1])[2].sort_index().unstack()
matrix.to_csv("matrix_final.tsv", index = False, header = False, sep = "\t")
```
The ```R``` script that plots the image using the ```image()``` function.
```r
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args[1]
args[2]
b <- as.numeric(args[2])
c <- args[1]
d <- ".jpg"
e <- paste(c, d)
x <- scan(args[1])
mat <- matrix(x, ncol = b, byrow = TRUE)
jpeg(file=e, width=1600, height=1600, res=300)
image(mat)
dev.off()
q()
```

## Plotting Coverage

Where the Jaccard index calculates an overall plot of similiarities between single experiments to cluster factors targeting similiar genomic regions or sites, another approach would be to plot the overall coverage of each experiment over *all* genomic regions covered in all experiments: This allows to identify genomic regions that are covered in most experiments or, to identify experiments covering special regions which are not covered by most other experiments. This is done with ```bedtools multiIntersectBed```.
```python
 SRX1094                                                                       SRX1094     SRX1095  +->   ...
+-----------------------+			  +-------------------------------------------------------------+
|chr1      199     2041 |                         |chr1      199     2041         1           0
|chr3      55911   59451|    +------------->      |chr1      2849    4822         0           1
|chr5      4280    10234|                         |chr1      32011   34331        0           1
+-----------------------+   multiIntersectBed()   |chr3      55911   59451        1           0
                                                  |chr5      4280    10234        1           0
 SRX1095                           +------->      |chr7      19331   20144        0           1
+-----------------------+          |              |                                     +
|chr1      2849    4822 |          |              |                        +            |
|chr1      32011   34331|   +------+              +	    matrix()       |            v
|chr7      19331   20144|                             +--------------------+
+-----------------------+                             |     image()                    ...
                                                      |
                                                      v
```
![alt text](./matrix_test.jpg)
```sh
multiIntersect(){

multiIntersectBed -i SRX* | cut -f 6- > "$FILE_NAME.matrix"
column_count=$(awk '{ FS = "\t" } ; { print NF}' $FILE_NAME.matrix | head -1)
Rscript --vanilla Jaccard_matrix.R $FILE_NAME.matrix $column_count

}

```

```r
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args[1]
args[2]
b <- as.numeric(args[2])
c <- args[1]
d <- ".jpg"
e <- paste(c, d)
x <- scan(args[1])
mat <- matrix(x, ncol = b, byrow = TRUE)
jpeg(file=e, width=1600, height=1600, res=300)
image(mat)
dev.off()
q()
```

<br/><br/>
## Sample Analysis
In the following, we will use the ChIP-Pipeline to analyze 35 ChIP-Seqs of **EPAS1** (Endothelial PAS domain-containing protein 1), as a sample analysis. Once started, the Pipeline ask to specify the input file and then requires no further interaction. This is the output of the Pipeline in the terminal:
The ```echo``` output resembles the respective steps in the Pipeline:
- Asking for the input file
- Checking the header
- Deciding whether the file should be split
- Extracting all SRX-IDs
- Removing excess data from mainfile
- Creating experiment files
- Intersecting files with G4s
- Calculating Jaccard indices
- Constructing Jaccard matrix
- Rendering matrix plot
- MultiIntersecting & rendering coverage plot

```sh
(base) Oth.ALL.05.EPAS1.AllCell.bed>./ChIP_Pipeline_final.sh 
Enter File name:                        
Oth.ALL.05.EPAS1.AllCell
Header looks good.
File is too small to split.
Extracting IDs...
Finished processing a total of 1 ID files.
Finished extracting IDs!
Formatting file...
File formatted.
Creating experiment files...
Done. A total of 35 experiment files was created.
Intersecting Experiment files with G4s...
Intersecting SRX114493.bed...
Intersecting SRX212355.bed...
Intersecting SRX212356.bed...
Intersecting SRX212357.bed...
Intersecting SRX212358.bed...
Intersecting SRX212363.bed...
Intersecting SRX212365.bed...
Intersecting SRX212364.bed...
Intersecting SRX212366.bed...
Intersecting SRX2346892.bed...
Intersecting SRX3051209.bed...
Intersecting SRX2584128.bed...
Intersecting SRX3346350.bed...
Intersecting SRX359885.bed...
Intersecting SRX3346354.bed...
Intersecting SRX4096727.bed...
Intersecting SRX4802301.bed...
Intersecting SRX4802302.bed...
Intersecting SRX4802309.bed...
Intersecting SRX4802308.bed...
Intersecting SRX4802315.bed...
Intersecting SRX4802322.bed...
Intersecting SRX4802316.bed...
Intersecting SRX4802323.bed...
Intersecting SRX4802335.bed...
Intersecting SRX4802337.bed...
Intersecting SRX4802342.bed...
Intersecting SRX4802344.bed...
Intersecting SRX4802349.bed...
Intersecting SRX4802350.bed...
Intersecting SRX4802363.bed...
Intersecting SRX968415.bed...
Intersecting SRX4802364.bed...
Intersecting SRX968419.bed...
Intersecting SRX968422.bed...
Calculating jaccard indices for 35 files...
Jaccard indices calculated.
-i used with no filenames on the command line, reading from STDIN.
-i used with no filenames on the command line, reading from STDIN.
Constructing matrix...
Matrix constructed.
Rendering matrix plot...
[1] "matrix_final.tsv"
[1] "35"
Read 1225 items
null device 
          1 
Matrix plot rendered.
Cleaning up...
Cleaning finished.
[1] "Oth.ALL.05.EPAS1.AllCell.matrix"
[1] "70"
Read 6047650 items
null device 
          1 
Removing temporary files...
rm: cannot remove 'temp': No such file or directory
Termporary files removed. Exiting.
```
Running the Pipeline took around 20 seconds from start to finish. Before exiting, all intermediary files such as temp files, jaccard matrices, intersect files a. o. are removed, as these are usually not required for further analysis and would only occupate storage space. However, these files can be kept by commenting out the corresponding line in the ``` cleaning()``` function, such as``` #rm intersected.* ```. We then can see 35 distinct SRX-files, containing the reads of each experiment. The ```Output_sorted.txt``` file contains a summary of the G4-intersect for each file, and sorts them according to the relative overlap, descending from the highest overlap. The file is structured as:
| SRX-ID        | Total reads       | Total G4 intersects  |Intersects / Reads |
| ------------- |:-------------:| :-----:|-----:| 
|SRX2346892.bed|90|190|2.1111111111111|


Adding the SRX-IDs' metadata, ```Output_sorted.txt``` looks like:
```sh
(base) Oth.ALL.05.EPAS1.AllCell.bed>head -10 Output_sorted.txt 
SRX2346892.bed	90	190	2.111111111111111	Cardiovascular	EPAS1_ChIPSeq
SRX968419.bed	202	314	1.5544554455445545	Kidney		ChIP-Seq of HIF-2a in 786-O with HIF-1a re-expression
SRX968415.bed	277	297	1.0722021660649819	Kidney		ChIP-Seq of HIF-2a in 786-O Vector Alone
SRX4802363.bed	2208	1548	0.7010869565217391	liver		HepG2_16hrs 0.5% O2_HIF-2a (PM9)_Rep 1
SRX3346354.bed	1432	950	0.6634078212290503	liver		ChIP-Seq_786-O_HIF2A
SRX4802309.bed	279	185	0.6630824372759857	Kidney		HKC8_0.5_16hr_HIF2a_rep2
SRX4802349.bed	1663	928	0.5580276608538786	Kidney		RCC4_Normoxia_HIF-2a (PM9)_Rep 1
SRX4802308.bed	204	112	0.5490196078431373	Kidney		HKC8_0.5_16hr_HIF2a_rep1
SRX4802350.bed	1510	827	0.547682119205298	Kidney		RCC4_Normoxia_HIF-2a (PM9)_Rep 2
SRX4802364.bed	1411	744	0.5272856130403969	Kidney		HepG2_16hrs 0.5% O2_HIF-2a (PM9)_Rep 2
SRX4802322.bed	237	105	0.4430379746835443	Kidney
SRX4802344.bed	5226	2266	0.43360122464600076	Kidney	
SRX4802316.bed	283	121	0.4275618374558304	Kidney
SRX968422.bed	3188	1355	0.42503136762860727	Kidney
SRX4802337.bed	4803	2027	0.4220278992296481	Kidney
SRX4802335.bed	3289	1288	0.3916083916083916	Kidney
SRX4802342.bed	6696	2521	0.37649342891278376	Kidney
SRX2584128.bed	1442	491	0.340499306518724	Epidermis
SRX3051209.bed	3335	1090	0.3268365817091454	Kidney
SRX3346350.bed	748	208	0.27807486631016043	Kidney
SRX212365.bed	1120	235	0.20982142857142858	Macrophage
SRX359885.bed	816	169	0.20710784313725492	Digestive tract
SRX212355.bed	1208	240	0.1986754966887417	Macrophages
SRX212363.bed	1154	228	0.1975736568457539	Macrophages
SRX212364.bed	1114	220	0.19748653500897667	Macrophages
SRX212357.bed	1273	222	0.1743912018853103	Macrophages
SRX114493.bed	895	153	0.17094972067039105	Kidney
SRX212358.bed	1383	232	0.16775126536514823	Macrophages
SRX212356.bed	1263	206	0.16310372129849565	Macrophages
SRX4096727.bed	3870	484	0.12506459948320414	Prostate
SRX4802315.bed	611	64	0.10474631751227496 	Kidney
SRX4802302.bed	666	60	0.09009009009009009	Kidney
SRX4802323.bed	637	44	0.06907378335949764	Kidney	
SRX4802301.bed	575	39	0.06782608695652174	Kidney
SRX212366.bed	1012	51	0.05039525691699605	Macrophages
```



