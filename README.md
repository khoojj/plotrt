# plotrt

I written this as a first try in py script writing to automate some of our work processes.
The script plots the number of reads vs sequencing time with fastq produced from Guppy basecalling after ONT sequencing as input.
Biopython is required to run the script.

Usage example:

For naming individual fastq files:

plotrt.py -i *.fastq -o out_dir_name

For naming folder containing fastq files:

plotrt.py -i /home/username/directory -o out_dir_name

210901  Version 1.0 released
