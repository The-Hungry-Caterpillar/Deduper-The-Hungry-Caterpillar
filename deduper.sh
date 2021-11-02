#!/bin/bash
#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=deduper
#SBATCH --output=deduper_%j.out
#SBATCH --error=deduper_%j.out
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# help functions, displays user options
help()
{
    echo "This program takes single-end sorted or unsorted sam files and removes PCR duplicates"
    echo " "
    echo "The following options are required:"
    echo "     -f     Input samfile name."
    echo " "
    echo "The following options are optional:"
    echo "     -s     Input any value if the input samfile is not already sorted. YOU MUST SORT THE SAMFILE IF NOT ALREADY SORTED!"
    echo "     -u     Input a text file of Unique Molecular Identifiers (UMIs)"
    echo "     -p     Input any value if you have paired end data."
    echo " "
}

# process input options
while getopts ":hf:s:u:p:" option
do
    case $option in 
        
        h) #displays help
            help
            exit;;

        f) #enter a file name
            file=$OPTARG;;

        s) #enter if samtools needs to sort or not
            sort=$OPTARG;;

        u) #enter UMI list
            umi=$OPTARG;;

        p) #paired end? too bad
            paired=$OPTARG;;

        \?) #displays invalid option
            echo "Error: Invalid option(s)"
            exit;;

    esac
done

#checks if -f was passed with a file, if not then exits script.
if ! [ -f "$file" ]; then
    echo "Error: You must pass a valid file to the -f option; see -h for help."
    exit 1
fi

if [ $paired ]; then
    echo "Sorry, we don't do paired end reads here."
fi

head -1000 sam.sam | grep '^@' > headers.txt

if [ $sort ]; then
    # module load samtools/1.5
    conda activate dedup
    samtools view -S -b $file > inter.bam
    rm $file
    samtools sort inter.bam -o sorted.bam
    rm inter.bam
    samtools view sorted.bam > $file
    rm sorted.bam
fi

if [ $umi ]; then
    python ogata_deduper.py -f $file -u $umi
else
    python ogata_deduper.py -f $file -u "NULL"
fi

rm headers.txt
