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
    echo " "
}

# process input options
while getopts ":hf:s:" option
do
    case $option in 
        
        h) #displays help
            help
            exit;;

        f) #enter a file name
            file=$OPTARG;;

        s) #enter if samtools needs to sort or not
            sort=$OPTARG;;

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


if [ $sort ]; then
    module load samtools/1.5
    samtools view -S -b $file > inter.bam
    rm $file
    samtools sort inter.bam -o sorted.bam
    samtools view sorted.bam > $file
fi

python dedup.py -f $file
