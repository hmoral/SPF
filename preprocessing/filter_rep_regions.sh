

## from gff to Bed file
sed 1d D2102046629.denovo.RepeatMasker.gff | awk '{print $1":"$4"-"$5}' > denovo_repetitive_regions.bed

sed -e s/:/\\t/g -e s/-/\\t/g denovo_repetitive_regions.bed > denovo_repetitive_regions.tab.bed

sed 1d D2102046629.known.RepeatMasker.gff | awk '{print $1":"$4"-"$5}' > known_repetitive_regions.bed

sed -e s/:/\\t/g -e s/-/\\t/g known_repetitive_regions.bed > known_repetitive_regions.tab.bed

## Merge both repetitive regions

bedtools sort -i <(cat known_repetitive_regions.tab.bed denovo_repetitive_regions.tab.bed) > repetitive_regions.both.tab.bed

bedtools merge -i repetitive_regions.both.tab.bed > repetitive_regions.both.merged.tab.bed

## To remove known repetitive regions from fasta file

module load bedtools

# Remove known regions from file
bedtools subtract -a Terpsiphone_corvina_B10K.tab.bed -b known_repetitive_regions.tab.bed > Nonrepetitive_Regions.known.bed

# Remove denovo regions from file
bedtools subtract -a Terpsiphone_corvina_B10K.tab.bed -b denovo_repetitive_regions.tab.bed > Nonrepetitive_Regions.denovo.bed

# Get bed file without known or de novo regions
bedtools subtract -a Nonrepetitive_Regions.known.bed -b denovo_repetitive_regions.tab.bed > Nonrepetitive_Regions.both.bed
