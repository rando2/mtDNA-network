# mtDNA-network
Here is the code that can be used to generate the haplotype networks for the farm foxes.
I had to modify the package pegas to make it consider deletions as evolutionarily significant (this requires that the sequences be the same lenghth, though, which was not a limitation of the original code).
It takes as input a spreadsheet where the frequencies of each haplotype in each population are noted.
It also needs a .fas file with the alignment of each of these sequences (obtained from NCBI, except the novel haploptype discovered at the ICG).
