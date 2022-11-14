# accuragen_coding_challenge
This project is designed for completion of Accuragen's coding challenge.

The goal of the coding challenge is to translate between transcript and genomic coordinates.
Two input files are required.
Input file 1 contains a set of transcripts, each one consisting of a name, chromosome, start position, and CIGAR string.
Input file 2 contains a set of queries of specific transcript coordinates to be mapped to corresponding genomic coordinates.

The function produces an output file of each query provided in input file 2 along with the corresponding genomic coordinates.

There is an optional argument for using reverse orientation instead of forward orientation for the transcript coordinates.
The limitation for this orientation is that the start position for each transcript only corresponds to the forward orientation,
so there is no way of knowing what the start position would be for the reverse orientation.
