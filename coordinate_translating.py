import argparse
import re


read_consuming_ops = ('M', 'I', 'X')
ref_consuming_ops = ('M', 'D', 'N', 'X')


def get_genomic_position(transcript_info, transcript_pos):
    """ Convert from a transcript coordinate to a corresponding genomic coordinate"""

    # Breakdown the CIGAR string into each component part, each segment consisting of a cigar length and cigar type
    parsed_cigar = re.findall(r'(\d+)([MXDIN])', transcript_info['cigar'])
    # Use the starting point of the transcript as the initial genomic position
    genomic_pos = transcript_info['start']
    # Keep track of how many coordinates within the transcript are accounted for
    accounted_transcript_length = 0

    # Iterate through the cigar string to keep incrementing the genomic position until the requested transcript position is reached
    for cigar_length, cigar_type in parsed_cigar:
        cigar_length = int(cigar_length)

        if cigar_type in ('M', 'X'):
            # For either matches ('M') or mismatches ('X'), increase the genomic position by the entire cigar length if the
            # requested transcript position is not yet reached, otherwise only increase by the remaining number of digits left
            # until the requested transcript position will be reached
            genomic_pos += min(cigar_length, transcript_pos - accounted_transcript_length)
            # The accounted transcript length can be incremented by the entire cigar length. Because the accounted transcript length
            # is only needed for knowing when to break out of the loop, it is okay to simply increment by entire cigar length even
            # if this could be overkill if the requested transcript position is before the end of this cigar segment
            accounted_transcript_length += cigar_length
        # For deletions ('D') or gaps ('N'), increase the genomic position without altering the transcript position
        elif cigar_type in ('D', 'N'):
            genomic_pos += cigar_length
        # For insertions ('I'), increase the accounted transcript length without altering the genomic position
        elif cigar_type == 'I':
            accounted_transcript_length += cigar_length

        # Once the requested transcript position is reached, there is no need to parse the remainder of the CIGAR string
        if transcript_pos < accounted_transcript_length:
            break
    return genomic_pos


def translate_coordinates(transcripts_input_p, queries_input_p, output_file_p):
    """ Given 2 input files, one containing information about various transcripts and another containing information
    regarding requested coordinates from a transcript, output the set of translated coordinates"""

    # Collect info from transcripts input into a dictionary of each transcript
    transcripts = {}
    for row in open(transcripts_input_p):
        cols = row.strip().split('\t')
        if len(cols) != 4:
            raise Exception('Invalid transcripts input provided. Expecting 4 columns, received {} columns.'.format(len(cols)))
        transcript_name, chrom, start, cigar = cols
        try:
            start = int(start)
        except ValueError:
            raise ValueError('Invalid start position of `{}` provided, must be an integer'.format(start))

        if transcript_name in transcripts:
            raise Exception('Duplicate transcripts provided, each transcript should be unique')

        transcripts[transcript_name] = {'chrom': chrom, 'start': start, 'cigar': cigar}

    output_file = open(output_file_p, 'w')
    # Read through queries input to translate each transcript position into its corresponding genomic position
    for row in open(queries_input_p):
        cols = row.strip().split('\t')
        if len(cols) != 2:
            raise Exception('Invalid queries input provided. Expecting 2 columns, received {} columns.'.format(len(cols)))
        transcript_name, transcript_pos = cols
        try:
            transcript_pos = int(transcript_pos)
        except ValueError:
            raise ValueError('Invalid transcript_pos position of `{}` provided, must be an integer'.format(transcript_pos))

        transcript_info = transcripts.get(transcript_name)
        if not transcript_info:
            raise Exception('Unknown transcript provided in queries input: `{}`'.format(transcript_name))

        # Retrieve the genomic position from the requested transcript position
        genomic_pos = get_genomic_position(transcript_info, transcript_pos)

        # Combine all of the outputs together and write to the output file
        output = cols + [transcript_info['chrom'], str(genomic_pos)]
        output_file.write('\t'.join(output) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Translate from transcript coordinates to genomic coordinates')
    parser.add_argument('-t', '--transcripts_input_p', help='Path of the transcripts input file.')
    parser.add_argument('-q', '--queries_input_p', help='Path of the queries input file.')
    parser.add_argument('-o', '--output_file_p', help='Path of the output file.')
    args = parser.parse_args()
    translate_coordinates(args.transcripts_input_p, args.queries_input_p, args.output_file_p)
