import argparse
from typing import TextIO

# Import in the sam file that is sorted
def get_args():
    parser = argparse.ArgumentParser(description='Pass in the sorted SAM file, if paired reads, and a file containing the umis')
    parser.add_argument('-f', '--file', help='Upload the sorted SAM file', required=True, type=str)
    parser.add_argument('-p', '--paired', help='Pass True or False for if they are paired or not. Default=False', required=True, default=False, type=bool)
    parser.add_argument('-u', '--umi', help='Specify path to UMI file that is separated by newlines.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Specify the output file')
    return parser.parse_args()

# Load argparse info into a variable
args = get_args()

################################################### Function Section

def get_important_information(read_file_line: str): 
    """
    This function takes in a SAM file line and spits out the important information for that is needed so you can check
    to see if that read is a duplicate or not.
    """
    # Split the line into a list
    full_line = read_file_line.strip('\n').split('\t')

    # store important variables that are needed to check
    qname, flag, rname, pos, cigar = full_line[0], full_line[1], full_line[2], full_line[3], full_line[5]

    # pull out the umi from the qname
    umi_qname = qname.split(':')[-1]

    # return the strand from the flag
    strand = check_strand(flag)
    
    return umi_qname, strand, rname, pos, cigar

def add_cigar_to_pos(cigar_variable, position):
    """
    Add the soft clipped ends to the beginning of the start position. This will be used for comparison with the positions in the 
    UMI dict
    After comparison it would remove them so it adds the original if the position isn't found in the dict.
    """
    # Check to see if there is an S in the cigar_variable
    soft_clipping = 'S' in cigar_variable
    # If there is soft clipping then do this
    if soft_clipping:
        # first get the amount that was soft clipped off
        cigar_clipped = cigar_variable.split('S')[0]
        # Now calculate the new position
        new_pos = int(position) - int(cigar_clipped)

    # If there is no soft clipping present then just spit back out the position
    else:
        new_pos = position 

    
    return new_pos

def check_strand(bitflag):
    """Check the bitflag to know if the strand is reversed or not"""
    if ((int(bitflag) & 16) == 16):
        return 'reverse'
    else:
        return 'forward'

def instantiate_duplicate_dictionary(umi_file: str) -> dict:
    """
    Create a dictionary that will use the UMI file and use the UMIs as keys and the value will be another dict that has pos, strand, chromosome, and BitFlag
    as the keys in the secondary dictionary. Then the values for those will be the ones that do not have duplicates up till that point.

    dict[UMI] = {
        pos = []
        strand = []
        chromosome = []
        quality_score = []
    }
    """
    # Create the dictionary that will be used
    temp_dict = dict()

    with open(umi_file, 'r') as file:

        # first read in the lines
        for umi_line in file:
            
            # Make the umis from the file the keys for this dict
            temp_dict[umi_line.strip('\n')] = {
                'pos': list(),
                'strand': list(),
                'chromosome': list(),
                'quality_score': list()
            }

    return temp_dict

def convert_phred(asci_char):
    """Convert the ASCII character to phred score 33"""
    return ord(asci_char) + 33

def duplicates(lst, item):
    """
    Get the duplicate indexes in a list
    """
    return [i for i, x in enumerate(lst) if x == item]


def check_duplicate(umi_qname, pos, strand, chrom, umi_dict):
    """
    Take the important variables and check them against the umi dictionary. First we check that the umi is present, then you go a level deeper
    and look at the chromosome. Then strand, and then pos. 

    If it is not in the dictionary then it will add it to the dictionary. As that is not a duplicate and needs to be stored
    so we can check to see if any new lines look like that read.
    """

    # First we need to check to see if the umi is within the dictionary keys
    check_umi = umi_qname in umi_dict
    
    # Get the lower level dictionary of that umi.
    if check_umi:    
        lower_level = umi_dict[umi_qname]

        # Check to see if the chrom is in the 'chromosome' key
        if chrom in lower_level['chromosome']:
            # Pull out all indexes where the matching chromosome is
            index_list = duplicates(lower_level['chromosome'], chrom)
            # Loop through index list
            for index in index_list:

                # pull out each pos and strand at that index
                pos_i = lower_level['pos'][index]
                strand_i = lower_level['strand'][index]
                
                # Check to see if the strands are the same
                if strand_i == strand:

                    # Check to see if the pos is the same
                    if pos_i == pos:
                        return True
                    
                    # Go to the next index if the position doesn't match
                    else:
                        next
                # Go to the next index if the strand doesn't match
                else:
                    next

    # umi isn't in the keys and doesn't need to be added to anything.
    else:
        return 'wrong'

    # If it gets through everything and it isn't within the dictionary then return False
    return False

########################################## Script Logic


# Create UMI dictionary
umi_dict = instantiate_duplicate_dictionary(args.umi)

output_file = open()

# Load in the file
with open(args.file, 'r') as sam_file:
    # Iter through each line
    for sam_line in sam_file:
        
        # Read in lines that are only read lines
        if r'@' not in sam_line:
            # Store the variables
            umi_qname, strand, rname, pos, cigar = get_important_information(sam_line)

            # Change the start position
            updated_pos = add_cigar_to_pos(cigar, pos)

            duplicate = check_duplicate(umi_qname, updated_pos, strand, rname, umi_dict)

            if duplicate:
                next
            else:
                umi_dict[umi_qname]['pos'].append(updated_pos)
                umi_dict[umi_qname]['chromosome'].append(rname)
                umi_dict[umi_qname]['strand'].append(strand)

                
                # Also, write to a file
                # f.wrtie(umid_dict) Blah
            print(duplicate)


            

