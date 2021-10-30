import argparse
import random
import numpy as np
import re


# Import in the sam file that is sorted
def get_args():
    parser = argparse.ArgumentParser(description='Pass in the sorted SAM file, if paired reads, and a file containing the umis')
    parser.add_argument('-f', '--file', help='Upload a sorted by chromosome and position SAM file', required=True, type=str)
    parser.add_argument('-p', '--paired', help='Pass True or False for if they are paired or not. Default=False', default=False, type=bool)
    parser.add_argument('-u', '--umi', help='Specify path to UMI file that is separated by newlines. If no UMI file is given then \
        default=random and the program will assume that UMIs are unknown and will generate their own from the reads.', default='random', type=str)
    parser.add_argument('-o', '--output', help='Specify the output file name, it will habe the .sam postfix', required=True, type=str)
    parser.add_argument('-ds', '--store_duplicates', help='Specify if you would like duplicates returned to a separate file. \
        Set to True. Default=False. If set True, you must specify output file name with argument -do', default=False, type=bool)
    parser.add_argument('-do', '--duplicate_output', help='Specify the file name that the duplicates should be written to.', default=None, type=str)
    parser.add_argument('-q', '--quality', help='Set to True if you want to have the best quality read returned instead \
        of the first duplicate. Default=False', default=False, type=bool)

    return parser.parse_args()

# Load argparse info into a variable
args = get_args()

# Check ds and do arguments
if args.store_duplicates == True:
    if args.duplicate_output == None:
        raise ValueError('You did not specify an output for the duplicates. Please use -do and specify and output')

################################################### Function Section

def get_important_information(read_file_line: str): 
    """
    This function takes in a SAM file line and spits out the important information for that is needed so you can check
    to see if that read is a duplicate or not.
    """
    # Split the line into a list
    full_line = read_file_line.strip('\n').split('\t')

    # store important variables that are needed to check
    qname, flag, rname, pos, cigar, quality_score = full_line[0], full_line[1], full_line[2], full_line[3], full_line[5], full_line[10]

    # pull out the umi from the qname
    umi_qname = qname.split(':')[-1]

    # return the strand from the flag
    strand = check_strand(flag)
    
    return read_file_line, umi_qname, strand, rname, pos, cigar, quality_score

def add_cigar_to_pos(cigar_variable, position, strand):
    """
    Add the soft clipped ends to the beginning of the start position. This will be used for comparison with the positions in the 
    UMI dict
    If the strand is forward then it would subtract the soft clipped end to the 5' of the forward start position.
    If the strand is reverse then it would add the soft clipped end to the 5' of the reverse start position.
    """
    # Check to see if there is an S in the cigar_variable
    soft_clipping = 'S' in cigar_variable
    
    # If there is soft clipping then do this
    if soft_clipping:

        # First break up the line:
        split_cigar = re.findall('\d+[A-Z]{1}', cigar_variable)

        # Then we check the orientation of the read
        if strand == 'forward':
            
            # If S is in the first position then we need to subtract it
            if 'S' in split_cigar[0]:

                # first get the amount that was soft clipped off
                cigar_clipped = cigar_variable.split('S')[0]

                # Now calculate the new position for the forward strand
                new_pos = int(position) - int(cigar_clipped)

            else:
                # Just return position as we don't care about soft clipping at the end of the read for forward strands
                new_pos = position
        # This means it is the reverse strand so we add the soft clipped ends back on
        else:
            
            # instantiate a sum
            sum = 0

            # First we have to remove the first S if it is in it.
            if 'S' in split_cigar[0]:
                # Remove that item the S
                split_cigar.remove(split_cigar[0])

                # Iter through the elements in the list and add them all up so I know how much to add
                for elem in split_cigar:
                    if 'I' in elem:
                        split_cigar.remove(split_cigar.index(elem))
                    else:
                        sum += int(re.sub('[A-Z]', '', elem))
            
            # If the first elem isnt an S then we just go through the same logic as above.
            else:
                for elem in split_cigar:
                    if 'I' in elem:
                        split_cigar.remove(split_cigar.index(elem))
                    else:
                        sum += int(re.sub('[A-Z]', '', elem))
            

            new_pos = int(position) - sum

    
    # If there is no soft clipping present then check what strand
    else:

         # First break up the line:
        split_cigar = re.findall('\d+[A-Z]{1}', cigar_variable)

        # If forward spit out the same position
        if strand == 'forward':
        
            new_pos = int(position)
        # If reverse and no S then just add 101 to the POS to get 5' since POS is left most.
        else:
            sum = 0

            if 'S' in split_cigar[0]:
                split_cigar.remove(0)
                for elem in split_cigar:
                    if 'I' in elem:
                        split_cigar.remove(split_cigar.index(elem))
                    else:
                        sum += int(re.sub('[A-Z]', '', elem))
            else:
                for elem in split_cigar:
                    if 'I' in elem:
                        split_cigar.remove(split_cigar.index(elem))
                    else:
                        sum += int(re.sub('[A-Z]', '', elem))
            

            new_pos = int(position) - sum

    
    return new_pos

def check_strand(bitflag):
    """Check the bitflag to know if the strand is reversed or not"""
    if ((int(bitflag) & 16) == 16):
        return 'reverse'
    else:
        return 'forward'

def instantiate_umi_set(umi_file: str) -> set:
    """
    Creates a set of unique umis that will be used for reference to throw out bad UMI reads. If no known UMIs then we will create
    not create this
    """
    # Create the dictionary that will be used
    temp_set = set()

    # If umi file isn't set/ no know UMIs
    if umi_file == 'random':
        
        return temp_set

    else:
        # UMi file open and store as a variable
        with open(umi_file, 'r') as file:

            # first read in the lines
            for umi_line in file:
                
                # Make the umis from the file the keys for this dict
                temp_set.add(umi_line.strip('\n'))

    return temp_set

def convert_phred(asci_char_str):
    """Convert the ASCII characters to phred score 33"""
    return [ord(x) + 33 for x in asci_char_str]


def store_or_check_read_against_dict(read_line, storing_dict, umi_set, output_file, duplicate_file):
    """Take in a read and checks to see if this read already exists in the dictionary. If it does not then it is stored into the dictionary.
    This will be the main chunk of code that will be running for our analysis.
    """

    # check storing dict length
    dict_length = len(storing_dict)

    # First we break the line into multiple parts
    full_line, umi_qname, strand, rname, pos, cigar, quality_score = get_important_information(read_line)

    ################# We need to check what chromosome we are on, if it is a new one then write all reads to the file and clear the dictionary.

    # if dict_length > 1 check last added chrom against the currently pulled one.
    if dict_length > 1:
        # Pull out the last entry in the dict
        last_entry = list(storing_dict.keys())[-1]
        # Check to see if we are still on the same chromosome
        if storing_dict[last_entry][2] != rname:
            # Iter through the dict and write them to the output file
            for key in storing_dict:
                output_file.write(storing_dict[key][0])
            # Clear the dictionary as we are on a new chromosome now
            storing_dict.clear()
        else:
            pass
    else:
        pass


################# Dictionary is sorted now we can check/store new entries into/against the dict

    # First we check to see if the umi has any values set
    if args.umi == 'random':
        
        # Look for Ns in the UMI
        if 'N' not in umi_qname:

            # Update the position if it needs to be changed from soft clipping
            updated_pos = add_cigar_to_pos(cigar, pos, strand)

            # Create the key needed to be checked and or put into the dict
            key_string = umi_qname + "-" + strand + "-" + str(updated_pos) + "-" + rname

            # Check to see if it is in the dict
            if key_string in storing_dict:
                
                # We need to check to see if the user wants higher quality scores
                if args.quality == True:
                    # Check the quality scores against each other
                    if np.mean(convert_phred(quality_score)) > np.mean(convert_phred(storing_dict[key_string][1])):

                        # check duplicate store option
                        if args.store_duplicates == True:
                            duplicate_file.write(storing_dict[key_string][0])
                        else:
                            None

                        # Replace the read with the better quality one.
                        storing_dict[key_string] = (full_line, quality_score, rname)

                    # The quality is not better so do not overwrite this dict entry
                    else:
                        # If the user wants duplicates stored then write to duplicate file
                        if args.store_duplicates == True:
                            duplicate_file.write(full_line)
                            return None
                        else:
                            return None
                
                # If not then don't do anything as we already have this entry present, but check stored_duplicate arg first
                else:
                    if args.store_duplicates == True:
                        duplicate_file.write(full_line)
                        return None
                    else:
                        return None
            
            # There is no entry in the dict yet, this is a new read. We have to add it.        
            else:
                storing_dict[key_string] = (full_line, quality_score, rname)
                return None

        # If it isn't then throw it out.
        else:
            return None


    # We have known UMIs we will follow this logic 
    else:
        # Check to see if the umi is in the set
        if umi_qname in umi_set:
            
            # Update the position if it needs to be changed from soft clipping
            updated_pos = add_cigar_to_pos(cigar, pos, strand)

            # Create the key needed to be checked and or put into the dict
            key_string = umi_qname + "-" + strand + "-" + str(updated_pos) + "-" + rname

            # Check to see if it is in the dict
            if key_string in storing_dict:
                
                # We need to check to see if the user wants higher quality scores
                if args.quality == True:
                    # Check the quality scores against each other
                    if np.mean(convert_phred(quality_score)) > np.mean(convert_phred(storing_dict[key_string][1])):

                        # check duplicate store option
                        if args.store_duplicates == True:
                            duplicate_file.write(storing_dict[key_string][0])
                        else:
                            None

                        # Replace the read with the better quality one.
                        storing_dict[key_string] = (full_line, quality_score, rname)

                    # The quality is not better so do not overwrite this dict entry
                    else:
                        # If the user wants duplicates stored then write to duplicate file
                        if args.store_duplicates == True:
                            duplicate_file.write(full_line)
                            return None
                        else:
                            return None
                
                # If not then don't do anything as we already have this entry present, but check stored_duplicate arg first
                else:
                    if args.store_duplicates == True:
                        duplicate_file.write(full_line)
                        return None
                    else:
                        return None
            
            # There is no entry in the dict yet, this is a new read. We have to add it.        
            else:
                storing_dict[key_string] = (full_line, quality_score, rname)
                return None

        # If it isn't then throw it out.
        else:
            return None

########################################## Script Logic

##### Create global variables

# Create a random run number for the ouptut
run_id = random.randint(1, 100000)

# If -ds is flagged True then create this file.
if args.store_duplicates == True:
    duplicate_file = open(f'output/duplicates/{args.duplicate_output}_{run_id}.sam',
    'w')
else:
    duplicate_file = None

# Set past_chrome variable 
past_chrome = 1

# Create UMI dictionary
umi_created_set = instantiate_umi_set(args.umi)

# Storing dict
read_dict = dict()



###### Run script

# Create the output_file to write for
output_file = open(f'output/{args.output}_{run_id}_deduped.sam', 'w')

# Load in the file
with open(args.file, 'r') as sam_file:
    # Iter through each line
    for sam_line in sam_file:
        
        # Read in lines that are only read lines
        if r'@' not in sam_line:
            
            # Run operation function
            store_or_check_read_against_dict(sam_line, read_dict, umi_created_set, output_file, duplicate_file)

        
        # Write all of the header lines to the output file already
        else:
            output_file.write(sam_line)


# Add the last 10 lines to the output, since the script doesn't trigger it
for read in read_dict:
    output_file.write(read_dict[read][0])

# Close the output file
output_file.close()

if args.store_duplicates == True:
    duplicate_file.close()


