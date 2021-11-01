import argparse
import random
import numpy as np


# Import in the sam file that is sorted
def get_args():
    parser = argparse.ArgumentParser(description='Pass in the sorted SAM file, if paired reads, and a file containing the umis')
    parser.add_argument('-f', '--file', help='Upload the sorted SAM file', required=True, type=str)
    parser.add_argument('-p', '--paired', help='Pass True or False for if they are paired or not. Default=False', default=False, type=bool)
    parser.add_argument('-u', '--umi', help='Specify path to UMI file that is separated by newlines.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Specify the output file name, it will habe the .sam postfix', required=True, type=str)
    parser.add_argument('-ds', '--store_duplicates', help='Specify if you would like duplicates returned to a separate file. \
        Type True or False. If set True, you must specify output file name with argument -do', default=False, type=bool)
    parser.add_argument('-do', '--duplicate_output', help='Specify the file to duplicate file name to write to.', default=None, type=str)

    return parser.parse_args()

# Load argparse info into a variable
args = get_args()

# Check ds and do arguments
if args.store_duplicates == True:
    if args.duplicate_ouptut == None:
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
        # first get the amount that was soft clipped off
        cigar_clipped = cigar_variable.split('S')[0]

        # Then we check the orientation of the read
        if strand == 'forward':
            # Now calculate the new position for the forward strand
            new_pos = int(position) - int(cigar_clipped)
        # This means it is the reverse strand so we add the soft clipped ends back on
        else:
            new_pos = int(position) + int(cigar_clipped)

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
        full_read = []
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
    return [ord(x) + 33 for x in asci_char]

def duplicates(lst, item):
    """
    Get the duplicate indexes in a list
    """
    return [i for i, x in enumerate(lst) if x == item]


def check_chromsome(current_chrom, chrom_in_sam, umi_dict, sliding_list, output):
    """If the chroms don't match then clear the lists in the dictionary as we have moved onto a new 
    chromosome and don't need to retain that info"""
    if current_chrom != chrom_in_sam:
        #print("Clearing Dict: Line 111")
        for key in umi_dict:
            umi_dict[key]['pos'].clear()
            umi_dict[key]['chromosome'].clear()
            umi_dict[key]['strand'].clear()
            umi_dict[key]['quality_score'].clear()

        # Write the rest of the sliding window reads to the file
        for read in sliding_list:
            output.write(read)

        # Clear the sliding_list:
        sliding_list.clear()

    return umi_dict, sliding_list
    

def quality_score_sliding_window(new_line, lst, output_write):
    """Store the most recent lines in memory so that they can be compared if a duplicate shows up.
    We can test the duplicate's quality score with the already stored read's quality score.
    As we move to the next read the oldest read is deleted from the list and written to the output file."""
    
    # First check length
    if len(lst) == 10:
        # Write oldest element to output file
        output_write.write(lst[0])
        # Remove just written element from list
        lst.pop(0)
        # Add the new line to the end of the list, making it the newest memeber
        lst.append(new_line)

    # If it isn't a length of ten then just add values    
    else:
        lst.append(new_line)

    return lst


def check_quality_score(dup_score, sliding_window, start_pos, full_dup_line):
    """There is a duplicate present. Check the duplicates quality vs the stored reads quality. If duplicate is better
    replace sliding windows entry with the duplicates."""

    # Find the same start position read in the sliding window
    index = 0
    for read in sliding_window:
        split_read = read.strip('\n').split('\t')
        
        # Check to see that the start positions match up 
        if int(start_pos) == int(add_cigar_to_pos(split_read[5], split_read[3], check_strand(split_read[1]))):
            
            # Check to see if the duplicates quality score is better or not.
            if dup_score > np.mean(convert_phred(read[10])):
                # Replace the read at that index with the better quality duplicate
                sliding_window[index] = full_dup_line

                # Return True bool for later processes
                return True
            # Go to the next read
            else:
                index+=1
        # Go to the next read
        else:
            index+=1
    # Return False bool for later processes
    return False




def check_duplicate(umi_qname_v, pos_v, strand_v, chrom_v, umi_dict_variable, current_chrom, window, out):
    """
    Take the important variables and check them against the umi dictionary. First we check that the umi is present, then you go a level deeper
    and look at the chromosome. Then strand, and then pos. 

    If it is not in the dictionary then it will add it to the dictionary. As that is not a duplicate and needs to be stored
    so we can check to see if any new lines look like that read.
    """
    # if chromosome is different then delete all the info within the dictionary since we are now moving onto a different chromosome
    umi_dict_actual, window = check_chromsome(current_chrom, chrom_v, umi_dict_variable, window, out)

    # First we need to check to see if the umi is within the dictionary keys
    check_umi = umi_qname_v in umi_dict_variable

    # Get the lower level dictionary of that umi.
    if check_umi:    
        lower_level = umi_dict_actual[umi_qname]

        # Check to see if the chrom is in the 'chromosome' key
        if chrom_v in lower_level['chromosome']:
            
            # Pull out all indexes where the matching chromosome is
            index_list = duplicates(lower_level['chromosome'], chrom_v)
            
            # Loop through index list
            for index in index_list:

                # pull out each pos and strand at that index
                pos_i = lower_level['pos'][index]
                strand_i = lower_level['strand'][index]
                
                # Check to see if the strands are the same
                if strand_i == strand_v:

                    # Check to see if the pos is the same
                    if pos_i == pos_v:
                        return True, index, window
                    
                    # Go to the next index if the position doesn't match
                    else:
                        next
                # Go to the next index if the strand doesn't match
                else:
                    next

    # umi isn't in the keys and doesn't need to be added to anything.
    else:
        return 'wrong', '_', window

    # If it gets through everything and it isn't within the dictionary then return False
    return False, '_', window

########################################## Script Logic

##### Create global variables

# Create a random run number for the ouptut
run_id = random.randint(1, 100000)

# If -ds is flagged True then create this file.
if args.store_duplicates == True:
    duplicate_file = open(f'output/duplicates/{args.duplicate_output}_{run_id}.sam',
    'w')

# Set past_chrome variable 
past_chrome = 1

# Create UMI dictionary
umi_dict = instantiate_duplicate_dictionary(args.umi)

# create sliding window list
sliding_window = list()


###### Run script


# Create the output_file to write for
output_file = open(f'output/{args.output}_{run_id}.sam', 'w')

# Load in the file
with open(args.file, 'r') as sam_file:
    # Iter through each line
    for sam_line in sam_file:
        
        # Read in lines that are only read lines
        if r'@' not in sam_line:
            # Store the variables
            full_line, umi_qname, strand, rname, pos, cigar, quality_score = get_important_information(sam_line)

            # Change the start position
            updated_pos = add_cigar_to_pos(cigar, pos)

            # This return a bool or a str result to check if the read is found within the dictionary.
            duplicate, dup_index, sliding_window = check_duplicate(umi_qname, updated_pos, strand, rname, umi_dict, past_chrome, sliding_window, output_file)
            #print(duplicate, "Line 199", sep='\t')
            
            # If duplicate is a string then the umi is not correct so just move on to the next read
            if type(duplicate) == str:
                #print('Throwing away crappy UMI read: line 202')
                past_chrome = rname
                next

            # If duplicate is True then move onto the next read
            elif duplicate:

                # We need to check the quality score of the read
                #stored_qual = umi_dict[umi_qname]['quality_score'][dup_index]
                dup_qual = np.mean(convert_phred(quality_score))

                qual_check = check_quality_score(dup_qual, sliding_window, updated_pos, full_line)

                #print('Moving on, this read is a duplicate: line 208')
                past_chrome = rname
                if args.store_duplicates == False:
                    next
                else:
                    duplicate_file.write(full_line)
            
            # If duplicate is False then it isn't a duplicate, so add to the dictionary
            else:
                #print('else: 214')
                umi_dict[umi_qname]['pos'].append(int(updated_pos))
                umi_dict[umi_qname]['chromosome'].append(rname)
                umi_dict[umi_qname]['strand'].append(strand)
                # Dont think I need to add this actually
                # umi_dict[umi_qname]['quality_score'].append(np.mean(convert_phred(quality_score)))

                
                # Also, write to a output file
                #print('Writing to File: 220')
                #output_file.write(full_line)
                sliding_window = quality_score_sliding_window(full_line, sliding_window, output_file)

                # Move to next line
                past_chrome = rname
                next
        
        # Write all of the header lines to the output file already
        else:
            output_file.write(sam_line)


# Add the last 10 lines to the output, since the script doesn't trigger it
for read in sliding_window:
    output_file.write(read)

# Close the output file
output_file.close()

if args.store_duplicates == True:
    duplicate_file.close()



# use umi-chr-pos-strand keys and just the read as a value
