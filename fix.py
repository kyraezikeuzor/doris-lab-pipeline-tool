# Antibuddy: Doris Lab IG Pipeline Tool
import subprocess
import os
from spinner import Spinner
from utils import check_tool_availability
from utils import create_file

'''
Requirements
- samtools-1.9
- seqtk
- fatotwobit
- blat
- blatSrc
- Folder with a bam file inside
- Strain specific target file
'''

# Global Variables
strain_specific_target_options = []
identifierStarter = "m84248"
IDENTIFIER_LIST = []

# Add identifiers from each 0000X.<strain-specific-target> file to the complete file
def append_list_of_identifiers(filename, identifiers):
    # Adding copied list of items from each mini ighg1 file to a giant IGHh1 file
    with open(filename, 'a') as file:  # Open the file in append mode
        for identifier in identifiers:
            file.write(identifier + '\n')  # Write each identifier on a new line
    
    print("\n✅ Added identifiers to " + filename + '\n')

# Capture global faToTwoBit script path
def set_file_header(header):
    global file_header
    file_header = header

# Function to extract identifiers
def extract_identifiers(filename):
    
    def get_sequence_score(line):
        if line.startswith(identifierStarter):
            recording = False
            startIndex = 0
            endIndex = 0
            sequenceScore = ''

            for a in range(1, len(line)-1):
                if line[a-1] == " " and line[a] != " ":
                    recording = True
                    startIndex = a
                    break

            if recording == True:
                for b in range(startIndex, len(line)-1):
                    if line[b+1] == " ":
                        endIndex = b
                        break

            for c in range(startIndex, endIndex+1):
                sequenceScore = sequenceScore + line[c]
 
            return int(sequenceScore)

    def get_max_sequence_score(filename):
        sequenceScores = []

        with open(filename, 'r') as file:
            for line in file:
                if line.startswith(identifierStarter):
                    sequenceScore = get_sequence_score(line)         

                    sequenceScores.append(int(sequenceScore))

        if len(sequenceScores) > 0:
            sequenceScores.sort()

        return sequenceScores[-1]

    # Define necessary variables
    maxSequenceScore = get_max_sequence_score(filename)
    identifiers = []

    # Creating identifiers array and adding only sequences from each line
    with open(filename, 'r') as file:
        recording = False
        for line in file:
            sequence = ''
            sequenceScoreForLine = get_sequence_score(line)

            if line.startswith(identifierStarter) and sequenceScoreForLine == maxSequenceScore:
                for i in range(0, len(line)-1):
                    if line[i+1] != " ":
                        sequence = sequence + line[i]
                    elif line[i+1] == " ":
                        sequence = sequence + line[i]
                        break
                    
                identifiers.append(sequence)
    
    # Print the command output
    print(f"✅ Extracted sequence identifiers from {filename} successfully.")
    return identifiers


# Use fasta identifiers to retrieve sequences in the original run that match them
def match_sequences(in_file, identifiers_file, out_file):
    # Command Usage: seqtk subseq STC654_Bio_Sample_8_flnc.fasta STC654_Bio_Sample_8_IGHG1.txt > IGHG1_seqs.fasta
    
    # Load spinner
    spinner = Spinner("Using fasta identifiers to retrieve matching sequences... ", speed=0.1)
    spinner.start()

    try:
        # Construct the command to be executed
        command = ["seqtk", "subseq", in_file, identifiers_file]
        
        # Run the command and redirect stdout to the output file
        with open(out_file, 'w') as f:
            result = subprocess.run(command, check=True, stdout=f, stderr=subprocess.PIPE)
    
        # Print the command output
        print("\n✅ Matched sequences successfully.")

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "❌ Error running Seqtk command:" + e.stderr.decode() + '\n')
    
    finally:
        # Stop spinner
        spinner.stop()
    

# Main processing function
def process_files_in_folder(folder_path, fa_file):
    global IDENTIFIER_LIST

    for filename in os.listdir(folder_path):

        # Skip hidden files like .DS_Store
        if filename.startswith('.'):
            continue
        
        file_path = os.path.join(folder_path, filename)

        print(filename)
        
        file_path = os.path.join(folder_path, filename)
        
        print(file_path)

        if os.path.isfile(file_path):
            # Extract the header (part before the period)
            header = filename.split('.')[0]
            set_file_header(header)
            
            # Extract identifiers and append them to the global list
            identifiers = extract_identifiers(file_path)
            IDENTIFIER_LIST.extend(identifiers)
            
    # Create a new file for each processed file
    name_of_combined_file = f"{file_header}.COMBINED.txt"
    full_combined_filename = os.path.join(folder_path, name_of_combined_file)

    create_file(full_combined_filename, "")

    append_list_of_identifiers(full_combined_filename, IDENTIFIER_LIST)

    # Create seqs file
    name_of_seqs_file = f"{file_header}-seqs.fasta"
    full_seqs_filename = os.path.join(folder_path, name_of_seqs_file)

    # Command Usage: seqtk subseq STC654_Bio_Sample_8_flnc.fasta STC654_Bio_Sample_8_IGHG1.txt > IGHG1_seqs.fasta
    match_sequences(fa_file, full_combined_filename, full_seqs_filename)




items_list = [
    {
        "file": "doris-lab/bio-sample-4/STC654-Bio-Sample-4-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-4/IGHA_A3",
            "doris-lab/bio-sample-4/IGHD1_A3",
            "doris-lab/bio-sample-4/IGHD2_A3",
            "doris-lab/bio-sample-4/IGHG1_A3",
            "doris-lab/bio-sample-4/IGHG2A_A3",
            "doris-lab/bio-sample-4/IGHG2B_A3",
            "doris-lab/bio-sample-4/IGHG2C_A3",
            "doris-lab/bio-sample-4/IGHG2C_A3long",
            "doris-lab/bio-sample-4/IGHA_B2",
            "doris-lab/bio-sample-4/IGHD_B2",
            "doris-lab/bio-sample-4/IGHG1_B2",
            "doris-lab/bio-sample-4/IGHG2A_B2",
            "doris-lab/bio-sample-4/IGHG2B_B2",
            "doris-lab/bio-sample-4/IGHG2C_B2",
            "doris-lab/bio-sample-4/IGHG2C_B2long"
        ]
    },
    {
        "file": "doris-lab/bio-sample-5/STC654-Bio-Sample-5-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-5/IGHA_A3",
            "doris-lab/bio-sample-5/IGHD1_A3",
            "doris-lab/bio-sample-5/IGHD2_A3",
            "doris-lab/bio-sample-5/IGHG1_A3",
            "doris-lab/bio-sample-5/IGHG2A_A3",
            "doris-lab/bio-sample-5/IGHG2B_A3",
            "doris-lab/bio-sample-5/IGHG2C_A3",
            "doris-lab/bio-sample-5/IGHG2C_A3long",
            "doris-lab/bio-sample-5/IGHA_B2",
            "doris-lab/bio-sample-5/IGHD_B2",
            "doris-lab/bio-sample-5/IGHG1_B2",
            "doris-lab/bio-sample-5/IGHG2A_B2",
            "doris-lab/bio-sample-5/IGHG2B_B2",
            "doris-lab/bio-sample-5/IGHG2C_B2",
            "doris-lab/bio-sample-5/IGHG2C_B2long"
        ]
    },
    {
        "file": "doris-lab/bio-sample-6/STC654-Bio-Sample-6-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-6/IGHA_A3",
            "doris-lab/bio-sample-6/IGHD1_A3",
            "doris-lab/bio-sample-6/IGHD2_A3",
            "doris-lab/bio-sample-6/IGHG1_A3",
            "doris-lab/bio-sample-6/IGHG2A_A3",
            "doris-lab/bio-sample-6/IGHG2B_A3",
            "doris-lab/bio-sample-6/IGHG2C_A3",
            "doris-lab/bio-sample-6/IGHG2C_A3long",
            "doris-lab/bio-sample-6/IGHA_B2",
            "doris-lab/bio-sample-6/IGHD_B2",
            "doris-lab/bio-sample-6/IGHG1_B2",
            "doris-lab/bio-sample-6/IGHG2A_B2",
            "doris-lab/bio-sample-6/IGHG2B_B2",
            "doris-lab/bio-sample-6/IGHG2C_B2",
            "doris-lab/bio-sample-6/IGHG2C_B2long"
        ]
    },
    {
        "file": "doris-lab/bio-sample-7/STC654-Bio-Sample-7-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-7/IGHA_A3",
            "doris-lab/bio-sample-7/IGHD1_A3",
            "doris-lab/bio-sample-7/IGHD2_A3",
            "doris-lab/bio-sample-7/IGHG1_A3",
            "doris-lab/bio-sample-7/IGHG2A_A3",
            "doris-lab/bio-sample-7/IGHG2B_A3",
            "doris-lab/bio-sample-7/IGHG2C_A3",
            "doris-lab/bio-sample-7/IGHG2C_A3long",
            "doris-lab/bio-sample-7/IGHA_B2",
            "doris-lab/bio-sample-7/IGHD_B2",
            "doris-lab/bio-sample-7/IGHG1_B2",
            "doris-lab/bio-sample-7/IGHG2A_B2",
            "doris-lab/bio-sample-7/IGHG2B_B2",
            "doris-lab/bio-sample-7/IGHG2C_B2",
            "doris-lab/bio-sample-7/IGHG2C_B2long"
        ]
    },
    {
        "file": "doris-lab/bio-sample-8/STC654-Bio-Sample-8-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-8/IGHA_A3",
            "doris-lab/bio-sample-8/IGHD1_A3",
            "doris-lab/bio-sample-8/IGHD2_A3",
            "doris-lab/bio-sample-8/IGHG1_A3",
            "doris-lab/bio-sample-8/IGHG2A_A3",
            "doris-lab/bio-sample-8/IGHG2B_A3",
            "doris-lab/bio-sample-8/IGHG2C_A3",
            "doris-lab/bio-sample-8/IGHG2C_A3long",
            "doris-lab/bio-sample-8/IGHA_B2",
            "doris-lab/bio-sample-8/IGHD_B2",
            "doris-lab/bio-sample-8/IGHG1_B2",
            "doris-lab/bio-sample-8/IGHG2A_B2",
            "doris-lab/bio-sample-8/IGHG2B_B2",
            "doris-lab/bio-sample-8/IGHG2C_B2",
            "doris-lab/bio-sample-8/IGHG2C_B2long"
        ]
    }
];

items_list = [
    {
        "file": "doris-lab/bio-sample-1/STC654-Bio-Sample-1-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-1/igha_a3",
        ]
    },
        {
        "file": "doris-lab/bio-sample-4/STC654-Bio-Sample-4-flnc.fa",
        "folders": [
            "doris-lab/bio-sample-4/IGHA_A3",
        ]
    }
]


def process(folder, fa_file):
    folder_path = folder
    fa_filename = fa_file

    # Run the processing function
    process_files_in_folder(folder_path, fa_filename)


def main():
    for item in items_list:
        for folder in item['folders']:
            process(folder, item['file'])


# Call main function
if __name__ == "__main__":
    main()