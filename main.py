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

# Get strain specific target options
def set_target_options():
    global strain_specific_target_options

    targets_file_folder = os.path.dirname(os.path.realpath(__file__))
    targets_file = os.path.join(targets_file_folder,'targets.txt')

    try:
        with open(targets_file, 'r') as file:
            strain_specific_target_options = [line.strip().upper() for line in file.readlines()]
    
    except FileNotFoundError:
        print("\n❌ The file targets.txt was not found. ")
    
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")

# Capture global faToTwoBit script path
def set_fatotwobit_script(script_path):
    global fatotwobit_script
    fatotwobit_script = script_path

# Capture global blat script path
def set_blat_script(script_path):
    global blat_script
    blat_script = script_path

# Convert bam file to fasta file
def convert_bam_to_fasta(in_file, out_file):
    # Command Usage: samtools fasta STC654_Bio_Sample_8_flnc.bam > STC654_Bio_Sample_8_flnc.fasta
    
    # Load spinner
    spinner = Spinner("Converting bam file to fasta file... ", speed=0.1)
    spinner.start()

    try:
        # Construct the command to be executed
        command = ["samtools", "fasta", "-0", out_file, in_file]
        
        # Run the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print the command output
        print('\n' + "✅ Converted bam to fasta file successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "❌ Error running Samtools command:" + e.stderr.decode() + '\n')
    
    finally:
        # Close spinner
        spinner.stop()
    
# Split fasta file into 10 mini fiels
def split_fasta_file(in_file, out_file):
    # Command Usage: seqkt split -n 10 STC654_Bio_Sample_8_flnc  STC654_Bio_Sample_8_flnc.fasta
    
    # Load spinner
    spinner = Spinner("Splitting fasta file into 10 chunks... ", speed=0.1)
    spinner.start()

    try:     
        # Construct the command to be executed
        command = ["seqtk", "split", "-n", "10", in_file, out_file]
        
        # Run the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print the command output
        print('\n' + "✅ Fasta files split successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "❌ Error running Seqtk command:" + e.stderr.decode() + '\n')

    finally:
        # Close spinner
        spinner.stop()

# Index each smaller fasta file with faToTwoBit
def index_fasta_file(in_file, out_file):
    # Command Usage: fatotwobit STC654_Bio_Sample_8_flnc.00001.fasta STC654_Bio_Sample_8_flnc.00001.2bit
    
    # Load spinner
    spinner = Spinner("Creating .2bit file from small fasta file... ", speed=0.1)
    spinner.start()

    # Create .2bit file
    create_file(out_file,'')

    try:           
        # Construct the command to be executed
        command = [fatotwobit_script, in_file, out_file]
        
        # Run the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print the command output
        print('\n' + "✅ .2bit from .fa file generated successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "\n❌ Error running faToTwoBit: " + e.stderr.decode() + '\n')

    finally:
        # Close spinner
        spinner.stop()

# Extract sequences of interest with blat
def extract_sequences_of_interest(database_file, query_file, output_file):
    # Command Usage: blat STC654_Bio_Sample_8_flnc.00001.2bit SHRA3_IGHG1.txt -out=blast STC654_Bio_Sample_8_IGHG1.00001.txt

    # Load spinner
    spinner = Spinner("Extracing sequences of interest with blat... ", speed=0.1)
    spinner.start()

    try:
        # Construct the command to be executed
        command = [blat_script, database_file, query_file, "-out=blast", output_file]
        
        # Run the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print the command output
        print('\n' + "✅ Extracted sequences with blat successfully.", result.stdout.decode())

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "❌ Error running BLAT command:" + e.stderr.decode() + '\n')

    finally:
        # Close spinner
        spinner.stop()

# Grab only the identifiers for reads
def extract_identifiers(filename):
    # This function goes into the strain target text file and only extracts the sequences with highest score

    # Define helper functions
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

    # Load spinner
    spinner = Spinner("Extracting sequences from identifer files... ", speed=0.1)
    spinner.start()

    # Creating identifiers array and adding only sequences from each line
    with open(filename, 'r') as file:
        recording = False
        for line in file:
            sequence = ''
            sequenceScoreForLine = get_sequence_score(line)

            if line.startswith("m84248") and sequenceScoreForLine == maxSequenceScore:
                for i in range(0, len(line)-1):
                    if line[i+1] != " ":
                        sequence = sequence + line[i]
                    elif line[i+1] == " ":
                        sequence = sequence + line[i]
                        break
                    
                identifiers.append(sequence)
    
    # Close spinner
    spinner.stop()

    # Print the command output
    print(f"✅ Extracted sequence idenfifiers from {filename} successfully.")

    return identifiers

# Add identifiers from each 0000X.<strain-specific-target> file to the complete file
def append_list_of_identifiers(filename, identifiers):
    # Adding copied list of items from each mini ighg1 file to a giant IGHh1 file
    with open(filename, 'a') as file:  # Open the file in append mode
        for identifier in identifiers:
            file.write(identifier + '\n')  # Write each identifier on a new line
    
    print("\n✅ Added identifiers to " + filename + '\n')

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
        
# Performs entire extraction process
def start():
    # --------------------------- STEP 1 PART 1: GET VARIABLES ---------------------------
    # --------------------------- GET BAM FILE ---------------------------
    
    # Get BAM file + Check if BAM file exists
    while True:
        step1_bam_file = input("STEP 1/5: Enter path to the .bam file: ")

        if os.path.exists(step1_bam_file):
            break
        else:
            print('\n' + f"‼️ Error: The file '{step1_bam_file}' does not exist." + '\n')
    
    print('\n')


    # --------------------------- GET STRAIN SPECIFIC TARGET ---------------------------
    # Get Strain Specific Target File + Check if Strain Specific Target Is one of our options
    while True: 
        step1_strain_specific_target = input(f"Options: {strain_specific_target_options} \nSTEP 2/5: Enter strain specific target (case-insensitive): ")
        step1_strain_specific_target_uppercase = step1_strain_specific_target.upper()

        if step1_strain_specific_target_uppercase in strain_specific_target_options:
            break
        else:
            print('\n' + f"‼️  {step1_strain_specific_target_uppercase} is not an available option. Options: ", strain_specific_target_options + '\n')

    print('\n')


    # --------------------------- GET STRAIN SPECIFIC TARGET FILE ---------------------------
    # Get Strain Specific Target File + Check if Strain Specific Target File Exists

    while True: 
        step1_strain_specific_target_file = input("STEP 3/5: Enter path to strain specific target file (ex: IGHG1_SHRA3.fasta): ")

        if os.path.exists(step1_strain_specific_target_file):
            break
        else:
            print('\n' + f"‼️ Error: The file '{step1_strain_specific_target_file}' does not exist." + '\n')

    print('\n')


    # --------------------------- GET BLAT SCRIPT PATH ---------------------------
    # Check if blat is available
    if not check_tool_availability("blat"):
        print("STEP 4/5: blat was not found systemwide.")
        
        while True:
            blat_script_path = input("STEP 4/5: Enter path to your blat script (e.g., /condaenv/bin/blat): ")    
            
            if os.path.exists(blat_script_path):
                set_blat_script(blat_script_path)
                break
            else:
                print('\n' + "‼️ Error: The provided path does not exist. Please ensure the path is correct and try again.")
    else:
        print("STEP 4/5: blat was found systemwide.")
    
    print('\n')


    # --------------------------- GET FATOTWOBIT SCRIPT PATH ---------------------------
    # Check if faTotTwoBit is available
    if not check_tool_availability("faToTwoBit"):
        print("STEP 5/5: faToTwoBit was not found systemwide.")

        while True:
            fatotwobit_script_path = input("STEP 4/4: Enter path to your faToTwoBit script (e.g., /condaenv/bin/faToTwoBit): ")    
            
            if os.path.exists(fatotwobit_script_path):
                set_fatotwobit_script(fatotwobit_script_path)
                break
            else:
                print('\n' + "‼️ Error: The provided path does not exist. Please ensure the path is correct and try again.")
        
    else:
        print("STEP 5/5: faToTwoBit was found systemwide.")



    # ---------------------------  STEP 1 PART 2: CREATE BAM AND STRAIN SPECIFIC TARGET FOLDER NAMES ---------------------------
    # Get the immediate directory the .bam file is in
    step1_bam_file_folder = os.path.dirname(step1_bam_file) 

    print(step1_bam_file_folder)

    # Create directory with strain specific target name
    # Ex. Output: doris-lab/bio-sample-4/IGHG1_A3
    step1_strain_specific_target_folder = os.path.join(step1_bam_file_folder, step1_strain_specific_target_uppercase)
    
    if not os.path.exists(step1_strain_specific_target_folder):
        os.makedirs(step1_strain_specific_target_folder)

    print(step1_strain_specific_target_folder)

    
    # Create a version of the bam file the user input without the file extension name
    step1_bam_file_without_extension = os.path.splitext(step1_bam_file)[0] # This is for eg: STC654_Bio_Sample_8_flnc
    
    # Create complete strain specific target txt file
    step1_bam_file_prefix = step1_bam_file_without_extension.replace("flnc", "")
    step1_sequences_txt_filename_prefix = step1_bam_file_prefix.replace(f"{step1_bam_file_folder}/", "")
    step1_sequences_txt_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}.COMPLETE.txt"
    step1_sequences_txt_file = os.path.join(step1_strain_specific_target_uppercase, step1_sequences_txt_filename)
    step1_sequences_txt_file = os.path.join(step1_bam_file_folder, step1_sequences_txt_file)

    print('\n')
    
    create_file(step1_sequences_txt_file, "")

    # --------------------------- STEP 1 PART 3: CONVERT BAM FILE TO FASTA FILE ---------------------------
    
    # Generate a new output filename with a .fasta extension
    step1_out_file_basename = os.path.basename(step1_bam_file)
    step1_out_file_fasta = os.path.splitext(step1_out_file_basename)[0] + '.fa'
    step1_out_file = os.path.join(os.path.dirname(step1_bam_file), step1_out_file_fasta)

    # Execute convert bam to fasta function
    convert_bam_to_fasta(step1_bam_file, step1_out_file)


    # --------------------------- STEP 2: SPLIT FASTA FILE INTO 10 MINI FILES ---------------------------
    split_fasta_file(step1_bam_file_without_extension, step1_out_file)

    print('\n')


    # --------------------------- STEP 3: ITERATE THROUGH FASTA FILES ---------------------------
    for i in range(1, 11):

        # --------------------------- STEP 3 PART 1: INDEX THROUGH EACH MINI FASTA FILE ---------------------------
        
        # Create variable file number to dynamically generate files
        step3_part1_file_number = f"{i:02}"
        
        # Set up files 
        step3_part1_fa_file = f"{step1_bam_file_without_extension}.000{step3_part1_file_number}.fa"
        step3_part1_2bit_file = f"{step1_bam_file_without_extension}.000{step3_part1_file_number}.2bit"
        
        # Execute index fasta file function
        index_fasta_file(step3_part1_fa_file, step3_part1_2bit_file)


        # --------------------------- STEP 3 PART 2: EXTRACT SEQUENCES OF INTEREST WITH BLAT --------------------------- 

        # Create output text file
        step3_part2_txt_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}.000{step3_part1_file_number}.txt"
        step3_part2_txt_file_with_folder = os.path.join(step1_strain_specific_target_uppercase, step3_part2_txt_filename)
        step3_part2_txt_file = os.path.join(step1_bam_file_folder, step3_part2_txt_file_with_folder)
        
        # Execute extract sequences function
        extract_sequences_of_interest(step3_part1_2bit_file, step1_strain_specific_target_file, step3_part2_txt_file)


        # --------------------------- STEP 3 PART 3: GRAB IDENTIFIERS FROM EACH TXT FILE --------------------------- 
        
        # Executive extract identifiers function
        step3_part3_identifiers = extract_identifiers(step3_part2_txt_file)


        # --------------------------- STEP 3 PART 4: ADD IDENTIFIERS TO LARGE TXT FILE --------------------------- 
        append_list_of_identifiers(step1_sequences_txt_file, step3_part3_identifiers)


    print('\n')

    # --------------------------- STEP 4: MATCH SEQUENCES TO ORIGINAL ---------------------------

    # Create output text file
    step4_seqs_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}_seqs.fa"
    step4_seqs_file_with_folder = os.path.join(step1_strain_specific_target_uppercase, step4_seqs_filename)
    step4_seqs_file = os.path.join(step1_bam_file_folder, step4_seqs_file_with_folder)

    print(step4_seqs_file)
    
    # Command Usage: seqtk subseq STC654_Bio_Sample_8_flnc.fasta STC654_Bio_Sample_8_IGHG1.txt > IGHG1_seqs.fasta
    match_sequences(step1_out_file, step1_sequences_txt_file, step4_seqs_file)

    print("\nProcess complete! 🎉")

# Performs sequence extraction process
def sequence():
    # --------------------------- STEP 1 PART 1: GET VARIABLES ---------------------------
    # --------------------------- GET BAM FILE ---------------------------
    
    # Get BAM file + Check if BAM file exists
    while True:
        step1_bam_file = input("STEP 1/4: Enter path to the .bam file: ")

        if os.path.exists(step1_bam_file):
            break
        else:
            print('\n' + f"‼️ Error: The file '{step1_bam_file}' does not exist." + '\n')
            step1_bam_file = input("STEP 1/4: Enter path to bam file: ")
    
    print('\n')

    # --------------------------- GET STRAIN SPECIFIC TARGET ---------------------------
    # Get Strain Specific Target File + Check if Strain Specific Target Is one of our options
    while True: 
        step1_strain_specific_target = input(f"Options: {strain_specific_target_options} \nSTEP 2/4: Enter strain specific target (case-insensitive): ")
        step1_strain_specific_target_uppercase = step1_strain_specific_target.upper()

        if step1_strain_specific_target_uppercase in strain_specific_target_options:
            break
        else:
            print('\n' + f"‼️  {step1_strain_specific_target_uppercase} is not an available option. Options: ", strain_specific_target_options + '\n')

    print('\n')


    # --------------------------- GET STRAIN SPECIFIC TARGET FILE ---------------------------
    # Get Strain Specific Target File + Check if Strain Specific Target File Exists

    while True: 
        step1_strain_specific_target_file = input("STEP 3/4: Enter path to strain specific target file (ex: IGHG1_SHRA3.fasta): ")

        if os.path.exists(step1_strain_specific_target_file):
            break
        else:
            print('\n' + f"‼️ Error: The file '{step1_strain_specific_target_file}' does not exist." + '\n')

    print('\n')

    # --------------------------- GET BLAT SCRIPT PATH ---------------------------
    # Check if blat is available
    if not check_tool_availability("blat"):
        print("STEP 4/4: blat was not found systemwide.") 
        
        while True:
            blat_script_path = input("STEP 4/4: Enter path to your blat script (e.g., /condaenv/bin/blat): ")    
            
            if os.path.exists(blat_script_path):
                set_blat_script(blat_script_path)
                break
            else:
                print("‼️ Error: The provided path does not exist. Please ensure the path is correct and try again.")
        
    else:
        print("STEP 4/4: blat was found systemwide.")
    
    print('\n')

    # ---------------------------  STEP 1 PART 2: CREATE BAM AND STRAIN SPECIFIC TARGET FOLDER NAMES ---------------------------
    # Find .fa file created from .bam file
    step1_out_file_basename = os.path.basename(step1_bam_file)
    step1_out_file_fasta = os.path.splitext(step1_out_file_basename)[0] + '.fa'
    step1_out_file = os.path.join(os.path.dirname(step1_bam_file), step1_out_file_fasta)


    # Get the immediate directory the .bam file is in
    step1_bam_file_folder = os.path.dirname(step1_bam_file) 

    print(step1_bam_file_folder)

    # Create directory with strain specific target name
    # Ex. Output: doris-lab/bio-sample-4/IGHG1_A3
    step1_strain_specific_target_folder = os.path.join(step1_bam_file_folder, step1_strain_specific_target_uppercase)
    
    if not os.path.exists(step1_strain_specific_target_folder):
        os.makedirs(step1_strain_specific_target_folder)

    print(step1_strain_specific_target_folder)

    print('\n')

    # Create a version of the bam file the user input without the file extension name
    step1_bam_file_without_extension = os.path.splitext(step1_bam_file)[0] # This is for eg: STC654_Bio_Sample_8_flnc
    
    # Create complete strain specific target txt file
    step1_bam_file_prefix = step1_bam_file_without_extension.replace("flnc", "")
    step1_sequences_txt_filename_prefix = step1_bam_file_prefix.replace(f"{step1_bam_file_folder}/", "")
    step1_sequences_txt_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}.COMPLETE.txt"
    step1_sequences_txt_file = os.path.join(step1_strain_specific_target_uppercase, step1_sequences_txt_filename)
    step1_sequences_txt_file = os.path.join(step1_bam_file_folder, step1_sequences_txt_file)


    # --------------------------- STEP 2 INDEX THROUGH FA FILES ---------------------------
    for i in range(1, 11):
        # --------------------------- STEP 2 PART 1: EXTRACT SEQUENCES OF INTEREST WITH BLAT --------------------------- 

        # Create variable file number to dynamically generate files
        step2_part1_file_number = f"{i:02}"

        # Set up files 
        step2_part1_fa_file = f"{step1_bam_file_without_extension}.000{step2_part1_file_number}.fa"
        step2_part1_2bit_file = f"{step1_bam_file_without_extension}.000{step2_part1_file_number}.2bit"

        # Create output text file
        step2_part1_txt_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}.000{step2_part1_file_number}.txt"
        step2_part1_txt_file_with_folder = os.path.join(step1_strain_specific_target_uppercase, step2_part1_txt_filename)
        step2_part1_txt_file = os.path.join(step1_bam_file_folder, step2_part1_txt_file_with_folder)
        
        # Execute extract sequences function
        extract_sequences_of_interest(step2_part1_2bit_file, step1_strain_specific_target_file, step2_part1_txt_file)

        # --------------------------- STEP 2 PART 2: GRAB IDENTIFIERS FROM EACH TXT FILE --------------------------- 
        
        # Executive extract identifiers function
        step2_part2_identifiers = extract_identifiers(step2_part1_txt_file)


        # --------------------------- STEP 2 PART 3: ADD IDENTIFIERS TO LARGE TXT FILE --------------------------- 
        append_list_of_identifiers(step1_sequences_txt_file, step2_part2_identifiers)

    
    print('\n')

    # --------------------------- STEP 3: MATCH SEQUENCES TO ORIGINAL ---------------------------

    # Create output text file
    step3_seqs_filename = f"{step1_sequences_txt_filename_prefix}{step1_strain_specific_target_uppercase}_seqs.fa"
    step3_seqs_file_with_folder = os.path.join(step1_strain_specific_target_uppercase, step3_seqs_filename)
    step3_seqs_file = os.path.join(step1_bam_file_folder, step3_seqs_file_with_folder)

    print(step1_out_file, step1_sequences_txt_file,step3_seqs_file)
    
    # Command Usage: seqtk subseq STC654_Bio_Sample_8_flnc.fasta STC654_Bio_Sample_8_IGHG1.txt > IGHG1_seqs.fasta
    match_sequences(step1_out_file, step1_sequences_txt_file, step3_seqs_file)


    print("\nProcess complete! 🎉")

# Doris Lab IG Pipeline Tool
def main():
    # --------------------------- CHECK IF TARGET OPTIONS EXISTS ---------------------------
    if len(strain_specific_target_options) > 0:
        # --------------------------- INTRODUCTION ---------------------------
        
        print("\nWelcome to a Doris Lab IG pipeline tool 🔬 \n")

        print(
        '''
    Task Options:

    'start':
    - Performs the entire extraction process from a single bam file.
    - Completes an extraction of identifers for one strain specific target.

    'sequence':
    - Performs only the extraction of identifiers for one strain specific target.
    - .2bit and .fa files for each split 1/10 fasta file must be already created from running 'start'.

        ''')

        task_options = ['start', 'sequence']

        while True:
            task = input("CHOOSE TASK: Enter task to perform (case-insensitive): ")
            task_lowercase = task.lower()

            if task_lowercase in task_options:
                break
            else:
                print('\n' + f"❌ Error: The option {task_lowercase} is not available."  + '\n')
        
        print('')


        # --------------------------- START ---------------------------
        
        if task_lowercase == 'start':
            start()
        elif task_lowercase == 'sequence':
            sequence()

    else: 
        print('\n' + "❌ Error: could not find file targets.txt file with strain specific target options.")

# Set strain specific target options global array
set_target_options()

# Call main function
if __name__ == "__main__":
    main()