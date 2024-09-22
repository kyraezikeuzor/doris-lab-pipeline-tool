# Antibuddy: Doris Lab IG Pipeline Tool
import subprocess
import os
from spinner import Spinner
from utils import check_tool_availability
from utils import create_file

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
        print("\n[Error] The file targets.txt was not found. ")
    
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

# Capture global blat script path
def set_igblast_script(script_path):
    global igblast_script
    igblast_script = script_path
     

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
        print('\n' + "[Success] Converted bam to fasta file successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "[Error] Error running Samtools command:" + e.stderr.decode() + '\n')
    
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
        print('\n' + "[Success] Fasta files split successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "[Error] Error running Seqtk command:" + e.stderr.decode() + '\n')

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
        print('\n' + "[Success] .2bit from .fa file generated successfully." + result.stdout.decode() + '\n')

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "[Error] Error running faToTwoBit: " + e.stderr.decode() + '\n')

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
        print('\n' + "[Success] Extracted sequences with blat successfully.", result.stdout.decode())

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "[Error] Error running BLAT command:" + e.stderr.decode() + '\n')

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
    print(f"[Success] Extracted sequence idenfifiers from {filename} successfully.")

    return identifiers

# Add identifiers from each 0000X.<strain-specific-target> file to the complete file
def append_list_of_identifiers(filename, identifiers):
    # Adding copied list of items from each mini ighg1 file to a giant IGHh1 file
    with open(filename, 'a') as file:  # Open the file in append mode
        for identifier in identifiers:
            file.write(identifier + '\n')  # Write each identifier on a new line
    
    print("\n[Success] Added identifiers to " + filename + '\n')

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
        print("\n[Success] Matched sequences successfully.")

    except subprocess.CalledProcessError as e:
        # Print the error message if the command fails
        print('\n' + "[Error] Error running Seqtk command:" + e.stderr.decode() + '\n')
    
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
        step2_part1_txt_filename_with_folder = os.path.join(step1_strain_specific_target_uppercase, step2_part1_txt_filename)
        step2_part1_txt_file = os.path.join(step1_bam_file_folder, step2_part1_txt_filename_with_folder)
        
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
    step3_seqs_filename_with_folder = os.path.join(step1_strain_specific_target_uppercase, step3_seqs_filename)
    step3_seqs_file = os.path.join(step1_bam_file_folder, step3_seqs_filename_with_folder)

    print(step1_out_file, step1_sequences_txt_file,step3_seqs_file)
    
    # Command Usage: seqtk subseq STC654_Bio_Sample_8_flnc.fasta STC654_Bio_Sample_8_IGHG1.txt > IGHG1_seqs.fasta
    match_sequences(step1_out_file, step1_sequences_txt_file, step3_seqs_file)

    print("\nProcess complete! [Success]")


def igblast(species="rat"):
    """
    Run IgBLAST with specified parameters.

    :param seqsfile: Path to the input sequence file.
    :param outputfile: Path to the output file.
    :param species: Species for germline gene database, default is 'rat'.
    """

    # Get seqs file + Check if seqs file exists
    while True:
        seqsfile = input("STEP 1/2: Enter path to the seqs file: ")

        if os.path.exists(seqsfile):
            break
        else:
            print('\n' + f"‼️ Error: The file '{seqsfile}' does not exist." + '\n')
            seqsfile = input("STEP 1/2: Enter path to seqs file: ")
    
    print('\n')

    # --------------------------- GET IGBLAST SCRIPT PATH ---------------------------
    # Check if blat is available
    if not check_tool_availability("igblast"):
        print("STEP 2/2: igblast was not found systemwide.")
        
        while True:
            igblast_script_path = input("STEP 2/2: Enter path to your igblast script (e.g., /condaenv/bin/igblast): ")    
            
            if os.path.exists(igblast_script_path):
                set_igblast_script(igblast_script_path)
                break
            else:
                print('\n' + "‼️ Error: The provided path does not exist. Please ensure the path is correct and try again.")
    else:
        print("STEP 2/2: igblast was found systemwide.")
    
    print('\n')

    # Generate a dynamic output file name based on the input file name
    input_base_name = os.path.splitext(seqsfile)[0]  # Remove the file extension
    outputfile = f"{input_base_name}.AIRR.original.tsv"  # Create the dynamic output file name



    # Define the IgBLAST parameters based on your input
    v_gene_database = f"IMGT_{species}_V"
    d_gene_database = f"IMGT_{species}_D"
    j_gene_database = f"IMGT_{species}_J"
    v_mismatch_penalty = -1
    min_d_gene_matches = 5
    d_mismatch_penalty = -2
    min_v_gene_length = 9
    min_j_gene_length = 0
    j_mismatch_penalty = -2
    num_v_gene_alignments = 1
    num_d_gene_alignments = 1
    num_j_gene_alignments = 1
    alignment_format = "AIRR"  # Using AIRR rearrangement tabular format

    # Construct the IgBLAST command
    igblast_command = [
        igblast_script,  # IgBLAST executable
        "-query", seqsfile,  # Input sequence file
        "-out", outputfile,  # Dynamic output file
        "-germline_db_V", v_gene_database,
        "-germline_db_D", d_gene_database,
        "-germline_db_J", j_gene_database,
        "-V_penalty", str(v_mismatch_penalty),
        "-min_D_match", str(min_d_gene_matches),
        "-D_penalty", str(d_mismatch_penalty),
        "-min_V_length", str(min_v_gene_length),
        "-min_J_length", str(min_j_gene_length),
        "-J_penalty", str(j_mismatch_penalty),
        "-num_alignments_V", str(num_v_gene_alignments),
        "-num_alignments_D", str(num_d_gene_alignments),
        "-num_alignments_J", str(num_j_gene_alignments),
        "-outfmt", "19",  # Correct AIRR tabular format code
        "-num_alignments", "10",  # Number of alignments to show
        "-evalue", "1",  # Expect value threshold for saving hits
        "-show_translation",  # Show amino acid translation
        "-domain_system", "imgt",  # V domain delineation system
        "-extend_align5end",  # Extend alignment at 5' end
        "-extend_align3end",  # Extend alignment at 3' end
        "-allow_vdj_overlap"  # Allow V(D)J overlap
    ]

    # Run the IgBLAST command
    try:
        subprocess.run(igblast_command, check=True)
        print(f"IgBLAST ran successfully. Output saved to {outputfile}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running IgBLAST: {e}")
            

import pandas as pd

def immuneref():
    input_file = input('Enter excel file name: ')  # Replace with your actual file path

    try:
        # Attempt to read the file as an .xls file using the xlrd engine
        df = pd.read_excel(input_file, engine='xlrd')
        print("File read successfully as .xls.")
    except Exception as e:
        print(f"Unexpected error reading the Excel file: {e}")
        return

    # Continue with processing if the file is read successfully
    print("Processing data...")

    # Step 1: Sort by "Stop codon" and remove rows where "Stop codon" is 'T'
    df_sorted = df.sort_values(by='stop_codon', ascending=False)
    df_filtered = df_sorted[df_sorted['stop_codon'] != 'T']

    print(df_filtered)

    # Step 2: Sort by "V alignment end" and remove rows where "V alignment end" < 290
    df_sorted_v_end = df_filtered.sort_values(by='v_alignment_end', ascending=True)
    df_filtered_v_end = df_sorted_v_end[df_sorted_v_end['v_alignment_end'] >= 290]

    print(df_filtered)

    # Step 3: Select specific columns and create a new DataFrame
    selected_columns = ['v_call', 'd_call', 'j_call', 'junction_aa']
    new_df = df_filtered_v_end[selected_columns]

    # Remove the thing he said to from the v_call, d_call, j_call columns
    for column in ['v_call', 'd_call', 'j_call']:
        new_df.loc[:, column] = new_df[column].str.replace(r'\*01', '', regex=True)
        new_df.loc[:, column] = new_df[column].str.replace(r'\*02', '', regex=True)
        new_df.loc[:, column] = new_df[column].str.replace(r'\*03', '', regex=True)

    # Ensure all columns used in agg are strings
    for column in ['v_call', 'd_call', 'j_call', 'junction_aa']:
        new_df.loc[:,column] = new_df[column].astype(str)

    # Combine columns into a single column
    new_df.loc[:,'combined'] = new_df[['v_call', 'd_call', 'j_call', 'junction_aa']].agg('-'.join, axis=1)

    # Count frequencies of each unique 'combined' sequence
    freq_df = new_df['combined'].value_counts().reset_index()
    freq_df.columns = ['combined', 'counts']

    # Calculate frequencies
    total_counts = freq_df['counts'].sum()
    freq_df['freqs'] = freq_df['counts'] / total_counts

    # Expand the 'combined' column back into separate columns
    expanded_df = freq_df['combined'].str.split('-', expand=True)
    
    # Check number of columns created
    num_columns = expanded_df.shape[1]
    expected_columns = ['v_call', 'd_call', 'j_call', 'junction_aa']
    
    # Adjust column names if the number of columns doesn't match expected
    if num_columns == len(expected_columns):
        expanded_df.columns = expected_columns
    else:
        print(f"Warning: Expected {len(expected_columns)} columns but got {num_columns}.")
        # Handle column mismatch (e.g., fill missing columns with NaN or adjust column names)
        expanded_df.columns = [f'column_{i}' for i in range(num_columns)]

    final_df = pd.concat([expanded_df, freq_df[['freqs']]], axis=1)

    # Determine the new file name and path
    directory = os.path.dirname(input_file)
    new_file_name = 'STC654-Bio-Sample-1-IGHG2B_A3_junc_aa_freqs.xlsx'
    output_file = os.path.join(directory, new_file_name)

    
    # Save the new DataFrame to the new Excel file
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        final_df.to_excel(writer, sheet_name='Filtered_Data', index=False)

    print(f"Filtered data has been saved to {output_file}")




# Doris Lab Ig Pipeline Tool
def main():
    # --------------------------- CHECK IF TARGET OPTIONS EXISTS ---------------------------
    if len(strain_specific_target_options) > 0:
        # --------------------------- INTRODUCTION ---------------------------
        
        print(
'''
<----- Doris Lab Immunoglobulin Pipeline ----->

Choose a task to get started:

'start':
- Performs the entire extraction process from a single bam file.
- Completes an extraction of identifers for one strain specific target.

'sequence':
- Performs only the extraction of identifiers for one strain specific target.
- .2bit and .fa files for each split 1/10 fasta file must be already created from running 'start'.

'igblast':
- Runs the IgBlast tool for a seqs file

''')

        task_options = ['start', 'sequence', 'igblast', 'immuneref']

        while True:
            task = input("CHOOSE TASK: Enter task to perform (case-insensitive): ")
            task_lowercase = task.lower()

            if task_lowercase in task_options:
                break
            else:
                print('\n' + f"[Error] Error: The option {task_lowercase} is not available."  + '\n')
        
        print('')


        # --------------------------- START ---------------------------
        
        if task_lowercase == 'start':
            start()
        elif task_lowercase == 'sequence':
            sequence()
        elif task_lowercase == 'igblast':
            igblast()
        elif task_lowercase == 'immuneref':
            immuneref()

    else: 
        print('\n' + "[Error] Error: could not find file targets.txt file with strain specific target options.")


# Set strain specific target options global array
set_target_options()

# Call main function
if __name__ == "__main__":
    main()