# IG-Buddy: Doris Lab IG Pipeline Tool
IGBuddy automates the data extraction and processing pipeline for the Doris Lab by converting and analyzing sequence data. This tool converts `.bam` files to `.fasta`, splits `.fasta` files, indexes them, and extracts sequences based on specific targets using various bioinformatics tools.

## Requirements

To run this tool, the following programs must be installed and accessible in your system's PATH:

- **samtools 1.9** - for converting `.bam` files to `.fasta`
- **seqtk** - for splitting `.fasta` files and retrieving specific sequences
- **fatotwobit** - for converting `.fasta` files to `.2bit`
- **blat** and **blatSrc** - for aligning sequences to targets

Ensure that these tools are available on your system before starting.

## Setup Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/username/Antibuddy.git
   cd Antibuddy
   ```
   
2. Prepare a folder structure in the following format:
   - Create a folder named `bio-sample-[number]`.
   - Place the `.bam` file for your sample within this folder.

3. Place your strain-specific target files in a designated location (or in the same directory as this tool for ease of use).

## Usage
### Set Up Environment
1. Run the setup functions to initialize the environment variables:

   ```python
   set_target_options()  # Loads strain-specific target options from `targets.txt`.
   set_fatotwobit_script("/path/to/fatotwobit")  # Set the path to the faToTwoBit script.
   set_blat_script("/path/to/blat")  # Set the path to the BLAT script.
   set_igblast_script("/path/to/igblast")  # Set the path to the IgBlast script.
   ```

2. Start the data processing pipeline with your `.bam` file by following these steps:

### Steps

1. **Convert `.bam` to `.fasta`:**
   ```python
   convert_bam_to_fasta("sample.bam", "sample.fasta")
   ```
   Converts `.bam` files to `.fasta` format for further processing.

2. **Split `.fasta` into smaller chunks:**
   ```python
   split_fasta_file("sample.fasta", "output_prefix")
   ```
   Splits the `.fasta` file into 10 smaller files for more manageable analysis.

3. **Index `.fasta` files using `faToTwoBit`:**
   ```python
   index_fasta_file("chunk_0001.fasta", "chunk_0001.2bit")
   ```
   Converts `.fasta` files to `.2bit` format for faster access by BLAT.

4. **Extract sequences of interest using BLAT:**
   ```python
   extract_sequences_of_interest("database.2bit", "query.txt", "output.txt")
   ```
   Uses BLAT to extract sequences that match specific targets, outputting them in BLAST format.

5. **Extract identifiers with the highest score:**
   ```python
   identifiers = extract_identifiers("target_file.txt")
   ```
   Retrieves identifiers of sequences with the highest score based on strain-specific targets.

6. **Append identifiers to a master file:**
   ```python
   append_list_of_identifiers("master_file.txt", identifiers)
   ```
   Adds identifiers to a specified master file for further analysis.

7. **Retrieve sequences that match specific identifiers:**
   ```python
   match_sequences("sample.fasta", "identifiers.txt", "matching_sequences.fasta")
   ```
   Extracts sequences from `.fasta` file based on a list of identifiers.

## Error Handling
Errors during file operations or external command executions will be caught and displayed, making it easier to troubleshoot issues such as missing files or incorrect paths.

## Contributing
If you'd like to contribute to this project, please fork the repository and use a feature branch. Pull requests are welcome.

## License
This project is licensed under the MIT License.
