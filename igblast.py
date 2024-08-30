import subprocess

def run_igblast(seqsfile, outputfile, species="rat"):
    """
    Run IgBLAST with specified parameters.

    :param seqsfile: Path to the input sequence file.
    :param outputfile: Path to the output file.
    :param species: Species for germline gene database, default is 'rat'.
    """
    # Define the IgBLAST parameters based on your input
    program = "blastn"
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
        "igblastn",  # IgBLAST executable
        "-query", seqsfile,  # Input sequence file
        "-out", outputfile,  # Output file
        "-germline_db_V", v_gene_database,
        "-germline_db_D", d_gene_database,
        "-germline_db_J", j_gene_database,
        "-penalty", str(v_mismatch_penalty),
        "-min_D_match", str(min_d_gene_matches),
        "-D_penalty", str(d_mismatch_penalty),
        "-min_V_length", str(min_v_gene_length),
        "-min_J_length", str(min_j_gene_length),
        "-J_penalty", str(j_mismatch_penalty),
        "-num_alignments_V", str(num_v_gene_alignments),
        "-num_alignments_D", str(num_d_gene_alignments),
        "-num_alignments_J", str(num_j_gene_alignments),
        "-outfmt", alignment_format,
        "-num_alignments", "10",  # Number of alignments to show
        "-evalue", "1",  # Expect value threshold for saving hits
        "-show_translation",  # Show amino acid translation
        "-domain_system", "imgt"  # V domain delineation system
    ]

    # Run the IgBLAST command
    try:
        subprocess.run(igblast_command, check=True)
        print(f"IgBLAST ran successfully. Output saved to {outputfile}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running IgBLAST: {e}")

# Example usage
run_igblast("path/to/your/sequences.fasta", "path/to/output/file.txt")

