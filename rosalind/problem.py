file_path = "rosalind_rna.txt"

with open(file_path) as f:
    lines = f.readlines()
    dna = lines[0]

    transcribed_rna = dna.replace("T", "U")

    print(transcribed_rna)
