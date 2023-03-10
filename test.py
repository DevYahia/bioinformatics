dataset = "GTGTCAGCCGTGAACAGCGATCTGGTTTTTGCTACGATTATGTAAAGCTGGAGGCTGTTCATTTCTTCGCTCGTGTTGTTAGGCCGGGATATAACACAACGCGACTATCCCTAAACGGTATCTTTCCGATCTGGTTCATGCATGGGCTGCCATATATAGCGCCCCCTGGTTCTTTCTAAGTCCCACTTGTGTTCCCTTTTACGATGTGTGCGTCGCGTCTTTCTGTAGAAAACTTCCATCGGGGCGCGAAGAAACGCCCAGTGGCGGGACGGTGTGGTCGTACGGATAGCATAAGGTCGGGAAACTACGTGGCCATTGTTACAGAATGTTATTAATGAAGCCAGCGACTCTAGCAACACTCCAAGCCCTCAGAGTGTAGTCATCAGAGTGTCACAGCTTAGGCAAGCTTGCCACCAGTGAAGCTCTTGCTTATTAGTAGCACTCGCGCGTGAGAGCCATAGCAGCTCGCTTCTACCCCAGAGATCAAGTGTCTTTGTTGAGATACGCAACCTCGTTTTGAGTTGAGCGGGTCGGCACAAAGTACTGCTTGCTGGAGGCGCGACCAAGGCTTCCATGGGGAACCACCCGGTCGACTCGCTAGTTTGCTGGTATCATGCGTCCAAAATAATACGATGGCCACGGGCCAGAAGATGCAAGATGCCTTTCGCCGGCGACGTACTAGATATAATCGTAGAAGGGAGCTCGGTCGGAAGATAGTCTTCTAAATCTTACCTATGCAGTTTGACTGAACTACCGTGGGGCTGCTCCCCTGATTTTCTGTGGAAACCTCATCTCCCTAACTTCTTTGACATGATTCTCCAACATAGAGCCTTTTACATGATGCCCCCACGGCCTATGGATCCTGTAAGATAGGTGAGGATACCCTACTGCTACTCGCCAGCCCATTATCAGGGGCGGAAATGTAGCTCAAGCAGCCTAAAAAGTTTAGGGCAGCGTTATCCATGGACCCGTTTGCTCCCATCCTAGGGC"


def count(string: str, letter: str) -> int:
    return string.count(letter)


print(count(dataset, "A"))
print(count(dataset, "C"))
print(count(dataset, "G"))
print(count(dataset, "T"))
