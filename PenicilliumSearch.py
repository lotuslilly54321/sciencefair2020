from Bio import SeqIO


#
#
#
def read_fasta_file(filename):
    data = ""
    with open(filename, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # identifier = record.id
            # description = record.description
            data = data + record.seq
    data.strip('Z')
    return data


#
#
#
def naive_approximate(protein, string):
    max_hamming_distance = 0
    occurrences = []
    for i in xrange(0, len(string) - len(protein) + 1):  # for	all	alignments
        nmm = 0
        for j in xrange(0, len(protein)):
            if string[i + j] != protein[j]:  # does	it	match?
                nmm += 1
                if nmm > max_hamming_distance:
                    break  # exceeded	maximum	distance
        if nmm <= max_hamming_distance:
            occurrences.append((i,nmm))
    return occurrences


chr1 = read_fasta_file('chr1.FASTA')
chr2 = read_fasta_file('chr2.FASTA')
chr3 = read_fasta_file('chr3.FASTA')
chr4 = read_fasta_file('chr4.FASTA')

chromosomes_dictionary = {
    "chromosome1": chr1,
    "chromosome2": chr2,
    "chromosome3": chr3,
    "chromosome4": chr4,
}

pcbAB = read_fasta_file('A6N339.fasta')
pcbC = read_fasta_file('P08703.fasta')
penD = read_fasta_file('P15802.fasta')

proteins_dictionary = {
    "pcbAB": pcbAB,
    "pcbC": pcbC,
    "penD": penD,
}

for chromosome_key, chromosome_value in chromosomes_dictionary.items():
    for protein_key, protein_value in proteins_dictionary.items():
        print('checking for protein {} in chromosome {}'.format(protein_key, chromosome_key))

