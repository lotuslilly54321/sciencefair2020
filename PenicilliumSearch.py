from Bio import SeqIO
import sys


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
    data = data.strip('Z')
    return data


#
#
#
def naive_approximate(gene, chromosome, max_hamming_distance):
    occurs = []
    for i in xrange(0, len(chromosome) - len(gene) + 1):  # for	all	alignments
        nmm = 0
        for j in xrange(0, len(gene)):
            if chromosome[i + j] != gene[j]:  # does	it	match?
                nmm += 1
                if nmm > max_hamming_distance:
                    break  # exceeded	maximum	distance
        if nmm <= max_hamming_distance:
            occurs.append((i, nmm))
    return occurs


chr1 = read_fasta_file('data/chr1.FASTA')
chr2 = read_fasta_file('data/chr2.FASTA')
chr3 = read_fasta_file('data/chr3.FASTA')
chr4 = read_fasta_file('data/chr4.FASTA')

chromosomes_dictionary = {
    "chromosome1": chr1,
    "chromosome2": chr2,
    "chromosome3": chr3,
    "chromosome4": chr4,
}

pcbAB = read_fasta_file('data/A6N339.fasta')
pcbC = read_fasta_file('data/P08703.fasta')
penDE = read_fasta_file('data/P15802.fasta')

genes_dictionary = {
    "pcbAB": pcbAB,
    "pcbC": pcbC,
    "penDE": penDE,
}


def protein_parts(protein):
    parts = []
    for i in range(0, len(protein), 10):
        start_index = i
        end_index = min(i + 10, len(protein))
        parts.append(str(protein[start_index:end_index]))
    return parts


# attempt to find using hamming distance on the entire protein

max_hamming_distance = 3

genes_found_dictionary = {}

for chromosome_key, chromosome_value in chromosomes_dictionary.items():
    for gene_key, gene_value in genes_dictionary.items():
        occurrences = naive_approximate(gene_value, chromosome_value, max_hamming_distance)
        for occurrence in occurrences:
            location, num_differences = occurrence
            genes_found_dictionary[gene_key] = location
            print('found gene {} in chromosome {} at location {} with {} differences'.format(
                gene_key, chromosome_key, location, num_differences))

# for any protein not matched, its possible that there was a shift (insertion or deletion)
# that caused the hamming distance algorithm to fail

genes_not_found_list = []
for gene_key, gene_value in genes_dictionary.items():
    if gene_key not in genes_found_dictionary:
        genes_not_found_list.append(gene_key)

print('Genes not found {}'.format(genes_not_found_list))
print('Attempting to run hamming distance on shorter substrings')

max_allowed_parts_mismatch_percent = 2

# attempt to run hamming distance on a shorter substrings
for chromosome_key, chromosome_value in chromosomes_dictionary.items():
    for gene_key, gene_value in genes_dictionary.items():
        if gene_key in genes_not_found_list:
            gene_part_values = protein_parts(gene_value)
            max_allowed_parts_mismatch = int(len(gene_part_values) * 1.0 * max_allowed_parts_mismatch_percent) / 100.0
            parts_mismatch_count = 0
            search_start_index = 0
            min_location = sys.maxint
            for gene_part_value in gene_part_values:
                occurrences = naive_approximate(gene_part_value, chromosome_value[search_start_index:],
                                                max_hamming_distance)
                # print('occurrences={}'.format(occurrences))
                min_occurrence_location = sys.maxint
                for occurrence in occurrences:
                    location, num_differences = occurrence
                    location = location + search_start_index
                    min_occurrence_location = min(min_occurrence_location, location)
                    # print('**debug** found part_value {} in chromosome {} at location {}'.format(gene_part_value,
                    #                                                                              chromosome_key,
                    #                                                                              location))
                search_start_index = max(search_start_index, min_occurrence_location)
                min_location = min(min_location, min_occurrence_location)
                if len(occurrences) == 0:
                    # print('**debug** part_value {} not found in chromosome {}'.format(
                    # gene_part_value,
                    # chromosome_key))
                    parts_mismatch_count += 1
            if parts_mismatch_count <= max_allowed_parts_mismatch:
                print('found gene {} in chromosome {} at location {} within mismatch percent of {}'.format(
                    gene_key, chromosome_key, min_location, max_allowed_parts_mismatch_percent))
