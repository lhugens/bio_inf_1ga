import re

def get_statistics(seq):
    letters = ["A", "C", "G", "T"]
    letters_dic = {"A":0, "C":0, "G":0, "T":0}

    for letter in seq:
        if letter in letters_dic:
            letters_dic[letter] = letters_dic.get(letter, 0) + 1

    for letter in reverse_complement(seq):
        if letter in letters_dic:
            letters_dic[letter] = letters_dic.get(letter, 0) + 1

    for key in letters:
        letters_dic[key] = letters_dic[key]/(len(seq)*2) * 100
     
    return letters_dic

def generate_codon_combinations():
    letters = ["A", "C", "G", "T"]

    codons = []
    
    for l1 in letters:
        for l2 in letters:
            for l3 in letters:
                codons.append(l1 + l2 + l3)

    return codons

def count_codon_occurances(seq, codon):
    codon_occs = 0
    
    for i in range(0, len(seq)):
        if seq[i:i + 3] == codon:
            codon_occs += 1

    return codon_occs

def reverse_complement(seq):
    res = ""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    for base in seq:
        if base not in complement:
            return "Error: Not a DNA sequence"
        else:
            res = complement[base] + res

    return res

def do_for_all(seq, func, params):
    rseq = reverse_complement(seq)
    res = []
    for i in range(3):        
        res.append(func(seq[i:], params))

    for i in range(3):    
        res.append(func(rseq[i:], params))

    return res

def count_all_codon_occurances(seq):
    codons = generate_codon_combinations()
    codon_occs = {}
    for codon in codons:
        codon_occs[codon] = sum(do_for_all(seq, count_codon_occurances, codon))
    return codon_occs

def step_A(seq):
    bps_dic = get_statistics(seq)
    codon_dic = count_all_codon_occurances(seq)
    least_common_codon = min(codon_dic, key=codon_dic.get)
    most_common_codon = max(codon_dic, key=codon_dic.get)

    print("1. ", len(seq))
    print("2. ", end=" ")
    for i, (key, value) in enumerate(bps_dic.items()):
        if i != 3:
            print(key, "=", value, "%,", end=" ")
        else:
            print(key, "=", value, "%")
    print("3. ", bps_dic["G"] + bps_dic["C"], "%")
    print("4. ", codon_dic["ATG"])
    print("5. ", codon_dic["TAA"] + codon_dic["TAG"] + codon_dic["TGA"])
    print("6.  Least common =", least_common_codon, ", Most common =", most_common_codon)

def find_all_occurences(seq, codon, step):
    codon_occs = []

    for i in range(0, len(seq), step):
        if seq[i:i + 3] == codon:
            codon_occs.append(i)
    return codon_occs

def get_all_orfs(seq, params=None):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    i = 0
    
    while i < len(seq):
        if seq[i:i + 3] == start_codon:
            for j in range(i, len(seq), 3):
                if seq[j:j + 3] in stop_codons:   
                    if j + 3 - i >= 149:
                        orfs.append((seq[i:j + 3], (i+1, j + 3)))
                    i = j + 3
                    break
        i += 3
     
    return orfs

def correct_orfs(orfs):
    info =  ""
    new_orfs = []
    for i in range(0, len(orfs)):
        for orf in orfs[i]:
            if i > 2:
                new_cord_first = int(orf[1][0] + i-3)
                new_cord_second = int(orf[1][1] + i-3)
                new_cord_first = 30000 - new_cord_first + 1
                new_cord_second = 30000 - new_cord_second + 1
            else:
                new_cord_first = int(orf[1][0] + i)
                new_cord_second = int(orf[1][1] + i)

            info = info + str(new_cord_first) + ", " + str(new_cord_second)

            new_orfs.append((orf[0], info))
            info = ""

    return new_orfs   

def step_B(seq):
    orfs = do_for_all(seq, get_all_orfs, None)
    orfs = correct_orfs(orfs)
    
    with open("all_potential_proteins.fasta", "w") as f:
        with open("orf_coordinates.txt", "w") as g:
            for i in range(0, len(orfs)):
                f.write(f'>Protein_ORF{i+1}\n{orfs[i][0]}\n')
                g.write(f'{orfs[i][1]}, ORF{i+1}\n')

def clean_genes(genes):
    new_genes = []
    for gene in genes:
        gene = gene.split("\t")
        if gene[2] == "exon":
            gene_id = gene[8].split(" ")[1][1:-2]
            # if gene[6] == "-":
            #     new_genes.append([gene[4], gene[3], gene_id])
            # else:
            new_genes.append([gene[3], gene[4], gene_id])
    return new_genes

def read_orfs():
    with open("orf_coordinates.txt") as f:
        lines = f.readlines()
        orfs = []
        for l in lines:
            l_split = l.split(", ")
            start = l_split[0]
            end = l_split[1]
            name = l_split[2][0:-1]
            if start > end:
                orfs.append([end, start, name])
            else:
                orfs.append([start, end, name])

    return orfs

def calculate_overlap(gene, orf):
    gene_start = int(gene[0])
    gene_end = int(gene[1])
    orf_start = int(orf[0])
    orf_end = int(orf[1])

    gene_length = gene_end - gene_start
    orf_length = orf_end - orf_start

    start_disalign = gene_start - orf_start
    end_disalign = gene_end - orf_end
    return 1 - (abs(start_disalign + end_disalign) / gene_length)

def step_C(seq, genes):
    genes = clean_genes(genes)
    orfs = read_orfs()
    for gene in genes:
        max_ov = 0
        max_ov_i = 0
        for i in range(0, len(orfs)):
            ov = calculate_overlap(gene, orfs[i])
            if max_ov < ov:
                max_ov = ov
                max_ov_i = i

        print("{:<5} {:<5} {:>5}%".format(orfs[max_ov_i][2], gene[2], round(max_ov * 100, 2)))   
    
def main():
    with open("./yeast/sequence_chr1.fasta") as seq_chr1_fasta:
        lines = seq_chr1_fasta.readlines()
        seq = ""
        for l in lines:
            if l[0] != ">":
                seq += l.replace("\n", "")
        seq = seq[:30000]

    with open("./yeast/genes_chr1.gtf") as genes_chr1_gtf:
        lines = genes_chr1_gtf.readlines()
        genes = []
        for l in lines:
            genes.append(l)

    step_A(seq)  
    step_B(seq)
    step_C(seq, genes)      

    return 0

if __name__ == "__main__":
    main()
