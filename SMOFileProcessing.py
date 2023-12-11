class Gene:
    def __init__(self, name):
        self.name = name
        # Simple list of expression levels - careful not to shuffle this
        self.fitness_levels = []


class GoTerm:
    def __init__(self, name):
        self.name = name
        # Dictionary from gene_name -> ([fitness_levels], annotation)
        self.gene_annotations = {}


def process_file():
    go_term_list = []
    gene_list = []
    num_characteristics = 0
    with open("screening_data.tsv", "r") as screening_data_file, open("labels.tsv", "r") as label_file:
        for i, fitness_file_line in enumerate(screening_data_file):
            # Read in and split the next line of each file
            split_fitness_file_line = fitness_file_line.split()
            num_characteristics = len(split_fitness_file_line) - 1
            split_label_file_line = label_file.readline().split()
            if i == 0:  # If it's the first line of each file
                for term in split_label_file_line:
                    go_term_list.append(GoTerm(term))
            else:
                gene = Gene(split_fitness_file_line.pop(0))  # First thing on the line is the name of the gene
                # Give the Gene all its fitness levels
                for fitness_level in split_fitness_file_line:
                    gene.fitness_levels.append(float(fitness_level))
                split_label_file_line.pop(0)

                # Assign the current gene its respective go term annotations
                for j, annotated in enumerate(split_label_file_line):
                    if annotated == "-1":
                        go_term_list[j].gene_annotations[gene.name] = (gene.fitness_levels, -1)
                    elif annotated == "1":
                        go_term_list[j].gene_annotations[gene.name] = (gene.fitness_levels, 1)
                gene_list.append(gene)
    print("Finished processing input file")
    return num_characteristics, gene_list, go_term_list
