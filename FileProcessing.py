class Gene:
    def __init__(self, name):
        self.name = name
        self.go_terms = {}  # Dictionary is go term -> boolean


class CellLine:
    def __init__(self, name):
        self.name = name
        self.expression_levels = {}  # Dictionary is gene name -> float


def process_file():
    go_term_list = []
    cell_line_list = []
    gene_list = []
    with open("screening_data.tsv", "r") as screening_data_file, open("labels.tsv", "r") as label_file:
        for i, expression_file_line in enumerate(screening_data_file):
            # Read in and split the next line of each file
            split_expression_file_line = expression_file_line.split()
            split_label_file_line = label_file.readline().split()
            if i == 0:  # If it's the first line of each file
                for term in split_label_file_line:
                    go_term_list.append(term)
                for cell_line_name in split_expression_file_line:
                    cell_line_list.append(CellLine(cell_line_name))
            else:
                gene = Gene(split_expression_file_line.pop(0))  # First thing on the line is the name of the gene
                for j, expression_level in enumerate(split_expression_file_line):
                    # For each expression level on the line, go into the corresponding column's cell line
                    # and place the (gene_name, expression_level) pair in it's dictionary
                    cell_line_list[j].expression_levels[gene.name] = float(expression_level)
                split_label_file_line.pop(0)
                # Assign the current gene its respective go term annotations
                for j, annotated in enumerate(split_label_file_line):
                    if annotated == "-1":
                        gene.go_terms[go_term_list[j].name] = -1
                    elif annotated == "1":
                        gene.go_terms[go_term_list[j].name] = 1
                    else:
                        gene.go_terms[go_term_list[j].name] = None
                gene_list.append(gene)
    return gene_list, go_term_list
