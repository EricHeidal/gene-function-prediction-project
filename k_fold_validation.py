import random

from sklearn import svm

from sklearn.model_selection import train_test_split  

class Gene:
    def __init__(self, name):
        self.name = name
        self.go_terms = {}  # Dictionary is go term -> boolean


class CellLine:
    def __init__(self, name):
        self.name = name
        self.expression_levels = {}  # Dictionary is gene name -> float


with open("screening_data.tsv", "r") as screening_data_file, open("labels.tsv", "r") as label_file:
    go_term_list = []
    cell_line_list = []
    gene_list = []
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
                    gene.go_terms[go_term_list[j]] = False
                elif annotated == "1":
                    gene.go_terms[go_term_list[j]] = True
                else:
                    gene.go_terms[go_term_list[j]] = None
            gene_list.append(gene)


# SVM Algorithm Exists

# K-Fold Cross Validation Algorithm


# Split the data set into k groups
# For each unique group
    # Take the group as a hold out or test data set
    # Take the remaining groups as training data
    # Fit a model on the training set and evaluate it on the test set
    # Retain the evaluation score and discard the model
    # Evaluation score = percent accuracy
# Summarize the skill of the model using the sample of model evaluation score

gene_list_copy = gene_list.copy()
# Shuffle the data set randomly
random.shuffle(gene_list_copy)
print("length of gene_list_copy ", len(gene_list_copy))

evalationScores = []
window = 1780



# Split the data set into k groups
# For each unique group
for k in range(10):

    print(gene_list_copy[k])
    #splitting data into training set and test set
    test = gene_list_copy[window - 1780: window]
    train = gene_list_copy[0 :window-1780] + gene_list_copy[window:]
    print("Size of test", len(test))
    print("Size of train" , len(train))

    #incrementing window for test data
    window += 1780

    #return evaluation score from svm method using test and training data


    #append evaluation score to array


#average evaulation scores in array?