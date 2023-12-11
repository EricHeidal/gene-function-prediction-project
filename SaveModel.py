import csv
import numpy
import math
from SMOFileProcessing import *


# Model is saved as csv with format "go_term, bias, coefficient1, coefficient2..."
def save_model(model_dictionary):
    with open("model.weights", 'w', encoding="utf-8") as f:
        for go_term in model_dictionary.keys():
            f.write(go_term + ",")
            for gene, coefficient in model_dictionary[go_term]:
                f.write(str(coefficient) + ",")
            f.write("\n")


def save_single_weight(file_name, go_term, w, b):
    # Appending to file so we don't overwrite data every time - if program crashes or we want to train in chunks
    # this will be useful
    with open(file_name, 'a') as f:
        f.write(go_term + ",")
        f.write(str(b) + ",")
        for weight in w:
            f.write(str(weight) + ",")
        f.write("\n")
        print("Saved a model")


# To predict, pass in a list of fitness levels
def predict(file_name, hundred_list, go_term_list, current_go_term_index):
    with open(file_name, 'r') as f:
        # For every weight vector in the file:
        confidence_matrix = []

        hundred_list_keys = list(hundred_list.keys())
        window = 20

        for i, row_string in enumerate(f):
            row = row_string.split(",")
            go_term = row.pop(0)
            bias = float(row.pop(0))
            empty_string = row.pop(len(row) - 1)
            # Dot-Product of weight vector and input vector
            predictions = []
            # splitting data into training set and test set
            test = hundred_list_keys[window - 20: window]
            test_genes = []
            for key in test:
                gene = Gene(key)
                gene.fitness_levels = hundred_list[key][0]
                test_genes.append(gene)

            for gene in test_genes:
                prediction = 0
                for j, coefficient in enumerate(row):
                    prediction += gene.fitness_levels[j] * float(coefficient)
                prediction -= bias
                if gene.name in go_term_list[current_go_term_index].gene_annotations:
                    actual_label = go_term_list[current_go_term_index].gene_annotations[gene.name][1]
                    predictions.append((prediction, actual_label))
            normalize(predictions)
            confidence_matrix.append(predictions)
            window += 20
            # sensitivity = true_positives / (true_positives + false_negatives)
            # specificity = true_negatives / (true_negatives + false_positives)
        return confidence_matrix
            # if prediction > 0:
            #     print("Gene is annotated with " + go_term + ": " + str(prediction))
            # else:
            #     print("Gene is not annotated with " + go_term + ": " + str(prediction))


def produce_predictions_matrix(file_name, gene_list, go_term_list):
    with open("moeller_predicted_labels.tsv", "w") as output:
        for i, go_term in enumerate(go_term_list):
            output.write(go_term.name)
            if i < len(go_term_list) - 1:
                output.write("\t")
            else:
                output.write("\n")
        models = []
        with open(file_name, 'r') as f:
            for row_string in f:
                row = row_string.split(",")
                go_term = row.pop(0)
                empty_string = row.pop(len(row) - 1)
                model = []
                for coefficient in row:
                    model.append(float(coefficient))
                models.append(model)

        for gene in gene_list:
            print(gene.name)
            predictions = []
            output.write(gene.name + "\t")

            for model in models:
                bias = model.pop(0)
                prediction = 0
                for j, coefficient in enumerate(model):
                    prediction += gene.fitness_levels[j] * coefficient
                prediction -= bias
                model.insert(0, bias)
                predictions.append(prediction)

            normalize(predictions)
            for i, label in enumerate(predictions):
                output.write(str(label))
                if i < len(predictions) - 1:
                    output.write("\t")
                else:
                    output.write("\n")


def normalize(array):
    max_val = float("-inf")
    min_val = float("inf")
    for pair in array:
        max_val = max(pair[0], max_val)
        min_val = min(pair[0], min_val)
    for i, pair in enumerate(array):
        if max_val - min_val != 0:
            array[i] = (pair[0] - min_val) / (max_val - min_val), pair[1]
        else:
            array[i] = 0, pair[1]
    array.sort(key=lambda x: x[0])


def produce_labels():
    num_characteristics, gene_list, go_term_list = process_file()
    produce_predictions_matrix("smo_full_model.weights", gene_list, go_term_list)

#produce_labels()
