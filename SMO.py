import numpy
from SaveModel import *
from SMOFileProcessing import *
#from RocGraph import *
import random


# The SMO Trainer will take in a training set for a single GO term
# Let input training set be a dictionary of gene_name -> ([fitness_levels], annotation)
class SMOTrainer:
    def __init__(self, training_set, num_features):
        # Dictionary of correct classifications/fitness for all ~17000 genes
        self.training_set = training_set
        # Lagrange Multipliers for each training gene
        self.alpha = {}
        for sample in training_set:
            self.alpha[sample] = 0.0
        # Weight vector
        self.w = [0.0] * num_features
        # Bias
        self.b = 0.0
        # Because our data will be imbalanced (i.e. having many more positive/negative annotations compared to the
        # opposite) I am fine tuning our misclassification penalty such that incorrectly classifying a member of the
        # less-frequent class results in a much larger penalty.
        num_negatives = 0
        num_positives = 0
        for sample in training_set:
            if training_set[sample][1] == -1:
                num_negatives += 1
            elif training_set[sample][1] == 1:
                num_positives += 1
        self.C = 0.5
        self.num_non_zero_non_C_alpha = 0
        self.error_cache = {}

    def train(self):
        # num_changed : the number of lagrange multipliers changed on each iteration
        num_changed = 0
        examine_all = True
        steps = 0
        while (num_changed > 0 or examine_all) and steps < 50000:
            steps += 1
            num_changed = 0
            if examine_all:
                for sample in self.training_set:
                    num_changed += self.examine_example(sample)
            else:
                for sample in self.training_set:
                    if self.alpha[sample] != 0 and self.alpha[sample] != self.C:
                        num_changed += self.examine_example(sample)
            if examine_all:
                examine_all = False
            elif num_changed == 0:
                examine_all = True
            if steps % 500 == 0:
                print(steps)
        # Training concludes by returning the newly-trained weight vector and the bias
        return self.w, self.b

    def examine_example(self, sample1):
        # y : the binary classification of a sample (+1 or -1)
        y1 = self.training_set[sample1][1]
        alpha1 = self.alpha[sample1]
        if 0 < alpha1 < self.C:
            e1 = self.error_cache[sample1]
        else:
            e1 = self.svm_predict(self.training_set[sample1][0]) - y1
        r1 = e1 * y1
        # Can decide on a tolerance. 10^-3 is recommended but we can play around with it if we have time
        tolerance = 0.001
        if (r1 < -tolerance and alpha1 < self.C) or (r1 > tolerance and alpha1 > 0):
            if self.num_non_zero_non_C_alpha > 1:
                # The selection of the other alpha we'd like to optimize is chosen by heuristic
                # Tier 1: First, check in error cache for maximum difference in error
                maximum = 0.0
                sample2 = None
                for i in self.error_cache.keys():
                    if 0 < self.alpha[i] < self.C:
                        e2 = self.error_cache[i]
                        temp = abs(e1 - e2)
                        if temp > maximum:
                            maximum = temp
                            sample2 = i
                if sample2 is not None and self.take_step(sample1, sample2) == 1:
                    return 1
            # Tier 2 - If we were unable to make positive progress with the above choice,
            # try any other non-zero, non-C alpha:
            length = len(self.training_set)
            total = 0
            i = random.randrange(length)
            training_set_list = list(self.training_set)
            while total < length:
                total += 1
                if 0 < self.alpha[training_set_list[i % length]] < self.C:
                    if self.take_step(sample1, training_set_list[i % length]) == 1:
                        return 1
                i += 1

            # Un-randomized loop for tier 2 selection
            # for sample1 in self.training_set:
            #     if self.alpha[sample1] != 0 and self.alpha[sample1] != self.C:
            #         if self.take_step(sample1, sample2, e2) == 1:
            #             return 1

            # Tier 3 - If all else fails, try any other sample until we eventually make progress towards optimization
            total = 0
            i = random.randrange(length)
            training_set_list = list(self.training_set)
            while total < length:
                total += 1
                if self.take_step(sample1, training_set_list[i % length]) == 1:
                    return 1
                i += 1

            # Un-randomized version of loop for tier 3 selection
            # for sample1 in self.training_set:
            #     if self.take_step(sample1, sample2, e2) == 1:
            #         return 1
        return 0

    def take_step(self, sample1, sample2):
        if sample1 == sample2:
            return 0
        alpha1_old = self.alpha[sample1]
        alpha2_old = self.alpha[sample2]
        # y : the binary classification of the sample (+1 or -1)
        y1 = self.training_set[sample1][1]
        y2 = self.training_set[sample2][1]

        l, h = self.compute_l_h(y1, y2, alpha1_old, alpha2_old)
        if l == h:
            return 0

        if 0 < self.alpha[sample1] < self.C:
            e1 = self.error_cache[sample1]
        else:
            e1 = self.svm_predict(self.training_set[sample1][0]) - y1
        if 0 < self.alpha[sample2] < self.C:
            e2 = self.error_cache[sample2]
        else:
            e2 = self.svm_predict(self.training_set[sample2][0]) - y2

        k11 = self.kernel(self.training_set[sample1][0], self.training_set[sample1][0])
        k12 = self.kernel(self.training_set[sample1][0], self.training_set[sample2][0])
        k22 = self.kernel(self.training_set[sample2][0], self.training_set[sample2][0])
        eta = (2 * k12) - k11 - k22
        # Can tweak value epsilon here for rounding errors - probably not super important to training, but it is
        # important to time. Platt recommends 10^-3
        epsilon = 0.001
        correct_classification = y1 * y2
        if eta < 0:
            alpha2_new = alpha2_old + (y2 * (e2 - e1) / eta)
            if alpha2_new < l:
                alpha2_new = l
            elif alpha2_new > h:
                alpha2_new = h
        else:
            f1 = eta / 2
            f2 = (y2 * (e1 - e2)) - (eta * alpha2_old)
            l_obj = f1 * l * l + f2 * l
            h_obj = f1 * h * h + f2 * h

            # f1 = y1 * (e1 + self.b) - alpha1_old * k11 - correct_classification * alpha2_old * k12
            # f2 = y2 * (e2 + self.b) - correct_classification * alpha1_old * k12 - alpha2_old * k22
            # l1 = alpha1_old + correct_classification * (alpha2_old - l)
            # h1 = alpha1_old + correct_classification * (alpha2_old * h)
            # # Evaluating objective function at L and H
            # l_obj = l1 * f1 + l * f2 + 0.5 * l1 * l1 * k11 + 0.5 * l * l * k22 + correct_classification * l * l1 * k12
            # h_obj = h1 * f1 + h * f2 + 0.5 * h1 * h1 * k11 + 0.5 * h * h * k22 + correct_classification * h * h1 * k12

            if l_obj > h_obj + epsilon:
                alpha2_new = l
            elif l_obj < h_obj - epsilon:
                alpha2_new = h
            else:
                alpha2_new = alpha2_old
        if abs(alpha2_new - alpha2_old) < epsilon * (alpha2_new + alpha2_old + epsilon):
            return 0
        alpha1_new = alpha1_old - correct_classification * (alpha2_new - alpha2_old)

        if alpha1_new < 0:
            alpha2_new += correct_classification * alpha1_new
            alpha1_new = 0

        elif alpha1_new > self.C:
            t = alpha1_new - self.C
            alpha2_new += correct_classification * t
            alpha1_new = self.C
        # Updating threshold
        if 0 < alpha1_new < self.C:
            new_b = self.b + e1 + y1 * (alpha1_new - alpha1_old) * k11 + y2 * (alpha2_new - alpha2_old) * k12
        else:
            if 0 < alpha2_new < self.C:
                new_b = self.b + e2 + y1 * (alpha1_new - alpha1_old) * k12 + y2 * (alpha2_new - alpha2_old) * k22
            else:
                b1 = e1 + y1 * (alpha1_new - alpha1_old) * k11 + y2 * (alpha2_new - alpha2_old) * k12 + self.b
                b2 = e2 + y1 * (alpha1_new - alpha1_old) * k12 + y2 * (alpha2_new - alpha2_old) * k22 + self.b
                new_b = (b1 + b2) / 2
        delta_b = new_b - self.b
        self.b = new_b
        # Perform vector addition to get a new weight vector
        # Old alphas subtracted from new alphas
        # REMOVE IF NOT USING LINEAR KERNEL
        self.w = numpy.array(self.w) + y1 * (alpha1_new - alpha1_old) * numpy.array(self.training_set[sample1][0]) \
            + y2 * (alpha2_new - alpha2_old) * numpy.array(self.training_set[sample2][0])

        t1 = y1 * (alpha1_new - alpha1_old)
        t2 = y2 * (alpha2_new - alpha2_old)
        for sample in self.training_set:
            if 0 < self.alpha[sample] < self.C:
                self.error_cache[sample] += t1 * self.kernel(self.training_set[sample][0], self.training_set[sample1][0]) + t2 * self.kernel(self.training_set[sample][0], self.training_set[sample2][0]) - delta_b
        self.error_cache[sample1] = 0.0
        self.error_cache[sample2] = 0.0

        self.alpha[sample1] = alpha1_new
        self.alpha[sample2] = alpha2_new
        # This is the only place where the values of alphas are changing, so update the number of alphas that are not 0
        # and are not C here
        if (alpha1_old == 0 or alpha1_old == self.C) and (alpha1_new != 0 and alpha1_new != self.C):
            self.num_non_zero_non_C_alpha += 1
        elif (alpha1_old != 0 and alpha1_old != self.C) and (alpha1_new == 0 or alpha1_new == self.C):
            self.num_non_zero_non_C_alpha -= 1

        if (alpha2_old == 0 or alpha2_old == self.C) and (alpha2_new != 0 and alpha2_new != self.C):
            self.num_non_zero_non_C_alpha += 1
        elif (alpha2_old != 0 and alpha2_old != self.C) and (alpha2_new == 0 or alpha2_new == self.C):
            self.num_non_zero_non_C_alpha -= 1
        return 1

    def svm_predict(self, x):
        return numpy.inner(self.w, x) - self.b

    @staticmethod
    # For now, I'm assuming our classes are linearly separable, so the kernel function is a simple dot-product
    def kernel(x, y):
        return numpy.inner(numpy.array(x), numpy.array(y))

    def compute_l_h(self, y1, y2, alpha1, alpha2):
        correct_classification = y1 * y2
        # C = max(self.C_positive, self.C_negative)
        # Since there are two C values, we need to figure out which is more constrictive and use that as the boundary
        # for alpha2
        if correct_classification == -1:
            if alpha1 - alpha2 > 0:
                return 0, self.C - (alpha1 - alpha2)
            else:
                return -(alpha1 - alpha2), self.C

        elif correct_classification == 1:
            if alpha1 + alpha2 > self.C:
                return alpha1 + alpha2 - self.C, self.C
            else:
                return 0, alpha1 + alpha2
        else:
            exit()


def main():
    start = int(input("Enter index of first GO term in range (start indexing at 0): "))
    end = int(input("Enter index of last GO term in range (highest index is 1213): "))
    if start > end:
        temp = start
        start = end
        end = temp
    file_name = "smo_model(" + str(start) + "-" + str(end) + ").weights"
    # Read in all fitness/annotation data
    print("Processing file")
    num_characteristics, gene_list, go_term_list = process_file()

    for i in range(start, end + 1):
        trainer = SMOTrainer(go_term_list[i].gene_annotations, num_characteristics)
        w, b = trainer.train()
        save_single_weight(file_name, go_term_list[i].name, w, b)

    file_name = "smo_model(0-0).weights"
    predict(file_name, gene_list, go_term_list)

# Compute sensitivity for each threshold
def sensitivityCompute(row):
    #TP / (TP + FN)
    sensArray = []
    #for each threshold (left value in array passed)
    for pair in row:
        tp = 0
        fn = 0
        threshold = pair[0]
        for pair2 in row:
            #if it is annotated
            if(pair2[1] == 1):
                #if confidence is above threshold then True Positive
                if(pair2[0] >= threshold):
                    tp += 1
                #False Negative
                else:
                    fn += 1
        #compute sensitivity
        if(tp == 0 and fn == 0):
            sens = 0
        else:  
            sens = tp / (tp + fn)
        #add to sensitivity array
        sensArray.append(sens)

    return sensArray

# Compute 1 - specificty for each threshold
def one_specificity(row):
    # 1 - (TN / (TN + FP))

    specArray = []
    #for each threshold (left value in array passed)
    for pair in row:
        tn = 0
        fp = 0
        threshold = pair[0]
        for pair2 in row:
            #if it is annotated
            if(pair2[1] == -1):
                #if confidence is above threshold then False positive
                if(pair2[0] >= threshold):
                    fp += 1
                #True negative
                else:
                    tn += 1
        #compute 1 - specificity
        if(tn == 0 and fp == 0):
            spec = 0
        else:  
            spec = 1 - (tn / (tn + fp))
        #add to specificity array
        specArray.append(spec)
    print(specArray)

    return specArray

def computeSensitivitySpecificity(values):
    print("Computing sensitivity and 1 - specificity")

    totalSpec = []
    totalSens = []

    #for each row of values
    for row in values:
        #compute sens and specificity for each row- these methods will each return an array
        sens = sensitivityCompute(row)
        spec = one_specificity(row)
        # add that array to an array of sensitivities and specificities for all iterations of k
        totalSpec.append(spec)
        totalSens.append(sens)

    #Since we need to average the columns, transpose the list
    #average all sensitivities and specificities for each column - this will yield 2 arrays which can be graphed
    averageSpec = []
    print("what is totalSpec", totalSpec)
    for column in zip(*totalSpec): 
        average = sum(column) / len(column)
        averageSpec.append(average)

    averageSens = []
    print("what is totalSens", totalSens)
    for column in zip(*totalSens): 
        average = sum(column) / len(column)
        averageSens.append(average)

    return averageSpec, averageSens



def crossValidation():
    start = int(input("Enter index of first GO term in range (start indexing at 0): "))
    end = int(input("Enter index of last GO term in range (highest index is 1213): "))
    if start > end:
        temp = start
        start = end
        end = temp

    # Read in all fitness/annotation data
    print("Processing file")
    num_characteristics, gene_list, go_term_list = process_file()

    #Execute the following code for each iteration of k-fold cross validation
    gene_list_copy = gene_list.copy()
    # Shuffle the data set randomly
    random.shuffle(gene_list_copy)
    print("length of gene_list_copy ", len(gene_list_copy))

    avgSpOverGo = []
    avgSeOverGo = []
    # call svm method using test and training data
    for current_go_term_index in range(start, end + 1):
        # Pick 100 genes, 50 annotated 50 unannotated
        hundred_list = {}
        unannotatedCount = 0
        annotatedCount = 0
        annotatedList = []
        for i, gene in enumerate(gene_list_copy):

            # if both counts are 50, break from loop
            if (unannotatedCount == 50 and annotatedCount == 50):
                break
            if gene_list_copy[i].name in go_term_list[current_go_term_index].gene_annotations:
                fitness_level = go_term_list[current_go_term_index].gene_annotations[gene_list_copy[i].name][0]
                annotation = go_term_list[current_go_term_index].gene_annotations[gene_list_copy[i].name][1]

                # if annotated count is < 50 and it is annotated, add to list
                if (annotatedCount < 50 and annotation == 1):
                    hundred_list[gene_list_copy[i].name] = (fitness_level, annotation)
                    annotatedList.append((gene_list_copy[i].name,fitness_level, annotation))
                    annotatedCount += 1

                # if unannotated count is < 50 and is unannotated , add to list
                if (unannotatedCount < 50 and annotation == -1):
                    hundred_list[gene_list_copy[i].name] = (fitness_level, annotation)
                    unannotatedCount += 1


        #The special case for the rare go terms that have less than 50 annotations
        #We will use oversampling to reach an annotatedCount of 50
        while(annotatedCount < 50):
            for annot in annotatedList:
                if(annotatedCount < 50):
                    hundred_list["dup" + str(annotatedCount) + annot[0]] = (annot[1],annot[2])
                    annotatedCount += 1
                else:
                    break


        # Split the data set into k groups
        # For each unique group
        window = 20
        for k in range(5):
            hundred_list_keys = list(hundred_list.keys())

            # splitting data into training set
            train = hundred_list_keys[0:window - 20] + hundred_list_keys[window:]

            train_dict = {}
            for key in train:
                train_dict[key] = hundred_list[key]
            file_name = "smo_model(go_term_index=" + str(current_go_term_index) + ").weights"

            trainer = SMOTrainer(train_dict, num_characteristics)
            w, b = trainer.train()
            save_single_weight(file_name, go_term_list[current_go_term_index].name, w, b)
            # incrementing window for test data
            window += 20

        file_name = "smo_model(go_term_index=" + str(current_go_term_index) + ").weights"
        array = predict(file_name, hundred_list, go_term_list, current_go_term_index)

        #Get final array of sensitivities and specifities for a single go term
        averageSpec, averageSens = computeSensitivitySpecificity(array)
        with open("go_term_" + str(current_go_term_index) + "_roc.txt", "a") as roc_file:
            roc_file.write(str((go_term_list[current_go_term_index].name, str(averageSens), str(averageSpec))))
            roc_file.write("\n")
            print(str(current_go_term_index) + " saved")


        #Calculate an average of these arrays for Oguz and add to arrays declared above loop
        #avgSp = sum(averageSpec) / len(averageSpec)
        #avgSe = sum(averageSens) / len(averageSens)
        #avgSeOverGo.append(avgSe)
        #avgSpOverGo.append(avgSp)

        #roc(averageSens, averageSpec, "Go term 1 sensitivity v. 1-specificity")

    #Calculate average sensitivty and specificity over all specified go terms
    #finalSpecAvg = sum(avgSpOverGo) / len(avgSpOverGo)
    #finalSensAvg = sum(avgSeOverGo) / len(avgSeOverGo)

    #print these vals
    #print("Final Spec average over all go terms", finalSpecAvg)
    #print("Final Sens average over all go terms", finalSensAvg)

#main()
crossValidation()
