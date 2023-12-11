import copy 
import csv

rowfirst_labels = []
rowfirst_screening = []

# Read in and split the next line of each file for labels
with open("labels.tsv", "r") as data:
    for i, dd_line in enumerate(data):
        parsed_data = dd_line.split()
        #first row is read in as it contains no integers(header row)
        if i == 0:
            parsed_data = ["gene\GO_Terms"] + parsed_data
            rowfirst_labels.append(parsed_data)
        else:
            for index1 in range(1, len(parsed_data)):
                #the first item in every row is a header, the rest are converted into floats before being read into the list
                if index1 != 0:
                    parsed_data[index1] = float(parsed_data[index1])
                    
            rowfirst_labels.append(parsed_data)

# Read in and split the next line of each file for screening
with open("screening_data.tsv", "r") as data:
    for i, dd_line in enumerate(data):
        parsed_data = dd_line.split()
        if i == 0:
            parsed_data = ["gene\cell_lines"] + parsed_data
            rowfirst_screening.append(parsed_data)
        else:
            for index1 in range(1, len(parsed_data)):
                if index1 != 0:
                    parsed_data[index1] = float(parsed_data[index1])
                    
            rowfirst_screening.append(parsed_data)

trained = []
skipped = []

#Function to check that a pair  has  not been used for training
def not_in_trained(cell):
    for index1 in range(0, len(trained)):
        if cell == trained[index1]:
            return False
    return True

#Finds a pair on a GO Term that hasn't been used for training
def find_col(target, col):
    idxc = 0
    idx1 = 0
    idx2 = 0
    for index1 in range(1, len(rowfirst_labels[0])):
        label_score= 0
        if col == -1:
            idxc = index1
        else:
            idxc = col
        for index2 in range(1, len(rowfirst_labels)):
            if rowfirst_labels[index2][idxc] == target:
                current_node = index2
                for index3 in range(1, len(rowfirst_labels)):
                    if index2 != index3 and rowfirst_labels[index3][idxc] == target and (not_in_trained((idxc, index2, index3))):
                        trained.append((idxc, index2, index3))
                        return (rowfirst_labels[index2][0], rowfirst_labels[index3][0], rowfirst_labels[0][idxc])               
    return -1


#Get the difference between two genes' expression values      
def diff(g1, g2):
    idx1 = 0
    idx2 = 0
    diff_d = {}
    diff_l = []
    i1 = 0
    i2 = 0
    for index1 in range(1, len(rowfirst_screening)):
        if rowfirst_screening[index1][0] == g1:
            i1 = index1
        if rowfirst_screening[index1][0] == g2:
            i2 = index1
    for index1 in range(1, len(rowfirst_screening[0])):
        diff = ()
        diff = (rowfirst_screening[0][index1], rowfirst_screening[i1][index1] - rowfirst_screening[i2][index1])
        if diff[1] < 0:
            diff = (diff[0], diff[1]*-1)
        diff_l.append(diff)
        diff_d[diff] = index1
    diff_l_r = []
    diff_l_r = sorted(diff_l, key=lambda x: x[1])
    return(diff_l_r)

#This function normalizes the difference between two genes' expression values
def c_builder(tbs, target):
    init_max = tbs[len(tbs)-1][1]
    for index1 in range(0, len(tbs)):
        tbs[index1] = (tbs[index1][0], (tbs[index1][1]/init_max)*target)
    return tbs


#Adjuster function to change weights during training
def adjuster(nlist, olist, target):
    i2 = 0
    i3 = 0
    test_list = copy.copy(olist)
    for index1 in range(1, 114):
        s_ind = "gene" + str(index1)
        for index2 in range(0, len(nlist)):
            if nlist[index2][0] == s_ind:
                i2 = index2
                break;
        for index3 in range(0, len(olist)):
            if olist[index3][0] == s_ind:
                i3 = index3
                break;

        if target != -1:
            adj = i2 - i3
            adj = adj*0.005
            if olist[i3][1] + adj >= 0:
                tp = olist[i3][1] + adj
                ttup = (olist[i3][0], tp)
                olist[i3] = ttup
            else:
                tx = (olist[i3][0], 0)
                olist[i3] = tx
        if target == -1:
            adj = i2 - i3
            adj = adj*0.005*-1
            if olist[i3][1] + adj <= 0:
                tp = olist[i3][1] + adj
                ttup = (olist[i3][0], tp)
                olist[i3] = ttup
            else:
                tx = (olist[i3][0], 0)
                olist[i3] = tx
    return olist

#Multiplier Function used to multiply expression values of a gene with both negative and positive weights 
def multip(gn, pweights, nweights):
    idx1 = 0
    ptotal = 0
    ntotal = 0
    for index1 in range(1, len(rowfirst_screening)):
        if rowfirst_screening[index1][0] == gn:
            idx1 = index1
            break
    if idx1 == 0:
        print("FAILED")
    for index1 in range(1, len(rowfirst_screening[0])):
        for index2 in range(0, len(pweights)):
            if pweights[index2][0] == rowfirst_screening[0][index1]:
                temp = pweights[index2][1]*rowfirst_screening[idx1][index1]
                ptotal = ptotal + temp
    for index1 in range(1, len(rowfirst_screening[0])):
        for index2 in range(0, len(nweights)):
            if nweights[index2][0] == rowfirst_screening[0][index1]:
                temp = nweights[index2][1]*rowfirst_screening[idx1][index1]
                
                ntotal = ntotal + temp
    x = (ptotal, ntotal)
    return x


Gol = {}
Goln = {}
PointP = {}
PointN = {}
Conf = {}

#SVM-Like ALgorithm to process the data and make continuous predictions
#t_matrix is the matrix that the algorithm will train and make predictions
#on, the target should be left as int "-1" and d is the number of columns
#that you want to analyze from the start. WARNING: a big number for d may
#increase the computation time significantly
#Example use of svm: svm(rowfirst_labels, -1,  3)

def svm(t_matrix, target, d):
    Counter = 0
    for index1 in range(1, len(t_matrix[0])):
        PointP[t_matrix[0][index1]] = []
        PointN[t_matrix[0][index1]] = []
        Gol[t_matrix[0][index1]] = []
        Goln[t_matrix[0][index1]] = []
        Conf[t_matrix[0][index1]] = []
    #Weight training and initialization
    for index_top in range(1, d):
        while(True):
            tp = find_col(target, index_top)
            if tp == -1:
                break
            lt = diff(tp[0], tp[1])
            ltf = c_builder(lt, target)
            Counter = Counter + 1
            if Gol[tp[2]] == []:
                Gol[tp[2]] = ltf
            else:
                old = Gol[tp[2]]
                Gol[tp[2]] = adjuster(ltf, Gol[tp[2]], target)
            if Counter >= 500:
                break
        target = target*-1
        while(True):
            tp = find_col(target, index_top)
            if tp == -1:
                break
            lt = diff(tp[0], tp[1])
            ltf = c_builder(lt, target)
            Counter = Counter + 1
            if Goln[tp[2]] == []:
                Goln[tp[2]] = ltf
            else:
                old = Goln[tp[2]]
                Goln[tp[2]] = adjuster(ltf, Goln[tp[2]], target)
            if Counter >= 2000:
                break
	#Initializing the point values for each GO Term
        for index1 in range(1, len(t_matrix[0])):
            cpos = 0
            cneg = 0
            vpos = 0
            vneg = 0
            for index2 in range(1, len(t_matrix)):
                if t_matrix[index2][index1] == 1 and Gol[t_matrix[0][index1]] != []:
                    cpos = cpos + 1
                    tmp = multip(t_matrix[index2][0], Gol[t_matrix[0][index1]], Goln[t_matrix[0][index1]])
                    vpos = vpos + tmp[0]
                if t_matrix[index2][index1] == -1 and Goln[t_matrix[0][index1]] != []:
                    cneg = cneg + 1
                    tmp = multip(t_matrix[index2][0], Gol[t_matrix[0][index1]], Goln[t_matrix[0][index1]])
                    vneg = vneg + tmp[1]
            if cpos + cneg != 0:
                PointP[t_matrix[0][index1]] = vpos/cpos
                PointN[t_matrix[0][index1]] = vneg/cneg
            

    #Make predictions about data in the already trained set and use
    #the information to readjust point values
    for index_top in range(1, d):
        correct = 0
        wrong = 0
        zero = 0
        P1 = 0
        CP1 = 0
        MIS1 = 0
        for index1 in range(1, len(t_matrix)):
            itr = 2
            itr1 = itr -1
            for index in range(1, itr):
                pred = 0
                curr = t_matrix[0][index_top]
                if t_matrix[index1][index_top] != 0:
                    temps = multip(t_matrix[index1][0], Gol[curr], Goln[curr])
                    if PointP[curr] != [] and PointN[curr] != []:
                        if abs(temps[0] - PointP[curr]) <= abs(temps[0] - PointN[curr]):
                            pred = -1
                        else:
                            if index >= itr1:
                                P1 = P1 + 1
                            pred = 1
                        if t_matrix[index1][index_top] == pred and index >= itr1:
                            correct = correct + 1
                            if pred == 1 and index >= itr1:
                                correct = correct + 1
                                CP1 = CP1 + 1
                        if t_matrix[index1][index_top] != pred:
                            if index >= itr1:
                                
                                wrong = wrong + 1
                            if pred == 1 and index >= itr1:
                                temp_v = (PointP[curr]-temps[0])/(100)
                                PointN[curr] = PointN[curr] + temp_v
                            if pred == -1 and index >= itr1:
                                MIS1 = MIS1 + 1
                                temp_v = (PointN[curr]-temps[1])/(100)
                                PointP[curr] = PointP[curr] + temp_v
                                
                elif index >= 5:
                    zero = zero + 1
        
        print("GO TERM: ", t_matrix[0][index_top])
        print("CORRECT: ", correct)
        print("WRONG: ", wrong)
        print("ZERO: ", zero)
        print("FP: ", P1-CP1)
        print("TP: ", CP1)
        print("ALL POSITIVES: ", MIS1+CP1)
    
    #Make predictions
    for index1 in range(1, d):
        curr = t_matrix[0][index1]
        max1 = -99999999999999999999999999999999999999
        min1 = 99999999999999999999999999999999999999
        bench = 0
        for index2 in range(1, len(t_matrix)):
            temps = multip(t_matrix[index2][0], Gol[curr], Goln[curr])
            Conf[t_matrix[0][index1]].append(temps[0])
            if temps[0] > max1:
                max1 = temps[0]
            if temps[0] < min1:
                min1 = temps[0]
            if PointP[curr] != [] and PointN[curr] != []:
                if abs(temps[0] - PointP[curr]) <= abs(temps[0] - PointN[curr]):
                    pred = -1
                else:
                    pred = 1
        if Gol[curr] > Goln[curr]:
            bench = max1
        else:
            bench = min1
        for index3 in range(0, len(Conf[t_matrix[0][index1]])):
            Conf[t_matrix[0][index1]][index3] = Conf[t_matrix[0][index1]][index3]/bench
            if Conf[t_matrix[0][index1]][index3] < 0:
                Conf[t_matrix[0][index1]][index3] = Conf[t_matrix[0][index1]][index3]*-1
    
    #write the prediction confidence scores in tab separated form
    with open('output.tsv', 'w', newline = '') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        first = []
        for index1 in range(0, len(t_matrix)):
            first.append(t_matrix[index1][0])
        tsv_output.writerow(first)
        second = []
        for index1 in range (1, 3):
            second = []
            for index2 in range(0, len(Conf[t_matrix[0][index1]])):
                if index2 == 0:
                    second.append(t_matrix[0][index1])
                else:
                    second.append(Conf[t_matrix[0][index1]][index2])
            tsv_output.writerow(second)
                            
svm(rowfirst_labels, -1,  6)


