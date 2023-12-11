import copy 

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

print(rowfirst_labels[0][0])
print(rowfirst_labels[0][1])
print(rowfirst_labels[0][2])
trained = []
skipped = []

def not_in_trained(cell):
    for index1 in range(0, len(trained)):
        if cell == trained[index1]:
            return False
    return True

def find_col(target):
    idxc = 0
    idx1 = 0
    idx2 = 0
    for index1 in range(1, len(rowfirst_labels[0])):
        label_score= 0
        idxc = index1
        for index2 in range(1, len(rowfirst_labels)):
            if rowfirst_labels[index2][index1] == target:
                current_node = index2
                for index3 in range(1, len(rowfirst_labels)):
                    if index2 != index3 and rowfirst_labels[index3][index1] == target and (not_in_trained((index1, index2, index3))):
                        trained.append((index1, index2, index3))
                        return (rowfirst_labels[index2][0], rowfirst_labels[index3][0], rowfirst_labels[0][index1])
                
    return -1

#for index in range(1, 5):
#    print(find_col(1))
#    print(find_col(-1))
        
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

def c_builder(tbs, target):
    init_max = tbs[len(tbs)-1][1]
    for index1 in range(0, len(tbs)):
        tbs[index1] = (tbs[index1][0], (tbs[index1][1]/init_max)*target)
    return tbs

Gol = {}

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
                #print("OLD: ", olist[i3])
                olist[i3] = ttup
                #print("NEW: ", olist[i3])
            else:
                tx = (olist[i3][0], 0)
                olist[i3] = tx
        if target == -1:
            adj = i2 - i3
            adj = adj*0.005*-1
            if olist[i3][1] + adj <= 0:
                tp = olist[i3][1] + adj
                ttup = (olist[i3][0], tp)
                #print("OLD: ", olist[i3])
                olist[i3] = ttup
                #print("NEW: ", olist[i3])
            else:
                tx = (olist[i3][0], 0)
                olist[i3] = tx
    return olist
def multip(gn, weights):
    idx1 = 0
    total = 0
    for index1 in range(1, len(rowfirst_screening)):
        if rowfirst_screening[index1][0] == gn:
            idx1 = index1
            break
    for index1 in range(1, len(rowfirst_screening[0])):
        for index2 in range(0, len(weights)):
            if weights[index2][0] == rowfirst_screening[0][index1]:
                temp = weights[index2][1]*rowfirst_screening[idx1][index1]
                
                total = total + temp

    return total
Goln = {}
def svm(t_matrix, p_matrix, target):
    Counter = 0
    for index1 in range(1, len(t_matrix[0])):
        Gol[t_matrix[0][index1]] = []
        Goln[t_matrix[0][index1]] = []
    while(True):
        tp = find_col(target)
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
        if Counter >= 1:
            break
    target = target*-1
    while(True):
        tp = find_col(target)
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
        if Counter >= 3:
            break
    print(t_matrix[606][0], t_matrix[1][0])
    val1 = multip(t_matrix[606][0], Gol['GO:0002576'])
    val2 = multip(t_matrix[1][0], Gol['GO:0002576'])
    #print(Gol['GO:0002576'])
    #print(Goln['GO:0002576'])
    correct = 0
    wrong = 0
    zero = 0
    for index1 in range(1, len(t_matrix)):
        pred = 0
        if t_matrix[index1][1] != 0:
            temps = multip(t_matrix[index1][0], Gol['GO:0002576'])
            if abs(temps - val1) < abs(temps - val1):
                pred = 1
            else:
                pred = -1
            if t_matrix[index1][1] == pred:
                correct = correct + 1
            if t_matrix[index1][1] != pred:
                wrong = wrong + 1
        else:
            zero = zero + 1
    print("CORRECT: ", correct)
    print("WRONG: ", wrong)
    print("ZERO: ", zero)

    print(Gol)
    print("FIN")

svm(rowfirst_labels, [], -1)
