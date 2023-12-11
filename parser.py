rowfirst_labels = []
rowfirst_screening = []

# Read in and split the next line of each file for labels
with open("labels.tsv", "r") as data:
    for i, dd_line in enumerate(data):
        parsed_data = dd_line.split()
        #first row is read in as it contains no integers(header row)
        if i == 0:
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
            rowfirst_screening.append(parsed_data)
        else:
            for index1 in range(1, len(parsed_data)):
                if index1 != 0:
                    parsed_data[index1] = float(parsed_data[index1])
                    
            rowfirst_screening.append(parsed_data)

            
