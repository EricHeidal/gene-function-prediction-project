import math
from scipy import integrate
from operator import itemgetter 

import plotly
plotly.tools.set_credentials_file(username='moell162', api_key='XWeDQVjHzfvEN01WsE7p')

import plotly.plotly as py
import plotly.graph_objs as go

# Create random data with numpy
import numpy as np



#Method used for graphing 5 go terms
def rocPrep(array): #array of tuples (go_term, sensNum, specNum, auc)
    #this is going to create the data object from the array passed in
    data = []
    for k in array:
        goTerm = go.Scatter(
            x = k[2], #specNum
            y = k[1], #sensNum
            name = k[0].name
        )
        data.append(goTerm)

    xaxis = [0, .2, .4, .6, .8, 1]
    yaxis = [0, .2, .4, .6, .8, 1]
    random = go.Scatter(
        x = xaxis,
        y = yaxis,
        name = "random"
    )

    data.append(random)

    return data

def roc(data, name):
    # Create a trace
    layout = go.Layout(
        title=go.layout.Title(
            text=name
            
        ),
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text='1-Specificity'
            )
        ),
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text='Sensitivity'
            )
        )
    )

    fig = go.Figure(data=data, layout=layout)
    py.plot(fig, filename=name)

#Finds best, worst and middle go term spec and sens scores and graphs them
def bestWorstMiddle(array): #array of tuples (go_term, sensNum, specNum)

    auc_list = []

    #Calculate AUC for every GoTerm
    for element in array:
        auc = integrate.simps(y=element[1], x=element[2])
        auc_list.append((element[0], element[1], element[2], auc))

    #sort the array (on the last element of each tuple which is auc)
    auc_list.sort(key=lambda tup: tup[3])

    #take the first 5 for the worst
    worst = auc_list[0:6]

    #take the last 5 for best
    best = auc_list[len(auc_list) - 5: len(auc_list)]

    mid = len(auc_list) / 2
    middle = auc_list[int(mid): int(mid) + 6]

    #pass it over to roc to create data that will be graphed, three graphs will be made
    data = rocPrep(worst)
    roc(data, "Worst 5 Go Terms based on AUC")

    data2 = rocPrep(middle)
    roc(data2, "Middle 5 Go Tems based on AUC")

    data3 = rocPrep(best)
    roc(data3, "Top 5 Go Tems based on AUC")

