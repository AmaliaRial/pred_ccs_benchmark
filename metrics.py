import pandas as pd
import numpy as np
import sys
import os


databases=['ccsbase','ccscompendium','cembio'] #this will depend on the databases I use
tools=['allccs','ccsbase','deepccs','darkchem']#need to add the rest (ccsPred2, HyperCCS, and finally ours once ready)
data= pd.DataFrame(columns=['Dataset','Tool','Mean_abs','SD_abs','Mean_perc','SD_perc','Outliers'])

for db in databases:
    for t in tools:
        df=pd.read_csv('./{}/{}/metrics{}.csv'.format(db,t,t),sep=';')
        data = pd.concat([data, df])
        
data.to_csv('./metrics/joined_metrics.csv', sep=',', index=False)