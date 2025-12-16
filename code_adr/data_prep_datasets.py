import pandas as pd
import numpy as np
import sys
import os

directory= sys.argv[1]
df=pd.read_csv('./{}/dataset_{}.csv'.format(directory, directory), sep=';')
classifications= pd.read_csv('./{}/classifications.csv'.format(directory),
sep=';')
df=pd.merge(df, classifications, on='SMILES')
df.to_csv('./{}/data_classified.csv'.format(directory), sep=';', index=False)
df2 = df[df['SMILES'].notna()]
df3 = df2.drop_duplicates('SMILES')

#---------------- CSS Base ---------------
adduct= np.repeat('All', df3.shape[0])
data_ccsbase=np.column_stack((adduct, df3.SMILES, df3.SMILES))
column_values_ccsbase = ['Adduct','Smiles','Name']
ccsbase = pd.DataFrame(data = data_ccsbase, columns = column_values_ccsbase)
if not os.path.exists('./{}/ccsbase'.format(directory)):
    os.makedirs('./{}/ccsbase'.format(directory))
ccsbase.to_csv('./{}/ccsbase/dataccsbase.csv'.format(directory), index=False)

#---------------- All CCS ----------------
number= np.arange(df3.shape[0])
data_allccs= np.column_stack((number, df3.SMILES))
allccs= pd.DataFrame(data = data_allccs)
if not os.path.exists('./{}/allccs'.format(directory)):
    os.makedirs('./{}/allccs'.format(directory))
allccs.to_csv('./{}/allccs/dataallccs.csv'.format(directory), index=False,
header=None)

#---------------- Dark Chem ----------------
data_darkchem= df3.SMILES
darkchem = pd.DataFrame(data = df3.SMILES)
darkchem= darkchem.rename(columns={'SMI': 'SMILES'})
if not os.path.exists('./{}/darkchem'.format(directory)):
    os.makedirs('./{}/darkchem'.format(directory))
darkchem.to_csv('./{}/darkchem/datadarkchem.tsv'.format(directory), sep='\t',index=False)

#---------------- Deep CCS ----------------
df4=df3.loc[df3.index.repeat(4)]
add=np.array(["M+H", "M+Na", "M-H", "M-2H"])
adducts=np.tile(add,len(df4)//4)
df4['Adducts']=adducts
data_deepccs=np.column_stack((df4.SMILES, df4.Adducts))
column_values_deepccs = ['SMILES','Adducts']
deepccs = pd.DataFrame(data = data_deepccs,
columns = column_values_deepccs)
if not os.path.exists('./{}/deepccs'.format(directory)):
    os.makedirs('./{}/deepccs'.format(directory))
deepccs.to_csv('./{}/deepccs/datadeepccs.csv'.format(directory), index=False)

#---------------- CCSPred2 ----------------
ccspred2 = pd.DataFrame(
    {'SMILES': df4.SMILES.values, 'Adduct': df4.Adducts.values}
)
if not os.path.exists('./{}/ccspred2'.format(directory)):
    os.makedirs('./{}/ccspred2'.format(directory))
ccspred2.to_csv('./{}/ccspred2/dataccspred2.csv'.format(directory), index=False)

#---------------- Hyper CCS ----------------
hyperccs = pd.DataFrame({'SMILES': df3.SMILES.values})
if not os.path.exists('./{}/hyperccs'.format(directory)):
    os.makedirs('./{}/hyperccs'.format(directory))
hyperccs.to_csv('./{}/hyperccs/datahyperccs.csv'.format(directory), index=False)

#---------------- Our model ----------------
our_model = df3.copy()
if not os.path.exists('./{}/our_model'.format(directory)):
    os.makedirs('./{}/our_model'.format(directory))
our_model.to_csv('./{}/our_model/dataourmodel.csv'.format(directory), sep=';', index=False)
