#Library import
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

#Plot format
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 



def open_filename(filename):
    """Function for opening the files with the right format
    created with the program Cytosim"""
    with open(filename, 'r') as file:
        lines = file.read().splitlines()

    colnames = ['frame'] + [i for i in lines[4].split(' ') if i not in ['','%']]
    
    lines = lines[5:]

    df = pd.DataFrame(columns=colnames)

    frame_counter = 0
    line_counter = 0
    for i in range(len(lines)):
        #print(i)
        if len(lines[i]) > 0 and '%' not in lines[i]: #miro si son dades o un separador
            data = [j for j in lines[i].split(' ') if j != '']
            df.loc[line_counter] = [int(frame_counter)]+data
            line_counter +=1 
        else:
            #es un separador i sempre que hi ha un separador hi ha 6 linies que thas de saltar
            frame_counter +=1/6
    return df

def position(data, frame):
    """Scatter pot of the myosin positions"""
    data_frame = data[data['frame']==frame]
    print(len(data_frame))
    fig = plt.figure(figsize = (6,6))
    sns.scatterplot(x="posX", y="posY", data=data_frame)    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    
   
   
runs = ['run0000', 'run0001', 'run0002','run0003']
cols = ['t = 0s', 't = 5s', 't = 10s']
rows = ['N = 100', 'N = 200', 'N = 300', 'N = 400']


fig, axes = plt.subplots(nrows=4, ncols=3,
                           sharex=True, sharey=True,
                           figsize=(6,8))
    
for i in range(len(runs)):
    run = runs[i]
    filename = run+'_myosin.csv'
    data = pd.read_csv(run+'_myosin.csv')
    data.loc[:,'posX'] = data.posX.astype(np.float)
    data.loc[:,'posY'] = data.posY.astype(np.float)
    ax=axes[i, 0]
    data_frame = data[data['frame']==0]
    sns.scatterplot(x="posX", y="posY", data=data_frame, ax = ax, alpha=0.5)    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    ax=axes[i, 1]
    data_frame = data[data['frame']==50]
    sns.scatterplot(x="posX", y="posY", data=data_frame, ax = ax, alpha=0.5)    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    ax=axes[i, 2]
    data_frame = data[data['frame']==100]
    sns.scatterplot(x="posX", y="posY", data=data_frame, ax = ax, alpha=0.5)    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    
for i, ax in enumerate(axes.flatten()[:3]):
    ax.set_title(cols[i], fontweight='bold')
    
for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel(row, rotation=90, size='large')
    
fig.savefig('myosin_cytosim.png', dpi=200)
  
