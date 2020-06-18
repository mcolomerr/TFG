#Script that process the Cytosim output files of myosin tracks
#and analysis the diffusion movement of the particles
#and the contraction of the network 
#@mcolomerrosell


#Importing libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.cluster import DBSCAN
from scipy import stats

import matplotlib.pyplot as plt
import matplotlib
import os
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

def open_filename(filename):
    with open(filename, 'r') as file:
        lines = file.read().splitlines()

    colnames = ['frame'] + [i for i in lines[4].split(' ') if i not in ['','%']]

    lines = lines[5:]
    df = pd.DataFrame(columns=colnames)

    frame_counter = 0
    line_counter = 0
    for i in range(len(lines)):
        if len(lines[i]) > 0 and '%' not in lines[i]: #miro si son dades o un separador
            data = [j for j in lines[i].split(' ') if j != '']
            df.loc[line_counter] = [int(frame_counter)]+data
            line_counter +=1 
        else:
            #es un separador i sempre que hi ha un separador hi ha 6 linies que thas de saltar
            frame_counter +=1/6
    return df

def open_filename(filename):
    exclude = [i for i, line in enumerate(open(filename)) if line.startswith('%')]
    data = pd.read_csv(filename,  delim_whitespace=True , skiprows = exclude[1:],
                       names=['class', 'identity', 'active', 'posX', 'posY', 'fiber1', 
                             'abscissa1', 'fiber2', 'abscissa2'], header=None)
    data = data.iloc[1:]
    data.reset_index()
    data.loc[:,'identity'] = data.identity.astype(np.float)
    nMolecules = max(data['identity'])
    data['frame'] = (data.index-1)//nMolecules
    data.loc[:,'frame'] = data.frame.astype(np.int)
    data['radius'] = np.sqrt(data['posX']**2+data['posY']**2)
    return data, nMolecules

def position(data, frame):
    """Intensity vs time plot"""
    data_frame = data[data['frame']==frame]
    print(len(data_frame))
    fig = plt.figure(figsize = (6,6))
    sns.scatterplot(x="posX", y="posY", data=data_frame)    
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    
def obtainTracks(dataSpots):
    """Obtain the x and y positions from a tracks datafarame,
    returns a vector x and a vector y with the coordinates"""
    tracks = list(dataSpots['identity'].unique())
    xTracks = []
    yTracks = []
    for trackID in tracks:
        xTrackID = []
        yTrackID = []
        xTracks.append(dataSpots[dataSpots['identity']==trackID]["posX"].values)
        yTracks.append(dataSpots[dataSpots['identity']==trackID]["posY"].values)
    return xTracks, yTracks

def obtainTracksFrame(dataSpots, frame):
    """Obtain the x and y positions from a tracks datafarame,
    returns a vector x and a vector y with the coordinates"""
    dataFrame = dataSpots[dataSpots['frame']==frame]
    tracks = list(dataFrame['identity'].unique())
    xTracks = []
    yTracks = []
    for trackID in tracks:
        xTrackID = []
        yTrackID = []
        xTracks.append(dataFrame[dataFrame['identity']==trackID]["posX"].values)
        yTracks.append(dataFrame[dataFrame['identity']==trackID]["posY"].values)
    return xTracks, yTracks
  
def msdPlot(df1, time_length, labelMolecule):
    num_tracks = df1['identity'].nunique()
    time_list = np.arange(1,time_length,1)
    x_pos = []
    y_pos = []
    df1['NewIndex'] = df1.groupby(['identity']).ngroup()
    df1.head()
    for m in range(num_tracks+1):
        x_pos.append([])
        y_pos.append([])
    for i in range(len(df1['posX'])):
        index = int(df1.iloc[i]['NewIndex'])
        x_pos[index].append(df1.iloc[i]["posX"]) 
        y_pos[index].append(df1.iloc[i]["posY"])
    
    #Computing the average msd
    average = []
    msd_time = []
    time = []
    time_list = np.arange(1,int(time_length/2),1)
    for s in range(num_tracks):
        if len(x_pos[s])>=time_length:
            msd_cell_i = []
            for m in time_list:
                displ0 = []
                for i in range((len(x_pos[s])-m)):
                    x1 = x_pos[s][i]
                    x2 = x_pos[s][i+m]
                    y1 = y_pos[s][i]
                    y2 = y_pos[s][i+m]
                    dis = ((x2-x1)**2+(y2-y1)**2)
                    displ0.append(dis)
                msd_cell_i.append(np.average(displ0))
            msd_time.append(msd_cell_i)
            time.append(time_list)
            
    #Plot MSD-t
    ave = []
    for t in range(len(time_list)):
        mean = []
        for m in range(len(msd_time)):
            mean.append(msd_time[m][t])
        ave.append(np.nanmean(mean))
    
    #Furth equation
    from scipy.optimize import curve_fit
    def f2(t, A, C): 
        return C*t**A
    A,C = curve_fit(f2, time_list, ave)[0] # your data x, y to fit
    pcov = curve_fit(f2, time_list, ave)[1] # your data x, y to fit
    errorA = np.sqrt(np.diag(pcov))[0]
    print(A, errorA)
    
    d = {'TIME':np.reshape(time,-1),'MSD':np.reshape(msd_time,-1)}
    df = pd.DataFrame(d)
    df['Nmolec'] = labelMolecule
    fig = plt.figure(3, figsize=(3,3))
    sns.lineplot(x="TIME", y="MSD", ci=None,data=df, label=labelMolecule)
    plt.xlim(0,time_length/2)
    plt.ylim(0,4)
    plt.legend(loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1))
    plt.xlabel('$\\tau $ (s)')
    plt.ylabel('MSD ($\mu m$)')
    fig.savefig('MSD_cytosim.png', bbox_inches='tight', dpi=600)
    return df

def contractility(data):
    data['sqRadius'] = data['radiusAverage']**2
    avRadius =  data['sqRadius'].groupby([data['frame']]).rolling(window=2).mean()
    effectiveRadius = np.sqrt(2*avRadius)
    return effectiveRadius



path = '/Volumes/SLN_ICFO_Mariona/cytosim/TFG_report10/'
runs = ['run0000', 'run0001','run0002','run0003','run0004']



#runs = ['run0000', 'run0001']

frames = [0,25,50,75,100]

fig, axes = plt.subplots(nrows=len(runs), ncols=len(frames),
                           sharex=True, sharey=True,
                           figsize=(2*(len(frames)),2*(len(runs))))

name = 'diffConc_'
cols = ['t = 0s', 't = 2.5s','t = 5s','t = 7.5s', 't = 10s']
rows = []
xlabel = []

#df = pd.DataFrame()

for i in range(len(runs)):
    filename = path+runs[i]+'/myosin.txt'
    data, nMolecules = open_filename(filename)
    data.loc[:,'posX'] = data.posX.astype(np.float)
    data.loc[:,'posY'] = data.posY.astype(np.float)
    for j in range(len(frames)):
        frame = frames[j]
        ax=axes[i, j]
        data_frame = data[data['frame']==frame]
        sns.scatterplot(x="posX", y="posY", data=data_frame, ax = ax, alpha=0.25)    
        plt.xlim(-5,5)
        plt.ylim(-5,5)
    print(data.tail())
    labelMolec = 'N = '+str(int(nMolecules))+'\n y ($\mu m$)'
    data['radiusAverage'] = data['radius'].groupby(data['identity']).transform('mean')
    dataBleb = data[data['radiusAverage']>2]
    dataNucleus = data[data['radiusAverage']<1.5]
    radius = contractility(data)
    print(radius)
    #data3 = data[data['frame']<201]
    #newDf = msdPlot(dataBleb, 100, labelMolec)
    #newDf = msdPlot(dataNucleus, 100, labelMolec)
    df = df.append(newDf)
    rows.append(labelMolec)
    xlabel.append('x ($\mu m$)')
    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.set_yticks([-4,0,4])
    ax.set_xticks([-4,0,4])
plt.subplots_adjust(wspace=0, hspace=0)
    
xlabel = 'x ($\mu m$)'
for i, ax in enumerate(axes.flatten()[:5]):
    ax.set_title(cols[i], fontweight='bold')
    
for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel(row, rotation=90, size='large')
    
for ax in axes[4,:]:
    ax.set_xlabel(xlabel, size='large')

fig.savefig(name+'myosin_cytosim.png', bbox_inches='tight', dpi=400)
