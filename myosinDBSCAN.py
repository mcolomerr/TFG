#Script that analyses the clustering of single particles
#using the DBSCAN clustering algorithm
#@mcolomerrosell


"""IMPORTING LIBRARIES"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from sklearn.cluster import DBSCAN
from sklearn.metrics.cluster import adjusted_rand_score

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

def obtainTracksFrame(dataSpots, frame):
    """Obtain the x and y positions from a tracks datafarame,
    returns a vector x and a vector y with the coordinates"""
    dataFrame = dataSpots[dataSpots['frame']==frame]
    tracks = list(dataFrame['id'].unique())
    xTracks = []
    yTracks = []
    for trackID in tracks:
        xTrackID = []
        yTrackID = []
        xTracks.append(dataFrame[dataFrame['id']==trackID]["POSITION_X"].values)
        yTracks.append(dataFrame[dataFrame['id']==trackID]["POSITION_Y"].values)
    return xTracks, yTracks

def DBSCAN_clustering(data, frame):
    xTracksFrame, yTracksFrame = obtainTracksFrame(data, frame)
    mergedTracks = np.array([(xTracksFrame[i][0], yTracksFrame[i][0]) for i in range(0, len(xTracksFrame))])
    clustering = DBSCAN(eps=0.4, min_samples=0).fit(mergedTracks)
    labels =clustering.labels_
    #print('nclusters: ', len(set(labels)))
    return mergedTracks, labels




#defining frame number

filename = 'thunderstorm_bleb_offstetBigger0-1_removeBorders.csv'
filename = 'thunderstorm_time.csv'



data = pd.read_csv(filename)
data['POSITION_X'] = data['x [nm]']/1000
data['POSITION_Y'] = data['y [nm]']/1000
dataHeatMap = pd.DataFrame() #creates a new dataframe that's empty

numParticles = []
for i in range(int(max(data['frame']))):
    frame = i+1
    xTracksFrame, yTracksFrame = obtainTracksFrame(data, frame)
    mergedTracks, labels = DBSCAN_clustering(data, frame)

    d = {x:list(labels).count(x) for x in list(labels)}
    dataHist = pd.DataFrame(d.items(), columns=['N', 'Size'])
    dataHist['frame'] = frame
    dataHeatMap = dataHeatMap.append(dataHist)
    numParticles.append(len(xTracksFrame))
    
plt.scatter(mergedTracks[:,0], mergedTracks[:,1], c= labels, cmap='jet')
plt.xlim(20,40)
plt.ylim(10,20)
    
"""
fig = plt.figure(figsize=(4,4))
plt.hist(d.values())
fig.savefig('histogram.png', dpi=200)
"""

fig = plt.figure(figsize=(6,4))
palette = sns.color_palette("GnBu_d")
data2 = dataHeatMap.groupby(["Size", "frame"])["N"].count().reset_index(name="count")
data2['time'] = 0.5*(data2['frame']-1)
dataAll = data2.pivot(index='Size', columns='time',values='count')
sns.heatmap(dataAll, cmap='YlGnBu')
plt.xlabel('t (s)')
fig.savefig('heatmap2.png', bbox_inches='tight', dpi=400)
