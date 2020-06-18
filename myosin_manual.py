#Script that analyses the manual tracks of myosins
#@mcolomerrosell

"""IMPORTING LIBRARIES"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib



fig, axes = plt.subplots(nrows=2, ncols=3,
                           sharex=True, sharey=True,
                           figsize=(9,6))
sns.set_palette("bright")
background_noise = 0
length_threshold = 20

tracksIDs = [['01','01'],['01','02'],['01','03'],['01','04'],['02','05'],['02','06']]

background = 0
for i in range(len(tracksIDs)):
    number = i
    video = tracksIDs[i][0]
    number = tracksIDs[i][1]
    path = '/Volumes/SLN_ICFO_Mariona/TIRF_data/BF_TRACKS'
    data_spots = pd.read_csv(path+'/'+'VIDEO'+video+'_SPOTS'+number+'.csv')
    data_tracks = pd.read_csv(path+'/'+'VIDEO'+video+'_TRACKS'+number+'.csv')
    data_background = pd.read_csv(path+'/'+'VIDEO'+video+'_BG'+number+'.csv')
    data_merge = pd.read_csv(path+'/'+'VIDEO'+video+'_MERGE'+number+'.csv')

    #data_spots = largeTracks(data_spots, data_tracks, length_threshold)
    num_tracks = data_spots['TRACK_ID'].nunique()
    tracks = list(data_spots['TRACK_ID'].unique())
    print('Number of tracks', num_tracks)
    #msdPlot(data_spots, length_threshold, num_tracks, tracks, 'hola','green')
    #intensityAnalysisSeabornBackground(data_spots, data_background,2, 'hola')
    #intensityAnalysisSeabornMerge(data_merge, data_background,2, 'hola')
    if background == 0:
        data_spots['LOCAL_BACKGROUND'] = 0
    elif background == 1:
        data_spots['LOCAL_BACKGROUND'] = float(data_background)
    else:
        data_background['FRAME'] = data_background.index+1
        def rowFunc(row):
            bg_value = float(data_background.loc[data_background['FRAME']==row['FRAME'],'Mean'])
            return bg_value
        data_spots['LOCAL_BACKGROUND'] = data_spots.apply(rowFunc, axis=1)
        

        
    data_merge['FRAME'] = data_merge.index+1
    
    if background == 0:
        data_merge['LOCAL_BACKGROUND'] = 0
    elif background == 1:
        data_merge['LOCAL_BACKGROUND'] = float(data_background)
    else:
        data_merge['LOCAL_BACKGROUND'] = data_background['Mean']

    data_merge['INTENSITY'] =  (data_merge['Mean'] - data_merge['LOCAL_BACKGROUND'])/1
    time_conversion = 0.5
    minimum_frame = min(data_spots['FRAME'])
    data_merge['FRAME'] = time_conversion*(data_merge['FRAME']-min(data_spots['FRAME']))
    data_spots['FRAME'] = time_conversion*(data_spots['FRAME']-min(data_spots['FRAME']))
    
    def rowFuncMerge(row):
            merge_value = float(data_merge.loc[data_merge['FRAME']==row['FRAME'],'Mean'])
            return merge_value
        
    data_spots['MERGE_INTENSITY'] = data_spots.apply(rowFuncMerge, axis=1)

    maximum_frame = max(data_spots['FRAME'])
    trackIDs = data_spots.TRACK_ID.unique()

    def typeTrack(c):
        if c['TRACK_ID'] == trackIDs[0]:
            return 'Spot A'
        else:
            return 'Spot B'
    
    data_spots['TRACK_NAME'] = data_spots.apply(typeTrack, axis=1)
    
    colorDictionary = {'Track 1': 'coral', 'Track 2': 'red'}

    #data_spots['LOCAL_BACKGROUND'] = data_background.loc[data_background['FRAME']==data_spots['FRAME'],'Mean']
    data_spots.head()
    data_spots['INTENSITY'] =  (data_spots['MEAN_INTENSITY'] - data_spots['LOCAL_BACKGROUND'])/1
    #plt.subplot(2, 3, number+1)
    ax=axes[i//3, i%3]

    #sns.set_palette("pastel")
    sns.set_palette("bright")
    sns.lineplot(x="FRAME", y="INTENSITY", hue="TRACK_NAME",palette=sns.color_palette('bright', n_colors=2), data=data_spots, estimator = None,  ax=ax)
    sns.lineplot(x="FRAME", y="MERGE_INTENSITY",data=data_spots, label= 'Merged spot',color='green',estimator = None,  ax=ax)
    ax.legend('', frameon=False)
    ax.set_xlim(0,20)
    ax.set_ylim(0,max(data_spots['MEAN_INTENSITY']))
    ax.set_ylim(0,3000)
    ax.set_xlabel('t (s)')
    ax.set_ylabel('I')
fig.legend(title='', labels=['Spot 1', 'Spot 2','Merged spot'],loc='upper center')
fig.savefig('Intensity_Single_Tracks', bbox_inches='tight', dpi=600)
