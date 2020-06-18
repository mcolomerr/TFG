#Script that analyses the output files of
#trackmate plugin
#@mcolomerrosell


"""IMPORTING LIBRARIES"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#Options
params = {'text.usetex' : True,
          'font.size' : 11,
          'font.family' : 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params) 

"""TRACKS SELECTION"""
def largeTracks(data_spots, data_tracks, length_threshold):
    """Threshold the Tracks for their length"""
    low_numberspots = []
    for i in range(len(data_tracks)):
        if data_tracks.NUMBER_SPOTS[i] <=length_threshold:
            label = int((data_tracks.Label[i].split("_"))[1])
            low_numberspots.append(label)
    indexNames = []
    for i in range(len(data_spots)):
        if data_spots.TRACK_ID[i] in low_numberspots:
            indexNames.append(i)
    data_spots.drop(indexNames , inplace=True)
    return data_spots
    
"""INTENSITY FEATURES"""
def intensityAnalysis(data_spots, length_threshold):
    """Intensity vs time plot"""
    
    intensity = []
    frames = []
    for m in range((num_tracks)):
        intensity.append([])
        frames.append([])

    for i in range(len(data_spots)):
        m = data_spots.iloc[i]['TRACK_ID']  
        index = tracks.index(m)
        intensity[index].append(data_spots.iloc[i]['MEAN_INTENSITY']-background_noise)
        frames[index].append(data_spots.iloc[i]['FRAME'])
    fig1 = plt.figure(1, figsize = (6,6))
    max_intensity = np.max(intensity)

    for m in range((num_tracks)):
        plt.plot(intensity[m], linewidth = 0.5)
    plt.xlabel('Time (frames)')
    plt.ylabel('Mean intensity')
    fig1.savefig('Intensity vs time')

    fig2 = plt.figure(2, figsize = (6,6))
    intensity_ave = []
    for i in range(length_threshold):
        ave = []
        for m in range(len(intensity)):
            ave.append(intensity[m][i])
        intensity_ave.append(np.average(ave))
    plt.plot(intensity_ave)
    plt.xlabel('Time (frames)')
    plt.ylabel('Mean intensity')
    plt.xlim([0,length_threshold])
    fig2.savefig('Intensity vs time - average')
    
"""INTENSITY FEATURES"""
def intensityAnalysisSeaborn(data_spots, length_threshold, background_noise, name):
    """Intensity vs time plot"""
    time_conversion = 0.5
    data_spots["Time zero"] = time_conversion*(data_spots.groupby("TRACK_ID").cumcount())
    data_spots['INTENSITY'] =  data_spots['MEAN_INTENSITY'] - background_noise 
    fig = plt.figure(1, figsize = (6,6))
    sns.lineplot(x="Time zero", y="INTENSITY", data=data_spots)
    plt.xlim(0,30)
    plt.ylim(0,750)
    plt.xlabel('Time (frames)')
    plt.ylabel('Mean intensity')
    fig.savefig('Intensity'+name, dpi=600)
    
    
def intensityAnalysisSeabornPercentage(data_spots, length_threshold, background_noise, name, number):
    """Intensity vs time plot"""
    time_conversion = 0.5
    data_spots["Time zero"] = time_conversion*(data_spots.groupby("TRACK_ID").cumcount())
    maximumIntensity = data_spots['MEAN_INTENSITY'].max() - background_noise
    data_spots['INTENSITY'] =  (data_spots['MEAN_INTENSITY'] - background_noise)/maximumIntensity
    #palette = sns.color_palette("GnBu_d", 7)
    palette = sns.dark_palette("orange",7)
    fig = plt.figure(1, figsize = (3,3))
    sns.lineplot(x="Time zero", y="INTENSITY", data=data_spots, color=palette[number])
    plt.xlim(0,30)
    plt.ylim(0,1)
    plt.xlabel('t (s)')
    plt.ylabel('I')
    fig.savefig('Intensity.png', bbox_inches='tight', dpi=600)
    
"""INTENSITY FEATURES"""
    
def position(data_spots, length_threshold, name):
    """Intensity vs time plot"""
    data_spots["Time zero"] = data_spots.groupby("TRACK_ID").cumcount()
    data_spots.groupby("TRACK_ID")["POSITION_X"].apply(lambda x:  x-x.iloc[0])
    data_spots['X-X0'] =  data_spots['POSITION_X'] - data_spots.groupby('TRACK_ID')['POSITION_X'].transform('first') 
    data_spots['Y-Y0'] =  data_spots['POSITION_Y'] - data_spots.groupby('TRACK_ID')['POSITION_Y'].transform('first') 
    fig = plt.figure(2, figsize = (6,6))
    sns.lineplot(x="X-X0", y="Y-Y0",hue="TRACK_ID" , data=data_spots, sort=False, lw=0.5)
    plt.xlim(-60,60)
    plt.ylim(-60,60)
    fig.savefig('Track position'+name, dpi=600)
    
def msdPlot(df1, time_length, num_tracks, tracks, name, color, number):
    time_conversion = 0.5
    pixel_conversion = 0.06
    time_list = np.arange(1,time_length,1)
    x_pos = []
    y_pos = []
    df1['NewIndex'] = df1.groupby(['TRACK_ID']).ngroup()
    df1.head()
    for m in range(num_tracks+1):
        x_pos.append([])
        y_pos.append([])
    for i in range(len(df1['POSITION_X'])):
        index = int(df1.iloc[i]['NewIndex'])
        x_pos[index].append(pixel_conversion*df1.iloc[i]["POSITION_X"]) 
        y_pos[index].append(pixel_conversion*df1.iloc[i]["POSITION_Y"])
        
    #Fitting the equation
    from scipy.optimize import curve_fit
    def f(t, D, P): # this is your 'straight line' y=f(x)
        return 4*D*(t-P*(1-np.exp(-t/P)))

    #Computing the average msd
    average = []
    msd_time = []
    time = []
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
            time.append(time_conversion*time_list)
    d = {'TIME':np.reshape(time,-1),'MSD':np.reshape(msd_time,-1)}
    df = pd.DataFrame(d)
                
    #Plot MSD-t
    ave = []
    for t in range(len(time_list)):
        mean = []
        for m in range(len(msd_time)):
            mean.append(msd_time[m][t])
        ave.append(np.nanmean(mean))
    
    #Equation 1
    from scipy.optimize import curve_fit
    def f1(t, A, C): 
        return C*t**A
    A,C = curve_fit(f1, time_list, ave)[0] # your data x, y to fit
    pcov = curve_fit(f1, time_list, ave)[1] # your data x, y to fit
    errorA = np.sqrt(np.diag(pcov))[0]
    errorC = np.sqrt(np.diag(pcov))[1]
    
    print('----------------------------')
    print('Formula 1: <msd> = D*t^a')
    print('Alpha= ', A)
    print('Error alpha= ', errorA)
    print('D1= ', C)
    print('Error D1= ', errorC)
    
    #Plotting the data
    fig = plt.figure(3, figsize=(3,3))
    palette = sns.color_palette("GnBu_d", 7)
    #palette = sns.dark_palette("orange",7)
    sns.lineplot(x="TIME", y="MSD", data=df, label = name, color=palette[number])
    plt.xlim(0,15)
    plt.ylim(0,5)
    plt.xlabel('$\\tau $ (s)')
    plt.ylabel('MSD ($\mu m$)')
    #plt.legend()
    plt.legend('',frameon=False)
    plt.xticks([0,5,10,15])
    fig.savefig('MSD'+name, bbox_inches='tight', dpi=600)
    
    return A, C
    
    
def box_plot(list_persistence, list_labels, feature):
    """This function creates different plots for study the response of cells to temperature.
    The parameters are:
    - list_files:list with the files that you want to analyse
    - number: when using the function more than one time in the same cell, you have to change the number
    - colors: list_colors that will be used"""
    #********************************************************************************
    #LIST OF SETTINGS THAT YOU MAY ADJUST
    time_length = 100 #minimum number of frames that the program will consider for analysing the data
    pixel_conversion = 0.3 #conversion between the pixels and micrometers
    time_step = 1 #number of time steps that we consider 
    velocity_threshold = 0
    msd_threshold = 0
    num_frames = 100 #minimum number of frames that the program will consider for analysing the data
    num_tracks = 30 #number of tracks that it will plot
    
    
    SMALL_SIZE = 14
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    font = {'family' : 'normal',
            'weight' : 'normal'}

    matplotlib.rc('font', **font)
    
    fig2 = plt.figure(figsize=(3,3))
    plt.boxplot(list_persistence, labels=list_labels)
    for i in range(len(list_labels)):
        y = list_persistence[i]
        x = np.random.normal(1+i, 0.04, size=len(y))
        plt.plot(x, y, 'r.', alpha=0.2)
    plt.ylabel(feature)
    plt.rc('xtick')     
    plt.rc('ytick')
    fig2.savefig('diff', bbox_inches='tight', dpi=600)

