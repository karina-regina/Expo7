from scipy.ndimage import rotate, shift
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
import scipy.stats as sp
import pandas as pd
import statsmodels.api as sm
import math
import cmath
import seaborn as sns
import openpyxl

#INPUT FOR DATA ANALYSIS

NumberOfTrials=1
np.random.seed(42)
print(pd.read_excel('LDA_Data.xlsx'))

#SET UP FOR DATA ANALYSIS

RobustnessCheck=0
exportdf={'Distance':[],'# Of Measurements':[],'Time Taken':[],'Error':[],'Failure':[],'Unusable':[],'Angle Of Plume':[],'Plume Origin':[],'Emission Rate':[]}

for experiment in range(NumberOfTrials):

    #INPUT FOR GAS PLUME

    #set from 1-6
    stability = 1                                                                                      
    #Origin Coordinates Of Plume
    angle= np.random.uniform(0, 2 * np.pi)
    stack_x = np.random.uniform(0, 51) * np.cos(angle)
    stack_y = np.random.uniform(0, 51) * np.sin(angle)
    stack_height = np.random.randint(1,11) 
    emission_rate = np.random.randint(1,21)                                                                 
    windspeed= np.random.randint(1,26)

    #INPUT FOR ALGORITHM

    #Input For Coordinates Of Initial Measurement Point Of Robot (Can Start Out With Multiple)
    Measurement_Points=[[0,0]] 

    #Wind Angle In Degrees
    wind_direction= np.pi/2 #np.random.uniform(0, 2*np.pi)

    #GRID & STEPSIZE CALCULATION 

    x_range=np.arange(-100,100,0.5)
    y_range=np.arange(-100,100,0.5)
    z_range=np.arange(0,50,0.5)
    #It Turns Each Into 3D-Arrays For Calculations In All 3 Axis
    X,Y,Z = np.meshgrid(x_range,y_range,z_range)

    #WIND RELATED CALCULATIONS

    #For Future (height_slice_direction= math.atan2(measured_y - stack_y, measured_x - stack_x) *(180/(math.pi)))
    #3D-Arrays For Distance Between Coordinate & Origin
    x_origin=X-stack_x
    y_origin=Y-stack_y

    wind_xcomponent=windspeed*math.sin(wind_direction-np.pi)
    wind_ycomponent=windspeed*math.cos(wind_direction-np.pi)

    dot_product=wind_xcomponent*x_origin+wind_ycomponent*y_origin
    magnitudes= windspeed*((x_origin**2)+(y_origin**2))**0.5
    subtended=np.arccos(dot_product/(magnitudes))
    hypotenuse=((x_origin**2)+(y_origin**2))**0.5
    #Downwind Is X Of Gaussian Plume Formula 
    downwind=np.cos(subtended)*hypotenuse
    #Turns It Into Complex Format To Make Correct Calculations Later On.
    downwindimag=np.vectorize(complex)(downwind.real, downwind.imag)

    #PASQUILL CONSTANTS & SIGMA CALCULATION

    #Pasquill Constants Found In Resources
    stability_class= {1 : (122.8,0.94470,24.1670,2.5334),
                    2: (90.673, 0.93198, 18.3330, 1.8096),
                    3: (61.141, 0.91465, 12.5, 1.0857),
                    4: (34.459, 0.86974, 8.3330, 0.72382),
                    5: (24.26, 0.83660, 6.25, 0.54287),
                    6: (15.209, 0.81558, 4.1667, 0.36191)
                    }
    P_a=stability_class[stability][0] 
    P_b=stability_class[stability][1]
    P_c=stability_class[stability][2]
    P_d=stability_class[stability][3]
    sig_z = P_a*(abs((downwind/1000))**P_b)
    sig_z[sig_z > 5000] = 5000
    theta=0.017453293*(P_c-P_d*np.log((downwindimag/1000)))
    sig_y=(465.11628*downwind/1000)*np.tan(theta)

    #CROSSWIND AND GAUSSIAN FORMULA

    crosswind=np.sin(subtended)*hypotenuse
    indix,indiy,indiz=np.where((downwind > 0))
    Concentration= np.zeros(np.shape(downwind))
    Concentration[indix,indiy,indiz] =((emission_rate/(2*math.pi*windspeed*sig_y[indix,indiy,indiz]*sig_z[indix,indiy,indiz]))*(math.e**(-crosswind[indix,indiy,indiz]**2/(2*sig_y[indix,indiy,indiz]**2))* (math.e**(-(Z[indix,indiy,indiz]-stack_height)**2/(2*sig_z[indix,indiy,indiz]**2))+ math.e**(-(Z[indix,indiy,indiz]+stack_height)**2/(2*sig_z[indix,indiy,indiz]**2)))))
    #Prevents Issues Regarding Algorithm And Measurements
    Concentration[Concentration < 3e-4] = 0
    Concentration2D=(Concentration[:,:, 0])

    exportdf['Angle Of Plume'].append(angle)
    exportdf['Plume Origin'].append(f'{stack_x,stack_y,stack_height}')
    exportdf['Emission Rate'].append(emission_rate)
    #SETUP FOR WHILE-LOOP
    
    exportdf['Error'].append(0)
    exportdf['Failure'].append(0)
    exportdf['Unusable'].append(0)
    i=0
    Total_Distance=0
    #This Is For When negative_condition Or positive_condition Have Invalid Values
    Error_String='.'
    #Redundancy Is A Method To Break The While-Loop From Inside A For-Loop As An Extra Measure
    Redundancy=1
    Concentration_At_Measurement = []
    #This Is Distance Between Gas Leak And Robot
    RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y])
    if RadiusCheck <= 2:
        exportdf['Unusable'][-1]+=1
        exportdf['Distance'].append(Total_Distance)
        exportdf['# Of Measurements'].append(len(Measurement_Points))
        exportdf['Time Taken'].append(len(Measurement_Points)*10+(Total_Distance/2000)*60)
        continue
    #Adds Concentration To Concentration_At_Measurement For Each Coordinate.]
    for XYCoordinate in Measurement_Points:
        Concentration_At_Measurement.append(Concentration2D[((round(100+XYCoordinate[1]))*2),((round(100+XYCoordinate[0]))*2)])

    if Concentration_At_Measurement[-1]>0:
        exportdf['Unusable'][-1]+=1
        exportdf['Distance'].append(Total_Distance)
        exportdf['# Of Measurements'].append(len(Measurement_Points))
        exportdf['Time Taken'].append(len(Measurement_Points)*10+(Total_Distance/2000)*60)
        continue

    #This Is Where Non-Iterative Calculations For The Algorithm Belong 

    octax = np.array([0,31.82,0, -31.82, -45, -31.82, 0, 31.82,45,75, 53.03, 0, -53.03, -75, -53.03, 0, 53.03 ])
    octay = np.array([0,31.82, 45, 31.82, 0, -31.82, -45, -31.82,0,0, 53.03, 75, 53.03, 0, -53.03, -75, -53.03])
    x= np.diff(octax)
    y= np.diff(octay)

    while RadiusCheck*Redundancy >= 2:   
        RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y]) 
        #This Is Where Iterative Calculations For The Algorithm Belong.

        if Concentration_At_Measurement[-1] <= 0:
            #In Case Of Overshooting
            if sum(Concentration_At_Measurement) > 0:
                negative_condition=[0,0]
                positive_condition=[0,0]
                move_direction=[-Measurement_Points[-1][0]-Measurement_Points[-(2%len(Concentration_At_Measurement))][0],-Measurement_Points[-1][1]-Measurement_Points[-(2%len(Concentration_At_Measurement))][1],positive_condition,negative_condition]
            else:
                #In Case Of Concentration Never Having Been Higher Than 0
                if i >= len(octax)-1:
                    exportdf['Failure'][-1]+=1
                    break
                negative_condition=[0,0]
                positive_condition=[0,0]
                move_direction=[x[i],y[i],positive_condition,negative_condition]
        #In Case Of Concentration Being Lower Than The Concentration At The Previous Measurement Point
        elif Concentration_At_Measurement[-1] < Concentration_At_Measurement[-(2%len(Concentration_At_Measurement))]:
            negative_condition=[0,0]
            positive_condition=[0,0]
            move_direction=[-Measurement_Points[-1][0]-Measurement_Points[-(2%len(Concentration_At_Measurement))][0],-Measurement_Points[-1][1]-Measurement_Points[-(2%len(Concentration_At_Measurement))][1],positive_condition,negative_condition]
        #In Case Of Concentration Being Higher Than The Concentration At The Previous Measurement Point    
        elif Concentration_At_Measurement[-1] >= Concentration_At_Measurement[-(2%len(Concentration_At_Measurement))]:
            negative_condition=[0,0]
            positive_condition=[0,0]
            break
            
            
        #REDUNDANCY CHECKS & NEXT POSITION CALCULATION

        Measurement_Points.append([])
        #For-loop Takes Into Account X & Y.
        for j in range(2):
            #If Suggested move_direction Is Infeasible, It Uses positive_condition
            if Measurement_Points[-2][j]+move_direction[j]>100:
                Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[2][j])
            #If Suggested move_direction Is Infeasible, It Uses negative_condition
            elif Measurement_Points[-2][j]+move_direction[j]<-100:
                Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[3][j])
            #Uses Normal move_direction Values     
            else:
                Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[j])  
            #Checks If Newest Measurement Point Is Invalid
            if abs(Measurement_Points[-1][j]) > 100:
                Redundancy=0
                Measurement_Points.pop(-1)
                Error_String=', Simulation Ended Due To An Error In Negative Or Positive Condition Calculation'
                exportdf['Error'][-1]+=1
                break
            #Rounds To The Nearest Half. (Due To Resolution Of Gas Plume Being In 0.5 By 0.5 Meters)
            Measurement_Points[-1][j]=(round(Measurement_Points[-1][j]*2)/2)
        if len(Measurement_Points)>=2:
            Total_Distance+=math.dist(Measurement_Points[i],Measurement_Points[i-1])
        if len(Measurement_Points)*10+(Total_Distance/2000)*60 > 300:
                exportdf['Failure'][-1]+=1
                break

        #UPDATES TO WHILE-LOOP

        i+=1
        Concentration_At_Measurement.append(Concentration2D[((round(100+Measurement_Points[i][1]))*2)][((round(100+Measurement_Points[i][0]))*2)])  
    print(f'Robot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {len(Measurement_Points)} Number Of Measurement Points, Until It Got To {[round(Measurement_Points[-1][0],ndigits=2),round(Measurement_Points[-1][1],ndigits=2)]} With A Concentration Of {Concentration_At_Measurement[-1]}. The Objective Coordinates Were {[stack_x, stack_y]}{Error_String}')
    exportdf['Distance'].append(Total_Distance)
    exportdf['# Of Measurements'].append(len(Measurement_Points))
    exportdf['Time Taken'].append(len(Measurement_Points)*10+(Total_Distance/2000)*60)
print(exportdf['Distance'],exportdf['# Of Measurements'],exportdf['Time Taken'],exportdf['Failure'],exportdf['Error'],exportdf['Unusable'])

pd.DataFrame(exportdf).to_excel('Test.xlsx')
