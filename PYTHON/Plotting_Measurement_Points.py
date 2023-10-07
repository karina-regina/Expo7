#This Python Script Plots The Route Of The Robot And Objective Coordinates.

#Right Now, The Route It Takes Is A Spiral. To Make Changes To The Algorithm, Go To Line 
#The Code Runs Until It Gets Within The 2 Meter Range.
#If It Does Get Within The Range, It Stops. (In Future, Change So That It Stops Once It Finalizes Its Guess)
#By Closing Plot, It Prints Total Distance Travelled, Number Of Measurements, Last Measurement Coordinates And Objective Coordinates.

#Input Of Initial Measurement_Points And Objective Coordinates Is At Line #21 And #24. 
#To Add An Algorithm, It Has Been Split Into 2 Locations. (Non-Iterative Calculations And Iterative Calculations)

#Cons Of Code
#It Is A Bit Long, There Are Ways To Shorten It.
#It Gets Laggy Once In the Hundreds Of Measurement Points.
#There Is No Algorithm Implemented Yet. (Just Random)
#Work Needs To Be Done To Make Sure Repetitious Measurement Points Aren't Taken
#positive_condition And negative_condition Would Be Unnecessary In Some Algorithms, But Necessary For Those That Don't Calculate Whether Their Movements Are Feasible (It's A Layer Of Redundancy)
#Calculations For positive_condition And negative_condition Have To Be A List With 2 Values, So They Have To Be Calculated Or Predetermined Regardless Of Whether Intended X And Y Calculations Are Always Feasible. 
#(Add More If You See Them) 

#Input For Coordinates Of Objective That Robot Is Heading Towards (Only Works With One)
Desired_Coordinates=[50,50]

#Input For Coordinates Of Initial Measurement Point Of Robot (Can Start Out With Multiple)
Measurement_Points=[[0,0]] 

#Not All Of These Libraries Are Necessary, 
#I Just Copy-Pasted Them From Some Of Our Other Code

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

#Setup For Interactive Plot
plt.ion()
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
ax.scatter(Desired_Coordinates[0], Desired_Coordinates[1], s=30, marker='D',label='Objective')
if len(Measurement_Points)>1:
     for l in range(len(Measurement_Points)):
          ax.scatter(Measurement_Points[l][0], Measurement_Points[l][1], s=10, color='b')
          ax.plot([Measurement_Points[l-1][0],Measurement_Points[l][0]], [Measurement_Points[l-1][1], Measurement_Points[l][1]],color='black')
annotation = plt.annotate('', xy=(0, 0),fontsize=8, xytext=(-130, 114))
plt.title('Simulation Of The Route Of Robot')
def update_annotation(new_text):
    annotation.set_text(new_text)

#Setup For While-Loop
i=0
Total_Distance=0
Error_String='.'
Redundancy=1

#This Is Where Non-Iterative Calculations For The Algorithm Belong 
#In This Case There Are None, Since It Is Random.

while math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), Desired_Coordinates)*Redundancy >= 2:   
    
    #This Is Where Iterative Calculations For The Algorithm Belong. 
    #move_direction Contains Movement In The X And Y Direction, 
    #But Also positive_condition And negative_condition. These Are Lists With 2 Values (X And Y).
    #They Act As A Contigency Plan In Case The First Values Lead To Outside The Feasible Range.
    #Make Sure positive_condition And negative_condition exist (Can Be Calculated In The Non-Iterative Section As Well)
    #Make Sure move_direction Is A List With 4 Elements, With The Last 2 Being Lists Of Length 2 And The First Two Being Integers Or Floats. Otherwise Code Breaks.
    negative_condition=[np.random.randint(0,20),np.random.randint(0,20)]
    positive_condition=[np.random.randint(-20,0),np.random.randint(-20,0)]
    move_direction=[np.random.randint(-10,10),np.random.randint(-10,10),positive_condition,negative_condition]
    
    #For-Loop And If Statements Make New Measurement Points, 
    #Makes Sure Robot Doesn't Pass Plot Boundaries And 
    #That Robot Doesn't Keep Getting Stuck In The Boundaries (X And Y)
    Measurement_Points.append([])
    for j in range(2):
        if Measurement_Points[-2][j]+move_direction[j]>100:
             Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[2][j])
        elif Measurement_Points[-2][j]+move_direction[j]<-100:
             Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[3][j])        
        else:
             Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[j])  
        if abs(Measurement_Points[-1][j]) > 100:
             Redundancy=0
             Measurement_Points.pop(-1)
             Error_String=', Simulation Ended Due To An Error In Negative Or Positive Condition Calculation'
             break
        
    #If Statement To Determine Whether This Calculation Of Total Distance Will Lead To Error
    if len(Measurement_Points)>=2:
        Total_Distance+=math.dist(Measurement_Points[-1],Measurement_Points[-2])

    #Updates To Plot And While-Loop
    update_annotation(f'Number Of Measurements Taken {i+1} \nCurrent Distance Between Last Measurement Point And Objective Is {round(math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), Desired_Coordinates),ndigits =2)} Meters\nTotal Distance Robot Has Travelled So Far {round(Total_Distance,ndigits=2)} Meters')
    ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, color='b')
    ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]],color='black')
    i+=1
    plt.draw()
    plt.pause(0.4)

print(f'Robot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {i} Number Of Measurement Points, Until It Got To {Measurement_Points[-1]}. The Objective Coordinates Were {Desired_Coordinates}{Error_String}')
plt.ioff()
plt.show()
