#This Python Script Plots The Route Of The Robot And Objective Coordinates.

#Right Now, The Route It Takes Depends On A Random Integer From -10 To 10. (Line #56, #59, #61) This Part Is Modifiable
#The Code Runs Until It Gets Within The 2 Meter Range.(Checked At Line #49) (Unlikely, Given How It Is Random)
#If It Does Get Within The Range, It Stops. (In Future, Change So That It Stops Once It Finalizes Its Guess)
#By Closing Plot, It Prints Total Distance Travelled, Number Of Measurements, Last Measurement Coordinates And Objective Coordinates.

#Input Of Code Is At Line #21 And #24. 
#To Add An Algorithm, Change Line #56, #59 And #61
#To Change Minimum Distance From Objective At Which Code Stops, Change Line #49

#Cons Of Code
#It Is A Bit Long, There Are Ways To Shorten It.
#Legend Stops At 23 Measurement Points. (More Won't Fit, A Fix Is To Call It M#i, So More Can Fit.)
#It Gets Laggy Once In the Hundreds Of Measurement Points.
#There Is No Algorithm Implemented Yet.
#Work Needs To Be Done To Make Sure Repetitious Measurement Points Aren't Taken
#(Add More If You See Them) 


#Input For Coordinates Of Objective That Robot Is Heading Towards
Desired_Coordinates=[50,50]

#Input For Coordinates Of Initial Measurement Point Of Robot
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

plt.ion()
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
ax.scatter(Desired_Coordinates[0], Desired_Coordinates[1], s=30, marker='D',label='Objective')
i=0
Total_Distance=0

while math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), Desired_Coordinates) >= 2:   
    #For-Loop And If Statements Make New Measurement Points, 
    #Makes Sure Robot Doesn't Pass Plot Boundaries And 
    #That Robot Doesn't Keep Getting Stuck In The Boundaries (X And Y)
    for j in range(2):
        if Measurement_Points[-1][j]>100:
             Measurement_Points[-1][j]=100
             Measurement_Points.append([Measurement_Points[-1][0]+np.random.randint(-20,0),Measurement_Points[-1][1]+np.random.randint(-20,0)])
        elif Measurement_Points[-1][j]<-100:
             Measurement_Points[-1][j]=-100
             Measurement_Points.append([Measurement_Points[-1][0]+np.random.randint(0,20),Measurement_Points[-1][1]+np.random.randint(0,20)])         
    else:
         Measurement_Points.append([Measurement_Points[-1][0]+np.random.randint(-10,10),Measurement_Points[-1][1]+np.random.randint(-10,10)])
    if len(Measurement_Points)>=2:
        Total_Distance+=math.dist(Measurement_Points[-1],Measurement_Points[-2])
     #This If-Statement Is To Not Add To Legend After 23 Measurement Points. (Lack Of Space)
    if i<23:
        ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, label=(f'Measurement #{i+1}'))
    else:
         ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10)
    ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]])
    i+=1
    plt.title(f'Route Of Robot ({i+1} Number Of Measurement Points Taken So Far)(Distance From Objective Is {math.dist(Measurement_Points[-1],Desired_Coordinates)})')
    ax.legend(fontsize='xx-small',bbox_to_anchor=(1.1, 0),ncol=6)
    plt.draw()
    plt.pause(0.4)
print(f'Robot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {i} Number Of Measurement Points, Until It Got To {Measurement_Points[-1]}. The Objective Coordinates Were {Desired_Coordinates}')
plt.ioff()
plt.show()
