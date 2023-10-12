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

#INPUT FOR GAS PLUME

stability = 1#set from 1-6                                                                                      

stack_x = 20#np.random.randint(-50,51)
stack_y = 5#np.random.randint(-50,51)
stack_height = 3#np.random.randint(1,11)      
emission_rate = 1#np.random.randint(1,21)                                                                 
windspeed= 4#np.random.randint(1,26)
wind_direction= -135#np.random.randint(0,360)

#INPUT FOR ALGORITHM

#Input For Coordinates Of Initial Measurement Point Of Robot (Can Start Out With Multiple)
Measurement_Points=[[0,0]] 

#GRID & STEPSIZE CALCULATION 

stability_str = ['Very unstable', 'Moderately unstable', 'Slightly unstable', 'Neutral', 'Moderately stable', 'Very stable'] # Possibly Useless
ssv= stability_str[stability-1] # Possibly Useless
x_range=np.arange(-100,100,0.5)
y_range=np.arange(-100,100,0.5)
z_range=np.arange(0,50,0.5)
X,Y,Z = np.meshgrid(x_range,y_range,z_range)

#WIND RELATED CALCULATIONS

#height_slice_direction= math.atan2(measured_y - stack_y, measured_x - stack_x) *(180/(math.pi))
x_origin=X-stack_x
y_origin=Y-stack_y

wind_xcomponent=windspeed*math.sin(math.radians(wind_direction-180))
wind_ycomponent=windspeed*math.cos(math.radians(wind_direction-180))

dot_product=wind_xcomponent*x_origin+wind_ycomponent*y_origin
magnitudes= windspeed*((x_origin**2)+(y_origin**2))**0.5
subtended=np.arccos(dot_product/(magnitudes+1e-5))
hypotenuse=((x_origin**2)+(y_origin**2))**0.5
downwind=np.cos(subtended)*hypotenuse
downwindimag=np.vectorize(complex)(downwind.real, downwind.imag)

#PASQUILL CONSTANTS & SIGMA CALCULATION

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
Concentration[indix,indiy,indiz] =1e6*((emission_rate/(2*math.pi*windspeed*sig_y[indix,indiy,indiz]*sig_z[indix,indiy,indiz]))*(math.e**(-crosswind[indix,indiy,indiz]**2/(2*sig_y[indix,indiy,indiz]**2))* (math.e**(-(Z[indix,indiy,indiz]-stack_height)**2/(2*sig_z[indix,indiy,indiz]**2))+ math.e**(-(Z[indix,indiy,indiz]+stack_height)**2/(2*sig_z[indix,indiy,indiz]**2)))))
Concentration[Concentration < 1e-2] = 0
Concentration2D=(Concentration[:,:, 0])

#Setup For Interactive Plot

plt.ion()
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
plt.pcolor(x_range, y_range, Concentration2D[:,:], shading=None, cmap='jet')
if len(Measurement_Points)>1:
     for l in range(len(Measurement_Points)):
          ax.scatter(Measurement_Points[l][0], Measurement_Points[l][1], s=5,alpha=0.5, color='white')
          ax.plot([Measurement_Points[l-1][0],Measurement_Points[l][0]], [Measurement_Points[l-1][1], Measurement_Points[l][1]],s=4,alpha=0.5,color='white')
ax.scatter(stack_x,stack_y,s=8,color='white',marker='x',label='Gas Leak Source')
annotation = plt.annotate('', xy=(0, 0),fontsize=8, xytext=(-130, 114))
plt.title('Simulation Of The Route Of Robot')
def update_annotation(new_text):
    annotation.set_text(new_text)
plt.legend()
plt.colorbar()

#Setup For While-Loop

i=0
Total_Distance=0
Error_String='.'
Redundancy=1
Concentration_At_Measurement = []
RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y])
for XYCoordinate in Measurement_Points:
     Concentration_At_Measurement.append(Concentration2D[((round(100+XYCoordinate[0]))*2),((round(100+XYCoordinate[1]))*2)])

#This Is Where Non-Iterative Calculations For The Algorithm Belong 

Spiral = np.linspace(0, 10 * np.pi, 25)
radius = np.linspace(0, 100, 25)  
x_coordinates = radius*np.cos(Spiral)
y_coordinates = radius*np.sin(Spiral)
x= np.diff(x_coordinates)
y= np.diff(y_coordinates)

while RadiusCheck*Redundancy >= 2:   
    
    #This Is Where Iterative Calculations For The Algorithm Belong.

    if Concentration_At_Measurement[i] <= 0:
        #In Case Of Overshooting
        if sum(Concentration_At_Measurement) > 0:
            negative_condition=[0,0]
            positive_condition=[0,0]
            print('wo')
            move_direction=[-Measurement_Points[i][0]-Measurement_Points[i-1][0],-Measurement_Points[i][1]-Measurement_Points[i-1][1],positive_condition,negative_condition]
        else:
            #In Case Of Concentration Never Having Been Higher Than 0
            if i >= len(Spiral)-1:
                break
            negative_condition=[0,0]
            positive_condition=[0,0]
            move_direction=[x[i],y[i],positive_condition,negative_condition]
    #In Case Of Concentration Being Lower Than The Concentration At The Previous Measurement Point
    elif Concentration_At_Measurement[i] < Concentration_At_Measurement[i-1]:
        negative_condition=[0,0]
        positive_condition=[0,0]
        move_direction=[-Measurement_Points[i][0]-Measurement_Points[i-1][0],-Measurement_Points[i][1]-Measurement_Points[i-1][1],positive_condition,negative_condition]
    #In Case Of Concentration Being Higher Than The Concentration At The Previous Measurement Point    
    elif Concentration_At_Measurement[i] >= Concentration_At_Measurement[i-1]:
        negative_condition=[0,0]
        positive_condition=[0,0]
        move_direction=[0,0,positive_condition,negative_condition]
        
    #Redundancy Checks & Next Position Calculation

    RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y]) 
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
        Measurement_Points[-1][j]=(round(Measurement_Points[-1][j]*2)/2)
    if len(Measurement_Points)>=2:
        Total_Distance+=math.dist(Measurement_Points[i],Measurement_Points[i-1])

    #Updates To Plot And While-Loop
    update_annotation(f'Current Distance Between Last Measurement Point And Objective Is {round(math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), [stack_x, stack_y]),ndigits =2)} Meters\nRobot Has Travelled {round(Total_Distance,ndigits=2)} Meters & Taken {i+1} Amount Of Measurements\nCurrent Concentration At Coordinate Is {round(Concentration_At_Measurement[-1],ndigits=2)}')
    ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, color='white',alpha=0.5)
    ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]],alpha=0.5,color='white')
    i+=1
    print(Measurement_Points[-1])
    Concentration_At_Measurement.append(Concentration2D[((round(100+Measurement_Points[i][0]))*2)][((round(100+Measurement_Points[i][1]))*2)])  
    plt.draw()
    plt.pause(0.4)
print(f'Robot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {i} Number Of Measurement Points, Until It Got To {[round(Measurement_Points[-1][0],ndigits=2),round(Measurement_Points[-1][1],ndigits=2)]} With A Concentration Of {Concentration_At_Measurement[-1]}. The Objective Coordinates Were {[stack_x, stack_y]}{Error_String}')
plt.ioff()
plt.show()
