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
from sympy import symbols, solve, Eq
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.optimize import curve_fit
np.random.seed(513)

#INPUT FOR GAS PLUME
#set from 1-6
stability = np.random.randint(1,7)                                                                                      
#Origin Coordinates Of Plume
angle= np.random.uniform(0, 2 * np.pi)
stack_x = np.random.uniform(0, 51) * np.cos(angle)
stack_y = np.random.uniform(0, 51) * np.sin(angle)
stack_height = np.random.randint(1,11) 
emission_rate = np.random.randint(1,21)                                                                 
windspeed= 1 #np.random.randint(1,16)
#Wind Angle In Degrees
wind_direction= np.random.uniform(0, 2*np.pi - 1e-5)

#INPUT FOR ALGORITHM

#Input For Coordinates Of Initial Measurement Point Of Robot (Can Start Out With Multiple)
Measurement_Points=[[0,0]] 

#GRID & STEPSIZE CALCULATION 

x_range=np.arange(-100,101,0.5)
y_range=np.arange(-100,101,0.5)
z_range=np.arange(0,51,0.5)
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
Concentration[Concentration < 3e-5] = 0
Concentration2D=(Concentration[:,:, 0])

#SETUP FOR INTERACTIVE PLOT

#plt.ion Allows For The Continuously Updated Plot 
plt.ion()
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
#Heatmap
#plt.contour(x_range, y_range, Concentration2D,cmap='summer',s=3)
plt.pcolor(x_range, y_range, Concentration2D, shading=None, cmap='jet')
#Checks For Length Of Measurement Points
if len(Measurement_Points)>1:
     for l in range(1,len(Measurement_Points)):
          ax.scatter(Measurement_Points[l][0], Measurement_Points[l][1], s=10,alpha=0.5, color='white')
          ax.plot([Measurement_Points[l-1][0],Measurement_Points[l][0]], [Measurement_Points[l-1][1], Measurement_Points[l][1]],alpha=0.5,color='white')
#Plots Origin Point Of Plume
ax.scatter(stack_x,stack_y,s=8,color='white',marker='x',label='Gas Leak Source')
#Setup For Info At Top-Left Of Plot
annotation = plt.annotate('', xy=(0, 0),fontsize=8, xytext=(-130, 114))
plt.title('Simulation Of The Route Of Robot')
#Setup For Updating Info At Top-Left Of Plot
def update_annotation(new_text):
    annotation.set_text(new_text)
plt.legend()
plt.colorbar()

#SETUP FOR WHILE-LOOP
i=0
Total_Distance=0
#This Is For When negative_condition Or positive_condition Have Invalid Values
Error_String='.'
#Redundancy Is A Method To Break The While-Loop From Inside A For-Loop As An Extra Measure
Redundancy=1
Concentration_At_Measurement = []

#Adds Concentration To Concentration_At_Measurement For Each Coordinate.]
for XYCoordinate in Measurement_Points:
     Concentration_At_Measurement.append(Concentration2D[((round(100+XYCoordinate[1]))*2),((round(100+XYCoordinate[0]))*2)])

#This Is Where Non-Iterative Calculations For The Algorithm Belong 
step1ss,step2ss=5,5
zeroCpoint=[]
phase1=True
phase3=False
phase4=False
phase5=False
begin_of_crosswind_index=0
correction=20
Meas_Count=0
wind_dir=(3*np.pi/2- wind_direction)-np.pi/2
octax = np.array([0,55.0, 47.631, 27.5, 0.0, -27.5, -47.631, -55.0, -47.631, -27.5, -0.0, 27.5, 47.631, 33.25, 20.731, -7.399, -29.957, -29.957, -7.399, 20.731, 77.75, 72.5, 57.458, 34.656, 7.174, -21.277, -46.855, -66.104, -76.426, -76.426, -66.104, -46.855, -21.277, 7.174, 34.656, 57.458, 72.5, 22.125, 6.837, -17.9, -17.9, 6.837, 88.875, 84.06, 70.135, 48.61, 21.818, -7.339, -35.701, -60.193, -78.163, -87.663, -87.663, -78.163, -60.193, -35.701, -7.339, 21.818, 48.61, 70.135, 84.06, 44.375, 35.9, 13.713, -13.713, -35.9, -44.375, -35.9, -13.713, 13.713, 35.9, 66.625, 60.027, 41.54, 14.825, -14.825, -41.54, -60.027, -66.625, -60.027, -41.54, -14.825, 14.825, 41.54, 60.027, 16.562, 0.0, -16.562, -0.0, 94.438, 89.815, 76.402, 55.509, 29.183, 0.0, -29.183, -55.509, -76.402, -89.815, -94.438, -89.815, -76.402, -55.509, -29.183, -0.0, 29.183, 55.509, 76.402, 89.815, 27.688, 13.844, -13.844, -27.688, -13.844, 13.844, 83.312, 78.288, 63.821, 41.656, 14.467, -14.467, -41.656, -63.821, -78.288, -83.312, -78.288, -63.821, -41.656, -14.467, 14.467, 41.656, 63.821, 78.288, 38.812, 27.445, 0.0, -27.445, -38.812, -27.445, -0.0, 27.445, 72.188, 66.693, 51.044, 27.625, 0.0, -27.625, -51.044, -66.693, -72.188, -66.693, -51.044, -27.625, -0.0, 27.625, 51.044, 66.693, 49.938, 42.01, 20.745, -7.107, -32.702, -47.915, -47.915, -32.702, -7.107, 20.745, 42.01, 61.062, 54.068, 34.687, 7.36, -21.653, -45.706, -59.288, -59.288, -45.706, -21.653, 7.36, 34.687, 54.068])
octay = np.array([0,0.0, 27.5, 47.631, 55.0, 47.631, 27.5, 0.0, -27.5, -47.631, -55.0, -47.631, -27.5, 0.0, 25.996, 32.416, 14.427, -14.427, -32.416, -25.996, 0.0, 28.087, 52.38, 69.599, 77.418, 74.782, 62.046, 40.93, 14.287, -14.287, -40.93, -62.046, -74.782, -77.418, -69.599, -52.38, -28.087, 0.0, 21.042, 13.005, -13.005, -21.042, 0.0, 28.858, 54.588, 74.403, 86.155, 88.571, 81.389, 65.387, 42.3, 14.628, -14.628, -42.3, -65.387, -81.389, -88.571, -86.155, -74.403, -54.588, -28.858, 0.0, 26.083, 42.203, 42.203, 26.083, 0.0, -26.083, -42.203, -42.203, -26.083, 0.0, 28.908, 52.09, 64.955, 64.955, 52.09, 28.908, -0.0, -28.908, -52.09, -64.955, -64.955, -52.09, -28.908, 0.0, 16.562, 0.0, -16.562, 0.0, 29.183, 55.509, 76.402, 89.815, 94.438, 89.815, 76.402, 55.509, 29.183, 0.0, -29.183, -55.509, -76.402, -89.815, -94.438, -89.815, -76.402, -55.509, -29.183, 0.0, 23.978, 23.978, 0.0, -23.978, -23.978, 0.0, 28.495, 53.552, 72.151, 82.047, 82.047, 72.151, 53.552, 28.495, 0.0, -28.495, -53.552, -72.151, -82.047, -82.047, -72.151, -53.552, -28.495, 0.0, 27.445, 38.812, 27.445, 0.0, -27.445, -38.812, -27.445, 0.0, 27.625, 51.044, 66.693, 72.188, 66.693, 51.044, 27.625, 0.0, -27.625, -51.044, -66.693, -72.188, -66.693, -51.044, -27.625, 0.0, 26.998, 45.425, 49.429, 37.74, 14.069, -14.069, -37.74, -49.429, -45.425, -26.998, 0.0, 28.377, 50.253, 60.617, 57.094, 40.492, 14.613, -14.613, -40.492, -57.094, -60.617, -50.253, -28.377])
x= np.diff(octax)
y= np.diff(octay)


while Redundancy == 1:   
    #This Is Where Iterative Calculations For The Algorithm Belong.

    if phase5 == True:
        move_direction=[-(step4ss+correction)*np.sin(wind_dir),(step4ss+correction)*np.cos(wind_dir),positive_condition,negative_condition]           
        Meas_Count+=1
        if Meas_Count>=4:
            break
    elif Concentration_At_Measurement[-1] <= 0:
        #In Case Of Overshooting
        if sum(Concentration_At_Measurement) > 0:
            negative_condition=[0,0]
            positive_condition=[0,0]
            zeroCpoint.append(Measurement_Points[-1])
            if phase1 == True:
                step1ss=-step1ss
                phase1=False
                move_direction=[Measurement_Points[np.where(np.array(Concentration_At_Measurement)>0)[0][0]][0]-Measurement_Points[-1][0]+step1ss*np.cos(wind_dir),Measurement_Points[np.where(np.array(Concentration_At_Measurement)>0)[0][0]][1]-Measurement_Points[-1][1]+step1ss*np.sin(wind_dir),positive_condition,negative_condition]
            elif phase4 == True:
                step4ss=4
                if correction == 0:
                    break
                move_direction=[-(step4ss+correction)*np.sin(wind_dir),(step4ss+correction)*np.cos(wind_dir),positive_condition,negative_condition]           
                correction=0
            else:
                phase3=True
                move_direction=[(zeroCpoint[0][0]-zeroCpoint[-1][0])/2 -10*np.sin(wind_dir),(zeroCpoint[0][1]-zeroCpoint[-1][1])/2 + 10*np.cos(wind_dir),positive_condition,negative_condition]
                end_of_crosswind_index=len(Measurement_Points)
        else:   
            #In Case Of Concentration Never Having Been Higher Than 0
            if i >= len(octax)-1:
                break
            negative_condition=[0,0]
            positive_condition=[0,0]
            move_direction=[x[i],y[i],positive_condition,negative_condition]
            begin_of_crosswind_index=len(Measurement_Points)
    #In Case Of Concentration Being Lower Than The Concentration At The Previous Measurement Point
    elif Concentration_At_Measurement[-1] < Concentration_At_Measurement[-(2%len(Concentration_At_Measurement))]:
        negative_condition=[0,0]
        positive_condition=[0,0]
        move_direction=[step1ss*np.cos(wind_dir),step1ss*np.sin(wind_dir),positive_condition,negative_condition]
        if phase4 == True:
            if correction == 0:
                phase5=True
            else:
                step4ss=4
            move_direction=[-(step4ss+correction)*np.sin(wind_dir),(step4ss+correction)*np.cos(wind_dir),positive_condition,negative_condition]           
            correction=0
    #In Case Of Concentration Being Higher Than The Concentration At The Previous Measurement Point    
    elif Concentration_At_Measurement[-1] >= Concentration_At_Measurement[-(2%len(Concentration_At_Measurement))]:
        negative_condition=[0,0]
        positive_condition=[0,0]
        move_direction=[(7*step1ss/4)*np.cos(wind_dir),(7*step1ss/4)*np.sin(wind_dir),positive_condition,negative_condition]
        if phase3 == True:
            move_direction=[20*np.sin(wind_dir),-20*np.cos(wind_dir),positive_condition,negative_condition]
            phase3=False
            phase4=True
            centerlineindex=len(Measurement_Points)
        elif phase4 == True:
            if correction == 20:
                step4ss=-4
                correction=0
            move_direction=[-(step4ss+correction)*np.sin(wind_dir),(step4ss+correction)*np.cos(wind_dir),positive_condition,negative_condition]

    #REDUNDANCY CHECKS & NEXT POSITION CALCULATION

    Measurement_Points.append([])
    #For-loop Takes Into Account X & Y.
    for j in range(2):
        #If Suggested move_direction Is Infeasible, It Uses positive_condition
        #For This Algorithm's Purposes, This Isn't Wanted
        if Measurement_Points[-2][j]+move_direction[j]>100:
            Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[2][j])
            if phase5 == True:
                Redundancy=0
        #If Suggested move_direction Is Infeasible, It Uses negative_condition
        #For This Algorithm's Purposes, This Isn't Wanted
        elif Measurement_Points[-2][j]+move_direction[j]<-100:
            Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[3][j])
            if phase5 == True:
                Redundancy=0
        #Uses Normal move_direction Values     
        else:
            Measurement_Points[-1].append(Measurement_Points[-2][j]+move_direction[j])  
        #Checks If Newest Measurement Point Is Invalid
        if abs(Measurement_Points[-1][j]) > 100:
            Redundancy=0
            Measurement_Points.pop(-1)
            Error_String=', Simulation Ended Due To An Error In Negative Or Positive Condition Calculation'
            break
        #Rounds To The Nearest Half. (Due To Resolution Of Gas Plume Being In 0.5 By 0.5 Meters)
        Measurement_Points[-1][j]=(round(Measurement_Points[-1][j]*2)/2)
    
    #UPDATES TO PLOT & WHILE-LOOP

    update_annotation(f'Current Distance Between Last Measurement Point And Objective Is {round(math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), [stack_x, stack_y]),ndigits =2)} Meters\nRobot Has Travelled {round(Total_Distance,ndigits=2)} Meters & Taken {i+1} Amount Of Measurements\nCurrent Concentration At Coordinate Is {round(Concentration_At_Measurement[-1],ndigits=2)}')
    ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, color='white',alpha=0.5)
    ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]],alpha=0.5,color='white')
    i+=1
    Concentration_At_Measurement.append(Concentration2D[((round(100+Measurement_Points[i][1]))*2)-1][((round(100+Measurement_Points[i][0]))*2)-1])  
    if len(Measurement_Points)>=2:
        Total_Distance+=math.dist(Measurement_Points[i],Measurement_Points[i-1])
        if Measurement_Points[-1] == Measurement_Points[-2]:
            Measurement_Points.pop(-1)
            Concentration_At_Measurement.pop(-1)
            break
    plt.draw()
    plt.pause(0.4)
update_annotation(f'Current Distance Between Last Measurement Point And Objective Is {round(math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), [stack_x, stack_y]),ndigits =2)} Meters\nRobot Has Travelled {round(Total_Distance,ndigits=2)} Meters & Taken {i+1} Amount Of Measurements\nCurrent Concentration At Coordinate Is {round(Concentration_At_Measurement[-1],ndigits=2)}')
ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, color='white',alpha=0.5)
ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]],alpha=0.5,color='white')
plt.draw()

#CALCULATING X1,X2,X3,X4

Relative_0_Coordinates= Measurement_Points[-1]
if math.dist(Measurement_Points[centerlineindex+2],[50 * np.cos(wind_dir),50 * np.sin(wind_dir)])<math.dist(Measurement_Points[centerlineindex+1],[50 * np.cos(wind_dir),50 * np.sin(wind_dir)]):
    Rotation_Constant=-1
else:
    Relative_0_Coordinates= Measurement_Points[centerlineindex+1]
    Rotation_Constant=1
min_C_ind,max_C_ind=Concentration_At_Measurement.index(min(Concentration_At_Measurement[centerlineindex:])),Concentration_At_Measurement.index(max(Concentration_At_Measurement[centerlineindex:]))
#X1
Crosswind_list=Measurement_Points[begin_of_crosswind_index:end_of_crosswind_index]
Crosswind_distance=math.dist(max(Crosswind_list),min(Crosswind_list))
#X2
lowmax_distance=(math.dist(Measurement_Points[min_C_ind],Measurement_Points[max_C_ind]))
#Rotating Measurement Points According To If We Moved Towards Origin or Away From It.
Rotated_MM=(Measurement_Points[centerlineindex+1:])[::Rotation_Constant]
x_manh_dist,y_manh_dist = np.array(Rotated_MM)[:,0],np.array(Rotated_MM)[:,1]
CvsX_X_Axis = [sum(np.sqrt(np.diff(x_manh_dist)**2 + np.diff(y_manh_dist)**2)[:i]) for i in range(len(np.sqrt(np.diff(x_manh_dist)**2 + np.diff(y_manh_dist)**2)) + 1)]
CvsX_C_Axis = (Concentration_At_Measurement[centerlineindex+1:])[::Rotation_Constant]
max_C_ind=CvsX_C_Axis.index(max(CvsX_C_Axis))

x_cols = ['Crosswind', 'Distance On Centerline', 'Negative Gradient','Positive Gradient']
#Handling Exceptions For When Negative Or Positive Slope Aren't Available Due To Area Bounds
if max_C_ind == 0:
    neg_slope=(CvsX_C_Axis[max_C_ind+1]-CvsX_C_Axis[max_C_ind])/(CvsX_X_Axis[max_C_ind+1]-CvsX_X_Axis[max_C_ind]+1e-100)
    x_cols.pop(-1)
    X_vars=[[Crosswind_distance,lowmax_distance,neg_slope]]
elif max_C_ind == len(CvsX_C_Axis)-1:
    pos_slope=(CvsX_C_Axis[max_C_ind]-CvsX_C_Axis[max_C_ind-1])/(CvsX_X_Axis[max_C_ind]-CvsX_X_Axis[max_C_ind-1]+1e-100)
    x_cols.pop(-2)
    X_vars=[[Crosswind_distance,lowmax_distance,pos_slope]]
else:
    neg_slope=(CvsX_C_Axis[max_C_ind+1]-CvsX_C_Axis[max_C_ind])/(CvsX_X_Axis[max_C_ind+1]-CvsX_X_Axis[max_C_ind]+1e-100)
    pos_slope=(CvsX_C_Axis[max_C_ind]-CvsX_C_Axis[max_C_ind-1])/(CvsX_X_Axis[max_C_ind]-CvsX_X_Axis[max_C_ind-1]+1e-100)
    X_vars=[[Crosswind_distance,lowmax_distance,neg_slope,pos_slope]]

#DEFINING LDA MODEL & DEPLOYING IT

LDA_Data=pd.read_excel('6000_LDA_Data.xlsx')
model = LinearDiscriminantAnalysis(solver='lsqr').fit(LDA_Data[x_cols], LDA_Data['Stability Class'])
stability_class_prediction = model.predict(pd.DataFrame(X_vars, columns=x_cols))[0]

#DEFINING CENTERLINE GAUSSIAN PLUME FORMULA

stability_values= {1 : (122.8,0.94470,24.1670,2.5334),
                2: (90.673, 0.93198, 18.3330, 1.8096),
                3: (61.141, 0.91465, 12.5, 1.0857),
                4: (34.459, 0.86974, 8.3330, 0.72382),
                5: (24.26, 0.83660, 6.25, 0.54287),
                6: (15.209, 0.81558, 4.1667, 0.36191)
                }
P_a,P_b,P_c,P_d=(stability_values[stability_class_prediction][0]),(stability_values[stability_class_prediction][1]),(stability_values[stability_class_prediction][2]),(stability_values[stability_class_prediction][3])

def cGP(x_or, Qu, x_leak, h):
    y_or = 0
    y_leak = 0
    z_or = 0
    
    x = x_or - x_leak
    y = y_or - y_leak

    a = P_a
    b = P_b
    c = P_c
    d = P_d

    theta = np.radians(c - d * np.log(x / 1000))
    sy = 465.11628 * (x / 1000) * np.tan(theta)
    sz = np.minimum(a * (x / 1000) ** b, 5000)

    c = (Qu / (2 * np.pi * sy * sz)) * (
            np.exp(- (z_or - h) ** 2 / (2 * sz ** 2)) +
            np.exp(- (z_or + h) ** 2 / (2 * sz ** 2))
    ) * np.exp(-y ** 2 / (2 * sy ** 2))

    return c
#Fixing Bounds For Solution Space

beta0 = [3.0, -20, 0.5]
lowerbound=-abs(math.dist([50 * np.cos(wind_dir),50 * np.sin(wind_dir)],Relative_0_Coordinates))
bounds=([1,lowerbound,1],[20,0,11])

#Getting Origin Prediction And Plotting
beta_fit, covmatrix = curve_fit(cGP, CvsX_X_Axis, CvsX_C_Axis, method='trf',bounds=bounds,maxfev=100000)
origin_prediction=[-(beta_fit[1])*np.sin(wind_dir)+Relative_0_Coordinates[0],(beta_fit[1])*np.cos(wind_dir)+Relative_0_Coordinates[1]]
ax.scatter(origin_prediction[0],origin_prediction[1], s=40, color='red',alpha=0.5,label='Prediction')
plt.legend()
plt.draw()

print(f'\nRobot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {i+1} Number Of Measurement Points.\nThe Algorithm Predicts The Origin {round(origin_prediction[0],2),round(origin_prediction[1],2),round(beta_fit[2],2)}. The True Origin Coordinates Were {[round(stack_x,2), round(stack_y,2),round(stack_height,2)]}.\nIt Predicted A Stability Class {stability_class_prediction}, The True Stability Was {stability}.\nIt Predicted An Emission Rate {round(beta_fit[0],2)}, With The True Emission Rate Being {round(emission_rate,2)}{Error_String}')
plt.ioff()
plt.show()

#PLOTTING CONCENTRATION VS. CENTERLINE

plt.title('Concentration Along Centerline')
plt.xlabel('Centerline Axis')
plt.ylabel('Concentration')
plt.scatter(CvsX_X_Axis,CvsX_C_Axis,color='red',label='Measurements')
plt.scatter(beta_fit[1],0,label='Origin',color='green')
XPLOT=np.linspace(-10,100,100)
plt.plot(XPLOT,cGP(XPLOT,*beta_fit),label='Fitted Parameters')
plt.legend()
plt.show()