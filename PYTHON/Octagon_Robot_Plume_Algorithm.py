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

#set from 1-6
stability = 1                                                                                      
#Origin Coordinates Of Plume
angle= np.random.uniform(0, 2 * np.pi)
stack_x = np.random.uniform(0, 51) * np.cos(angle)
stack_y = np.random.uniform(0, 51) * np.sin(angle)
# stack_x = -100
# stack_y = 100
stack_height = 3 #np.random.randint(1,11) 
emission_rate = 10 #np.random.randint(1,21)                                                                 
windspeed= 4.7 #np.random.randint(1,26)
#Wind Angle In Degrees
wind_direction= 45 #np.random.randint(0,360)

#INPUT FOR ALGORITHM

#Input For Coordinates Of Initial Measurement Point Of Robot (Can Start Out With Multiple)
Measurement_Points=[[0,0]] 

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

wind_xcomponent=windspeed*math.sin(math.radians(wind_direction-180))
wind_ycomponent=windspeed*math.cos(math.radians(wind_direction-180))

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
     for l in range(len(Measurement_Points)):
          ax.scatter(Measurement_Points[l][0], Measurement_Points[l][1], s=5,alpha=0.5, color='white')
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
#This Is Distance Between Gas Leak And Robot
RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y])
#Adds Concentration To Concentration_At_Measurement For Each Coordinate.]
for XYCoordinate in Measurement_Points:
     Concentration_At_Measurement.append(Concentration2D[((round(100+XYCoordinate[1]))*2),((round(100+XYCoordinate[0]))*2)])

#This Is Where Non-Iterative Calculations For The Algorithm Belong 

# octax = np.array([0,31.82,0, -31.82, -45, -31.82, 0, 31.82,45,75, 53.03, 0, -53.03, -75, -53.03, 0, 53.03 ])
# octay = np.array([0,31.82, 45, 31.82, 0, -31.82, -45, -31.82,0,0, 53.03, 75, 53.03, 0, -53.03, -75, -53.03])
#small to big 
# octax = np.array([0,55.0, 47.63139720814413, 27.500000000000007, 3.3677786976552213e-15, -27.49999999999999, -47.63139720814413, -55.0, -47.63139720814412, -27.500000000000025, -1.0103336092965664e-14, 27.500000000000007, 47.63139720814411, 33.25, 20.73103591180289, -7.398821054047452, -29.95721485775543, -29.957214857755435, -7.39882105404743, 20.731035911802906, 77.75, 72.49971583618867, 57.45794331390624, 34.656157161625856, 7.173864948271731, -21.277297478104426, -46.85484297848718, -66.10438230297748, -76.42615850042337, -76.42615850042337, -66.1043823029775, -46.854842978487135, -21.277297478104398, 7.173864948271764, 34.65615716162587, 57.4579433139063, 72.49971583618868, 22.125, 6.837001000545713, -17.89950100054571, -17.899501000545712, 6.837001000545707, 44.375, 35.90012912538829, 13.712629125388293, -13.712629125388288, -35.90012912538829, -44.375, -35.90012912538829, -13.712629125388299, 13.712629125388283, 35.90012912538829, 66.625, 60.02705082399868, 41.540008048838125, 14.82545722483945, -14.825457224839443, -41.54000804883814, -60.027050823998685, -66.625, -60.02705082399866, -41.54000804883808, -14.8254572248394, 14.825457224839493, 41.54000804883815, 60.02705082399869, 88.875, 84.0595073561439, 70.13486277260448, 48.610017553130696, 21.81752266963853, -7.33923932885353, -35.70068086603264, -60.19339967823772, -78.16322963847671, -87.66286083991695, -87.66286083991696, -78.16322963847674, -60.19339967823777, -35.7006808660327, -7.339239328853572, 21.817522669638493, 48.61001755313066, 70.13486277260442, 84.05950735614387, 16.5625, 1.0141606305439018e-15, -16.5625, -3.0424818916317053e-15, 27.6875, 13.843750000000004, -13.843749999999995, -27.6875, -13.843750000000012, 13.843750000000004, 38.8125, 27.444581944802877, 2.376580194595332e-15, -27.444581944802874, -38.8125, -27.44458194480288, -7.129740583785997e-15, 27.44458194480287, 49.9375, 42.01009829575711, 20.744787211781706, -7.106847236272181, -32.70210790139267, -47.914680369873956, -47.91468036987397, -32.702107901392694, -7.106847236272181, 20.744787211781723, 42.01009829575711, 61.0625, 54.06815856644913, 34.6874535972712, 7.360271038090661, -21.653060915034835, -45.705937560197846, -59.28813472657829, -59.28813472657831, -45.705937560197896, -21.6530609150349, 7.360271038090565, 34.6874535972711, 54.068158566449085, 72.1875, 66.6925537531585, 51.0442707669039, 27.624960273854924, 4.420209540672478e-15, -27.624960273854914, -51.04427076690389, -66.6925537531585, -72.1875, -66.69255375315852, -51.044270766903914, -27.6249602738549, -1.3260628622017433e-14, 27.624960273854935, 51.044270766903885, 66.69255375315852, 83.3125, 78.28814146922599, 63.821077667349854, 41.65625000000001, 14.46706380187614, -14.467063801876131, -41.65624999999998, -63.82107766734985, -78.28814146922599, -83.3125, -78.28814146922599, -63.821077667349854, -41.656250000000036, -14.467063801876133, 14.467063801876103, 41.65625000000001, 63.82107766734984, 78.28814146922599, 94.4375, 89.81539975762357, 76.4015424062841, 55.50896976337043, 29.1827924062841, 5.782629104723908e-15, -29.18279240628409, -55.50896976337042, -76.4015424062841, -89.81539975762357, -94.4375, -89.81539975762357, -76.4015424062841, -55.508969763370445, -29.18279240628411, -1.7347887314171724e-14, 29.18279240628408, 55.50896976337041, 76.4015424062841, 89.81539975762357])
# octay = np.array([0,0.0, 27.499999999999996, 47.63139720814412, 55.0, 47.63139720814413, 27.499999999999996, 6.735557395310443e-15, -27.500000000000007, -47.631397208144115, -55.0, -47.63139720814412, -27.500000000000025, 0.0, 25.99589679206199, 32.41635308004564, 14.426634325658812, -14.426634325658803, -32.41635308004564, -25.995896792061977, 0.0, 28.08653954605114, 52.37983629351982, 69.5989459028561, 77.41833220693893, 74.78194375668669, 62.04583942103863, 40.930100663714434, 14.286525010238346, -14.28652501023836, -40.93010066371441, -62.04583942103866, -74.7819437566867, -77.41833220693893, -69.59894590285609, -52.379836293519766, -28.086539546051082, 0.0, 21.04212542303027, 13.00474870697097, -13.004748706970966, -21.042125423030274, 0.0, 26.082970570478494, 42.20313291059744, 42.203132910597446, 26.0829705704785, 5.43437017121638e-15, -26.08297057047849, -42.20313291059744, -42.203132910597446, -26.082970570478505, 0.0, 28.90750411870731, 52.089522519432485, 64.954572149114, 64.954572149114, 52.08952251943247, 28.907504118707287, -2.1428234306941183e-14, -28.907504118707354, -52.08952251943251, -64.954572149114, -64.95457214911399, -52.08952251943246, -28.90750411870727, 0.0, 28.85766532556624, 54.58815484029423, 74.40317075558222, 86.15544863535798, 88.57144681596779, 81.38935440646823, 65.38746256107459, 42.29982455616992, 14.628344211200238, -14.628344211200215, -42.29982455616987, -65.38746256107454, -81.38935440646821, -88.57144681596779, -86.155448635358, -74.40317075558224, -54.588154840294315, -28.85766532556634, 0.0, 16.5625, 2.0283212610878037e-15, -16.5625, 0.0, 23.978078367281643, 23.978078367281647, 3.390740825139234e-15, -23.97807836728164, -23.978078367281643, 0.0, 27.444581944802877, 38.8125, 27.444581944802877, 4.753160389190664e-15, -27.444581944802874, -38.8125, -27.44458194480288, 0.0, 26.998250821688902, 45.42474776801626, 49.429208253929076, 37.740244369315775, 14.069019557268913, -14.06901955726888, -37.740244369315754, -49.429208253929076, -45.42474776801625, -26.9982508216889, 0.0, 28.377158692922617, 50.2534523111314, 60.61728562461242, 57.09442931897814, 40.49192731882857, 14.613212750559017, -14.61321275055895, -40.49192731882852, -57.094429318978115, -60.617285624612435, -50.25345231113147, -28.377158692922702, 0.0, 27.624960273854917, 51.0442707669039, 66.6925537531585, 72.1875, 66.6925537531585, 51.0442707669039, 27.624960273854928, 8.840419081344956e-15, -27.62496027385491, -51.04427076690389, -66.69255375315852, -72.1875, -66.6925537531585, -51.044270766903914, -27.624960273854903, 0.0, 28.494553190819776, 53.5522427320098, 72.15074145279104, 82.04679592282957, 82.04679592282957, 72.15074145279105, 53.55224273200982, 28.49455319081979, 1.0202838645396386e-14, -28.49455319081977, -53.5522427320098, -72.15074145279102, -82.04679592282957, -82.04679592282959, -72.15074145279104, -53.55224273200983, -28.494553190819765, 0.0, 29.182792406284094, 55.50896976337043, 76.4015424062841, 89.81539975762357, 94.4375, 89.81539975762357, 76.4015424062841, 55.508969763370445, 29.182792406284104, 1.1565258209447816e-14, -29.182792406284083, -55.50896976337042, -76.4015424062841, -89.81539975762357, -94.4375, -89.81539975762357, -76.40154240628411, -55.50896976337045, -29.182792406284122])

#small big medium to medium big? 
octax = np.array([0,55.0, 47.631, 27.5, 0.0, -27.5, -47.631, -55.0, -47.631, -27.5, -0.0, 27.5, 47.631, 33.25, 20.731, -7.399, -29.957, -29.957, -7.399, 20.731, 77.75, 72.5, 57.458, 34.656, 7.174, -21.277, -46.855, -66.104, -76.426, -76.426, -66.104, -46.855, -21.277, 7.174, 34.656, 57.458, 72.5, 22.125, 6.837, -17.9, -17.9, 6.837, 88.875, 84.06, 70.135, 48.61, 21.818, -7.339, -35.701, -60.193, -78.163, -87.663, -87.663, -78.163, -60.193, -35.701, -7.339, 21.818, 48.61, 70.135, 84.06, 44.375, 35.9, 13.713, -13.713, -35.9, -44.375, -35.9, -13.713, 13.713, 35.9, 66.625, 60.027, 41.54, 14.825, -14.825, -41.54, -60.027, -66.625, -60.027, -41.54, -14.825, 14.825, 41.54, 60.027, 16.562, 0.0, -16.562, -0.0, 94.438, 89.815, 76.402, 55.509, 29.183, 0.0, -29.183, -55.509, -76.402, -89.815, -94.438, -89.815, -76.402, -55.509, -29.183, -0.0, 29.183, 55.509, 76.402, 89.815, 27.688, 13.844, -13.844, -27.688, -13.844, 13.844, 83.312, 78.288, 63.821, 41.656, 14.467, -14.467, -41.656, -63.821, -78.288, -83.312, -78.288, -63.821, -41.656, -14.467, 14.467, 41.656, 63.821, 78.288, 38.812, 27.445, 0.0, -27.445, -38.812, -27.445, -0.0, 27.445, 72.188, 66.693, 51.044, 27.625, 0.0, -27.625, -51.044, -66.693, -72.188, -66.693, -51.044, -27.625, -0.0, 27.625, 51.044, 66.693, 49.938, 42.01, 20.745, -7.107, -32.702, -47.915, -47.915, -32.702, -7.107, 20.745, 42.01, 61.062, 54.068, 34.687, 7.36, -21.653, -45.706, -59.288, -59.288, -45.706, -21.653, 7.36, 34.687, 54.068])
octay = np.array([0,0.0, 27.5, 47.631, 55.0, 47.631, 27.5, 0.0, -27.5, -47.631, -55.0, -47.631, -27.5, 0.0, 25.996, 32.416, 14.427, -14.427, -32.416, -25.996, 0.0, 28.087, 52.38, 69.599, 77.418, 74.782, 62.046, 40.93, 14.287, -14.287, -40.93, -62.046, -74.782, -77.418, -69.599, -52.38, -28.087, 0.0, 21.042, 13.005, -13.005, -21.042, 0.0, 28.858, 54.588, 74.403, 86.155, 88.571, 81.389, 65.387, 42.3, 14.628, -14.628, -42.3, -65.387, -81.389, -88.571, -86.155, -74.403, -54.588, -28.858, 0.0, 26.083, 42.203, 42.203, 26.083, 0.0, -26.083, -42.203, -42.203, -26.083, 0.0, 28.908, 52.09, 64.955, 64.955, 52.09, 28.908, -0.0, -28.908, -52.09, -64.955, -64.955, -52.09, -28.908, 0.0, 16.562, 0.0, -16.562, 0.0, 29.183, 55.509, 76.402, 89.815, 94.438, 89.815, 76.402, 55.509, 29.183, 0.0, -29.183, -55.509, -76.402, -89.815, -94.438, -89.815, -76.402, -55.509, -29.183, 0.0, 23.978, 23.978, 0.0, -23.978, -23.978, 0.0, 28.495, 53.552, 72.151, 82.047, 82.047, 72.151, 53.552, 28.495, 0.0, -28.495, -53.552, -72.151, -82.047, -82.047, -72.151, -53.552, -28.495, 0.0, 27.445, 38.812, 27.445, 0.0, -27.445, -38.812, -27.445, 0.0, 27.625, 51.044, 66.693, 72.188, 66.693, 51.044, 27.625, 0.0, -27.625, -51.044, -66.693, -72.188, -66.693, -51.044, -27.625, 0.0, 26.998, 45.425, 49.429, 37.74, 14.069, -14.069, -37.74, -49.429, -45.425, -26.998, 0.0, 28.377, 50.253, 60.617, 57.094, 40.492, 14.613, -14.613, -40.492, -57.094, -60.617, -50.253, -28.377])
x= np.diff(octax)
y= np.diff(octay)

while RadiusCheck*Redundancy >= 2:   
    RadiusCheck=math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]),[stack_x,stack_y]) 
    print(Measurement_Points[-1], '\nc', Concentration_At_Measurement[-1])
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
        move_direction=[0,0,positive_condition,negative_condition]
        
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
             break
        #Rounds To The Nearest Half. (Due To Resolution Of Gas Plume Being In 0.5 By 0.5 Meters)
        Measurement_Points[-1][j]=(round(Measurement_Points[-1][j]*2)/2)
    if len(Measurement_Points)>=2:
        Total_Distance+=math.dist(Measurement_Points[i],Measurement_Points[i-1])

    #UPDATES TO PLOT & WHILE-LOOP

    update_annotation(f'Current Distance Between Last Measurement Point And Objective Is {round(math.dist((Measurement_Points[-1][0], Measurement_Points[-1][1]), [stack_x, stack_y]),ndigits =2)} Meters\nRobot Has Travelled {round(Total_Distance,ndigits=2)} Meters & Taken {i+1} Amount Of Measurements\nCurrent Concentration At Coordinate Is {round(Concentration_At_Measurement[-1],ndigits=5)}')
    ax.scatter(Measurement_Points[i][0], Measurement_Points[i][1], s=10, color='white',alpha=0.5)
    ax.plot([Measurement_Points[i-1][0],Measurement_Points[i][0]], [Measurement_Points[i-1][1], Measurement_Points[i][1]],alpha=0.5,color='white')
    i+=1
    Concentration_At_Measurement.append(Concentration2D[((round(100+Measurement_Points[i][1]))*2)][((round(100+Measurement_Points[i][0]))*2)])  
    plt.draw()
    plt.pause(0.4)
print(f'Robot Travelled {round(Total_Distance,ndigits=2)} Meters And Took {i} Number Of Measurement Points, Until It Got To {[round(Measurement_Points[-1][0],ndigits=2),round(Measurement_Points[-1][1],ndigits=2)]} With A Concentration Of {Concentration_At_Measurement[-1]}. The Objective Coordinates Were {[stack_x, stack_y]}{Error_String}')
plt.ioff()
plt.show()
