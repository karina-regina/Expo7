# import libraries that are used further in the code (please do not touch these four lines)
import numpy as np                              
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sympy as sym

# Define the function for the Guassian plume model
def func(data, f, g, h, Q):                             # f is here the x coördinate of the leakage, g is the y coördinate of the leakage, h is the height of the leakage and Q is the mass emission rate in a certain point
    v, w, y, z, U, stability = data                     # These are the data that are given by the robot
    list_a = []                                         # v is here the x coördinate of the datapoint, w is the y coördinate of the data point
    list_b = []                                         # y & z are the y and z as in the formula and are always 0 in our case
    list_c = []                                         # U is the windspeed and stability is the pasquil stability category which is a value from 1-6, where 1 = A and 6 = F
    list_d = []                                         # these four lists for a,b,c and d will create the values for these variables for a datapoint, when the stability category is given
    for i in range(len(stability)):
        stabilityvalues= {1 : (122.8,0.94470,24.1670,2.5334),
                  2: (90.673, 0.93198, 18.3330, 1.8096),
                  3: (61.141, 0.91465, 12.5, 1.0857),
                  4: (34.459, 0.86974, 8.3330, 0.72382),
                  5: (24.26, 0.83660, 6.25, 0.54287),
                  6: (15.209, 0.81558, 4.1667, 0.36191)
                 }
        list_a.append(stabilityvalues[stability[i]][0])
        list_b.append(stabilityvalues[stability[i]][1])                           # y & z are the y and z as in the formula and are always 0 in our case
        list_c.append(stabilityvalues[stability[i]][2])                                         # U is the windspeed and stability is the pasquil stability category which is a value from 1-6, where 1 = A and 6 = F
        list_d.append(stabilityvalues[stability[i]][3])

    distance = np.sqrt((v - f)**2 + (w - g)**2)                                             # since our y value is always zero, we fill in the distance for x in the formula (read user manual to understand exactly why). this is the distance between the leakage and the datapoint
    theta = 0.017453293 * (list_c - list_d * np.log(distance/1000))                         # here you can see the calculation of: theta
    sigma_y = 465.11628 * (distance / 1000) * np.tan(theta)                                 # sigma_y
    sigma_z = list_a * (distance / 1000)**list_b                                            # sigma_Z
    first_fracture = Q / (2 * np.pi * U * sigma_y * sigma_z)                                # the first fracture of the formula
    exp_arg_z1 = -(z - h)**2 / (2 * sigma_z**2)                                             # the first exponential of the formula
    exp_arg_z2 = -(z + h)**2 / (2 * sigma_z**2)                                             # the second exponential of the formula
    exp_arg_y = -(y)**2 / (2 * sigma_y**2)                                                  # the last exponential of the formula
    eq = first_fracture * (np.exp(exp_arg_z1) + np.exp(exp_arg_z2)) * np.exp(exp_arg_y)     # here the fracture and exponentials come together to create the entire equation

    return eq

# data from robot
# left side of the equation                                     here are the values for C that are retrieved from the robot
C = np.array([1.2906*(10**-5), 1.7442*(10**-22), 3.9516*(10**-7), 2.9693*(10**-4), 2.5659*(10**-4), 5.3601*(10**-7), 7.7545*(10**-4), 3.1532*(10**-5)]) 

# right side of the equation                                                These are all values that should be retrieved from the robot
v = np.array([-14.0, 10.0, 34.0, 88.0, 64.0, 46.0, 88.0, 54.4])             # values for v (the x coördinate of the datapoint)   
w = np.array([5.0, 6.0, 22.0, 7.0, -38.0, 12.0, 10.0, 20.0])                # values for w (the y coördinate of the datapoint)
U = np.array([3.1, 5.0, 1.9, 6.2, 4.7, 6.9, 2.4, 5.0])                      # values for the U (the winspeed measured at the datapoint)
y = np.zeros(np.size(v))                                                    # y is always 0 in our case
z = np.zeros(np.size(v))                                                    # z is always 0 in our case
stability = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0])              # stability per datapoint
x = v, w, y, z, U, stability                                                # collection of all the data, so it can be used in the curvefit function

# loop that uses ceveral initial guesses and picks the best one
solution = [0,0,0,0]                                                        # an empty solution that will be replaced with new solutions by running the code
diff_gues = 100000                                                          # chose a high number, such that the solution will be changed the first time
for i in range(21):                                                         # for loops over f,g,h and Q with estimations of the solution ("initial guesses")
    f_gues = i
    for j in range(2):
        g_gues = j
        for p in range(8):
            h_guess = p
            for q in range(2):
                q_gues = q
                p0 = [f_gues, g_gues, h_guess, q_gues]                                                              # create the initial guess
                popt, pcov = curve_fit(func, x, C, p0, maxfev = 100000, method = 'dogbox')                          # the function that finds the solutions for f,g,h and Q
                diff_guess = abs(p0[0] - popt[0])+abs(p0[1] - popt[1])+abs(p0[2] - popt[2])+abs(p0[3] - popt[3])    # calculates the difference between the guess and the solution to determine how good the solution is (why? read user manual)
                if diff_guess < diff_gues:                                                                          # if this solution is better then the previous best solution, then change the solution to this one
                    diff_gues = diff_guess
                    solution = popt

# show the best found values for f, g, h and Q
print('f =',solution[0],',','g =',solution[1],',','h =',solution[2],',','Q =',solution[3])                          # this prints the output




# The code has ended here!
# this is just some information to get to understand the code better and to show the results
# The values that where used for f, g, h and Q in this example that the model should find are:
# f = 20 g = 0 h = 7 Q = 1
# what the model found with 3 datapoints = f = 0.0 , g = 0.0 , h = 3.0 , Q = 0.0
# what the model found with 4 datapoints = f = 17.0 , g = 0.0 , h = 7.0 , Q = 1.0
# what the model found with 5 datapoints = f = 19.0 , g = 0.0 , h = 7.0 , Q = 1.0
# what the model found with 6 datapoints = f = 19.0 , g = 0.0 , h = 7.0 , Q = 1.0
# what the model found with 7 datapoints = f = 20.0 , g = 0.0 , h = 7.0 , Q = 1.0
# what the model found with 8 datapoints = f = 20.0 , g = 0.0 , h = 7.0 , Q = 1.0