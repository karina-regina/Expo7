import numpy as np                              
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sympy as sym

#input from MatLab
x_sample = np.array([-14.0, 10.0, 34.0, 88.0, 64.0, 46.0, 88.0, 54.4])   
y_sample = np.array([5.0, 6.0, 22.0, 7.0, -38.0, 12.0, 10.0, 20.0])
windspeed = np.array([3.1, 5.0, 1.9, 6.2, 4.7, 6.9, 2.4, 5.0])
y = np.zeros(np.size(x_sample))
z = np.zeros(np.size(y_sample))
stability_class = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
x_sample, y_sample, y, z, windspeed, stability_class = data

# Define the function for the Guassian plume model
def func(data, x_leakage, y_leakage, z_leakage, emission_rate):
    list_a, list_b, list_c, list_d = ([] for i in range(4))
    for i in range(len(stability_class)):
        stability_values= {1 : (122.8,0.94470,24.1670,2.5334),
                  2: (90.673, 0.93198, 18.3330, 1.8096),
                  3: (61.141, 0.91465, 12.5, 1.0857),
                  4: (34.459, 0.86974, 8.3330, 0.72382),
                  5: (24.26, 0.83660, 6.25, 0.54287),
                  6: (15.209, 0.81558, 4.1667, 0.36191)
                 }
        list_a.append(stability_values[stability_class[i]][0])
        list_b.append(stability_values[stability_class[i]][1])
        list_c.append(stability_values[stability_class[i]][2])
        list_d.append(stability_values[stability_class[i]][3])

    distance = np.sqrt((x_sample - x_leakage)**2 + (y_sample - y_leakage)**2)                                             # since our y value is always zero, we fill in the distance for x in the formula (read user manual to understand exactly why). this is the distance between the leakage and the datapoint
    theta = 0.017453293 * (list_c - list_d * np.log(distance/1000))                         # here you can see the calculation of: theta
    sigma_y = 465.11628 * (distance / 1000) * np.tan(theta)                                 # sigma_y
    sigma_z = list_a * (distance / 1000)**list_b                                            # sigma_Z
    first_fraction = Q / (2 * np.pi * U * sigma_y * sigma_z)                                # the first fracture of the formula
    exp_arg_z1 = -(z - h)**2 / (2 * sigma_z**2)                                             # the first exponential of the formula
    exp_arg_z2 = -(z + h)**2 / (2 * sigma_z**2)                                             # the second exponential of the formula
    exp_arg_y = -(y)**2 / (2 * sigma_y**2)                                                  # the last exponential of the formula
    eq = first_fraction * (np.exp(exp_arg_z1) + np.exp(exp_arg_z2)) * np.exp(exp_arg_y)     # here the fracture and exponentials come together to create the entire equation

    return eq

# data from robot
# left side of the equation                                     here are the values for C that are retrieved from the robot
C = np.array([1.2906*(10**-5), 1.7442*(10**-22), 3.9516*(10**-7), 2.9693*(10**-4), 2.5659*(10**-4), 5.3601*(10**-7), 7.7545*(10**-4), 3.1532*(10**-5)]) 

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
print('x_leakage =',solution[0],',','y_leakage=',solution[1],',','z_leakage =',solution[2],',','emmision_rate =',solution[3])


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