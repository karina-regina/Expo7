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

LDA_Data=pd.read_excel('6000_LDA_Data.xlsx')
x_cols = ['Crosswind', 'Distance On Centerline', 'Negative Gradient','Positive Gradient','Stability Class']


#print(LDA_Data[x_cols].cov())
model = LinearDiscriminantAnalysis(solver='lsqr').fit(LDA_Data[x_cols], LDA_Data['Stability Class'])
stability_class_prediction = model.predict(LDA_Data[x_cols])
print(len(np.where(abs(stability_class_prediction-LDA_Data['Stability Class'])<=2)[0]),len(stability_class_prediction))
