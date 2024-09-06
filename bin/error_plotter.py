
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



log = np.genfromtxt("Error/erros_vectorial_abs.csv",delimiter=',')

dH = log[1:,1]
error = log[1:,3]
print("Calculating error decay")
order = (np.log(error[-1]) - np.log(error[0])) / (np.log(dH[-1]) - np.log(dH[0]))

plt.plot(np.log(dH),np.log(error),marker='*',color='black')
plt.title("Error \u0394t = \u0394h- Conv. = {:.4f}".format(order))

plt.savefig("error_log.png")
