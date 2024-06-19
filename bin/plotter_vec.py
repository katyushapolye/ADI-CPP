import sys
import numpy as np
import matplotlib.pyplot as plt
import csv
import re

from os import listdir

def readCSV(filename):
    data = []

    # Read the CSV file
    with open(filename, 'r') as csvfile:
        # Create a CSV reader object
        csvreader = csv.reader(csvfile, delimiter=',')

        # Iterate over each row in the CSV
        for row in csvreader:
            # Append the row to the data list after splitting each element on comma
            data.append([elem for elem in row])
    x_values = []
    y_values = []

    # Split each string and extract values
    for row in data:
        for entry in row:
            x_values.append(float(entry.split(';')[0]))
            y_values.append(float(entry.split(';')[1]))

    x_values = np.array(x_values)            
    y_values = np.array(y_values)
    




    return x_values,y_values

N = int(sys.argv[1])

p = 1

for files in listdir('VectorFields'):
    filename = files.split('.')[0]
    c = int(filename.split("_")[1])
    u,v =  readCSV("VectorFields/VectorField_{}.csv".format(c))
    U = u.reshape((N,N));
    V = v.reshape((N,N))


    norm = np.sqrt(U**2,V**2)



    x = np.linspace(-1, 1, N)
    y = np.linspace(-1, 1, N)
    X, Y = np.meshgrid(x, y)




    plt.figure(figsize=(8,8))
    plt.quiver(X,Y,U,V,norm,cmap='viridis',scale=30)
    plt.title("Taylor-Green Vortex - IT {}".format(c));
    print("Plotting {:.2f}% ".format((p/ len(listdir('VectorFields')))*100),end='\r',)



    plt.savefig('VectorFrames/VectorFrame_{}.png'.format(c));
    plt.close() #saves mem

    import gc
    gc.collect()
    p+=1

