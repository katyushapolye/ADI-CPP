import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import csv
import re
import mpl_toolkits

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
error = []
it = []
for files in listdir('Fields/VectorFields'):
    filename = files.split('.')[0]
    c = int(filename.split("_")[1])
    u,v =  readCSV("Fields/VectorFields/VectorField_{}.csv".format(c))



    U = u.reshape((N,N));
    V = v.reshape((N,N))

    Ue,Ve = readCSV("Fields/VectorFieldsExact/VectorField_{}.csv".format(c))
    Ue = Ue.reshape((N,N))
    Ve = Ve.reshape((N,N))


    Pxe, Pye = readCSV("Fields/GradPressure/GradPressureField_{}.csv".format(c))



    error.append(np.max(np.abs((Ue - U))))
    #error.append(np.max(np.abs((Ve - V))))

    
    it.append(c)



    mag = np.sqrt(U**2,V**2)
    normalization = Normalize(vmin=0,vmax=0.5)


    x = np.linspace(0, 1, N)
    y = np.linspace(0 ,1, N)
    X, Y = np.meshgrid(x, y)

"""

    plt.figure(figsize=(8,8))
    plt.quiver(X,Y,U,V,mag,cmap='viridis',scale=5,norm=normalization)
    plt.title("Taylor-Green Vortex Aprox {}x{} - IT {}".format(N,N,c));

    plt.savefig('Frames/VectorFrames/VectorFrame_{}.png'.format(c));
    plt.close() #saves mem


    plt.figure(figsize=(8,8))
    plt.quiver(X,Y,Ue,Ve,mag,cmap='viridis',scale=5,norm=normalization)
    plt.title("Taylor-Green Vortex Exact {}x{} - IT {}".format(N,N,c));

    plt.savefig('Frames/VectorFramesExact/VectorFrame_{}.png'.format(c));
    plt.close() #saves mem


    #Pressure plot
    #plt.figure(figsize=(8,8))
    #plt.quiver(X,Y,Pxe,Pye,mag,cmap='viridis',scale=10,norm=normalization)
    #plt.title("Pressure Gradient {}x{} - IT {}".format(N,N,c));
    #plt.savefig('Frames/GradPressureFrames/VectorFrame_{}.png'.format(c));
    #plt.close() #saves mem

    print("Plotting {:.2f}% ".format((p/ len(listdir('Fields/VectorFields')))*100),end='\r',)
    


    p+=1
   

    import gc
    gc.collect()
"""



print("\nError:  " + str(max(error)))

file = open("Error/erros_vectorial_abs.csv",'a')
file.write("{},{},{},{}\n".format(N,1/N,1/N,max(error)))
file.close()




