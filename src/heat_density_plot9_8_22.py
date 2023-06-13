# OUTDATED 10/6/22
# Example run:
# python heat_density_plot9_8_22.py -p '../data/2020/10_27/geno_complexity10_27FMNccons.csv' -x 'robustness' -y 'complexity'  -d true -s "test.png"
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import optparse

def get_commands():
       #Reading commands from the line with optparse
       parser = optparse.OptionParser()
       #Core functionality
       parser.add_option("-p", "--p", dest = "csv_path", help = "Path to the csv file that will be used as a dataframe")
       parser.add_option("-x", "--x_column", dest = "x_column", help = "Name of the column of the CSV that is to be put on the x-axis")
       parser.add_option("-y", "--y_column", dest = "y_column", help = "Name of the column of the CSV that is to be put on the y-axis")
       parser.add_option("-s", "--save_path", dest = "save_path", help = "The desired path and file name of the graph to be saved")

       #Additional arguments that help aesthetically
       parser.add_option("--x_label", dest = "x_label", help = "Label for the x-axis")
       parser.add_option("--y_label", dest = "y_label", help = "Label for the y-axis")
       parser.add_option("-f", "--fontsize", dest = "fontsize", help = "Font size for labels")
       parser.add_option("-d", "--display", dest = "display", help = "Boolean value which chooses whether or not to execute the show command")
       (options, arguments) = parser.parse_args()
       return options

arguments = get_commands()
csv_path = arguments.csv_path
save_path = arguments.save_path
x_column = arguments.x_column
y_column = arguments.y_column
y_lab = arguments.y_label
x_lab = arguments.x_label
fontsize = arguments.fontsize

if (not y_lab):
       y_lab = y_column
if (not x_lab):
       x_lab = x_column
if (not fontsize):
       fontsize = 15

data = pd.read_csv(csv_path,comment='#')
y = data.loc[:, y_column]
x = data.loc[:, x_column]

# Calculate the point density
xy = np.vstack([x, y])
z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots()
plt.scatter(x, y, c=z, s=10, edgecolor='none',cmap = 'Reds')
plt.tick_params(labelsize=fontsize)
plt.xticks(fontsize = 15,family = 'Times New Roman')
plt.yticks(fontsize = 15,family = 'Times New Roman')
plt.xlabel(x_column,fontsize = 20,family = 'Times New Roman')
plt.ylabel(y_column,fontsize = 20,family = 'Times New Roman')
cb = plt.colorbar(shrink = 0.5)
cb.ax.tick_params(labelsize=15)
for l in cb.ax.yaxis.get_ticklabels():
       l.set_family('Times New Roman')
plt.figtext(0.76,0.73,'density',size = 20,family = 'Times New Roman')
plt.savefig(save_path)     
print("[+] Plot using " + csv_path + " with x = " + x_column + " and y = " + y_column + " saved to " + save_path)

if(arguments.display):
       plt.show()
