## Bug on setting the axis labels 10/5/22  NOT CORRECT
## Relies upon DataFrames, PyCall, CSV, PyPlot and SciPy
using PyCall, PyPlot, SciPy
## Data needs at least two named columns
## Example function call:
## heatmap( "../data/counts/count_outputs_ch_4funcs_3inputs_7gates_4lb_cmplxC.csv", "../data/9_4_22/heat_map_plot.png",  "lgints7_4", "k_complexity" )
function heatmap(dataloc::String, saveloc::String, xvar::String, yvar::String; ylab::String="y", xlab::String="x", labfontsize::Int=20)
  hdf = read_dataframe( dataloc )
  heatmap( hdf.xvar, hdf.yvar, ylab, xlab, labfontsize, csvfile=saveloc )
end

function heatmap(x::Vector, y::Vector, ylab::String="y", xlab::String="x", labfontsize::Int=20; csvfile::String="" )
    #load data in using Dataframe casting
    #Importing python numpy module for vstack function
    np = pyimport("numpy")

    #Creating heat density plot with matlab plotting method
    #x = data[:,xvar]
    #y = data[:,yvar]

    #Calculating point density

    #Stacks the vectors vertically
    xy = np.vstack([x, y])

    #Computes weights for each point
    z = SciPy.stats.gaussian_kde(xy)(xy)

    #Customizing the plot
    fig, ax = PyPlot.subplots()
    PyPlot.scatter(x, y, c = z, s=10, edgecolor="none", cmap = "Reds")

    PyPlot.tick_params(labelsize=15)
    PyPlot.xticks(size = 15, family = "Times New Roman")
    PyPlot.yticks(size = 15, family = "Times New Roman")
    PyPlot.xlabel(xlab, 20, family = "Times New Roman")
    PyPlot.ylabel(ylab, 20, family = "Times New Roman")
    cb = PyPlot.colorbar(shrink = .5)
    cb.ax.tick_params(labelsize = 15)
    for i in cb.ax.yaxis.get_ticklabels()
        i.set_family("Times New Roman")
    end

    PyPlot.figtext(0.76, 0.73, "density", size = 20, family = "Times New Roman")

    PyPlot.savefig(saveloc)
    println("Saved a heatmap density plot to " * saveloc)
end
