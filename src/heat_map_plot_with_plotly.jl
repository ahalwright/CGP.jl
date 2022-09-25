using PlotlyJS
using DataFrames
using CSV
using PyCall
using SciPy

export heatmap_density_plot

# Example:  heatmap_density_plot( pdf.kcomp, pdf.log_redund, "../data/9_23_22/heat_map.png",xlab="log redundancy",ylab="K complexity",toShowScale=true)
# Example:  heatmap_density_plot( pdf.kcomp, pdf.log_redund, "../data/9_23_22/heat_map.png",xlab="log redundancy",ylab="K complexity",toShowScale=true, html_loc="../data/9_23_22/heat_map.html" )
# The procedure for generating the pdf dataframe is given at the end of notes/8_15_22.txt.

function heatmap_density_plot(xvar::Vector{}, yvar::Vector{}, pngsaveloc::String;
     ylab::String = "y",xlab::String = "x", colorsc::String = "Reds", toShowScale::Bool = true,
      background_color::String = "rgba(0,0,0,0)", to_html::Bool=false, html_loc::String="")

    np = pyimport("numpy")

    #Stacks the vectors vertically
    xy = np.vstack([xvar, yvar])

    #Computes weights for each point
    z = SciPy.stats.gaussian_kde(xy)(xy)
    
    #Determines most stylistic elements of the plot
    layout = Layout(
    xaxis_title=xlab,
    yaxis=attr(
        showline=true, 
        linewidth=2, 
        linecolor="black", 
        mirror=true, 
        showgrid=false, 
        zeroline = false,
        ticks="outside",
        tickwidth=1
    ),
    xaxis=attr(
        showline=true, 
        linewidth=2, 
        linecolor="black", 
        mirror=true, 
        showgrid=false, 
        zeroline=false,
        ticks="outside",
        tickwidth=1
    ),
    yaxis_title=ylab,
    colorbar_title="Density",
    font=attr(
        family="Times New Roman",
        size=18,
    ),
    plot_bgcolor=background_color
    )

    #plots the points and creats a trace to be added to a plot
    trace1 = scattergl(
    x=xvar,
    y=yvar,
    mode="markers",
    marker=attr(color = z, colorscale=colorsc, showscale=toShowScale)
    )

    pconfig=PlotConfig(responsive=true, displayModeBar=false, scrollZoom=true)
    #Constructs the plot
    fig = plot(trace1, layout, config=pconfig)

    #exports a png of the plot 
    savefig(fig, pngsaveloc)

    if length(html_loc) > 0
        open(html_loc, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end 
        println("[+] HTML Plot succesfully created and saved to " * html_loc)
    end
    print("[+] Plot succesfully created and saved to " * pngsaveloc)
end


