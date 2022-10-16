using PlotlyJS
using DataFrames
using CSV
using PyCall
using SciPy

export heatmap_density_plot

# Example: heatmap_scatter_plot(  pdf.kcomp, pdf.log_redund, xlabl="K complexity", ylabl="log redundancy", pngfile="../data/10_5_22/log_redund_vs_Kcomp_4x1_10gts5lb4funcs_10_5.png" )

function heatmap_scatter_plot( xvar::Vector{}, yvar::Vector{}; xlabl::String="x", ylabl::String="y", pngfile::String="plot.html",
  colorsc::String="Reds", toShowScale::Bool=true, background_color::String="rgba(0,0,0,0)", html_loc::String="" )
  to_html = length(html_loc) > 0 ? true : false
  heatmap_density_plot( xvar, yvar, pngfile, ylab=ylabl, xlab=xlabl, colorsc=colorsc, toShowScale=toShowScale,
      background_color=background_color, to_html=to_html, html_loc=html_loc )
end

function heatmap_density_plot(xvar::Vector{}, yvar::Vector{}, pngsaveloc::String;
      ylab::String = "y",xlab::String = "x", colorsc::String = "Reds", toShowScale::Bool = false,
      background_color::String = "rgba(0,0,0,0)", to_html::Bool=false, html_loc::String="plot.html" )

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
    fig = PlotlyJS.plot(trace1, layout, config=pconfig  )

    if to_html
        open(html_loc, "w") do io
            PlotlyBase.to_html(io, fig.plot)
        end 
        println("[+] HTML Plot succesfully created and saved to " * html_loc)
    end
    print("[+] Plot succesfully created and saved to " * pngsaveloc)

    #exports a png of the plot 
    PlotlyJS.savefig(fig, pngsaveloc)
end
#=
function heatmap_density_plot(csv_path::String, pngsaveloc::String, y_variable_name::String, x_variable_name::String;
    ylab::String = "y",xlab::String = "x", colorsc::String = "Reds", toShowScale::Bool = false,
     background_color::String = "rgba(0,0,0,0)", to_html::Bool=false, html_loc::String="plot.html")

   np = pyimport("numpy")

   #Reads in the 
   data = DataFrame(CSV.File(csv_path))
   xvar = data[!, x_variable_name]
   yvar = data[!, y_variable_name]
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

   if to_html
       open(html_loc, "w") do io
           PlotlyBase.to_html(io, fig.plot)
       end 
       println("[+] HTML Plot succesfully created and saved to " * html_loc)
   end
   print("[+] Plot succesfully created and saved to " * pngsaveloc)
end
=#
