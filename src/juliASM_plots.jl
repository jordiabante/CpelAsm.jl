"""
plot_single_histogram(x,title)

Function to generate a histogram of a vectors x.

# Examples
```julia-repl
julia>p=juliASM.plot_single_histogram(x,title)
```
"""
function plot_single_histogram(x,title; xlabel=L"Value",ylabel=L"Frequency")
    xmin = 0.9*minimum(reduce(vcat,x))
    xmax = 1.1*maximum(reduce(vcat,x))
    p = histogram(x, normalize=true, line=(1,0.2), xlabel=xlabel,
      ylabel=ylabel, fillcolor=:blue, fillalpha=0.2, xtickfont=font(6,
      "Arial"), title=title, titlefontsize=12, xlims=(xmin,xmax),
      nbins=range(xmin,xmax,length=25)
  )
  return p
end # end plot_single_histogram
"""
plot_histogram(x, labels)

Function to generate an overlayed histogram of vectors x[1], x[2], ....

# Examples
```julia-repl
julia>p=juliASM.plot_histogram(x, labels)
```
"""
function plot_histogram(x, labels; xlabel=L"Value",ylabel=L"Frequency")
    xmin = 0 #0.9 * minimum(reduce(vcat,x))
    xmax = 0.5 #1.1 * maximum(reduce(vcat,x))
    # x_colors = [:blue :red :green :orange :pink :]
    # x_colors = x_colors[1:length(labels)]
    p = histogram(x, label=labels, normalize=true, line=(1,0.2), xlabel=xlabel,
        ylabel=ylabel, fillcolor=:match, fillalpha=0.25, xlims=(xmin,xmax),
        xtickfont=font(6, "Arial"), nbins=range(xmin,xmax,length=25)
        # title=title, titlefontsize=12,
  )
  return p
end # end plot_histogram
# """
# plot_boxplot(x, labels)
#
# Function to generate an boxplot of vectors x[1], x[2], ....
#
# # Examples
# ```julia-repl
# julia>p=juliASM.plot_histogram(x, labels)
# ```
# """
# function plot_boxplot(x, labels; xlabel=L"Label",ylabel=L"Values")
#     ymin = 0 #0.9 * minimum(reduce(vcat,x))
#     ymax = 0.5 #1.1 * maximum(reduce(vcat,x))
#     # x_colors = [:blue :red :green :orange :pink :]
#     # x_colors = x_colors[1:length(labels)]
#     p = boxplot(labels, x, leg=false, line=(1,0.2), xlabel=xlabel, ylabel=ylabel,
#         fillcolor=:match, ylims=(ymin,ymax), xtickfont=font(6, "Arial"),
#         fillalpha=0.25
#         # title=title, titlefontsize=12,
#   )
#   return p
# end # end plot_boxplot
