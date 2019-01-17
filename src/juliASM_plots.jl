"""
plot_histogram(x, labels)

Function to generate an overlayed histogram of vectors x[1], x[2], ....

# Examples
```julia-repl
julia>p=JuliASM.plot_histogram(x, labels)
```
"""
function plot_histogram(x, labels; xlabel=L"Value",ylabel=L"Frequency")
    xmin = 0.9 * minimum(reduce(vcat,x))
    xmax = 1.1 * maximum(reduce(vcat,x))
    # x_colors = [:blue :red :green :orange :pink :]
    # x_colors = x_colors[1:length(labels)]
    p = histogram(x, label=labels, normalize=true, line=(1,0.2), xlabel=xlabel,
        ylabel=ylabel, fillcolor=:match, fillalpha=0.25, xlims=(xmin,xmax),
        xtickfont=font(6, "Arial"), nbins=range(xmin,xmax,length=25)
        # title=title, titlefontsize=12,
        )
    return p
end # end plot_histogram
"""
plot_heatmap(x)

Function to generate an heatmap of log-likelihood function for a gixen
observation vector x.

# Examples
```julia-repl
julia>p=JuliASM.plot_heatmap(x)
```
"""
function plot_heatmap(x::Array{Array{Int64,1},1},eta::Array{Float64,1}; maxVal=4.0)
    # Ranges
    par1=-maxVal:0.1:maxVal
    par2=-maxVal:1:maxVal

    # Create (minus) log-likelihood function
    Loglike = JuliASM.create_Llkhd(x)

    # Matrix of (minus) log-likelihood values
    minY = Inf
    min_ij = [NaN,NaN]
    Y = zeros(Float64, length(par1), length(par1))
    for i in 1:length(par1), j in 1:length(par1)
        Y[i,j] = Loglike([par1[i],par1[j]])
        if Y[i,j]<minY
            minY=Y[i,j]
            min_ij=[i,j]
        end
    end

    # Plot: check keys(PlotUtils._gradients) for
    p = heatmap(Y, fillalpha=0.25, xlabel=L"\beta", ylabel=L"\alpha",
                title=latexstring("-\\log P[x_{obs};\\eta];~\\eta_0=$eta"),
                xticks=(1:10:length(par1),par2),
                yticks=(1:10:length(par1),par2),
                fillcolor=:thermal)
    annotate!(min_ij[2], min_ij[1], text(".", 25, :red, :center))

    # Return plot object
    return p
end # end plot_heatmap
