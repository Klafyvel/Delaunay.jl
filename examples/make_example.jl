using Plots
using DelaunayMaps

m  = Map(20, 20, 3);
plot_map(m, window=:fit)
savefig("../images/example_3_points.png")
m  = Map(20, 20, 10);
plot_map(m, show_nums=false, window=:fit)
savefig("../images/example_10_points.png")
m  = Map(20, 20, 100);
plot_map(m, show_nums=false, window=:fit)
savefig("../images/example.png")
