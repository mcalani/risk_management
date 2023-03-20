using Gadfly, RDatasets
iris = dataset("datasets", "iris")
#plot(data::AbstractDataFrame, elements::Element...; mapping...)
p = plot(iris, x=:SepalLength, y=:SepalWidth, Geom.point);
img = SVG("iris_plot.svg", 14cm, 8cm)
draw(img, p)
plot(iris, x=:SepalLength, y=:SepalWidth)

function get_to_it(d)
  ppoint = plot(d, x=:SepalLength, y=:SepalWidth, Geom.point)
  pline = plot(d, x=:SepalLength, y=:SepalWidth, Geom.line)
  ppoint, pline
end
ps = get_to_it(iris)
map(display, ps)

gasoline = dataset("Ecdat", "Gasoline")
plot(gasoline, x=:Year, y=:LGasPCar, color=:Country, Geom.point, Geom.line)
