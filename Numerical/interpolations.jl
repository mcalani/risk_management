using Interpolations
using Plots

u=rand(10)
x=collect(1:1:10)

plot(u)
intp=LinearInterpolation(x,u)

u[9]
intp(9) #funcion continua

Interp = interpolate((collect(grid),), v0, Gridded(Linear()) )
Interp = extrapolate(Interp,Interpolations.Flat())

int_f= interpolate((x,),u, BSpline(Linear()) ,OnGrid())
