require(png)
require(INLA)


img = readPNG("data/frame1.png")
i1[1:10,] = i1[6:15,]
f1 = 1+0.5*matrix(runif(30*30),30,30)


i2 = array(0,dim=c(30,30))
i2[,2:30] = i1[,1:29]
i2[,6:30] = i1[,1:25]
f2 = 1+0.5*matrix(runif(30*30),30,30)

Y = i1*f1 + i2*f2

fl = flow(Y,f1,f2)

image(fl$u)
image(fl$v)

fl = flow(-grad.t(i1,i2),grad.x(i1,i2),grad.y(i1,i2))

image(fl$u)
image(fl$v)


field = "spdex"
property = "mean"
colx = fl$result$summary.ran[[field]][[property]]
rgbx = rgb(normalize(colx),0,0)
plot(fl$mesh,col=rgbx,rgl=TRUE)


hist(fl$result$summary.ran[["spdex"]][["mean"]])
hist(fl$result$summary.ran[["spdey"]][["mean"]])


i1 = matrix(0,14,34)
i2 = matrix(0,14,34)
i1[3:12,3:32]=matrix(c(1,0,0),10,30,byrow=TRUE)
i2[3:12,3:32]=matrix(c(0.7,0.3,0),10,30,byrow=TRUE)





i1 = readPNG("data/frame1.png")
i2 = matrix(0.1+runif(900),30,30)

fl = fit(i1,grad.x(i1,i2))


field = "spdex"
property = "mean"
colx = fl$result$summary.ran[[field]][[property]]
rgbx = rgb(normalize(colx),0,0)
plot(fl$mesh,col=rgbx,rgl=TRUE)
hist(colx)

property = "mean"
coly = fl$result$summary.ran[[field]][[property]]
rgby = rgb(normalize(coly),0,0)
plot(fl$mesh,col=rgby,rgl=TRUE)
