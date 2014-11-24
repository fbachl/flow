require(png)
require(INLA)

i1 = readPNG("data/frame1.png")
i2 = array(0,dim=c(30,30))
i2[,2:30] = i1[,1:29]

fl = flow(i1,i2)



image(fl$u)
image(fl$v)


field = "spdey"
property = "sd"
colx = fl$result$summary.ran[[field]][[property]]
rgbx = rgb(normalize(colx),0,0)
plot(fl$mesh,col=rgbx,rgl=TRUE)


