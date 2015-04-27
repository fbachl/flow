require(png)
require(INLA)
require(pracma)

# # gaussian blob
# i1 = readPNG("data/frame1.png")
# i1 = cbind(i1,i1)
# i2 = array(0,dim=c(30,60))
# #i2[,2:30] = i1[,1:29]
# #i2[1:29,] = i1[2:30,]
# i2[2:30,1:59] = i1[1:29,2:60]

# Moving rectangle
i1 = array(0,dim=c(30,30))
i1[10:20,10:20] = 1

i2 = array(0,dim=c(30,30)); i2[,1:29] = i1[,2:30] # movement to the left
i2 = array(0,dim=c(30,30)); i2[,2:30] = i1[,1:29] # movement to the right
i2 = array(0,dim=c(30,30)); i2[2:30,] = i1[1:29,] # movement down
i2 = array(0,dim=c(30,30)); i2[1:29,] = i1[2:30,] # movement up
i2 = array(0,dim=c(30,30)); i2[1:29,1:29] = i1[2:30,2:30] # movement up left

mimage(i1)
mimage(i2)

# Incompressible flow
fl = flow(i1,i2,lik.fix = 1, compressible = FALSE)
plot(fl,scale=0.5)

# Compressible flow
cfl = flow(i1,i2,lik.fix = 1, compressible = TRUE)
plot(cfl,scale = 0.5)

# Magnitude
mean(sqrt(fl$u^2+fl$v^2))
mean(sqrt(cfl$u^2+cfl$v^2))

# Gradient components
mimage(fl$gx)
mimage(fl$gy)
mimage(fl$gt)

# Compare components
mimage(cbind(fl$u,fl$v))
image(fl$u) ; range(fl$u)
image(fl$v) ; range(fl$v)

# Plot standard deviation
field = "spdey"
property = "sd"
model = fl
colx = fl$result$summary.ran[[field]][[property]]
rgbx = rgb(normalize(colx),0,0)
plot(model$mesh,col=rgbx,rgl=TRUE)


