# Data mesh:

mesh.nrow = 30
mesh.ncol = 30
lattice = inla.mesh.lattice(x=seq(1,mesh.ncol,length.out=mesh.ncol),y=seq(1,mesh.nrow,length.out=mesh.nrow))
data.mesh = inla.mesh.create(lattice=lattice,extend=list(n=5),boundary=lattice$segm)

# Flow mesh
flow.mesh = data.mesh

# Data

im1 = data.mesh$loc[,1] > 10 & data.mesh$loc[,1] < 20 & data.mesh$loc[,2] > 10 & data.mesh$loc[,2] < 20
im2 = data.mesh$loc[,1] > 11 & data.mesh$loc[,1] < 21 & data.mesh$loc[,2] > 10 & data.mesh$loc[,2] < 20 # right
im2 = data.mesh$loc[,1] > 9 & data.mesh$loc[,1] < 19 & data.mesh$loc[,2] > 10 & data.mesh$loc[,2] < 20 # left
im2 = data.mesh$loc[,1] > 10 & data.mesh$loc[,1] < 20 & data.mesh$loc[,2] > 9 & data.mesh$loc[,2] < 19 # down

plot(data.mesh, rgl = TRUE, col = im1)
plot(data.mesh, rgl = TRUE, col = im2)

# Bend space
data.mesh$loc[,2] = (data.mesh$loc[,2]-15) * cos(2*(data.mesh$loc[,1]-15)/30)
flow.mesh = data.mesh

# Run flow estimation

fl = flow2(im1, im2,
           lik.fix = TRUE,
           data.mesh,
           flow.mesh,
           compressible = FALSE,
           u.spde.args = "default",
           v.spde.args = "default")

plot.flow(fl)
plot(flow.mesh, rgl = TRUE, col = fl$gx)
plot(flow.mesh, rgl = TRUE, col = fl$gy)
plot(flow.mesh, rgl = TRUE, col = fl$gt)
plot(flow.mesh, rgl = TRUE, col = fl$u)
plot(flow.mesh, rgl = TRUE, col = fl$v)
