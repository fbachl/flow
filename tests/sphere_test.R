require(INLA)
require(pracma)

# Data mesh:
coords = expand.grid(lon=seq(-180,180,by=10),lat=seq(-90,90,by=10))
lattice = inla.mesh.lattice(x=coords$lon,y=coords$lat,units="longlat")
data.mesh = inla.mesh.create(lattice=lattice)
plot(data.mesh)

# Flow mesh:
coords = expand.grid(lon=seq(-180,180,by=10),lat=seq(-80,80,by=10))
lattice = inla.mesh.lattice(x=coords$lon,y=coords$lat,units="longlat")
flow.mesh = inla.mesh.create(lattice=lattice)
plot(flow.mesh)


# Generate data (a square)

  eloc = data.mesh$loc; colnames(eloc) = c("x","y","z")
  gloc = euc.to.geo(eloc, R = 1)

  im1 = 1 * (gloc[,"lat"] < 20 & gloc[,"lat"] > -20 & gloc[,"lon"] < 20 & gloc[,"lon"] > -20)
  im2 =  1 * (gloc[,"lat"] < 20 & gloc[,"lat"] > -20 & gloc[,"lon"] < 25 & gloc[,"lon"] > -10) # right
  im2 =  1 * (gloc[,"lat"] < 20 & gloc[,"lat"] > -20 & gloc[,"lon"] < 10 & gloc[,"lon"] > -25) # left
  im2 =  1 * (gloc[,"lat"] < 25 & gloc[,"lat"] > -10 & gloc[,"lon"] < 20 & gloc[,"lon"] > -20) # up
  im2 =  1 * (gloc[,"lat"] < 5 & gloc[,"lat"] > -25 & gloc[,"lon"] < 20 & gloc[,"lon"] > -20) # down

  rgl.sphere(data.mesh, im1)
  rgl.sphere(data.mesh, im2)

# Compute flow

  fl = flow(im1, im2,
             lik.fix = TRUE,
             data.mesh = data.mesh,
             flow.mesh = flow.mesh,
             compressible = FALSE,
             u.spde.args = "default",
             v.spde.args = "default",
             delta = 10)

# Plot flow

  plot(fl,im1)
  plot(fl,im2)
