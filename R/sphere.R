sphere.grad = function(data.mesh, flow.mesh, val1, val2, phi = 1) {

  #
  # Gradient along latitude direction
  #

  # convert to geographic coordinates
  eloc = flow.mesh$loc; colnames(eloc) = c("x","y","z")
  gloc = euc.to.geo(eloc, R = 1)

  # Geographic coordinates of gradient1 points
  g1.gloc = gloc
  g1.gloc$lat = g1.gloc$lat + phi
  msk = ( g1.gloc$lat > 90 ) & ( g1.gloc$lon > 0)
  g1.gloc[msk,"lon"] = g1.gloc$lon[msk] - 180
  msk = ( g1.gloc$lat > 90 ) & ( g1.gloc$lon < 0)
  g1.gloc[msk,"lon"] = g1.gloc$lon[msk] + 180

  # Euclidean coordinates of gradent1 points
  g1.eloc = geo.to.euc(g1.gloc,R=1)

  # Read values at gradient1 points and compute gradient
  A = inla.spde.make.A(data.mesh,loc = as.matrix(g1.eloc))
  B = inla.spde.make.A(data.mesh,loc = as.matrix(eloc))
  lgrad1 = as.vector(A%*%val1 - B%*%val1)
  lgrad2 = as.vector(A%*%val2 - B%*%val2)
  lgrad = 0.5 * ( lgrad1 + lgrad2 )


  #
  # Gradient orthogonal to latitude direction (gradient2)
  #

  cx = normalize(pracma::cross(as.matrix(eloc),as.matrix(g1.eloc)))
  g2.eloc = cos(pi*phi/180)*eloc + sin(pi*phi/180)*cx
  g2.gloc = euc.to.geo(g2.eloc, R=1)

  A = inla.spde.make.A(data.mesh,loc = as.matrix(g2.eloc))
  B = inla.spde.make.A(data.mesh,loc = as.matrix(eloc))

  ograd1 = as.vector(A%*%val1 - B%*%val1)
  ograd2 = as.vector(A%*%val2 - B%*%val2)
  ograd = 0.5 * ( ograd1 + ograd2 )

  if(any(!(abs(norm(g2.eloc)==1) - 1)<0.01)) {
    stop("WTF")
  }

  #
  # temporal gradient
  #

  A = inla.spde.make.A(data.mesh,loc = as.matrix(eloc))

  gradt = as.vector(A%*%val2 - A%*%val1)

  return(list(eloc = eloc,
              grad1 = lgrad,
              grad2 = ograd,
              eloc1 = g1.eloc,
              eloc2 = g2.eloc,
              gloc1 = g1.gloc,
              gloc2 = g2.gloc,
              gradt = gradt))
}




norm = function(x){sqrt(apply(x*x,MARGIN=1,sum))}
normalize = function(x) {t(t(x)/norm(x))}


euc.to.geo = function(euc,R=6371) {
  lat = 180*asin(euc[,"z"] / R)/pi
  lon = 180*atan2(euc[,"y"], euc[,"x"])/pi
  return(data.frame(lat=lat,lon=lon))
}

geo.to.euc = function(geo,R=6371) {
  euc = data.frame(x = R*cos(0.5*pi*geo$lat/90)*cos(0.5*pi*geo$lon/90),
                   y = R*cos(0.5*pi*geo$lat/90)*sin(0.5*pi*geo$lon/90),
                   z = R*sin(0.5*pi*geo$lat/90) )
  return(euc)
}

rot.x = function(eloc, phi) {
  M = matrix(c(1,0,0,
             0, cos(phi), -sin(phi),
             0, sin(phi), cos(phi)),nrow=3,byrow=TRUE)
  reloc = as.matrix(eloc) %*% M
  colnames(reloc) = colnames(eloc)
  return(reloc)
}

rot.y = function(eloc, phi) {
  M = matrix(c(cos(phi), 0, sin(phi),
               0, 1, 0,
               -sin(phi), 0, cos(phi)),nrow=3,byrow=TRUE)
  reloc = as.matrix(eloc) %*% M
  colnames(reloc) = colnames(eloc)
  return(reloc)
}

rot.z = function(eloc, phi) {
  M = matrix(c(cos(phi), -sin(phi), 0,
               sin(phi), cos(phi), 0,
               0, 0, 1),nrow=3,byrow=TRUE)
  reloc = as.matrix(eloc) %*% M
  colnames(reloc) = colnames(eloc)
  return(reloc)
}

