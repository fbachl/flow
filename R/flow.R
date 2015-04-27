flow = function(i1,i2,lik.fix=0,compressible=TRUE){


  gx = grad.x(i1,i2)
  gy = grad.y(i1,i2)
  gt = -grad.t(i1,i2)

  imv = as.vector((i1+i2)/2)
  gxv = as.vector(gx)
  gyv = as.vector(gy)
  y = as.vector(gt)

  npix = length(gx)
  nrow = dim(gx)[1]
  ncol = dim(gx)[2]

  loc.x = matrix(seq(1, nrow, len=nrow), nrow, ncol)
  loc.y = matrix(seq(1, ncol, len=ncol), nrow, ncol, byrow=TRUE)



  # Define lattice and mesh
  #lattice = inla.mesh.lattice(x=1:ncol,y=1:nrow)
  mesh.nrow = ncol
  mesh.ncol = nrow
  lattice = inla.mesh.lattice(x=seq(1,ncol,length.out=mesh.ncol),y=seq(1,nrow,length.out=mesh.nrow))
  mesh = inla.mesh.create(lattice=lattice,extend=list(n=5),boundary=lattice$segm)

  sigma0 = 1
  kappa0 = 1e-3
  tau0 = 1/(4*pi*kappa0^2*sigma0^2)^0.5

  spde.args = list(alpha=2, constr=FALSE,
                   B.tau = cbind(log(tau0),0),
                   B.kappa = cbind(log(kappa0),1))
  #spde.args = list(alpha=2)
  spde.mdlx = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))
  spde.mdly = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))


  formula = y ~ f(spdex, model=spde.mdlx) + f(spdey, model=spde.mdly) -1

  A = inla.spde.make.A(mesh, loc=cbind(as.vector(loc.x),y=as.vector(loc.y)))
  Mx = A*matrix(gxv,npix,mesh$n)
  My = A*matrix(gyv,npix,mesh$n)

  if (compressible) {
    # Locations used to construct divergence
    delta = 1

    ln.x = loc.x - delta
    ln.y = loc.y
    ln.A = inla.spde.make.A(mesh, loc=cbind(as.vector(ln.x),y=as.vector(ln.y)))


    rn.x = loc.x + delta
    rn.y = loc.y
    rn.A = inla.spde.make.A(mesh, loc=cbind(as.vector(rn.x),y=as.vector(rn.y)))

    tn.x = loc.x
    tn.y = loc.y + delta
    tn.A = inla.spde.make.A(mesh, loc=cbind(as.vector(tn.x),y=as.vector(tn.y)))

    bn.x = loc.x
    bn.y = loc.y - delta
    bn.A = inla.spde.make.A(mesh, loc=cbind(as.vector(bn.x),y=as.vector(bn.y)))

    Mx = Mx + 0.5 * ( ln.A * matrix(imv,npix,mesh$n) - rn.A * matrix(imv,npix,mesh$n) )
    My = My + 0.5 * ( bn.A * matrix(imv,npix,mesh$n) - tn.A * matrix(imv,npix,mesh$n) )

    warning("This might not be correct. imv should maybe replaced by image evaluated at the nb.x etc.")
    warning("We also did not check yet if all the directions are correct!")
  }


  #effects = list(y=y,spdex = as.vector(gx))
  stk = inla.stack(data=list(y=y),A=list(Mx,My),tag="fpp",
                   effects=list(spdex=1:spde.mdlx$n.spde,spdey=1:spde.mdly$n.spde))

  #stk = inla.stack(data=list(y=y),A=list(Mx),tag="fpp",
  #                 effects=list(spdex=1:spde.mdlx$n.spde))

  lik.hyper = list(prec = list(prior = "loggamma",
                               param = c(1.0,0.01),
                               initial = 1e1,
                               fixed=lik.fix))

  result = inla(formula,family ="gaussian",data=inla.stack.data(stk),
                control.family = list(hyper=lik.hyper),
                control.predictor=list(A=inla.stack.A(stk)),verbose=TRUE)

  spdex.mean = result$summary.ran[["spdex"]][["mean"]]
  spdex.mean = matrix(spdex.mean[mesh$idx$lattice],mesh.ncol,mesh.nrow)

  spdey.mean = result$summary.ran[["spdey"]][["mean"]]
  spdey.mean = matrix(spdey.mean[mesh$idx$lattice],mesh.ncol,mesh.nrow)

  result = list(result=result,mesh=mesh,u=spdex.mean,v=spdey.mean,gx=gx,gy=gy,gt=gt)
  class(result) = c("flow","list")

  return(result)
}


grad.x = function(i1,i2){
  gx = array(0,dim(i1))
  gx[1:dim(i1)[1]-1,] = 0.5 * ( i1[2:dim(i1)[1],]-i1[1:(dim(i1)[1]-1),] + i2[2:dim(i1)[1],]-i2[1:(dim(i1)[1]-1),] )
  return(gx)
}

grad.y = function(i1,i2){
  gy = array(0,dim(i1))
  gy[,1:dim(i1)[2]-1] = 0.5 * ( i1[,2:dim(i1)[2]]-i1[,1:(dim(i1)[2]-1)] + i2[,2:dim(i1)[2]]-i2[,1:(dim(i1)[2]-1)] )
  return(gy)
}

grad.t = function(i1,i2) {i2-i1}

normalize = function(x) {
  return( (x-min(x))/(max(x)-min(x)))
}
