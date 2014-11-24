flow = function(i1,i2,i3){


  gx = i2#grad.x(i1,i2)
  gy = i3# grad.y(i1,i2)
  gt = i1#-grad.t(i1,i2)

  gxv = as.vector(gx)
  gyv = as.vector(gy)
  y = as.vector(gt)

  npix = length(gx)
  nrow = dim(gx)[1]
  ncol = dim(gx)[2]

  loc.x = matrix(seq(1, nrow, len=nrow), nrow, ncol)
  loc.y = matrix(seq(1, ncol, len=ncol), nrow, ncol, byrow=TRUE)


  # Define lattice and mesh
  lattice = inla.mesh.lattice(x=1:ncol,y=1:nrow)
  mesh = inla.mesh.create(lattice=lattice,extend=list(n=5),boundary=lattice$segm)

  #spde.args = list(alpha=2,prior.variance.nominal=10,theta.prior.prec=0.01)
  spde.args = list(alpha=2)
  spde.mdlx = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))
  spde.mdly = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))


  formula = y ~ f(spdex, model=spde.mdlx) + f(spdey, model=spde.mdly) -1

  A = inla.spde.make.A(mesh, loc=cbind(as.vector(loc.x),y=as.vector(loc.y)))
  Mx = A*matrix(gxv,npix,npix)
  My = A*matrix(gyv,npix,npix)


  #effects = list(y=y,spdex = as.vector(gx))
  stk = inla.stack(data=list(y=y),A=list(Mx,My),tag="fpp",
                   effects=list(spdex=1:spde.mdlx$n.spde,spdey=1:spde.mdly$n.spde))

  #stk = inla.stack(data=list(y=y),A=list(Mx),tag="fpp",
  #                 effects=list(spdex=1:spde.mdlx$n.spde))


  result = inla(formula,family ="gaussian",data=inla.stack.data(stk),
                control.predictor=list(A=inla.stack.A(stk)),verbose=TRUE)

  spdex.mean = result$summary.ran[["spdex"]][["mean"]]
  spdex.mean = matrix(spdex.mean[mesh$idx$lattice],nrow,ncol)

  spdey.mean = result$summary.ran[["spdey"]][["mean"]]
  spdey.mean = matrix(spdey.mean[mesh$idx$lattice],nrow,ncol)



  return(list(result=result,mesh=mesh,ux=spdex.mean,v=spdey.mean))
}

grad.x = function(i1,i2){
  gx = array(0,dim(i1))
  # gx[,1:dim(i1)[2]-1] = i1[,2:dim(i1)[2]]-i1[,1:dim(i1)[2]-1]
  gx[,1:dim(i1)[2]-1] = 0.5 * ( i1[,2:dim(i1)[2]]-i1[,1:dim(i1)[2]-1] + i2[,2:dim(i1)[2]]-i2[,1:dim(i1)[2]-1] )
  return(gx)
}

grad.y = function(i1,i2){
  gy = array(0,dim(i1))
  # gy[1:dim(i1)[2]-1,] = i1[2:dim(i1)[2],]-i1[1:dim(i1)[2]-1,]
  gy[1:dim(i1)[2]-1,] = 0.5 * ( i1[2:dim(i1)[2],]-i1[1:dim(i1)[2]-1,] + i2[2:dim(i1)[2],]-i2[1:dim(i1)[2]-1,] )
  return(gy)
}
grad.t = function(i1,i2) {i2-i1}
#
# ## a graph on a file
# cat("3 1 1 2 2 1 1 3 0\n", file="g.dat")
# g = inla.read.graph("g.dat")
# ## writing an inla.graph-object to file
# g.file = inla.write.graph(g, mode="binary")
# ## re-reading it from that file
# gg = inla.read.graph(g.file)

normalize = function(x) {
  return( (x-min(x))/(max(x)-min(x)))
}
