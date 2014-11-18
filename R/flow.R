flow = function(i1,i2){
  gx = grad.x(i1,i2)
  gy = grad.y(i1,i2)
  gt = -grad.t(i1,i2)

  gxv = as.vector(gx)
  gyv = as.vector(gy)
  y = as.vector(gt)

  npix = length(gx)
  nrow = dim(i1)[1]
  ncol = dim(i1)[2]

  loc.x = matrix(seq(1, ncol, len=ncol), nrow, ncol, byrow=TRUE)
  loc.y = matrix(seq(1, nrow, len=nrow), nrow, ncol)


  # Define lattice and mesh
  lattice = inla.mesh.lattice(x=1:ncol,y=1:nrow)
  mesh = inla.mesh.create(lattice=lattice,extend=list(n=5),boundary=lattice$segm)

  spde.args = list(alpha=2,prior.variance.nominal=10,theta.prior.prec=0.01)
  spde.mdlx = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))
  spde.mdly = do.call(inla.spde2.matern,c(list(mesh=mesh),spde.args))


  formula = y ~ f(spde, model=spde.mdlx) -1
  #+ f(gxy, model=spde.mdly) -1

  A = inla.spde.make.A(mesh, loc=as.matrix(x=as.vector(loc.x),y=as.vector(loc.y)))

  #effects = list(y=y,spdex = as.vector(gx))
  stk = inla.stack(data=list(y=y),A=list(A),tag="fpp",
                   effects=list(spde=1:spde.mdlx$n.spde))

  result = inla(formula,family ="gaussian",data=inla.stack.data(stk),
                control.predictor=list(A=inla.stack.A(stk)),verbose=TRUE)

  field = "spde"
  property = "mean"
  col = result$summary.ran[[field]][[property]]
  col = rgb(normalize(col),0,0)
  plot(mesh,col=col,rgl=TRUE)
  return(result)
}

grad.x = function(i1,i2){
  gx = array(0,dim(i1))
  gx[,1:dim(i1)[2]-1] = i1[,2:dim(i1)[2]]-i1[,1:dim(i1)[2]-1]
  return(gx)
}

grad.y = function(i1,i2){return(grad.x(t(i1),t(i2)))}

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
