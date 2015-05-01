flow = function(i1, i2,
                 lik.fix = TRUE,
                 data.mesh,
                 flow.mesh,
                 compressible=TRUE,
                 u.spde.args = "default",
                 v.spde.args = "default",
                 delta = 1){
  #
  # Step 1: Determine locations for estimation of gradient
  #
  if (flow.mesh$manifold == "R2"){
    g.loc = flow.mesh$loc[,1:2]
    gx.loc = cbind(flow.mesh$loc[,1] + delta, flow.mesh$loc[,2])
    gy.loc = cbind(flow.mesh$loc[,1], flow.mesh$loc[,2] + delta)
    gt.loc = flow.mesh$loc

    Agx = inla.spde.make.A(data.mesh, gx.loc)
    Agy = inla.spde.make.A(data.mesh, gy.loc)
    Ag = inla.spde.make.A(data.mesh, g.loc)

    gx = as.vector( 0.5 * ( ( Agx %*% i1 - Ag %*% i1 ) + ( Agx %*% i2 - Ag %*% i2 ) ) )
    gy = as.vector( 0.5 * ( ( Agy %*% i1 - Ag %*% i1 ) + ( Agy %*% i2 - Ag %*% i2 ) ) )
    gt = as.vector( Ag %*% i2 - Ag %*% i1 )

    grad = NULL

  } else if ( flow.mesh$manifold == "S2") {


    grad = sphere.grad(data.mesh, flow.mesh, i1, i2, phi = delta)

    plot(flow.mesh, rgl = TRUE, col = grad$gradt)
    plot(flow.mesh, rgl = TRUE, col = grad$grad1)
    plot(flow.mesh, rgl = TRUE, col = grad$grad2)

    if(any(!(abs(norm(grad$eloc2)==1) - 1)<0.01)) {
      stop("WTF")
    }

    g.loc = grad$loc
    gx.loc = grad$g1.eloc
    gy.loc = grad$g2.eloc
    gt.loc = grad$loc

    gx = grad$grad1
    gy = grad$grad2
    gt = grad$gradt

  } else {
    stop("Not implemented")
  }

  #
  # Default SPDE parameterization
  #

  sigma0 = 1
  kappa0 = 1e-3
  tau0 = 1/(4*pi*kappa0^2*sigma0^2)^0.5

  if (is.character(u.spde.args) & ( u.spde.args == "default") ) {
    u.spde.args = list(alpha=2, constr=FALSE,
                       B.tau = cbind(log(tau0),0),
                       B.kappa = cbind(log(kappa0),1))
  }

  if (is.character(v.spde.args) & ( v.spde.args == "default") ) {
    v.spde.args = list(alpha=2, constr=FALSE,
                       B.tau = cbind(log(tau0),0),
                       B.kappa = cbind(log(kappa0),1))
  }

  spde.mdlx = do.call(inla.spde2.matern,c(list(mesh=flow.mesh),u.spde.args))
  spde.mdly = do.call(inla.spde2.matern,c(list(mesh=flow.mesh),v.spde.args))

  formula = y ~ f(spdex, model=spde.mdlx) + f(spdey, model=spde.mdly) -1

  A = inla.spde.make.A(flow.mesh, loc = flow.mesh$loc)

  Mx = A * matrix(gx, length(gx), flow.mesh$n)
  My = A * matrix(gy, length(gy), flow.mesh$n)

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


  stk = inla.stack(data = list(y = -gt),
                   A = list(Mx,My),
                   tag = "flow.tag",
                   effects = list(spdex = 1:spde.mdlx$n.spde,
                                  spdey = 1:spde.mdly$n.spde)
  )

  lik.hyper = list(prec = list(prior = "loggamma",
                               param = c(1.0,0.01),
                               initial = 1e1,
                               fixed=lik.fix))

  result = inla(formula,family ="gaussian",data=inla.stack.data(stk),
                control.family = list(hyper=lik.hyper),
                control.predictor=list(A=inla.stack.A(stk)),verbose=TRUE)

  spdex.mean = result$summary.ran[["spdex"]][["mean"]]
  spdey.mean = result$summary.ran[["spdey"]][["mean"]]

  result = list(result = result,
                data.mesh = data.mesh,
                flow.mesh = flow.mesh,
                u=spdex.mean,
                v=spdey.mean,
                gx = gx,
                gy = gy,
                gt = gt,
                g.loc = g.loc,
                gx.loc = gx.loc,
                gy.loc = gy.loc,
                gt.loc = gt.loc,
                grad = grad)

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

