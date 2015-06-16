plot.flow = function(fl, col, scale=1, manifold = fl$flow.mesh$manifold) {

if ( manifold == "R2" ){
  plot.inla.mesh(fl$data.mesh,edge.color=rgb(0.7,0.7,0.7))

  vx = cbind(fl$u+fl$ub,fl$u+fl$ub) * (fl$gx.loc - fl$g.loc)
  vy = cbind(fl$v+fl$vb,fl$v+fl$vb) * (fl$gy.loc - fl$g.loc)

  arrows(flow.mesh$loc[,1],
         flow.mesh$loc[,2],
         flow.mesh$loc[,1] + vx[,1] + vy[,1],
         flow.mesh$loc[,2] + vx[,2] + vy[,2] ,
         length = 0.1, col = "blue")
  }
  else if ( manifold == "S2"){
    u = fl$u + fl$ub
    v = fl$v + fl$vb

    line.matrix = function(from,to){
      dfr = data.frame(from,to, NA,NA,NA,NA,NA,NA)
      return( matrix(t(as.matrix(dfr)),ncol=3,byrow=TRUE) )
    }

    rgl.sphere(fl$data.mesh,col)
#     rgl.points(fl$grad$eloc,col="black")
#     rgl.points(fl$grad$eloc1,col="red")
#     rgl.points(fl$grad$eloc2,col="blue")

    scale = 0.5 * 1/max(max(abs(u)),max(abs(v)))

#     rgl.lines(1.01*line.matrix(fl$grad$eloc,fl$grad$eloc1), col = "red", lwd = 2)
#     rgl.lines(1.01*line.matrix(fl$grad$eloc,fl$grad$eloc2), col = "green", lwd = 2)

    # flow arrows
    arrows = line.matrix(fl$grad$eloc,fl$grad$eloc + scale * ((u) * (fl$grad$eloc1 - fl$grad$eloc) + (v) * (fl$grad$eloc2 - fl$grad$eloc)))
    rgl.lines(1.01*arrows, col = "blue", lwd = 2)
  }
}

rgl.sphere = function(mesh, col){

  plot(mesh,rgl=TRUE,col)
  rgl.lines(rbind(data.frame(x=0,y=0,z=0),1.2*data.frame(x=1,y=0,z=0)), col = "red", lwd = 4)
  rgl.lines(rbind(data.frame(x=0,y=0,z=0),1.2*data.frame(x=0,y=1,z=0)), col = "green", lwd = 4)
  rgl.lines(rbind(data.frame(x=0,y=0,z=0),1.2*data.frame(x=0,y=0,z=1)), col = "blue", lwd = 4)

  rgl.viewpoint(fov = 120)

}



mimage = function(x){
  image(t(x[rev(1:dim(x)[1]),]))
}

whiten = function(x) {
  return( (x-min(x))/(max(x)-min(x)))
}
