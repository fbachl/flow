plot.flow = function(fl,scale=1) {
  plot.inla.mesh(fl$data.mesh,edge.color=rgb(0.7,0.7,0.7))

  vx = cbind(fl$u,fl$u) * (fl$gx.loc - fl$g.loc)
  vy = cbind(fl$v,fl$v) * (fl$gy.loc - fl$g.loc)

  arrows(flow.mesh$loc[,1],
         flow.mesh$loc[,2],
         flow.mesh$loc[,1] + vx[,1] + vy[,1],
         flow.mesh$loc[,2] + vx[,2] + vy[,2] ,
         length = 0.1, col = "blue")

}

mimage = function(x){
  image(t(x[rev(1:dim(x)[1]),]))
}

whiten = function(x) {
  return( (x-min(x))/(max(x)-min(x)))
}
