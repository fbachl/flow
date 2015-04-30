plot.flow = function(fl,scale=1) {
  grid = expand.grid(1:dim(fl$u)[1],1:dim(fl$u)[2])
  plot(x=NULL,y=NULL,xlim=c(1,dim(fl$u)[1]),ylim=c(1,dim(fl$u)[1]))

  len = 0.7*whiten(as.vector(sqrt(fl$u^2+fl$v^2)))

  invisible(mapply(arrows,
         as.vector(grid$Var1),
         as.vector(grid$Var2),
         grid$Var1+scale*as.vector(fl$v),
         grid$Var2+scale*as.vector(-fl$u),
         length=0.1*len))
}

mimage = function(x){
  image(t(x[rev(1:dim(x)[1]),]))
}

whiten = function(x) {
  return( (x-min(x))/(max(x)-min(x)))
}
