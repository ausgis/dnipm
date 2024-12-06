tbl2mat = \(.tbl,x = 1, y = 2, z = 3){
  .tbl = .tbl[,c(x,y,z)]
  img = terra::rast(.tbl)
  zs = terra::as.matrix(img,wide = TRUE)

  return()
}
