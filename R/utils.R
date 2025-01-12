.xyz_vec2mat = \(z,coords){
  xyz_m = cbind(coords,z)
  names(xyz_m) = c('x','y','z')
  xyz_m = tibble::as_tibble(xyz_m)
  return(sdsfun::tbl_xyz2mat(xyz_m))
}
