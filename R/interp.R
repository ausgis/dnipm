dnipm_pd = \(z,coords,predicts,
             method = "lagrange"){
  coords = as.matrix(coords)
  predicts = as.matrix(predicts)
  if (method == "lagrange") {
    res = LagrangeInterp(predicts,coords,z)
  }
  return(res)
}
