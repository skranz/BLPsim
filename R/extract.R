

blp_coef = function(blp, which=c("linear","rc")[1]) {
  restore.point("coef_blp")
  rc = which=="rc"
  theta = if (rc) blp$theta_rc else blp$theta_lin
  if (!rc) {
    names = rownames(theta)
    theta = as.vector(theta)
    names(theta) = names
  }

  theta
}

