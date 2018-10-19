#' Prepare a data set with product data for estimation of a simple logit model
#' without random coefficients as initially described in Berry (1994)
berry_logit_data = function(productData, lhs="log_s_log_s0",s=NULL, s0=NULL, market_identifier=NULL) {
  restore.point("simple_linear_blp_data")

  if (! lhs %in% colnames(productData)) {
    if (is.null(s0)) s0 = "share_outside_product"

    if (!s0 %in% colnames(productData)) {
      productData = productData %>%
        group_by(!!(as.name(market_identifier))) %>%
        mutate(!!(as.name(s0)) := 1- sum(( !!(as.name(s)) ))) %>%
        ungroup()
    }
    productData[[lhs]] = log(productData[[s]]) - log(productData[[s0]])
  }
  productData
}


#' Compute elasticities for the simple logit model without
#' random coefficients initially described in Berry (1994)
#'
#' Let s_i be the market share of firm i, p_i the price,
#' and beta_b the coeficient of price in the utility function.
#'
#' The own price derivative is
#'
#' ds_i / dp_i = -beta_p * s_i * (1-s_i)
#'
#' The own price elasticity is
#'
#' ds_i / dp_i * (p_i / s_i) = -beta_p * (1- s_i) * p_i
#'
#' The cross price derivative with respect to the price of a good j is given by
#'
#' ds_i / dp_j = beta_p * s_i * s_j
#'
#' The cross price elasticity is
#'
#' ds_i / dp_j * (p_j / s_i) = beta_p * s_j * p_j
#' @export
berry_logit_get_elasticities = function(productData, estimation, variable, product_id, market_identifier,s=as.character(formula(estimation)[[2]]),  products=NULL, markets=NULL,simplify=TRUE) {
  restore.point("get_elasticities_simple_berry_logit")

  if (!is.null(products)) {
    productData = productData[productData[[product_id]] %in% products, ,drop=FALSE]
  }

  if (!is.null(markets)) {
    productData = productData[productData[[market_identifier]] %in% markets, ,drop=FALSE]
  } else {
    markets = unique(productData[[market_identifier]])
  }
  if (length(markets)>1) {
    li = lapply(markets, function(market) {
      get_elasticities_simple_berry_logit(productData, estimation, variable,product_id, market_identifier,s, markets=market)
    })
    return(li)
    stop("currently only works for a single market")
  }

  coef = coef(estimation)[[variable]]

  prods = productData[[product_id]]
  J = NROW(prods)

  shares = productData[[s]]
  xval = productData[[variable]]



  sj.mat = matrix(shares,J,J,byrow = TRUE)
  xj.mat = matrix(xval, J,J,byrow=TRUE)

  e = coef*sj.mat*xj.mat
  own = -coef * (1-shares) * xval
  diag(e) = own

  colnames(e) = rownames(e) = products


  if (!simplify) return(list(e))
  return(e)

}

#' Estimate a simple linear logit demand model without
#' random coefficients as initially described in Berry (1994)
#'
#' This can be used as comparison to the BLP estimate with random
#' coefficients. Alternative, one could derive an initial guess for
#' the delta of the BLP model
#'
#' @param blp_formula The formula for the BLP model
#' @param productData The product data of the BLP model
#' @param formula Optionally the direct formula to estimate the simple logit model via instrumental variables regression (using ivreg). If NULL this formula is created from blp_formula.
#' @param market_identifier The variable that identifies a market
#' @param s variable name that contains market shares of all products. Note that these markets shares should have been computed by dividing the sold number of products by the POTENTIAL market size. This means the shares should generally not add up to 1 in a market. The remaining share is the assumed market share of the 'buy-nothing' option. If NULL we use the LHS of blp_formula
#' @param s0 the variable that contains the share of customers that don't buy a product in the market. If NULL, we compute s0 automatically using the variable s.
#' @return The result of an ivreg regression.
#' @export
berry_logit_estimate = function(blp_formula, productData, formula=NULL,market_identifier=NULL, s=NULL, s0=NULL, lhs="log_s_log_s0") {
  restore.point("simple_linear_blp")

  if (is.null(s))
    s = as.character(blp_formula[[2]])

  dat = simple_berry_logit_data(productData,lhs,s,s0, market_identifier)

  if (is.null(formula)) {
    formula = blp_to_simple_formula(blp_formula, lhs=as.name(lhs))
  }

  library(AER)
  ivreg(formula, data=dat)

}

blp_to_berry_logit_formula = function(blp_formula, lhs=NULL) {
  F = Formula(blp_formula)

  if (!is.null(lhs)) {
    lhs = formula(F,lhs=1, rhs=0)[[2]]
  }
  ols = formula(F, lhs=0, rhs=c(1), collapse=TRUE)[[2]]
  exo = formula(F, lhs=0, rhs=2)[[2]]
  inst = formula(F, lhs=0, rhs=4)[[2]]

  form = as.formula(substitute(lhs ~ ols | exo + inst, list(lhs=lhs, ols=ols, exo=exo, inst=inst)))
  form
}
