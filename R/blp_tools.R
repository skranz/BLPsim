
#' Similar to get_elasticities, but gets corresponding derivatives
#'
#' E.g. if the variable is the price, we compute the matrix
#' of cross-price derivatives for the selected products for
#' a particular market.
#'
#' This matrix is crucial to compute Bertrand equilibria.
#' @export
blp_get_derivatives = function(blp_data, blp_estimation, variable, products=NULL, market, e=NULL) {
  restore.point("get_blp_derivatives")

  if (is.null(e)) {
    e = get_elasticities(blp_data, blp_estimation, variable , products, market)
  }


  # We need to extract the relevant rows from blp_data
  par = blp_data$parameters
  if (!is.null(products)) {
    rows = which(par$market_id_char_in == market & par$product_id %in% products)
  } else {
    rows = which(par$market_id_char_in == market)
  }

  if (is.null(products)) {
    products = rownames(e)
  }


  # We need the values of the market share and the variable
  # for the given rows
  s = blp_data$data$shares[rows]

  # We assume the variable is in Xlin
  x = blp_data$data$X_lin[rows, variable]

  # We need to reorder according to products
  prods = par$product_id[rows]
  order = match(products,prods)
  #cbind(prods[order], products)
  s = s[order]
  x = x[order]

  s.mat = matrix(s,NROW(e),NCOL(e), byrow=FALSE)
  x.mat = matrix(x,NROW(e),NCOL(e), byrow=TRUE)

  deriv = e * s.mat / x.mat
  deriv
}

#' Estimate the simple logit model without random coefficients to get an initial
#' guess of the representative utilities delta that can be passed to BLP_data.
#'
#' @param blp_formula The formula for the BLP model
#' @param productData The product data of the BLP model
#' @param market_identifier The variable that identifies a market
#' @param s variable name that contains market shares of all products. Note that these markets shares should have been computed by dividing the sold number of products by the POTENTIAL market size. This means the shares should generally not add up to 1 in a market. The remaining share is the assumed market share of the 'buy-nothing' option. If NULL we use the LHS of blp_formula
#' @param s0 the variable that contains the share of customers that don't buy a product in the market. If NULL, we compute s0 automatically using the variable s.
#'
#' @return A vector of delta values with as many rows as prductData
#'
#' @export
blp_delta_guess = function(blp_formula, productData,market_identifier=NULL,s=NULL, s0=NULL) {
  iv = estimate_simple_berry(blp_formula,productData,s=s,s0=s0, market_identifier=market_identifier)
  fitted(iv)
}

#' Helper function to generate a formula for
blp_make_formula = function(dep, lin=NULL, exo=NULL, rc=NULL, inst=NULL) {
  restore.point("blp.make.formula")

  if (is.null(lin)) lin = "0"
  if (is.null(exo)) exo = "0"
  if (is.null(rc)) rc = "0"
  if (is.null(inst)) inst = "0"


  str = paste0(dep ,"~ ", paste0(lin, collapse="+"),
    " | ", paste0(exo, collapse="+"),
    " | ", paste0(rc, collapse="+"),
    " | ", paste0(inst, collapse="+")
  )
  as.formula(str)
}

#' Compute BLP instruments and add to product data
#'
#' BLP suggest to instrument prices with the values of other product
#' characteristics.
#'
#' If the characteristics of competing products change between markets,
#' one may react with the price for the own product. For example, if ceteris
#' paribus other products get more attractive attributes, it can be optimal
#' to reduce the prices of the own product to mitigate market share losses.
#'
#' The validity of these instrument depends on the stronger assumption that
#' other firm's observed attributes are not correlated with the own product's
#' unobserved attributes.
#'
#' Consider a BLP model for the demand of cars and the attribute horse power.
#' Concrentely BLP suggest to compute to instrumental variables
#' related to this attribute:
#'
#' 1. The sum of horse power of all other cars from the same firm
#'
#' 2. The sum of horse power of all cars from a different firms
#'
#' This function allows to quickly compute such instrumental variables
#'
#' @param productData The data frame that contains the product data
#' @param vars The variable names of all the product attributes for which instruments shall be computed
#' @param market_identifier The variable name that identifies a market
#' @param firm_var The variable that identifies a firm. If NULL only a single instrumental variable per attribute is computed that consists of the sum of the attribute for all other products.
#' @param prefix a prefix to add to the new instrumental variables. If NULL a default will be chosen. If firm_var is not NULL prefix must be a vector of two strings. The first is the prefix for the intra-firm sums and the second for the inter-firm sums.
#' @param add_cols if TRUE the new instrumental variables will be added as new columns to product data. Otherwise only a data frame with the new instruments will be returned.
#'
#' @return A data frame that contains the new instruments as columns. The data frame also has an attribute "instruments" that contains the names of all generated instrumental variables.

blp_create_instruments = function(productData, vars, market_identifier, firm_var=NULL, prefix = NULL, add_cols=TRUE) {
  restore.point("blp_create_instruments")

  verb = if (add_cols) "mutate" else "transmute"
  if (is.null(firm_var)) {
    if (is.null(prefix)) prefix = "other_sum_"

    dat = group_by_at(productData,market_identifier)

    new_vars = paste0(prefix,vars)
    code = paste0(verb,"(dat,",
      paste0(new_vars," = sum(",vars,") - ", vars, collapse=","),
    ")")
    res = eval(parse(text=code))
  } else {
    if (is.null(prefix)) {
      prefix = c(
        paste0("other_sum_in_",firm_var,"_"),
        paste0("other_sum_outside_",firm_var,"_")
      )
    }

    # Intra group other sums
    dat = group_by_at(productData,c(market_identifier,firm_var))
    new_vars1 = paste0(prefix[1],vars)
    code = paste0(verb,"(dat,",
      paste0(new_vars1," = sum(",vars,") - ", vars, collapse=","),
    ")")
    res = eval(parse(text=code))

    # Inter group other sums
    res = group_by_at(res,market_identifier)
    new_vars2 = paste0(prefix[2],vars)
    code = paste0("mutate(res,",
      paste0(new_vars2," = sum(",new_vars1,") - ", new_vars1, collapse=","),
    ")")
    res = eval(parse(text=code))
    new_vars = c(new_vars1, new_vars2)
  }

  if (!add_col) {
    res = res[, new_vars]
  }
  attr(res,"instruments") = union(attr(res,"instruments"),new_vars)

  return(res)

}

#' Compute the sum of all other variables
other_sum = function(x) {
  sum(x)-x
}
