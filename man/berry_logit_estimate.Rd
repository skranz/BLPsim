% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/berry_logit.R
\name{berry_logit_estimate}
\alias{berry_logit_estimate}
\title{Estimate a simple linear logit demand model without
random coefficients as initially described in Berry (1994)}
\usage{
berry_logit_estimate(blp_formula, productData, formula = NULL,
  market_identifier = NULL, s = NULL, s0 = NULL, lhs = "log_s_log_s0")
}
\arguments{
\item{blp_formula}{The formula for the BLP model}

\item{productData}{The product data of the BLP model}

\item{formula}{Optionally the direct formula to estimate the simple logit model via instrumental variables regression (using ivreg). If NULL this formula is created from blp_formula.}

\item{market_identifier}{The variable that identifies a market}

\item{s}{variable name that contains market shares of all products. Note that these markets shares should have been computed by dividing the sold number of products by the POTENTIAL market size. This means the shares should generally not add up to 1 in a market. The remaining share is the assumed market share of the 'buy-nothing' option. If NULL we use the LHS of blp_formula}

\item{s0}{the variable that contains the share of customers that don't buy a product in the market. If NULL, we compute s0 automatically using the variable s.}
}
\value{
The result of an ivreg regression.
}
\description{
This can be used as comparison to the BLP estimate with random
coefficients. Alternative, one could derive an initial guess for
the delta of the BLP model
}
