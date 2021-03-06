% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blp_tools.R
\name{blp_delta_guess}
\alias{blp_delta_guess}
\title{Estimate the simple logit model without random coefficients to get an initial
guess of the representative utilities delta that can be passed to BLP_data.}
\usage{
blp_delta_guess(blp_formula, productData, market_identifier = NULL,
  s = NULL, s0 = NULL)
}
\arguments{
\item{blp_formula}{The formula for the BLP model}

\item{productData}{The product data of the BLP model}

\item{market_identifier}{The variable that identifies a market}

\item{s}{variable name that contains market shares of all products. Note that these markets shares should have been computed by dividing the sold number of products by the POTENTIAL market size. This means the shares should generally not add up to 1 in a market. The remaining share is the assumed market share of the 'buy-nothing' option. If NULL we use the LHS of blp_formula}

\item{s0}{the variable that contains the share of customers that don't buy a product in the market. If NULL, we compute s0 automatically using the variable s.}
}
\value{
A vector of delta values with as many rows as prductData
}
\description{
Estimate the simple logit model without random coefficients to get an initial
guess of the representative utilities delta that can be passed to BLP_data.
}
