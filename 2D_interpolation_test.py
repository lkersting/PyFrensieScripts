#! /usr/bin/env python
import PyFrensie.Data.Native as Native
import PyFrensie.Utility as Utility
import PyFrensie.Utility.Prng as Prng
import PyFrensie.Utility.Distribution as Distribution
import PyFrensie.MonteCarlo.Collision as Collision
import PyTrilinos.Teuchos as Teuchos
import numpy
import matplotlib.pyplot as plt
from collections import Counter


def lin_interp( alpha, y_0, y_1):
  return y_0 + (y_1 - y_0)*alpha

def log_interp( alpha, y_0, y_1):
  return y_0*pow(y_1/y_0,alpha)

def unit_base_pdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = e_loss_0_max - e_loss_0_min

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = e_loss_1_max - e_loss_1_min

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = e_loss_max - e_loss_min

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper pdf values
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )

  # Calculate the pdf
  return lin_interp( E_alpha, pdf_0*e_loss_0_L, pdf_1*e_loss_1_L )/e_loss_L

def unit_base_pdf_log( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
    # Get lower min, max and length
    e_loss_0_min = e_losses_0[0]
    e_loss_0_max = e_losses_0[len(e_losses_0)-1]
    e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

    # Get upper min, max and length
    e_loss_1_min = e_losses_1[0]
    e_loss_1_max = e_losses_1[len(e_losses_1)-1]
    e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

    # Get E' min, max and length
    e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
    e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
    e_loss_L = numpy.log(e_loss_max/e_loss_min)

    # Calculate eta for E'
    eta = numpy.log( e_loss/e_loss_min )/e_loss_L

    # Calculate the lower and upper E' values
    if eta == 0.0:
        e_loss_0 = e_loss_0_min
        e_loss_1 = e_loss_1_min
    elif eta >= 1.0:
        e_loss_0 = e_loss_0_max
        e_loss_1 = e_loss_1_max
    else:
        e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
        e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

    # Calculate the lower and upper pdf values
    pdf_0 = dist_0.evaluatePDF( e_loss_0 )
    pdf_1 = dist_1.evaluatePDF( e_loss_1 )

    # Calculate the pdf
    return log_interp( E_alpha, pdf_0*e_loss_0_L, pdf_1*e_loss_1_L )/e_loss_L

def unit_base_cdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = (e_loss_0_max - e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = (e_loss_1_max - e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = (e_loss_max - e_loss_min)

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = lin_interp( E_alpha, cdf_0, cdf_1 )

  return cdf

def unit_base_cdf_log( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = numpy.log(e_loss_max/e_loss_min)

  # Calculate eta for E'
  eta = numpy.log( e_loss/e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
    e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = log_interp( E_alpha, cdf_0, cdf_1 )

  return cdf

def unit_base2_pdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = e_loss_0_max - e_loss_0_min

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = e_loss_1_max - e_loss_1_min

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = e_loss_max - e_loss_min

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper pdf values
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )

  # Calculate the pdf
  return lin_interp( E_alpha, pdf_0*e_loss_0_L, pdf_1*e_loss_1_L )/e_loss_L

def unit_base2_pdf_log( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
    # Get lower min, max and length
    e_loss_0_min = e_losses_0[0]
    e_loss_0_max = e_losses_0[len(e_losses_0)-1]
    e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

    # Get upper min, max and length
    e_loss_1_min = e_losses_1[0]
    e_loss_1_max = e_losses_1[len(e_losses_1)-1]
    e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

    # Get E' min, max and length
    e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
    e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
    e_loss_L = numpy.log(e_loss_max/e_loss_min)

    # Calculate eta for E'
    eta = numpy.log( e_loss/e_loss_min )/e_loss_L

    # Calculate the lower and upper E' values
    if eta == 0.0:
        e_loss_0 = e_loss_0_min
        e_loss_1 = e_loss_1_min
    elif eta >= 1.0:
        e_loss_0 = e_loss_0_max
        e_loss_1 = e_loss_1_max
    else:
        e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
        e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

    # Calculate the lower and upper pdf values
    pdf_0 = dist_0.evaluatePDF( e_loss_0 )
    pdf_1 = dist_1.evaluatePDF( e_loss_1 )

    # Calculate the pdf
    return log_interp( E_alpha, pdf_0*e_loss_0_L*e_loss_0_L, pdf_1*e_loss_1_L*e_loss_1_L )/(e_loss_L*e_loss_L)

def unit_base2_cdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = (e_loss_0_max - e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = (e_loss_1_max - e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = (e_loss_max - e_loss_min)

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = lin_interp( E_alpha, cdf_0, cdf_1 )

  return cdf

def unit_base2_cdf_log( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = numpy.log(e_loss_max/e_loss_min)

  # Calculate eta for E'
  eta = numpy.log( e_loss/e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
    e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = log_interp( E_alpha, cdf_0, cdf_1 )

  return cdf

def correlated_unit_base_pdf_lin( cdf, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower min, max, length, sample and pdf
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = e_loss_0_max - e_loss_0_min
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )
  product_0 = pdf_0*e_loss_0_L

  # Get upper min, max, length, sample and pdf
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = e_loss_1_max - e_loss_1_min
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )
  product_1 = pdf_1*e_loss_1_L

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = e_loss_max - e_loss_min

  numerator = product_0*product_1
  if numerator == 0.0:
    return 0.0
  else:
    denominator = lin_interp( E_alpha, product_1, product_0)*e_loss_L
    return numerator/denominator

def correlated_unit_base_pdf_log( cdf, e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
    if cdf > 1.0:
        cdf = 1.0 - 1e-15
    if cdf == 0.0:
      return correlated_unit_base_pdf_lin( cdf, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha)

    # Get lower min, max, length, sample and pdf
    e_loss_0_min = e_losses_0[0]
    e_loss_0_max = e_losses_0[len(e_losses_0)-1]
    e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)
    e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
    pdf_0 = dist_0.evaluatePDF( e_loss_0 )
    product_0 = (e_loss_0 - e_loss_0_min)*pdf_0

    # Get upper min, max, length, sample and pdf
    e_loss_1_min = e_losses_1[0]
    e_loss_1_max = e_losses_1[len(e_losses_1)-1]
    e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)
    e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
    pdf_1 = dist_1.evaluatePDF( e_loss_1 )
    product_1 = (e_loss_1 - e_loss_1_min)*pdf_1

    # Get E' min, max and length
    e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
    e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
    e_loss_L = numpy.log(e_loss_max/e_loss_min)

    numerator = product_0*product_1
    denominator = (e_loss-e_loss_min)*lin_interp( E_alpha, product_1, product_0 )
    print pdf_0, pdf_1, numerator/denominator
    return numerator/denominator

def sample_correlated_unit_base_lin( cdf, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower min, max, length, sample and eta
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = e_loss_0_max - e_loss_0_min
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
  eta_0 = ( e_loss_0 - e_loss_0_min )/e_loss_0_L

  # Get upper min, max, length, sample and eta
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = (e_loss_1_max - e_loss_1_min)
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
  eta_1 = ( e_loss_1 - e_loss_1_min )/e_loss_1_L

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = (e_loss_max - e_loss_min)

  # Calculate eta for E'
  if eta_0 == eta_1:
    eta = eta_0
  else:
    eta = lin_interp( E_alpha, eta_0, eta_1 )

  # Calculate the sampled value
  sample = e_loss_min + e_loss_L*eta

  return sample

def sample_correlated_unit_base_log( cdf, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower min, max, length, sample and eta
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
  eta_0 = numpy.log( e_loss_0/e_loss_0_min )/e_loss_0_L

  # Get upper min, max, length, sample and eta
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
  eta_1 = numpy.log( e_loss_1/e_loss_1_min )/e_loss_1_L

  # Get E' min, max and length
  e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = numpy.log(e_loss_max/e_loss_min)

  # Calculate eta for E'
  if eta_0 == eta_1:
    eta = eta_0
  else:
    eta = log_interp( E_alpha, eta_0, eta_1 )

  # Calculate the sampled value
  sample = numpy.exp( numpy.log(e_loss_min) + e_loss_L*eta )

  return sample

def correlated_direct_pdf_lin( cdf, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower sample and pdf
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )

  # Get upper sample and pdf
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )

  numerator = pdf_0*pdf_1
  if numerator == 0.0:
    return 0.0
  else:
    denominator = lin_interp( E_alpha, pdf_1, pdf_0)
    return numerator/denominator

def correlated_direct_pdf_log( cdf, e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
    if cdf > 1.0:
      cdf = 1.0 - 1e-15

    # Get lower sample and pdf
    e_loss_0 = dist_0.sampleWithRandomNumber( cdf )
    pdf_0 = dist_0.evaluatePDF( e_loss_0 )
    product_0 = e_loss_0*pdf_0

    # Get upper sample and pdf
    e_loss_1 = dist_1.sampleWithRandomNumber( cdf )
    pdf_1 = dist_1.evaluatePDF( e_loss_1 )
    product_1 = e_loss_1*pdf_1

    numerator = product_0*product_1
    if numerator == 0.0:
        return 0.0
    else:
        denominator = e_loss*lin_interp( E_alpha, product_1, product_0 )
        return numerator/denominator

def sample_direct_lin( cdf, dist_0, dist_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower sample
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )

  # Get upper sample
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )

  # Calculate the sampled value
  sample = lin_interp( E_alpha, e_loss_0, e_loss_1 )

  return sample

def sample_direct_log( cdf, dist_0, dist_1, E_alpha):
  if cdf > 1.0:
    cdf = 1.0 - 1e-15

  # Get lower sample
  e_loss_0 = dist_0.sampleWithRandomNumber( cdf )

  # Get upper sample
  e_loss_1 = dist_1.sampleWithRandomNumber( cdf )

  # Calculate the sampled value
  sample = log_interp( E_alpha, e_loss_0, e_loss_1 )

  return sample

def construct_processed_cdf_distribution( pdfs, cdfs, e_losses ):

  # Process the e_losses
  process_cdf_e_losses = numpy.zeros( shape=( len( e_losses ) ) )
  for i in range(0, len(e_losses)):
    process_cdf_e_losses[i] = numpy.log(e_losses[i])

  # construct the distribution
  process_cdf_dist = Distribution.TabularCDFDistribution_LogLin(process_cdf_e_losses, cdfs, True)

  # Get the processed cdf dist cdf values
  process_cdf_cdfs = numpy.zeros( shape=( len( cdfs ) ) )
  process_cdf_pdfs = numpy.zeros( shape=( len( cdfs ) ) )
  for i in range(0, len(process_cdf_cdfs)):
    process_cdf_cdfs[i] = process_cdf_dist.evaluateCDF( process_cdf_e_losses[i] )
    process_cdf_pdfs[i] = process_cdf_dist.evaluatePDF( process_cdf_e_losses[i] )

  return process_cdf_e_losses, process_cdf_pdfs, process_cdf_cdfs, process_cdf_dist

def recast_distributions( dist_0, dist_1, cdfs_0, cdfs_1 ):

  # Combine unique cdf values from dist_0 and dist_1
  combined_cdf = numpy.unique( numpy.concatenate((cdfs_0, cdfs_1),0))

  # Reevaluate the e_loss and pdf values on new grid
  combined_e_loss_0 = numpy.zeros( shape=( len( combined_cdf ) ) )
  combined_e_loss_1 = numpy.zeros( shape=( len( combined_cdf ) ) )
  combined_pdf_0 = numpy.zeros( shape=( len( combined_cdf ) ) )
  combined_pdf_1 = numpy.zeros( shape=( len( combined_cdf ) ) )
  for i in range(0, len(combined_cdf)):
    # Sample new energy loss values
    combined_e_loss_0[i] = dist_0.sampleWithRandomNumber( combined_cdf[i] )
    combined_e_loss_1[i] = dist_1.sampleWithRandomNumber( combined_cdf[i] )

    # Evaluate new PDF values
    combined_pdf_0[i] = dist_0.evaluatePDF( combined_e_loss_0[i] )
    combined_pdf_1[i] = dist_1.evaluatePDF( combined_e_loss_1[i] )

  return combined_cdf, combined_pdf_0, combined_pdf_1, combined_e_loss_0, combined_e_loss_1

def sub_unit_base_pdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower index, min, max and length
  index_0 = 0
  for i in range(0, len(e_losses_0) ):
    if e_losses_0[i] < e_loss:
      index_0 = i
    else:
      break
  e_loss_0_min = e_losses_0[index_0]
  e_loss_0_max = e_losses_0[index_0+1]
  e_loss_0_L = e_loss_0_max - e_loss_0_min

  # Get upper index, min, max and length
  index_1 = 0
  for i in range(0, len(e_losses_1) ):
    if e_losses_1[i] < e_loss:
      index_1 = i
    else:
      break
  e_loss_1_min = e_losses_0[index_1]
  e_loss_1_max = e_losses_0[index_1+1]
  e_loss_1_L = e_loss_1_max - e_loss_1_min

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = e_loss_max - e_loss_min

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper pdf values
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )

  # Calculate the pdf
  return lin_interp( E_alpha, pdf_0*e_loss_0_L, pdf_1*e_loss_1_L )/e_loss_L

def sub_unit_base_pdf_log( e_loss, dist_0, dist_1, e_losses, e_losses_0, e_losses_1, E_alpha):
  # Get the energy loss index
  index = 0
  for i in range(0, len(e_losses) ):
    if e_losses[i] < e_loss:
      index = i
    else:
      break

  # Get lower min, max and length
  e_loss_0_min = e_losses_0[index]
  e_loss_0_max = e_losses_0[index+1]
  e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

  # Get upper  min, max and length
  e_loss_1_min = e_losses_1[index]
  e_loss_1_max = e_losses_1[index+1]
  e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = e_losses[index]
  e_loss_max = e_losses[index+1]
  e_loss_L = numpy.log(e_loss_max/e_loss_min)

  # Calculate eta for E'
  eta = numpy.log( e_loss/e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
    e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

  # Calculate the lower and upper pdf values
  pdf_0 = dist_0.evaluatePDF( e_loss_0 )
  pdf_1 = dist_1.evaluatePDF( e_loss_1 )

  # Calculate the pdf
  return log_interp( E_alpha, pdf_0*e_loss_0_L, pdf_1*e_loss_1_L )/e_loss_L

def sub_unit_base_cdf_lin( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = (e_loss_0_max - e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = (e_loss_1_max - e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = lin_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = lin_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = (e_loss_max - e_loss_min)

  # Calculate eta for E'
  eta = ( e_loss - e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = e_loss_0_min + e_loss_0_L*eta
    e_loss_1 = e_loss_1_min + e_loss_1_L*eta

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = lin_interp( E_alpha, cdf_0, cdf_1 )

  return cdf

def sub_unit_base_cdf_log( e_loss, dist_0, dist_1, e_losses_0, e_losses_1, E_alpha):
  # Get lower min, max and length
  e_loss_0_min = e_losses_0[0]
  e_loss_0_max = e_losses_0[len(e_losses_0)-1]
  e_loss_0_L = numpy.log(e_loss_0_max/e_loss_0_min)

  # Get upper min, max and length
  e_loss_1_min = e_losses_1[0]
  e_loss_1_max = e_losses_1[len(e_losses_1)-1]
  e_loss_1_L = numpy.log(e_loss_1_max/e_loss_1_min)

  # Get E' min, max and length
  e_loss_min = log_interp( E_alpha, e_loss_0_min, e_loss_1_min )
  e_loss_max = log_interp( E_alpha, e_loss_0_max, e_loss_1_max )
  e_loss_L = numpy.log(e_loss_max/e_loss_min)

  # Calculate eta for E'
  eta = numpy.log( e_loss/e_loss_min )/e_loss_L

  # Calculate the lower and upper E' values
  if eta == 0.0:
    e_loss_0 = e_loss_0_min
    e_loss_1 = e_loss_1_min
  elif eta >= 1.0:
    e_loss_0 = e_loss_0_max
    e_loss_1 = e_loss_1_max
  else:
    e_loss_0 = numpy.exp( numpy.log(e_loss_0_min) + e_loss_0_L*eta )
    e_loss_1 = numpy.exp( numpy.log(e_loss_1_min) + e_loss_1_L*eta )

  # Calculate the lower and upper cdf values
  cdf_0 = dist_0.evaluateCDF( e_loss_0 )
  cdf_1 = dist_1.evaluateCDF( e_loss_1 )

  # Calculate the pdf
  cdf = log_interp( E_alpha, cdf_0, cdf_1 )

  return cdf



Utility.initFrensiePrng()

#datadir = '/home/software/mcnpdata/'
datadir = '/home/lkersting/frensie/src/packages/test_files/'

source = Teuchos.FileInputSource( datadir + '/cross_sections.xml' )
xml_obj = source.getObject()
cs_list = Teuchos.XMLParameterListReader().toParameterList( xml_obj )

# -------------------------------------------------------------------------- ##
#  Brem Data
# -------------------------------------------------------------------------- ##
# Possible Elements ['H-Native', 'Al-Native', 'Pb-Native']
elements = ['Al-Native']
# Possible Interpolation Policies ["LogLogLog", "LinLinLin", "LinLinLog"]
interps = ["LogLogLog"]
# Possible Interpolation Schemes ["Unit-base", "Processed CDF Unit-base", "Unit-base CDF", "Correlated Unit-base", "Corresponding Energies", "Cumulative Points"]
schemes = ["Unit-base CDF", "Correlated Unit-base", "Corresponding Energies"]
# Show Relative difference in schemes (True/False)
show_difference = False
# Possible energies [1e-2, 1e-1, 1.0, 15.7, 20.0]
energies = [0.0173]
# # Step length between plot points
step = 0.01
length = int(1.0/step)

plot_number = 1
for z in elements:
    print "\n----------------------------"
    print "-----", z, "Tests -----"
    print "----------------------------"
    data_list = cs_list.get( z )
    file_name = datadir + data_list.get( 'electroatomic_file_path' )
    native_data = Native.ElectronPhotonRelaxationDataContainer( file_name )
    energy_grid = native_data.getBremsstrahlungEnergyGrid()

    # Plot the given energy
    for energy in energies:
      print "Energy = ", energy
      print "----------------------------"

      # Find energy in energy grid (lower bin index)
      index = 0
      for i in range(0,len(energy_grid)):
        if energy_grid[i] < energy:
          index = i

      # Get lower and upper energy and the interpolation alpha
      energy_0 = energy_grid[index]
      energy_1 = energy_grid[index+1]
      lin_E_alpha = ( energy - energy_0 )/( energy_1 - energy_0 )
      log_E_alpha = numpy.log( energy/energy_0 )/numpy.log( energy_1/energy_0 )

      # Get lower energy bin data
      e_losses_0 = native_data.getBremsstrahlungPhotonEnergy(energy_0)
      pdfs_0 = native_data.getBremsstrahlungPhotonPDF(energy_0)
      dist_0 = Distribution.TabularDistribution_LinLin( e_losses_0, pdfs_0 )
      log_dist_0 = Distribution.TabularDistribution_LogLog( e_losses_0, pdfs_0 )

      # Get the lower cdf values
      cdfs_0 = numpy.zeros( shape=( len( pdfs_0 ) ) )
      log_pdfs_0 = numpy.zeros( shape=( len( pdfs_0 ) ) )
      log_cdfs_0 = numpy.zeros( shape=( len( pdfs_0 ) ) )
      for i in range(0,len(cdfs_0)):
        pdfs_0[i] = dist_0.evaluatePDF( e_losses_0[i] )
        log_pdfs_0[i] = dist_0.evaluatePDF( e_losses_0[i] )
        cdfs_0[i] = dist_0.evaluateCDF( e_losses_0[i] )
        log_cdfs_0[i] = dist_0.evaluateCDF( e_losses_0[i] )

      # Process lower energy bin data
      process_cdf_e_losses_0, process_cdf_pdfs_0, process_cdf_cdfs_0, process_cdf_dist_0 = construct_processed_cdf_distribution(pdfs_0, cdfs_0, e_losses_0)


      e_losses_1 = native_data.getBremsstrahlungPhotonEnergy(energy_1)
      pdfs_1 = native_data.getBremsstrahlungPhotonPDF(energy_1)
      dist_1 = Distribution.TabularDistribution_LinLin( e_losses_1, pdfs_1 )
      log_dist_1 = Distribution.TabularDistribution_LogLog( e_losses_1, pdfs_1 )

      # Get the upper cdf values
      cdfs_1 = numpy.zeros( shape=( len( pdfs_1 ) ) )
      log_pdfs_1 = numpy.zeros( shape=( len( pdfs_1 ) ) )
      log_cdfs_1 = numpy.zeros( shape=( len( pdfs_1 ) ) )
      for i in range(0, len(cdfs_1)):
        pdfs_1[i] = dist_1.evaluatePDF( e_losses_1[i] )
        log_pdfs_1[i] = dist_1.evaluatePDF( e_losses_1[i] )
        cdfs_1[i] = dist_1.evaluateCDF( e_losses_1[i] )
        log_cdfs_1[i] = dist_1.evaluateCDF( e_losses_1[i] )

      # Process upper energy bin data
      process_cdf_e_losses_1, process_cdf_pdfs_1, process_cdf_cdfs_1, process_cdf_dist_1 = construct_processed_cdf_distribution(pdfs_1, cdfs_1, e_losses_1)

      process_percent_values_0 = (e_losses_0 - e_losses_0[0])/(energy_0 - e_losses_0[0])
      process_percent_values_1 = (e_losses_1 - e_losses_1[0])/(energy_1 - e_losses_1[0])

      cdfs_0, pdfs_0, pdfs_1, e_losses_0, e_losses_1 = recast_distributions( dist_0, dist_1, cdfs_0, cdfs_1 )
      cdfs_1 = cdfs_0

      # Select distribution parameters
      for interp in interps:
        print "Interp = ", interp
        print "----------------------------"

        title = 'Bremsstrahlung Energy Loss PDF at ' + str(energy) + ' MeV'
        fig = plt.figure(num=plot_number, figsize=(10,5))
        if show_difference:
          plt.subplot2grid((2,1),(0, 0), colspan=5)
        else:
          plt.subplot2grid((1,1),(0, 0), colspan=5)
        plt.xlabel('Unit-base Energy Loss')
        plt.ylabel('PDF')
        plt.title( title )

        percent_values_0 = (e_losses_0 - e_losses_0[0])/(energy_0 - e_losses_0[0])
        percent_values_1 = (e_losses_1 - e_losses_1[0])/(energy_1 - e_losses_1[0])
        label0 = str(energy_0) +' MeV'
        log_label0 = "log "+label0
        process_cdf_label0 = "processed cdf "+label0

        label1 = str(energy_1) +' MeV'
        log_label1 = "log "+label1
        process_cdf_label1 = "processed cdf "+label1

        plt.plot( percent_values_0, pdfs_0, label=label0)
        plt.plot( percent_values_1, pdfs_1, label=label1)

        #plt.plot( process_percent_values_0, log_pdfs_0, label=log_label0)
        #plt.plot( process_percent_values_0, process_cdf_pdfs_0, label=process_cdf_label0)

        #plt.plot( process_percent_values_1, log_pdfs_1, label=log_label1)
        #plt.plot( process_percent_values_1, process_cdf_pdfs_1, label=process_cdf_label1)

        plot_number = plot_number + 1

        title = 'Bremsstrahlung Energy Loss CDF at ' + str(energy) + ' MeV'
        fig = plt.figure(num=plot_number, figsize=(10,5))
        if show_difference:
          plt.subplot2grid((2,1),(0, 0), colspan=5)
        else:
          plt.subplot2grid((1,1),(0, 0), colspan=5)
        plt.xlabel('Unit-base Energy Loss')
        plt.ylabel('CDF')
        plt.title( title )

        plt.plot( percent_values_0, cdfs_0, label=label0)
        #plt.plot( process_percent_values_0, log_cdfs_0, label=log_label0)
        #plt.plot( process_percent_values_0, process_cdf_cdfs_0, label=process_cdf_label0)

        plt.plot( percent_values_1, cdfs_1, label=label1)
        # plt.plot( process_percent_values_1, log_cdfs_1, label=log_label1)
        # plt.plot( process_percent_values_1, process_cdf_cdfs_1, label=process_cdf_label1)
        plt.xscale('log')
        plt.yscale('log')
        # plt.xlim(x_min,1.0)
        plt.legend( loc=3)

        # Plot all schemes on one graph
        pdfs = numpy.zeros(shape=( len(schemes), length ) )
        cdfs = numpy.zeros(shape=( len(schemes), length ) )
        e_losses = numpy.zeros(shape=( len(schemes) , length ) )

        combined_length = len( cdfs_0 )
        combined_cdfs = numpy.zeros(shape=( len(schemes), combined_length ) )
        combined_e_losses = numpy.zeros(shape=( len(schemes) , combined_length ) )

        e_loss_min = 0.0
        x_min = 1.0
        y_min = min(cdfs_0[1], cdfs_1[1])
        for n in range(0, len(schemes) ):
          print schemes[n]

          if schemes[n] == "Unit-base":
            if interp == "LogLogLog":
              e_loss_min = log_interp( log_E_alpha, e_losses_0[0],e_losses_1[0] )
              e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              for i in range(0, length):
                # Calculate the pdf
                pdfs[n,i] = unit_base_pdf_log( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)

                # Calculate the CDF
                if i == 0:
                  cdfs[n,i] = 0.0
                else:
                  cdfs[n,i] = unit_base_cdf_log( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)
                  cdfs[n,i] = (i+1.0)/(length + 1.0)
                  if e_losses[n,i] > 0.1*energy:
                    cdfs[n,i] = 0.5

            elif interp == "LinLinLin":
              e_loss_min = lin_interp( lin_E_alpha, e_losses_0[0],e_losses_1[0] )
              e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              for i in range(0, length):
                # Calculate the pdf
                pdfs[n,i] = unit_base_pdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

                # Calculate the CDF
                if i == 0:
                  cdfs[n,i] = 0.0
                else:
                  cdfs[n,i] = unit_base_cdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

            # Normalize to 1
            pdfs[n] = pdfs[n]/cdfs[n,length-1]
            cdfs[n] = cdfs[n]/cdfs[n,length-1]
            print cdfs[n]

          if schemes[n] == "Unit-base CDF":
            if interp == "LogLogLog":
              e_loss_min = log_interp( log_E_alpha, e_losses_0[0],e_losses_1[0] )
              e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              for i in range(0, length):
                # Calculate the pdf
                pdfs[n,i] = unit_base2_pdf_log( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)

                # Calculate the CDF
                if i == 0:
                  cdfs[n,i] = 0.0
                else:
                  cdfs[n,i] = unit_base2_cdf_log( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)

            elif interp == "LinLinLin":
              e_loss_min = lin_interp( lin_E_alpha, e_losses_0[0],e_losses_1[0] )
              e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              for i in range(0, length):
                # Calculate the pdf
                pdfs[n,i] = unit_base2_pdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

                # Calculate the CDF
                if i == 0:
                  cdfs[n,i] = 0.0
                else:
                  cdfs[n,i] = unit_base2_cdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

            # Normalize to 1
            pdfs[n] = pdfs[n]/cdfs[n,length-1]
            cdfs[n] = cdfs[n]/cdfs[n,length-1]

          if schemes[n] == "Processed CDF Unit-base":
            if interp == "LogLogLog":
              e_loss_min = log_interp( log_E_alpha, e_losses_0[0], e_losses_1[0] )
              e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              for i in range(0, length):
                # Calculate the processed e loss
                process_e_loss = numpy.log(e_losses[n,i])

                # Calculate the pdf
                pdfs[n,i] = unit_base2_pdf_lin( process_e_loss, process_cdf_dist_0, process_cdf_dist_1, process_cdf_e_losses_0, process_cdf_e_losses_1, log_E_alpha)

                # Calculate the CDF
                if i == 0:
                  cdfs[n,i] = 0.0
                else:
                  cdfs[n,i] = unit_base2_cdf_lin( process_e_loss, process_cdf_dist_0, process_cdf_dist_1, process_cdf_e_losses_0, process_cdf_e_losses_1, log_E_alpha)
              print cdfs[n]
            elif interp == "LinLinLin":
              print "The processed unit-base is naturally in logarithmic space!"
              # e_loss_min = lin_interp( lin_E_alpha, e_losses_0[0],e_losses_1[0] )
              # e_losses[n] = numpy.logspace(numpy.log10(e_loss_min), numpy.log10(energy), num=length)

              # for i in range(0, length):
              #   # Calculate the pdf
              #   pdfs[n,i] = unit_base2_pdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

              #   # Calculate the CDF
              #   if i == 0:
              #     cdfs[n,i] = 0.0
              #   else:
              #     cdfs[n,i] = unit_base2_cdf_lin( e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha)

            # Normalize to 1
            pdfs[n] = pdfs[n]/cdfs[n,length-1]
            cdfs[n] = cdfs[n]/cdfs[n,length-1]

          elif schemes[n] == "Correlated Unit-base":
            if n > 0:
              cdfs[n] = cdfs[n-1]
            else:
              cdfs[n] = numpy.logspace(numpy.log10(1e-12), numpy.log10(1.0), num=length)
            combined_cdfs[n] = cdfs_0
            if interp == "LogLogLog":
              for i in range(0, length):
                e_losses[n,i] = sample_correlated_unit_base_log( cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha )
                pdfs[n,i] = correlated_unit_base_pdf_log( cdfs[n,i], e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)
              for i in range(0, combined_length):
                combined_e_losses[n,i] = sample_correlated_unit_base_log( combined_cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha )
            elif interp == "LinLinLin":
              for i in range(0, length):
                e_losses[n,i] = sample_correlated_unit_base_lin( cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha )
                pdfs[n,i] = correlated_unit_base_pdf_lin( cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha )
              for i in range(0, combined_length):
                combined_e_losses[n,i] = sample_correlated_unit_base_lin( combined_cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha )

          elif schemes[n] == "Corresponding Energies":
            if n > 0:
              cdfs[n] = cdfs[n-1]
            else:
              cdfs[n] = numpy.logspace(numpy.log10(1e-12), numpy.log10(1.0), num=length)
            combined_cdfs[n] = cdfs_0
            if interp == "LogLogLog":
              for i in range(0, length):
                e_losses[n,i] = sample_direct_log( cdfs[n,i], dist_0, dist_1, log_E_alpha )
                pdfs[n,i] = correlated_direct_pdf_log( cdfs[n,i], e_losses[n,i], dist_0, dist_1, e_losses_0, e_losses_1, log_E_alpha)
              for i in range(0, combined_length):
                combined_e_losses[n,i] = sample_direct_log( combined_cdfs[n,i], dist_0, dist_1, log_E_alpha )
            elif interp == "LinLinLin":
              for i in range(0, length):
                e_losses[n,i] = sample_direct_lin( cdfs[n,i], dist_0, dist_1, lin_E_alpha )
                pdfs[n,i] = correlated_direct_pdf_lin( cdfs[n,i], dist_0, dist_1, e_losses_0, e_losses_1, lin_E_alpha )
              for i in range(0, combined_length):
                combined_e_losses[n,i] = sample_direct_lin( combined_cdfs[n,i], dist_0, dist_1, lin_E_alpha )

          elif schemes[n] == "Cumulative Points":

            combined_cdfs[n] = cdfs_0

            new_length = combined_length*2-2
            new_pdfs = numpy.zeros( shape=( new_length ) )
            new_e_losses = numpy.zeros( shape=( new_length ) )

            if interp == "LogLogLog":
              # Calculate the first point on the energy grid
              combined_e_losses[n,0] = log_interp( log_E_alpha, e_losses_0[0],e_losses_1[0] )
              new_e_losses[0] = combined_e_losses[n,0]
              new_e_losses[1] = new_e_losses[0] + 1e-12

              # Calculate new energy loss grid
              for i in range(1, combined_length-1):
                combined_e_losses[n,i] = log_interp( log_E_alpha, e_losses_0[i], e_losses_1[i] )
                j = 2*i
                new_e_losses[j] = combined_e_losses[n,i]
                new_e_losses[j+1] = new_e_losses[j] + 1e-12

              # Calculate the last point on the energy grid
              combined_e_losses[n,combined_length-1] = log_interp( log_E_alpha, e_losses_0[combined_length-1],e_losses_1[combined_length-1] )
              new_e_losses[new_length-1] = combined_e_losses[n,combined_length-1]

              # Calculate new PDF values
              for i in range(0, len( new_pdfs )):
                new_pdfs[i] = sub_unit_base_pdf_log( new_e_losses[i], dist_0, dist_1, combined_e_losses[n], e_losses_0, e_losses_1, log_E_alpha)

          plt.figure(plot_number - 1)
          label = interp + " " + schemes[n]

          # Normalize energy loss to percent
          if schemes[n] == "Cumulative Points":
            percent_values = (new_e_losses - new_e_losses[0])/(energy - new_e_losses[0])
            if new_pdfs.all() > 0.0:
              plt.plot( percent_values, new_pdfs, label=label)
              x_min = min( percent_values[2], x_min)
          else:
            percent_values = (e_losses[n] - e_losses[n,0])/(energy - e_losses[n,0])
            if pdfs[n].all() > 0.0:
              plt.plot( percent_values, pdfs[n], label=label)
              x_min = min( percent_values[1], x_min)

          plt.xscale('log')
          plt.yscale('log')
          plt.xlim(x_min,1.0)
          plt.legend( loc=3)

          plt.figure( plot_number )
          if schemes[n] == "Cumulative Points":
            percent_values = (combined_e_losses[n] - combined_e_losses[n,0])/(energy - combined_e_losses[n,0])
            plt.plot( percent_values, combined_cdfs[n], label=label)
          else:
            plt.plot( percent_values, cdfs[n], label=label)

          plt.figure( plot_number )
          plt.xscale('log')
          plt.yscale('log')
          plt.xlim(x_min,1.0)
          plt.ylim(y_min, 1.0)
          plt.legend( loc=4)

        # Plot Differences in Unit-Base and Unit-Base CDF
        if "Unit-base" in schemes and "Unit-base CDF" in schemes:
          plot_number = plot_number + 1
          fig1 = plt.figure(num=plot_number, figsize=(10,5))
          ax1 =plt.subplot2grid((2,1),(0, 0), colspan=5)
          plt.ylabel('Absolute Difference')
          plt.title( 'Differences' )

          i_0 = schemes.index("Unit-base")
          i_1 = schemes.index("Unit-base CDF")
          print i_0, i_1
          label = schemes[i_0] + ' vs ' + schemes[i_1]
          # Calculate difference between e_loss/cdfs
          cdf_diff = numpy.zeros(shape=length)
          pdf_diff = numpy.zeros(shape=length)
          for i in range(0,length):
            cdf_diff[i] = abs( cdfs[i_0,i] - cdfs[i_1,i] )
            pdf_diff[i] = abs( pdfs[i_0,i] - pdfs[i_1,i] )
          print cdf_diff
          print pdf_diff
          plt.plot( e_losses[i_0], cdf_diff, label="CDF Diff")
          plt.plot( e_losses[i_0], pdf_diff, label="PDF Diff")

          plt.xscale('log')
          plt.yscale('log')
          #plt.xlim(y_min,1.0)
          plt.legend( loc=4)

          plt.subplot2grid((2,1),(1, 0), colspan=5 )
          plt.ylabel('Absolute Relative Difference')
          plt.xlabel('Energy Loss')

          # Calculate rel. difference between e_loss/cdfs
          index = 0
          for i in range(0,length):
            if cdf_diff[i]*cdfs[i_0,i] > 0:
              cdf_diff[i] = cdf_diff[i]/cdfs[i_0,i]
              if cdf_diff[i] < 1e-15:
                cdf_diff[i] = 0.0
                index = index + 1
            else:
              index = index + 1

            if pdf_diff[i]*pdfs[i_0,i] > 0:
              pdf_diff[i] = pdf_diff[i]/pdfs[i_0,i]

          if index < length:
            plt.plot( e_losses[i_0], cdf_diff, label="CDF Rel. Diff")
          plt.plot( e_losses[i_0], pdf_diff, label="PDF Rel. Diff")
          print pdf_diff

          plt.xscale('log')
          plt.yscale('log')
          #plt.xlim(y_min,1.0)
          plt.legend( loc=4)

        # Plot differences
        if len(schemes) > 1:
          # Plot differences in energy loss
          plot_number = plot_number + 1
          fig1 = plt.figure(num=plot_number, figsize=(10,5))
          ax1 =plt.subplot2grid((2,1),(0, 0), colspan=5)
          plt.ylabel('Absolute Difference')
          plt.title( 'Differences in Energy Loss' )

          # Set up which plots to compare
          index = []; data = []; combined_data = []; index_found = False
          for n in range(0,len(schemes)):
            if cdfs[n,length-1] > 0.0 and combined_cdfs[n,combined_length-1] > 0.0 and not index_found:
              index = n
              index_found = True
            elif cdfs[n,length-1] > 0.0:
              data.append(n)
            elif combined_cdfs[n,combined_length-1] > 0.0:
              combined_data.append(n)

          if not index_found and len(data) > 1:
            index = data.pop(0)
            index_found = True
          elif not index_found and len(combined_data) > 1:
            index = combined_data.pop(0)
            index_found = True

          # Plot Differences
          for n in data:
            label = schemes[n] + ' vs ' + schemes[index]
            # Calculate difference between e_loss/cdfs
            for i in range(0,length):
              e_losses[n,i] = abs(e_losses[n,i] - e_losses[index,i])

            plt.plot( cdfs[n], e_losses[n], label=label)

          for n in combined_data:
            # Calculate difference between e_loss/cdfs
            for i in range(0,combined_length):
              combined_e_losses[n,i] = abs(combined_e_losses[n,i] - combined_e_losses[index,i])

            plt.plot( combined_cdfs[n], combined_e_losses[n], label=label)

          plt.xscale('log')
          plt.yscale('log')
          plt.xlim(y_min,1.0)
          plt.legend( loc=4)

          plt.subplot2grid((2,1),(1, 0), colspan=5 )
          plt.ylabel('Absolute Relative Difference')
          plt.xlabel('CDF')

          # Plot Relative Differences
          for n in data:
            # Calculate rel. difference between e_loss/cdfs
            for i in range(0,length):
              e_losses[n,i] = e_losses[n,i]/e_losses[index,i]

            label = schemes[n] + ' vs ' + schemes[index]
            plt.plot( cdfs[n], e_losses[n], label=label)

          for n in combined_data:
            # Calculate rel. difference between e_loss/cdfs
            for i in range(0,combined_length):
              combined_e_losses[n,i] = combined_e_losses[n,i]/combined_e_losses[index,i]

            label = schemes[n] + ' vs ' + schemes[index]
            plt.plot( combined_cdfs[n], combined_e_losses[n], label=label)

          plt.xscale('log')
          plt.yscale('log')
          plt.xlim(y_min,1.0)
          #plt.autoscale()
          #plt.xlim(cdfs[0],cdfs[len(cdfs)-1])
          #plt.ylim(y_min,y_max)
          plt.legend( loc=4)

        plot_number = plot_number + 1


plt.show()