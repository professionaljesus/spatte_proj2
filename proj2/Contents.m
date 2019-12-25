% Spatial Statistics with Image Analysis.
%
% general files:
%  fmsn20path - Set path to fmsn20-subdirectories.
%
% classification:
%  dirichletrnd - Sample from a Dirichlet distribution
%  gibbs_mu_sigma - Sample the posterior mean and covariance given y
%  kmeans - Classify data using the K-means algorithm.
%  manualclass - Mark image pixels as belonging to different classes.
%  normmix_classify - Classify data in a Gaussian mixture model.
%  normmix_em - Estimate parameters in a Gaussian mixture model.
%  normmix_gibbs - Sample parameters in a Gaussian mixture model using Gibbs.
%  normmix_initial - Computes initial estimates for K classes
%  normmix_kmeans - Use a single K-means to get a random inital estimate
%  normmix_posterior - Compute class probs. for a Gaussian mixture model.
%
% DMRF:
%  gibbs_alpha_beta - Samples a new set of alpha,beta given z using a MH-gibbs step
%  mrf_gaussian_est - Weighted estimate of Gaussian parameters
%  mrf_gaussian_post - Posterior alpha-parameters in a MRF with Gaussian data
%  mrf_icm - Estimate the MAP field from an MRF model
%  mrf_negLogPL - Computes the negative log pseudo-likelihood of a MRF.
%  mrf_ple - Pseudo-likelihood estimation of MRF parameters
%  mrf_sim - Simulate a samples for a MRF
%
% fields:
%  covest_ls - estimates a Matérn covariance with the Least Squares method.
%  covest_ml - estimates a Matérn covariance with the Maximum Likelihood method.
%  covest_nonparametric - nonparametric covariance estimator
%  distance_matrix - calculates the distance matrix for one or two sets of locations
%  matern_covariance - calculates Matérn covariances
%
% GMRF:
%  calc_gmrf_props - Helper function for calculating pointwise GMRF properties.
%  gmrf_mcmc_skeleton - Simulate GMRF parameters with MCMC, and perform calculations
%  gmrf_negloglike_Po_skeleton - Calculate the GMRF data likelihood, non-Gaussian observations
%  gmrf_negloglike_skeleton - Calculate the GMRF data likelihood, x integrated out
%  gmrf_param_hessian - Calculates the hessian for a negated log-likelihood
%  gmrf_param_map - Finds the MAP estimate in a simple field model.
%  gmrfprec - Constructs a precision matrix for a GMRF on a regular grid
%  gmrf_taylor_Po_skeleton - Taylor expansion of the conditional for non-Gaussian observations
%  igmrfprec - Constructs a precision matrix for a 1:st or 2:nd order IGMRF
%  igmrfprec_sphere - Constructs a precision matrix for spherical IGMRF
%  matern_prec_matrices - Calculate matrices need to build Matérn precisions
%  reorder_trigraph - Reorder nodes in a triangular graph for nice Cholesky
%
% graphics:
%  coordTrans - Transforms coordinates between long/lat, 3D and lambert.
%  globe_plot - Plot a field defined on a spherical grid 
%  hist2 - 2D histogram.
%  landsatimage - Make an RGB image matrix from LANDSAT data.
%  rgbimage - Make an RGB image from several weight images.
%  trisphere2anglegrid - Map a spherical field to a flat projection
%
% misc:
%  colstack - Stack columns of an image 
%  defstruct - Fill struct with default values if needed.
%  fillholes - Fill NaN-holes in an image.
%  helmert - Construct a Helmert (sub-)matrix.
%  icolstack - Invert column stacking of an image 
%  indicshape - Compute an indicator image for a deformed template.
%  ivec - The inverse of vec, anti-vectorise landmarks.
%  lanread - Read Landsat data file type .lan
%  noise01 - Modify binary image with noise.
%  pca - Principal component transformation of data.
%  polyimage - Computes an aliased indicator image for a polygon.
%  preshape - Compute the preshapes of landmark objects
%  reparameterise - Reparamatises a shape to obtain equidistant landmarks.
%  rotmat - Compute a 2D rotation matrix
%  simplespline - Compute a spline interpolation of a sequence of landmarks.
%  sphereHarmonics - Creates spherical harmonic functions.
%  triSphere - Triangulates a sphere.
%  vec - Vectorise a landmark matrix.
%
% shape:
%  gmrf_snake - Estimate a closed shape using a GMRF-snake model.
%  mark - Mark landmarks in an image.
%  procrustes_align - Perform Procrustes alignment of landmark data.
%  procrustes_dist - Compute Procrustes distances.
%  procrustes_mean - Estimate the Procrustes mean.
%  shape_tangent_inv - Compute pre-shapes from shape tangent coordinates.
%  shape_tangent - Compute vectorised shape space tangent coordinates.
%  snake_neg_loglike - Calculate negative log-likelihood for a snake.
%
% warp:
%  tps_prep - Precompute data for use in  tps_pull  and tps_push
%  tps_pull - Computes a pull warp
%  tps_push - Computes a push warp
%  tps_warp0 - Deform an image using a TPS warp.
%  tps_warp0_prep - Precompute data for use in  tps_warp0
%  tps_warp1 - Deform an image using an inverse TPS warp.
%  tps_warp1_prep - Precompute data for use in  tps_warp1
%
% $Date: 2014-12-10 12:13:39 +0100 (ons, 10 dec 2014) $

% $Id: Contents.m 4837 2014-12-10 11:13:39Z johanl $
