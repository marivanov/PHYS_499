#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 14:06:06 2017

@author: mivanov
"""
import numpy as np

def bayes_linear(x, y, x_err, y_err, p1, nWalkers=10, nBurn=100, nSample=1000,
                 conf_interval=[15.9, 84.1], verbose=False,
                 return_samples=False):
    '''
    Fit a line with errors in both variables using MCMC.
    Original version of this function is Erik Rosolowsky's:
    https://github.com/low-sky/py-low-sky/blob/master/BayesLinear.py
    Parameters
    ----------
    x : `~numpy.ndarray`
        x data.
    y : `~numpy.ndarray`
        y data.
    x_err : `~numpy.ndarray`
        x errors.
    y_err : `~numpy.ndarray`
        y errors.
    nWalkers : int, optional
        Number of walkers in the sampler (>2 required). Defaults to 10.
    nBurn : int, optional
        Number of steps to burn chain in for. Default is 100.
    nSample : int, optional
        Number of steps to sample chain with. Default is 1000.
    conf_interval : list, optional
        Upper and lower percentiles to estimate the bounds on the parameters.
        Defaults to the 1-sigma intervals (34.1% about the median).
    verbose : bool, optional
        Plot the resulting fit.
    return_samples : bool, optional
        Returns the entire chain of samples, when enabled.
    Returns
    -------
    params : `~numpy.ndarray`
        Fit parameters (slope, intercept)
    errors : `~numpy.ndarray`
        Confidence interval defined by values given in `conf_interval`
        (slope, intercept).
    samples : `~numpy.ndarray`
        Samples from the chain. Returned only when `return_samples` is enabled.
    '''

    try:
        import emcee
    except ImportError:
        raise ImportError("emcee must be installed to use Bayesian fitting.")

    def _logprob(p, x, y, x_err, y_err):
        theta, b = p[0], p[1]
        if np.abs(theta - np.pi / 4) > np.pi / 4:
            return -np.inf
        Delta = (np.cos(theta) * y - np.sin(theta) * x - b * np.cos(theta))**2
        Sigma = (np.sin(theta))**2 * x_err**2 + (np.cos(theta))**2 * y_err**2
        lp = -0.5 * np.nansum(Delta / Sigma) - 0.5 * np.nansum(np.log(Sigma))

        return lp

    ndim = 2
    p0 = np.zeros((nWalkers, ndim))
    p0[:, 0] = np.arctan(p1[0]) + np.random.randn(nWalkers) * 0.1
    p0[:, 1] = np.random.randn(nWalkers) + p1[1]
    sampler = emcee.EnsembleSampler(nWalkers, ndim, _logprob,
                                    args=[x, y, x_err, y_err])
    pos, prob, state = sampler.run_mcmc(p0, nBurn)
    sampler.reset()
    sampler.run_mcmc(pos, nSample)

    slopes = np.tan(sampler.flatchain[:, 0])
    intercepts = sampler.flatchain[:, 1]

    slope = np.median(slopes)
    intercept = np.median(intercepts)

    params = np.array([slope, intercept])

    # Use the percentiles given in conf_interval
    error_intervals = np.empty((2, 2))
    error_intervals[0] = np.percentile(slopes, conf_interval)
    error_intervals[1] = np.percentile(intercepts, conf_interval)

    return params, error_intervals, np.vstack([slopes, intercepts])