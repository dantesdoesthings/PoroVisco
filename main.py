import numpy as np
import scipy.stats as sp
from scipy.optimize import least_squares
import pandas as pd
import math


def main():
    # Parameters
    radius = 12.7e-3
    input_file_name = 'C:\\dev\\repos\\personal\\EricaMatlab\\reference\\datasets\\0.2per alg gel em 08112022 111747_tdf.CSV'
    output_file_name = 'test.txt'
    save_figure = False
    skip_rows = list(range(52)) + [53, 104]  # Which rows to skip: First 52 are other data, the 2 others are junk
    col_names_original = ['Elapsed Time ', 'Load 2 ', 'Disp     ']
    col_names_rename = {
        'Elapsed Time ': 'elapsed',
        'Load 2 ': 'load',
        'Disp     ': 'disp'
    }
    col_names_new = ['elapsed', 'load', 'disp']
    start_index = 50  # The blank row in the input data, "index" in the old matlab code

    # Pull in the data
    data = pd.read_csv(input_file_name, skiprows=skip_rows, usecols=col_names_original)
    # Rename and reorder for easier use and access
    data.rename(columns=col_names_rename, inplace=True)
    data = data[col_names_new]
    # Adjust the data as needed
    data['load'] *= -1
    data['disp'] *= 1e-3
    # Start calculations #####
    # Separate the ramp and relax times
    ramp_time = data['elapsed'][start_index] - data['elapsed'][0]
    hold_time = data['elapsed'].iloc[-1] - data['elapsed'][start_index]
    # Run the calcs
    plot_flag = 1
    i = 1  # Loop through all the files, TODO: Remove
    legend = 'test'  # TODO: Fix/remove
    result = poro_visco_elastic_model(data,
                                      radius,
                                      ramp_time,
                                      hold_time,
                                      plot_flag,
                                      i,
                                      legend,
                                      start_index)
    # print(data, ramp_time)
    # print(data)


def poro_visco_elastic_model(test_data: pd.DataFrame,
                             radius,
                             ramp_time,
                             hold_time,
                             plot_flag,
                             i,
                             leg,
                             dwellid):
    # TODO: All of the x_ values seem to be like 3 actual values, then a list afterwards.
    #  Seperate into appropriate things if needed after verifying functionality
    elapsed = test_data['elapsed'].values
    load = test_data['load'].values
    disp = test_data['disp'].values
    # Stats from the input
    hmax = abs(disp[dwellid])
    pmax = load[dwellid]
    # Parameters
    m = 1.5
    av = 8 * radius**0.5 / 3
    ap = (radius * hmax)**0.5
    # Poroelastic Initialization
    p0_guess = load[1]
    p_inf_guess = load[-1]
    g_guess = 3*p0_guess/(16 * hmax * ap)
    v_guess = 1 - (p0_guess / (2 * p_inf_guess))
    kappa_guess = 1e-14
    log_kappa_guess = np.log10(kappa_guess)
    d_guess = (2 * (1 - v_guess) / (1 - 2 * v_guess)) * g_guess * kappa_guess
    log_d_guess = np.log10(d_guess)

    # Variable bounds
    if v_guess < 0:
        vmin = 0.5 / v_guess
        vmax = -1 / v_guess
    else:
        vmin = -1/v_guess
        vmax = 0.5/v_guess

    lower_bound1 = [0, 0, vmin]
    upper_bound1 = [1, 1, vmax]
    x01 = [1, 1, 1]
    xg1 = [log_d_guess, p_inf_guess, v_guess]

    # Viscoelastic initialization
    num_params = 2
    c0_guess = load[-1] / (hmax**m * av)
    c_guess = c0_guess * np.ones(num_params)
    g_guess = (c0_guess + np.sum(c_guess)) / 2
    t_guess = [ramp_time, hold_time]
    x02 = np.ones([1, 2 * num_params + 1])
    xg2 = [c0_guess, c_guess, t_guess]
    lower_bound2 = [0, 0, 0, 0, 0]
    upper_bound2 = [1, 1, 1, 1, 1]

    # Opt run
    x0 = [x01, x02]
    x_guess = [xg1, xg2]
    lower_bound = [lower_bound1, lower_bound2]
    upper_bound = [upper_bound1, upper_bound2]
    # options.optim = optimset('MaxFunEvals', 500000, 'Display', 'none', 'TolX', 1E-30);
    # X = lsqnonlin( @ OBJPoroVisco_mri, X0, LB, UB, options.optim, DT, Xguess, rampTime, hmax, av, m, R);
    x = least_squares(poro_visco_optimization, x0, bounds=[lower_bound, upper_bound],
                      args=[test_data, x_guess, ramp_time, hmax, av, m, radius],
                      max_nfev=500000, xtol=1e-30)

    # Poroelastic parameter estimation
    log_d_id = np.real(x[0] * log_d_guess)
    d_id = 10**log_d_id
    p_inf = x[1] * p_inf_guess
    v = x[2] * v_guess
    p0 = p_inf * (2 * (1 - v))
    g = 3 * p0 / (16 * hmax * ap)
    kappa = d_id * (1 - 2 * v) / (2 * (1 - v) * g)
    kappa = 0.89 * 1e-3 * kappa
    e = 2 * g * 1.5

    # Viscoelastic parameter estimation
    x_id = np.real(x[3:] * x_guess[3:])
    c0 = x_id[0]
    c = x_id[1: num_params + 1]
    t = x_id[num_params + 1:]
    rcf = t / ramp_time * (np.exp(ramp_time / t) - 1)
    c1 = c[0] / rcf[0]
    c2 = c[1] / rcf[1]
    g0 = (c0 + np.sum(c)) / 2
    e0 = 2 * g0 * 1.5
    g_inf = c0 / 2
    e_inf = 2 * g_inf * 1.5
    tau1 = t

    # Load relaxation data
    num_points = len(elapsed)
    tau = [d_id * elapsed[k] / (radius * hmax) for k in range(len(num_points))]
    g = [0.491 * np.exp(-0.908 * np.sqrt(tau[k])) + 0.509 * np.exp(-1.679 * (tau[k])) for k in range(num_points)]
    fitP = [g[k] * (p0 - p_inf) + p_inf for k in range(num_points)]
    fitV = [c0 + c1 * np.exp(-1 * elapsed[k] / t[0]) + c2 * np.exp(-1 * elapsed[k] / t[1]) for k in range(num_points)]
    fitV = [fitV[k] * av * hmax ** 1.5 for k in range(num_points)]
    fit = np.array([fitV[k] * fitP[k] / load[-1] for k in range(num_points)])

    s = list(range(num_points))
    rsq_visco = coefficient_determination(load, fitV)
    rsq_poro = coefficient_determination(load, fitP)
    rsq_both = coefficient_determination(load, fit)

    pp1 = [g0 / 1000, e0 / 1000, g_inf / 1000, e_inf / 1000, g_inf / g0, tau[0], tau[1], rsq_visco]
    pp2 = [g / 1000, e / 1000, v, kappa, d_id, rsq_poro, rsq_both]
    pp = [pp1, pp2]
    dp = [test_data, fit]

    # PLOT

def coefficient_determination(experimental_data: np.ndarray,
                              fit_data: np.ndarray,
                              num_regressors: int = 0):
    num_points = len(experimental_data)
    if len(fit_data) != num_points:
        raise ValueError('Experimental data and fit data have unequal lengths.')

    # Apparently there's a "numpy.corrcoef" function that probably does most of this
    # you can look it up later to simplify if you want
    sse = np.sum((experimental_data - fit_data) ** 2)
    sst = np.sum((experimental_data - np.mean(experimental_data)) ** 2)
    r_square = 1 - sse/sst
    adjusted_r_square = 1 - (1 - r_square) * (num_points - 1) / (num_points - num_regressors - 1)
    return r_square, adjusted_r_square


def poro_visco_optimization(x,
                            exp_data,
                            x_guess,
                            rise_time,
                            hmax,
                            a,
                            m,
                            radius):
    # Setup
    x1 = x[:3]
    x2 = x[3:]
    xg1 = x_guess[:3]
    xg2 = x_guess[3:]
    # Poroelastic
    log_d = x1[0] * xg1[0]
    p_inf = x1[1] * xg1[1]
    v = x1[2] * xg1[2]
    d = 10**log_d
    p0 = p_inf*(2*(1-v))
    # Viscoelastic
    num_param = (len(xg2) - 1) / 2
    c0 = x2[0] * xg2[0]
    c = x2[1:num_param + 1] * xg2[1: num_param + 1]
    t = x2[num_param + 1:] * xg2[num_param + 1:]
    rcf = t / rise_time * (np.exp(rise_time / t) - 1)  # Check this in particular since the original used ./ and .* and exp()
    c1 = c[0] / rcf[0]
    c2 = c[1] / rcf[1]

    tau = [d * exp_data[k] / (radius * hmax) for k in range(len(exp_data))]
    g = [0.491 * np.exp(-0.908 * np.sqrt(tau[k])) + 0.509 * np.exp(-1.679 * (tau[k])) for k in range(len(exp_data))]
    fitP = [g[k] * (p0 - p_inf) + p_inf for k in range(len(exp_data))]
    fitV = [c0 + c1 * np.exp(-1 * exp_data[k][0] / t[0]) + c2 * np.exp(-1 * exp_data[k][0] / t[1]) for k in range(len(exp_data))]
    fitV = [fitV[k] * a * hmax ** 1.5 for k in range(len(exp_data))]
    fit = np.array([fitV[k] * fitP[k] / exp_data[-1][1] for k in range(len(exp_data))])
    err_vect = fit - exp_data[:][1] / np.mean(exp_data[:][1])  # This is probably getting the column wrong

    g0 = (c0 + np.sum(c)) / 2
    ap = np.sqrt(radius * hmax)
    g = 3 * p0 / (16 * hmax * ap)
    if (g - g0) / g0 > 0.05:
        err_vect *= 1e6
    return err_vect


if __name__ == '__main__':
    main()
