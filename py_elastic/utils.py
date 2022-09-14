import numpy as np
import scipy.stats as sp
import pandas as pd
from scipy.optimize import least_squares, curve_fit
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def read_input(input_file_name: str):
    skip_rows = list(range(52)) + list(range(53, 105))  # Which rows to skip ADJUSTED TO FOR 9.17.20 (Removed)
    col_names_original = ['Elapsed Time ', 'Load 2 ', 'Disp     ']
    col_names_rename = {
        'Elapsed Time ': 'elapsed',
        'Load 2 ': 'load',
        'Disp     ': 'disp'
    }

    # Pull in the data
    data = pd.read_csv(input_file_name, skiprows=skip_rows, usecols=col_names_original)
    # Rename and reorder for easier use and access
    data.rename(columns=col_names_rename, inplace=True)
    # Adjust the data as needed
    data['load'] *= -1
    data['disp'] *= 1e-3
    return data['elapsed'].values, data['load'].values, data['disp'].values



def poro_visco_elastic_model(elapsed: np.ndarray,
                             load: np.ndarray,
                             disp: np.ndarray,
                             radius,
                             ramp_time,
                             hold_time):
    # Stats from the input
    hmax = abs(disp[0])
    # Parameters
    m = 1.5
    av = 8 * radius**0.5 / 3
    ap = (radius * hmax)**0.5
    # Poroelastic Initialization
    p0_guess = load[0]
    p_inf_guess = load[-1]
    g_guess = 3*p0_guess/(16 * hmax * ap)
    v_guess = 1 - (p0_guess / (2 * p_inf_guess))
    kappa_guess = 1e-14
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
    t_guess = [ramp_time, hold_time]
    x02 = np.ones([2 * num_params + 1])
    xg2 = np.concatenate([[c0_guess], c_guess, t_guess])
    lower_bound2 = [0, 0, 0, 0, 0]
    upper_bound2 = [1, 1, 1, 1, 1]

    # Opt run
    x_guess = np.concatenate([xg1, xg2])
    lower_bound = np.concatenate([lower_bound1, lower_bound2])
    upper_bound = np.concatenate([upper_bound1, upper_bound2])
    x0 = (lower_bound + upper_bound)/2 #Redefined upper bounds and lower bounds to be averaged (remove x0 bound error)
    # print(type(x0), type(elapsed), type(load), type(disp), type(x_guess), type(ramp_time), type(hmax), type(av), type(m), type(radius))
    # options.optim = optimset('MaxFunEvals', 500000, 'Display', 'none', 'TolX', 1E-30);
    # X = lsqnonlin( @ OBJPoroVisco_mri, X0, LB, UB, options.optim, DT, Xguess, rampTime, hmax, av, m, R);
    #ERICA ERROR LOG: "Each lower bound must be less than 0". Zeros in data from lack of force, change in data file to -0.0001
    opt_result = least_squares(poro_visco_optimization, x0, bounds=[lower_bound, upper_bound],
                               args=[elapsed, load, x_guess, ramp_time, hmax, av, m, radius],
                               max_nfev=500000, xtol=1e-30)
    exp_decay_result, _ = curve_fit(exp_decay_function, elapsed, load, bounds=(0, [np.inf, np.inf, 2 * load[-1]]))
    print(exp_decay_result)
    x = opt_result.x

    # Poroelastic parameter estimation
    log_d_id = np.real(x[0] * log_d_guess)
    d_id = 10**log_d_id
    p_inf = x[1] * p_inf_guess
    v = x[2] * v_guess
    p0 = p_inf * (2 * (1 - v))
    G = 3 * p0 / (16 * hmax * ap)
    kappa = d_id * (1 - 2 * v) / (2 * (1 - v) * G)
    kappa = 0.89 * 1e-3 * kappa
    e = 2 * G * 1.5

    # Viscoelastic parameter estimation
    x_id = np.real(x[3:] * x_guess[3:])
    c0 = x_id[0]
    c = x_id[1: num_params + 1]
    t = x_id[num_params + 1:]
    rcf = t / ramp_time * (np.exp(ramp_time / t) - 1)
    c1 = c[0] / rcf[0]
    c2 = c[1] / rcf[1]
    g0 = (c0 + np.sum(c)) / 2
    e0 = exp_decay_function(np.array([0]), *exp_decay_result)
    e_inf = exp_decay_result[2]
    g_inf = c0 / 2
    tau1 = t

    # Load relaxation data
    num_points = len(elapsed)
    tau = [d_id * elapsed[k] / (radius * hmax) for k in range(num_points)]
    g = [0.491 * np.exp(-0.908 * np.sqrt(tau[k])) + 0.509 * np.exp(-1.679 * (tau[k])) for k in range(num_points)]
    fitP = [g[k] * (p0 - p_inf) + p_inf for k in range(num_points)]
    fitV = [c0 + c1 * np.exp(-1 * elapsed[k] / t[0]) + c2 * np.exp(-1 * elapsed[k] / t[1]) for k in range(num_points)]
    fitV = [fitV[k] * av * hmax ** 1.5 for k in range(num_points)]
    fit = np.array([fitV[k] * fitP[k] / load[-1] for k in range(num_points)])

    tau = np.array(tau)
    g = np.array(g)
    fitP = np.array(fitP)
    fitV = np.array(fitV)

    rsq_visco, _ = coefficient_determination(load, fitV)
    rsq_poro, _ = coefficient_determination(load, fitP)
    rsq_both, _ = coefficient_determination(load, fit)

    pp1 = [g0 / 1000, e0 / 1000, g_inf / 1000, e_inf / 1000, g_inf / g0, tau[0], tau[1], rsq_visco]
    pp2 = [G / 1000, e / 1000, v, kappa, d_id, rsq_poro, rsq_both]

    # PLOT
    return pp1, pp2, fit, exp_decay_result


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
                            elapsed,
                            load,
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
    num_param = (len(xg2) - 1) // 2
    c0 = x2[0] * xg2[0]
    c = x2[1:num_param + 1] * xg2[1: num_param + 1]
    t = x2[num_param + 1:] * xg2[num_param + 1:]
    rcf = t / rise_time * (np.exp(rise_time / t) - 1)  # Check this in particular since the original used ./ and .* and exp()
    c1 = c[0] / rcf[0]
    c2 = c[1] / rcf[1]

    num_points = len(elapsed)
    tau = [d * elapsed[k] / (radius * hmax) for k in range(num_points)]
    g = [0.491 * np.exp(-0.908 * np.sqrt(tau[k])) + 0.509 * np.exp(-1.679 * (tau[k])) for k in range(num_points)]
    fitP = [g[k] * (p0 - p_inf) + p_inf for k in range(num_points)]
    fitV = [c0 + c1 * np.exp(-1 * elapsed[k] / t[0]) + c2 * np.exp(-1 * elapsed[k] / t[1]) for k in range(num_points)]
    fitV = [fitV[k] * a * hmax ** 1.5 for k in range(num_points)]
    fit = np.array([fitV[k] * fitP[k] / load[-1] for k in range(num_points)])
    err_vect = (fit - load) / np.mean(load)  # This is probably getting the column wrong

    g0 = (c0 + np.sum(c)) / 2
    ap = np.sqrt(radius * hmax)
    g = 3 * p0 / (16 * hmax * ap)
    if (g - g0) / g0 > 0.05:
        err_vect *= 1e6
    return err_vect


def exp_decay_function(x: np.ndarray, a: float, b: float, c: float):
    result = a * np.exp(-b * x) + c
    return result


def plot_results(time_vals: np.ndarray,
                 point_data: np.ndarray,
                 exp_decay_result: np.ndarray,
                 show_plot: bool = False,
                 save_plot: bool = False,
                 plot_path: str = 'ResultPlot.png',
                 exp_fit=None):
    fig = plt.figure()  # type: Figure
    p1 = fig.add_subplot()
    p1.set_ylabel('load')
    p1.set_xlabel('elapsed time')
    p1.plot(time_vals, point_data, 'o', color='#004466', label='original data')
    p1.plot(time_vals, exp_decay_function(time_vals, *exp_decay_result), '-r', label='exponential decay fit')
    p1.legend()
    if exp_fit:
        y_data = np.exp(exp_fit.slope * time_vals + exp_fit.intercept)
        p1.plot(time_vals, y_data, '-b', label='Testing')
    if save_plot:
        plt.savefig(plot_path)
    if show_plot:
        plt.show()


def run_analysis_for_file(
        input_file_name: str,
        output_file_name: str,
        show_figure: bool = False,
        save_figure: bool = False,
        figure_file_name: str = None,
        radius: float = 5.5,
        indent_depth: float = 1.21
):
    # Start calculations
    # Separate the ramp and relax times
    start_index = 50  # The blank row in the input data, "index" in the old matlab code
    elapsed, load, disp = read_input(input_file_name)
    ramp_time = elapsed[start_index] - elapsed[0]
    hold_time = elapsed[-1] - elapsed[start_index]
    # Run the calcs
    legend = ['G0 (kPa)',
              'E0 (kPa)',
              'Ginf (kPa)',
              'Einf (kPa)',
              'Ginf/G0',
              'T1',
              'T2',
              'R-SQ, VISCO',
              'G (kPa)',
              'E (kPa)',
              'nu',
              'k (m^2)',
              'D (m^2/s)',
              'R-SQ, PORO',
              'R-SQ, PVE']
    f_0, f_inf, e_0, e_inf, exp_decay = find_e_vals(
        elapsed, load, disp, radius, indent_depth)

    # Gather and save the results
    result_df1 = pd.DataFrame({
        'F0': f_0,
        'E0 (kPa)': e_0,
        'F_inf': f_inf,
        'E_inf (kPa)': e_inf,
    }, index=[0])
    result_df1.to_csv(output_file_name, index=None)
    plot_results(elapsed, load, exp_decay, show_figure, save_figure, figure_file_name)


def find_e_vals(elapsed: np.ndarray,
                load: np.ndarray,
                disp: np.ndarray,
                radius: float,
                indent_depth: float):
    # Find an exponential decay fit
    exp_decay_result, _ = curve_fit(exp_decay_function, elapsed, load, bounds=(0, [np.inf, np.inf, 2 * load[-1]]))
    # Convert radius and indent_depth from mm to meters
    radius *= 1e-3
    indent_depth *= 1e-3
    # Find F_0 as the approx displacement from the fit at time=0
    f_0 = exp_decay_function(elapsed[0], *exp_decay_result)
    # Find F_inf from the same fit
    f_inf = exp_decay_result[2]
    # Use calc_e to find the E values
    e_0 = calc_e(f_0, radius, indent_depth)
    e_inf = calc_e(f_inf, radius, indent_depth)
    return f_0, f_inf, e_0, e_inf, exp_decay_result


def calc_e(f, r, h, poisson_ratio=0.5):
    e_star = 3 * f / (4 * r ** 0.5 * h ** 1.5)
    e = (1 - poisson_ratio ** 2) * e_star
    # Convert to KPa from Pa
    e *= 1e-3
    return e


def serialize(obj):

    # Serialize a numpy array
    if isinstance(obj, np.ndarray):
        return obj.tolist()

    # Serialize a numpy matrix
    if isinstance(obj, np.matrix):
        return obj.tolist()

    # Attempt to serialize a number
    try:
        return float(obj)
    except Exception:
        pass

    return obj.__dict__
