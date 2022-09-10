import numpy as np
import scipy.stats as sp
import pandas as pd
from scipy.optimize import least_squares, curve_fit
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def read_input(input_file_name: str):
    skip_rows = list(range(52)) + list(range(53, 105))  # Which rows to skip
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
