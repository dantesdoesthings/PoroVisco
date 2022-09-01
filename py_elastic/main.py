import pandas as pd

from utils import read_input, poro_visco_elastic_model


def main():
    # Parameters
    radius = 12.7e-3
    input_file_name = '../reference/datasets/0.2per alg gel em 08112022 111747_tdf.CSV'
    output_file_name = 'test.csv'
    save_figure = False
    # Start calculations #####
    # Separate the ramp and relax times
    start_index = 50  # The blank row in the input data, "index" in the old matlab code
    elapsed, load, disp = read_input(input_file_name)
    ramp_time = elapsed[start_index] - elapsed[0]
    hold_time = elapsed[-1] - elapsed[start_index]
    # Run the calcs
    plot_flag = 1
    i = 1  # Loop through all the files, TODO: Remove, deal with how the old code ran multiple files at once
    legend = ['G0 (kPa)','E0 (kPa)','Ginf (kPa)','Einf (kPa)','Ginf/G0','T1','T2','R-SQ, VISCO','G (kPa)','E (kPa)','nu','k (m^2)','D (m^2/s)','R-SQ, PORO', 'R-SQ, PVE']
    pp1, pp2, fit = poro_visco_elastic_model(
        elapsed, load, disp,
        radius, ramp_time, hold_time, plot_flag, i, legend, start_index)
    result_df1 = pd.DataFrame({
        'G0 (kPa)': pp1[0],
        'E0 (kPa)': pp1[1],
        'Ginf (kPa)': pp1[2],
        'Einf (kPa)': pp1[3],
        'Ginf/G0': pp1[4],
        'T1': pp1[5],
        'T2': pp1[6],
        'R-SQ, VISCO': pp1[7],
        'G (kPa)': pp2[0],
        'E (kPa)': pp2[1],
        'nu': pp2[2],
        'k (m^2)': pp2[3],
        'D (m^2/s)': pp2[4],
        'R-SQ, PORO': pp2[5],
        'R-SQ, PVE': pp2[6]
    }, index=[0])
    result_df1.to_csv(output_file_name, index=None)


if __name__ == '__main__':
    main()
