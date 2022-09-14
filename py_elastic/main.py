import glob
from pathlib import Path
import os

import pandas as pd

from utils import run_analysis_for_file


def main():
    file_or_dir = 'file'   # Set this to 'file' or to 'dir' for file or full directory of CSV files.

    # Parameters to set
    radius = 5.5  # Radius, r, in millimeters
    indent_depth = 1.21  # Indent depth, h, in millimeters

    if file_or_dir == 'file':
        # Settings
        input_file_name = Path('../reference/datasets/agae1per1-08232022-011721_tdf.csv')
        output_file_name = 'test.csv'
        figure_file_name = 'test.png'
        show_figure = True
        save_figure = True

        run_analysis_for_file(input_file_name,
                              output_file_name,
                              show_figure,
                              save_figure,
                              figure_file_name,
                              radius,
                              indent_depth)
    elif file_or_dir == 'dir':
        input_dir = Path('C:/dev/repos/personal/EricaMatlab/reference/datasets')
        output_dir = Path('C:/dev/repos/personal/EricaMatlab/reference/test_output')
        figure_file_name_base = 'fig'
        show_figure = True
        save_figure = True
        if not output_dir.exists():
            try:
                os.mkdir(output_dir)
            except Exception as e:
                print(f'ERROR: The output directory specified, {output_dir} does not exist '
                      f'and python does not have permissions to create it.')
                raise e
        input_files = glob.glob('*.csv', root_dir=input_dir)
        for input_file in input_files:
            base_name = input_file.lower().strip('.csv')
            in_path = input_dir / input_file
            out_path = output_dir / (base_name + '_result.csv')
            fig_path = output_dir / (figure_file_name_base + '_' + base_name + '.png')
            print(f'Running analysis for {in_path}')
            try:
                run_analysis_for_file(in_path,
                                      out_path,
                                      show_figure,
                                      save_figure,
                                      fig_path,
                                      radius,
                                      indent_depth)
                print(f'Results saved to {out_path}')
            except pd.errors.EmptyDataError as e:
                print(f'Unexpected error occurred for this file: {e}')

    else:
        raise ValueError('the file_or_dir variable should be set to \'file\' or \'dir\'. '
                         f'Instead it has value {file_or_dir=}')


if __name__ == '__main__':
    main()

print("Success!")
