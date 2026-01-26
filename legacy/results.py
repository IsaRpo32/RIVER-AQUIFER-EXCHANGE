"""
Results Module (Legacy)

NOTE: This is legacy code preserved for reference.
The active functionality has been moved to src/modflow_utils.py
"""

import os
import flopy
import warnings
import numpy as np

from .geology import Geology
from .hydrology import Hydrology

from flopy.utils.sfroutputfile import SfrFile


class Results(Geology, Hydrology):
    def __init__(self):
        super().__init__()
        self.figsize = (15, 5)

    def check_data(self):
        """Check whether there are properties missing or not"""
        missing = [attr for attr in ["topography", "basic_river", "kinoshita_river"]
                  if not hasattr(self, attr)]

        if len(missing) > 0:
            for attr in missing:
                warnings.warn(f"Attribute '{attr}' is missing. If it doesn't worry, omit this message.")

    def model_error(self, verbose=True):
        lst_file_path = os.path.join(self.path_results, self.model.name + ".list")
        mfl = flopy.utils.MfListBudget(lst_file_path)

        incrementaldf, cumulativedf = mfl.get_dataframes()
        if verbose:
            print("Cumulative error STEADY:", str(cumulativedf.iloc[0,-1])+ "%", "in L**3")
            print("Incremental error STEADY:", str(incrementaldf.iloc[0,-1])+ "%", "in L**3/T")
            print("Cumulative error TRANSIENT:", str(cumulativedf.iloc[1,-1])+ "%", "in L**3")
            print("Incremental error TRANSIENT:", str(incrementaldf.iloc[1,-1])+ "%", "in L**3/T")

            return None

        return incrementaldf, cumulativedf

    @staticmethod
    def adjust_coordinates(i_riv, j_riv, orientation="v"):
        """Helper function to adjust river coordinates based on orientation"""
        if orientation == 'h':
            idx_sort = np.lexsort((i_riv, j_riv))
            i_riv, j_riv = j_riv[idx_sort], i_riv[idx_sort]

        idx_place, idx_unique_counts = np.unique(i_riv, return_counts=True, return_index=True)[1:]
        idx_start_repeats = np.ones(i_riv.shape)
        if len(idx_place) > 0:
            idx_start_repeats[idx_place-1] = 0

        for n in range(len(idx_start_repeats) - 1):
            if idx_start_repeats[n] == 1 and idx_start_repeats[n+1] == 1:
                continue

        return i_riv, j_riv

    def extract_sfr_data(self, reach_numbers, labels, orientation='v'):
        """
        Extract SFR data for specific reaches without plotting.

        Parameters
        ----------
        reach_numbers : int or list
            Reach number(s) to extract (0-indexed)
        labels : list
            SFR variables to extract ['Qin', 'Qaquifer', 'Qout', 'stage', 'depth', 'gradient']
        orientation : str
            River orientation ('v' or 'h')

        Returns
        -------
        dict
            Dictionary with data structured by reach and variable
        """
        sfr_df = SfrFile(f"{self.model.model_ws}\\{self.model.name}.sfr.out").df

        if isinstance(reach_numbers, int):
            reach_numbers = [reach_numbers]

        file_reaches = [rn + 1 for rn in reach_numbers]
        filtered_df = sfr_df[sfr_df['reach'].isin(file_reaches)]

        data_dict = {}

        for file_rn, orig_rn in zip(file_reaches, reach_numbers):
            reach_df = filtered_df[filtered_df['reach'] == file_rn]

            data_dict[orig_rn] = {
                'reach_id': orig_rn,
                'file_reach_id': file_rn,
                'n_steps': len(reach_df),
                'time_steps': np.arange(len(reach_df)),
                'data': {}
            }

            for label in labels:
                if label in reach_df.columns:
                    values = reach_df[label].to_numpy()
                    data_dict[orig_rn]['data'][label] = {
                        'values': values,
                        'max_value': values.max(),
                        'min_value': values.min(),
                        'max_index': np.argmax(values),
                        'min_index': np.argmin(values),
                        'mean': values.mean(),
                        'std': values.std()
                    }
                else:
                    print(f"Warning: Variable '{label}' not found in SFR data")

        return data_dict
