import os

from calculations.normalization import read_edi_files, mtedi, normalize_rho
from ui.DataWidget import MTComponent


class Normalization:
    def __init__(self, edis: list[mtedi.Edi], mt_points: int, period: float, file_paths: list = None):
        self.edis = edis
        self.mt_points = mt_points
        self.period = period
        self.file_paths = file_paths

        periods = self.edis[0].Z.res_xy
        if periods[0] >= period:
            self.period_index = 0
        elif periods[-1] <= period:
            self.period_index = len(periods) - 1
        else:
            self.period_index = 0
            for i, x in enumerate(periods):
                if x <= period:
                    self.period_index = i
                else:
                    break

        self.result_edis = normalize_rho(self.edis, self.period_index, self.mt_points)
        self.data_widget = MTComponent(self, self.result_edis)

    def save_results(self, dir_path):
        names = [os.path.splitext(os.path.basename(i))[0] for i in self.file_paths]

        for i, edi in enumerate(self.result_edis):
            edi.write_edi_file(f'{dir_path}/{names[i]}.edi')

    def return_parameters(self):
        return self.edis, self.mt_points, self.period
