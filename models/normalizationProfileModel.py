import os

from calculations.normalization import read_edi_files, mtedi, normalize_rho
from models.normalizationModel import Normalization
from ui.DataWidget import MTComponent


class NormalizationProfileModel:
    """
    Класс профиля

    Инициализация -> дается список путей к файлам EDI, внутри это читается и сохраняется
    """

    def __init__(self, file_paths: list[str]):
        self.file_paths = file_paths
        self.edis = read_edi_files(self.file_paths)

        self.normalizations = {}
        self.data_widget = None
        self.create_widget()

    def add_normalization(self, period, mt_points):
        if len(self.normalizations):
            norm_id = list(self.normalizations.keys())[-1] + 1
        else:
            norm_id = 1
        norm = Normalization(self.edis, mt_points, period, file_paths=self.file_paths)
        self.normalizations[norm_id] = norm

        return norm

    def add_edi(self, file_path: list[str]):
        self.file_paths.extend(file_path)
        self.edis.extend(read_edi_files(file_path))

    def delete_edi(self, file_name):
        for file in self.file_paths:
            if file_name == os.path.basename(file):
                self.file_paths.remove(file)
                break

    def delete_normalization(self, norm_id):
        del self.normalizations[norm_id]

    def create_widget(self):
        self.data_widget = MTComponent(self, self.edis)
