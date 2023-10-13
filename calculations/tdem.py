import numpy as np
import SimPEG.electromagnetics.time_domain as tdem
from SimPEG import maps


def forward_tdem_1d(rho_array, h_array, times, I, loop_radius, h=1):
    """
    :param rho_array: массив сопротивлений
    :param h_array: массив мощностей слоев
    :param times: массив времен становлений
    :param I: сила тока
    :param loop_radius: радиус петли
    :param h: высота петли над поверхностью
    """

    mu0 = 4 * np.pi * 1e-7
    time_eps = 1e-9  # сдвиг для расчета производной магнитной индукции

    # simpeg периодически криво считает dB/dT, поэтому производную считать будем вручную
    times_for_calc = np.empty(len(times) * 2)
    times_for_calc[::2] = times - time_eps / 2
    times_for_calc[1::2] = times + time_eps / 2

    sigma_array = 1 / rho_array

    source_location = np.array([0.0, 0.0, h])
    source_orientation = "z"

    receiver_location = np.array([0.0, 0.0, h])
    receiver_orientation = "z"

    receiver_list = [
        tdem.receivers.PointMagneticFluxDensity(
            receiver_location, times_for_calc, orientation=receiver_orientation
        )
    ]

    source_list = [
        tdem.sources.CircularLoop(
            receiver_list=receiver_list,
            location=source_location,
            waveform=tdem.sources.StepOffWaveform(0),
            current=I,
            radius=loop_radius,
        )
    ]

    survey = tdem.Survey(source_list)

    simulation = tdem.Simulation1DLayered(
        survey=survey,
        thicknesses=h_array,
        sigmaMap=maps.IdentityMap(nP=len(sigma_array)),
    )

    B = simulation.dpred(sigma_array)
    dBdT = -(B[1::2] - B[::2]) / (times_for_calc[1::2] - times_for_calc[::2])

    S = np.pi * loop_radius ** 2
    dUq = S * dBdT

    rho_tau = mu0 / np.pi / times * (S * S * mu0 * I / 20 / times / dUq) ** (
                2 / 3)

    return rho_tau


def export_tdem(file_name: str, tdem_data: tuple):
    ro_t, kernel_abs, w_list, kernel_ifft_abs, T, dHdt_V = tdem_data
    result = ['\t'.join(['Ro(t)', 'Efi(w)', 'w', 'Efi(t)', 't(c)', 'dH/dt(В)'])]
    for i in range(len(w_list)):
        result.append([ro_t[i], kernel_abs[i], w_list[i], kernel_ifft_abs[i], T[i], dHdt_V[i]])

    for i in range(1, len(result)):
        result[i] = [str(x) for x in result[i]]
        result[i] = '\t'.join(result[i])

    result_text = '\n'.join(result)
    with open(file_name, 'w') as file:
        file.write(result_text)
    return

