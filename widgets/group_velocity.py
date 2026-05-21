# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "attoworld==2026.1.2",
#     "lightwaveexplorer==2026.1",
#     "marimo>=0.23.6",
#     "matplotlib==3.10.9",
#     "numpy==2.4.6",
#     "scipy==1.17.1",
# ]
# [tool.marimo.display]
# theme = "dark"
# ///

import marimo

__generated_with = "0.23.6"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(r"""
    ## Pulse propagation toy
    The set of sliders below let you control the propagation of a pulse through fused silica. Try to adjust the parameters to get a feeling for how the pulse stretches and gets delayed, depending on its frequency and bandwidth, and on the thickness of the glass. Some of the relevant material properties, the dispersion curves, are plotted underneath it.

    All of the waves have equal amplitude, and the bandwidth is finite: what function describes the pulse envelope (for thickness 0)?

    Bonus question: why does the pulse repeat if the number of frequencies is low, or the length of the time window is too long?
    """)
    return


@app.cell
def _(mo):
    thickness_slider = mo.ui.slider(
        start=0.0,
        stop=300,
        step=1,
        value=160.0,
        label="Thickness (\u03bcm)",
        show_value=True,
    )
    frequency_slider = mo.ui.slider(
        start=101,
        stop=600,
        step=1,
        value=300,
        label="Frequency (THz)",
        show_value=True,
    )
    bandwidth_slider = mo.ui.slider(
        start=10,
        stop=200,
        step=5,
        value=120,
        label="Bandwidth (THz)",
        show_value=True,
    )
    time_length_slider = mo.ui.slider(
        start=10,
        stop=100,
        step=5,
        value=40,
        label="Time range (fs)",
        show_value=True,
    )
    N_frequencies_slider = mo.ui.slider(
        start=1,
        stop=32,
        step=1,
        value=12,
        label="Number of frequencies",
        show_value=True,
    )
    moving_reference_frame_checkbox = mo.ui.checkbox(
        value=True, label="Moving reference frame"
    )
    gauss_check = mo.ui.checkbox(value=False, label="Gaussian amplitudes")

    mo.output.append(thickness_slider)
    mo.output.append(frequency_slider)
    mo.output.append(bandwidth_slider)
    mo.output.append(time_length_slider)
    mo.output.append(N_frequencies_slider)
    mo.output.append(moving_reference_frame_checkbox)
    mo.output.append(gauss_check)
    return (
        N_frequencies_slider,
        bandwidth_slider,
        frequency_slider,
        gauss_check,
        moving_reference_frame_checkbox,
        thickness_slider,
        time_length_slider,
    )


@app.cell
def _(
    N_frequencies_slider,
    aw,
    bandwidth_slider,
    constants,
    frequency_slider,
    gauss_check,
    lwe,
    moving_reference_frame_checkbox,
    np,
    plt,
    sellmeier_coefficients,
    sig,
    thickness_slider,
    time_length_slider,
):
    def plot_group_velocity(
        t_length: float = time_length_slider.value * 1e-15, N_grid: int = 1024
    ):
        f_low = frequency_slider.value - 0.5 * bandwidth_slider.value
        f_high = f_low + bandwidth_slider.value
        freqs = np.linspace(
            f_low * 1e12,
            (f_low + bandwidth_slider.value) * 1e12,
            N_frequencies_slider.value,
        )
        freq_offsets = 2 * np.pi * (freqs - 1e12 * frequency_slider.value)

        central_freq = np.mean(freqs)
        lams = 1e6 * constants.speed_of_light / freqs
        ns = np.real(lwe.sellmeier(lams, sellmeier_coefficients, 0))
        if moving_reference_frame_checkbox.value:
            n0 = np.real(
                lwe.sellmeier(
                    np.array([1e6 * constants.speed_of_light / central_freq]),
                    sellmeier_coefficients,
                    0,
                )
            )
        else:
            n0 = 0.0

        _t = np.linspace(-t_length, t_length, N_grid)
        waves = np.zeros((_t.shape[0], freqs.shape[0]))
        t0 = 0
        spacing = 12 / N_frequencies_slider.value
        _colors = []

        def color_curve(f, f0, sigma):
            return np.exp(-((f0 - f) ** 2) / (2 * sigma**2))

        for _i in range(freqs.shape[0]):
            _colors.append(
                (
                    color_curve(
                        f=freqs[_i] / 1e12, f0=f_low, sigma=0.5 * bandwidth_slider.value
                    ),
                    color_curve(
                        f=freqs[_i] / 1e12,
                        f0=(f_low + f_high) / 2,
                        sigma=0.3 * bandwidth_slider.value,
                    ),
                    color_curve(
                        f=freqs[_i] / 1e12,
                        f0=f_high,
                        sigma=0.5 * bandwidth_slider.value,
                    ),
                )
            )
            _phase = (
                thickness_slider.value
                * 1e-6
                * 2
                * np.pi
                * freqs[_i]
                * (ns[_i] - n0)
                / constants.speed_of_light
            )
            if gauss_check.value:
                _amplitude = np.exp(
                    -(freq_offsets[_i] ** 2) / (4e24 * bandwidth_slider.value**2)
                )
            else:
                _amplitude = 1
            waves[:, _i] = _amplitude * np.cos(2 * np.pi * freqs[_i] * _t - _phase)
            (_lineplot,) = plt.plot(
                1e15 * _t, waves[:, _i] + spacing * _i, color=_colors[_i]
            )

        field = 2 * np.sum(waves, axis=1) / freqs.shape[0]
        plt.plot(1e15 * _t, field - 2 * spacing, color="magenta")
        plt.plot(1e15 * _t, np.abs(sig.hilbert(field)) - 2 * spacing, color="gray")

        plt.xlim(-1e15 * t_length / 2, 1e15 * t_length / 2)
        plt.xlabel("Time (fs)")
        plt.yticks([])
        return plt.gcf()

    plot_group_velocity()
    aw.plot.showmo()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Material properties

    Here are some of the relevant material properties for fused silica, as defined in the lecture notes.

    First, we'll plot the refractive index, and the group index, followed by the GVD.

    You can change the material by replacing the Sellmeier coefficients below, in the format used by Lightwave Explorer (i.e. copy-and-paste from CrystalDatabase.txt).
    """)
    return


@app.cell
def _(mo):
    sellmeier_input = mo.ui.text_area(
        label="Sellmeier coefficients:",
        value="1 0 0.6961663 -0.00467914825849 0 0.4079426 -0.013512063074 0 0.8974794 -97.9340025379 0 0 1 0 0 0 0 0 0 0 0 0",
        full_width=True,
        rows=1,
    )
    mo.output.append(sellmeier_input)
    eqn_type_input = mo.ui.number(start=0, stop=2, value=0, label="Equation type:")
    mo.output.append(eqn_type_input)
    return eqn_type_input, sellmeier_input


@app.cell
def _(
    aw,
    constants,
    eqn_type,
    get_group_index,
    get_gvd,
    lwe,
    np,
    plt,
    sellmeier_coefficients,
):
    freq = np.linspace(150e12, 1000e12, 1024)
    wavelength = 1e6 * constants.speed_of_light / freq
    group_index = np.zeros(freq.shape)
    gvd = np.zeros(freq.shape)
    for _i in range(freq.shape[0]):
        group_index[_i] = get_group_index(
            sellmeier_coefficients, eqn_type, frequency=freq[_i]
        )
        gvd[_i] = get_gvd(sellmeier_coefficients, eqn_type, frequency=freq[_i])
    plt.plot(
        wavelength,
        np.real(lwe.sellmeier(wavelength, sellmeier_coefficients, eqn_type)),
        label="Refractive index",
    )
    plt.plot(wavelength, group_index, label="Group index")
    plt.xlabel("Wavelength (\u03bcm)")
    plt.ylabel("Index")
    plt.legend()
    aw.plot.showmo()
    return gvd, wavelength


@app.cell
def _(aw, gvd, plt, wavelength):
    plt.plot(wavelength, 1e27 * gvd)
    plt.xlabel("Wavelength (\u03bcm)")
    plt.ylabel("GVD (fs$^2$/mm)")
    aw.plot.showmo()
    return


@app.cell
def _(eqn_type_input, np, sellmeier_input):
    sellmeier_coefficients = np.fromstring(sellmeier_input.value, sep=" ")
    eqn_type = int(eqn_type_input.value)
    return eqn_type, sellmeier_coefficients


@app.cell
def _(aw, constants, lwe, np):
    def get_gvd(
        sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float = 375e12
    ):
        positions = np.array([-3, -2, -1, 0, 1, 2, 3])
        df = 1e-3 * frequency
        frequencies = frequency + df * positions
        wavelengths = 1e6 * constants.speed_of_light / frequencies
        stencil = aw.numeric.fornberg_stencil(
            2, np.array([-3, -2, -1, 0, 1, 2, 3], dtype=float)
        )
        n = np.real(lwe.sellmeier(wavelengths, sellmeier_coeffs, equation_type))
        k = n * 2 * np.pi * frequencies / constants.speed_of_light
        d2kdw2 = np.sum(stencil * k) / (4 * (np.pi * df) ** 2)
        return d2kdw2

    def get_group_velocity(
        sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float = 375e12
    ):
        positions = np.array([-3, -2, -1, 0, 1, 2, 3])
        df = 1e-3 * frequency
        frequencies = frequency + df * positions
        wavelengths = 1e6 * constants.speed_of_light / frequencies
        stencil = aw.numeric.fornberg_stencil(
            1, np.array([-3, -2, -1, 0, 1, 2, 3], dtype=float)
        )
        n = np.real(lwe.sellmeier(wavelengths, sellmeier_coeffs, equation_type))
        k = n * 2 * np.pi * frequencies / constants.speed_of_light
        dkdw = np.sum(stencil * k) / (2 * np.pi * df)
        return 1.0 / dkdw

    def get_group_index(
        sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float = 375e12
    ):
        return constants.speed_of_light / get_group_velocity(
            sellmeier_coeffs, equation_type, frequency
        )

    return get_group_index, get_gvd


@app.cell
async def _():
    import sys

    import LightwaveExplorer as lwe
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np

    is_in_web_notebook = sys.platform == "emscripten"
    if is_in_web_notebook:
        import zipfile

        import micropip

        await micropip.install(
            "https://nickkarpowicz.github.io/wheels/attoworld-2026.1.2-cp312-cp312-emscripten_3_1_58_wasm32.whl"
        )
    import attoworld as aw

    aw.plot.set_style("nick_dark", font_size=14)
    from scipy import constants
    from scipy import signal as sig

    return aw, constants, lwe, mo, np, plt, sig


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
