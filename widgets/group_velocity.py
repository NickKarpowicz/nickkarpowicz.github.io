import marimo

__generated_with = "0.13.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Pulse propagation toy
    The set of sliders below let you control the propagation of a pulse through fused silica. Try to adjust the parameters to get a feeling for how the pulse stretches and gets delayed, depending on its frequency and bandwidth, and on the thickness of the glass. Some of the relevant material properties, the dispersion curves, are plotted underneath it.

    All of the waves have equal amplitude, and the bandwidth is finite: what function describes the pulse envelope (for thickness 0)?

    Bonus question: why does the pulse repeat if the number of frequencies is low, or the length of the time window is too long?
    """
    )
    return


@app.cell
def _(mo):
    thickness_slider = mo.ui.slider(start=0.0,stop=300,step=5,value=160.0, label="Thickness (\u03bcm)", show_value=True)
    frequency_slider = mo.ui.slider(start=101,stop=600,step=1,value=300, label="Frequency (THz)", show_value=True)
    bandwidth_slider = mo.ui.slider(start=10,stop=200, step=5, value=120, label="Bandwidth (THz)", show_value=True)
    time_length_slider = mo.ui.slider(start=10, stop=100, step=5, value=40, label="Time range (fs)", show_value=True)
    N_frequencies_slider = mo.ui.slider(start=3, stop=20, step=1, value=12, label="Number of frequencies", show_value=True)
    mo.output.append(thickness_slider)
    mo.output.append(frequency_slider)
    mo.output.append(bandwidth_slider)
    mo.output.append(time_length_slider)
    mo.output.append(N_frequencies_slider)
    return (
        N_frequencies_slider,
        bandwidth_slider,
        frequency_slider,
        thickness_slider,
        time_length_slider,
    )


@app.cell
def _(
    N_frequencies_slider,
    bandwidth_slider,
    constants,
    frequency_slider,
    get_group_index,
    lwe,
    np,
    plt,
    sellmeier_coefficients,
    showmo,
    sig,
    thickness_slider,
    time_length_slider,
):
    def plot_group_velocity(thickness: float=0.0, t_length: float = time_length_slider.value*1e-15, N_grid: int=1024):
        f_low = frequency_slider.value - 0.5*bandwidth_slider.value
        f_high = (f_low+bandwidth_slider.value)
        freqs = np.linspace(f_low*1e12, (f_low+bandwidth_slider.value)*1e12, N_frequencies_slider.value)
        lams = 1e6*constants.speed_of_light/freqs
        ns = np.real(lwe.sellmeier(lams,sellmeier_coefficients,0))
        _t = np.linspace(-t_length,t_length,N_grid)
        waves = np.zeros((_t.shape[0],freqs.shape[0]))
        t0 = ns[0] * thickness / constants.speed_of_light
        spacing = 12/N_frequencies_slider.value
        _colors=[]
        def color_curve(f,f0,sigma):
            return np.exp(-(f0-f)**2/(2*sigma**2))
        for _i in range(freqs.shape[0]):
            _colors.append((
                color_curve(f=freqs[_i]/1e12,f0=1.05*f_low,sigma=0.4*bandwidth_slider.value),
                color_curve(f=freqs[_i]/1e12,f0=(f_low+f_high)/2,sigma=0.35*bandwidth_slider.value),
                color_curve(f=freqs[_i]/1e12,f0=0.95*f_high,sigma=0.4*bandwidth_slider.value)))
            _phase = -thickness * 2*np.pi*freqs[_i] * (ns[_i]-ns[0])/constants.speed_of_light
            waves[:,_i] = np.cos(2*np.pi*freqs[_i]*_t + _phase)
            _lineplot, = plt.plot(1e15*_t,waves[:,_i]+spacing*_i,color=_colors[_i],linewidth=1.25)
            _markerplot, = plt.plot(1e15*(ns[_i]-ns[0]) * thickness/constants.speed_of_light, spacing*_i + 1,'d',color=_colors[_i])
        n_group = get_group_index(sellmeier_coefficients,0,frequency=np.mean(freqs))
        plt.plot(1e15*_t,2*np.sum(waves,axis=1)/freqs.shape[0] - 2*spacing,color='black', linewidth=2)
        plt.plot(1e15*_t,2*np.abs(sig.hilbert(np.sum(waves,axis=1)))/freqs.shape[0] - 2*spacing,color='gray')
        plt.plot(1e15*(n_group-ns[0]) * thickness/constants.speed_of_light,-2*spacing + 2,'d',color='black')
        plt.xlim(-1e15*t_length/2,1e15*t_length/2)
        plt.xlabel("Time (fs)")
        plt.yticks([])
        return plt.gcf()

    plot_group_velocity(thickness=1e-6*thickness_slider.value)
    showmo()
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    ## Material properties

    Here are some of the relevant material properties for fused silica, as defined in the lecture notes.

    First, we'll plot the refractive index, and the group index, followed by the GVD.

    You can change the material by replacing the Sellmeier coefficients below, in the format used by Lightwave Explorer (i.e. copy-and-paste from CrystalDatabase.txt).
    """
    )
    return


@app.cell
def _(mo):
    sellmeier_input = mo.ui.text_area(label="Sellmeier coefficients:", value="1 0 0.6961663 -0.00467914825849 0 0.4079426 -0.013512063074 0 0.8974794 -97.9340025379 0 0 1 0 0 0 0 0 0 0 0 0",full_width=True, rows=1)
    mo.output.append(sellmeier_input)
    eqn_type_input = mo.ui.number(start=0,stop=2,value=0, label="Equation type:")
    mo.output.append(eqn_type_input)
    return eqn_type_input, sellmeier_input


@app.cell
def _(
    constants,
    eqn_type,
    get_group_index,
    get_gvd,
    lwe,
    np,
    plt,
    sellmeier_coefficients,
    showmo,
):
    freq = np.linspace(150e12, 1000e12, 1024)
    wavelength = 1e6*constants.speed_of_light/freq
    group_index = np.zeros(freq.shape)
    gvd = np.zeros(freq.shape)
    for _i in range(freq.shape[0]):
        group_index[_i] = get_group_index(sellmeier_coefficients,eqn_type,frequency=freq[_i])
        gvd[_i] = get_gvd(sellmeier_coefficients,eqn_type,frequency=freq[_i])
    plt.plot(wavelength,np.real(lwe.sellmeier(wavelength, sellmeier_coefficients, eqn_type)), label="Refractive index")
    plt.plot(wavelength, group_index, label="Group index")
    plt.xlabel("Wavelength (\u03bcm)")
    plt.ylabel("Index")
    plt.legend()
    #plt.savefig("group_index.pdf",bbox_inches='tight')
    showmo()
    return gvd, wavelength


@app.cell
def _(gvd, plt, showmo, wavelength):
    plt.plot(wavelength, 1e27*gvd)
    plt.xlabel("Wavelength (\u03bcm)")
    plt.ylabel("GVD (fs$^2$/mm)")
    #plt.savefig("gdd.pdf",bbox_inches='tight')
    showmo()
    return


@app.cell
def _():
    return


@app.cell
def _(constants, fornberg_stencil, lwe, np):
    def get_gvd(sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float=375e12):
        positions = np.array([-3, -2, -1, 0, 1, 2, 3])
        df = 1e-3 * frequency
        frequencies = frequency + df * positions
        wavelengths = 1e6 * constants.speed_of_light/frequencies
        stencil = fornberg_stencil(2, positions)
        n = np.real(lwe.sellmeier(wavelengths, sellmeier_coeffs, equation_type))
        k = n * 2 * np.pi * frequencies / constants.speed_of_light
        d2kdw2 = np.sum(stencil * k)/(4 * (np.pi * df)**2)
        return d2kdw2

    def get_group_velocity(sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float=375e12):
        positions = np.array([-3, -2, -1, 0, 1, 2, 3])
        df = 1e-3 * frequency
        frequencies = frequency + df * positions
        wavelengths = 1e6 * constants.speed_of_light/frequencies
        stencil = fornberg_stencil(1, positions)
        n = np.real(lwe.sellmeier(wavelengths, sellmeier_coeffs, equation_type))
        k = n * 2 * np.pi * frequencies / constants.speed_of_light
        dkdw = np.sum(stencil * k)/(2 * np.pi * df)
        return 1.0/dkdw

    def get_group_index(sellmeier_coeffs: np.ndarray, equation_type: int, frequency: float=375e12):
        return constants.speed_of_light/get_group_velocity(sellmeier_coeffs, equation_type, frequency)
    return get_group_index, get_gvd


@app.cell
def _():

    #plot_group_velocity(thickness=0e-6).savefig("group_velocity_0.pdf",bbox_inches='tight')
    return


@app.cell
def _():
    #plot_group_velocity(thickness=160e-6).savefig("group_velocity_1.pdf",bbox_inches='tight')
    return


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import LightwaveExplorer as lwe
    from matplotlib import rcParams, cycler
    import matplotlib.pyplot as plt
    from scipy import constants
    import scipy.signal as sig
    import io
    def light_plot():
        """
        Use a light style for matplotlib plots.
        """
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Verdana', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']
        rcParams['axes.prop_cycle'] = cycler(color=["blue", "magenta", "orange", "purple", "indigo", "cyan", "black"])
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Verdana', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']
    light_plot()

    def fornberg_stencil(order: int, positions: np.ndarray, position_out: float = 0.0) -> np.ndarray:
        """
        Generate a finite difference stencil using the algorithm described by B. Fornberg
        in Mathematics of Computation 51, 699-706 (1988).

        Args:
            order (int): the order of the derivative
            positions (np.ndarray): the positions at which the functions will be evaluated in the stencil.
        Returns:
            np.ndarray: the finite difference stencil with weights corresponding to the positions in the positions input array

        Examples:

            >>> stencil = fornberg_stencil(1, [-1,0,1])
            >>> print(stencil)
            [-0.5 0. 0.5]
        """
        # Contributed by Nick Karpowicz
        number_of_positions = len(positions)
        delta = np.zeros((number_of_positions, number_of_positions, order + 1), dtype=float)
        delta[0, 0, 0] = 1.0
        c1 = 1.0

        for n in range(1, number_of_positions):
            c2 = 1.0
            for v in range(n):  # v from 0 to n-1
                c3 = positions[n] - positions[v]
                c2 *= c3
                if n <= order:
                    delta[n-1, v, n] = 0.0
                for m in range(min(n, order) + 1):
                    if m == 0:
                        last_element = 0.0
                    else:
                        last_element = m * delta[n-1, v, m-1]
                    delta[n, v, m] = ((positions[n] - position_out) * delta[n-1, v, m] - last_element) / c3
            for m in range(min(n, order) + 1):
                if m == 0:
                    first_element = 0.0
                else:
                    first_element = m * delta[n-1, n-1, m-1]
                delta[n, n, m] = (c1 / c2) * (first_element - (positions[n-1] - position_out) * delta[n-1, n-1, m])
            c1 = c2

        return delta[-1, :, -1].squeeze()

    def showmo():
        """
        Helper function to plot as an svg to have vector plots in marimo notebooks
        """
        # Contributed by Nick Karpowicz
        svg_buffer = io.StringIO()
        plt.savefig(svg_buffer, format='svg')
        return mo.Html(svg_buffer.getvalue())
    return constants, fornberg_stencil, lwe, mo, np, plt, showmo, sig


@app.cell
def _(eqn_type_input, np, sellmeier_input):
    sellmeier_coefficients = np.fromstring(sellmeier_input.value, sep=" ")
    eqn_type = int(eqn_type_input.value)
    return eqn_type, sellmeier_coefficients


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
