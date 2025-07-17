import marimo

__generated_with = "0.14.11"
app = marimo.App()


@app.cell
def _(mo):
    mo.md(
        r"""
    # Solving the TDSE of a 2-level system
    This is a very simple PDE to solve, but it gives you an impression of how the CN scheme works, and can be expanded to more levels, or using a different basis (e.g. space or momentum), if you change the dipole operator. Bonus question - *what should you change the dipole matrix to, if you want to represent the wavefunction in space?*
    """
    )
    return


@app.cell
def _():
    import numpy as np
    import matplotlib.pyplot as plt
    import marimo as mo
    import io
    import marimo as mo
    from cycler import cycler

    # I'm copying a few functions from the attoworld module so I don't have to include it, since there's no wasm version
    font_size=12
    plt.rcParams.update(
        {
            "font.size": font_size,
            "xtick.labelsize": 0.9 * font_size,
            "ytick.labelsize": 0.9 * font_size,
            "legend.fontsize": 0.9 * font_size,
            "figure.autolayout": True,
        }
    )
    plt.rcParams.update(
        {
            "font.sans-serif": [
                "Helvetica",
                "Nimbus Sans L",
                "Arial",
                "Verdana",
                "Nimbus Sans L",
                "DejaVu Sans",
                "Liberation Sans",
                "Bitstream Vera Sans",
                "sans-serif",
            ],
            "font.family": "sans-serif",
            "axes.prop_cycle": cycler(
                color=["cyan", "magenta", "orange", "blueviolet", "lime"]
            ),
            "figure.facecolor": "#171717",
            "figure.edgecolor": "#171717",
            "savefig.facecolor": "#171717",
            "savefig.edgecolor": "#171717",
            "axes.facecolor": "black",
            "text.color": "white",
            "axes.edgecolor": "white",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "ytick.color": "white",
            "grid.color": "white",
            "lines.color": "white",
            "image.cmap": "magma",
        }
    )
    def showmo():
        """Display a plot as an svg in marimo notebooks."""
        svg_buffer = io.StringIO()
        plt.savefig(svg_buffer, format="svg")
        return mo.output.append(mo.Html(svg_buffer.getvalue()))
    return mo, np, plt, showmo


@app.cell
def _(mo):
    mo.output.append(mo.md("## Simulation parameter sliders (atomic units)"))
    time_slider = mo.ui.slider(start=100, stop=3000, value=1600, step=100, label="time range", show_value=True)
    mo.output.append(time_slider)

    dt_slider = mo.ui.slider(start=0.01, stop=1, value=0.2, step = 0.01, label="dt", show_value=True)
    mo.output.append(dt_slider)

    frequency_slider=mo.ui.slider(start=0.0, stop=1.0, value=0.5, step = 0.01, label="omega", show_value=True)
    mo.output.append(frequency_slider)

    pulse_duration_slider = mo.ui.slider(start=10.0, stop = 1000.0, value=200.0, step = 10.0, label="pulse duration", show_value=True)
    mo.output.append(pulse_duration_slider)

    field_strength_slider = mo.ui.slider(start=0.0, stop=1.0, value=0.01, step=0.001, label="field strength", show_value=True)
    mo.output.append(field_strength_slider)

    dipole_moment_slider = mo.ui.slider(start=0.0, stop=1.0, value = 0.1, step=0.01, label="dipole element", show_value=True)
    mo.output.append(dipole_moment_slider)

    energy_difference_slider = mo.ui.slider(start=0.0, stop=2.0, value=0.5, step=0.01, label="energy difference", show_value=True)
    mo.output.append(energy_difference_slider)
    return (
        dipole_moment_slider,
        dt_slider,
        energy_difference_slider,
        field_strength_slider,
        frequency_slider,
        pulse_duration_slider,
        time_slider,
    )


@app.cell
def _(dipole_moment_slider, energy_difference_slider, np):
    #express energy and dipole operators (parts of the Hamiltonian) as 2x2 matrices
    #note that I'm using atomic units throughout, so hbar=1, and e=1
    energy = np.array([
        [0.1,0],
        [0,0.1+energy_difference_slider.value]])
    dipole = np.array([
        [0,dipole_moment_slider.value],
        [dipole_moment_slider.value,0]])

    #we'll use the identity matrix in writing out the Crank-Nicolson propagator
    Identity = np.eye(2)

    #Perform a step using the Crank-Nicolson approach
    #h is the time step
    #Note that this isn't the most efficient way to do this, but I hope it's transparent
    #enough to understand
    def cnStep(psi, h, initialField, nextField):
        #we express the CN scheme like this:
        #Ubackward psi_{n+1} = Uforward psi_n
        Uforward = Identity + (1j*h/2) * (energy + initialField*dipole)
        Ubackward = Identity - (1j*h/2) * (energy + nextField*dipole)

        #We solve the linear system of equations to get the new value of psi
        return np.linalg.solve(Ubackward, Uforward@psi)

    #polarization is the expectation values of the dipole operator
    def get_polarization(psi):
        return -np.real(np.dot(np.conj(psi.T),dipole@psi))
    return cnStep, get_polarization


@app.cell
def _(cnStep, get_polarization, np, plt, showmo):
    def generate_field(Nsteps, dt, maxField, duration, frequency):
        t = dt * np.linspace(1, Nsteps + 1, Nsteps + 1, dtype=float)
        t -= np.mean(t)
        pulse = (
            maxField 
            * np.exp(-(t**2) / (2 * duration**2)) 
            * np.cos(frequency * t))
        return t, pulse


    def run_simulation(Nsteps, dt, pulse):
        # initial state: everything in the low energy level
        psi0 = np.array([1.0, 0], dtype=complex)

        # we'll store the history of psi, as well as the polarization:
        psi = np.zeros((2, Nsteps), dtype=complex)
        polarization = np.zeros(Nsteps, dtype=float)

        for i in range(Nsteps):
            # calculate the polarization
            polarization[i] = get_polarization(psi0)

            # store the state of the wavefunction
            psi[:, i] = psi0

            # propagate to the next time step
            psi0 = cnStep(psi0, dt, pulse[i], pulse[i + 1])

        return polarization, psi


    def plot_simulation(
        Nsteps: int,
        dt: float,
        field_strength: float,
        pulse_duration: float,
        frequency: float,
    ):
        t, pulse = generate_field(
            Nsteps, dt, field_strength, pulse_duration, frequency
        )
        polarization, psi = run_simulation(Nsteps, dt, pulse)

        tplot = t[0:Nsteps]
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(tplot, polarization, label="Field")
        ax[0].plot(tplot, pulse[0:Nsteps], label="Polarization")
        ax[0].set_ylabel("Amplitude")
        ax[0].set_xlabel("t")
        ax[0].legend()

        ax[1].plot(tplot, np.abs(psi[0, :]) ** 2, label="Ground level")
        ax[1].plot(tplot, np.abs(psi[1, :]) ** 2, label="Excited level")
        ax[1].set_xlabel("t")
        ax[1].set_ylabel("Probability")
        ax[1].legend()
        showmo()
    return (plot_simulation,)


@app.cell
def _(
    dt_slider,
    field_strength_slider,
    frequency_slider,
    plot_simulation,
    pulse_duration_slider,
    time_slider,
):
    plot_simulation(
        Nsteps=int(time_slider.value / dt_slider.value),
        dt=dt_slider.value,
        field_strength=field_strength_slider.value,
        pulse_duration=pulse_duration_slider.value,
        frequency=frequency_slider.value,
    )
    return


if __name__ == "__main__":
    app.run()
