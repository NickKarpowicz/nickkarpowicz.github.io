import marimo
# /// script
# [tool.marimo.display]
# theme = "dark"
# ///
__generated_with = "0.11.24"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import LightwaveExplorer as lwe
    import io
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rcParams, cycler
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Verdana', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']
    rcParams['axes.prop_cycle'] = cycler(color=["cyan", "magenta", "orange", "purple"]) 
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'Verdana', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']
    rcParams['figure.facecolor'] = 'black'
    rcParams['figure.edgecolor'] = 'black'
    rcParams['savefig.facecolor'] = 'black'
    rcParams['savefig.edgecolor'] = 'black'
    rcParams['axes.facecolor'] = 'black'
    rcParams['text.color'] = 'white'
    rcParams['axes.edgecolor'] = 'white'
    rcParams['axes.labelcolor'] = 'white'
    rcParams['xtick.color'] = 'white'
    rcParams['ytick.color'] = 'white'
    rcParams['grid.color'] = 'white'
    rcParams['lines.color'] = 'white'


    def showmo():
        """
        Helper function to plot as an svg and have it display in marimo in vector form
        """
        svg_buffer = io.StringIO()
        plt.savefig(svg_buffer, format='svg')
        svg_buffer.seek(0)
        svg_data = svg_buffer.getvalue()
        return mo.Html(svg_data)
    return cycler, io, lwe, mo, np, plt, rcParams, showmo


@app.cell
def _(mo):
    f0 = mo.ui.slider(start=50, stop=1000, step=10, value=400, label="Frequency (THz)", show_value=True)
    f0
    return (f0,)


@app.cell
def _(mo):
    bandwidth = mo.ui.slider(start=2, stop=300, value=80, label="Bandwidth (THz)", show_value=True)
    bandwidth
    return (bandwidth,)


@app.cell
def _(mo):
    Ns = mo.ui.slider(start=2, stop=16, step=2, value=2, label="SG order", show_value=True)
    Ns
    return (Ns,)


@app.cell
def _(mo):
    cep = mo.ui.slider(start=0.0, stop=2.0, step=0.1, value=0, label="CEP/pi", show_value=True)
    cep
    return (cep,)


@app.cell
def _(mo):
    tau = mo.ui.slider(start=-20, stop=20, step=1, value=0, label="Delay (fs)", show_value=True)
    tau
    return (tau,)


@app.cell
def _(mo):
    phi2 = mo.ui.slider(start=-100, stop=100, step=1, value=0, label="GDD (fs^2)", show_value=True)
    phi2
    return (phi2,)


@app.cell
def _(mo):
    phi3 = mo.ui.slider(start=-300, stop=300, step=1, value=0, label="TOD (fs^3)", show_value=True)
    phi3
    return (phi3,)


@app.cell
def _(mo):
    grid_length = mo.ui.slider(start=10, stop=1000, step=1, value=160, label="Grid length (fs)", show_value=True)
    grid_length
    return (grid_length,)


@app.cell
def _(mo):
    grid_dt = mo.ui.slider(start=0.05, stop = 1.0, step=0.01, value=0.25, label="dt (fs)", show_value=True)
    grid_dt
    return (grid_dt,)


@app.cell
def _(grid_dt, grid_length, np):
    dt = grid_dt.value * 1e-15
    Nt = int(grid_length.value*1e-15/dt)
    t = np.linspace(-0.5e-15*grid_length.value,0.5e-15*grid_length.value,Nt)
    w = 2*np.pi*np.fft.fftfreq(Nt,dt)
    return Nt, dt, t, w


@app.cell
def _(
    Ns,
    bandwidth,
    cep,
    f0,
    grid_length,
    lwe,
    np,
    phi2,
    phi3,
    plt,
    showmo,
    t,
    tau,
    w,
):

    ws = (w-2*np.pi*f0.value*1e12)
    phi = np.pi*cep.value + (-0.5e-15*grid_length.value + tau.value*1e-15)*w + 0.5*phi2.value*1e-30*ws**2 + (1.0/6)*phi3.value*1e-45*ws**3
    Ew = np.exp(-ws**Ns.value/(2*np.pi*bandwidth.value*1e12)**Ns.value -1.0j * phi)
    Ew[w<0]=0
    Et = lwe.norma(np.fft.ifft(Ew))
    fwhm_value = 1e15*lwe.fwhm(t,np.abs(Et)**2)
    plt.plot(1e15*t,np.real(Et))
    plt.xlabel("Time (fs)")
    plt.title(f"FWHM duration: {fwhm_value: .2f} fs")
    showmo()
    return Et, Ew, fwhm_value, phi, ws


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
