import marimo

__generated_with = "0.11.26"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Adding a crystal
        This worksheet guides you through all the things you need to make an entry for the crystal database, then produces a text block you can paste inside.

        To get started, just some imports:
        """
    )
    return


@app.cell
def _():
    import numpy as np
    import LightwaveExplorer as lwe
    import marimo as mo
    return lwe, mo, np


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Next, let's write down the info we'll need to make the database entry
        This is a list of what you need to know about your crystal to describe it to LWE.

        This basic form will assume we can get our refractive index directly from refractiveindex.info; I'll also provide a more advanced example where we can fit to a different set of data.
        """
    )
    return


@app.cell
def _(mo):
    #First, the name of the crystal. I'm putting L afterwards to show that the refractive index is going to be stored in complex valued Lorentzians. A G means that the absorption features are Gaussians.
    name=mo.ui.text(value="My new crystal", label="Crystal name")
    name
    return (name,)


@app.cell
def _(mo):
    crystal_options = options=["Isotropic","Uniaxial","Biaxial"]
    crystal_type = mo.ui.dropdown(options=crystal_options,value="Isotropic", label="Crystal type")
    crystal_type
    return crystal_options, crystal_type, options


@app.cell
def _(mo):
    sellmeier_eqn = mo.ui.dropdown(options=["Multielement fit", "Lorentzians", "Gaussians"], value="Lorentzians", label="Target equation")
    sellmeier_eqn
    return (sellmeier_eqn,)


@app.cell
def _(mo):
    sellmeier_ref = mo.ui.text(value="Journal of ... 52, 1023 (2028).",label="Refractive index reference", full_width=True)
    sellmeier_ref
    return (sellmeier_ref,)


@app.cell
def _(mo):
    has_chi2 = mo.ui.checkbox(value=True,label="Has chi(2)")
    has_chi2
    return (has_chi2,)


@app.cell
def _(mo):
    chi2_tensor_text = mo.ui.text_area(value="0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0", label="d tensor")
    chi2_tensor_text
    return (chi2_tensor_text,)


@app.cell
def _(mo):
    d_tensor_ref = mo.ui.text(value="Journal of ... 52, 1023 (2028).",label="d tensor reference", full_width=True)
    d_tensor_ref
    return (d_tensor_ref,)


@app.cell
def _(chi2_tensor_text, np):
    chi2 = np.fromstring(chi2_tensor_text.value, sep=' ')
    if chi2.size == 18:
        chi2 = np.reshape(chi2,(3,6))
    print(chi2)
    return (chi2,)


@app.cell
def _(mo):
    freq_element_0 = mo.ui.text(value="375e12")
    freq_element_1 = mo.ui.text(value="750e12")
    nonlinear_frequencies = mo.ui.array([freq_element_0]*2 + [freq_element_1], label="chi(2) measurement frequencies")
    nonlinear_frequencies
    return freq_element_0, freq_element_1, nonlinear_frequencies


@app.cell
def _(mo):
    chi3_options = ["none", "full tensor", "single value (assume isotropic)", "single value from n_2 (assume isotropic)"]
    chi3_type = mo.ui.dropdown(options=chi3_options, value=chi3_options[2], label="chi(3) type")
    chi3_type
    return chi3_options, chi3_type


@app.cell
def _(chi3_options, chi3_type, mo):
    chi3 = 0.0
    if chi3_type.value == chi3_options[1]:
        chi3_entry = mo.ui.text_area(label="chi(3) tensor")
    elif chi3_type.value == chi3_options[2]:
        chi3_entry = mo.ui.text(label="chi(3)")
    elif chi3_type.value == chi3_options[3]:
        chi3_entry = mo.ui.array([mo.ui.text(label="n2"), mo.ui.text(label="n")],label="n2 and n")

    chi3_entry
    return chi3, chi3_entry


@app.cell
def _(mo):
    chi3_ref = mo.ui.text(value="Journal of ... 52, 1023 (2028).",label="Chi(3) reference", full_width=True)
    chi3_ref
    return (chi3_ref,)


@app.cell
def _(freq_element_0, mo):
    chi3_frequencies = mo.ui.array([freq_element_0]*4, label="chi(3) measurement frequencies")
    chi3_frequencies
    return (chi3_frequencies,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Generate a block of text that can be pasted inside the CrystalDatabase.txt file""")
    return


@app.cell
def _(
    chi2,
    chi3,
    chi3_frequencies,
    chi3_ref,
    chi3_type,
    d_tensor_ref,
    has_chi2,
    mo,
    name,
    nonlinear_frequencies,
    sellmeier_eqn,
    sellmeier_ref,
):
    with mo.capture_stdout() as buffer:
        print("Name:")
        print(name.value)
        print("Type:")
        print(chi3_type.value)
        print("Sellmeier equation:")
        print(sellmeier_eqn.value)
    
    
        print("Sellmeier reference:")
        print(sellmeier_ref.value)
        print("chi2 type:")
        print(has_chi2.value)
        print("d:")
        print('\n'.join(' '.join(map(str, row)) for row in chi2))
        print("d reference:")
        print(d_tensor_ref.value)
        print("chi3 type:")
        print(chi3_type.value)
        print("chi3:")
        if chi3_type.value == 1:
            print('\n'.join(' '.join(map(str, row)) for row in chi3))
        else:
            print(chi3)
            print(0)
            print(0)
        print("chi3 reference:")
        print(chi3_ref.value)
        print("Spectral file:")
        print("None")
        print("Nonlinear reference frequencies:")
        print(f"{float(nonlinear_frequencies[0].value):.8g} {float(nonlinear_frequencies[1].value):.8g} {float(nonlinear_frequencies[2].value):.8g} {float(chi3_frequencies[0].value):.8g} {float(chi3_frequencies[1].value):.8g} {float(chi3_frequencies[2].value):.8g} {float(chi3_frequencies[3].value):.8g}")
        print("~~~crystal end~~~\n")

    output_box = mo.ui.text_area(value=buffer.getvalue(),rows=35)
    output_box


    return buffer, output_box


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
