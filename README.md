# SpecView

## Overview
SpecView is an interactive dashboard for viewing and analyzing spectral images, such as those collected by cathodoluminescence or EELS techniques. It is runnable in a Jupyter Notebook, which allows easy access to the underlying data. For example, SpecView enables the user to collect an average spectra over a specific ROI by drawing a box. The resulting spectra is accessable as a class attribute for further in-notebook manipulation.

## Installation
The following instructions assume that Anaconda and Jupyter Lab have already been installed. Jupyter Notebook likely works too, but the code has not been tested there. Other environment managers besides Anaconda should be fine too, as long as they can create an environment from a YML file. The YML file contains the package dependencies.

Here are step-by-step instructions for installing SpecView:

1. Clone this repo, or alternatively download it as a zip and extract.
2. Open up Anaconda Prompt and navigate to the specview directory.
3. Create a new environment by entering: `conda env create --file environment.yml`
4. Activate the environment: `conda activate specview`
5. Open JupyterLab by entering `jupyter lab`.
6. Open `Example.ipynb` and explore.


