# Time_Spectral_Granger_Causality_EEG
Granger Causality (GC) pipeline for evaluating channels predicitability from data from [Xu et al 2023] described here https://www.pnas.org/doi/abs/10.1073/pnas.2216268120. For downloading the data refer to the Zenodo link here and download the entire [here](https://zenodo.org/record/7803212#.ZC3Cb-zML0q) and download the entire dataset. For any inquiry about data request or extra detail in code execution please don't hesitate to reach professor Jimo Borjigin, PhD [here](mailto:borjigin@umich.edu).

## Requirements

Before running be sure you have installed **Matlab** version **R2025a** or more recent as well as

1. Current **EEGlab** version downloaded from the SCCN [download website](https://sccn.ucsd.edu/eeglab/download.php)
2. **MVGC v1.0 version** from the [stable Github link](https://github.com/lcbarnett/MVGC1) supporting the GC calculation used in this repo.

Before running any command from Matlab prompt localized in your local machine using **cd** or from Matlab file **navigation bar**. Follow this command:

```matlab
cd /path/to/Time_Spectral_Granger_Causality_EEG
```

## 1.Run EEG preprocessing

For execution a Matlab preprocessing follow this command:

```matlab
time_frequency_gc_single_edf('<edf_filename_local>', '<output_file_suffix>')
```

For instance it is possible to evaluate the preprocessing with average rerefence using this command

```matlab
time_frequency_gc_single_edf('JF_20250225', '_average')
```

Or with Fz rereference with the following command

```matlab
time_frequency_gc_single_edf('JF_20250225', '_Fz_reref')
```

**Depending on the throughput quality and the amount of cores your of your machine processor this preprocessing can take 15-20mins for each patient .edf file**

The function **time_frequency_gc_single_edf** performs:
 
  1) **Data Loading**
     
       Supports large EDF files via chunked reading, or direct loading using biosig depending on the size of the file.
     
  2) **Re-referencing**
     
       Average reference or Fz reference.
     
  3) **Resampling**
     
       Downsample to 256 Hz.
  
  4) **Filtering**
     
    a) Notch filter at 60 Hz.
    b) Bandpass filter (0.1–100 Hz).
    
  5)  **Artifact Reduction**
     
    a) ASR-based cleaning (light EMG suppression).
    b) ICA decomposition (runica).
    c) Automatic IC rejection (heuristic-based) adding light severity suppressing for avoiding GC matrices singularity.
    
  6) **ECG Artifact Removal**
     
     Lagged linear regression using ECG reference channel
     
  7) **Channel Localization**
     
     Standard 10–20 montage assignment
     
  8) **GC Analysis** (optional) **can be evaluated later**
     
     Time-domain and spectral GC via **MVGC toolbox**

After executing this preprocessing command you will observe the following **.mat** files. The filename string will contain the suffix you added as input before the .mat extension.
In this case the suffix string is empty. This for the sake of comparison between different preprocessing options and for facilitating the exploration of resulting files across the **result_plots**, **preprocessed_save**, and **GC_estimation** folders.

<img width="709" height="398" alt="image" src="https://github.com/user-attachments/assets/fa0e2606-ddea-4234-88b7-626b82c7ffeb" />


Each file contains an EEGlab structure for each Near-Dead-Evennt (NDE) stages described in [Xu et al 2023] between **S1-S11** stages: having S1 as baseline stage before removing the ventilator, and from S2-S11 all the subsequent stages without the ventilator, switching multiple time the peacemaker activation. The file denoted with the preffix **JF_20250225** contains information with all the stages **S1-S11** separated for an adeaquate load using [**edfread.m**](https://www.mathworks.com/help/signal/ref/edfread.html). For loading large .edf files take into account that **edfreadm.m** can represent large loading throuput times. It is recommended a file or time-length segmentation before or as **edfread.m** input parameters.

## 2. GC Estimation using MVGC toolbox

After you get the preprocessed files in the **preprocessed_save**, you can proceed with the GC estimation following the logic/workflow in the function called **reading_eeg_saved_MVGC**. 
### GC definition following [Barret et al 2013](https://arxiv.org/pdf/1606.08644)

For defininng the GC between two EEG channels $x$ and $\hat{x}$, we must first establish the linear interaction between both channels in time as follows:

$$ 
    \hat{x(t)} = \sum_{k=0}^{p} A_{k} x(t-k) + \lambda \sigma(t)
$$

where:
$A_{k}$ are autoregressive coefficients or **VAR** matrix. 
- $p$ is the model maximum order.
- $\sigma(t)$ is the innovation (residual noise) or **SIGMA** matrix.
- $\lambda$ is a scaling factor.

In general way the **MVGC** toolbox evaluate whether $x$ Granger-causes $\hat{x}$. After $A_{k}$ and $sigma$ are estimatied, GC score is defined as taking into account a general model from all the possible contributions from any channel $y$:

$$
F_{x \to \hat{x}} = \ln \frac{\text{var}(\epsilon_x^{\text{reduced}})}{\text{var}(\epsilon_x^{\text{full}})}
$$

where:
- $\epsilon_x^{\text{full}}$ is the residual of the full model
- $\epsilon_x^{\text{reduced}}$ is the residual of the reduced model

A higher value of \( F_{y \to x} \) indicates stronger predictive influence from \( y \) to \( x \).

### Autocovariance Representation

The autocovariance function is defined as:

$$
\Gamma(\tau) = \operatorname{cov}(x(t), x(t - \tau)), \quad \tau = 0,1,2,\dots
$$

and the residual covariance matrix is denoted as:

$$
\Sigma
$$

These quantities are used internally in the MVGC toolbox for stable VAR estimation and spectral GC computation.

Before running any GC estimation take into account the following limitation of this toolbox:

1) EEG time-series are **highly colinear** and can generate singularities and numerical instabilities if model overfitting is not controlled propely. Measures, such as, Spectral Radius larger than one, or $A$. or $\Sigma$ prediction singularities are used by the toolbox to control matrix colinear ill-condition of each estimation. This information is shown during the estimation using **disp**.
2) The number of channels included in the analysis is also a very important factor to avoid overfitting and numerical instabilities on $A$ or $\Sigma$ estimations.


