# Time_Spectral_Granger_Causality_EEG
Granger Causality (GC) pipeline for evaluating channels predictability from data from [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120). For downloading the data refer to the Zenodo link [here](https://zenodo.org/record/7803212#.ZC3Cb-zML0q) and download the entire dataset. For any inquiry about data request or extra detail in code execution please don't hesitate to reach professor Jimo Borjigin, PhD [here](mailto:borjigin@umich.edu). Preliminary topoplot results showing spectral GC for patient **JF_20250225** are shown in the following Figure.

<img width="1980" height="1024" alt="image" src="https://github.com/user-attachments/assets/e911f1eb-835d-4dc8-9952-53a0f66be33e" />

An increased GC is observed in $\gamma_{1}$ and $\gamma_{2}$ supported by high causality measures on $\text{FC} \rightarrow \text{TPO}$ areas in **S4** and **S6** stages. **P4** and **O2** are causality hubs receiving most of the $\text{FC}$ connections from the entire  $\text{FC}$ cortex. Similar GC measures are concentrated in $\theta$ and $\alpha$ bands for **S3** and **S4** being consistent with the right PAC and NTSE activations reported in [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120).

## Requirements

Before running be sure you have installed **Matlab** version **R2025a** or more recent, as well as

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

**Depending on the throughput quality and the amount of cores your of your machine processor this preprocessing can take 15-20mins for each patient .edf file om average laptop**

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
In this case the suffix string is empty. This is done for evaluation a comparison between different preprocessing options and for facilitating the exploration of resulting files across the **result_plots**, **preprocessed_save**, and **GC_estimation** folders.

<img width="709" height="398" alt="image" src="https://github.com/user-attachments/assets/fa0e2606-ddea-4234-88b7-626b82c7ffeb" />


Each file contains an EEGlab structure for each Near-Dead-Event (NDE) stages described in [Xu et al 2023] between **S1-S11** stages: having S1 as baseline stage before removing the ventilator, and from S2-S11 all the subsequent stages without the ventilator, switching multiple time the peacemaker activation. The file denoted with the preffix **JF_20250225** contains information with all the stages **S1-S11** separated for an adeaquate load using [**edfread.m**](https://www.mathworks.com/help/signal/ref/edfread.html). For loading large .edf files take into account that **edfreadm.m** can represent large loading throuput times. It is recommended a file or time-length segmentation before or as **edfread.m** input parameters.

## 2. GC Estimation using MVGC toolbox

The broad pipeline of the GC evaluation is shown in the following Figure

<img width="1245" height="359" alt="image" src="https://github.com/user-attachments/assets/2606beaf-f96b-4fc3-abcd-22b79b41d554" />

After the preprocessed files are located in the **preprocessed_save**, it is possible to proceed with the GC estimation following the logic/workflow in the function called **reading_eeg_saved_MVGC**. 
### GC definition following [Barret et al 2013](https://arxiv.org/pdf/1606.08644)

For defining the GC between two EEG channels $x$ and $\hat{x}$, we must first establish the linear interaction between both channels in time as follows:

$$ 
    \hat{x(t)} = \sum_{k=0}^{p} A_{k} x(t-k) + \lambda \sigma(t)
$$

where:
$A_{k}$ are autoregressive coefficients or **VAR** matrix. 
- $p$ is the model maximum order.
- $\sigma(t)$ is the innovation (residual noise) or **SIGMA** matrix.
- $\lambda$ is a scaling factor.

In general way the **MVGC** toolbox evaluate whether $x$ Granger-causes $\hat{x}$. After $A_{k}$ and $sigma$ are estimatied, GC score is defined as taking into account a general model from all the possible contributions from any channel $y$. For extend the channel to channel GC to a full channel extended model, we can rewrite the interaction equation following  this:

$$
\hat{x}(t) = \sum_{k=1}^{p} A_k x(t-k) + \sum_{k=1}^{p} B_k y(t-k)  + \sigma(t)
$$

Thus, the GC estimatiion for each channel $y$ is defined as:

$$
F_{y \to \hat{x}} = \ln \frac{\text{var}(\epsilon_x^{\text{reduced}})}{\text{var}(\epsilon_x^{\text{full}})}
$$


where:
- $\epsilon_x^{\text{full}}$ is the residual of the full model
- $\epsilon_x^{\text{reduced}}$ is the residual of the reduced model

A higher value of $F_{y \to x}$ indicates stronger predictive influence from any $y$ to $\hat{x}$.

### Autocovariance Representation

The autocovariance function used for inferring $A_{k}$ and $sigma$ is defined as

$$
\Gamma(\tau) = \text{Cov}(x(t), x(t - \tau)), \quad
\Gamma(\tau) = \mathbb{E}(x(t),x(t - \tau)), \quad \tau = 0,1,2,\dots
$$

and the residual covariance matrix is denoted as $Sigma$.

These quantities are used internally in the MVGC toolbox for stable **VAR** estimation and spectral GC computation. With the following considerations.

---

### Autocovariance Representation

The autocovariance function used for inferring the VAR coefficients $A_k$ and the residual covariance $\Sigma$ is defined as:

$$
\Gamma(\tau) = \text{Cov}(x(t), x(t - \tau)) = \mathbb{E}[x(t)x(t - \tau)^\top], \quad \tau = 0,1,2,\dots
$$

where $\Gamma(\tau) \in \mathbb{R}^{n \times n}$ is the lag-$\tau$ covariance matrix for an $n$-dimensional signal. $tau$ is the maximum order of the Autocovariance model $p$.

---

### Toeplitz Structure

The sequence of autocovariance matrices defines a block Toeplitz matrix:

$$
\mathbf{\Gamma}_p =
\begin{bmatrix}
\Gamma(0) & \Gamma(1) & \Gamma(2) & \cdots & \Gamma(p-1) \\
\Gamma(1)^\top & \Gamma(0) & \Gamma(1) & \cdots & \Gamma(p-2) \\
\Gamma(2)^\top & \Gamma(1)^\top & \Gamma(0) & \cdots & \Gamma(p-3) \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
\Gamma(p-1)^\top & \Gamma(p-2)^\top & \Gamma(p-3)^\top & \cdots & \Gamma(0)
\end{bmatrix}
$$

This structure arises from the stationarity assumption of the VAR process.

---

### Yule–Walker Equations

The **VAR(p)** model parameters are estimated by solving the multivariate **Yule–Walker** equations:

$$
\Gamma(\tau) = \sum_{k=1}^{p} A_k \Gamma(\tau - k), \quad \tau = 1,2,\dots,p
$$

This linear system can be written compactly as:

```math
\Gamma_1 =
\begin{bmatrix}
\Gamma(1) & \Gamma(2) & \cdots & \Gamma(p)
\end{bmatrix}
=
\begin{bmatrix}
A_1 & A_2 & \cdots & A_p
\end{bmatrix}
\Gamma_p
```
---

### Residual Covariance Estimation

Once the coefficients $A_k$ are estimated, the residual covariance matrix is given by:

$$
\Sigma = \Gamma(0) - \sum_{k=1}^{p} A_k \Gamma(k)^\top
$$

The matrix $\Sigma$ represents the covariance of the innovation process and is critical for computing Granger causality.

An evaluation of numerical stability is reflecting to keep the Spectral Radius $\rho$ of $\Gamma$ lower than one in both **VAR** and **SIGMA** estimations
following this mathematical formulation:

```math
\mathcal{A} =
\begin{bmatrix}
A_1 & A_2 & \cdots & A_{p-1} & A_p \\
I   & 0   & \cdots & 0       & 0 \\
0   & I   & \cdots & 0       & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0   & 0   & \cdots & I       & 0
\end{bmatrix}
```

with

```math
\rho(\mathcal{A}) < 1
```

assuring a finite spectral rank of $\mathcal{A}$.

---

### Connection to MVGC

The MVGC toolbox uses this autocovariance sequence to:

- estimate stable VAR models
- compute time-domain and spectral GC
- ensure positive-definiteness of $\Sigma$
- verify model stability $\rho < 1$

Before running any GC estimation take into account the following limitation of this toolbox:

1) EEG time-series are **highly collinear** and can generate singularities and numerical instabilities if model overfitting is not controlled properly. Measures, such as, Spectral Radius $\rho$ smaller than one, or singularities in the coefficient $A$ or residual covariance $\Sigma$  are used by the toolbox to control matrix collinear ill-conditions on each estimation. This information is shown during the estimation using **disp**.
2) The number of channels included in the analysis is also a very important factor **to avoid overfitting and numerical instabilities** on $A$ or $\Sigma$ estimations. For this particular analysis we tested the GC between this group channels out of [Xu et al 2023 dataset](https://zenodo.org/record/7803212#.ZC3Cb-zML0q) including frontal, central, and Temporal-Parietal-Occipital (TPO) regions ['T3','T4','T5','T6','P3','P4','O1','O2','C3','C4','F3','F4','F7','F8'].
3) $p$ representing the maximum order of the estimation is also really important to avoid model overfitting or numerical instabilities. Extending an hyperparameter tuning is still an open plausible option in this case.

Thus, the function **reading_eeg_saved_MVGC**` implements the full pipeline for loading preprocessed EEG stages from **preprocessed_save** folder and performing GC analysis using the MVGC toolbox.

### Function Overview

The pipeline performs the following steps:

1. Load preprocessed EEG segments (S1–S11)
2. Visualize stacked EEG signals from those segments
3. Apply GC estimation using MVGC on each segments
4. Generate frequency-specific GC measures for the different EEG band specified in [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120), such as, $\delta$, $\theta$, $\alpha$, $\beta$, $\gamma_{1}$, and $\gamma_{2}$.

Use the following the following Matlab command for executing the 

```matlab
reading_eeg_saved_MVGC('JF_20250225', '_remove_midline')
```

The full GC estimation workflow consist in integrating data loading, visualization, and MVGC-based inference into a single execution pipeline.

At a high level, the function:

- Iterates over pre-segmented EEG stages (S1–S11) stored as `.mat` files  
- Reconstructs EEGLAB structures ('EEG_k') for each stage  
- Aggregates all stages for global visualization (stacked plots)  
- Performs stage-specific GC analysis by calling the core routine:
  
```matlab
  MVGC_application(X_sub_eeg.EEG_k, ...)
```

After the processes finishes you will obtain the following plots for the stacked EEG segments visualization

<img width="1102" height="545" alt="image" src="https://github.com/user-attachments/assets/e40a95aa-bae6-46fd-ba1a-017b99fc3979" />

And some examples of GC matrix in different frequency bands for instance $\gamma_{1}$ can be seen here

<img width="1156" height="572" alt="image" src="https://github.com/user-attachments/assets/b22174e2-1c1c-485f-a3b1-67f72ca41fee" />

<img width="1154" height="570" alt="image" src="https://github.com/user-attachments/assets/b72b5c8f-9f17-48ec-9976-562409d9f485" />

These preliminary results shows an **increased causality in $\text{FC} \rightarrow \text{TPO}$ areas in S4 and S6** for patient **JF_20250225**.

The complete execution of the GC calculation for all the $\delta$, $\theta$, $\alpha$, $\beta$, $\gamma_{1}$, and $\gamma_{2}$ of patient  **JF_20250225** takes **$\approx$ 5h, 36mins, and 17 seconds on an avarage laptop**.
