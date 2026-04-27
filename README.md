# Time_Spectral_Granger_Causality_EEG
Granger Causality (GC) pipeline for evaluating channels predictability from data from [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120). For downloading the data refer to the Zenodo link [here](https://zenodo.org/record/7803212#.ZC3Cb-zML0q) and download the entire dataset. For any inquiry about data request or extra detail in code execution please don't hesitate to reach professor Jimo Borjigin, PhD [here](mailto:borjigin@umich.edu). 

Preliminary topoplot results with the spectral GC meaaures for patient **JF_20250225** are shown in the following Figures. We implemented a sensitivity analysis across different unipolar re-reference on **Fz**, **Cz**, and **Pz** channels for testing re-reference variability on the middle line, as well as average reference. 

[Faes et al 2011](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.051112) and [Trongnetrpunya et al 2015](https://www.frontiersin.org/journals/systems-neuroscience/articles/10.3389/fnsys.2015.00189/full) suggest to use a multipolar reference after the application of ASR and ICA artifact removal, to reduce sensitivity in GC estimations. All the GC pairwise directions over the half of the quartile of the amplitude per frequency are statisticially significant after FDR correction $t(2,135) \geq 5.47, p<0.001$. 

These plots included the GC measures for **average** reference 

<img width="1810" height="753" alt="image" src="https://github.com/user-attachments/assets/842e2bc9-28da-401f-83da-14d8f388384e" />

<br>
<br>

for **Fz** reference


<img width="1809" height="758" alt="image" src="https://github.com/user-attachments/assets/f07c0048-412d-4588-abf2-0cf69c7dd514" />

<br>
<br>

for **Cz** reference

<img width="1812" height="764" alt="image" src="https://github.com/user-attachments/assets/b5f1ae92-8cfb-4437-bf7c-944d276be58b" />

<br>
<br>

and for **Pz** reference


<img width="1812" height="762" alt="image" src="https://github.com/user-attachments/assets/d3292a57-af6d-424d-a5b2-6814b9e2ffd6" />

<br>
<br>

An increased GC is observed in $\gamma_{1}$ and $\gamma_{2}$ supported by high causality measures on $\text{FC} \rightarrow \text{TPO}$ areas in **S5** and **S7** stages. **TPO** channels are causality hubs receiving most of the $\text{FC}$ connections being concentrated in the left for **S2** and **S3** stages and in the right. This GC measures increased are concentrated in $\gamma_{1}$ and $\gamma_{2}$ are consistent with the bilateral Phase Amplitude Coupling (PAC) and NTSE activations reported in [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120). 


GC contributions between S2-S11 stages are shown in the following Figure for $\gamma_{1}$


<img width="1724" height="734" alt="image" src="https://github.com/user-attachments/assets/0f7a2dda-4fd7-4fd2-9a0a-ac05e082eceb" />

<br>
<br>

and for $\gamma_{2}$ here


<img width="1728" height="737" alt="image" src="https://github.com/user-attachments/assets/add816e2-e40f-479b-b009-c67de68d5572" />

<br>
<br>

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
time_frequency_gc_single_edf_removing_ECG_definitive('<edf_filename_local>', '<output_file_suffix>', '<cell_preprocessing_activation>')
```

The **cell_preprocessing_activation** must contain three values **{ASR_activator, ICA_activator, ECG_removal_activator}**, if any of the values in the cell is one. The corresponding prepreocessing phase will be activated (e.g. ASR, ICA, or ECG removal  stage-lagged regressor) and a corresponding suffix will be added as an identified in the **result_plots** and **preprocessed_save** folders.

For instance it is possible to evaluate the preprocessing with average rerefence using this command

```matlab
time_frequency_gc_single_edf_removing_ECG_definitive('JF_20250225', '_average_reref', {1, 1, 1})
```

Activating all the preprocessing knobs as Independent Components Analysis [**ICA**](https://sccn.ucsd.edu/~arno/eeglab/auto/runica.html), [**ASR**](https://eeglab.org/plugins/clean_rawdata/), and [**ECG removal**](https://github.com/meiyor/Time_Spectral_Granger_Causality_EEG/blob/main/utils/remove_ecg_multiple_stages_lagged_regression.m).

Or with Fz rereference with the following command

```matlab
time_frequency_gc_single_edf_removing_ECG_definitive('JF_20250225', '_Fz_reref', {1, 1, 1})
```

**Depending on the throughput quality and the amount of cores your of your machine processor this preprocessing can take 15-20mins for each patient .edf file om average laptop**

The function **time_frequency_gc_single_edf_removing_ECG_definitive** performs the following preprocessin in that precise order to avoid GC sensitivity across different channel re-references and preprocessing selection:
 
  1) **Data Loading**
     
       Supports large EDF files via chunked reading, or direct loading using biosig depending on the size of the file.
     
  2) **Resampling**
     
       Downsample to 256 Hz.
  
  3) **Filtering**
     
    a) Notch filter at 60 Hz.
    b) Bandpass filter (0.1–100 Hz).

  4) **Channel Localization**
     Standard 10–20 montage assignment for now
     
  5) **ECG Artifact Removal**
     
     Stage-lagged linear regression using ECG detrended and filtered signal between 1-100Hz. OLS regressor is used with a window of 80ms and including the derivative components of the ECG filtered signal. No Ridge regularization is needed.
   
  6)  **Artifact Reduction**
     
    a) ASR-based cleaning (light EMG suppression).
    b) ICA decomposition (runica).
    c) Automatic IC rejection (heuristic-based) adding light severity suppressing for avoiding GC matrices singularity - even using pseudoinverse approaches.

  7) **Re-referencing**
       Average reference or Fz reference.
     
     
  8) **GC Analysis** (optional) **can be evaluated later**
     
     Time-domain and spectral GC via **MVGC toolbox**, explained in the following section.

After executing this preprocessing command you will observe the following **.mat** files. The filename string will contain the suffix you added as input before the .mat extension.
In this case the suffix string is empty. This is done for evaluation a comparison between different preprocessing options and for facilitating the exploration of resulting files across the **result_plots**, **preprocessed_save**, and **GC_estimation** folders.

<img width="709" height="398" alt="image" src="https://github.com/user-attachments/assets/fa0e2606-ddea-4234-88b7-626b82c7ffeb" />


Each file contains an EEGlab structure for each Near-Dead-Event (NDE) stages described in [Xu et al 2023] between **S1-S11** stages: having S1 as baseline stage before removing the ventilator, and from S2-S11 all the subsequent stages without the ventilator, switching multiple time the peacemaker activation. The file denoted with the preffix **JF_20250225** contains information with all the stages **S1-S11** separated for an adeaquate load using [**edfread.m**](https://www.mathworks.com/help/signal/ref/edfread.html). For loading large .edf files take into account that **edfreadm.m** can represent large loading throuput times. It is recommended a file or time-length segmentation before or as **edfread.m** input parameters.

## 2. GC Estimation using MVGC toolbox

The broad pipeline of the GC evaluation is shown in the following Figure.

<img width="1818" height="618" alt="image" src="https://github.com/user-attachments/assets/c741f201-fa12-41a1-b31e-80066fa40a18" />


After the preprocessed files are located in the **preprocessed_save**, it is possible to proceed with the GC estimation following the logic/workflow in the function called **reading_eeg_saved_MVGC**. 
### GC definition following [Barnett et al 2013](https://arxiv.org/pdf/1606.08644)

For defining the GC between two EEG channels $x$ and $\hat{x}$, we must first establish the linear interaction between both channels in time as follows:

$$ 
    \hat{x(t)} = \sum_{k=0}^{p} A_{k} x(t-k) + \lambda \sigma(t)
$$

where:
$A_{k}$ are autoregressive coefficients or **VAR** matrix. 
- $p$ is the model maximum order.
- $\sigma(t)$ is the innovation (residual noise) or **SIGMA** matrix.
- $\lambda$ is a scaling factor.

The **MVGC** toolbox evaluate whether $x$ Granger-causes $\hat{x}$. After $A_{k}$ and $\sigma$ are estimatied, GC score is defined as taking into account a general model from all the possible contributions from any channel $y$. For extend the channel to channel GC to a full channel extended model, we can rewrite the interaction equation following this formulation:

$$
\hat{x}(t) = \sum_{k=1}^{p} A_k x(t-k) + \sum_{k=1}^{p} B_k y(t-k)  + \sigma(t)
$$

GC estimatiion for each channel $y$ is defined as:

$$
F_{y \to \hat{x}} = \ln \frac{\text{var}(\epsilon_x^{\text{reduced}})}{\text{var}(\epsilon_x^{\text{full}})}
$$


where:
- $\epsilon_x^{\text{full}}$ is the residual of the full model
- $\epsilon_x^{\text{reduced}}$ is the residual of the reduced model

A higher value of $F_{y \to x}$ indicates stronger predictive influence from any $y$ to $\hat{x}$.

### Autocovariance Representation

The autocovariance function used for inferring $A_{k}$ and $\sigma$ is defined as

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

where $\Gamma(\tau) \in \mathbb{R}^{n \times n}$ is the lag-$\tau$ covariance matrix for an $n$-dimensional signal. $\tau$ is the maximum order of the Autocovariance model $p$.

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
- verify model stability $\rho(\mathcal{A}) < 1$

Before running any GC estimation take into account the following limitation of this toolbox:

1) EEG time-series are **highly collinear** and can generate singularities and numerical instabilities if model overfitting is not controlled properly. Measures, such as, Spectral Radius $\rho$ smaller than one, or singularities in the coefficient $A$ or residual covariance $\Sigma$  are used by the toolbox to control matrix collinear ill-conditions on each estimation. This information is shown during the estimation using **disp**.
2) The number of channels included in the analysis is also a very important factor **to avoid overfitting and numerical instabilities** on $A$ or $\Sigma$ estimations. For this particular analysis we tested the GC between this group channels out of [Xu et al 2023 dataset](https://zenodo.org/record/7803212#.ZC3Cb-zML0q) including frontal, central, and Temporal-Parietal-Occipital (TPO) regions ['T3','T4','T5','T6','P3','P4','O1','O2','C3','C4','F3','F4','F7','F8'].
3) $p$ representing the maximum order of the estimation is also really important to avoid model overfitting or numerical instabilities. Extending an hyperparameter tuning is still an open plausible option in this case.

The function **reading_eeg_saved_MVGC** implements the full pipeline for loading preprocessed EEG stages from **preprocessed_save** folder and performing GC analysis using the MVGC toolbox.

### Frequency Decomposition using Geweke Spectral Granger Causality

After time-domain GC measure is obtained, the MVGC toolbox calculates the directional predictability as a function of frequency using the formulation introduced by [Geweke](https://www.tandfonline.com/doi/pdf/10.1080/01621459.1982.10477803?casa_token=8wvugzOgOysAAAAA:Uc3dtknMgdeaTNUygova2AcwXDcfnceRnpuAH4fJmFxLNmuunewjYTW04SGkHMsK2sf67sPGaTAmww). This allows evaluating whether one channel predicts another preferentially the specific EEG bands evaluated here, such as, $\delta$, $\theta$, $\alpha$, $\beta$, $\gamma_{1}$, $\gamma_{2}$). This process consists in:


Calculating the transfer function from the GC in time and translate it in the frequency domain as:

```math
H(f)=
\left(
I-\sum_{k=1}^{p} A_k e^{-i2\pi fk}
\right)^{-1}
```

where:

- $H(f)$ is the frequency response of the multivariate system  
- $I$ is the identity matrix  
- $f$ is frequency in Hz, normalized internally by sampling rate of the input signal, in this case 256 Hz.

The spectral density matrix is then calculated as:

$$
S(f)=H(f)\Sigma H(f)^{*}
$$

where $\cdot^{*}$ is the conjugate transpose.

For the directional influence from channel $j$ to channel $i$, Geweke spectral GC is defined as:

```math
f_{j \rightarrow i}(f)
=
\ln
\frac
{S_{ii}(f)}
{S_{ii}(f)-H_{ij}(f)\Sigma_{jj}H_{ij}^{*}(f)}
```

where:

- $S_{ii}(f)$ is the autospectrum of target channel $i$
- $H_{ij}(f)$ is the transfer component between the spectrum from source $j$ to target $i$
- $\Sigma_{jj}$ is the innovation variance of source channel $j$

These spectral components are then calcuated following the Geweke formulations with

$$
S_{ii}(f)=\left[S(f)\right]_{ii}
$$

the $i$-th diagonal element of the spectral density matrix

$$
S(f)=H(f)\Sigma H(f)^{*}
$$

representing the autospectral power of target channel $i$ at frequency $f$.

$$
H_{ij}(f)=\left[H(f)\right]_{ij}
$$

and with $(i,j)$-th element of the spectral transfer matrix

```math
H(f)=\left(I-\sum_{k=1}^{p}A_k e^{-i2\pi fk}\right)^{-1}
```

representing the frequency-domain propagation from source channel $j$ to target channel $i$.

$$
\Sigma_{jj}=\left[\Sigma\right]_{jj}
$$

the $j$-th diagonal element of the innovation covariance matrix

$$
\Sigma=\text{Cov}(\varepsilon_t)
$$

Summarizing these spectral measures, a larger value of $f_{j \rightarrow i}(f)$ indicates stronger predictive influence from channel $j$ to channel $i$ at frequency $f$.

The time-domain GC is then recovered by integrating spectral GC across frequencies:

```math
F_{j \rightarrow i}
=
\frac{1}{2\pi}
\int_{-\pi}^{\pi}
f_{j \rightarrow i}(\omega)\, d\omega
```

This decomposition is particularly useful for the EEG representation evaluated here because directional interactions often emerge selectively in oscillatory bands rather than uniformly across the spectrum. In this repository, band-wise averages are reported for:

- $\delta$ : 0–4 Hz  
- $\theta$ : 4–8 Hz  
- $\alpha$ : 8–13 Hz  
- $\beta$ : 13–25 Hz  
- $\gamma_{1}$ : 25–55 Hz  
- $\gamma_{2}$ : 80–100 Hz

### Function Overview

The pipeline performs the following steps:

1. Load preprocessed EEG segments (S1–S11)
2. Visualize stacked EEG signals from those segments
3. Apply GC estimation using MVGC on each segments
4. Generate frequency-specific GC measures for the different EEG band specified in [Xu et al 2023](https://www.pnas.org/doi/abs/10.1073/pnas.2216268120), such as, $\delta$, $\theta$, $\alpha$, $\beta$, $\gamma_{1}$, and $\gamma_{2}$ as described above.

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
