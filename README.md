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

After you executed this preprocessing command you will observe the folllowing **.mat** files. The filename string will contain the suffix you added as input before the .mat extension.
In this case the suffix string is empty. This for the sake of comparison between different preprocessing options and for facilitating the exploration of resulting files across the **result_plots**, **preprocessed_save**, and **GC_estimation** folders.

<img width="709" height="398" alt="image" src="https://github.com/user-attachments/assets/fa0e2606-ddea-4234-88b7-626b82c7ffeb" />


Each file contains the EEGlab structure of each Near-Dead-End () stages proposed in [Xu et al 2023] between **S1-S11**, being S1 the baseline stage before removing the ventilator, and from S2-S11 all the subsequent stages without the ventilator, switching multiple time the peacemaker activation. 


