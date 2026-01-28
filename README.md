# Individualised Temporal Interference: ITI
This repository holds code for the registered report titled "A comparison of individualized and generic Temporal Interference Stimulation of reward circuits". 
This project aims to test whether delivering Temporal Interference stimulation (TIS) to a personalised vs generic target more effectively modulates brain function, task performance, and the relationship between the two. This aim is achieved through analyses of fMRI data collected while participants play the Iowa Gambling Task and concurrently recieve sham, personalised, or generic TIS. 

This repository will contain both code and data. Below you will first find a description of the code in this repository, then a descriptive figure of the pipeline built by the code, and at the end, a description of the data. 

For any questions concerning this project, please email danielle.kurtin18@imperial.ac.uk. 

## Code
The code in this repository will analyse task performance, and identify which neural metric best relates to task performance. It will then use participant's neuroimaging data to inform the construction of whole-brain oscillatory models (i.e., virtual brain twins) of participant's brain function while playing the IGT. Then, using the models, we can identify which brain region should be stimulated to shift the target neural metric in the direction that relates to improved task performance. 
The scripts composing this pipeline are provided below. For each section of the pipeline, a description of each step is provided, followed by the name of the script that runs that step. 

**Compute brain-behaviour relationships** 
<ol>
<li> Analyse task performance - PVL, followed by BehavSummaryAndConcatenateRuns.m </li>
<li> Organise inputs: create a data frame of extracted regional timeseries per subject and session - CreateDF.m </li> 
<li> Compute measures of between-network connectivity for every trial of the Iowa Gambling Task (IGT) - ComputeTrialMetric.m </li> 
<li> Identify which connectivity metric best relates to IGT performance - ComputeBrainBehav_RandomForest.m </li> 
</ol>

**Build virtual brain twins**
<ol>
<li> Prepare structural connectivity data - ComputeSC.m </li>
<li> Compute the target frequency for each region - ComputeHopfFreq.m </li>
<li> Create the virtual brain twin - Ceff.m </li> 
<li> Identify the target region for stimulation - Perturb.m </li> 
</ol>

**Coming soon:**
Scripts to pre- and post-process the neuroimaging data

<img width="1654" height="2339" alt="Fig1" src="https://github.com/user-attachments/assets/e437ffdd-30d5-40ad-be0f-71a8ab955101" />

## Data
**Coming soon**
Extracted timeseries and other pilot data 
