%% Automated Signal Processing Script
% This script automates the process of loading data, plotting the initial figure, and loading trials.

%% Variable Name Changes

%% GUI Input Variables
% These variables come directly from GUI elements.

% File Selection & Data Acquisition
% FileNameEditField -> str_FullFileName -> fileName
% SampleFqHzEditField -> f_target -> targetSampleRate

% Channel & Signal Parameters
% ChannelEditField -> sel_Ch -> selectedChannel
% ThresholdEditField -> Th -> threshold

% Time Parameters
% InitialTimeEditField -> ti -> initialTime
% FinalTimeEditField -> tf -> finalTime
% InterStimTimeEditField -> tis -> interStimulusTime

% Trial Selection
% AllTrialsCheckBox -> all_trials -> allTrials
% SingleTrialEditField -> single_trial -> singleTrial
% SelectedTrialsEditField -> v_Trials -> selectedTrialIndices

% Time-Frequency Analysis Parameters
% F1 -> ps_MinFreqHz -> minFreqHz
% F2 -> ps_MaxFreqHz -> maxFreqHz
% CyclesEditField -> ps_StDevCycles -> stDevCycles
% NormCheckBox -> norm -> normalizeFlag
% LogCheckBox -> log -> logTransformFlag
% ResampleFr -> f_target -> targetSampleRate

% Compatibility Variables
% StimNumEditField -> Stim_num -> stimulusNumber
% NoStimCheckBox -> GF -> noStimulusFlag

%% Step-Based Processing Variables
% These variables are derived from processing steps in the script.

% Data Loading
% m_Data -> dataMatrix
% srate -> actualSampleRate
% h -> headerInfo
% v_Time -> timeVector
% num_Ch -> numChannels

% Trial Extraction & Processing
% m_Trials1 -> trialMatrix
% v_ini -> initialVector

% Trial Assembly
% m_Trials1_w -> trialWindowMatrix
% m_Trials1_wt -> trialTimeWindowMatrix
% v_Timew -> trialTimeVector
% m_Trials1_freq -> trialFrequencyMatrix

% Trial Selection & Filtering
% m_Trials1_selected -> selectedTrialsMatrix
% m_Trials1_freq_sel -> selectedTrialsFrequencyMatrix

% Time-Frequency Processing
% srt -> sampleRateReduction

clc; clear; close all;

%% Define File and Parameters

% File Information
fileName = '17308005.abf';

% Channel and Signal Parameters
selectedChannel = 1;
threshold = 0.8;

% Time Parameters
initialTime = -0.05;
finalTime = Inf;

% Sampling Rate Parameters
targetSampleRate = 1000;

% Trial Selection
allTrials = 1;
singleTrial = 1;
selectedTrialIndices = [2, 5, 10, 14, 21, 30];

% Compatibility Variables (Do Not Modify)
noStimulusFlag = 0;
interStimulusTime = 0.1;
stimulusNumber = 1;

% Time-Frequency Analysis Parameters
minFreqHz = 2;
maxFreqHz = 256;
stDevCycles = 2;
normalizeFlag = 0;
logTransformFlag = 1;

%% Load Signal Data
[dataMatrix, actualSampleRate, headerInfo, timeVector, numChannels] = f_LoadSignal(fileName);

%% Plot Initial Figure
f_Plot_Ch(timeVector, dataMatrix, headerInfo, numChannels);

%% Load Trials
[trialMatrix, initialVector] = f_Cols_Trials(dataMatrix, selectedChannel, threshold, stimulusNumber, initialTime, numChannels, actualSampleRate);
[trialWindowMatrix, timeEnd, trialFrequencyMatrix] = f_Assemble_Trials(trialMatrix, stimulusNumber, interStimulusTime, actualSampleRate);
[trialTimeWindowMatrix, trialTimeVector] = f_Assemble_Time(trialWindowMatrix, initialVector, initialTime, finalTime, timeEnd, actualSampleRate);

[selectedTrialsMatrix, selectedTrialsFrequencyMatrix] = f_Select_Trials(trialTimeWindowMatrix, selectedTrialIndices, trialTimeVector, allTrials, trialFrequencyMatrix, noStimulusFlag, headerInfo, selectedChannel);

%% Optionally plot a single trial
% f_plot_single(1, headerInfo, singleTrial, selectedChannel);

%% Resample Data
sampleRateReduction = round(actualSampleRate / targetSampleRate);

%% Time-Frequency Average Analysis
% f_tf_average(selectedTrialsMatrix, minFreqHz, maxFreqHz, sampleRateReduction, actualSampleRate, trialTimeVector, stDevCycles, normalizeFlag, logTransformFlag, headerInfo, selectedChannel);

%% Time-Frequency Histogram Analysis
f_tf_histogram(selectedTrialsMatrix, minFreqHz, maxFreqHz, sampleRateReduction, actualSampleRate, stDevCycles, normalizeFlag, headerInfo, selectedChannel);