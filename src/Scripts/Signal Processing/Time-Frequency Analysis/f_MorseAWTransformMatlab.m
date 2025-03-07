% Function: f_MorseAWTransformMatlab.m
%
% Description:
% This function calculates a Time-Frequency Transform using an Analytic
% Morse Wavelet.
% The Sample Rate of the input signal is considered in order to compute
% the transform.
%
% Parameters:
% pv_Signal(*): Signal to process
% ps_SampleRate(*): Sample rate
% ps_MinFreqHz: Min frequency (in Hz) to process from
% ps_MaxFreqHz: Max frequency (in Hz) to process to
% ps_FreqSeg: Number of segments used to calculate the size of the
% resulting matrix in the frequency direction
%
% ps_NumOfCycles: The number of cycles you want to include for the
% transform at every frequency (the wavelet duration). In other words,
% this parameter corresponds to the number of oscillations at the peak
% frequency that fit within the central window of the time-domain wavelet.
%
% ps_Magnitudes: Set to 1 (default) if the magnitudes of the coefficients
% must be returned; 0 for analytic values (complex values).
%
% ps_SquaredMag: Set to 1 if the magnitudes of the coefficients divided by
% the squared of the corresponding scale must by power to 2
%
% ps_MakeBandAve: Set to 1 if instead of returning a matrix with all values
% in the time-frequency map, the function returns just a vector with the
% average along all the frequency scales for each time moment.
%
% ps_Phases: Set to 1 if the phases of the coefficients
% must be returned; 0 for analytic values (complex values).
%
% ps_TimeStep: Time step between values that are going to be kept in the
% output matrix. Each time moment is the average of the previous values
% according to the size of the window defined by this parameter.
%
% ps_DyadicScale: Set to 1 if the scale must be dyadic or 0 (default) if
% not. Dyadic scale means that the scale vector is constructed by taking
% the log2 of ps_MinFreqHz and ps_MaxFreqHz and create a linear vector
%
% ps_Beta: Beta parameter from the Generalized Morse Analytic Wavelet.
% This parameter is adjusted automatically when the ps_NumOfCycles
% parameter is non empty.
%
% ps_Gamma: Gamma parameter from the Generalized Morse Analytic Wavelet.
% This parameter is set to 3 by default (see reference below).
%
% Lilly and Olhede (2012), Generalized Morse wavelets as a superfamily
% of analytic wavelets. IEEE Trans. Sig. Proc., 60 (11), 6036--6041.
%
% (*) Required parameters
%
% Outputs:
% m_MorseWT: Matrix containing the scalogram. Time in rows, frequency in
% colums. Frequencies in descending order
% v_TimeAxis: Array containing the time axis values (second units)
% v_FreqAxis: Array containing the frequency axis values in descending
% order (Hz units)
%
% Author: Mario Valderrama
%
function [m_MorseWT, v_TimeAxis, v_FreqAxis] = ...
        f_MorseAWTransformMatlab( ...
        pv_Signal, ...
        ps_SampleRate, ...
        ps_MinFreqHz, ...
        ps_MaxFreqHz, ...
        ps_FreqSeg, ...
        ps_NumOfCycles, ...
        ps_Magnitudes, ...
        ps_SquaredMag, ...
        ps_MakeBandAve, ...
        ps_Phases, ...
        ps_TimeStep, ...
        ps_DyadicScale, ...
        ps_Beta, ...
        ps_Gamma)

    if nargin < 2
        return;
    end

    if ~exist('ps_MinFreqHz', 'var') || isempty(ps_MinFreqHz) || ...
            ps_MinFreqHz == 0
        ps_MinFreqHz = 0.1;
    end

    if ~exist('ps_MaxFreqHz', 'var') || isempty(ps_MaxFreqHz) || ...
            ps_MaxFreqHz > ps_SampleRate / 2;
        ps_MaxFreqHz = ps_SampleRate / 2;
    end

    if ~exist('ps_FreqSeg', 'var') || isempty(ps_FreqSeg) || ...
            ps_FreqSeg <= 0
        ps_FreqSeg = round(ps_MaxFreqHz - ps_MinFreqHz);
    end

    if ~exist('ps_NumOfCycles', 'var') || isempty(ps_NumOfCycles)
        ps_NumOfCycles = 3;
    end

    if ~exist('ps_Magnitudes', 'var') || isempty(ps_Magnitudes)
        ps_Magnitudes = 1;
    end

    if ~exist('ps_SquaredMag', 'var') || isempty(ps_SquaredMag)
        ps_SquaredMag = 0;
    end

    if ~exist('ps_MakeBandAve', 'var') || isempty(ps_MakeBandAve)
        ps_MakeBandAve = 0;
    end

    if ~exist('ps_Phases', 'var') || isempty(ps_Phases)
        ps_Phases = 0;
    end

    if ~exist('ps_TimeStep', 'var')
        ps_TimeStep = [];
    end

    if ~exist('ps_DyadicScale', 'var') || isempty(ps_DyadicScale)
        ps_DyadicScale = 0;
    end

    if ~exist('ps_Beta', 'var')
        ps_Beta = [];
    end

    if ~exist('ps_Gamma', 'var') || isempty(ps_Gamma)
        ps_Gamma = 3;
    end

    pv_Signal = pv_Signal(:);

    %     display('[f_MorseAWTransformMatlab] - Computing wavelet transform...');
    %     tic

    if ps_DyadicScale
        s_FreqStep = (log2(ps_MaxFreqHz) - log2(ps_MinFreqHz)) / (ps_FreqSeg - 1);
        v_FreqAxis = log2(ps_MinFreqHz):s_FreqStep:log2(ps_MaxFreqHz);
        v_FreqAxis = 2 .^ v_FreqAxis;
    else
        s_FreqStep = (ps_MaxFreqHz - ps_MinFreqHz) / (ps_FreqSeg - 1);
        v_FreqAxis = ps_MinFreqHz:s_FreqStep:ps_MaxFreqHz;
    end

    v_FreqAxis = fliplr(v_FreqAxis);

    s_IsEven = 0;

    if mod(numel(pv_Signal), 2) == 0
        pv_Signal = pv_Signal(1:end - 1);
        s_IsEven = 1;
    end

    v_TimeAxis = (0:numel(pv_Signal) - 1) ./ ps_SampleRate;
    s_Len = numel(v_TimeAxis);
    s_SqrtLen = sqrt(s_Len);
    s_HalfLen = floor(s_Len / 2) + 1;

    v_WAxis = (2 .* pi ./ s_Len) .* ...
        (0:(s_Len - 1));
    v_WAxis = v_WAxis .* ps_SampleRate;
    v_WAxisHalf = v_WAxis(1:s_HalfLen);

    if isempty(ps_TimeStep)
        s_SampAve = 1;
    else
        s_SampAve = round(ps_TimeStep * ps_SampleRate);

        if s_SampAve < 1
            s_SampAve = 1;
        end

    end

    v_SampAveFilt = [];

    if s_SampAve > 1
        v_IndSamp = 1:s_SampAve:numel(v_TimeAxis);
        v_TimeAxis = v_TimeAxis(v_IndSamp);
        v_SampAveFilt = ones(s_SampAve, 1);
    end

    v_InputSignalFFT = fft(pv_Signal, numel(pv_Signal));

    clear m_MorseWT
    m_MorseWT = zeros(numel(v_FreqAxis), numel(v_TimeAxis));

    if isempty(ps_Beta)
        ps_Beta = (ps_NumOfCycles * pi) ^ 2 / ps_Gamma;
    end

    s_CenterPeakRad = (ps_Beta / ps_Gamma) ^ (1 / ps_Gamma);

    s_FreqInd = 0;

    for s_FreqCounter = v_FreqAxis

        s_Scale = s_CenterPeakRad / (2 * pi * s_FreqCounter);

        v_WAxisHalfAux = v_WAxisHalf * s_Scale;

        clear v_WinFFT
        v_WinFFT = zeros(s_Len, 1);
        v_WinFFT(1:s_HalfLen) = 2 .* exp(((ps_Beta / ps_Gamma) * ...
            ((1 + log(ps_Gamma)) - log(ps_Beta))) + ...
            ps_Beta .* log(v_WAxisHalfAux) - v_WAxisHalfAux .^ ps_Gamma);
        %         v_WinFFT = v_WinFFT.* s_SqrtLen./ norm(v_WinFFT, 2);

        s_FreqInd = s_FreqInd + 1;

        if s_SampAve > 1
            clear v_MorseTemp
            v_MorseTemp = zeros(s_Len + (s_SampAve - 1), 1);
            v_MorseTemp(s_SampAve:end) = ifft(v_InputSignalFFT .* v_WinFFT) ./ ...
                sqrt(s_StDevSec);

            if ps_Magnitudes
                v_MorseTemp = abs(v_MorseTemp);
            end

            if ps_SquaredMag
                v_MorseTemp = v_MorseTemp .^ 2;
            end

            v_MorseTemp(1:(s_SampAve - 1)) = ...
                flipud(v_MorseTemp(s_SampAve + 1:2 * s_SampAve - 1));
            v_MorseTemp = filter(v_SampAveFilt, 1, v_MorseTemp) ./ s_SampAve;
            v_MorseTemp = v_MorseTemp(s_SampAve:end);

            m_MorseWT(s_FreqInd, :) = v_MorseTemp(v_IndSamp);
        else
            m_MorseWT(s_FreqInd, :) = ifft(v_InputSignalFFT .* v_WinFFT);
        end

    end

    clear v_WinFFT v_MorseTemp v_SampAveFilt

    %     toc

    if s_SampAve > 1
        return;
    end

    if s_IsEven
        m_MorseWT = horzcat(m_MorseWT, m_MorseWT(:, end));
        v_TimeAxis = horzcat(v_TimeAxis, v_TimeAxis(end) + ...
            (1 / ps_SampleRate));
    end

    if ps_Phases
        m_MorseWT = angle(m_MorseWT);
        return;
    end

    if ps_Magnitudes ~= 1
        return;
    end

    m_MorseWT = abs(m_MorseWT);

    if ps_SquaredMag
        m_MorseWT = m_MorseWT .^ 2;
    end

    if ps_MakeBandAve
        m_MorseWT = mean(m_MorseWT, 2);
        m_MorseWT = flipud(m_MorseWT);
        v_TimeAxis = [];
        v_FreqAxis = fliplr(v_FreqAxis);
    end

    return;
