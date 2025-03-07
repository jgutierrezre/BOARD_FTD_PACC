% Function: f_ImageMatrix
%
% Descrption:
% This function makes the image of the data contained in a matrix. Each
% cell in the matrix is the center of a rectangular pixel.
%
% Parameters:
% pm_Data(*): Matrix to display as image
% pv_XAxis: Vector containing the bounds of x axis (see "help imagesc")
% pv_YAxis: Vector containing the bounds of y axis (see "help imagesc")
% pv_Limits: Two-element vector that limits the range of data values in the
% matrix (see "help imagesc")
% pstr_Color: String containing the name of the desired colormap (see
% "help colormap")
% ps_ColorLevels: Number of levels inside the colormap
% ps_NewFigure: Set to 1 to force the creation of a new figure
% ps_InvertImage: Set to 1 to invert upside down the resulting image.
% Default 1.
% ps_SetColorMap: Set to 1 to set the colormap passed above. Default 1.
% ps_NonEquAxis: Set to 1 if elements in X and/or Y arrays are not equally
% spaced.
%
% (*) Required parameters
%
function f_ImageMatrix( ...
        pm_Data, ...
        pv_XAxis, ...
        pv_YAxis, ...
        pv_Limits, ...
        pstr_Color, ...
        ps_ColorLevels, ...
        ps_NewFigure, ...
        ps_InvertImage, ...
        ps_SetColorMap, ...
        ps_NonEquAxis)

    if nargin < 1
        error('[f_ImageMatrix] - ERROR: Bad number of parameters')
    end

    v_Limits = [];
    str_Color = 'jet';
    s_ColorLevels = 256;
    s_NewFigure = 0;
    s_InvertImage = 1;

    if nargin >= 4 && ~isempty(pv_Limits)
        v_Limits = pv_Limits;
    end

    if nargin >= 5 && ~isempty(pstr_Color)
        str_Color = pstr_Color;
    end

    if nargin >= 6 && ~isempty(ps_ColorLevels)
        s_ColorLevels = ps_ColorLevels;
    end

    if nargin >= 7 && ~isempty(ps_NewFigure)
        s_NewFigure = ps_NewFigure;
    end

    if nargin >= 8 && ~isempty(ps_InvertImage)
        s_InvertImage = ps_InvertImage;
    end

    if ~exist('ps_SetColorMap', 'var') || isempty(ps_SetColorMap)
        ps_SetColorMap = 1;
    end

    if ~exist('ps_NonEquAxis', 'var') || isempty(ps_NonEquAxis)
        ps_NonEquAxis = 0;
    end

    if s_NewFigure
        figure
    end

    if ps_SetColorMap
        str_ColorMap = sprintf('colormap(%s(%d))', str_Color, s_ColorLevels);
        eval(str_ColorMap)
    end

    if ps_NonEquAxis
        pv_XAxis = pv_XAxis(:);
        v_XAxisAux = zeros(numel(pv_XAxis) + 1, 1);
        v_XAxisAux(2:end - 1) = diff(pv_XAxis) / 2;
        v_XAxisAux(1) = v_XAxisAux(2);
        v_XAxisAux(end) = v_XAxisAux(end - 1);
        v_XAxisAux(1:end - 1) = pv_XAxis - v_XAxisAux(1:end - 1);
        v_XAxisAux(end) = pv_XAxis(end) + v_XAxisAux(end);

        s_Min = min(v_XAxisAux);
        s_Max = max(v_XAxisAux);
        s_Dis = s_Max - s_Min;
        v_XAxisAux = (v_XAxisAux - s_Min) ./ s_Dis;

        s_Min = min(pv_XAxis);
        s_Max = max(pv_XAxis);
        s_Dis = s_Max - s_Min;
        v_XAxisAux = v_XAxisAux .* s_Dis;
        v_XAxisAux = v_XAxisAux + s_Min;

        pv_YAxis = pv_YAxis(:);
        v_YAxisAux = zeros(numel(pv_YAxis) + 1, 1);
        v_YAxisAux(2:end - 1) = diff(pv_YAxis) / 2;
        v_YAxisAux(1) = v_YAxisAux(2);
        v_YAxisAux(end) = v_YAxisAux(end - 1);
        v_YAxisAux(1:end - 1) = pv_YAxis - v_YAxisAux(1:end - 1);
        v_YAxisAux(end) = pv_YAxis(end) + v_YAxisAux(end);

        s_Min = min(v_YAxisAux);
        s_Max = max(v_YAxisAux);
        s_Dis = s_Max - s_Min;
        v_YAxisAux = (v_YAxisAux - s_Min) ./ s_Dis;

        s_Min = min(pv_YAxis);
        s_Max = max(pv_YAxis);
        s_Dis = s_Max - s_Min;
        v_YAxisAux = v_YAxisAux .* s_Dis;
        v_YAxisAux = v_YAxisAux + s_Min;

        m_DataAux = zeros(size(pm_Data, 1) + 1, size(pm_Data, 2) + 1);
        m_DataAux(1:end - 1, 1:end - 1) = pm_Data;
        m_DataAux(end, :) = min(pm_Data(:));
        m_DataAux(:, end) = min(pm_Data(:));

        pcolor(v_XAxisAux, v_YAxisAux, m_DataAux);
        shading(gca, 'flat');

        if s_InvertImage
            set(gca, 'YDir', 'reverse');
        end

        if ~isempty(v_Limits)
            set(gca, 'CLim', v_Limits);
        end

    else

        if ~isempty(v_Limits)
            imagesc(pv_XAxis, pv_YAxis, pm_Data, v_Limits);
        else
            imagesc(pv_XAxis, pv_YAxis, pm_Data);
        end

        if s_InvertImage
            axis('xy');
        end

    end

    return;
