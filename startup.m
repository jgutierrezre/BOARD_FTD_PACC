% startup.m
% --------------------------------
% Automatically sets up paths when MATLAB starts in this project

disp('Initializing MATLAB environment for BOARD_FTD_PACC...');

% Get the directory of this file (i.e., the project root)
projectRoot = fileparts(mfilename('fullpath'));

% Define the base directory to scan for scripts
srcDir = fullfile(projectRoot, 'src');

% Ensure the base directory exists
if exist(srcDir, 'dir')
    % Recursively add all subdirectories to the MATLAB path
    addpath(genpath(srcDir));
    disp(['Added all subdirectories under: ', srcDir]);
else
    warning(['Source directory "', srcDir, '" not found. No paths added.']);
end

disp('MATLAB project paths initialized successfully.');