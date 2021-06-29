

set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);

[filepath, ~, ~] = fileparts(mfilename('fullpath'));
cd(filepath);
filepath = filepath + "/..";
addpath(genpath(fullfile(filepath)));
cd ../;
cd components-library;