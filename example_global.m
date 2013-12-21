%% Reproducing Figure 4b) in "In Defense of 3D-Label Stereo"
% N.B. that the result may differ slightly because the first assigment is random

clear all; close all;
% Example usage of the code.
% Generate a disparty map object (dm) and iterativly perform binary fusion

%% Settings
% Image pair from
% http://vision.middlebury.edu/stereo/
images{1} = double(imread('data/teddy/im2.png'));
images{2} = double(imread('data/teddy/im6.png'));

% Settings from imrender's download_stereo.m
% Syntax for P matrix is defined in imrender/ojw/download_stereo.m
% Note it's not identical for every sequence in middlebury.
P = repmat([eye(3) zeros(3, 1)], [1 1 2]);
P(1,end,end) = -0.25;
disp_range = [0 59];
disparity_factor = 4;

%% Setup object
root = fileparts(which(mfilename));
addpath([root filesep 'imrender' filesep 'vgg']);
addpath([root filesep 'imrender' filesep 'ojw']);

options = ojw_default_options('cvpr08');
options.smoothness_kernel =1; % Default
dm  = dispmap_globalstereo(images,P, disp_range, disparity_factor, options);

%% Generate and fuse SegPln proposals
segplns = dm.segpln();

for iter = 1:length(segplns)
	dm.binary_fusion(segplns{iter});
	dm.display_current_dispmap();

	xlabel(sprintf('SegPln %d/%d \n', iter, length(segplns)))
	drawnow();
end