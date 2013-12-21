%% Reproducing Figure 4 in "Simultaneous Fusion Moves for 3D-Label Stereo"
% "Surface view of a toy example on the Baby 2 sequence using l1
% regularization. Fusing only the 14 SegPln proposals".
% N.B. that the result may differ slightly because the first assigment is random
clear all; close all;

% Image pair from
% http://vision.middlebury.edu/stereo/'
images{1} = double(imread('data/baby2/im2.png'));
images{2} = double(imread('data/baby2/im6.png'));

% Settings from imrender's download_stereo.m
% Syntax for P matrix is defined in imrender/ojw/download_stereo.m
% Note it's not identical for every sequence in middlebury.
P = repmat([eye(3) zeros(3, 1)], [1 1 2]);
P(1,end,end) = -0.25;
disp_range = [0 85];
disparity_factor = 3;

%% Setup object
root = fileparts(which(mfilename));
addpath([root filesep 'imrender' filesep 'vgg']);
addpath([root filesep 'imrender' filesep 'ojw']);

options = ojw_default_options('cvpr08');
options.smoothness_kernel =1; % Default
dm  = dispmap_globalstereo(images,P, disp_range, disparity_factor, options);

%% Generate SegPln proposals
segplns = dm.segpln();

% Region of intrest
roi.x = 1:350;
roi.y = 50:360;

%% Binary fusion
show_steps = true;
iterations = dm.binary_fuse_until_convergence(segplns, show_steps);

% Display
figure(1);
clf; hold on;
subplot(1,2,1);
dm.display_surfaces(roi)
view(0,270);camorbit(0,5);camlight('headlight')
xlabel('Iterative fusion');

%% Simultaneous fusion 
dm.restart();
dm.maxiter = 3000;
dm.max_relgap = 1e-5;
[e,lb, iterations] = dm.simultaneous_fusion(segplns);

% Display
subplot(1,2,2); 
dm.display_surfaces(roi);
view(0,270);camorbit(0,5);camlight('headlight')
xlabel('Simultaneous fusion');