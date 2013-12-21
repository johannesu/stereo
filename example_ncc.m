% Shows of the use the code with a NCC data (unary) term
clear all; close all;
if (matlabpool('size') == 0)
	matlabpool open
end

%% Settings
% Image pair from
%http://vision.middlebury.edu/stereo/'
images{1} = double(imread('data/teddy/im2.png'));
images{2} = double(imread('data/teddy/im6.png'));

disparities = 0:1:50;
tol = 8*(disparities(2)-disparities(1));
unary_weight = 40;
kernel = 1;

%% Setup object
dm  = dispmap_ncc(images,disparities, kernel, unary_weight, tol);

%%
% Generate new proposals by samlping taking the
% best NCC disparity at each voxel ad fitting planes
radius = 5;
id = 0;
proposal_cell ={};
for x = 10:50:dm.sz(2)
	for y = 10:50:dm.sz(1)
		id = id+1;
		proposal_cell{id} =dm.generate_new_plane_RANSAC(x,y, radius); 
	end
end

% Add a few Fronto parallo
for d = 0:10:max(disparities)
   proposal = zeros(size(dm.assignment));
   proposal(3,:) = 1;
   proposal(4,:) = -d;

   proposal_cell{end+1} = proposal;
end


%% Iterative binary fusion
for iter = 1:length(proposal_cell)
  dm.binary_fusion(proposal_cell{iter});
	dm.display_current_dispmap();
	drawnow();
end

% Display
single_energy = dm.energy();
figure(1);
dm.display_current_dispmap;
xlabel('Iterative fusion');

%% Simultaneous fusion
dm.restart();
dm.simultaneous_fusion(proposal_cell);
sim_energy = dm.energy();

figure(2);
dm.display_current_dispmap
xlabel('Simultaneous fusion');
