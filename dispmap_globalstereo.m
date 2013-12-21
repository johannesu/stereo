% Johannes Ulén 2013
%
% Uses the unary term from
% Woodford, Oliver, et al.
% "Global stereo reconstruction under second-order smoothness priors."
% Pattern Analysis and Machine Intelligence, IEEE Transactions on 31.12 (2009): 2115-2128.
% Alot of the code below is modification of the source code from
% http://www.robots.ox.ac.uk/~ojw/software.htm.

classdef dispmap_globalstereo < dispmap_super
	
	properties (SetAccess = public)
		tol;
	end
	
	properties (SetAccess = protected)
		options;
		ephoto;
		P;
		disp_range;
		disparity_factor;
		disps;
		d_min;
		d_step;
		start_disparity;
	end
	
	methods
		function self = dispmap_globalstereo(images,P, disp_range, disparity_factor, options)
			% Superclass constructor
			kernel = options.smoothness_kernel;
			self = self@dispmap_super(images, kernel);
			self.tol = options.disp_thresh;
			
			% Using the orignal code from woodford to get settings segment image et cetera.
			root = fileparts(which(mfilename));
			addpath([root filesep 'imrender' filesep 'vgg']);
			addpath([root filesep 'imrender' filesep 'ojw']);
			
			if max(abs(P([1:6 9])-[1 0 0 0 1 0 1])) > 1e-12
				error('First image must be reference image');
			end
			self.P = permute(P(:,:,:), [2 1 3]);
			
			self.disp_range = disp_range;
			self.disparity_factor  = disparity_factor;
			
			self.disps = disp_range(1)*disparity_factor:disp_range(2)*disparity_factor;
			self.disps = sort(self.disps,'descend');
			
			self.d_min = self.disps(end);
			self.d_step =  self.disps(1) - self.d_min;
			self.options = options;
			
			preprocess(self);
			self.start_disparity = rand(self.sz) * self.d_step + self.d_min;
			init_solution(self);
		end
		
		function proposals = segpln(self)
			images = self.images;
			disps = self.disps;
			options = self.options;
			
			% Modification of ojw_segpln.m
			% which keeps the normal direction of the approximated planes
			P = permute(self.P(:,:,:), [2 1 3]);
			R = uint8(images{1});
			sz = size(R);
			if numel(sz) < 3
				sz(3) = 1;
			end
			
			Rvec = reshape(double(R),[],sz(3));
			[X Y] = meshgrid(1:sz(2), 1:sz(1));
			WC = ones(sz(1)*sz(2), 3);
			WC(:,1) = X(:);
			WC(:,2) = Y(:);
			X = sz - 2 * options.window;
			corr = zeros(X(1), X(2), numel(disps));
			filt = fspecial('average', [1 1+2*options.window]);
			
			% For each image...
			for a = 1:numel(images)
				% Project the points
				X = WC * P(:,1:3,a)';
				P_ = P(:,4,a);
				
				% For each disparity...
				for b = 1:numel(disps)
					% Vary image coordinates according to disparity
					d = disps(b) * P_;
					Z = 1 ./ (X(:,3) + d(3));
					
					% Look up the colours
					Y = squeeze(vgg_interp2(images{a}, (X(:,1) + d(1)) .* Z, (X(:,2) + d(2)) .* Z, 'linear', -1000));
					
					% Calculate the RSSD
					% (MATLAB wont nest anon. functions)
					Y = self.ephoto(Y-Rvec);
					Y = conv2(filt, filt', reshape(Y, sz(1:2)), 'valid');
					corr(:,:,b) = corr(:,:,b) + Y;
				end
			end
			
			% Normalize
			X = self.ephoto(-1000 -Rvec) * numel(images);
			corr = (X(1) - corr) / X(1);
			clear X Y Z h1 h2
			
			% Extract highest scoring matches (winner takes all)
			[info.corr corr] = max(corr, [], 3);
			corr = disps(corr);
			corr(info.corr<0.07) = 0;
			corr = padarray(corr, options.window([1 1]), 'symmetric'); % Return to original size
			
			% Generate image segmentations
			if size(R, 3) == 1
				R = repmat(R, [1 1 3]);
			end
			segment_params = [1 1.5 10 100];
			mults = [1:7 3 5 8 12 24 50 100];
			nMaps = numel(mults);
			info.segments = zeros(sz(1), sz(2), nMaps, 'uint32');
			for b = 1:nMaps
				sp = segment_params * mults(b);
				if b < 8
					% Segment the image using mean shift
					info.segments(:,:,b) = vgg_segment_ms(R, sp(1), sp(2), sp(3));
				else
					% Segment the image using Felzenszwalb's method
					info.segments(:,:,b) = vgg_segment_gb(R, 0, sp(4), sp(3), 1);
				end
			end
			
			clear A
			
			% World coordinates for plane fitting
			[X Y] = meshgrid(1:sz(2), 1:sz(1));
			WC = zeros(sz(2)*sz(1), 3);
			WC(:,3) = 1 ./ corr(:);
			WC(:,2) = WC(:,3) .* Y(:);
			WC(:,1) = WC(:,3) .* X(:);
			clear X Y
			
			% Switch off annoying warnings
			warning_state = warning('query', 'all');
			warning off MATLAB:divideByZero
			warning off MATLAB:singularMatrix
			warning off MATLAB:nearlySingularMatrix
			warning off MATLAB:illConditionedMatrix
			warning off MATLAB:rankDeficientMatrix
			
			% Generate piecewise-planar disparity maps
			
			% Initilize
			proposals = cell(nMaps,1);
			
			% 0 disparity
			prop = zeros(4,sz(1)*sz(2));
			prop(3,:) = 1;
			for i = 1:nMaps
				proposals{i} = prop;
			end
			
			rt = 0.1; %2 * min(abs(diff(Z(:))));
			for b = 1:nMaps
				for a = 1:max(max(info.segments(:,:,b)))
					% Choose a segment
					M = info.segments(:,:,b) == a;
					N = WC(M,:);
					N = N(N(:,3)~=0,:);
					
					local_WC_points = N;
					if size(N, 1) > 3
						% Ransac to weed out outliers
						M_ = rplane(self, N, rt);
						local_WC_points = N(M_,:);
					end
					
					if size(local_WC_points,1) > 2
						% Find least squares plane from inliers
						N_ = local_WC_points \ repmat(-1, [size(local_WC_points, 1) 1]);
						
						[Y X] = ind2sub(sz, find(M));
						D = -(X * N_(1) + Y * N_(2) + N_(3));
						
						p = [repmat([N_(1) N_(2) 1 N_(3)], [numel(X) 1])];
						proposals{b}(:,M(:)) = p';
					end
				end
			end
			
			% Reset warnings
			warning(warning_state);
			
			for b = 1:nMaps
				proposals{b}(isnan(proposals{b})) = 1e-100;
				proposals{b}(isinf(proposals{b})) = 1e-100;
			end
		end
		
		function display(self)
			fprintf('Disparity map with unary term from GlobalStereo \n');
			display@dispmap_super(self);
		end
		
		function restart(self)
			% Reset to inital solution
			init_solution(self);
		end
		function set.tol(self, tol)
			if (tol < 0)
				error('Tolerance weight must be positive');
			end
			
			self.tol = tol;
			update_energy(self);
		end
		
		function display_surfaces(self, roi)
			% "Settings"
			im0 = self.images{1};
			sz = size(self.images{1});
			tol = 1000;
			assignment = self.assignment;
			
			% Choosen a subset
			if nargin == 2
				assert(min(roi.x) > 0);
				assert(max(roi.x) <= size(im0,1));
				
				assert(min(roi.y) > 0);
				assert(max(roi.y) <= size(im0,2));
				
				[xx,yy] = meshgrid(roi.y, roi.x);
				im0 = im0(roi.x,roi.y,:);
				inds = sub2ind([sz(1) sz(2)], yy,xx);
				%
				assignment = assignment(:,inds(:));
			else
				
				[xx,yy] = meshgrid(1:size(im0,2),1:size(im0,1));
			end
			impoints = [xx(:)'; yy(:)'];
	
			f = 3740;
			imnr = 0;
			b = 40*(imnr-1);
			dmin = 100;
						
			disps = -(sum(assignment(1:2,:).*impoints) + assignment(4,:))./assignment(3,:);
			
			% Triangualte
			nodes = zeros(size(im0(:,:,1)));
			nodes(:) = 1:length(nodes(:));
			
			corner1 = [nodes(1:end-1,1:end-1)];
			corner2 = [nodes(1:end-1,2:end)];
			corner3 = [nodes(2:end,1:end-1)];
			
			corner1 = corner1(:);
			corner2 = corner2(:);
			corner3 = corner3(:);
			
			disp1 = disps(corner1);
			disp1prim = -(sum(assignment(1:2,corner2).*impoints(:,corner1)) + assignment(4,corner2))./assignment(3,corner2);
			disp2 = disps(corner2);
			disp2prim = -(sum(assignment(1:2,corner1).*impoints(:,corner2)) + assignment(4,corner1))./assignment(3,corner1);
			
			cutedges = (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			disp1 = disps(corner3);
			disp1prim = -(sum(assignment(1:2,corner2).*impoints(:,corner3)) + assignment(4,corner2))./assignment(3,corner2);
			disp2 = disps(corner2);
			disp2prim = -(sum(assignment(1:2,corner3).*impoints(:,corner2)) + assignment(4,corner3))./assignment(3,corner3);
			
			cutedges = cutedges | (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			disp1 = disps(corner3);
			disp1prim = -(sum(assignment(1:2,corner1).*impoints(:,corner3)) + assignment(4,corner1))./assignment(3,corner1);
			disp2 = disps(corner1);
			disp2prim = -(sum(assignment(1:2,corner3).*impoints(:,corner1)) + assignment(4,corner3))./assignment(3,corner3);
			
			cutedges = cutedges | (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			tri = [corner1(~cutedges) corner2(~cutedges) corner3(~cutedges)];
			imgray = rgb2gray(im0);
			
			color = imgray(1:end-1,1:end-1);
			tricolor = color(~cutedges);
						
			corner1 = [nodes(1:end-1,2:end)];
			corner3 = [nodes(2:end,1:end-1)];
			corner2 = [nodes(2:end,2:end)];
			
			corner1 = corner1(:);
			corner2 = corner2(:);
			corner3 = corner3(:);
			
			disp1 = disps(corner1);
			disp1prim = -(sum(assignment(1:2,corner2).*impoints(:,corner1)) + assignment(4,corner2))./assignment(3,corner2);
			disp2 = disps(corner2);
			disp2prim = -(sum(assignment(1:2,corner1).*impoints(:,corner2)) + assignment(4,corner1))./assignment(3,corner1);
			
			cutedges = (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			disp1 = disps(corner3);
			disp1prim = -(sum(assignment(1:2,corner2).*impoints(:,corner3)) + assignment(4,corner2))./assignment(3,corner2);
			disp2 = disps(corner2);
			disp2prim = -(sum(assignment(1:2,corner3).*impoints(:,corner2)) + assignment(4,corner3))./assignment(3,corner3);
			
			cutedges = cutedges | (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			disp1 = disps(corner3);
			disp1prim = -(sum(assignment(1:2,corner1).*impoints(:,corner3)) + assignment(4,corner1))./assignment(3,corner1);
			disp2 = disps(corner1);
			disp2prim = -(sum(assignment(1:2,corner3).*impoints(:,corner1)) + assignment(4,corner3))./assignment(3,corner3);
			
			cutedges = cutedges | (disp1-disp1prim) > tol | (disp2-disp2prim) > tol;
			
			color = imgray(1:end-1,1:end-1);
			tricolor = [tricolor color(~cutedges)];
			
			tri = [tri; [corner1(~cutedges) corner2(~cutedges) corner3(~cutedges)]];
			
			depth = f*b./(-disps+dmin);
			points3D = [impoints/10;disps];
			patch('Faces',tri,'Vertices',points3D','FaceColor','b','Edgecolor','none');		
			axis equal
			title(sprintf('Energy %g \n', energy(self)));
		end
	end
	
	methods (Access = protected)
		function disps = disparitymap_from_assignment(self, assignment, points)
			% Overloaded
			if nargin < 3
				points = self.get_points();
			end
			
			% Rescale
			disps = disparitymap_from_assignment@dispmap_super(self, assignment, points);
			disps = (disps- self.d_min)/self.d_step;
		end
		
		function init_solution(self)
			assignment = zeros(4, prod(self.sz(1:2)) );
			assignment(3,:) = 1;
			assignment(4,:) = -self.start_disparity(:);
			
			self.assignment  = assignment;
		end
		
		function U = unary_cost(self,assignment)
			disp = self.d_step*(reshape(disparitymap_from_assignment(self,assignment),self.sz(1:2)) + self.d_min);
			
			tp = prod(self.sz(1:2));
			X = repmat(1:self.sz(2), [self.sz(1) 1]);
			Y = repmat((1:self.sz(1))', [1 self.sz(2)]);
			WC = [X(:) Y(:) ones(tp, 1) disp(:)];
			
			a = 2;
			T = WC * self.P(:,:,a);
			N = 1 ./ T(:,3);
			T(:,1) = T(:,1) .* N;
			T(:,2) = T(:,2) .* N;
			
			oobv = -1000;
			R = double(reshape(self.images{1},self.sz(1)*self.sz(2),[]));
			M = vgg_interp2(self.images{a}, T(:,1), T(:,2), 'linear', oobv);
			M = squeeze(M) - R;
			
			U = self.ephoto(M);		
		end
		
		function preprocess(self)
			Rorig = uint8(self.images{1});
			sz = size(Rorig);
			sz = sz(1:2);
			
			% Only planar terms
			self.options.planar = 0;
			colors = size(Rorig,3);
			
			if colors == 1
				Rorig = repmat(Rorig, [1 1 3]);
			end
			
			% Segment the image using mean shift
			info.segment = vgg_segment_ms(Rorig, self.options.seg_params(1), ...
				self.options.seg_params(2), self.options.seg_params(3));
			self.improve = (self.options.improve > 0);
			
			% Setting up equvivalent smoothness
			T = reshape(uint32(1:prod(sz)), sz);
			num_in = numel(self.images);
			
			% Find smoothness edges which don't cross segmentation boundaries
			nbh = [self.neighborhood.ind1 self.neighborhood.ind2]';
			EW = reshape(~any(diff(int32(info.segment(nbh))), 1), 1, []);
			EW = EW * self.options.lambda_h + ~EW * self.options.lambda_l;
			EW = EW * (num_in / ((self.options.connect==8) + 1));
			
			self.ephoto = @(F) log(2) - log(exp(sum(F .^ 2, 2)*(-1/(self.options.col_thresh*colors)))+1);
			
			self.smooth_weights = EW;
			
			
			if (self.smoothness_kernel == 2)
				self.smooth_weights = self.smooth_weights/self.tol;
				self.tol = self.tol^2;
			end
		end
		
		% LO-RANSAC functions
		function inls = rplane(self, pts, th)
			MAX_SAM = 500;
			conf = .95;
			
			len = size(pts, 1);
			max_i = 3;
			max_sam = MAX_SAM;
			no_sam = 0;
			div = repmat(-1, [3 1]);
			inls = false(len, 1);
			
			while no_sam < max_sam
				no_sam = no_sam + 1;
				sam = randperm(len);
				sam = sam(1:3);
				
				%%% compute a distance of all points to a plane given by
				%%% pts(:,sam) to dist
				N = pts(sam,:) \ div;
				dist = abs((pts * N) + 1);
				v = dist < th;
				no_i = sum(v);
				
				if max_i < no_i
					% Re-estimate plane and inliers
					N = pts(v,:) \ repmat(-1, [no_i 1]);
					dist = abs((pts * N) + 1);
					v = dist < th;
					
					if sum(v) > sum(inls)
						inls = v;
						max_i = no_i;
						max_sam = min([max_sam,nsamples(self, sum(inls), len, 3, conf)]);
					end
				end
			end
		end
		function SampleCnt = nsamples(self, ni, ptNum, pf, conf)
			q  = prod ([(ni-pf+1) : ni] ./ [(ptNum-pf+1) : ptNum]);
			
			if (1 -q) < eps
				SampleCnt = 1;
			else
				SampleCnt  = log(1 - conf) / log(1 - q);
			end
			
			if SampleCnt < 1
				SampleCnt = 1;
			end
		end
		function p = FitPlaneToPoints(self, points)
			p = zeros(4,1);
			points0 = points-repmat(mean(points,2),[1 size(points,2)]);
			
			A = zeros(3,3);
			for i = 1:size(points,2);
				A = A + points0(1:3,i)*points0(1:3,i)';
			end
			[V,D] = eig(A);
			
			p(1:3) = V(:,1);
			p(4) = -(p(1:3)'*mean(points(1:3,:),2));
		end
	end
end