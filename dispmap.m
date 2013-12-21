% Johannes Ulén and Carl Olsson 2013
%
classdef dispmap < handle
	properties
		smoothness_kernel
		tol;
		unary_weight;
		assignment;

		% For TRW-S
		maxiter = 1000;
	end
	properties (SetAccess = private)
		sz;
		images;
		smooth;
		neighborhood;
		disparities; % Possible disparities
		ncc;
		stored_energy; % Keep track on energy
	end
	methods
		%% Constructor
		function self = dispmap(images,disparities, kernel, unary_weight,tol)
			self.disparities= disparities;
			self.images = images;
			self.sz = size(images{1});
			self.sz = self.sz(1:2);

			self.tol = tol;
			self.unary_weight = unary_weight;
			self.smoothness_kernel = kernel;

			% Initilization
			construct_neighborhood(self);

			% Pre-calculate NCC at each disparity level and save it in NCC
			compute_ncc(self,2);

			% Intitlize a solution which is based on the NCC volume
			init_solution(self);
		end

		%% Set functions
		function set.tol(self, tol)
			if (tol < 0)
				error('Tolerance weight must be positive');
			end

			self.tol = tol;
			update_energy(self);
		end
		function set.unary_weight(self, weight)
			if (weight < 0)
				error('Unary weight must be positive');
			end

			self.unary_weight = weight;
			update_energy(self);
		end
		function set.smoothness_kernel(self, kernel)
			switch (kernel)
				case 1
					self.smooth = @(F) min(sqrt((F).^2), self.tol);
				case 2
					self.smooth = @(F) min((F).^2,self.tol);
				otherwise
					error('Unkown kernel type');
			end

			self.smoothness_kernel = kernel;
			update_energy(self);
		end
		function set.assignment(self, assignment)
			self.assignment  = assignment;
			update_energy(self);
		end

		%% Normal methods
        function restart(self)
           % Reset to inital solution 
           init_solution(self);
        end
		function E = energy(self)
			E = self.stored_energy;
		end

		function [e,lb] = binary_fusion(self, proposal)
			% Perform one binary fusion

			% Input check
			if (~all(size(proposal) == size(self.assignment)))
				error('Binary fusion: Proposals is of wrong size');
			end

			% Unary cost for all combinations
			[E00,E01, E10,E11] = pairwise_cost(self, self.assignment, proposal);

			% Unary cost for the two choses
			points = get_points(self);
			U0 = self.unary_weight*unary_cost(self, self.assignment);
			U1 = self.unary_weight*unary_cost(self, proposal);
			connectivity = uint32([self.neighborhood.ind1'; self.neighborhood.ind2']);

			%
			% Call Roof duality solver
			[labelling, e, lb] = rd(U0, U1, E00, E01, E10, E11, connectivity);

			% Set
			self.assignment(:,labelling == 1) = proposal(:,labelling == 1);
		end
		function [e,lb, iterations] = multi_fusion(self,proposal_cell)
			assert(isa(proposal_cell,'cell'));

			% Keep the pixel choose.
			proposal_cell{end+1} = self.assignment;

			%Data cost
			points = get_points(self);
			unary = zeros(length(proposal_cell),size(points,2));

			parfor i = 1:length(proposal_cell);
				unary(i,:) = self.unary_weight*unary_cost(self, proposal_cell{i});
			end

			connectivity = uint32([self.neighborhood.ind1'; self.neighborhood.ind2']);
           
			% Calculate all p's and q's
			ind1 = self.neighborhood.ind1;
			ind2 = self.neighborhood.ind2;	
	
            alphas = ones(numel(ind1),1);
            kernel = int32(self.smoothness_kernel);
            
			points = get_points(self);
			q = zeros(length(proposal_cell),length(ind1));
			qprim = zeros(length(proposal_cell),length(ind2));

			parfor p_id = 1:length(proposal_cell)
				q(p_id,:) = disparitymap_from_assignment(self, proposal_cell{p_id}(:,ind2), points(:,ind2));
				qprim(p_id,:) = disparitymap_from_assignment(self, proposal_cell{p_id}(:,ind1), points(:,ind2));
            end
            
			options_struct.maxiter = self.maxiter;
			[L,e,lb, iterations] = trws(kernel, unary, connectivity, q, qprim, alphas, self.tol, options_struct);

			%Extract solution
			assignments = zeros(size(proposal_cell{1}));
			for i = 1:length(proposal_cell)
				assignments(:,L == i) = proposal_cell{i}(:,L==i);
            end
            
            self.assignment = assignments;
		end

		function im = current_dispmap(self)
			im = reshape(disparitymap_from_assignment(self, self.assignment), [self.sz]);
		end
		function display_current_dispmap(self)
			imagesc(self.current_dispmap());

			colormap gray(256);
			caxis([min(self.disparities) max(self.disparities)]);
			title(sprintf('Solution energy: %g \n', self.energy()));
			axis equal;
		end
		function display(self)
			fprintf('-- \n');
			fprintf('Disparity object \n');
			fprintf('-- \n');
			fprintf('Energy of current solution: %g. \n', self.energy());
			fprintf('Image pair size: (%d,%d) \n', self.sz(1), self.sz(2));
			fprintf('Disparity levels %d in range [%g, %g]. \n', numel(self.disparities), max(self.disparities), min(self.disparities));
			fprintf('Settings: \n');
			fprintf('Unary weight         : %g \n', self.unary_weight);
			fprintf('Tolerance            : %g \n', self.tol);
			fprintf('Smoothness kernel    : %d \n', self.smoothness_kernel)
			fprintf('Maximum iterations    : %d \n', self.maxiter)
			fprintf('-- \n');

			display_current_dispmap(self);
		end

	end

	methods (Access = private)
		%% NCC methods
		function compute_ncc(self, patchsize)
			% Start multi-core
			if (matlabpool('size') == 0)
				matlabpool open
			end

			% Readability
			im0 = self.images{1};
			im1 = self.images{2};
			d   = self.disparities;

			ncc = zeros(size(im0,1),size(im0,2),length(d));

			%compute ncc-mean and ncc-norm of right image
			meanpatch = ones(2*patchsize+1)./sum(sum(ones(2*patchsize+1)))/3;
			patch = ones(2*patchsize+1);
			Rright = double(im0(:,:,1));
			Gright = double(im0(:,:,2));
			Bright = double(im0(:,:,3));
			mean_right = conv2(Rright,meanpatch,'same')+conv2(Gright,meanpatch,'same')+conv2(Bright,meanpatch,'same');

			term1R = conv2(Rright.^2,patch,'same');
			term1G = conv2(Gright.^2,patch,'same');
			term1B = conv2(Bright.^2,patch,'same');

			term2R = mean_right.*conv2(Rright,patch,'same');
			term2G = mean_right.*conv2(Gright,patch,'same');
			term2B = mean_right.*conv2(Bright,patch,'same');

			term4 = sum(patch(:))*3*mean_right.^2;
			norm_right = sqrt(term1R+term1G+term1B-2*(term2R+term2G+term2B)+term4);

			parfor i = 1:length(d);
				H = eye(3);
				H(1,3) = d(i);
				tform = maketform('projective',H');
				imtr = imtransform(im1,tform,'xdata',[1 size(im0,2)],'ydata',[1 size(im0,1)],'size',size(im0));
				bnd_im = imtransform(ones(size(im1(:,:,1))),tform,'xdata',[1 size(im0,2)],'ydata',[1 size(im0,1)],'size',size(im0(:,:,1)),'fill',0);

				Rtr = double(imtr(:,:,1));
				Gtr = double(imtr(:,:,2));
				Btr = double(imtr(:,:,3));

				mean_tr = conv2(Rtr,meanpatch,'same')+conv2(Gtr,meanpatch,'same')+conv2(Btr,meanpatch,'same');

				%Compute ncc-mean and ncc- norm of transformed-left image
				term1R = conv2(Rtr.^2,patch,'same');
				term1G = conv2(Gtr.^2,patch,'same');
				term1B = conv2(Btr.^2,patch,'same');

				term2R = mean_tr.*conv2(Rtr,patch,'same');
				term2G = mean_tr.*conv2(Gtr,patch,'same');
				term2B = mean_tr.*conv2(Btr,patch,'same');

				term4 = sum(patch(:))*3*mean_tr.^2;
				norm_tr = sqrt(term1R+term1G+term1B-2*(term2R+term2G+term2B)+term4);

				%Compute ncc
				term1R = conv2(Rright.*Rtr,patch,'same');
				term1G = conv2(Gright.*Gtr,patch,'same');
				term1B = conv2(Bright.*Btr,patch,'same');

				term2R = mean_right.*conv2(Rtr,patch,'same');
				term2G = mean_right.*conv2(Gtr,patch,'same');
				term2B = mean_right.*conv2(Btr,patch,'same');

				term3R = mean_tr.*conv2(Rright,patch,'same');
				term3G = mean_tr.*conv2(Gright,patch,'same');
				term3B = mean_tr.*conv2(Bright,patch,'same');

				term4 = sum(patch(:))*3*mean_tr.*mean_right;
				ncci = term1R+term1G+term1B-(term2R+term2G+term2B)-(term3R+term3G+term3B)+term4;
				ncci = ncci./norm_right./norm_tr;

				ncci(~isfinite(ncci)) = 0;
				ncci(~(bnd_im>=1-1e-8)) = 0;
				ncc(:,:,i) = real(ncci);
			end

			% parfor does not allow to write to self.ncc directly
			% Data is not copied with this operation
			self.ncc = ncc;
		end
		function init_solution(self)
			d = self.disparities;
			[y2,t2] = max(self.ncc,[],3);
			ncc_size = size(self.ncc);

			okdepth = (t2 < ncc_size(3) & t2 > 1);

			d2 = d(t2);
			t1 = zeros(size(t2));
			t1(okdepth) = t2(okdepth)-1;
			t1(~okdepth) = t2(~okdepth);
			d1 = d(t1);

			t3 = zeros(size(t2));
			t3(okdepth) = t2(okdepth)+1;
			t3(~okdepth) = t2(~okdepth);
			d3 = d(t3);

			% sample at t1 and t2
			[col,row] = meshgrid(1:size(self.ncc,2), 1:size(self.ncc,1));
			y1 = self.ncc(row+(col-1)*ncc_size(1)+(t1-1).*ncc_size(1)*ncc_size(2));
			y3 = self.ncc(row+(col-1)*ncc_size(1)+(t3-1).*ncc_size(1)*ncc_size(2));

			% interpolate local minima
			a = y1./(d1-d2)./(d1-d3);
			b = y2./(d2-d1)./(d2-d3);
			c = y3./(d3-d1)./(d3-d2);

			%polynom: r*d^2+p*d+q
			r = a+b+c;
			p = -(a.*(d2+d3)+b.*(d1+d3)+c.*(d1+d2));
			q = a.*d2.*d3 + b.*d1.*d3 + c.*d1.*d2;

			% local maxima
			best_disp = -p./r/2;
			best_disp(~okdepth) = d2(~okdepth);

			assignment = zeros(4, prod(self.sz(1:2)) );
			assignment(3,:) = 1;
			assignment(4,:) = -best_disp(:);
            
            self.assignment = assignment;
		end
		function nccs = sample_ncc_from_disp(self, new_pixel_disps)
			largeval = 1e6;
			new_pixel_disps = reshape(new_pixel_disps,self.sz(1:2));
			d = self.disparities;

			closest_depth_ind = ones(size(self.ncc(:,:,1)));
			smallest_dist = abs(new_pixel_disps-d(1));
			y2 = ones(size(self.ncc(:,:,1)));
			for i = 1:length(d);
				new_dists = abs(new_pixel_disps - d(i));
				closest_depth_ind(new_dists <= smallest_dist) = i;
				ncci = self.ncc(:,:,i);
				y2(new_dists <= smallest_dist) = ncci(new_dists <= smallest_dist);
				smallest_dist(new_dists <= smallest_dist) = new_dists(new_dists <= smallest_dist);
			end

			t2 = closest_depth_ind;

			okdepth = (t2 < size(self.ncc,3) & t2 > 1);
			good_disp = new_pixel_disps <= max(d) & new_pixel_disps >= min(d);

			d2 = d(t2);
			t1 = zeros(size(t2));
			t1(okdepth) = t2(okdepth)-1;
			t1(~okdepth) = t2(~okdepth);
			d1 = d(t1);

			t3 = zeros(size(t2));
			t3(okdepth) = t2(okdepth)+1;
			t3(~okdepth) = t2(~okdepth);
			d3 = d(t3);

			% sample at t1 and t2
			[col,row] = meshgrid(1:size(self.ncc,2), 1:size(self.ncc,1));
			y1 = self.ncc(row+(col-1)*size(self.ncc,1)+(t1-1).*size(self.ncc,1)*size(self.ncc,2));
			y3 = self.ncc(row+(col-1)*size(self.ncc,1)+(t3-1).*size(self.ncc,1)*size(self.ncc,2));

			% local minimum
			a = y1./(d1-d2)./(d1-d3);
			b = y2./(d2-d1)./(d2-d3);
			c = y3./(d3-d1)./(d3-d2);

			%polynom: r*d^2+p*d+q
			r = a+b+c;
			p = -(a.*(d2+d3)+b.*(d1+d3)+c.*(d1+d2));
			q = a.*d2.*d3 + b.*d1.*d3 + c.*d1.*d2;

			nccs = r.*new_pixel_disps.^2+p.*new_pixel_disps+q;
			ncci = self.ncc(:,:,1);
			nccs(t2 == 1) = ncci(t2 == 1);
			ncci = self.ncc(:,:,end);
			nccs(t2 == length(d)) = ncci(t2 == length(d));
			nccs(~good_disp) = -largeval;
		end
		function U = unary_cost(self, assignment)
			disps = disparitymap_from_assignment(self,assignment);

			%sample ncc at depths
			nccs = sample_ncc_from_disp(self, disps);

			% Unary cost
			U = 1-nccs(:);
		end
		% General methods

		function [E00,E01, E10,E11] = pairwise_cost(self, assignment, proposals)
			% Smoothcost
			points = get_points(self);

			ind1 = self.neighborhood.ind1;
			ind2 = self.neighborhood.ind2;

			q = disparitymap_from_assignment(self, assignment(:,ind2), points(:,ind2));
			qprim = disparitymap_from_assignment(self, assignment(:,ind1), points(:,ind2));

			E00 = self.smooth(q-qprim);

			if nargout > 1
				%new points
				new_q = disparitymap_from_assignment(self, proposals(:,ind2), points(:,ind2));
				new_qprim = disparitymap_from_assignment(self, proposals(:,ind1), points(:,ind2));

				%Smoothness cost p = alpha, q = alpha => V(p,q)
				E11 = self.smooth(new_q - new_qprim);

				%Smoothness cost p = alpha (=0), q = beta (=1)
				E10 = self.smooth(q - new_qprim);

				%Smoothness cost p = beta (=1), q = alpha (=0)
				E01 = self.smooth(new_q - qprim);
			end
		end
		function update_energy(self)
			% Return inf is no assigment is set yet
			if (isempty(self.assignment))
				self.stored_energy = inf;
			else
				% Call everytime assigments is update
				U = unary_cost(self, self.assignment);
				P = pairwise_cost(self, self.assignment);

				self.stored_energy = self.unary_weight * sum(U(:)) + sum(P(:));
			end
		end
		function points = get_points(self)
			[xx,yy] = meshgrid(1:self.sz(2),1:self.sz(1));
			points = [xx(:)'; yy(:)'];
		end
		function construct_neighborhood(self)
			%Generate neighborhood
			nodenr = zeros(self.sz(1), self.sz(2) );
			nodenr(:) = 1:length(nodenr(:));

			%Vertical edges
			start = nodenr(1:end-1,:);
			finish = nodenr(2:end,:);
			self.neighborhood.ind1 = [start(:); finish(:)];
			self.neighborhood.ind2 = [finish(:); start(:)];

			%Horisotal edges
			start = nodenr(:,1:end-1);
			finish = nodenr(:,2:end);
			self.neighborhood.ind1 = [self.neighborhood.ind1; start(:); finish(:)];
			self.neighborhood.ind2 = [self.neighborhood.ind2; finish(:); start(:)];

			%remove edgesl from or to background node
			keep_edges = (self.neighborhood.ind2 ~= 0) & (self.neighborhood.ind1 ~= 0);
			self.neighborhood.ind1 = self.neighborhood.ind1(keep_edges);
			self.neighborhood.ind2 = self.neighborhood.ind2(keep_edges);

			self.neighborhood.nodenr = nodenr;
		end
		function set_disparity(self, assignment)
			self.assignment(1:2,:) = 0;
			self.assignment(3,:) = 1;
			self.assignment(4,:) = -assignment(:);
		end
		function init_assigments(self)
			self.assignment = zeros(4, prod(self.sz(1:2)) );
			self.assignment(3,:) = 1;
		end

		% Remark:
		% It is possible to make this implementation with assigments(3,i) = 1 \forall i and hence removing it.
		% However keeping all four variables for each voxel leads to easier to read code.

		% Convert given assigment to a disparty map.
		function disps = disparitymap_from_assignment(self, assignment, points)
			if nargin < 3
				points = self.get_points();
            end
            
            if (any(any(assignment(3,:) == 0)))
                error('Infinite disparity');
            end

			disps = -(sum(assignment(1:2,:).*points)+assignment(4,:))./assignment(3,:);
		end

	end
end