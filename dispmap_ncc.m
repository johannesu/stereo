% Johannes Ulén and Carl Olsson 2013
%
% Example on how to create a disparity map
% Inherit all functinally except the intilization and unary cost from dispmap_super
classdef dispmap_ncc < dispmap_super

	properties (SetAccess = protected)
		ncc;
		disparities; % Disparities to sample NCC at
		unary_weight;
		smooth;
		tol;
	end

	methods
		function self = dispmap_ncc(images,disparities, kernel, unary_weight,tol)
			 % Superclass constructor
			 self = self@dispmap_super(images, kernel);
			 self.disparities = disparities;
			 self.unary_weight = unary_weight;
			 self.smoothness_kernel = kernel;
			 self.tol = tol;
			 % Special for normalized cross correlation
			 compute_ncc(self,2);

			 % Intitlize a solution which is based on the NCC volume
			 init_solution(self);
		end

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

		function proposal = generate_new_plane_RANSAC(self, x, y, r)
			% Generate a plane by samlping the NCC volume around x,y
			points = get_points(self);
			best_disp = best_disp_from_ncc(self);

			% Find points within raidus r
			dist = points -  ...
				[repmat(x,[1 size(points,2)]); repmat(y,[1 size(points,2)])];

			ids = sqrt(dist(1,:).^2 + dist(2,:).^2) < r;

			% Approximated from NCC sampling
			points3D = [points(:,ids); best_disp(ids)];

			% Fit a plane to points inside radius r
			p = fit_plane_to_points(self, points3D);

			proposal =repmat(p, [1 size(self.assignment,2)]);
		end
		function p = fit_plane_to_points(self, points)
			c= mean(points,2);
			n = size(points,2);

			cost_func = -(points-repmat(c,[1 n]))';
			p = zeros(4,1);

			if (self.smoothness_kernel == 1)
				w = ones(size(cost_func,1),1);
				% Iterative reweighted least squares
				% To calculate normal direction of plane
				for irls_iteration = 1:20
					w = repmat(w,[1 3]);
					% Solve using SVD
					[~,~,V] = svd(w.*cost_func,'econ');
					p(1:3) = V(:,end);
					w = sqrt(abs(cost_func*V(:,end)));
				end
			elseif (self.smoothness_kernel == 2)
				[~,~,V] = svd(cost_func,'econ');
				p(1:3) = V(:,end);
			end

			p(4) = -(p(1:3)'*mean(points(1:3,:),2));
			p = p/p(3);
		end
		function display(self)
				fprintf('Disparity map with normalized cross correlation unary term \n');
				fprintf('Disparity levels %d in range [%g, %g]. \n', numel(self.disparities), max(self.disparities), min(self.disparities));
				display@dispmap_super(self);
				fprintf('Unary weight         : %g \n', self.unary_weight);
			  fprintf('Tolerance            : %g \n', self.tol);
		end
    function restart(self)
       % Reset to inital solution
       init_solution(self);
    end
	end
	methods (Access = protected)
		%% NCC methods
		function U = unary_cost(self, assignment)
			disps = disparitymap_from_assignment(self,assignment);

			%sample ncc at depths
			nccs = sample_ncc_from_disp(self, disps);

			% Unary cost
			U = self.unary_weight*(1-nccs(:));
		end
		function compute_ncc(self, patchsize)
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
				% Move image according to disparit d(i);
				bnd_im = zeros(self.sz);
				bnd_im(:,round(d(i)+1):end) = 1;
				y_span = ceil(d(i)+1):self.sz(2);

				[X,Y] = meshgrid(linspace(1,self.sz(2)-d(i),numel(y_span)), ...
					1:self.sz(1));
				imtr = zeros([self.sz 3]);
				for dim = 1:3
					imtr(:,y_span,dim) = interp2(im1(:,:,dim),X,Y);
				end

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
			best_disp = best_disp_from_ncc(self);

			assignment = zeros(4, prod(self.sz(1:2)) );
			assignment(3,:) = 1;
			assignment(4,:) = -best_disp(:);

			self.assignment = assignment;
		end
		function best_disp = best_disp_from_ncc(self)
			% Best assigments based on ncc volume
			% all normal direction initilized to be fronto-parallel
			d = self.disparities;
			[y2,t2] = max(self.ncc,[],3);
			ncc_size = size(self.ncc);

			okdepth = (t2 < ncc_size(3) & t2 > 1);
			[r,p,q, d2] = interpolate_ncc(self, t2, y2, okdepth);

			% local maxima
			best_disp = -p./r/2;
			best_disp(~okdepth) = d2(~okdepth);
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
			[r,p,q] = interpolate_ncc(self, t2, y2, okdepth);

			nccs = r.*new_pixel_disps.^2+p.*new_pixel_disps+q;
			ncci = self.ncc(:,:,1);
			nccs(t2 == 1) = ncci(t2 == 1);
			ncci = self.ncc(:,:,end);
			nccs(t2 == length(d)) = ncci(t2 == length(d));
			nccs(~good_disp) = -largeval;
		end
		function [r,p,q, d2] = interpolate_ncc(self, t2, y2, okdepth)
			d2 = self.disparities(t2);
			t1 = zeros(size(t2));
			t1(okdepth) = t2(okdepth)-1;
			t1(~okdepth) = t2(~okdepth);
			d1 = self.disparities(t1);

			t3 = zeros(size(t2));
			t3(okdepth) = t2(okdepth)+1;
			t3(~okdepth) = t2(~okdepth);
			d3 = self.disparities(t3);

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
		end
	end
end