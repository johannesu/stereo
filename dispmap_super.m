% Johannes Ulén and Carl Olsson 2013
%
classdef dispmap_super < handle
	properties
		smoothness_kernel
		assignment;
		
		 % For TRW-S or Exhaustive binary fusion
		maxiter = 1000;
		max_relgap = 1e-4; 
		
		% Use RD-Improve
		improve = false; 
	end
	properties (SetAccess = protected)
		sz;
		images;
		neighborhood;
		stored_energy; % Keep track on energy
		smooth_weights;
	end
	methods
		%% Constructor
		function self = dispmap_super(images,kernel)
			self.images = images;
			self.sz = size(images{1});
			self.sz = self.sz(1:2);
			
			self.smoothness_kernel = kernel;

			% Initilization
			construct_neighborhood(self);
			
			% Default all smooth weights are 1.
			self.smooth_weights = ones(1,numel(self.neighborhood.ind1));
		end
	
		%% Set functions
		function set.max_relgap(self, max_relgap)
			if (max_relgap < 0)
				error('Maximum relative gap must be non-negative');
			end
			
			self.max_relgap = max_relgap;
		end
		function set.improve(self, improve)
				self.improve = logical(improve);
		end
		function set.assignment(self, assignment)
			self.assignment  = assignment;
			update_energy(self);
		end	
		function set.smoothness_kernel(self, kernel)
			self.smoothness_kernel = kernel;
			update_energy(self);
		end
		function E = energy(self)
			E = self.stored_energy;
		end

		function [e,lb, num_unlabelled] = binary_fusion(self, proposal)
			% Perform one binary fusion

			% Input check
			if (~all(size(proposal) == size(self.assignment)))
				error('Binary fusion: Proposals is of wrong size');
			end

			% Unary cost for all combinations
			[E00,E01, E10,E11] = all_pairwise_costs(self, self.assignment, proposal);
    			
			% Unary cost for the two choses
			U0 = unary_cost(self, self.assignment);
			U1 = unary_cost(self, proposal);
			connectivity = uint32([self.neighborhood.ind1'; self.neighborhood.ind2']);

			% Call Roof duality solver
			rd_options.improve = self.improve;
			[labelling, e, lb, num_unlabelled] =  ...
				rd(U0, U1, E00, E01, E10, E11, connectivity, rd_options);
                         
			% Set
			self.assignment(:,labelling == 1) = proposal(:,labelling == 1);
		end
		function number_of_iterations = binary_fuse_until_convergence(self, proposal_cell, show_steps)
			% Fuse binary proposals until either max iterations is reached or no proposals can improve the result
			
			if nargin < 3
				show_steps = false;
			end
			
			if ~(isa(proposal_cell,'cell'));
				error('Input proposals should be given in cell array.')
			end
			
			% Iteration order
			number_of_random_ids = self.maxiter*5;
			ids = [1:length(proposal_cell)' randi([1 length(proposal_cell)],number_of_random_ids,1)'];
			ids([diff(ids) == 0]) = 0; % Removing repating proposal number
			ids(ids < 1) = [];
			ids(ids > length(proposal_cell)) = [];
			
			propsals_fun = @(n) proposal_cell{ids(n)};
			iter = 0;
			E = energy(self);
			
			visited_propsals = zeros(numel(proposal_cell),1);
			fused_propsals_in_order = {};
			
			for iter = 1:self.maxiter
				% If random sets runs out
				if iter > number_of_random_ids
					ids = [ids ids];
				end
				
				iter = iter+1;
				
				% No proposals have fused since last visist
				% hence no need to run this again.
				if ( visited_propsals(ids(iter)))
					continue;
				end
				
				P = propsals_fun(iter);
				
				% Fuse
				binary_fusion(self, P);
				E(end+1) = energy(self);
				
				if (show_steps)
					display_current_dispmap(self);
					drawnow();
				end
				
				% If we tried ALL propsals without improvement we quit
				if (iter > 1)
					if (E(end-1) ~= E(end))
						visited_propsals(:) = 0;
					else
						visited_propsals(ids(iter)) = 1;
					end
				else
					visited_propsals(ids(1)) = 1;
				end
				
				if (all(visited_propsals))
					break;
				end
			end
			
			number_of_iterations = length(E);
		end
		function [e,lb, iterations] = simultaneous_fusion(self,proposal_cell)
			if ~(isa(proposal_cell,'cell'));
				error('Input proposals should be given in cell array.')
			end
			% Keep the pixel choose.
			proposal_cell{end+1} = self.assignment;

			%Data cost
			points = get_points(self);
			unary = zeros(length(proposal_cell),size(points,2));

			parfor i = 1:length(proposal_cell);
				unary(i,:) = unary_cost(self, proposal_cell{i});
			end

			connectivity = uint32([self.neighborhood.ind1'; self.neighborhood.ind2']);

			% Calculate all p's and q's
			ind1 = self.neighborhood.ind1;
			ind2 = self.neighborhood.ind2;

			kernel = int32(self.smoothness_kernel);

			points = get_points(self);
			q = zeros(length(proposal_cell),length(ind1));
			qprim = zeros(length(proposal_cell),length(ind2));
			
			parfor p_id = 1:length(proposal_cell)
				q(p_id,:) = disparitymap_from_assignment(self, proposal_cell{p_id}(:,ind2), points(:,ind2));
				qprim(p_id,:) = disparitymap_from_assignment(self, proposal_cell{p_id}(:,ind1), points(:,ind2));
			end


			options_struct.maxiter = self.maxiter;
			options_struct.max_relgap = self.max_relgap;
			alphas = self.smooth_weights(:);
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
			title(sprintf('Solution energy: %g \n', self.energy()));
			axis equal;
		end
		function display(self)
			fprintf('Energy of current solution: %g. \n', self.energy());
			fprintf('Image pair size: (%d,%d) \n', self.sz(1), self.sz(2));
			fprintf('Settings: \n');
			fprintf('Smoothness kernel    : %d \n', self.smoothness_kernel)
			fprintf('Maximum iterations   : %d \n', self.maxiter)
			fprintf('Max relative gap     : %g \n', self.max_relgap);
			fprintf('RD-Improve			  : %d \n', int32(self.improve));
			display_current_dispmap(self);
		end
	end

	methods (Access = protected)
		function U = unary_cost(self, assignment)
			error('Overload unary_cost \n');
		end
		function c = pairwise_cost(self, p,q)
			switch (self.smoothness_kernel)
				case 1
					c = self.smooth_weights.* min( abs( (p-q)), self.tol);
				case 2
					c= self.smooth_weights.* min( ((p-q)).^2, self.tol);
				otherwise
					error('Unkown kernel type');
			end
		end
		function [E00,E01, E10,E11] = all_pairwise_costs(self, assignment, proposals)
			% Smoothcost
			points = get_points(self);

			ind1 = self.neighborhood.ind1;
			ind2 = self.neighborhood.ind2;

			q = disparitymap_from_assignment(self, assignment(:,ind2), points(:,ind2));
			qprim = disparitymap_from_assignment(self, assignment(:,ind1), points(:,ind2));

			E00 = pairwise_cost(self,q,qprim);

			if nargout > 1
				%new points
				new_q = disparitymap_from_assignment(self, proposals(:,ind2), points(:,ind2));
				new_qprim = disparitymap_from_assignment(self, proposals(:,ind1), points(:,ind2));

				%Smoothness cost p = alpha, q = alpha => V(p,q)
				E11 =  pairwise_cost(self,new_q, new_qprim);

				%Smoothness cost p = alpha (=0), q = beta (=1)
				E10 =  pairwise_cost(self,q, new_qprim);

				%Smoothness cost p = beta (=1), q = alpha (=0)
				E01 =  pairwise_cost(self,new_q, qprim);
			end
		end
		function update_energy(self)
			% Return inf if no assigment is set yet
			if (isempty(self.assignment))
				self.stored_energy = inf;
			else
				% Call everytime assigments is update
				U = unary_cost(self, self.assignment);
				P = all_pairwise_costs(self, self.assignment);

				self.stored_energy =  sum(U(:)) + sum(P(:));
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
		% However keeping all four variables for each voxel leads more straightforward code.

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