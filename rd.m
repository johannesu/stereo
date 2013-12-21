% Wrapper function for the roof duality solver

function [solution, energy, lower_bound,num_unlabelled] = rd(U0,U1, E00, E01, E10, E11, connectivity, options)

assert(min(connectivity(:) > 0));
assert( max(connectivity(:)) <= numel(U0) );

% Compile if need be
cpp_file = ['cpp' filesep 'rd_mex.cpp'];
out_file = 'rd_mex';
extra_arguments = {};
sources = {['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO.cpp'], ...
			['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_extra.cpp'], ...
			['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_maxflow.cpp'], ...
			['cpp' filesep 'QPBO-v1.3.src' filesep 'QPBO_postprocessing.cpp']};
        
compile(cpp_file, out_file, sources, extra_arguments);

% Solve
% Change from matlab from base 1 to base 0.
[solution, energy, lower_bound, num_unlabelled] = rd_mex(U0,U1, E00, E01, E10, E11, connectivity-1, options);