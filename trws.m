% Wrapper function for TRW-S solver
function [solution, energy, lower_bound, iterations] =  ...
			trws(kernel, unary, connectivity, q, qprim, alphas, tol, options)

assert(min(connectivity(:) > 0));
assert( max(connectivity(:) <= numel(unary)))
kernel = int32(kernel);

if (any(isnan(q(:))))
    error('q contains NaN');
end

if (any(isnan(qprim(:))))
    error('qprim contains NaN');
end

% Compile if need be
my_name = mfilename('fullpath');
my_path = fileparts(my_name);

cpp_file = ['cpp' filesep 'trws_mex.cpp'];
out_file = 'trws_mex';
extra_arguments = {['-I"' my_path '"']};
sources = {['cpp' filesep 'trw-s' filesep 'MRFEnergy.cpp'], ...
						['cpp' filesep 'trw-s' filesep 'minimize.cpp'], ...
						['cpp' filesep 'trw-s' filesep 'ordering.cpp'], ...
						['cpp' filesep 'trw-s' filesep 'treeProbabilities.cpp']};

compile(cpp_file, out_file, sources, extra_arguments)

% Solve
% Change from matlab from base 1 to base 0.
[solution, energy, lower_bound, iterations] = trws_mex(kernel, unary, connectivity-1, q, qprim, alphas, tol, options);