function [U, Uall] = laplacian_smooth_edges(V, F, L_method, b, ...
    lambda, method, S, max_iter)
%LAPLACIAN_SMOOTH_EDGES Smooth edge-based scalar quantities using Laplacian
%smoothing based on Crouzeix-Raviart non-conforming finite elements
%   
%   INPUT PARAMETERS:
%
%       - V:            #VxVD vertex coordinate list
%
%       - F:            #Fx3 face connectivity list
%
%       - L_method:     Method for Laplacian:
%                           - 'uniform'
%                           - 'crouzeix-raviart ' (default)
%
%       - b:            #Bx2 vector of indices of fixed edges
%
%       - lambda:       Diffusion speed parameter (0.1)
%
%       - method:       Solution method:
%                           - 'implicit' (default)
%                           - 'explicit'
%
%       - S:            #ExSD edge-based scalar field to smooth 
%                       (default is edge lengths)
%
%       - max_iter:     Maximum number of smoothing iterations (1000)
%
%   OUTPUT PARAMETERS:
%
%       - U:            #ExSD smoothed edge-based scalar field
%
%       - Uall:         #ExSDxiters list of smoothed values for each
%                       iteration
%
%   by Dillon Cislo 2024/02/20

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------
validateattributes(V, {'numeric'}, {'2d', 'finite', 'real'});
numV = size(V,1); dimV = size(V,2);
assert((dimV == 2) || (dimV == 3), ['Input triangulation must be ' ...
    '2D or 3D']);

validateattributes(F, {'numeric'}, {'2d', 'integer', 'positive', ...
    'finite', 'real', 'ncols', 3, '<=', numV});
% numF = size(F,1);

% Construct edge list. NOTE: This method appears to produce identical
% results to MATLAB's 'edges(triangulation(F,V))'
allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
[E, ~, ~] = unique(sort(allE, 2),'rows');
numE = size(E,1);

if (~exist('L_method', 'var') || isempty(L_method))
    L_method = 'crouzeix-raviart';
else
    validateattributes(L_method, {'char'}, {'vector'});
    assert(ismember(L_method, {'uniform', 'crouzeix-raviart'}), ...
        'Invalid Laplacian construction method');
end

if (~exist('b', 'var') || isempty(b))
    b = [];
else
    validateattibutes(b, {'numeric'}, {'vector', 'integer', ...
        'positive', 'finite', 'real', '<=', numE});
end

if (~exist('lambda', 'var') || isempty(lambda))
    lambda = 0.1;
else
    validateattributes(lambda, {'numeric'}, {'scalar', ...
        'nonnegative', 'finite', 'real'});
end

if (~exist('method', 'var') || isempty(method))
    method = 'implicit';
else
    validateattributes(method, {'char'}, {'vector'});
    assert(ismember(method, {'implicit', 'explicit'}), ...
        'Invalid diffusion update solution method');
end

if (~exist('S', 'var') || isempty(S))
    S = edge_lengths(V, E);
else
    validateattributes(S, {'numeric'}, {'2d', 'finite', 'real', ...
        'nrows', numE});
end
% dimS = size(S,2);

if (~exist('max_iter', 'var') || isempty(max_iter))
    max_iter = 1000;
else
    validateattributes(max_iter, {'numeric'}, {'scalar', 'integer', ...
        'nonnegative', 'finite', 'real'});
end

% Handle the trivial case
if ((lambda == 0) || (max_iter == 0))
    U = S; Uall = S;
    return;
end

% TODO: This should be exposed to the user
h = avgedge(V,F);
if(~exist('tol','var'))
    tol = 0.001;
end

%--------------------------------------------------------------------------
% RUN SMOOTHING
%--------------------------------------------------------------------------

I = speye(numE, numE);

% Right now we only compute the Laplacian once regardless of method. Is
% there a way to construct the Crouzeix-Raviart Laplacian using intrinsic
% edge lengths?
L = -crouzeix_raviart_cotmatrix(V, F);
if strcmpi(L_method, 'uniform')
    error('Uniform edge Laplacian smoothing is not working yet');
    L(abs(L(:)) > 0) = 1; % This is really a lazy way to do this
end

% Factorization and symmtery flag used by 'min_quad_with_fixed'
P = [];
% sym = [];

iter = 0;
U = S;
U_prev = S;
if nargout >= 2
    Uall = [];
end

while( iter < max_iter && (iter == 0 || max(abs(U(:)-U_prev(:)))>tol*h) )
    
    U_prev = U;
    
    switch method
        
        case 'implicit'
            Q = (I-lambda*L);
            % could prefactor Q for 'uniform' case
            for d = 1:size(S,2)
                [U(:,d),P] = min_quad_with_fixed(Q*0.5,-U(:,d),b,S(b,d),[],[],P);
            end
            
        case 'explicit'
            Q = (I+lambda*L);
            U = Q * U;
            % enforce boundary
            U(b,:) = S(b,:);
            
        otherwise
            error(['' method ' is not a supported smoothing method']);
            
    end
    
    if (nargout >= 2), Uall = cat(3, Uall, U); end
    
    iter = iter + 1;
    
end

end

