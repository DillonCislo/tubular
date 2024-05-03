function mu = bc_metric_intrinsic(F, V2D, L)
%BC_METRIC_INTRINSIC Calculate the Beltrami coefficient of a metric
%represented by a mesh triangulation using only the edge lengths of the
%mapping (i.e. no explicit mapping is required)
%
%   INPUT PARAMETERS
%
%       - F:        #Fx3 face connectivity list
%       - V2D:      #Vx2 2D vertex coordinate list
%       - L:        #Ex1 edge length list
%
%   OUTPUT PARAMETERS:
%
%       - mu:    #Fx1 complex Beltrami coefficient
%
%   by Dillon Cislo 08/31/2021

% NOTE: Further input validation is performed in 'inducedMetricIntrinsic'
if (nargin < 1), error('Please suppply face connectivity list'); end
if (nargin < 2), error('Please supply initial vertex coordinate list'); end
if (nargin < 3), error('Please supply edge length list'); end

% Calculate the metric tensor on faces
g = inducedMetricIntrinsic(F, V2D, L);

% Calculate the Beltrami coefficient
mu = cellfun( @(x) ...
    (x(1,1)-x(2,2)+2i*x(1,2)) ./ (trace(x)+2*sqrt(det(x))), g );

end
