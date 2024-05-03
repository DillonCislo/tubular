function K = gaussianCurvatureIntrinsic(F, L)
%GAUSSIANCURVATUREINTRINSIC Calculate the Gaussian curvature of a surface
%using only the intrinsic representation of the metric via the edge
%lengths, (i.e. no explicit representation of the surface is required)
%
%   INPUT PARAMETERS:
%
%       - F:        #Fx3 face connectivity list
%       - L:        #Ex1 edge length list
%
%   OUTPUT PARAMETERS:
%
%       - K:        #Vx1 Gaussian curvature on mesh vertices
%
%   by Dillon Cislo 08/30/2021

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply edge length list'); end

validateattributes( F, {'numeric'}, ...
    {'ncols', 3, 'positive', 'integer', 'finite', 'real'} );
validateattributes( L, {'numeric'}, ...
    {'vector', 'positive', 'finite', 'real'} );

if (size(L,2) ~= 1), L = L.'; end

% The sorted list of vertex ID's in the mesh
allVIDx = unique(F(:));
assert( isequal(allVIDx, (1:numel(allVIDx)).'), ...
    'Face connectivity list contains unused vertices' );

% A 'fake' vertex coordinate list
V = repmat(allVIDx, 1, 3) .* ones(numel(allVIDx), 3);

% A MATLAB-style representation of the triangulation
TR = triangulation(F, V);

% The edge connectivity list - it is assumed that the edge length list is
% sorted according to the ordering calculated by the 'edges' function
E = TR.edges;
if (numel(L) ~= size(E,1))
    error('Improperly sized edge length list');
end

% Vertex IDs on the boundary
bdyIDx = freeBoundary(TR);
bdyIDx = bdyIDx(:,1);

% Vertex IDs in the bulk
bulkIDx = allVIDx(~ismember(allVIDx, bdyIDx));

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, sort(E,2), 'rows' );
[~, e2IDx] = ismember( e2IDx, sort(E,2), 'rows' );
[~, e3IDx] = ismember( e3IDx, sort(E,2), 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

%--------------------------------------------------------------------------
% Calculate Gaussian Curvature
%--------------------------------------------------------------------------

% Edge lengths defined on faces
L_F = L(feIDx);

% Calculate internal angles from triangulation edge lengths ---------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);

% The internal angles
intAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );
intAng = acos(intAng);

% Combine sum internal angles around vertex 1-ring to calculate Gaussian
% curvatures --------------------------------------------------------------
K = full(sparse( F(:), 1, intAng(:), numel(allVIDx), 1 ));
K(bulkIDx) = 2 * pi - K(bulkIDx);
K(bdyIDx) = pi - K(bdyIDx);

end
