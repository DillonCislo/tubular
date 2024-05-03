function [compErr, allGMPC, RE] = discreteFormCompatibility(F, L, Phi)
%DISCRETEFORMCOMPATIBILITY Determine if a set of edge lengths and bending
%angles are compatible, in the sense that they may correspond to the
%discrete fundamental forms of an immersible triangulation of a 2D surface.
% See "Linear Surface Reconstruction from Discrete Fundamental Forms on 
% Triangle Meshes" by Wang, Liu, and Tong (2012).
%NOTE: ASSUMES THAT FACES ARE CCW ORIENTED
%
%   INPUT PARAMETERS:
%
%       - F:        #Fx3 face connectivity list
%       - L:        #Ex1 edge length list
%       - Phi:      #Ex1 list of bending angles, i.e. angles between
%                   adjacent face normal vectors
%
%   OUTPUT PARAMETERS:
%
%       - compErr:  #Vx1 list of compatibility conditions.
%                   Boundary vertex errors are set to NaN
%
%       - allGMPC:  #Vx1 cell array of matrix compatibility conditions
%                   Boundary vertex matrices are set to nan(3)
%
%       - RE:       #Ex2 cell array of frame field rotation matrices
%
%   by Dillon Cislo 01/15/2021

%------------------------------------------------------------------------
% INPUT PROCESSING
%------------------------------------------------------------------------
if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply edge length list'); end
if (nargin < 3), error('Please supply bending angle list'); end

validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive', 'real'});

% A set of 'phantom' vertices
V = zeros(max(F(:)), 3);

% A MATLAB-style representation of the input triangulation
TR = triangulation(F, V);
E = TR.edges;

validateattributes(L, {'numeric'}, ...
    {'vector', 'positive', 'real', 'numel', size(E,1)});
validateattributes(Phi, {'numeric'}, ...
    {'vector', 'real', 'numel', size(E,1)});

if (size(L,2) ~= 1), L = L.'; end
if (size(Phi, 2) ~= 1), Phi = Phi.'; end

% Mesh Topology Processing ----------------------------------------------

% Edge-face correspondence tool
try

    % Attempt to use NES edge ordering
    efIDx = edgeFaceCorrespondence(F);

catch

    warning('Unable to use NES edge ordering');
    efIDx = TR.edgeAttachments(E);
    efIDx = cell2mat(cellfun(@(x) repmat(x, 1, mod(numel(x),2)+1), ...
        efIDx, 'UniformOutput', false));
    % efIDx = sort(efIDx, 2); % Arbitrary choice

end

% Face-edge correspondence tool
e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, sort(E, 2), 'rows' );
[~, e2IDx] = ismember( e2IDx, sort(E, 2), 'rows' );
[~, e3IDx] = ismember( e3IDx, sort(E, 2), 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

% Edge-in-face correspondence tool
EinFIDx = full(sparse( feIDx, ...
    2-( efIDx(feIDx,1) == repmat((1:size(F,1)).', 3, 1) ), ...
    repmat([1 2 3], size(F,1), 1), size(E,1), 2 ));
EinFIDx(EinFIDx(:,2) == 0, 2) = EinFIDx(EinFIDx(:,2) == 0, 1);

% A CCW sorted list of faces attached to each vertex
vA = TR.vertexAttachments;

% Edges on the boundary (logical vector)
bdyEIDx = efIDx(:,1)-efIDx(:,2) == 0;

% Vertex IDs on the boundary/bulk (logical vector)
bdyIDx = TR.freeBoundary;
bdyIDx = bdyIDx(:,1);
bulkIDx = ~ismember((1:size(V,1)).', bdyIDx);
bdyIDx = ~bulkIDx;

% A CCW sorted list of edges attached to each vertex and an indicator for
% whether the face ordering follows the ordering in the edge-face
% correspondence tool

% OLD NON-VECTORIZED METHOD
% vAE = cell(size(vA,1), 1);
% for v = 1:size(vAE,1)
% 
%     if bdyIDx(v)
% 
%         vAE{v} = nan(1);
%         continue;
% 
%     end
% 
%     % The number of faces in the current vertex 1-ring
%     numF1R = numel(vA{v});
% 
%     % The faces in the vertex 1-ring and the next face in the CCW ordering
%     F1R = [vA{v}.', circshift(vA{v}.', -1, 1)];
% 
%     E1R = zeros(numF1R, 2);
%     for j = 1:numF1R
% 
%         [goodEdge, E1R(j,1)] = ismember(F1R(j,:), efIDx, 'rows');
%         if goodEdge
%             E1R(j,2) = 1;
%         else
%             [goodEdge, E1R(j,1)] = ismember(F1R(j,:), efIDx(:, [2 1]), 'rows');
%             E1R(j,2) = 2;
%             assert(goodEdge, 'Vertex-edge attachment list invalid');
%         end
% 
%     end
% 
%     vAE{v} = E1R;
% 
% end

vAE = cellfun(@(x) [x.' circshift(x.', -1, 1)], vA, 'Uni', false);
vAE = vertcat(vAE{:});

[~, vAE1] = ismember( vAE, efIDx, 'rows' );
[isE2, vAE2] = ismember( vAE, fliplr(efIDx), 'rows' );

vAE = vAE1;
vAE(isE2) = vAE2(isE2);
vAE = [ vAE, (1+isE2)];
vAE = mat2cell(vAE, cellfun(@numel, vA));

% Mesh Geometry Processing ----------------------------------------------

% Re-cast edge lengths onto faces
L_F = L(feIDx);
L2_F = L_F.^2;

% Some convenience variables to vectorize the cosine law calculation
Gi2 = L2_F; Gj2 = circshift(L2_F, [0 -1]); Gk2 = circshift(L2_F, [0 -2]);
Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]); 

% The internal angles of each face
intAng = ( Gj2 + Gk2 - Gi2 ) ./ ( 2 .* Gj .* Gk );
intAng = acos(intAng);

%--------------------------------------------------------------------------
% CALCULATE LOCAL COMPATIBILITY CONDITIONS
%--------------------------------------------------------------------------
% This calculation is performed by considering the parallel transport of
% an orthonormal frame (defined on mesh faces) around a CCW loop
% surrounding each bulk vertex. The fundamental forms are compatible of the
% frame is unchanged after transportation around the loop. For simplicity,
% we always set the first axis of the frame in each face to align with the
% first edge in that same face

% Pre-Compute Rotation Matrices on Each Mesh Edge -------------------------

% The rotation angle to align the first axis in the 2D orthonormal frame of
% each face with one of its edges
theta_FE = [ 0 .* intAng(:,1), pi - intAng(:,3), pi + intAng(:,2) ];

RotZ = @(x) [cos(x), -sin(x), 0; sin(x), cos(x), 0; 0 0 1 ];
RotX = @(x) [1 0 0; 0, cos(x), -sin(x); 0, sin(x), cos(x)];
RotZPi = [-1 0 0; 0 -1 0; 0 0 1];

RE = cell(size(E,1), 2);
for i = 1:size(E,1)

    % Ignore boundary edges
    if bdyEIDx(i)
        
        RE{i,1} = nan(3);
        RE{i,2} = nan(3);
        continue;

    end

    RE{i,1} = ...
        RotZ( theta_FE(efIDx(i,1), EinFIDx(i,1)) ) * ...
        RotX( Phi(i) ) * RotZPi * ...
        RotZ( -theta_FE(efIDx(i,2), EinFIDx(i,2)) );
    
    RE{i,2} = RE{i,1}.';

end

% Generate Compatibility Matrix for Each Mesh Vertex --------------------

allGMPC = cell(size(V,1), 1);
for v = 1:size(V,1)

    if bdyIDx(v)

        allGMPC{v} = nan(3);
        continue;

    end

    % The edge in the vertex 1-ring
    E1R = vAE{v};

    R = eye(3);
    for i = 1:size(E1R, 1)

        R = R * RE{E1R(i,1), E1R(i,2)};

    end

    allGMPC{v} = R;

end

% Generate Scalar Error Value -------------------------------------------

compErr = cellfun(@(x) norm(x-eye(3), 'fro'), allGMPC);

end