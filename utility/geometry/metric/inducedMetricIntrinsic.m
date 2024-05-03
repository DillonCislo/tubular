function g = inducedMetricIntrinsic(F, V2D, L3D)
%INDUCEDMETRICINTRINSIC Calculates the discrete induced metric of a surface
%represented by a mesh triangulation using only the edge lengths of the 3D
%surface (i.e. no explicit 3D parameterization is required)
%
%   INPUT PARAMETERS
%
%       - F:        #Fx3 face connectivity list
%       - V2D:      #Vx2 2D vertex coordinate list
%       - L3D:      #Ex1 3D edge length list
%
%   OUTPUT PARAMETERS
%
%       - g:    #Fx1 cell array. Each entry is a 2x2 matrix representing
%               the induced metric on a face
%
%   by Dillon Cislo 05/12/2020

% Validate input triangulation
validateattributes( V2D, {'numeric'}, ...
    {'2d', 'ncols', 2, 'finite', 'real'} );
validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'finite', ...
    'real', 'positive', '<=', size(V2D,1)} );

TR2D = triangulation(F, V2D);
E = TR2D.edges;

validateattributes(L3D, {'numeric'}, ...
    {'vector', 'numel', size(E,1), 'finite', 'real'} );

% Construct topological structure tools
e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, sort(E, 2), 'rows' );
[~, e2IDx] = ismember( e2IDx, sort(E, 2), 'rows' );
[~, e3IDx] = ismember( e3IDx, sort(E, 2), 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

% The 3D edge lengths
L3D_F = L3D(feIDx);
L1 = L3D_F(:,1);
L2 = L3D_F(:,2);
L3 = L3D_F(:,3);

% The 2D directed edge vectors
e12D = V2D(F(:,3), :) - V2D(F(:,2), :);
e22D = V2D(F(:,1), :) - V2D(F(:,3), :);
e32D = V2D(F(:,2), :) - V2D(F(:,1), :);

% Augmented edge unit vectors (used for cross products)
e12D = [ e12D, zeros(size(e12D,1), 1) ];
e22D = [ e22D, zeros(size(e22D,1), 1) ];
e32D = [ e32D, zeros(size(e32D,1), 1) ];

% The face areas/face unit normals of the pullback mesh
% ( +Z if the face is CCW oriented, -Z if the face is CW oriented )
fN = cross( e12D, e22D, 2 );
fA = sqrt( sum( fN.^2, 2 ) ) ./ 2;

fN = fN ./ (2 .* fA);

% The outward pointing in-plane edge normals
t1 = cross(e12D, fN, 2 ); t1 = t1(:, [1 2]);
t2 = cross(e22D, fN, 2 ); t2 = t2(:, [1 2]);
t3 = cross(e32D, fN, 2 ); t3 = t3(:, [1 2]);

% Construct the induced metric
g = cell( size(F,1), 1 );
for f = 1:size(F,1)
    
    g{f} = -( ...
        ( L1(f).^2-L2(f).^2-L3(f).^2 ) .* kron( t1(f,:), t1(f,:)' ) + ...
        ( L2(f).^2-L3(f).^2-L1(f).^2 ) .* kron( t2(f,:), t2(f,:)' ) + ...
        ( L3(f).^2-L1(f).^2-L2(f).^2 ) .* kron( t3(f,:), t3(f,:)' ) ...
        ) ./ ( 8 .* fA(f).^2 );
    
end

end
