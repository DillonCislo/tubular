function [F, V] = sliceMeshAlongPath(F0, V0, P0)
%SLICEMESHALONGPATH Slices a mesh triangulation along a path specified by
%by an ordered list of points. The output is a manifold mesh triangulation
%which is forced to include all of the edges defined by the input path
%
%   INPUT PARAMETERS:
%
%       - F0:       #F0x3 input face connectivity list
%
%       - V0:       #V0xD input vertex coordinate list
%
%       - P0:       #P0xD input path point coordinate list
%
%   OUTPUT PARAMETERS:
%
%       - F:        #Fx3 output face connectivity list
%
%       - V:        #VxD output vertex coordinate list
%
%   by Dillon Cislo 2024/02/05

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------
validateattributes(V0, {'numeric'}, {'finite', 'real', '2d'});
inputDim = size(V0,2);
assert(ismember(inputDim, [2 3]), 'Input vertices must be 2D or 3D');
numV0 = size(V0, 1);

validateattributes(F0, {'numeric'}, {'integer', 'positive', 'finite', ...
    'real', '2d', 'ncols', 3, '<=', numV0});

validateattributes(P0, {'numeric'}, ...
    {'finite', 'real', '2d', 'ncols', inputDim});

if (inputDim == 2)
    V0 = [V0 zeros(numV0, 1)];
    P0 = [P0 zeros(size(P0,1), 1)];
end

assert( ~any(is_vertex_nonmanifold(F0)), ...
    'Input mesh contains nonmanifold vertices');
assert( isempty(nonmanifold_edges(F0)), ...
    'Input mesh contains nonmanifold edges');

E0 = edges(triangulation(F0, V0));

e1IDx = sort( [ F0(:,3), F0(:,2) ], 2 );
e2IDx = sort( [ F0(:,1), F0(:,3) ], 2 );
e3IDx = sort( [ F0(:,2), F0(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, E0, 'rows' );
[~, e2IDx] = ismember( e2IDx, E0, 'rows' );
[~, e3IDx] = ismember( e3IDx, E0, 'rows' );

feIDx0 = [ e1IDx e2IDx e3IDx ];

% Snap the path point to the input mesh
[~, fID, P] = point_mesh_squared_distance(P0, V0, F0);
BC = barycentric_coordinates(P, ...
    V0(F0(fID,1),:), V0(F0(fID,2),:), V0(F0(fID,3),:));
% newEdges = [ (1:(size(P,1)-1)).', (2:size(P,1)).'];

% The distance threshold used to detect duplicated points
minDist = 1e-7;
BC(abs(BC) < 1e-10) = 0;

% Remove duplicate points in the input path
duplPoints = find_duplicate_rows(P, minDist);
if ~isempty(duplPoints)
    duplPoints = {duplPoints.idx}.';
    rmIDx = cell2mat(cellfun(@(x) ...
        [ x(1) * ones(numel(x)-1, 1), x(2:end) ], ...
        duplPoints, 'Uni', false));
    P(rmIDx(:,2), :) = [];
    BC(rmIDx(:,2), :) = [];
    fID(rmIDx, :) = [];
    % newEdges = changem(newEdges, rmIDx(:,1), rmIDx(:,2));
end

% Update the new prospective edges assuming that the points will be
% appended to the end of the vertex coordinate list
% newEdges = newEdges + numV0;

% Remove points that are already present in the vertex coordinate list
% [nnIDx, nnDists] = knnsearch(V0, P);
[~, nnDists] = knnsearch(V0, P);
rmIDx = (nnDists < minDist);
if any(rmIDx)
    % oldIDx = (1:size(P,1)).';
    % oldIDx = oldIDx(rmIDx);
    % newIDx = nnIDx(rmIDx);
    P(rmIDx, :) = [];
    BC(rmIDx, :) = [];
    fID(rmIDx, :) = [];
    % newEdges = changem(newEdges, newIDx, oldIDx);
end

assert(isempty(find_duplicate_rows([V0; P], minDist)), ...
    'Processed input still contains duplicate rows');

%--------------------------------------------------------------------------
% SLICE MESH ALONG PATH
%--------------------------------------------------------------------------

% Split Edges To Include New Vertices -------------------------------------

onEdge = sum(BC == 0, 2) == 1;
splitFIDx = fID(onEdge);
[I, splitEFIDx] = find(BC(onEdge, :) == 0);
splitEFIDx(I) = splitEFIDx;
splitEIDx = sub2ind(size(F0), splitFIDx, splitEFIDx);
splitEIDx = feIDx0(splitEIDx);
splitEMP = (V0(E0(splitEIDx, 1), :) + V0(E0(splitEIDx, 2), :))/2;

splitEdges = zeros(numel(splitEFIDx), 2);
splitEdges(splitEFIDx == 1, :) = F0(splitFIDx(splitEFIDx == 1), [2 3]);
splitEdges(splitEFIDx == 2, :) = F0(splitFIDx(splitEFIDx == 2), [3 1]);
splitEdges(splitEFIDx == 3, :) = F0(splitFIDx(splitEFIDx == 3), [1 2]);

[V, F] = split_edges(V0, F0, splitEdges);
assert(size(V,1) == (numV0+size(P,1)), 'Extra points added');
pIDx = knnsearch(splitEMP, V((numV0+1):end, :));
V((numV0+1):end, :) = P(pIDx, :);

assert( ~any(is_vertex_nonmanifold(F)), ...
    'Input mesh contains nonmanifold vertices');
assert( isempty(nonmanifold_edges(F)), ...
    'Input mesh contains nonmanifold edges');

% Flip Edges Until All Required Edges Are Contained in the Mesh -----------

[pathIDx, pathDists] = knnsearch(V, P0);
assert( all(pathDists < minDist), ...
    'Points were not addded to the mesh properly');
pathEdges = [pathIDx(1:(end-1)), pathIDx(2:end)];

count = 0;
maxIter = 100;
while true
    
    TR = triangulation(F,V);
    E = TR.edges;
    
    edgeInMesh = ismember(sort(pathEdges, 2), E, 'rows');
    if all(edgeInMesh)
        if (inputDim == 2), V = V(:, [1 2]); end
        return;
    end
    
    % YOU NEED TO ACTUALLY IMPLEMENT THIS!!!
    
    count = count+1;
    if (count > maxIter)
        error('Edge flips failed after %d iterations', maxIter);
    end
    
end

end


function [ dupl, C ] = find_duplicate_rows(A, distThresh)

[C, ia, ~] = unique( A, 'rows' );

if (size(A,1) == size(C,1))
    % disp('There are no duplicate rows!');
    dupl = {};
    return;
end

rep_idx = setdiff(1:size(A,1), ia);
rep_val = unique( A(rep_idx,:), 'rows');

dupl_val = cell( size(rep_val,1), 1 );
dupl_idx = cell( size(rep_val,1), 1 );

for i = 1:size(rep_val,1)
    
    dupl_val{i} = rep_val(i,:);
    
    diffA = A - repmat( rep_val(i,:), size(A,1), 1 );
    dupl_idx{i} = find( sqrt(sum(diffA .* conj(diffA), 2)) < distThresh );
    
end

dupl = struct( 'val', dupl_val, 'idx', dupl_idx );

end

