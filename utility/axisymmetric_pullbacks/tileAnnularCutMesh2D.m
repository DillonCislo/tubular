function [ TF, TV2D, TQ, TVQ, tileNum ] = tileAnnularCutMesh2D( cutMesh, tileCount, options )
%TILEANNULARCUTMESH2D This function vertically tiles the orbifold pullback 
% of an annular cutMesh and returns the parameters of a single 
% triangulation that covers the pullback multiple times.
%
% Note: as of 2021, this required that the pathPairs(:, 2) are the last
% indices of the vertices in cutMesh.u. Now in 2022 this is no longer the
% case.
%
% Note: if the tiles appear on the "wrong side" of the original tile, so
% that the order of veritces is nonmonotonic with their y coordinate, set
% cutMesh.pathPairs = cutMesh.pathPairs(:, [2,1]) and try again.
%
%   INPUT PARAMETERS:
%       - cutMesh:          A struct defining the cut 3D annulus with
%                           fields (See 'cylinderCutMesh.m'):
%                               f : #faces x 3 int array
%                               u : #vertices x 2 float array 
%                           and optional fields
%                               vn : #vertices x D numeric array
%                               quality : #faces x QDim numeric array
%
%       - tileCount:        The vertical tiling parameters.
%                           tileCount(1) tiles above the basic tile.
%                           tileCount(2) tiles below the basic tile.
%  
%       - options:          struct with fields 
%                           enforceUniformShift: (bool) check that all path
%                               pairs are equidistant in pullback space (u)
%                           preview: (bool) inspect progress as we go
%                           forceShift: (numeric) override shift from mesh
%                           path pairs and use this value instead
%
%   OUTPUT PARAMETERS:
%       - TF:               #Fx3 face connectivity list of the combined
%                           triangulation.
%       - TV2D:             #Vx2 2D pullback coordinate list of the
%                           combined triangulation.
%       - TQ:               #FxQDim face quality array for tiled mesh
%       - TVQ:              #VxQDim vertex quality array for tiled mesh
%       - tileNum:          #Vx1 int, tile index of each vertex
%
%  SEE ALSO
%   tileAnnularCutMesh.m --> for meshes with 3D push-forward
%  
% by Dillon Cislo / NPMitchell 2019-2022
% 2022 updates: allows more flexible indexing of pathPairs
%
%==========================================================================
% THE GEOMETRY OF THE CUT MESH:
%
% The field 'cutMesh.pathPairs' is a (#CV)x2 array that holds the
% correspondences between the original vertex IDs of the cut vertices and
% the vertex IDs of their duplicates.  
%
% Its columns are given by: [ (S1->E1), (S2->E2) ]
%
% Where (S1,E1) are the IDs of the original start and end points of the cut
% path and (S2,E2) are the IDs of the duplicated start and end points. The
% final output of the embedding procedure will have the following geometry
%
%           (4)
%    (S1)--------(E1)
%     |            |
%     |            |
% (1) |            | (3)
%     |            |
%     |            |
%    (S2)--------(E2)
%           (2)
%
% Where segments (2)+(4) are identified to create the topological annulus.
% The template for the output region in the (u,v) plane is the domain:
% [0 1] X [0 1]
%
%==========================================================================
%
% Example usage                           
% -------------     
%                              19--20-21 <---- topSeamIDx for positive shifts
%                               | \| \|
%                              16--17-18
%                               | \| \|
%           13--14-15 <----    13--14-15 <---- topSeamIDx
%            | \| \|            | \| \|
%           10--11-12          10--11-12
%            | \| \|            | \| \|
% 7--8--9    7--8--9  <----     7--8--9  <---- topSeamIDx
% | \| \|    | \| \|            | \| \|
% 4--5--6    4--5--6            4--5--6
% | \| \|    | \| \|            | \| \|
% 1--2--3    1--2--3  <----     1--2--3  <---- 

% 7--8--9    7--8--9            7--8--9         for negative shifts
% | \| \|    | \| \|            | \| \|
% 4--5--6    4--5--6            4--5--6
% | \| \|    | \| \|            | \| \|
% 1--2--3    1--2--3  <----     1--2--3  <---- topSeamIDx
%            | \| \|            | \| \|
%           13--14-15          13--14-15 
%            | \| \|            | \| \|
%           10--11-12 <----    10--11-12 <---- topSeamIDx
%                               | \| \|
%                              19--20-21 
%                               | \| \|
%                              16--17-18 <---- topSeamIDx
% 
% cutMesh = struct() ;                       
% x = [ 0,1,2, 0,1,2, 0,1,2];               
% y = [0, 0, 0, 1, 1, 1, 2, 2, 2] ;      
% cutMesh.u = [x(:), y(:)] ;
% cutMesh.pathPairs = [7,1; 8,2; 9,3] ;
% cutMesh.pathPairs = [1,7; 2,8; 3,9] ;
% cutMesh.f = [1,2,4; 2,5,4; 2,3,5; 3,6,5; ...
%              4,5,7; 5,8,7; 5,6,8; 6,9,8] ;
% tileAnnularCutMesh2D(cutMesh, [0, 1])
%
%                              19--20-21 <----
%                               | \| \|
%                              16--17-18
%                               | \| \|
%           14--10-15 <----    13--14-15 <----
%            | \| \|            | \| \|
%           11--12-13          10--11-12
%            | \| \|            | \| \|
% 7--2--9    7--2--9  <----     7--8--9  <----
% | \| \|    | \| \|            | \| \|
% 4--5--6    4--5--6            4--5--6
% | \| \|    | \| \|            | \| \|
% 1--8--3    1--8--3  <----     1--2--3  <---- topSeamIDx
%
%
% 7--2--9    7--2--9            7--2--9         for negative shifts
% | \| \|    | \| \|            | \| \|
% 4--5--6    4--5--6            4--5--6
% | \| \|    | \| \|            | \| \|
% 1--8--3    1--8--3  <----     1--8--3  <---- topSeamIDx
%            | \| \|            | \| \|
%           12--13-14          12--13-14 
%            | \| \|            | \| \|
%           10--15-11 <----    10--15-11 <---- topSeamIDx
%                               | \| \|
%                              18--19-20 
%                               | \| \|
%                              16--21-17 <---- topSeamIDx
% 
% cutMesh = struct() ;                       
% x = [ 0,1,2, 0,1,2, 0,1,2];               
% y = [ 0,2,0, 1,1,1, 2,0,2] ;      
% cutMesh.u = [x(:), y(:)] ;
% cutMesh.pathPairs = [7,1; 8,2; 9,3] ;
% cutMesh.pathPairs = [1,7; 8,2; 3,9] ;
% cutMesh.f = [1,8,4; 8,5,4; 8,3,5; 3,6,5; ...
%              4,5,7; 5,2,7; 5,6,2; 6,9,2] ;
% options = struct('enforceUniformShift', true, 'preview', true) ;
% tileAnnularCutMesh2D(cutMesh, [0, 1], options)

% Default tiling creates three stacked tiles
if nargin < 2
    tileCount = [1 1];
end

% Input options (optional)
enforceUniformShift = false ;
forceShift = NaN ;
preview = false ;
if nargin > 2
    if isfield(options, 'enforceUniformShift')
        enforceUniformShift = options.enforceUniformShift ;
    end
    if isfield(options, 'preview')
        preview = options.preview ;
    end
    if isfield(options, 'shift')
        forceShift = options.shift ;
    elseif isfield(options, 'forceShift')
        forceShift = options.forceShift ;
    end
end

% Output handling prep
if nargout > 2
    compute_face_quality = isfield(cutMesh, 'quality') ;
else
    compute_face_quality = false ;
end
if nargout > 3
    compute_vertex_quality = isfield(cutMesh, 'vertex_quality') ;
else
    compute_vertex_quality = false ;
end



% Verify input cut mesh
if ~isfield( cutMesh, 'u' )
    error("Cut mesh must contain 2D pullback coordinates as field 'u'");
end

pathPairs = cutMesh.pathPairs;

% The vertex IDs defining the current top seam of the combined tile
topSeamIDx = pathPairs(:,1);

% Find the location of the basic bottom seam vertices in the basic face
% connectivity list as a #(BF)x3 logical array
bottomSeamLoc = ismember( cutMesh.f, pathPairs(:,2) );

% Combined triangulation is just the basic tile to start
TF = cutMesh.f;
TV2D = cutMesh.u;
if compute_face_quality 
    TQ = cutMesh.quality ;
    TQ0 = TQ ;
else
    TQ = [] ;
end
if compute_vertex_quality
    TVQ = cutMesh.vertex_quality ;
    TVQ0 = TVQ ;
end

% Find the vertical shift between tiles (should just be 1 or 2*pi)
if isnan(forceShift)
    shift = mean(cutMesh.u( pathPairs(1,1), 2 ) - cutMesh.u( pathPairs(1,2), 2 ) );
else
    shift = forceShift ;
end

if enforceUniformShift
    shifts = cutMesh.u( pathPairs(:,1), 2 ) - cutMesh.u( pathPairs(:,2), 2 );
    assert(all(abs(shifts - shift) < 1e-7))
end
disp(['Tiling 2D mesh with spacing dY = ' num2str(shift)])

% Due to the structure of the cutMesh generation process it is easiest to
% add all new tiles to the top of the basic tile and then shift to reflect
% the desired tiles below the basic tile
tileNum = ones(size(cutMesh.u,1),1) ;
for i = 1:sum(abs(tileCount))
    
    % The parameters of the basic tile
    face = cutMesh.f(:);
    V2D = cutMesh.u;
    
    % Sew a shifted basic tile to the top of the current combined
    % triangulation -------------------------------------------------------
    
    % The shifted pullback vertices
    V2D(:,2) = V2D(:,2) + i * shift;
    if compute_vertex_quality
        TVQi = TVQ0 ;
        TVQi( pathPairs(:, 2), : ) = [] ;
    end
    
    % Remove the bottom seam from the vertex lists
    V2D( pathPairs(:,2), : ) = [];
    unreferenced = pathPairs(:,2) + size(TV2D, 1) ;
        
    % Part 1: Replace basic bottom seam IDs with the current top seam IDs
    % This replaces the BottomSeamLocations pathPairs(:,2) with top indices 
    for j = 1:length(topSeamIDx)
        face( ismember(face, pathPairs(j,2) ) ) = topSeamIDx(j);
        % Note these updates are protected in part 2 since
        % the replaced elements are in bottomSeamLocations.
        % Note these updates are protected in part 3 since unreferenced
        % vertices are necessarily higher index than any element of 
        % pathPairs(:).
    end
    
    % Part 2: Update the basic tile face list to reflect the fact that 
    % vertices will be added at the end of the combined vertex coordinate 
    % list
    face( ~bottomSeamLoc(:) ) = face( ~bottomSeamLoc(:) ) + size(TV2D,1);
           
    % Part 3: Update face list
    % check that the unreferenced vertices are not actually in the face
    % list
    assert(~any(ismember(unreferenced, face(:))))
    face0 = face ;
    for pp = 1:length(unreferenced)
        pt2rm = unreferenced(pp) ;
        face( face0 > pt2rm ) = face( face0 > pt2rm ) - 1 ;
    end
    
    % Reshape the face connectivity list
    face = reshape( face, size( cutMesh.f ) );
    % maxval = size(TV2D, 1) + size(V2D, 1) ;
    assert(max(face(:)) == size(TV2D, 1) + size(V2D, 1))
    
    % Update the current top seam vertex IDs
    % Offset by #(eliminated indices < kept indices pathPairs(:, 1))
    for pp = 1:length(pathPairs(:, 1))
        offsets(pp) = sum(pathPairs(:, 2) < pathPairs(pp, 1));    
    end
    topSeamIDx = pathPairs(:,1) + size(TV2D,1) - offsets(:) ;
    if preview
        disp(['topSeamIDx -> [' num2str(topSeamIDx') ']'])
    end
    
    % Update combined lists
    TF = [ TF; face ];
    TV2D = [ TV2D; V2D ];
    
    % % Update vertices to remove unreferenced (tiled vertices)
    % if i == 1
    %     unreferenced2 = setdiff(1:size(TV3D, 1), TF(:)) ;
    %     assert(all(sort(unreferenced(:)) == sort(unreferenced2(:)))) ;
    % end
    % 
    % [ TF, TV3D, oldVertexIDx, newVertexIDx] = ...
    %     remove_vertex_from_mesh( TF, TV3D, unreferenced) ;
    % TV2D = TV2D(oldVertexIDx, :) ;
    
    % Face quality is simply concatenated since physical face list is
    % unaltered, simply translated and reindexed along seams
    if compute_face_quality
        TQ = [TQ; TQ0 ] ;
    end
    if compute_vertex_quality 
        TVQ = [TVQ; TVQi] ;
    end
    tileNum = [tileNum ; ones(size(V2D,1),1) + i ] ;
    
    % Final checks for this addition
    assert(length(unique(TF(:))) == size(TV2D, 1))
    % Inspect output
    if preview
        clf
        trisurf(triangulation(TF, cat(2, TV2D, (1:size(TV2D,1))')), 1:size(TV2D,1)) ;
        hold on;
        plot3(TV2D(:, 1), TV2D(:, 2), (1:size(TV2D, 1))', '.')
    end
    
end

% Shift the coordinates of the combined triangulation to reflect the
% desired number of tiles below the basic tile
TV2D(:,2) = TV2D(:,2) - abs(tileCount(2)) * shift;

% Final checks
assert(length(unique(TF(:))) == size(TV2D, 1))
% Inspect output
if preview
    clf
    trisurf(triangulation(TF, cat(2, TV2D, (1:size(TV2D,1))')), 1:size(TV2D,1)) ;
    hold on;
    plot3(TV2D(:, 1), TV2D(:, 2), (1:size(TV2D, 1))', '.')
    % Look for unreferenced vertices
    missing = setdiff(1:max(TF(:)), TF(:)) ;
    assert(isempty(missing))
end

end

