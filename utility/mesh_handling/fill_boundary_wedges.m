function F = fill_boundary_wedges(FOld, V, varargin)
%FILL_BOUNDARY_WEDGES Iteratively fill missing 'wedges' in a mesh with
%boundary. At each iteration, for each boundary component, this method
%identifies next-nearest-neighbor boundary vertices (x-o-x) that can be
%joined by an edge. If the incident angle opposite this new edge is below a
%threshold and the edge length is below a threshold the edge is accepted.
%Included optional checking for mesh self-intersections. Note that vertices
%are not removed or added.
%
%   INPUT PARAMETERS:
%
%       - FOld:     #FOx3 input face connectivity list
%       - V:        #Vx3 vertex coordinate list
%
%   OPTIONAL INPUT PARAMETERS (Name, Value)-pairs:
%
%       - ('MaxAngle', maxAngle = 120): Maximum angle (in degrees!) above
%       which the candidate edge opposite this angle is rejected.
%
%       - ('MaxLength', maxLength = 2): Maximum length as a fraction of the
%       mean edge length in the original mesh above which a candidate edge
%       is rejected.
%
%       - ('MaxIterations', maxIter = 35): Maximum number of iterations to
%       run the fill process. Set to 'Inf' to run until all wedges are
%       filled.
%
%       - ('CheckQuality', checkQuality = false): Whether or not
%       to check for mesh quality. Output mesh should have no
%       self-intersections, no nonmanifold edges/vertices, and the normal
%       vectors of the new faces should remain properly oriented with the
%       original face normals. NOTE: This method assumes that the input
%       mesh is already of passing quality.
%
%       - ('Verbose', verbose = false): Whether or not to produce verbose
%       output.
%
%   OUTPUT PARAMETERS:
%
%       - F:    #Fx3 output face connectivity list
%
%   by Dillon Cislo 2023/12/15

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------
validateattributes(V, {'numeric'}, {'2d', 'real', 'finite', 'nonnan'});
if (size(V,2) ~= 3)
    if (size(V,2) == 2)
        V = [V, zeros(size(V,1), 1)];
    else
        error('Vertices must be 2D or 3D');
    end
end

VN = per_vertex_normals(V, FOld, 'Weighting', 'angle');

validateattributes(FOld, {'numeric'}, ...
    {'2d', 'ncols', 3, 'real', 'positive', 'integer', '<=', size(V,1)});
EOld = edges(triangulation(FOld, V));
LOld = V(EOld(:,2),:)-V(EOld(:,1),:);
LOld = sqrt(sum(LOld.^2, 2));
meanLength = mean(LOld);

% Set default parameters
maxAngle = 120;
maxLength = 2;
maxIter = 35;
checkQuality = false;
verbose = false;

for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double'), continue; end
    if isa(varargin{i}, 'logical'), continue; end
    
    if strcmpi(varargin{i}, 'MaxAngle')
        maxAngle = varargin{i+1};
        validateattributes(maxAngle, {'numeric'}, ...
            {'scalar', 'real', 'positive'});
    end
    
    if strcmpi(varargin{i}, 'MaxLength')
        maxLength = varargin{i+1};
        validateattributes(maxLength, {'numeric'}, ...
            {'scalar', 'real', 'positive'});
    end
    
    if strcmpi(varargin{i}, 'MaxIterations')
        maxIter = varargin{i+1};
        if isinf(maxIter)
            validateattributes(maxIter, {'numeric'}, ...
            {'scalar', 'real', 'positive'});
        else
            validateattributes(maxIter, {'numeric'}, ...
            {'scalar', 'real', 'positive', 'integer'});
        end
    end
    
    if strcmpi(varargin{i}, 'CheckQuality')
        checkQuality = varargin{i+1};
        validateattributes(checkQuality, {'logical'}, {'scalar'});
    end
    
    if strcmpi(varargin{i}, 'Verbose')
        verbose = varargin{i+1};
        validateattributes(verbose, {'logical'}, {'scalar'});
    end
    
end

% Convert length scale into absolute length
maxLength = maxLength * meanLength;

% Convert maximum angle to radians
maxAngle = deg2rad(wrapTo360(maxAngle));

% Check to see if self-intersection code is available
if (exist('mesh_self_intersection_3d', 'file') ~= 3)
    if checkQuality
        warning('Cannot find self-intersection code! Ignoring command...');
    end
    checkQuality = false;
end

%--------------------------------------------------------------------------
% Fill In Wedges
%--------------------------------------------------------------------------

F = FOld;
allBdys = DiscreteRicciFlow.compute_boundaries(F);
if isempty(allBdys), return; end

% Make all boundary lists column vectrors
allBdys = cellfun(@(x) x.', allBdys, 'Uni', false);

% Check orientation of each boundary element
allCOM = barycenter(V, F);
TR = triangulation(F, V);
allFN = TR.faceNormal;
isCCW = false(numel(allBdys), 1);
for bidx = 1:numel(allBdys)
    
    bdy = allBdys{bidx}; % The current boundary list
    nbdy = circshift(bdy, [-1 0]); % The next vertex
    emp = (V(bdy, :) + V(nbdy, :))./2; % Boundary edge midpoints
    
    % The IDs of the attached faces
    fID = edgeAttachments(TR, bdy, nbdy);
    assert(all(cellfun(@numel, fID, 'Uni', true) == 1), ...
        'Invalid boundary edges');
    fID = cell2mat(fID);
    
    FN = allFN(fID, :);
    COM = allCOM(fID, :);
    
    orientationSign = normalizerow(cross(V(nbdy,:)-emp, COM-emp, 2));
    orientationSign = sign(dot(orientationSign, FN, 2));
    
    %*******************************************************************
    % This test works produces flipped results for outer vs. inner
    % boundaries. Checked for straight-forward multiply connected
    % topologies, but I'm not sure this is 100% robust for all topologies
    %*******************************************************************
    if all(orientationSign > 0)
        if (bidx == 1)
            isCCW(bidx) = true;
        else
            isCCW(bidx) = false;
        end
    elseif all(orientationSign < 0)
        if (bidx == 1)
            isCCW(bidx) = false;
        else
            isCCW(bidx) = true;
        end
    else
        error('Invalid mesh. Non-consistently ordered faces?');
    end
    
end

for iter = 1:maxIter
    
    if verbose
        fprintf('Iteration %d\n', iter);
    end
    
    curE = sort(edges(triangulation(F, V)), 2);
    anyBoundaryAltered = false;
    
    for bidx = 1:numel(allBdys)
        
        bdy = allBdys{bidx}; % The current boundary list
        if isempty(bdy), continue; end
        
        pbdy = circshift(bdy, [1 0]); % The previous vertex
        nbdy = circshift(bdy, [-1 0]); % The next vertex

        % Calculate external angle of edge triplets
        ep = normalizerow(V(pbdy, :) - V(bdy, :));
        en = normalizerow(V(nbdy, :) - V(bdy, :));
        
        z = cross(ep, en, 2);
        zhat = normalizerow(z);
        
        isConvex = dot(zhat, VN(bdy,:), 2);
        if isCCW(bidx)
            if (bidx == 1)
                isConvex = (isConvex < 0);
            else
                isConvex = (isConvex > 0);
            end
        else
            if (bidx == 1)
                isConvex = (isConvex > 0);
            else
                isConvex = (isConvex < 0);
            end
        end
        
        % This is always going to yield the smallest angle between edges
        extAng = wrapTo2Pi(2 .* atan2(dot(z, zhat, 2), 1+dot(ep, en, 2)));
        extAng(isConvex) = 2 * pi - extAng(isConvex);
        
        % Calculate candidate boundary edge lengths
        newL = V(nbdy, :) - V(pbdy, :);
        newL = sqrt(sum(newL.^2, 2));
        
        % Determine which edges satisfy the thresholds
        edgeIDx = (1:numel(bdy)).';
        goodEdges = (newL <= maxLength) & (extAng <= maxAngle);
        goodEdges = goodEdges & ... % Ensure no duplicate edges
            ~ismember(sort([pbdy, nbdy], 2), curE, 'rows');
        edgeIDx(~goodEdges) = [];
        
        if isempty(edgeIDx), continue; end
        
        if checkQuality
            
            minID = -1;
            while (~isempty(edgeIDx) && (minID < 0))
                
                % The new edge is the one with the smallest exterior angle
                [~, testIDLoc] = min(extAng(edgeIDx));
                testID = edgeIDx(testIDLoc);
                
                if isCCW(bidx)
                    if (bidx == 1)
                        testF = [F; pbdy(testID), nbdy(testID), bdy(testID) ];
                    else
                        testF = [F; pbdy(testID), bdy(testID), nbdy(testID) ];
                    end
                else
                    if (bidx == 1)
                        testF = [F; pbdy(testID), bdy(testID), nbdy(testID) ];
                    else
                        testF = [F; pbdy(testID), nbdy(testID), bdy(testID) ];
                    end
                end
                testTR = triangulation(testF, V);
                
                neighborIDx = testTR.neighbors(size(testF,1));
                neighborIDx = unique(neighborIDx(~isnan(neighborIDx)));
                neighborIDx = neighborIDx(:);
                testFN = testTR.faceNormal(size(testF,1));
                nnFN = testTR.faceNormal(neighborIDx);
                
                wrongDirection = ...
                    dot(nnFN, repmat(testFN, size(nnFN,1), 1), 2);
                wrongDirection = any(wrongDirection <= 0);
                
                [intersects, ~] = mesh_self_intersection_3d(testF, V);
                
                nonmanifold = any(is_vertex_nonmanifold(testF)) || ...
                    ~isempty(nonmanifold_edges(testF));
                
                % Only accept a new edge if it does not cause mesh
                % self-intersections
                if (intersects || nonmanifold || wrongDirection)
                    edgeIDx(testIDLoc) = [];
                else
                    minID = testID;
                    break;  
                end

            end
            
            if (minID < 0), continue; end
            
        else
            
            % The new edge is the one with the smallest exterior angle
            [~, minID] = min(extAng(edgeIDx));
            minID = edgeIDx(minID);
            
        end
        
        if isCCW(bidx)
            if (bidx == 1)
                F = [F; pbdy(minID), nbdy(minID), bdy(minID) ];
            else
                F = [F; pbdy(minID), bdy(minID), nbdy(minID) ];
            end
        else
            if (bidx == 1)
                F = [F; pbdy(minID), bdy(minID), nbdy(minID) ];
            else
                F = [F; pbdy(minID), nbdy(minID), bdy(minID) ];
            end
        end
        
        bdy(minID) = [];
        allBdys{bidx} = bdy;
        anyBoundaryAltered = true;

    end
    
    if ~anyBoundaryAltered, break; end
    
end

assert(isequal(unique(F, 'rows'), sortrows(F)), ...
    'Output contains duplicate faces');

end

