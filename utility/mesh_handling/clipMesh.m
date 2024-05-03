function [newFace, oldFaceIDx, newVertex, oldVertexIDx] = ...
    clipMesh(face, vertex, rmAllTypes)
% CLIPMESH Remove dangling triangles, ears, nonmanifold edges, nonmanifold
% vertices and unreferenced vertices from a mesh triangulation
% 
% INPUT PARAMETERS:
%
%   face:         #Fx3 input face connectivity list
%
%   vertex:       #VxD input vertex coordinate list
%
%   rmAllTypes:   If true, this algorithm will remove dangling triangles,
%                 ears, nonmanifold edges, nonmanifold vertices, 
%                 and unreferenced vertices.  If false, it will only
%                 remove unreferenced vertices.
%
% OUTPUT PARAMETERS:
%
%   newFace:      #F'x3 output face connectivity list
%
%   oldFaceIDx:   #F'x1 indices of new faces in old triangulation, i.e.
%                 newFace = face(oldFaceIDx, :);
%
%   newVertex:    #V'xD output vertex coordinate list
%
%   oldVertexIDx: #V'x1 indices of new vertices in old triangulation, i.e.
%                 newVertex = vertex(oldVertexIDx, :);
% 
% By Dillon Cislo 02/04/2020

%--------------------------------------------------------------------------
% Validate Inputs
%--------------------------------------------------------------------------

validateattributes( vertex, {'numeric'}, ...
    {'2d', 'finite', 'real', 'nonnan'} );
validateattributes( face, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real', 'integer', 'positive'} );

assert(all(ismember(face, (1:size(vertex,1)).'), 'all'), ...
    'Face connectivity list contains unreferenced vertices!');

if nargin < 3
    rmAllTypes = true;
end

%--------------------------------------------------------------------------
% Remove Undesired Vertices
%--------------------------------------------------------------------------

% Lone points are defined as vertices belonging to ears (boundary triangles
% sharing only one edge with the mesh) or vertices belonging to triangles
% sharing only one vertex with the edge of the mesh or vertices that are
% not referenced by the input triangulation
num_lone_points = Inf;

% Make temporary copies of the input parameters
f = face;
v = vertex;

while( num_lone_points ~= 0 )
    
    % Determine the number of faces attached to each vertex
    vertexFaceCount = full(sparse(f, 1, 1, size(v,1), 1));
    
    % The number and ID of vertices to remove
    if rmAllTypes
        
        % Remove all vertices in less than one face
        lone_points = vertexFaceCount < 2;
        
        % Remove non-manifold vertices
        lone_points = lone_points | is_vertex_nonmanifold(f);
        
        % Remove vertices comprising non-manifold edges
        e = edges(triangulation(f,v));
        lone_points = lone_points | ...
            ismember((1:size(v,1)).', unique(e(nonmanifold_edges(f), :)));
        
    else
        
        % Remove only unreferenced vertices
        lone_points = vertexFaceCount < 1;
        
    end
    
    num_lone_points = sum(lone_points);
    lone_points = find(lone_points);    
    
    % Remove the lone vertices at the specified index values
    newVertex = v;
    newVertex(lone_points, :) = [];
    
    % Find the new index for each of the new vertices   
    [~, newVertexIDx] = ismember( v, newVertex, 'rows');
    
    % Find any faces that contained the lone vertices and remove them
    newFace = f;
    lone_faces = any( ismember( newFace, lone_points ), 2 );
    newFace( lone_faces, : ) = [];
    
    % Now update the vertex indices in the connectivity list
    newFace = newVertexIDx(newFace);
    
    % Update the temporary copies
    f = newFace;
    v = newVertex;
    
end

% if isempty(newVertex)
%     warning('All vertices have been removed!');
% end
% 
% if isempty(newFace)
%     warning('All faces have been removed!');
% end

% Find the indices of the new vertices in the old triangulation
% NOTE: This assumes there are no duplicate vertices
[~, oldVertexIDx] = ismember( newVertex, vertex, 'rows' );

% Find the indices of the new faces in the old triangulation
% NOTE: this assumes there are no duplicate faces
oldFaceIDx = ismember(sort(newFace,2), sort(face,2), 'rows');

end


    
