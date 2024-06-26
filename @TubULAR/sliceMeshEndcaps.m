function sliceMeshEndcaps(tubi, opts, methodOpts)
% SLICEMESHENDCAPS()
% Create cut meshes with endcaps removed 
% Note that a later step involves cutMesh, which slices along AP.
% Noah Mitchell 2019
% This version relies on Gabriel Peyre's toolbox called
% toolbox_fast_marching/
%
% Parameters
% ----------
% tubi : TubULAR class instance
%   class instance for which we cut off endcaps of mesh
% opts : struct with optional fields
%   adist_thres
%   pdist_thres
%   adist_thres2
%   pdist_thres2
%   aCapMethod
%   pCapMethod
%   aOffset
%   pOffset
%   aOffset2
%   pOffset2
%   aOffsetRate
%   pOffsetRate
%   aOffset2Rate
%   pOffset2Rate
%   aDistRate
%   pDistRate
%   tref
%   custom_aidx
%   custom_pidx
% methodOpts : struct with optional fields
%   tref : int (timestamp, not index)
%       reference time stamp for dorsal definition, after and before which 
%       we point match to find the "dorsal" boundary vertex on each endcap
%   quickScan : bool
%       first iterate through all points with large dt between frames to
%       check that settings are good, then go back and compute each frame 
%       in order
%   overwrite
%   save_figs
%   preview
%   timePoints
%
% Prerequisites
% -------------
% alignMeshesAPDV() after identifying APD pints
% extractCenterlineSeries() 
% 
%
% First run extract_centerline.m or extractCenterlineSeries() before 
% running this code.
% Run this code only after identifying anterior (A), posterior (P), and 
% dorsal anterior (D) points via computeAPDpoints().

%% Parameters from tubi
timePoints = tubi.xp.fileMeta.timePoints ;

%% Parameters
tref = tubi.xp.fileMeta.timePoints(1) ;     % which timepoint to use as reference, when
                                          % dorsal points on the boundary
                                          % are defined. Subsequent and
                                          % previous timepoints have dorsal
                                          % points point-matched in a way
                                          % which traces back to the
                                          % definition at tref
adist_thres = 20 ;  % distance threshold for cutting off anterior in pix
pdist_thres = 15 ;  % distance threshold for cutting off posterior in pix
adist_thres2 = adist_thres ;  % Second threshold on anterior distance
pdist_thres2 = pdist_thres ;  % Second threshold on posterior distance
overwrite = false ;  % recompute sliced endcaps
save_figs = true ;  % save images of cntrline, etc, along the way
preview = false ;  % display intermediate results
posterior_phi_subtractCOM = true ;  % For posterior pole determination 
                                    % of dorsal point, subtract the center
                                    % of mass                              
aCapMethod = 'geodesic' ;  % [geodesic,ball,xvalue,yvalue,zvalue] method for determining vertices to remove on anterior, could be 'geodisic', 'ball', 'cone', 'yvalue', etc
pCapMethod = 'geodesic' ;  % [geodesic,ball,xvalue,yvalue,zvalue] method for determining vertices to remove on posterior
aOffset = [0, 0, 0] ;  % Offset in pixels in APDV coord sys (will be rotated to pixel space)
pOffset = [0, 0, 0] ; 
aOffset2 = [0, 0, 0] ;  % Offset in pixels in APDV coord sys (will be rotated to pixel space)
pOffset2 = [0, 0, 0] ; 
aOffsetRate = [0, 0, 0] ;  % Offset rate (per timepoint) in pixels in APDV coord sys (will be rotated to pixel space)
pOffsetRate = [0, 0, 0] ; 
aOffset2Rate = [0, 0, 0] ;  % Offset rate (per timepoint) in pixels in APDV coord sys (will be rotated to pixel space)
pOffset2Rate = [0, 0, 0] ; 
aDistRate = [0, Inf] ;  % [ramp,duration; ramp2,duration2; etc] 
                        % Offset rate (per timepoint) in slope or pixels in 
                        % APDV coord sys, second number is number of
                        % timepoints to have this ramp
pDistRate = [0, Inf] ; 
custom_aidx = [];   % User supplied vertex IDs that will be removed to make the anterior cap
custom_pidx = [];   % User supplied vertex IDs that will be removed to make the posterior cap
quickScan = true ;  % first iterate through all points with large dt between frames to
                    %       check that settings are good, then go back and compute each frame 
                    %       in order


if isfield(opts, 'adist_thres')
    adist_thres = opts.adist_thres ;
end
if isfield(opts, 'pdist_thres')
    pdist_thres = opts.pdist_thres ;
end
if isfield(opts, 'adist_thres2')
    adist_thres2 = opts.adist_thres2 ;
end
if isfield(opts, 'pdist_thres2')
    pdist_thres2 = opts.pdist_thres2 ;
end
if isfield(opts, 'aCapMethod')
    aCapMethod = opts.aCapMethod ;
end
if isfield(opts, 'pCapMethod')
    pCapMethod = opts.pCapMethod ;
end
if isfield(opts, 'aOffset')
    aOffset = opts.aOffset ;
end
if isfield(opts, 'pOffset')
    pOffset = opts.pOffset ;
end
if isfield(opts, 'aOffset2')
    aOffset2 = opts.aOffset2 ;
end
if isfield(opts, 'pOffset2')
    pOffset2 = opts.pOffset2 ;
end
if isfield(opts, 'aOffsetRate')
    aOffsetRate = opts.aOffsetRate ;
end
if isfield(opts, 'pOffsetRate')
    pOffsetRate = opts.pOffsetRate ;
end
if isfield(opts, 'aOffset2Rate')
    aOffset2Rate = opts.aOffset2Rate ;
end
if isfield(opts, 'pOffset2Rate')
    pOffset2Rate = opts.pOffset2Rate ;
end
if isfield(opts, 'aDistRate')
    aDistRate = opts.aDistRate ;
end
if isfield(opts, 'pDistRate')
    pDistRate = opts.pDistRate ;
end
if isfield(opts, 'tref')
    tref = opts.tref ;
end
if isfield(opts, 'custom_aidx')
    custom_aidx = opts.custom_aidx;
end
if isfield(opts, 'custom_pidx')
    custom_pidx = opts.custom_pidx;
end
if isfield(methodOpts, 'overwrite')
    overwrite = methodOpts.overwrite ;
elseif isfield(opts, 'overwrite')
    overwrite = opts.overwrite ;
end
if isfield(methodOpts, 'save_figs')
    save_figs = methodOpts.save_figs ;
end
if isfield(methodOpts, 'preview')
    preview = methodOpts.preview;
end
if isfield(methodOpts, 'timePoints')
    timePoints = methodOpts.timePoints;
end
if isfield(methodOpts, 'quickScan')
    quickScan = methodOpts.quickScan;
end

% Determine if valid custom input had been supplied
useCustomPts = iscell(custom_aidx) && iscell(custom_pidx) &&...
    (numel(custom_aidx) == numel(timePoints)) && ...
    (numel(custom_pidx) == numel(timePoints));
if useCustomPts
    disp('Slicing endcaps based on user input');
end

% Invert rotation for offset
aOffXYZ = tubi.APDV2dxyz(aOffset) / tubi.ssfactor ;
pOffXYZ = tubi.APDV2dxyz(pOffset) / tubi.ssfactor ;
aOffRateXYZ = tubi.APDV2dxyz(aOffsetRate) / tubi.ssfactor ;
pOffRateXYZ = tubi.APDV2dxyz(pOffsetRate) / tubi.ssfactor ;
aOff2XYZ = tubi.APDV2dxyz(aOffset2) / tubi.ssfactor ;
pOff2XYZ = tubi.APDV2dxyz(pOffset2) / tubi.ssfactor ;
aOffRate2XYZ = tubi.APDV2dxyz(aOffset2Rate) / tubi.ssfactor ;
pOffRate2XYZ = tubi.APDV2dxyz(pOffset2Rate) / tubi.ssfactor ;
apDir = tubi.APDV2dxyz([1, 0,0]) ;
apDir = apDir ./ sqrt(sum(apDir .^ 2 )) ;

% figure parameters
xwidth = 16 ; % cm
ywidth = 10 ; % cm

% subsampling factor for the h5s used to train for mesh/acom/pcom/dcom
ssfactor = tubi.ssfactor ; 
[~, ~, xyzlim_um] = tubi.getXYZLims() ;

% Name output directory
outdir = tubi.dir.cylinderMesh ;
figoutdir = fullfile(outdir, 'images') ;
if ~exist(figoutdir, 'dir')
    mkdir(figoutdir) ;
end

%% Load AP coms
[apt_sm, ppt_sm] = tubi.getAPpointsSm ;
trefIDx = tubi.xp.tIdx(tref) ;

%% Iterate through each mesh

% Set up timepoint sequence to quickly scan through all timepoints (every
% 50), then more finely (every 20), then more finely (every 10), then do
% every timepoint. Since the pointmatching is serial, each subsequent pass
% recomputes meshes for previously done timepoints.
% ONLY DO THIS IS WE ARE THEN GOING TO OVERWRITE THE QUICK SCAN.
if (overwrite && ~useCustomPts && quickScan)
    if length(timePoints) > 49
        todo = trefIDx:50:length(timePoints); % first preview how slices will look
    else
        todo = [] ;
    end
    if length(timePoints) > 19
        todo2 = trefIDx:20:length(timePoints) ;
    else
        todo2 = [];
    end
    if length(timePoints) > 9
        todo3 = trefIDx:10:length(timePoints) ;
    else
        todo3 = [] ;
    end
else
    todo = [] ; 
    todo2 = [];
    todo3 = []; 
end
% NOTE: begin with tref, advance to end, then return to tref and go
% backwards
if useCustomPts
    todo4 = 1:length(timePoints);
else
todo4 = [trefIDx:length(timePoints), fliplr(1:(trefIDx-1)) ] ;
end
todo = [todo, todo2, todo3, todo4] ;
for ii=todo
    tt = timePoints(ii) ;
    Dt = tt - tref ;
    disp(['tt = ' num2str(tt)])
    
    apt = apt_sm(ii, :) ;
    ppt = ppt_sm(ii, :) ;
    
    %% Name the output mesh filename
    name = sprintf(tubi.fileBase.mesh, tt) ;
    meshfn = sprintfm(tubi.fullFileBase.mesh, tt) ;
    outfn = sprintfm(tubi.fullFileBase.cylinderMesh, tt) ;
    keepfn = sprintfm(tubi.fullFileBase.cylinderKeep, tt) ; 
    boundaryfn = sprintfm(tubi.fullFileBase.apBoundary, tt) ;
    
    outapd_boundaryfn = tubi.fileName.apBoundaryDorsalPts ;  % not a base name
    
    %% Read the mesh
    disp(['Loading mesh ' name]) ;
    
    % Compute the endcaps if not already saved
    if overwrite || ~exist(outfn, 'file') || ~exist(keepfn, 'file')
        
        if exist(outfn, 'file')
            disp(['Overwriting cylinder mesh: ' outfn])
        else
            disp(['Cylinder mesh not on disk: ' outfn])
        end
        if exist(keepfn, 'file')
            disp(['Overwriting keep indices: ' keepfn])
        else
            disp(['Keep indices not on disk: ' keepfn])
        end
        disp(['Computing endcaps for ' name])
        
        % Load raw mesh
        mesh = read_ply_mod(meshfn);
        
        if eulerCharacteristic(mesh) ~= 2
            disp(['WARNING: input mesh is not topological sphere: Euler Characteristic = ' ...
                num2str(eulerCharacteristic(mesh))])
            disp('...attempting to continue...')
        end
        
        if useCustomPts
            
            rmIDx = [custom_aidx{ii}; custom_pidx{ii}];
            
            % Remove custom vertices from mesh
            [faces, vtx, ~, ~] = ...
                remove_vertex_from_mesh(mesh.f, mesh.v, rmIDx);
            
            % Keep only the largest remaning connected component
            [faces, vtx, ~, ~] = ...
                remove_isolated_mesh_components(faces, vtx);

            % Remove any unreferences vertices
            [faces, vtx, ~, ~] = ...
                remove_unreferenced_vertices_from_mesh(faces, vtx);
            
            % Determine which vertices remain in the sliced mesh
            keep = knnsearch(mesh.v, vtx);
            vn = mesh.vn(keep, :);
            
            % Save the data in units of pixels (same as original mesh)
            disp(['Saving cylinder mesh to ' outfn])
            plywrite_with_normals(outfn, faces, vtx, vn)
            
            % Save the indices to keep when cutting off endcaps
            save(keepfn, 'keep') ;
            
            % Subsample the mes to match acom, pcom
            vtx = vtx / ssfactor;
            
            aOff = aOffXYZ + aOffRateXYZ * Dt ;
            acomOff = apt + aOff ;
            
            pOff = pOffXYZ + pOffRateXYZ * (tt - timePoints(1)) ;
            pcomOff = ppt + pOff ;
            
        else
            
            mesh = read_ply_mod(meshfn);
            % subsample the mesh to match acom, pcom
            vtx = mesh.v / ssfactor ;

            % Remove unreferenced vertices
            [ mesh.f, vtx, ~, oldVertexIDx] = remove_unreferenced_vertices_from_mesh( mesh.f, vtx ) ;
            fv = struct('f', mesh.f, 'v', vtx, 'vn', mesh.vn(oldVertexIDx, :)) ;
            disp('loaded closed mesh.')
            
            %% Remove anterior endcap
            % Measure distance to the posterior
            % Strategy: remove within distance of acom
            aOff = aOffXYZ + aOffRateXYZ * Dt ;
            acomOff = apt + aOff ;
            acom2 = apt + aOff2XYZ + aOffRate2XYZ * Dt ;
            if numel(aDistRate) > 6
                error('Code for more than three ramp rates here')
            end
            
            % Adjust parameters based on "rates", which ramp endcap
            % parameters based on which timepoint is under consideration.
            if abs(aDistRate(1)) > 0 || numel(aDistRate) > 2
                if numel(aDistRate) > 2
                    % second column is duration of ramp rate
                    if Dt > aDistRate(1, 2)
                        % Apply entire first ramp
                        adist0 = adist_thres + aDistRate(1,2) * aDistRate(1,1) ;
                        
                        if (Dt - aDistRate(1,2)) < aDistRate(2, 2)
                            adist0 = adist0 + (Dt - aDistRate(1,2)) * aDistRate(2,1) ;
                            disp(['Ramping with second rate, adist = ' num2str(adist0)])
                        else
                            adist0 = adist0 + aDistRate(2,2) * aDistRate(2, 1) ;
                        end
                        
                        if Dt > aDistRate(2, 2)
                            try
                                adist0 = adist0 + (Dt - aDistRate(1,3)) * aDistRate(3,1) ;
                                disp(['Done ramping, adist = ' num2str(adist0)])
                            catch
                                disp(['Done ramping, adist = ' num2str(adist0)])
                            end
                        end
                    else
                        adist0 = adist_thres + Dt * aDistRate(1,1) ;
                        disp(['Ramping with first rate, adist = ' num2str(adist0)])
                    end
                else
                    disp('Ramping with first value')
                    adist0 = adist_thres + aDistRate(1) * Dt ;
                end
            else
                adist0 = adist_thres ;
            end
            
            % Now remove the endcap
            disp(['Removing first/anterior endcap using method: ' aCapMethod])
            if strcmpi(aCapMethod, 'geodesic')    
                % Options for the marching
                [~, nearestID_a] = min(vecnorm(vtx - acomOff, 2, 2)) ;
                startpts = nearestID_a ;
                fmoptions = struct() ;
                fmoptions.end_points = [];   

                %   output is the distance function to the set of starting points.
                [adist_g] = perform_fast_marching_mesh(vtx', fv.f, startpts, fmoptions) ;
                pts_to_remove = find(adist_g < adist0) ;
            elseif strcmpi(aCapMethod, 'ball')
                adist2 = sum((vtx - acomOff) .^ 2, 2);
                pts_to_remove = find(adist2 < adist0^2) ;
            elseif strcmpi(aCapMethod, 'cone')
                % reference: https://stackoverflow.com/questions/12826117/how-can-i-detect-if-a-point-is-inside-a-cone-or-not-in-3d-space
                apex = acomOff ;
                % project vertices onto direction of cone (here is ap)
                coneDir = -apDir ;
                cone_dist = sum((vtx - apex) .* (coneDir .* ones(size(vtx))), 2) ;
                orth_dist = vecnorm(vtx - apex - cone_dist * coneDir, 2, 2) ;
                cone_radius = cone_dist * adist0 ;
                insideCone = orth_dist < cone_radius ;
                pts_to_remove = find(insideCone);
            elseif strcmpi(aCapMethod, 'ballCone') || strcmpi(aCapMethod, 'coneBall')
                % reference: https://stackoverflow.com/questions/12826117/how-can-i-detect-if-a-point-is-inside-a-cone-or-not-in-3d-space
                apex = acomOff ;
                % project vertices onto direction of cone (here is ap)s
                coneDir = -apDir ;
                cone_dist = sum((vtx - apex) .* (coneDir .* ones(size(vtx))), 2) ;
                orth_dist = vecnorm(vtx - apex - cone_dist * coneDir, 2, 2) ;
                cone_radius = cone_dist * adist0 ;
                insideCone = orth_dist < cone_radius ;
                % Also threshold on absolute distance from COMOffset
                adist2 = sum((vtx - acom2) .^ 2, 2);
                pts_to_remove = find(insideCone & (adist2 < adist_thres2^2));
            elseif strcmpi(aCapMethod, 'xvalue')
                adist = vtx(:, 1) - acomOff(1) ;
                pts_to_remove = find(adist < adist0) ;
            elseif strcmpi(aCapMethod, 'yvalue')
                adist = vtx(:, 2) - acomOff(2) ;
                pts_to_remove = find(adist < adist0) ;
            elseif strcmpi(aCapMethod, 'zvalue')
                adist = vtx(:, 3) - acomOff(3) ;
                pts_to_remove = find(adist < adist0) ;
            end
            
            % Check it
            if preview
                clf
                trisurf(triangulation(mesh.f,vtx), ...
                    'edgecolor', 'none', 'facealpha', 0.1); hold on;
                scatter3(apt(1), apt(2), apt(3), 30, 'filled')
                scatter3(apt(1)+aOff(1), apt(2)+aOff(2), ...
                    apt(3)+aOff(3), 30, 'filled')
                scatter3(vtx(pts_to_remove, 1), vtx(pts_to_remove,2), ...
                    vtx(pts_to_remove, 3), 'markeredgealpha', 0.1)
                axis equal
                title('Points to remove for anterior face')
                pause(1)
                close all
            end
            
            % Make sure that we are removing a connected component
            if numel(pts_to_remove) > 2
                % Form a mesh from the piece(s) to be removed
                allpts = 1:size(vtx,1) ;
                all_but_acut = setdiff(allpts', pts_to_remove) ;
                [ acutfaces, acutvtx, ~] = remove_vertex_from_mesh( fv.f, fv.v, all_but_acut ) ;
                [ ~, ~, connected_indices, npieces ] = ...
                    remove_isolated_mesh_components( acutfaces, acutvtx ) ;
                if any(npieces > 1)
                    disp('Ensuring that only a single component is removed')
                    pts_to_remove = pts_to_remove(connected_indices) ;
                end

                % Check removed vertices
                % clf; trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none'); hold on; scatter3(fv.v(pts_to_remove, 1), fv.v(pts_to_remove, 2), fv.v(pts_to_remove, 3), 100)

                % Check removed component
                %trimesh(acutfaces, acutvtx(:, 1), acutvtx(:, 2), acutvtx(:, 3), npieces)

                % Could use Face-based connected component check:
                % Remove all faces... etc. Decided not to implement this
            else
                disp('Fewer than three vertices are being removed for anterior/proximal component. We assume this will be a connected component')
            end
            
            % Remove it
            [faces, vtx, keep_acut] = remove_vertex_from_mesh(fv.f, fv.v, pts_to_remove ) ;
            vn = fv.vn(keep_acut, :) ;
            nfaces_removed = size(fv.f, 1) - size(faces, 1) ;
            disp(['Removed ' num2str(length(pts_to_remove)) ' vertices and ' num2str(nfaces_removed) ' faces with acut'])
            clear nfaces_removed
            
            % CHECK that the output is correct topology
            % MATLAB-style triangulation
            meshTri = triangulation( faces, vtx );
            % The #Ex2 edge connectivity list of the mesh
            edgeTri = edges( meshTri );
            % Check that the input mesh is a topological cylinder
            if ( length(vtx) - length(edgeTri) + length(faces) ) == 1
                disp('Remaining mesh is topological disk')
            else
                disp('Remaining mesh after removing anterior cap is not disk! Removing smaller pieces...')
                % Check it
                % trisurf(triangulation(faces, vtx), 'edgecolor', 'none')
                [ faces, vtx, ~, npieces ] = ...
                    remove_isolated_mesh_components( faces, vtx ) ;
                meshTri = triangulation( faces, vtx );
                % The #Ex2 edge connectivity list of the mesh
                edgeTri = edges( meshTri );
                
                if ( length(vtx) - length(edgeTri) + length(faces) ) ~= 1
                    error('Topology after removing anterior cap is not disk')
                end
            end
            
            % Inspect
            % if preview
            %     fig = figure;
            %     trimesh(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3))
            %     view(30,145)
            %     waitfor(fig)
            % end
            
            %% Remove the posterior part as well
            
            % Measure distance to the posterior
            % Offset the posterior center of mass by the offset rate x time
            pOff = pOffXYZ + pOffRateXYZ * (tt - timePoints(1)) ;
            pcomOff = ppt + pOff ;
            disp(['Removing second/posterior endcap using method: ' aCapMethod])
            
            if strcmpi(pCapMethod, 'geodesic')    
                % Options for the marching
                fmoptions = struct() ;
                fmoptions.end_points = [];   
                
                [~, nearestID_p] = min(vecnorm(vtx - pcomOff, 2, 2)) ;

                %   output is the distance function to the set of starting points.
                pdist_g = perform_fast_marching_mesh(vtx, faces, nearestID_p, fmoptions);
                pdist2 = pdist_g.^2 ;
            elseif strcmpi(pCapMethod, 'ball')
                pdist2 = sum((vtx - pcomOff) .^ 2, 2);
            elseif strcmpi(pCapMethod, 'xvalue')
                pdist2 = (vtx(:, 1) - pcomOff(1)).^2 ;
            elseif strcmpi(pCapMethod, 'yvalue')
                pdist2 = (vtx(:, 2) - pcomOff(2)).^2 ;
            elseif strcmpi(pCapMethod, 'zvalue')
                pdist2 = (vtx(:, 3) - pcomOff(3)).^2 ;
            end
            
            % STRATEGY 1: within distance of pcom
            pdist_thres_ii = pdist_thres ;
            pcut_done = false ;
            while ~pcut_done
                disp(['finding points within ' num2str(pdist_thres_ii) ' of pcom'])
                pts_to_remove = find(pdist2 < pdist_thres_ii^2) ;
                
                loopIdx = 1 ;
                while numel(pts_to_remove) == 0 && loopIdx < 10
                    disp('Increasing radius by 10% since no points enclosed in pdist')
                    pts_to_remove = find(pdist2 < (pdist_thres_ii * (1 + 0.1 * loopIdx))^2) ;
                    loopIdx = loopIdx + 1 ;
                end
                
                if numel(pts_to_remove) == 0
                    error('Cannot find any piece of the mesh near p point even after increasing threshold ball radius by 10% 10 times')
                end
                
                % Remove the posterior cap and check if result is topological
                % cylinder
                [faces_postpcut, vtx_postpcut, keep_pcut] = ...
                    remove_vertex_from_mesh(faces, vtx, pts_to_remove ) ;
                
                % Check if result is cylinder topologically
                % MATLAB-style triangulation
                meshTri = triangulation( faces_postpcut, vtx_postpcut );
                % The #Ex2 edge connectivity list of the mesh
                edgeTri = edges( meshTri );
                % Check that the input mesh is a topological cylinder
                
                if ( length(vtx_postpcut) - length(edgeTri) + length(faces_postpcut) ) == 0
                    disp(['Remaining vertices are cylinder after pcut'])
                else
                    disp(['Selecting largest component of pcut only'])
                    % Make sure that we are removing a connected component
                    % form a mesh from the piece(s) to be removed
                    allpts = linspace(1, length(vtx), length(vtx)) ;
                    all_but_pcut = setdiff(allpts, pts_to_remove) ;
                    [ pcutfaces, pcutvtx, ~] = remove_vertex_from_mesh( faces, vtx, all_but_pcut ) ;
                    
                    [ pcutfacesLargest, pcutvtxLargest, connected_indices, npieces ] = ...
                        remove_isolated_mesh_components( pcutfaces, pcutvtx ) ;
                    % If there were more than one piece selected, remove the bigger one only
                    if any(npieces > 1)
                        pts_to_remove = pts_to_remove(connected_indices) ;
                    else
                        disp('Posterior cut includes only one region')
                    end
                    
                    % Remove the posterior cap
                    [faces_postpcut, vtx_postpcut, keep_pcut] = ...
                        remove_vertex_from_mesh(faces, vtx, pts_to_remove ) ;
                    
                    % % check it
                    % if preview
                    %     fig = figure ;
                    %     scatter3(vtx(:, 1), vtx(:, 2), vtx(:, 3))
                    %     hold on;
                    %     scatter3(vtx(pts_to_remove, 1), vtx(pts_to_remove, 2), vtx(pts_to_remove, 3),  'filled')
                    %     plot3(pcom(1), pcom(2), pcom(3), 's')
                    %     waitfor(fig)
                    % end
                end
                disp(['Removed ' num2str(length(pts_to_remove)) ' vertices with pcut'])
                % nremain = length(keep_pcut) ;
                % nbefore = length(keep_acut) ;
                % try
                %     assert(nbefore - nremain == length(pts_to_remove))
                % catch
                %     error('Pre vs post posterior cut vertices have different count')
                % end
                
                % Repeat if no vertices are removed
                if length(pts_to_remove) > 0
                    % Check also that result is cylinder topologically
                    % MATLAB-style triangulation
                    meshTri = triangulation( faces_postpcut, vtx_postpcut );
                    % The #Ex2 edge connectivity list of the mesh
                    edgeTri = edges( meshTri );
                    % Check that the input mesh is a topological cylinder
                    if ( length(vtx_postpcut) - length(edgeTri) + length(faces_postpcut) ) == 0
                        pcut_done = true ;
                    else
                        % Try simply isolating largest component -- this
                        % happens sometimes when there is an unreferenced
                        % vertex.
                        [ faces_postpcut, vtx_postpcut, ~, npieces ] = ...
                            remove_isolated_mesh_components( faces_postpcut, vtx_postpcut ) ;
                        
                        if ( length(vtx_postpcut) - length(edgeTri) + length(faces_postpcut) ) == 0
                            pcut_done = true ;
                        else
                            close all
                            trisurf(triangulation(pcutfaces, pcutvtx), 'edgecolor', 'k');
                            hold on;
                            trisurf(triangulation(faces_postpcut, vtx_postpcut), 'edgecolor', 'none')
                            axis equal
                            waitfor(gcf)
                            disp('BAD PCUT: NOT TOPOLOGICAL CYLINDER! Trying again with larger threshold')
                            pdist_thres_ii = pdist_thres_ii * 1.02 ;
                        end
                    end
                else
                    % repeat with larger threshold
                    pdist_thres_ii = pdist_thres_ii * 1.02 ;
                    if pdist_thres_ii > (max(vtx(:)) - min(vtx(:)))
                        error('Removing entire sample to attain correct topology. Address this.')
                    end
                end
            end
            faces = faces_postpcut ;
            vtx = vtx_postpcut ;
            vn = vn(keep_pcut, :) ;
            
            %% figure out which indices were kept
            keep = keep_acut(keep_pcut) ;
            
            %% Check that the remainder is a single connected component
            % currently have faces, vtx.
            [ faces, vtx, keep_final_pass, npieces ] = ...
                remove_isolated_mesh_components( faces, vtx ) ;
            nremain = length(keep_final_pass) ;
            nbefore = length(keep) ;
            disp(['Removed ' num2str(nbefore - nremain) ' vertices with final pass'])
            
            if any(npieces > 1)
                % Update keep here to reflect final pass
                keep = keep(keep_final_pass) ;
                vn = vn(keep_final_pass, :) ;
            end
            
            %% Save the data in units of pixels (same as original mesh)
            disp(['Saving cylinder mesh to ' outfn])
            plywrite_with_normals(outfn, faces, vtx * ssfactor, vn)
            
            % Save the indices to keep when cutting off endcaps
            save(keepfn, 'keep') ;
            
        end
    
    else
        
        meshcut = read_ply_mod(outfn) ;
        faces = meshcut.f ;
        vtx = meshcut.v / ssfactor ; 
        vn = meshcut.vn ;
        
    end
        
    %% Compute dorsal vertex on anterior and posterior free boundaries
    TR = triangulation(faces, vtx) ;
    boundary = freeBoundary(TR) ; 
    nn = length(boundary(:)) ;
    bb = zeros(nn, 1) ;
    bb(1) = boundary(1) ;
    dmyk = 2;
    for kk = 1:length(boundary)
        seg = boundary(kk, :) ;
        if ~any(bb == seg(1))
            bb(dmyk) = seg(1) ;
            dmyk = dmyk + 1;
        end
        if ~any(bb == seg(2))
            bb(dmyk) = seg(2) ;
            dmyk = dmyk + 1 ;
        end
    end
    bb = bb(bb > 0) ;
    
    % Check the boundary segments
    % for kk = 1:length(boundary)
    %     row = boundary(kk, :) ;
    %     plot3(vtx(row, 1), vtx(row, 2), vtx(row, 3), '-')
    %     hold on
    % end
    % Check the slimmed boundary
    % plot3(vtx(bb, 1), vtx(bb, 2), vtx(bb, 3), '-')

    % figure out if boundary is anterior or posterior
    % Note that this assumes an elongated structure so that anterior
    % boundary is closer to acom than pcom and vice versa
    adist = sum((vtx(bb, :) - apt) .^ 2, 2) ;
    pdist = sum((vtx(bb, :) - ppt) .^ 2, 2) ;
    ab = bb(adist < pdist) ;
    pb = bb(adist > pdist) ;
    % check them
    if preview
        fig = figure('visible', 'on');
        trisurf(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3),...
            'edgecolor', 'none', 'facealpha', 0.1)
        hold on;
        plot3(vtx(ab, 1), vtx(ab, 2), vtx(ab, 3), 'b.-')
        plot3(vtx(pb, 1), vtx(pb, 2), vtx(pb, 3), 'r.-')
        plot3(apt(1), apt(2), apt(3), 'o')
        plot3(ppt(1), ppt(2), ppt(3), 'o')
        axis equal
        pause(5)
        close all
    end
    
    % Save the anterior and posterior boundary indices (after endcap cut)
    save(boundaryfn, 'ab', 'pb')
    
    % choose anterior Dorsal as point with smallest phi for reference TP
    % (closest to zero or 2pi) if this is ref TP. Otherwise, match "prev"
    % Note that previous dorsal point could be that of the NEXT timepoint
    % if tt < tref. 
    if Dt == 0 
        % Choose dorsal point based on angle in yz plane. 
        % Since this is the first timepoint, taking angle with wrt x axis
        % is same as wrt AP axis given that AP axis is a straight line
        % presently.        
        % transform to APDV coords
        vrs = tubi.xyz2APDV(vtx * ssfactor) ;
        % subtract pi/2 to make dorsal be zero
        a_phipi = mod(atan2(vrs(ab, 3), vrs(ab, 2)) - pi * 0.5, 2*pi) ;
        % For posterior, subtract COM of the boundary in z -- this is a bit
        % of a hack, but sometimes the boundary is not well centered around
        % the posterior pole.
        if posterior_phi_subtractCOM
            p_phipi = mod(atan2(vrs(pb, 3) - mean(vrs(pb, 2)), vrs(pb, 2))...
                            - pi * 0.5, 2*pi) ;
        else
            p_phipi = mod(atan2(vrs(pb, 3), vrs(pb, 2)) - pi * 0.5, 2*pi) ;
        end
        % Now make the range go from -pi to pi with dorsal being zero
        a_phipi(a_phipi > pi) = a_phipi(a_phipi > pi) - 2 * pi ;
        p_phipi(p_phipi > pi) = p_phipi(p_phipi > pi) - 2 * pi ;
        [~, adb_tmp] = min(abs(a_phipi)) ;
        [~, pdb_tmp] = min(abs(p_phipi)) ;
        
        % anterior/posterior dorsal pt on boundary ['antr dorsal boundary']
        adb = ab(adb_tmp) ;
        pdb = pb(pdb_tmp) ;
        
        % check it
        if preview
            clf; set(gcf, 'visible', 'on')
            trisurf(faces, vtx(:, 1), vtx(:, 2), vtx(:, 3), 0*vtx(:, 1),...
                'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            scatter3(vtx(ab, 1), vtx(ab, 2), vtx(ab, 3), 10, a_phipi)
            scatter3(vtx(pb, 1), vtx(pb, 2), vtx(pb, 3), 10, p_phipi)
            plot3(vtx(adb, 1), vtx(adb, 2), vtx(adb, 3), 'ks')
            plot3(vtx(pdb, 1), vtx(pdb, 2), vtx(pdb, 3), 'k^')
            caxis([-pi, pi])
            colorbar()
            axis equal
            title('Identification of dorsal points through phi [pix/ssfactor]')
            pause(5)
            close all 
            
            % Plot in RS coords
            clf; set(gcf, 'visible', 'on')
            trisurf(faces, vrs(:, 1), vrs(:, 2), vrs(:, 3), 0*vrs(:, 1),...
                'edgecolor', 'none', 'facealpha', 0.1)
            hold on;
            scatter3(vrs(ab, 1), vrs(ab, 2), vrs(ab, 3), 10, a_phipi)
            scatter3(vrs(pb, 1), vrs(pb, 2), vrs(pb, 3), 10, p_phipi)
            plot3(vrs(adb, 1), vrs(adb, 2), vrs(adb, 3), 'ks')
            plot3(vrs(pdb, 1), vrs(pdb, 2), vrs(pdb, 3), 'k^')
            caxis([-pi, pi])
            colorbar()
            axis equal
            title('Identification of dorsal points through phi [microns, APDV]')
            pause(5)
            close all 
        end
        previous_avtx = vtx(adb, :) ;
        previous_pvtx = vtx(pdb, :) ;
        backward_avtx = vtx(adb, :) ;
        backward_pvtx = vtx(pdb, :) ;
    elseif Dt > 0 
        % CURRENT TIME IS A TIMEPOINT AFTER TREF
        % transform to APDV coords in case topology is wrong (for
        % inspection in case of error)
        vrs = tubi.xyz2APDV(vtx * ssfactor) ;
        
        % Match previous timepoint
        ka = dsearchn(vtx(ab, :), previous_avtx) ;
        kp = dsearchn(vtx(pb, :), previous_pvtx) ;
        adb = ab(ka) ;
        pdb = pb(kp) ;
        % UPDATE "PREVIOUS" to be CURRENT
        previous_avtx = vtx(adb, :) ;
        previous_pvtx = vtx(pdb, :) ;
    elseif Dt < 0 
        % CURRENT TIME IS A TIMEPOINT BEFORE TREF
        % transform to APDV coords in case topology is wrong (for
        % inspection in case of error)
        vrs = tubi.xyz2APDV(vtx * ssfactor) ;
        
        % Match previous timepoint
        ka = dsearchn(vtx(ab, :), backward_avtx) ;
        kp = dsearchn(vtx(pb, :), backward_pvtx) ;
        adb = ab(ka) ;
        pdb = pb(kp) ;        
        % UPDATE "PREVIOUS" (which here is in future) to be CURRENT
        backward_avtx = vtx(adb, :) ;
        backward_pvtx = vtx(pdb, :) ;
    end
        
    %% Save the anterior dorsal and posterior dorsal vertex points
    try 
        h5create(outapd_boundaryfn, ['/' name '/adorsal'], size(adb)) ;
    catch
        disp('adorsal pt already exists --> overwriting ') ;
    end
    try
        h5create(outapd_boundaryfn, ['/' name '/pdorsal'], size(pdb)) ;
    catch
        disp('pdorsal pt already exists --> overwriting') ;
    end
    h5write(outapd_boundaryfn, ['/' name '/adorsal'], adb) ;
    h5write(outapd_boundaryfn, ['/' name '/pdorsal'], pdb) ;
    
    %% Plot the result
    figfn = fullfile(figoutdir, [name '.png']) ;
    if save_figs && (~exist(figfn, 'file') || overwrite)
        disp(['Saving figure for frame ' num2str(ii)])
        fig = figure('Visible', 'Off');
        vrs = tubi.xyz2APDV(vtx * ssfactor) ;
        tmp = trisurf(faces, vrs(:, 1), vrs(:, 2), vrs(:, 3), ...
            vrs(:, 3), 'edgecolor', 'none', 'FaceAlpha', 0.4) ;
        hold on
        plot3(vrs(adb, 1), vrs(adb, 2), vrs(adb, 3), 's')
        plot3(vrs(pdb, 1), vrs(pdb, 2), vrs(pdb, 3), '^')
        acomrs = tubi.xyz2APDV(apt * ssfactor) ;
        pcomrs = tubi.xyz2APDV(ppt * ssfactor) ;
        acomOffrs = tubi.xyz2APDV(acomOff * ssfactor) ;
        pcomOffrs = tubi.xyz2APDV(pcomOff * ssfactor) ;
        plot3(acomrs(1), acomrs(2), acomrs(3), 'ks')
        plot3(pcomrs(1), pcomrs(2), pcomrs(3), 'k^')
        plot3(acomOffrs(1), acomOffrs(2), acomOffrs(3), 'ko')
        plot3(pcomOffrs(1), pcomOffrs(2), pcomOffrs(3), 'ko')
        axis equal
        xlim(xyzlim_um(1, :) + [-5, 5])
        ylim(xyzlim_um(2, :) + [-5, 5])
        zlim(xyzlim_um(3, :) + [-5, 5])
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex')
        titlestr = 'Mesh with cylindrical cuts, $t=$' ;
        timestr = sprintf('%03d', tt * tubi.timeInterval) ;
        titlestr = [titlestr timestr ' ' tubi.timeUnits] ;
        title(titlestr, 'Interpreter', 'Latex')
        legend({'surface', 'cutPath start', 'cutPath end', ...
            'anterior', 'posterior', ...
            'anterior point with offset', 'posterior point with offset'})
        %view(50, 145)
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=16cm
        
        % view(0,0)
        saveas(fig, figfn)
        
        view(2)
        figfn = fullfile(figoutdir, ['xy_' name '.png']) ;
        saveas(fig, figfn)
       
        % view(90,0)
        % figfn = fullfile(figoutdir, ['yz_' name '.png']) ;
        % saveas(fig, figfn)

        if preview
            set(fig, 'Visible', 'On') ;
            waitfor(fig) ;
        else
            close(fig)
        end
    end
    

    %% Verify output triangulation topology -------------------------------------

    % MATLAB-style triangulation
    meshTri = triangulation( faces, vrs );
    % The #Ex2 edge connectivity list of the mesh
    edgeTri = edges( meshTri );
    % Check that the input mesh is a topological cylinder
    eulerChar = ( length(vrs) - length(edgeTri) + length(faces) ) ;
    numBdy = numel(DiscreteRicciFlow.compute_boundaries(faces));
    if ((eulerChar ~= 0) || (numBdy ~= 2))
        trisurf(meshTri, 'FaceColor', 'none')
        axis equal
        title(['Input mesh is NOT a topological cylinder: TP=' num2str(tt)]);
        disp(['Euler characteristic = ', num2str(eulerChar)]);
        disp(['Number of boundaries = ', num2str(numBdy)]);
        error( 'Input mesh is NOT a topological cylinder!' );
    end

end

disp('done')
