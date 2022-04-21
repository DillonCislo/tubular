classdef TubULAR < handle
    % Tube-like sUrface Lagrangian Analysis Resource (TubULAR) class
    %
    % Coordinate Systems
    % ------------------    
    % uv :  (conformal map onto unit square)
    %       Conformally mapping the cylinderCutMesh onto the unit square 
    %       in the plane results in the instantaneous uv coordinate system. 
    %       Corners of the unit square are taken directly from the cutMesh,
    %       so are liable to include some overall twist.
    % sphi : (proper length x rectified azimuthal coordinate)
    %       quasi-axisymmetric system in which first coordinate is 
    %       proper length along surface and second is a rectified azimuthal
    %       coordinate. Rectification means that the surface is rotated as
    %       v->phi(s), where each s coordinate is offset by a "rotation"
    %       about the surface along a direction that is perpendicular to s
    %       in pullback space. This rectification may be based on surface
    %       positions in R^3 (geometric) or based on intensity motion in
    %       pullback space (material/Lagrangian) inferred through 
    %       phasecorrelation of tissue strips around discretized s values.
    % uvprime : (conformal map
    %       [same as uvprime_sm, since uvprime is currently computed via
    %       sphi_sm coordinates]
    % r/spr/sphir :  (proper length x rectified azimuthal coordinate)
    %       same as sphi but with aspect ratio relaxed to minimize isoareal
    %       energy cost --> makes triangles more similar in area by scaling
    %       the s axis (longitudinal axis) by a scalar factor.
    %
    % iLastik's coordinate systems : I recommend using cxyz as the axis 
    %       when outputting from iLastik. For our setup, this will mean
    %       that the meshes are mirrored with respect to the lab frame, so
    %       tubi.flipy = true.
    % 
    % PIV measurements fall into two classes: 
    %   - 'piv': principal surface-Lagrangian-frame PIV (sp_sme or up_sme)
    %       --> note that the designation of coordinate system is not
    %       explicitly specified in the filenames:
    %       tubi.dir.mesh/gridCoords_nU0100_nV0100/piv/piv3d, etc
    %   - 'piv_uvp_sme': PIV in coordSys for quasiconformal measurements 
    %       --> note that these are less Lagrangian than
    %       sp_sme or up_sme, so they are treated as independent from the 
    %       principal pipeline in which we measure velocities in a
    %       surface-Lagrangian frame (ie sp_sme or up_sme).
    %
    % Properties
    % ----------
    % imSize        : 2x1 int, size of pullback images to create
    % xyzlim        : 3x2 float, mesh limits in full resolution pixels, in data space
	% xyzlim_um     : 3x2 float, mesh limits in lab APDV frame in microns
    % resolution    : float, resolution of pixels in um
    % rot           : 3x3 float, APDV rotation matrix
    % trans         : 3x1 float, APDV translation 
    % a_fixed       : float, aspect ratio for fixed-width pullbacks
    % phiMethod     : str, '3dcurves' or 'texture'
    % flipy         : bool, APDV coord system is mirrored XZ wrt raw data
    % nV            : int, sampling number along circumferential axis
    % nU            : int, sampling number along longitudinal axis
    % uvexten       : str, naming extension with nU and nV like '_nU0100_nV0100'
    % t0            : int, reference timePoint in the experiment
    % velocityAverage : struct with fields
    %   vsmM : (#timePoints-1) x (nX*nY) x 3 float array
    %       3d velocities at PIV evaluation coordinates in um/dt rs
    %   vfsmM : (#timePoints-1) x (2*nU*(nV-1)) x 3 float array
    %       3d velocities at face barycenters in um/dt rs
    %   vnsmM : (#timePoints-1) x (nX*nY) float array
    %       normal velocity at PIV evaluation coordinates in um/dt rs
    %   vvsmM : (#timePoints-1) x (nU*nV) x 3 float array
    %       3d velocities at (1x resolution) mesh vertices in um/min rs
    %   v2dsmM : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in pixels/ min
    %   v2dsmMum : (#timePoints-1) x (nX*nY) x 2 float array
    %       2d velocities at PIV evaluation coordinates in scaled pix/min, but 
    %       proportional to um/min (scaled by dilation of map)
    %
    % APDV coordinate system and Centerline specification
    % ---------------------------------------------------
    % QuapSlap allows the user to designate the AP axis, a DV axis, and a
    % centerline for the segmented object. 
    % To designate APDV coordinate system, use iLastik training on the
    % timepoint t0 (which is an attribute of tubi). In particular, train on
    % anterior (A), posterior (P), background (B), and dorsal (D) location 
    % in different iLastik channels. Small blobs of high probability near
    % the anterior end for A, somewhere along the posterior end for P, and 
    % anywhere along dorsal for D are best. Then the centers of mass of
    % thresholded contiguous high probability are computed for each to
    % define the coordinate system. By default, the ilastik results are
    % read as ch1=A, ch2=P, ch3=B, ch4=D, but anteriorChannel, 
    % posteriorChannel, and dorsalChannel specify the iLastik
    % training channel that is used for each specification.
    % Name the h5 file output from iLastik as 
    % ..._Probabilities_APDVcoords.h5
    % For example, dorsal for the gut was chosen at the fused site where 
    % additional 48YGAL4-expressing muscle-like cells form a seam.
    % Posterior is at the rear of the yolk, where the endoderm closes, for 
    % apical surface training. Anterior is at the junction of the midgut 
    % with the foregut.
    % Separately, define the AP points for centerline extraction. For most gut
    % data, the posterior point is different in centerline than it is for AP
    % centerline specification, so we use a different ILP to train for A and P
    % for all timepoints of dynamic data. 
    %
    % 
    properties
        xp                      % ImSAnE experiment class instance or struct with fields
                                %   expMeta : struct with fields
                                %   fileMeta : struct with fields
        dynamic                 % true if multiple timepoints, false if fixed
        timeInterval = 1        % increment in time between timepoints with 
                                % indices differing by 1. For example, if
                                % timePoints are [0,1,2,4] and these are 
                                % [1,1,2] minutes apart, then timeInterval 
                                % is 1. 
        timeUnits = 'min'       % units of the timeInterval (ex 'min')
        spaceUnits = '$\mu$m'   % units of the embedding space (ex '$\mu$m')
        imSize                  % size of pullback images to create (default is [a_ratio * 1000, 1000])
        dir                     % str, directory where QuapSlap data lives
        fileName                % fileName
        fileBase                % fileNames to be populated by timestamp ('...%06d...mat')
        fullFileBase            % full path of filenames (like fullfile(tubi.dir.X, tubi.fileBase.X))
        ssfactor                % subsampling factor for probabilities 
        APDV = struct(...
            'resolution', [], ...       % resolution of data in spaceUnits / pixel
            'rot', [], ...              % rotation matrix to transform data into APDV frame (rot*v+trans)*resolution
            'trans', [])                % translation vector to transform data into APDV frame (rot*v+trans)*resolution
        flipy                           % whether data is mirror image of lab frame coordinates
        nV                              % sampling number along circumferential axis
        nU                              % sampling number along longitudinal axis
        uvexten                         % naming extension with nU and nV like '_nU0100_nV0100'
        t0                              % reference time in the experiment
        normalShift = 0                 % shift to apply to meshes in pixel space along normal direction
        a_fixed = 1                     % aspect ratio for fixed geometry pullback meshes
        phiMethod = '3dcurves'          % method for determining Phi map in pullback mesh creation, with 
                                        % the full map from embedding to pullback being [M'=(Phi)o()o()]. 
                                        % This string specifier must be '3dcurves' (geometric phi stabilization) 
                                        % or 'texture' (optical flow phi stabilization)
        endcapOptions                   % struct with fields adist, pdist, and tref.
                                        %   adist : distance around anterior point A which is removed from tubular mesh (sliced off)
                                        %   pdist : distance around anterior point P which is removed from tubular mesh (sliced off) 
                                        %   tref : reference time used to define th point on the endcap
                                        %       at which we cut the cylinder mesh into a cylinderCutMesh (a topological disk/square). 
                                        %       This "dorsal" point for other timepoints are identified by pointmatching.
                                        %   Additional fields allowed : 
        plotting = struct('preview', false, ... % display intermediate results
            'save_ims', true, ...       % save images
            'xyzlim_um_buff', [], ...   % xyzlimits in um in RS coord sys with buffer
            'xyzlim_raw', [], ...       % xyzlimits in pixels
            'xyzlim_pix', [], ...       % xyzlimits in pixels RS
            'xyzlim_um', [], ...        % xyzlimits in um in RS coord sys
            'colors', [])               % color cycle for tubi
        apdvPts = struct('anteriorPts', [], ...
            'posteriorPts', [], ... 
            'antPts_sm', [], ...
            'postPts_sm', [], ... 
            'dorsalPts', [], ... 
            'antPts_rs', [], ... 
            'postPts_rs', [], ... 
            'dorsPts_rs', [])
        apdvOptions
        currentTime
        currentMesh = struct('rawMesh', [], ... % original mesh found by surface detection
            'alignedMesh', [], ...      % APDV rotated and scaled mesh (raw mesh in APDV coordinates)     
            'cylinderMesh', [], ...     % original mesh with endcaps cut off
            'cylinderMeshClean', [], ...    % cylinder mesh with "ears" removed (ears give difficulty in mapping to the plane)
            'cutMesh', [], ...          % cylinder mesh with a seam given by cutPath
            'cutPath', [], ...          % vertex indices of the cutMesh along which the periodic seam is cut
            'uvcutMesh', [], ...        % rectilinear cutMesh in (u,v) from Dirichlet map result to rectangle 
            'spcutMesh', [], ...        % rectilinear cutMesh in (s,phi) 'surface Lagrangian' parameterization
            'spcutMeshSm', [], ...      % rectilinear cutMesh in (s,phi) smoothed in time
            'spcutMeshSmRS', [], ...    % rectilinear cutMesh in (s,phi) smoothed in time with rotated scaled embedding
            'spcutMeshSmRSC', [], ...   % rectilinear cutMesh as closed cylinder (topological annulus), in (s,phi) smoothed, with rotated scaled embedding
            'ricciMesh', [])            % ricci flow result pullback mesh, topological annulus          
        currentCline = struct('mss', [], ...
            'mcline', [], ...
            'avgpts', []) ;
        data = struct('adjustlow', 0, ...
            'adjusthigh', 0, ...
            'axisOrder', [1 2 3], ...
            'ilastikOutputAxisOrder', 'cxyz') % options for scaling and transposing image intensity data
        currentData = struct('IV', [], ...
            'adjustlow', 0, ...
            'adjusthigh', 0 )           % image intensity data in 3d and scaling
        currentVelocity = struct('piv3d', struct()) ;     
        piv = struct( ...
            'imCoords', 'sp_sme', ...   % image coord system for measuring PIV / optical flow) ;
            'Lx', [], ...               % width of image, in pixels (x coordinate)
            'Ly', [], ...               % height of image, in pixels (y coordinate)
            'raw', struct(), ...        % raw PIV results from disk/PIVLab
            'smoothed', struct(), ...   % smoothed PIV results after gaussian blur
            'smoothing_sigma', 1 ) ;    % sigma of gaussian smoothing on PIV, in units of PIV sampling grid pixels
        velocityAverage = struct(...
            'v3d', [], ...              % 3D velocities in embedding space [pix/dt]
            'v2d', [], ...              % 2D tangential velocities in pullback
            'v2dum', [], ...            % 2D tangential velocity scaled by speed in true embedding space
            'vn', [], ...               % normal velocity in spaceUnits per timeInterval timeUnits
            'vf', [], ...               % velocity vielf on face barycenters after Lagrangian avg
            'vv', []) ;                 % velocity field on vertices after Lagrangian avg
        cleanCntrlines                  % centerlines in embedding space after temporal averaging
        smoothing = struct(...
            'lambda', 0.00, ...             % diffusion const for field smoothing on mesh
            'lambda_mesh', 0.00, ...        % diffusion const for vertex smoothing of mesh itself
            'nmodes', 7, ...                % number of low freq modes to keep per DV hoop
            'zwidth', 1) ;                  % half-width of tripulse filter applied along zeta/z/s/u direction in pullback space, in units of du/dz/ds/dzeta
        pathlines = struct('t0', [], ...    % timestamp (not an index) at which pathlines form regular grid in space
            'refMesh', [], ...              % reference mesh for pathline advection
            'piv', [], ...                  % Lagrangian pathlines from piv coords
            'vertices', [], ...             % Lagrangian pathlines from mesh vertices
            'vertices3d', [], ...           % Lagrangian pathlines from mesh vertices
            'faces', [], ...                % Lagrangian pathlines from mesh face barycenters
            'beltrami', ...                 % beltrami coefficient evaluated along pathlines
                struct('mu_material', [], ...
                'mu_material_filtered', [], ...
                'mu_material_vertices', [], ...
                'fitlerOptions', []));     
       currentStrain = struct(...
            'pathline', ...                 % strain from pathlines
            struct('t0Pathlines', [], ...   % t=0 timepoint for pathlines in question
            'strain', [], ...               % strain evaluated along pathlines
            'beltrami', [])) ;              % beltrami coefficient for pathlines
    end
    
    % Some methods are hidden from public view. These are used internally
    % to the class.
    methods (Hidden)
        function tubi = TubULAR(xp, opts)
            tubi.initializeTubULAR(xp, opts)
        end
        initializeTubULAR(tubi, xp, opts)
        plotSPCutMeshSmSeriesUtility(tubi, coordsys, options)
        plotMetricKinematicsTimePoint(tubi, tp, options)
        [XX, YY] = pullbackPathlines(tubi, x0, y0, t0, options) 
        plotAverageVelocitiesTimePoint(tubi, tp, options)
        plotPathlineVelocitiesTimePoint(tubi, tp, options)
        plotStrainRateTimePoint(tubi, tp, options) 
        plotPathlineStrainRateTimePoint(tubi, tp, options)
        plotPathlineStrainTimePoint(tubi, tp, options)
    end
    
    % Public methods, accessible from outside the class and reliant on 
    % properties of the class instance
    methods
        function setTime(tubi, tt)
            % Set the current time of the dataset and clear current data
            % which was associated with the previously considered time
            %
            % Parameters
            % ----------
            % tt : int or float
            %   timePoint to set to be current, from available times in
            %   tubi.xp.fileMeta.timePoints
            %
            if tt ~= tubi.currentTime
                tubi.clearTime() ;
            end
            tubi.currentTime = tt ;
            if ~isa(tubi.xp, 'struct')
                % disp('WARNING: tubi.xp is interpreted as a struct rather than imsane Experiment() class instance')
                tubi.xp.setTime(tt) ;
            end
        end
        
        function [tpTrue, timestr ] = trueTime(tubi, tp, div60)
            if nargin < 3
                div60 = false ;
            end
            tpTrue = (tp - tubi.t0set()) * tubi.timeInterval ;
            if div60
                tpTrue = tpTrue / 60 ;
            end
            
            if nargout > 1
                if div60
                    if contains(tubi.timeUnits, 'min')
                        timestr = [num2str(tpTrue) ' hr'] ;
                    elseif contains(tubi.timeUnits, 'sec')
                        timestr = [num2str(tpTrue) ' min'] ;
                    end
                else
                    timestr = [sprintf('%03d', tpTrue) ' ' tubi.timeUnits] ;
                end
            end
        end
        
        function clearTime(tubi)
            % clear current timepoint's data for tubi instance
            tubi.currentTime = [] ;
            tubi.currentMesh.rawMesh = [] ;
            tubi.currentMesh.rawMeshRS = [] ;
            tubi.currentMesh.cylinderMesh = [] ;
            tubi.currentMesh.cylinderMeshClean = [] ;
            tubi.currentMesh.cutMesh = [] ;
            tubi.currentMesh.cutPath = [] ;
            tubi.currentMesh.spcutMesh = [] ;
            tubi.currentMesh.cutMesh = [] ;
            tubi.currentMesh.spcutMesh = [] ;
            tubi.currentMesh.spcutMeshSm = [] ;
            tubi.currentMesh.spcutMeshSmRS = [] ;
            tubi.currentMesh.spcutMeshSmRSC = [] ;
            tubi.currentMesh.uvpcutMesh = [] ;
            tubi.currentData.IV = [] ;
            tubi.currentData.adjustlow = 0 ;
            tubi.currentData.adjusthigh = 0 ;
            tubi.currentVelocity.piv3d = struct() ;
        end
        
        function t0 = t0set(tubi, t0)
            % t0set(tubi, t0) Set time offset from file or manually 
            if nargin < 2
                if exist(tubi.fileName.t0, 'file')
                    % Note that fold_onset is in units of timepoints, not 
                    % indices into timepoints
                    dlmread(tubi.fileName.t0, ',', 1, 0) ;
                else
                    tubi.t0 = tubi.xp.fileMeta.timePoints(1) ;
                end
            else
                tubi.t0 = t0 ;
            end
            t0 = tubi.t0 ;
        end
        
        function makeMIPs(tubi, dim, pages, timePoints, adjustIV)
            if nargin < 5
                adjustIV = false ;
            end
            if nargin < 4 
                timePoints = tubi.xp.fileMeta.timePoints;
            elseif isempty(timePoints)
                timePoints = tubi.xp.fileMeta.timePoints;
            end
            if ~iscell(pages)
                pages = {pages} ;
            end
            % create mip directories if needed
            for qq = 1:length(pages)
                outdir = sprintf(tubi.dir.mip, dim, ...
                            min(pages{qq}), max(pages{qq})) ;
                if ~exist(outdir, 'dir')
                    mkdir(outdir)
                end
            end
            
            % make the mips
            for tp = timePoints
                for qq = 1:length(pages)
                    im = tubi.mip(tp, dim, pages{qq}, adjustIV) ;
                    imfn = sprintf(tubi.fullFileBase.mip, dim, ...
                        min(pages{qq}), max(pages{qq}), tp) ;
                    imwrite(im, imfn,'tiff','Compression','none')
                end
            end
            disp('done making mips')
        end
        
        function im = mip(tubi, tp, dim, pages, adjustIV)
            if nargin < 5
                adjustIV = false ;
            end
            tubi.setTime(tp)
            tubi.getCurrentData(adjustIV)
            for qq = 1:length(tubi.currentData.IV)
                if dim == 1
                    im = squeeze(max(tubi.currentData.IV{qq}(pages, :, :), [], dim)) ;
                elseif dim == 2
                elseif dim == 3
                else
                    error('dim > 3 not understood')
                end
            end
        end
        
        prepareIlastik(tubi, opts)
        
        function getMeshes(tubi, detectOptions)
            % Obtain mesh surfaces of volumetric data (like ImSAnE's surface
            % detection methods), here using integralDetector methods for
            % activecontouring
            if isfield(detectOptions, 'pressure')
                pressure = detectOptions.pressure ;
            else
                pressure = 0 ;
            end
            if isfield(detectOptions, 'tension')
                tension = detectOptions.tension ;
            else
                tension = 0 ;
            end
            if isfield(detectOptions, 'pre_pressure')
                pre_pressure = detectOptions.pre_pressure ;
            else
                pre_pressure = 0 ;
            end
            if isfield(detectOptions, 'pre_tension')
                pre_tension = detectOptions.pre_tension ;
            else
                pre_tension = 0 ;
            end
            if isfield(detectOptions, 'post_pressure')
                post_pressure = detectOptions.post_pressure ;
            else
                post_pressure = 0 ;
            end
            if isfield(detectOptions, 'post_tension')
                post_tension = detectOptions.post_tension ;
            else
                post_tension = 0 ;
            end
            opts = detectOptions ;
            
            
            
            % ilastik internally swaps axes. 1:x, 2:y, 3:z, 4:class
            % strategy: put into xyzc format, then pop last index
            if strcmp(opts.ilastikaxisorder, 'xyzc')
                pred = file ;
            elseif strcmp(opts.ilastikaxisorder, 'yxzc')
                % to convert yxzc to xyzc, put x=2 y=1 z=3 c=4
                pred = permute(file,[2,1,3,4]);
            elseif strcmp(opts.ilastikaxisorder, 'zyxc')
                % to convert yxzc to xyzc, put x=3 y=2 z=1 c=4
                pred = permute(file,[3,2,1,4]);
            elseif strcmp(opts.ilastikaxisorder, 'yzcx')
                % to convert yxzc to xyzc, put x=4 y=1 z=2 c=3
                pred = permute(file,[4,1,2,3]);
            elseif strcmp(opts.ilastikaxisorder, 'cxyz')
                % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
                pred = permute(file,[2,3,4,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cyxz')
                % to convert yxzc to xyzc, put x=2 y=3 z=4 c=1
                pred = permute(file,[3,2,4,1]);
            elseif strcmp(opts.ilastikaxisorder, 'czyx')
                % to convert yxzc to xyzc, put x=1>4 y=2>3 z=3>2 c=4>1
                pred = permute(file,[4,3,2,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cyzx')
                % to convert cyzx to xyzc put x=1>4 y=2>2 z=3>3 c=4>1
                pred = permute(file,[4,2,3,1]);
            elseif strcmp(opts.ilastikaxisorder, 'cxzy')
                % to convert cxzy to xyzc put x=1>2 y=2>4 z=3>3 c=4>1
                pred = permute(file,[2,4,3,1]);
            else
                error('Have not coded for this axisorder. Do so here')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Identify the surface using the loaded probabilities here
            % Convert the current image to a level set using morphological snakes
            ssfactor = opts.ssfactor ;
            niter = opts.niter ;
            niter0 = opts.niter0 ;
            mslsDir = opts.mslsDir ;
            pressure = opts.pressure ;
            tension = opts.tension ;
            exit_thres = opts.exit_thres ;
            ofn_ply = opts.ofn_ply ;
            ofn_ls = opts.ofn_ls ;
            channel = opts.channel ;
            ms_scriptDir = opts.ms_scriptDir ;
            % tpstamp = fileMeta.timePoints(first_tp);
            timepoint = opts.timepoint ;
            ofn_smoothply = opts.ofn_smoothply ;
            mlxprogram = opts.mlxprogram;
            init_ls_fn = opts.init_ls_fn ;
            % Run MS on a series of timepoints in a parent directory
            run_full_dataset = opts.run_full_dataset ;
            % Radius of initial guess if init_ls_fn does not exist or is
            % not supplied
            radius_guess = opts.radius_guess ;
            save = opts.save ;
            center_guess = opts.center_guess ;
            plot_mesh3d = opts.plot_mesh3d ;
            dtype = opts.dtype ; 
            mask = opts.mask ;
            use_pointcloud = opts.mesh_from_pointcloud ;
            dataset_prob_searchstr = opts.prob_searchstr ;
            ilastikaxisorder = opts.ilastikaxisorder ;
            smooth_with_matlab = opts.smooth_with_matlab ;
            % since python flips axes wrt MATLAB, flip them here
            morphsnakesaxisorder = fliplr(ilastikaxisorder) ;
            
            % Create the output dir if it doesn't exist
            if ~exist(mslsDir, 'dir')
                mkdir(mslsDir)
            end
            
            outputMesh = fullfile(mslsDir, msls_mesh_outfn) ;
            
            % Check if previous time point's level set exists to use as a seed
            % First look for supplied fn from detectOptions.
            % If not supplied (ie init_ls_fn is none or empty string, then
            % seek previous timepoint output from MS algorithm.
            
            disp(['init_ls_fn = ', init_ls_fn])
            disp(['ofn_ls = ', ofn_ls])
            if strcmp(init_ls_fn, 'none') || strcmp(init_ls_fn, '')
                % User has NOT supplied fn from detectOptions
                init_ls_fn = [ofn_ls, ...
                    num2str(timepoint - 1, '%06d' ) '.' dtype] ;
            end
            
            disp([ 'initial level set fn = ', init_ls_fn])
            if exist(init_ls_fn, 'file')
                % It does exist, and the given name is the RELATIVE path.    
                % Use it as a seed (initial level set) 
                disp('running using initial level set')
                init_ls = load(init_ls_fn) ;
            elseif exist(fullfile(mslsDir, init_ls_fn), 'file')
                % It does exist, and given name is the relative path
                % without the extension. 
                % Use it as a seed (initial level set)
                disp('running using initial level set')
                init_ls = load(fullfile(mslsDir, init_ls_fn)) ;
            elseif exist(fullfile(mslsDir, [ init_ls_fn '.h5']), 'file')
                % It does exist, and given name is the FULL path 
                % without the extension.
                disp('running using initial level set')
                init_ls = load(fullfile(mslsDir, init_ls_fn)) ;                
            else
                % The guess for the initial levelset does NOT exist, so use
                % a sphere for the guess.
                disp(['Using default sphere for init_ls -- no such file on disk: ' fullfile(mslsDir, [ init_ls_fn '.h5'])])
                init_ls = [] ;
            end
            
            % Flip axis order LR of output mesh to return to MATLAB
            % orientation? NO, not helpful to do "-permute_mesh 'zyx'"
            % since already did morphsnakesaxisorder = fliplr() earlier
            % command = [command ' -adjust_for_MATLAB_indexing'] ;
            
            disp(['does outputMesh exist: ', num2str(exist(outputMesh, 'file'))])
            disp(outputMesh)
            % error('here')
            if ~exist(outputMesh, 'file')

                if use_dataset_command
                    % User has elected to run as a dataset, so pass a directory
                    % with _Probabilities.h5 files to run on.
                    error('handle here')
                    msls_mesh_outfn = ofn_ply;
                    ls_outfn = ofn_ls;
                else
                    % We are running MS on a single file. Give the filename of
                    % the ilatik output run on filename.h5
                    prob_infn = [opts.fileName, '_Probabilities.h5'] ;
                    msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d'), '.ply'];
                    ls_outfn = [ofn_ls, num2str(timepoint, '%06d'), '.', dtype];

                    BW = activecontour(data, init_ls, niter,...
                        'SmoothFactor', tension, 'ContractionBias', pressure) ;
                    
                    % Convert BW to mesh
                    mesh = isosurface(BW, 0.5) ;
                end

            else
                disp(['output PLY already exists: ', msls_mesh_outfn])
            end
            
            %% Clean up mesh file for this timepoint using MeshLab --------

            % Here use the boundary mesh from marching cubes to make a
            % smooth mesh
            % Check if we need to smooth the full dataset of meshes
            if use_dataset_command
                disp('Using dataset command for morphsnakes...')
                % find all ms_...ply files in mslsDir, and smooth them all
                files_to_smooth = dir(fullfile(mslsDir, [ofn_ply '*.ply'])) ;
                lsfns_to_smooth = dir(fullfile(mslsDir, [ls_outfn '*' dtype])) ;
                
                if abs(length(lsfns_to_smooth) - length(files_to_smooth))>1 
                    error('The number of output levelsets does not equal the number of output unsmoothed PLYs. These must match.')
                end
                for i=1:length(files_to_smooth)
                    msls_mesh_outfn = files_to_smooth(i).name ;
                    infile = fullfile( mslsDir, msls_mesh_outfn );
                    % Note that LS file is outputLs ;
                    split_fn = strsplit(msls_mesh_outfn, ofn_ply) ;
                    extension_outfn = split_fn{2} ;
                    base_outfn_for_pointcloud = ofn_ply ;
                    mesh_outfn = [ofn_smoothply, extension_outfn] ;
                    outputMesh = fullfile(mslsDir, mesh_outfn);
                    
                    disp(['outputMesh = ', outputMesh])
                    %bad = so_bad
                    if ~exist( outputMesh, 'file')
                        % Smooth with either meshlab or matlab
                        if smooth_with_matlab < 0
                            % USE MESHLAB, not matlab
                            if use_pointcloud
                                % Use the pointcloud from the level set rather than the
                                % boundary mesh from marching cubes
                                %----------------------------------------------------------------------
                                % Extract the implicit level set as a 3D binary array
                                %----------------------------------------------------------------------

                                % The file name of the current time point
                                ls_outfn_ii = lsfns_to_smooth(i).name ;
                                % The 3D binay array
                                bwLS = h5read( ls_outfn_ii, '/implicit_levelset' );

                                % Extract the (x,y,z)-locations of the level set boundary (in pixel
                                % space)
                                bwBdyIDx = bwperim( bwLS );

                                clear bwBdy
                                [ bwBdy(:,1), bwBdy(:,2), bwBdy(:,3) ] = ind2sub( size(bwLS), ...
                                    find(bwBdyIDx) );

                                %----------------------------------------------------------------------
                                % Create output mesh
                                %----------------------------------------------------------------------

                                % Write the points to a .obj file as a point cloud for ouput to Meshlab
                                clear OBJ
                                OBJ.vertices = bwBdy;
                                OBJ.objects(1).type='f';
                                OBJ.objects(1).data.vertices=[];

                                pointcloud_fn = [base_outfn_for_pointcloud '.obj'] ;
                                disp(['Writing point cloud ' pointcloud_fn]);
                                write_wobj(OBJ, pointcloud_fn );

                                % Run the meshlab script
                                command = ['meshlabserver -i ' pointCloudFileName, ...
                                    ' -o ' outputMesh, ' -s ' mlxprogram ' -om vn'] ;
                                disp(['running ' command])
                                system( command );
                            else
                                % Use the marching cubes mesh surface to smooth
                                command = ['meshlabserver -i ' infile ' -o ' outputMesh, ...
                                    ' -s ' mlxprogram ' -om vn'];
                                % Either copy the command to the clipboard
                                clipboard('copy', command);
                                % or else run it on the system
                                disp(['running ' command])
                                system(command)
                            end
                        elseif smooth_with_matlab == 0
                            disp('No smoothing, with either matlab or meshlab')
                            mesh = read_ply_mod(infile) ;
                            disp('Compute normals...')
                            mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                            plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                        elseif smooth_with_matlab > 0 
                            % Smooth with MATLAB
                            disp(['Smoothing with MATLAB using lambda = given value of ' num2str(smooth_with_matlab)])
                            mesh = read_ply_mod(infile) ;

                            % Check that this behaves the way we want
                            % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                            % size(mesh.vn)
                            % size(mesh.v)
                            % assert(size(mesh.vn, 1) == size(mesh.v, 1))

                            newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
                            disp('Compute normals...')
                            mesh.v = newV ;
                            mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                            disp(['Saving smoothed mesh to ' outputMesh])
                            plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                        end
                    else
                        disp(['t=', num2str(timepoint) ': smoothed mesh file found...'])
                        disp([' --> file to smooth was ' files_to_smooth(i).name])
                    end

                end
            else
                disp('Using individual timepoint command for morphsnakes')
                mesh_outfn = [ofn_smoothply, num2str(timepoint, '%06d'), '.ply'];
                outputMesh = fullfile(mslsDir, mesh_outfn) ;
                if ~exist( outputMesh, 'file')
                    msls_mesh_outfn = [ofn_ply, num2str(timepoint, '%06d' ), '.ply'];
                    infile = fullfile( mslsDir, msls_mesh_outfn );
                    if smooth_with_matlab < 0
                        % Build meshlab command to smooth meshes
                        command = ['meshlabserver -i ' infile ' -o ' outputMesh, ...
                            ' -s ' mlxprogram ' -om vn'];
                        % Either copy the command to the clipboard
                        clipboard('copy', command);
                        % or else run it on the system
                        disp(['running ' command])
                        system(command)
                    elseif smooth_with_matlab == 0
                        disp('No smoothing, with either matlab or meshlab')
                         mesh = read_ply_mod(infile) ;
                         disp('Compute normals...')
                         mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                         plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)      
                    elseif smooth_with_matlab > 0 
                        disp(['Smoothing with MATLAB using lambda = given value of ' num2str(smooth_with_matlab)])
                        mesh = read_ply_mod(infile) ;
                        
                        % Check that this behaves the way we want
                        % mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                        % size(mesh.vn)
                        % size(mesh.v)
                        % assert(size(mesh.vn, 1) == size(mesh.v, 1))
                        
                        newV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], smooth_with_matlab) ;
                        disp('Compute normals...')
                        mesh.v = newV ;
                        mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
                        disp(['Saving smoothed mesh to ' outputMesh])
                        plywrite_with_normals(outputMesh, mesh.f, mesh.v, mesh.vn)
                    end
                else
                    disp(['t=', num2str(timepoint) ': smoothed mesh file found, loading...'])    
                end    
            end
            
        end
        
        % Obtain global coordinate system for viewing (rotated coord sys)
        [apt,ppt,dpt] = computeAPDVCoords(tubi, opts)
        
        function [apts_sm, ppts_sm] = getAPpointsSm(tubi) 
            % Load the anterior and posterior 'centers of mass' ie the
            % endpoints of the object's centerline
            try
                apts_sm = h5read(tubi.fileName.apdv, '/apts_sm') ;
                ppts_sm = h5read(tubi.fileName.apdv, '/ppts_sm') ;
                assert(size(apts_sm, 1) == length(tubi.xp.fileMeta.timePoints))
            catch
                disp(['Warning: did not find APpointsSm on disk in ' tubi.fileName.apdv ', recomputing...'])
                try
                    opts = load(tubi.fileName.apdvOptions) ;
                catch
                    opts = struct() ;
                    opts.apdvOpts = struct() ;
                end
                [apts_sm, ppts_sm] = tubi.computeAPDpoints(opts.apdvOpts) ;
            end
        end
        
        function [rot, trans] = getRotTrans(tubi)
            % Load the translation to put anterior to origin and AP axis
            % along x axis 
            if ~isempty(tubi.APDV.trans)
                % return from self
                trans = tubi.APDV.trans ;
            else
                % load from disk
                trans = importdata(tubi.fileName.trans) ;
                tubi.APDV.trans = trans ;
            end
            % Load the rotation from XYZ to APDV coordinates
            if ~isempty(tubi.APDV.rot)
                rot = tubi.APDV.rot ;
            else
                % Load the rotation matrix
                rot = importdata(tubi.fileName.rot) ;
                tubi.APDV.rot = rot ;
            end
        end
        
        function [xyzlim_raw, xyzlim_pix, xyzlim_um, xyzlim_um_buff] = ...
                getXYZLims(tubi)
            %[raw, pix, um, um_buff] = GETXYZLIMS(tubi)
            % Grab each xyzlim from self, otherwise load from disk
            % full resolution pix
            if ~isempty(tubi.plotting.xyzlim_raw)
                xyzlim_raw = tubi.plotting.xyzlim_raw ;
            else
                try
                    xyzlim_raw = dlmread(tubi.fileName.xyzlim_raw, ',', 1, 0) ; 
                    tubi.plotting.xyzlim_raw = xyzlim_raw ;
                catch
                    [tubi.plotting.xyzlim_raw, tubi.plotting.xyzlim_pix, ...
                        tubi.plotting.xyzlim_um, ...
                        tubi.plotting.xyzlim_um_buff] = ...
                        tubi.measureXYZLims() ;
                    xyzlim_raw = tubi.plotting.xyzlim_raw ;
                end
            end
            % rotated scaled in full resolution pix
            if ~isempty(tubi.plotting.xyzlim_pix)
                xyzlim_pix = tubi.plotting.xyzlim_pix ;
            else
                try
                    xyzlim_pix = dlmread(tubi.fileName.xyzlim_pix, ',', 1, 0) ; 
                    tubi.plotting.xyzlim_pix = xyzlim_pix ;
                catch
                    [~, tubi.plotting.xyzlim_pix, ...
                        tubi.plotting.xyzlim_um, ...
                        tubi.plotting.xyzlim_um_buff] = ...
                        tubi.measureXYZLims() ;
                    xyzlim_pix = tubi.plotting.xyzlim_pix ;
                end
            end
            % rotated scaled APDV in micron
            if ~isempty(tubi.plotting.xyzlim_um)
                xyzlim_um = tubi.plotting.xyzlim_um ;
            else
                try
                    xyzlim_um = dlmread(tubi.fileName.xyzlim_um, ',', 1, 0) ;
                    tubi.plotting.xyzlim_um = xyzlim_um ;
                catch
                    [~, ~, tubi.plotting.xyzlim_um, ...
                        tubi.plotting.xyzlim_um_buff] = ...
                        tubi.measureXYZLims() ;
                    xyzlim_um = tubi.plotting.xyzlim_um ;
                end
            end
            % rotated scaled APDV in micron, with padding
            if ~isempty(tubi.plotting.xyzlim_um_buff)
                xyzlim_um_buff = tubi.plotting.xyzlim_um_buff ;
            else
                try
                    xyzlim_um_buff = dlmread(tubi.fileName.xyzlim_um_buff, ',', 1, 0) ;
                    tubi.plotting.xyzlim_um_buff = xyzlim_um_buff ;
                catch
                    [~, ~, ~, tubi.plotting.xyzlim_um_buff] = ...
                        tubi.measureXYZLims() ;
                    xyzlim_um_buff = tubi.plotting.xyzlim_um_buff ;
                end
            end
        end
        
        function data = loadBioFormats(tubi, fullFileName)
            % Load a tiff from disk using bioFormats importer
            r = bfGetReader(fullFileName);
            r.setSeries(tubi.xp.fileMeta.series-1);
            nChannelsUsed = numel(tubi.xp.expMeta.channelsUsed);
            if tubi.xp.fileMeta.swapZT == 0
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];
            else
                stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeT(), r.getSizeZ()];
            end
            debugMsg(2, ['stack size (xyzt) ' num2str(stackSize) '\n']);

            xSize = stackSize(1);
            ySize = stackSize(2);
            zSize = stackSize(3);
            
            % number of channels
            nTimePts = stackSize(4);
            
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
            for i = 1:r.getImageCount()

                ZCTidx = r.getZCTCoords(i-1) + 1;
                
                % in the fused embryo data coming out of the python script,
                % Z and T are swaped. In general this isn't the case, thus
                % introduce a file metaField swapZT
                if tubi.xp.fileMeta.swapZT == 0
                    zidx = ZCTidx(1);
                    tidx = ZCTidx(3);
                else 
                    zidx = ZCTidx(3);
                    tidx = ZCTidx(1);
                end
                cidx = ZCTidx(2);

                % see above: if there is only one timepoint all the planes
                % should be read, if there are multiple timepoints, only
                % the correct time should be read
                if nTimePts == 1 || (nTimePts > 1 && this.currentTime == tidx-1)
                    
                    debugMsg(1,'.');
                    if rem(i,80) == 0
                        debugMsg(1,'\n');
                    end

                    dataCidx = find(tubi.xp.expMeta.channelsUsed == cidx);
                    if ~isempty(dataCidx)
                        data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
                    else
                        disp('skipping channel and z plane')
                    end
                end
            end
        end
        
        function setDataLimits(tubi, tp, adjustlow_pctile, adjusthigh_pctile)
            % Use timepoint (tp) to obtain hard values for intensity limits
            % so that data is rescaled to fixed limits instead of
            % percentile. This is useful to avoid flickering of overall
            % intensity in data in which a few voxels vary a lot in
            % intensity.
            tubi.xp.loadTime(tp);
            tubi.xp.rescaleStackToUnitAspect();
            IV = tubi.xp.stack.image.apply() ;
            try
                assert(adjusthigh_pctile > 0 && adjustlow_pctile < 100)
                assert(adjusthigh_pctile > adjustlow_pctile)
            catch
                error('adjustment values must be 0<=val<=100 and increasing')
            end
            adjustlow = prctile(IV{1}(:), adjustlow_pctile) ;
            adjusthigh = prctile(IV{1}(:), adjusthigh_pctile) ;
            tubi.data.adjustlow = adjustlow ;
            tubi.data.adjusthigh = adjusthigh ;
        end
        
        IV = getCurrentData(tubi, adjustIV, varargin)
        
        function IV = adjustIV(tubi, IV, adjustlow, adjusthigh, forceValues)
            % Parameters
            % ----------
            % tubi : QuapSlap class instance
            % IV : 3D numeric (optional, loads currentData if empty)  
            %   data to adjust
            % adjustlow : numeric or #channels x 1 numeric
            %   minimum intensity or percent intensity to adjust data
            % adjusthigh : numeric or #channels x 1 numeric
            %   maximum intensity or percent intensity to adjust data
            % forceValues : enforce the adjustlow/high to be intensity
            % values, not percentile, even if < 100. 
            if nargin > 2 
                adjustlow = tubi.data.adjustlow ;
                adjusthigh = tubi.data.adjusthigh ;
            end
            if nargin < 2 
                if ~isempty(tubi.currentData.IV) 
                    IV = tubi.currentData.IV ;
                end
            end
            if nargin < 5
                forceValues = false ;
            end
            
            % If only one value of intensity limit is supplied, duplicate 
            % for each channel of IV            
            if numel(adjustlow) == 1 && length(IV) > 1
                adjustlow = adjustlow *  ones(size(IV)) ;
            end
            if numel(adjusthigh) == 1 && length(IV) > 1
                adjusthigh = adjusthigh *  ones(size(IV)) ;
            end

            % custom image intensity adjustment
            if all(adjustlow == 0) && all(adjusthigh == 0) && ~forceValues
                disp('Using default limits for imadjustn')
                for ii = 1:length(IV)
                    IV{ii} = imadjustn(IV{ii});
                end
            elseif all(adjustlow < 100) && all(adjusthigh < 100) && ~forceValues
                disp('Taking custom limits for imadjustn as prctile')
                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(prctile( IVii(:) , adjustlow(ii) )) ;
                    vhi = double(prctile( IVii(:) , adjusthigh(ii))) ;
                    disp(['  raw --> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    
                    % Now as fraction for imadjustn
                    vlo = vlo / double(max(IVii(:))) ;
                    vhi = vhi / double(max(IVii(:))) ;
                    disp(['  frac--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    IV{ii} = imadjustn(IVii, [double(vlo); double(vhi)]) ;
                end
            else
                % adjusthigh is > 100 or forceValues, so interpret as an intensity value
                disp('Taking custom limits for imadjustn as direct intensity limit values')

                for ii = 1:length(IV)
                    IVii = IV{ii} ;
                    vlo = double(adjustlow(ii)) ;
                    vhi = double(adjusthigh(ii)) ;
                    disp(['--> ', num2str(vlo), ', ', num2str(vhi), ...
                        ' for ', num2str(adjustlow(ii)), '/', num2str(adjusthigh(ii))])
                    tmp = (double(IVii) - vlo) / (vhi - vlo) ;
                    tmp(tmp > (vhi - vlo)) = 1.0 ;
                    IV{ii} = uint16(2^16 * tmp) ;
                    % cast(tmp, class(IVii)) ;  
                end
            end
            if nargout > 0
                disp('Attributing to self.currentData.IV')
                tubi.currentData.IV = IV ;
                tubi.currentData.adjustlow = adjustlow ;
                tubi.currentData.adjustlow = adjusthigh ;
            else
                disp('WARNING: returning IV instead of attributing to self')
            end
        end
        
        % Get velocity
        function getCurrentVelocity(tubi, varargin)
            if isempty(tubi.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            if isempty(varargin) 
                do_all = true ;
            else
                do_all = false ;
            end
            
            no_piv3d = isempty(fieldnames(tubi.currentVelocity.piv3d)) ;
            if (do_all || contains(varargin, 'piv3d')) && no_piv3d
                % Load 3D data for piv results
                piv3dfn = tubi.fullFileBase.piv3d ;
                load(sprintf(piv3dfn, tubi.currentTime), 'piv3dstruct') ;
                tubi.currentVelocity.piv3d = piv3dstruct ;
            end
        end
        
        % APDV methods
        [apts_sm, ppts_sm] = computeAPDpoints(tubi, opts)
        function ars = xyz2APDV(tubi, a)
            %ars = xyz2APDV(tubi, a)
            %   Transform 3d coords from XYZ data space to APDV coord sys
            [rot, trans] = tubi.getRotTrans() ;
            ars = ((rot * a')' + trans) * tubi.APDV.resolution ;
            if tubi.flipy
                ars(:, 2) = - ars(:, 2) ;
            end
        end
        
        function axyz = APDV2xyz(tubi, a)
            %ars = xyz2APDV(tubi, a)
            %   Transform 3d coords from APDV coord sys to XYZ data space
            [rot, trans] = tubi.getRotTrans() ;
            if tubi.flipy
                a(:, 2) = - a(:, 2) ;
            end
            invRot = tubi.invertRotation(rot) ;
            preRot = a / tubi.APDV.resolution - trans ; 
            axyz = (invRot * preRot')' ;
            % Note: ars = ((rot * axyz')' + trans) * tubi.APDV.resolution ;
        end
        
        function daxyz = APDV2dxyz(tubi, a)
            %ars = xyz2APDV(tubi, a)
            %   Transform 3d vectors from APDV coord sys to XYZ data space
            [rot, trans] = tubi.getRotTrans() ;
            if tubi.flipy
                a(:, 2) = - a(:, 2) ;
            end
            invRot = tubi.invertRotation(rot) ;
            preRot = a / tubi.APDV.resolution ; 
            daxyz = (invRot * preRot')' ;
            % Note: ars = ((rot * axyz')' + trans) * tubi.APDV.resolution ;
        end
        
        function dars = dx2APDV(tubi, da)
            %dars = dx2APDV(tubi, da)
            %   Transform 3d difference vector from XYZ data space to APDV 
            %   coord sys
            [rot, ~] = tubi.getRotTrans() ;
            dars = ((rot * da')') * tubi.APDV.resolution ;
            if tubi.flipy
                dars(:, 2) = - dars(:, 2) ;
            end
        end
        function setAPDVOptions(tubi, apdvOpts)
            tubi.apdvOptions = apdvOpts ;
        end
        function apdvOptions = loadAPDVOptions(tubi)
            load(tubi.fileName.apdvOptions, 'apdvOptions')
            tubi.apdvOptions = apdvOptions ;
        end     
        function apdvOptions = saveAPDVOptions(tubi, apdvOptions)
            if nargin < 2
                apdvOptions = tubi.APDVs.apdvOptions ;
            end
            save(tubi.fileName.apdvOptions, 'apdvOptions')
        end
        [rot, trans, xyzlim_raw, xyzlim, xyzlim_um, xyzlim_um_buff] = ...
            alignMeshesAPDV(tubi, alignAPDVOpts) 
        
        % Plot aligned meshes for pretty presentation 
        % -- for ex, for optogenetic figures
        plotAlignedMeshesPretty(tubi, options)

        % Load raw mesh or alignedMesh (rotated & scaled to APDV)
        function rawMesh = loadCurrentRawMesh(tubi)
            meshfn = sprintf(tubi.fullFileBase.mesh, tubi.currentTime) ;
            rawMesh = read_ply_mod(meshfn) ;
            tubi.currentMesh.rawMesh = rawMesh ;
        end
        function rawMesh = getCurrentRawMesh(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.rawMesh)
                tubi.loadCurrentRawMesh() ;
            end
            if nargout > 0
                rawMesh = tubi.currentMesh.rawMesh ;
            end
        end
        function alignedMesh = loadCurrentAlignedMesh(tubi)
            meshfn = sprintf(tubi.fullFileBase.alignedMesh, tubi.currentTime) ;
            alignedMesh = read_ply_mod(meshfn) ;
            tubi.currentMesh.alignedMesh = alignedMesh ;
        end
        function amesh = getCurrentAlignedMesh(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.spcutMeshSm)
                tubi.loadCurrentAlignedMesh() ;
            end
            if nargout > 0
                amesh = tubi.currentMesh.alignedMesh ;
            end
        end
        
        % Masked Data
        generateMaskedData(tubi, options)
        alignMaskedDataAPDV(tubi)
        
        % TexturePatch and demo mesh visualization
        % tubi.visualizeMeshEvolution() gives evolution sequence with
        % orthogonalviews of data planes -- aka morphsnakes demo
        visualizeMeshEvolution(tubi, options)
        plotSeriesOnSurfaceTexturePatch(tubi, overwrite, metadat, ...
                                        TexturePatchOptions)
        
        % Surface Area and Volume over time
        measureSurfaceAreaVolume(tubi, options)
        
        % Centerlines & cylinderMesh
        extractCenterlineSeries(tubi, cntrlineOpts)
        function setEndcapOptions(tubi, endcapOpts)
            tubi.endcapOptions = endcapOpts ;
        end        
        function loadEndcapOptions(tubi)
            tmp = load(tubi.fileName.endcapOptions, 'endcapOptions');
            tubi.endcapOptions = tmp.endcapOptions ;
        end        
        function saveEndcapOptions(tubi)
            endcapOptions = tubi.endcapOptions ;
            save(tubi.fileName.endcapOptions, 'endcapOptions')
        end
        sliceMeshEndcaps(tubi, endcapOpts, methodOpts)
        generateCleanCntrlines(tubi, idOptions)
        function getCleanCntrlines(tubi)
            if isempty(tubi.cleanCntrlines)
                try
                    tmp = load(tubi.fileName.cleanCntrlines, 'cntrlines') ;
                    tubi.cleanCntrlines = tmp.cntrlines ;
                    disp('Loaded clean centerlines from disk')
                catch
                    disp('No clean centerlines on disk, generating...')
                    tubi.cleanCntrlines = tubi.generateCleanCntrlines() ;
                end
            end
        end
        function loadCurrentCylinderMeshlean(tubi)
            cylmeshfn = ...
                sprintf( tubi.fullFileBase.cylinderMesh, tubi.currentTime ) ;
            tubi.currentMesh.cylinderMesh = read_ply_mod( cylmeshfn );
        end
        function loadCurrentCylinderMeshClean(tubi)
            cylmeshfn = ...
                sprintf( tubi.fullFileBase.cylinderMeshClean, tubi.currentTime ) ;
            disp(['Loading cylinderMeshClean ' cylmeshfn])
            tubi.currentMesh.cylinderMeshClean = read_ply_mod( cylmeshfn );
        end
        
        % RawRicci meshes -- can be used to form cutMeshes
        [rawRicciMesh, rawRicciMu] = ...
            generateRawRicciMeshTimePoint(tubi, tp, options)
        
        % cutMesh
        generateCurrentCutMesh(tubi, options)
        plotCutPath(tubi, cutMesh, cutPath)
        function loadCurrentCutMesh(tubi)
            if isempty(tubi.currentTime)
                error('No currentTime set. Use QuapSlap.setTime()')
            end
            cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tubi.currentTime) ;
            cutPfn = sprintf(tubi.fullFileBase.cutPath, tubi.currentTime) ;
            tmp = load(cutMeshfn, 'cutMesh') ;
            tmp.cutMesh.v = tmp.cutMesh.v + tmp.cutMesh.vn * tubi.normalShift ;
            tubi.currentMesh.cutMesh = tmp.cutMesh ;
            try
                tubi.currentMesh.cutPath = dlmread(cutPfn, ',', 1, 0) ;
            catch
                debugMsg(1, 'Could not load cutPath, cutMesh is limited\n')
                % Wait, isn't cutP a field of cutMesh?
                tmp.cutMesh.cutP
                error('check this here --> is cutP a field?')
            end
        end
        
        % spcutMesh
        generateCurrentSPCutMesh(tubi, cutMesh, overwrite)
        function spcutMesh = getCurrentSPCutMesh(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.spcutMesh)
                tubi.loadCurrentSPCutMesh() ;
            end
            if nargout > 0
                spcutMesh = tubi.currentMesh.spcutMesh ;
            end
        end
        function spcutMesh = loadCurrentSPCutMesh(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            spcutMeshfn = sprintf(tubi.fullFileBase.spcutMesh, tubi.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMesh') ;
            tubi.currentMesh.spcutMesh = tmp.spcutMesh ;
            if nargout > 0
                spcutMesh = tubi.currentMesh.spcutMesh ;
            end
        end        
        function spcutMeshSm = getCurrentSPCutMeshSm(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.spcutMeshSm)
                tubi.loadCurrentSPCutMeshSm() ;
            end
            if nargout > 0
                spcutMeshSm = tubi.currentMesh.spcutMeshSm ;
            end
        end
        function loadCurrentSPCutMeshSm(tubi)
            spcutMeshfn = sprintf(tubi.fullFileBase.spcutMeshSm, tubi.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSm') ;
            tubi.currentMesh.spcutMeshSm = tmp.spcutMeshSm ;
        end
        function spcutMeshSmRS = getCurrentSPCutMeshSmRS(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.spcutMeshSmRS)
                tubi.loadCurrentSPCutMeshSmRS() ;
            end
            if nargout > 0
                spcutMeshSmRS = tubi.currentMesh.spcutMeshSmRS ;
            end
        end
        function mesh = loadCurrentSPCutMeshSmRS(tubi)
            spcutMeshfn = sprintf(tubi.fullFileBase.spcutMeshSmRS, tubi.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSmRS') ;
            tubi.currentMesh.spcutMeshSmRS = tmp.spcutMeshSmRS ;
            if nargout > 0
                mesh = tubi.currentMesh.spcutMeshSmRS ;
            end
        end
        function spcutMeshSmRSC = getCurrentSPCutMeshSmRSC(tubi)
            if isempty(tubi.currentTime)
                error('First set currentTime')
            end
            if isempty(tubi.currentMesh.spcutMeshSmRSC)
                tubi.loadCurrentSPCutMeshSmRSC() ;
            end
            if nargout > 0
                spcutMeshSmRSC = tubi.currentMesh.spcutMeshSmRSC ;
            end
        end
        function loadCurrentSPCutMeshSmRSC(tubi)
            spcutMeshfn = sprintf(tubi.fullFileBase.spcutMeshSmRSC, tubi.currentTime) ;
            tmp = load(spcutMeshfn, 'spcutMeshSmRSC') ;
            tubi.currentMesh.spcutMeshSmRSC = tmp.spcutMeshSmRSC ;
        end
        
        
        function [xyzrs, fieldfaces, tri, baryc] = uv2APDV(tubi, uv, coordSys, umax, vmax)
            % Convert (u,v) pullback coordinates to xyz coordinates in the
            % APDV frame (rotated and scaled frame aligned with AP axis
            % along x and DV axis along z)
            if nargin < 4
                umax = 1 ;
            end
            if nargin < 5 
                vmax = 1 ;
            end
            if isempty(coordSys)
                coordSys = 'spsmrs' ;
            end
            if strcmpi(coordSys, 'spsmrs') 
                mesh = tubi.getCurrentSPCutMeshSmRS() ;
                mesh.u(:, 1) = mesh.u(:, 1) * umax / max(mesh.u(:, 1)) ;
            else
                error('handle coordSys here')
            end
            
            [xyzrs, fieldfaces, tri, baryc] = interpolate2Dpts_3Dmesh(mesh.f, mesh.u,...
                                        mesh.v, uv) ;
            % If the loaded coordinates are not already in APDV, but 
            % instead in the data coordSys, convert to APDV
            if ~contains(coordSys, 'rs')
                xyzrs = xyz2APDV(xyz) ;
            end                       
        end
        
        % Radii in a coordsys
        function [radii, radius_cline] = getRadii(tubi, options)
            % Load radii from disk in specified coordinate system
            %
            % Parameters
            % ----------
            % options : struct with fields
            %   coordSys: str specifier for how to return radii
            %       'spcutMeshSmRSC' --> return nU x (nV-1) array for radii
            %       of vertices in spcutMeshSmRSC embedding
            %       <add other options here>
            if nargin < 2
                options = struct() ;
            end
            if isfield(options, 'coordSys')
                coordSys = options.coordSys ;
            else
                if length(tubi.xp.fileMeta.timePoints) > 1
                    coordSys = 'spcutMeshSmRSC' ;
                else
                    coordSys = 'spcutMesh' ;
                end
            end
            if strcmpi(coordSys, 'spcutMeshSmRSC')
                tubi.loadCurrentSPCutMeshSmRSC() ;
                radii = tubi.currentMesh.spcutMeshSmRSC.radius_um ;
                radius_cline = mean(tubi.currentMesh.spcutMeshSmRSC.radius_um, 2) ;
            elseif strcmpi(coordSys, 'spcutMesh')
                tubi.loadCurrentSPCutMesh() 
                radii = tubi.currentMesh.spcutMesh.radii_from_mean_uniform_rs ;
                radius_cline = mean(tubi.currentMesh.spcutMesh.radii_from_mean_uniform_rs, 2) ;
            else
                error(['Have not coded for this coord sys yet:' coordSys])
            end
        end
        
        % t0_for_phi0 (uvprime cutMesh)
        function mesh = getCurrentUVPrimeCutMesh(tubi)
            if isempty(tubi.currentMesh.uvpcutMesh)
                tubi.loadCurrentSPCutMeshSm() ;
            end
            mesh = tubi.currentMesh.uvpcutMesh ;
        end
        function loadCurrentUVPrimeCutMesh(tubi)
            uvpcutMeshfn = sprintf(tubi.fullFileBase.uvpcutMesh, tubi.currentTime) ;
            tmp = load(uvpcutMeshfn, 'uvpcutMesh') ;
            tubi.currentMesh.uvpcutMesh = tmp.uvpcutMesh ;
        end
        measureUVPrimePathlines(tubi, options)
        % Note: measureBeltramiCoefficient() allows uvprime
        % coordSys.
        
        % Ricci flow for (r=log(rho), phi) coordinate system
        function mesh = getCurrentRicciMesh(tubi)
            if isempty(tubi.currentMesh.ricciMesh)
                try
                    tubi.loadCurrentRicciMesh() ;
                    mesh = tubi.currentMesh.ricciMesh ;
                catch
                    disp('Ricci mesh not on disk! Generating...')
                    mesh = tubi.generateRicciMeshTimePoint(tubi.currentTime) ;
                end
            else
                mesh = tubi.currentMesh.ricciMesh ;
            end
        end
        function ricciMesh = loadCurrentRicciMesh(tubi, maxIter)
            if nargin < 2
                maxIter = 200 ;
            end
            ricciMeshfn = sprintf(tubi.fullFileBase.ricciMesh, ...
                maxIter, tubi.currentTime) ;
            load(ricciMeshfn, 'ricciMesh') ;
            tubi.currentMesh.ricciMesh = ricciMesh ;
        end
        measureRPhiPathlines(tubi, options)
        % Note: measureBeltramiCoefficient() allows ricci coordSys
        measureBeltramiCoefficient(tubi, options)
        function beltrami = getubieltramiCoefficient(tubi, options)
            if nargin < 2
                options = struct() ;
            end
            if isempty(tubi.pathlines.beltrami)
                tubi.pathlines.beltrami = tubi.loadBeltramiCoefficient(options) ;
            end
            if nargout > 0 
                beltrami = tubi.pathlines.beltrami ;
            end
        end
        function beltrami = loadBeltramiCoefficient(tubi, options)    
            if isfield(options, 't0Pathlines')
                t0Pathlines = options.t0Pathlines ;
            else
                t0Pathlines = tubi.t0set() ;
            end
            fn = sprintf(tubi.fileName.pathlines.quasiconformal, t0Pathlines) ;
            tubi.pathlines.beltrami = load(fn)  ;
            if nargout > 0 
                beltrami = tubi.pathlines.beltrami ;
            end
        end
        measureBeltramiCoefficientPullbackToPullback(tubi, options)
        
        % Radial indentation for pathlines
        function indentation = getPathlineIndentation(tubi, options)
            if nargin < 2
                options = struct() ;
            end
            indentation = measurePathlineIndentation(tubi, options) ;
        end
        function indentation = measurePathlineIndentation(tubi, options)
            overwrite = false ;
            t0p = tubi.t0set() ;
            if isfield(options, 'overwrite')
                overwrite = options.overwrite ;
            end
            if isfield(options, 't0Pathline')
                t0p = options.t0Pathline ;
            end
            indentFn = sprintf(tubi.fileName.pathlines.indentation, t0p) ;
            if ~exist(indentFn, 'file') || overwrite 
                radFn = sprintf(tubi.fileName.pathlines.radius, t0p) ;
                if ~exist(radFn, 'file')
                    disp(['pathline radii not on disk: ' radFn])
                    tubi.measurePullbackPathlines(options) ;
                end
                
                load(radFn, 'vRadiusPathlines')
                rad = vRadiusPathlines.radii ;
                nU = size(rad, 2) ;
                indentation = 0 * rad ;
                rad0 = rad(vRadiusPathlines.tIdx0, :, :) ;
                for tidx = 1:length(tubi.xp.fileMeta.timePoints)
                    indentation(tidx, :, :) = -(rad(tidx, :, :) - rad0) ./ rad0 ;
                end
                save(indentFn, 'indentation')
                
                % Plot the indentation as a kymograph
                close all
                figfn = fullfile(sprintf(tubi.dir.pathlines.data, ...
                    t0p), 'indentation_kymograph.png') ;
                set(gcf, 'visible', 'off')
                indentAP = mean(indentation, 3) ;
                uspace = linspace(0, 1, nU) ;
                imagesc(uspace, tubi.xp.fileMeta.timePoints, indentAP)
                xlabel('ap position, $u''/L$', 'interpreter', 'latex')
                ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
                caxis([-max(abs(indentAP(:))), max(abs(indentAP(:)))])
                colormap blueblackred
                cb = colorbar() ;
                ylabel(cb, 'indentation $\delta r/r_0$', 'interpreter', 'latex')
                saveas(gcf, figfn)
                
                % Plot in 3d
                % load reference mesh and pathline vertices in 3d
                load(sprintf(tubi.fileName.pathlines.refMesh, t0p), ...
                    'refMesh') ;
                load(sprintf(tubi.fileName.pathlines.v3d, t0p), 'v3dPathlines') ;
                indentDir = sprintf(tubi.dir.pathlines.indentation, t0p) ;
                if ~exist(indentDir, 'dir')
                    mkdir(indentDir) 
                end
                [~,~,~,xyzlim] = tubi.getXYZLims() ;
                for tidx = 1:size(rad, 1)
                    tp = tubi.xp.fileMeta.timePoints(tidx) ;
                    fn = fullfile(indentDir, 'indentation_%06d.png') ;
                    if ~exist(fn, 'file') || overwrite
                        close all 
                        fig = figure('visible', 'off') ;
                        opts = struct() ;
                        opts.fig = fig ;
                        opts.ax = gca ;
                        xx = v3dPathlines.vXrs(tidx, :) ;
                        yy = v3dPathlines.vYrs(tidx, :) ;
                        zz = v3dPathlines.vZrs(tidx, :) ;
                        v3d = [ xx(:), yy(:), zz(:) ] ;
                        indent = indentation(tidx,:,:) ;
                        opts.sscale = 0.5 ;
                        opts.axisOff = false ;
                        opts.label = 'constriction, $\delta r/r_0$' ;
                        opts.ax_position = [0.1141, 0.1100, 0.6803, 0.8150] ;
                        scalarFieldOnSurface(refMesh.f, v3d, indent(:), opts) ;
                        view(0, 0)
                        axis equal
                        axis off
                        xlim(xyzlim(1, :))
                        ylim(xyzlim(2, :))
                        zlim(xyzlim(3, :))
                        sgtitle(['constriction, $t=$', ...
                            sprintf('%03d', (tp-t0p)*tubi.timeInterval), ...
                            ' ', tubi.timeUnits ], 'interpreter', 'latex')
                        saveas(gcf, sprintf(fn, tp)) ;
                        close all
                    end
                end
                
            else
                load(indentFn, 'indentation')
            end
        end
        
        % Radial indentation for UVprime pathlines
        function indentation = measureUVPrimePathlineIndentation(tubi, options)
            overwrite = false ;
            t0p = tubi.t0set() ;
            if isfield(options, 'overwrite')
                overwrite = options.overwrite ;
            end
            if isfield(options, 't0Pathline')
                t0p = options.t0Pathline ;
            end
            indentFn = sprintf(tubi.fileName.pathlines_uvprime.indentation, t0p) ;
            if ~exist(indentFn, 'file') || overwrite
                radFn = sprintf(tubi.fileName.pathlines_uvprime.radius, t0p) ;
                if ~exist(radFn, 'file')
                    disp(['pathline radii not on disk: ' radFn])
                    tubi.measureUVPrimePathlines(options) ;
                end
                
                load(radFn, 'vRadiusPathlines')
                rad = vRadiusPathlines.radii ;
                nU = size(rad, 2) ;
                indentation = 0 * rad ;
                rad0 = rad(vRadiusPathlines.tIdx0, :, :) ;
                for tidx = 1:length(tubi.xp.fileMeta.timePoints)
                    indentation(tidx, :, :) = -(rad(tidx, :, :) - rad0) ./ rad0 ;
                end
                save(indentFn, 'indentation')
                
                % Plot the indentation as a kymograph
                close all
                figfn = fullfile(sprintf(tubi.dir.pathlines_uvprime.data, ...
                    t0p), 'indentation_kymograph.png') ;
                set(gcf, 'visible', 'off')
                indentAP = mean(indentation, 3) ;
                uspace = linspace(0, 1, nU) ;
                imagesc(uspace, tubi.xp.fileMeta.timePoints, indentAP)
                xlabel('ap position, $u''/L$', 'interpreter', 'latex')
                ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
                caxis([-max(abs(indentAP(:))), max(abs(indentAP(:)))])
                colormap blueblackred
                cb = colorbar() ;
                ylabel(cb, 'indentation $\delta r/r_0$', 'interpreter', 'latex')
                saveas(gcf, figfn)
                
                % Plot in 3d
                % load reference mesh and pathline vertices in 3d
                load(sprintf(tubi.fileName.pathlines_uvprime.refMesh, t0p), ...
                    'refMesh') ;
                load(sprintf(tubi.fileName.pathlines_uvprime.v3d, t0p), 'v3dPathlines') ;
                indentDir = sprintf(tubi.dir.pathlines_uvprime.indentation, t0p) ;
                if ~exist(indentDir, 'dir')
                    mkdir(indentDir) 
                end
                [~,~,~,xyzlim] = tubi.getXYZLims() ;
                for tidx = 1:size(rad, 1)
                    tp = tubi.xp.fileMeta.timePoints(tidx) ;
                    fn = fullfile(indentDir, 'indentation_%06d.png') ;
                    if ~exist(fn, 'file') || overwrite
                        close all 
                        fig = figure('visible', 'off') ;
                        opts = struct() ;
                        opts.fig = fig ;
                        opts.ax = gca ;
                        xx = v3dPathlines.vXrs(tidx, :) ;
                        yy = v3dPathlines.vYrs(tidx, :) ;
                        zz = v3dPathlines.vZrs(tidx, :) ;
                        v3d = [ xx(:), yy(:), zz(:) ] ;
                        indent = indentation(tidx,:,:) ;
                        opts.sscale = 0.5 ;
                        opts.axisOff = false ;
                        opts.label = 'constriction, $\delta r/r_0$' ;
                        opts.ax_position = [0.1141, 0.1100, 0.6803, 0.8150] ;
                        scalarFieldOnSurface(refMesh.f, v3d, indent(:), opts) ;
                        view(0, 0)
                        axis equal
                        axis off
                        xlim(xyzlim(1, :))
                        ylim(xyzlim(2, :))
                        zlim(xyzlim(3, :))
                        sgtitle(['constriction, $t=$', sprintf('%03d', tp), ...
                            ' ', tubi.timeUnits ], 'interpreter', 'latex')
                        saveas(gcf, sprintf(fn, tp)) ;
                        close all
                    end
                end
                
            else
                load(indentFn, 'indentation')
            end
        end
        
        
        % Pullbacks
        generateCurrentPullbacks(tubi, cutMesh, spcutMesh, spcutMeshSm, pbOptions)
        
        function doubleCoverPullbackImages(tubi, options)
            % options : struct with fields
            %   coordsys : ('sp', 'uv', 'up')
            %       coordinate system to make double cover 
            %   overwrite : bool, default=false
            %       whether to overwrite current images on disk
            %   histeq : bool, default=true
            %       perform histogram equilization during pullback
            %       extension
            %   ntiles : int, default=50 
            %       The number of bins in each dimension for histogram equilization for a
            %       square original image. That is, the extended image will have (a_fixed *
            %       ntiles, 2 * ntiles) bins in (x,y).
            %   a_fixed : float, default=tubi.a_fixed
            %       The aspect ratio of the pullback image: Lx / Ly       
            
            if nargin > 1
                % unpack options
                if isfield(options, 'coordsys')
                    coordsys = options.coordsys ;
                    options = rmfield(options, 'coordsys') ;
                    if strcmp(coordsys, 'sp')
                        imDir = tubi.dir.im_sp ;
                        imDir_e = tubi.dir.im_spe ;
                        fn0 = tubi.fileBase.im_sp ;
                        ofn = tubi.fileBase.im_spe ;
                    elseif strcmp(coordsys, 'spsm') || strcmp(coordsys, 'sp_sm')
                        imDir = tubi.dir.im_sp_sm ;
                        imDir_e = tubi.dir.im_sp_sme ;
                        fn0 = tubi.fileBase.im_sp_sm ;
                        ofn = tubi.fileBase.im_sp_sme ;
                    elseif strcmp(coordsys, 'spsm2') || ...
                            strcmp(coordsys, 'sp_sm2') || ...
                            strcmp(coordsys, 'spsmLUT') || ...
                            strcmp(coordsys, 'sp_smLUT')
                        % equalize the histogram in patches of the image
                        options.histeq = true ;
                        imDir = tubi.dir.im_sp_sm ;
                        imDir_e = tubi.dir.im_sp_smeLUT ;
                        fn0 = tubi.fileBase.im_sp_sm ;
                        ofn = tubi.fileBase.im_sp_smeLUT ;
                    elseif strcmp(coordsys, 'rsm') || strcmp(coordsys, 'r_sm')
                        imDir = tubi.dir.im_r_sm ;
                        imDir_e = tubi.dir.im_r_sme ;
                        fn0 = tubi.fileBase.im_r_sm ;
                        ofn = tubi.fileBase.im_r_sme ;
                    elseif strcmp(coordsys, 'uv')
                        imDir = tubi.dir.im_uv ;
                        imDir_e = tubi.dir.im_uve ;
                        fn0 = tubi.fileBase.im_uv ;
                        ofn = tubi.fileBase.im_uv_sme ;
                    elseif strcmp(coordsys, 'up')
                        imDir = tubi.dir.im_up ;
                        imDir_e = tubi.dir.im_upe ;
                        fn0 = tubi.fileBase.im_up ;
                        ofn = tubi.fileBase.im_up_e ;
                    elseif strcmp(coordsys, 'uvprime')
                        imDir = tubi.dir.im_uvprime ;
                        imDir_e = tubi.dir.im_uvprime_e ;
                        fn0 = tubi.fileBase.im_uvprime ;
                        ofn = tubi.fileBase.im_uvprime_e ;
                    elseif strcmp(coordsys, 'r') || strcmp(coordsys, 'spr') || strcmp(coordsys, 'sphir') || strcmp(coordsys, 'rsp')
                        imDir = tubi.dir.im_r ;
                        imDir_e = tubi.dir.im_re ;
                        fn0 = tubi.fileBase.im_r ;
                        ofn = tubi.fileBase.im_re ;
                    else
                        error(['did not recognize coordsys: ' coordsys])
                    end
                else
                    % Default value of coordsys = 'sp' ;
                    imDir = tubi.dir.im_sp ;
                    imDir_e = tubi.dir.im_spe ;
                    fn0 = tubi.fileBase.im_sp ;
                    ofn = tubi.fileBase.im_spe ;
                end
                
                % Customization of the paths
                if isfield(options, 'subdir')
                    imDir = fullfile(imDir, options.subdir) ;
                end
                if isfield(options, 'fn0')
                    fn0 = options.fn0 ;
                end
                
                % pack options if missing fields
                if ~isfield(options, 'histeq')
                    % equalize the histogram in patches of the image
                    options.histeq = false ;
                end
                if ~isfield(options, 'a_fixed')
                    % Assign the aspect ratio for histogram equilization
                    options.a_fixed = tubi.a_fixed ;
                end
                if ~isfield(options, 'ntiles')
                    % Number of tiles per circumference and per unit ap
                    % length, so that if aspect ratio is two, there will be
                    % 2*ntiles samplings for histogram equilization along
                    % the ap axis
                    options.ntiles = 50 ;
                end
                if ~isfield(options, 'overwrite')
                    options.overwrite = false ;
                end
            else
                % Default options
                % Default value of coordsys = 'sp' ;
                options = struct() ;
                imDir = tubi.dir.im_sp ;
                imDir_e = tubi.dir.im_spe ;
                fn0 = tubi.fileBase.im_sp ;
                ofn = tubi.fileBase.im_spe ;
                options.histeq = false ;
                options.a_fixed = tubi.a_fixed ;
                options.ntiles = ntiles ;
            end
            options.outFnBase = ofn ;
            extendImages(imDir, imDir_e, fn0, tubi.xp.fileMeta.timePoints, options)
            disp(['done ensuring extended tiffs for ' imDir ' in ' imDir_e])
        end
        
        % measure mean hoop centerline, Length, & Writhe of curve
        function [mss, mcline, avgpts] = getCurrentClineDVhoop(tubi)
            clineDVhoopBase = tubi.fullFileBase.clineDVhoop ;
            % load the centerline for this timepoint
            if isempty(tubi.currentTime)
                error('Please first set currentTime')
            end
            fn = sprintf(clineDVhoopBase, tubi.currentTime) ;
            disp(['Loading DVhoop centerline from ' fn])
            load(fn, 'mss', 'mcline', 'avgpts')
            
            % Is this needed here? Yes!
            if tubi.flipy
                mcline(:, 2) = -mcline(:, 2) ;
                avgpts(:, 2) = -avgpts(:, 2) ;
            end
            tubi.currentCline.mss = mss ;
            tubi.currentCline.mcline = mcline ;
            tubi.currentCline.avgpts = avgpts ;
        end

        [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
            measureWrithe(tubi, options)
        function Length_t = measureLength(tubi, options)
            if nargin < 2
                options = struct() ;
            end
            [~,~,~, Length_t, ~] = tubi.measureWrithe(options) ;
        end
                
        % Smooth meshes in time
        [v3dsmM, nsmM] = smoothDynamicSPhiMeshes(tubi, options) ;
        function plotSPCutMeshSm(tubi, options) 
            plotSPCutMeshSmSeriesUtility(tubi, 'spcutMeshSm', options)
        end
        function plotSPCutMeshSmRS(tubi, options) 
            plotSPCutMeshSmSeriesUtility(tubi, 'spcutMeshSmRS', options)
        end
        function plotSPCutMeshSmRSC(tubi, options) 
            plotSPCutMeshSmSeriesUtility(tubi, 'spcutMeshSmRSC', options)
        end
        function [v3dsmM, nsmM] = loadSPCutMeshSm(tubi) 
            timePoints = tubi.xp.fileMeta.timePoints ;
            v3dsmM = zeros(length(timePoints), tubi.nU*tubi.nV, 3);
            nsmM = zeros(length(timePoints), tubi.nU*tubi.nV, 3) ;
            % Load each mesh into v3dsmM and nsmM    
            for qq = 1:length(timePoints)
                load(sprintf(tubi.fullFileBase.spcutMeshSm, ...
                    timePoints(qq)), 'spcutMeshSm') ;
                v3dsmM(qq, :, :) = spcutMeshSm.v ;
                nsmM(qq, :, :) = spcutMeshSm.vn ;
            end
        end
        
        % Mean & Gaussian curvature videos
        measureCurvatures(tubi, options)
        
        % Cell segmentation
        generateCellSegmentation2D(tubi, options)
        processCorrectedCellSegmentation2D(tubi, options)
        generateCellSegmentation3D(tubi, options)
        generateCellSegmentationPathlines3D(tubi, options)
        function seg2d = getCurrentSegmentation2D(tubi, options)
            % Obtain the cell segmentation in 3D pullback space
            if isempty(tubi.currentSegmentation.seg2d)
                try
                    tubi.currentSegmentation.seg2d = load(sprintf(tubi.fullFileBase.segmentation2d, tubi.currentTime)) ;
                catch
                    options.timePoints = [tubi.currentTime] ;
                    tubi.generateCellSegmentation2D(options) ;
                    tubi.currentSegmentation.seg2d = load(sprintf(tubi.fullFileBase.segmentation2d, tubi.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg2d = tubi.currentSegmentation.seg2d ;
            end
        end
        function seg2dCorr = getCurrentSegmentation2DCorrected(tubi, options)
            % Decide on coordinate system for corrected pullback binary 
            coordSys = 'spsme';
            if nargin > 1
                if isfield(options, 'coordSys')
                    coordSys = options.coordSys ;
                end
            end
                
            % Obtain the cell segmentation in 3D pullback space
            if isempty(tubi.currentSegmentation.seg2dCorrected)
                try
                    tubi.currentSegmentation.seg2dCorrected = ...
                        load(sprintf(tubi.fullFileBase.segmentation2dCorrected, coordSys, tubi.currentTime)) ;
                catch
                    options.timePoints = [tubi.currentTime] ;
                    tubi.processCorrectedCellSegmentation2D(options) ;
                    tubi.currentSegmentation.seg2dCorrected = ...
                        load(sprintf(tubi.fullFileBase.segmentation2dCorrected, coordSys, tubi.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg2dCorr = tubi.currentSegmentation.seg2dCorrected ;
                assert(strcmpi(coordSys, seg2dCorr.coordSys))
            end
        end
        function seg3d = getCurrentSegmentation3D(tubi, options)
            % Obtain the cell segmentation in 3D pushforward space
            if isempty(tubi.currentSegmentation.seg3d)
                try
                    tubi.loadCurrentSegmentation3D() ;
                catch
                    options.timePoints = [tubi.currentTime] ;
                    tubi.generateCellSegmentation3D(options) ;
                    tubi.loadCurrentSegmentation3D() ;
                end    
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg3d = tubi.currentSegmentation.seg3d ;
            end
        end
        function seg3dCorr = getCurrentSegmentation3DCorrected(tubi, options)
            % Obtain the cell segmentation in 3D pullback space
            if isempty(tubi.currentSegmentation.seg3dCorrected)
                try
                    tubi.currentSegmentation.seg3dCorrected = ...
                        load(sprintf(tubi.fullFileBase.segmentation3dCorrected, tubi.currentTime)) ;
                catch
                    options.timePoints = [tubi.currentTime] ;
                    options.corrected = true ;
                    tubi.processCorrectedCellSegmentation3D(options) ;
                    tubi.currentSegmentation.seg3dCorrected = ...
                        load(sprintf(tubi.fullFileBase.segmentation3dCorrected, tubi.currentTime)) ;
                end
            end
            % if requested, return segmentation as output
            if nargout > 0
                seg3dCorr = tubi.currentSegmentation.seg3dCorrected ;
            end
        end
        function seg3d = loadCurrentSegmentation3D(tubi) 
            tubi.currentSegmentation.seg3d = ...
                load(sprintf(tubi.fullFileBase.segmentation3d, tubi.currentTime)) ;
            if nargout > 0
                seg3d = tubi.currentSegmentation.seg3d ;
            end
        end
        function seg3dCorr = loadCurrentSegmentation3DCorrected(tubi) 
            tubi.currentSegmentation.seg3dCorrected = ...
                load(sprintf(tubi.fullFileBase.segmentation3dCorrected, tubi.currentTime)) ;
            if nargout > 0
                seg3dCorr = tubi.currentSegmentation.seg3dCorrected ;
            end
        end
        % Note: anisotropy is stored in seg3d.quality
        % no need for: measureCellAnisotropy(tubi, options)
        plotSegmentationStatisticsLobes(tubi, options)
        estimateIntercalationRate(tubi, options)
        
        % density of cells -- nuclei or membrane based
        measureCellDensity(tubi, nuclei_or_membrane, options)
        function loadCurrentCellDensity(tubi)
            if tubi.currentData
                disp('Loading from self')
            else
                disp('Loading from disk')
            end
        end
        plotCellDensity(tubi, options)
        plotCellDensityKymograph(tubi, options)
        
        % spcutMeshStack 
        %   --> similar to spcutMeshSmStack, but no smoothing
        %   --> this is useful for fixed data
        generateSPCutMeshStack(tubi, spcutMeshStackOptions)
        
        % spcutMeshSmStack
        generateSPCutMeshSmStack(tubi, spcutMeshSmStackOptions)
        measureThickness(tubi, thicknessOptions)
        phi0_fit = fitPhiOffsetsViaTexture(tubi, uspace_ds_umax, vspace,...
            phi0_init, phi0TextureOpts)
       
        % uvprime cutMeshSm
        generateUVPrimeCutMeshes(tubi, options)
        
        % ricci Mesh (truly conformal mapped mesh)
        [ricciMesh, ricciMu] = generateRicciMeshTimePoint(tubi, tp, options) 
        
        % spcutMeshSm coordinate system demo
        coordinateSystemDemo(tubi)
        
        % flow measurements
        function getPIV(tubi, options)
            % Load PIV results and store in tubi.piv if not already loaded
            if isempty(fieldnames(tubi.piv.raw)) || isempty(tubi.piv.Lx) ...
                    || isempty(tubi.piv.Lx) 
                % Load raw PIV results
                if nargin > 1
                    tubi.loadPIV(options)
                else
                    tubi.loadPIV() 
                end
            end
            
            if isempty(fieldnames(tubi.piv.smoothed)) 
                % Additionally smooth the piv output by sigma
                if tubi.piv.smoothing_sigma > 0
                    piv = tubi.piv.raw ;
                    disp(['Smoothing piv output with sigma=' num2str(tubi.piv.smoothing_sigma)])
                    for tidx = 1:length(tubi.xp.fileMeta.timePoints)-1
                        velx = piv.u_filtered{tidx} ;
                        vely = piv.v_filtered{tidx} ;
                        piv.u_filtered{tidx} = imgaussfilt(velx, tubi.piv.smoothing_sigma) ;
                        piv.v_filtered{tidx} = imgaussfilt(vely, tubi.piv.smoothing_sigma) ;
                    end
                    disp('done smoothing')
                    tubi.piv.smoothed = piv ;
                end
            end
            
        end
        function piv = loadPIV(tubi, options)
            % Load PIV results from disk and store in tubi.piv
            if ~isempty(fieldnames(tubi.piv.raw)) && isempty(tubi.piv.Lx) ...
                    && isempty(tubi.piv.Lx) 
                disp("WARNING: Overwriting tubi.piv with piv from disk")
            end
            tubi.piv.raw = load(tubi.fileName.pivRaw.raw) ;  
            timePoints = tubi.xp.fileMeta.timePoints ;
            if strcmp(tubi.piv.imCoords, 'sp_sme')
                im0 = imread(sprintf(tubi.fullFileBase.im_sp_sme, ...
                    timePoints(1))) ;
                % for now assume all images are the same size
                tubi.piv.Lx = size(im0, 1) * ones(length(timePoints), 1) ;
                tubi.piv.Ly = size(im0, 2) * ones(length(timePoints), 1) ;
            else
                error(['Unrecognized imCoords: ' tubi.piv.imCoords])
            end
            if nargout > 0
                piv = tubi.piv ;
            end
        end
        measurePIV2d(tubi, options)
        measurePIV3d(tubi, options)
        
        %% Metric
        plotMetric(tubi, options) 
        
        %% Pathlines
        function [p2d, p3d] = samplePullbackPathlines(tubi, XY0, options)
            %[p2d, p3d] = samplePullbackPathlines(tubi, XY0, options)
            % start at XY0, folow flow using barycentric coordinates of PIV
            % pullback pathlines. This is the same as advecting along a
            % fixed Langrangian coordinate. 
            
            % default options
            t0 = tubi.t0set() ;
            preview = false ;
            forceShift = NaN ;
            
            % unpack options
            if nargin < 3
                options = struct() ;
            end
            if isfield(options, 't0')
                t0 = options.t0 ;
            elseif isfield(options, 't0Pathlines')
                t0 = options.t0Pathlines ;
            end
            if isfield(options, 'preview')
                preview = options.preview ;
            end
            
            tidx0 = tubi.xp.tIdx(t0) ;
            timePoints = tubi.xp.fileMeta.timePoints ;
            
            % get image size if we push to 3d or forceShift
            if strcmp(tubi.piv.imCoords, 'sp_sme')
                im = imread(sprintf(tubi.fullFileBase.im_sp_sme, t0)) ;
                forceShift = -size(im, 1) * 0.5 ;
            else
                error('Handle coordSys here')
            end
            
            pathlines = tubi.getPullbackPathlines(t0, 'vertexPathlines', ...
                'vertexPathlines3D') ;
            XX = squeeze(pathlines.vertices.vX(tidx0, :, :)) ;
            YY = squeeze(pathlines.vertices.vY(tidx0, :, :)) ;
            nX = length(unique(XX(:))) ;
            nY = length(unique(YY(:))) ;
            % pivfaces = defineFacesRectilinearGrid(XYgrid, nX, nY) ;
            faces = pathlines.refMesh.f ;
            XY = [XX(:), YY(:)] ;
            cmesh = struct('f', faces, 'u', XY, ...
                'pathPairs', pathlines.refMesh.pathPairs) ;
            tileOpts = struct() ;
            tileOpts.forceShift = forceShift ;
            [tiledFaces, tiledXY] = tileAnnularCutMesh2D(cmesh, ...
                [1, 1], tileOpts) ;
            
            % Note: the following is an exerpt from barycentricMap2d(faces, v2d, vmap, uv)
            % tr0 = triangulation(pivfaces, XY) ;
            % [fieldfaces, baryc0] = pointLocation(tr0, XY0) ; 
            
            p2d = zeros(length(timePoints), size(XY0, 1), 2) ;
            if nargout > 1
                p3d = zeros(length(timePoints), size(XY0, 1), 3) ;
            end
            for tidx = 1:length(timePoints)
                tp = timePoints(tidx) ;
                if mod(tidx, 10) == 0
                    disp(['t=' num2str(tp)])
                end
                % X1 = squeeze(pathlines.piv.XX(tidx, :, :))' ;
                % Y1 = squeeze(pathlines.piv.YY(tidx, :, :))' ;
                X1 = squeeze(pathlines.vertices.vX(tidx, :, :)) ;
                Y1 = squeeze(pathlines.vertices.vY(tidx, :, :)) ;
                XY1 = [X1(:), Y1(:)] ;
                
                
                if nargout > 1
                    % tile current mesh in 2d and 3d
                    v3d1x = squeeze(pathlines.vertices3d.vXrs(tidx, :, :)) ;
                    v3d1y = squeeze(pathlines.vertices3d.vYrs(tidx, :, :)) ;
                    v3d1z = squeeze(pathlines.vertices3d.vZrs(tidx, :, :)) ;
                    
                    cmesh = struct('f', faces, 'u', XY1, ...
                        'v', [v3d1x(:), v3d1y(:), v3d1z(:)], ...
                        'pathPairs', pathlines.refMesh.pathPairs) ;
                    [~, XY1, XYZ1] = tileAnnularCutMesh(cmesh, [1, 1], tileOpts) ;
                    p2d(tidx, :, :) = barycentricMap2d(tiledFaces, tiledXY, XY1, XY0) ;
                    p3d(tidx, :, :) = barycentricMap2d(tiledFaces, tiledXY, XYZ1, XY0) ;
                    assert(~any(isnan(p2d(:))))
                                        
                    if preview && mod(tidx, 10) == 0
                        clf
                        subplot(1, 2, 1)
                        plot(p2d(tidx, :, 1), p2d(tidx, :, 2), '.')
                        axis equal
                        title(['t=' num2str(tp)])
                        subplot(1, 2, 2)
                        plot3(p3d(tidx, :, 1), p3d(tidx, :, 2), ...
                            p3d(tidx, :, 3), '.')
                        axis equal
                        pause(0.1)
                    end
                else
                    % tile current mesh
                    cmesh = struct('f', faces, 'u', XY1, ...
                        'pathPairs', pathlines.refMesh.pathPairs) ;
                    [~, XY1] = tileAnnularCutMesh2D(cmesh, [1, 1]) ;
                    p2d(tidx, :, :) = barycentricMap2d(...
                        tiledFaces, tiledXY, XY1, XY0) ;
                    assert(~any(isnan(p2d(:))))
                end
            end
            
            if nargout > 1
                tubi.clearTime() ;
            end
        end
        
        measurePullbackPathlines(tubi, options)
        function pathlines = getPullbackPathlines(tubi, t0, varargin)
            % Discern if we must load pathlines or if already loaded
            doneLoading = false ;
            if nargin > 1 
                if tubi.pathlines.t0 ~= t0
                    % The timestamp at which pathlines form grid that is 
                    % requested is different than the one that is loaded,
                    % if any are indeed already loaded. Therefore, we load
                    % anew
                    if nargin > 2
                        % pass varargin along to load method
                        tubi.loadPullbackPathlines(t0, varargin)
                    else
                        tubi.loadPullbackPathlines(t0)
                    end
                    doneLoading = true ;
                end
            end
            if ~doneLoading
                % No t0 supplied or t0 is the same as what is stored
                % in tubi.pathlines.t0, if any is already stored (ie if any
                % pathlines are already loaded)
                if isempty(tubi.pathlines.t0)
                    % no pathlines loaded. Load here
                    if nargin > 2
                        tubi.loadPullbackPathlines(t0, varargin)
                    elseif nargin > 1
                        tubi.loadPullbackPathlines(t0)
                    else
                        tubi.loadPullbackPathlines()
                    end
                else
                    % There are pathlines loaded already. Which are
                    % requested here in varargin? First check if varargin 
                    % is empty or not  
                    if nargin > 2            
                        if any(contains(lower(varargin), 'pivpathlines')) || ...
                                any(contains(lower(varargin), 'piv'))
                            if isempty(tubi.pathlines.piv)
                                disp('Loading pivPathlines') 
                                tubi.loadPullbackPathlines(t0, 'pivPathlines')
                            end          
                        end
                        if any(contains(lower(varargin), 'vertexpathlines')) || ...
                                any(contains(lower(varargin), 'vertex')) || ...
                                any(contains(lower(varargin), 'vertices'))
                            if isempty(tubi.pathlines.vertices)
                                disp('Loading vertexPathlines') 
                                tubi.loadPullbackPathlines(t0, 'vertexPathlines')
                            end            
                        end            
                        if any(contains(lower(varargin), 'facepathlines')) || ...
                                any(contains(lower(varargin), 'face')) || ...
                                any(contains(lower(varargin), 'faces'))
                            disp('Loading facePathlines') 
                            if isempty(tubi.pathlines.faces)
                                tubi.loadPullbackPathlines(t0, 'facePathlines')
                            end
                        end
                        if any(contains(lower(varargin), 'vertexpathlines3d')) || ...
                                any(contains(lower(varargin), 'vertex3d')) || ...
                                any(contains(lower(varargin), 'vertices3d'))
                            disp('Loading vertexPathlines3d') 
                            if isempty(tubi.pathlines.vertices3d)
                                tubi.loadPullbackPathlines(t0, 'vertexPathlines3d')
                            end
                        end
                    else
                        % varargin is not supplied, so load all three if
                        % not already loaded
                        % First grab t0
                        if nargin < 2
                            t0 = tubi.t0set() ;
                        end
                        if isempty(tubi.pathlines.piv)
                            disp('Loading pivPathlines') 
                            tubi.loadPullbackPathlines(t0, 'pivPathlines')
                        end            
                        if isempty(tubi.pathlines.vertices)
                            disp('Loading vertexPathlines') 
                            tubi.loadPullbackPathlines(t0, 'vertexPathlines')
                        end            
                        if isempty(tubi.pathlines.faces)
                            disp('Loading facePathlines') 
                            tubi.loadPullbackPathlines(t0, 'facePathlines')
                        end  
                        if isempty(tubi.pathlines.vertices3d)
                            disp('Loading vertexPathlines3d') 
                            tubi.loadPullbackPathlines(t0, 'vertexPathlines3d')
                        end
                    end
                end
            end
            if nargout > 0
                pathlines = tubi.pathlines ;
            end
        end
        function loadPullbackPathlines(tubi, t0, varargin)
            if nargin < 2
                t0 = tubi.t0set() ;
            elseif isempty(t0)
                t0 = tubi.t0set() ;
            else
                try
                    assert(isnumeric(t0))
                catch
                    error('t0 supplied must be numeric')
                end
            end
            % assign t0 as the pathline t0
            tubi.pathlines.t0 = t0 ;
            tmp = load(sprintf(tubi.fileName.pathlines.refMesh, t0), ...
                'refMesh') ;
            tubi.pathlines.refMesh = tmp.refMesh ;
            
            if nargin < 3
                varargin = {'pivPathlines', 'vertexPathlines', ...
                    'facePathlines'} ;
            end
            if any(contains(varargin, 'pivPathlines'))
                load(sprintf(tubi.fileName.pathlines.XY, t0), 'pivPathlines')
                tubi.pathlines.piv = pivPathlines ;
            end
            if any(contains(varargin, 'vertexPathlines'))
                load(sprintf(tubi.fileName.pathlines.vXY, t0), 'vertexPathlines')
                tubi.pathlines.vertices = vertexPathlines ;
            end
            if any(contains(varargin, 'facePathlines'))
                load(sprintf(tubi.fileName.pathlines.fXY, t0), 'facePathlines')
                tubi.pathlines.faces = facePathlines ;
            end
            if any(contains(varargin, 'vertexPathlines3d'))
                tmp = load(sprintf(tubi.fileName.pathlines.v3d, t0)) ;
                tubi.pathlines.vertices3d = tmp.v3dPathlines ;
                tubi.pathlines.vertices3d.smoothing_sigma = tmp.smoothing_sigma ;
            end
        end
        
        %% Velocities -- loading Raw / noAveraging
        function loadVelocityRaw(tubi, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(tubi.fileName.pivRaw.v3d, 'vsmM') ;
                tubi.velocityAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(tubi.fileName.pivRaw.v2dum, 'v2dsmMum') ;
                tubi.velocityAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(tubi.fileName.pivRaw.vn, 'vnsmM') ;
                tubi.velocityAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(tubi.fileName.pivRaw.vf, 'vfsmM') ;
                tubi.velocityAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(tubi.fileName.pivRaw.vv, 'vvsmM') ;
                tubi.velocityAverage.vv = vvsmM ;
            end
        end
        function getVelocityRaw(tubi, varargin)
            % todo: check if all varargin are already loaded
            loadVelocityRaw(tubi, varargin{:})
        end
        
        %% Velocities -- Lagrangian Averaging
        timeAverageVelocities(tubi, samplingResolution, options)
        function loadVelocityAverage(tubi, varargin)
            % Load and pack into struct
            if isempty(varargin)
                varargin = {'v3d', 'v2dum', 'v2d', 'vn', 'vf', 'vv'};
            end
            if any(strcmp(varargin, 'v3d'))
                load(tubi.fileName.pivAvg.v3d, 'vsmM') ;
                tubi.velocityAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(tubi.fileName.pivAvg.v2dum, 'v2dsmMum') ;
                tubi.velocityAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(tubi.fileName.pivAvg.vn, 'vnsmM') ;
                tubi.velocityAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(tubi.fileName.pivAvg.vf, 'vfsmM') ;
                tubi.velocityAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'vv'))
                load(tubi.fileName.pivAvg.vv, 'vvsmM') ;
                tubi.velocityAverage.vv = vvsmM ;
            end
        end
        function getVelocityAverage(tubi, varargin)
            % todo: check if all varargin are already loaded
            loadVelocityAverage(tubi, varargin{:})
        end
        plotTimeAvgVelocities(tubi, options)
        helmholtzHodge(tubi, options)
        measurePathlineVelocities(tubi, options)
        plotPathlineVelocities(tubi, options)
        
        %% Velocities -- simple/surface-Lagrangian averaging
        timeAverageVelocitiesSimple(tubi, samplingResolution, options)
        function loadVelocitySimpleAverage(tubi, varargin)
            % Load and pack into struct
            if any(strcmp(varargin, 'v3d'))
                load(tubi.fileName.pivSimAvg.v3d, 'vsmM') ;
                tubi.velocitySimpleAverage.v3d = vsmM ;
            end
            if any(strcmp(varargin, 'v2dum'))
                load(tubi.fileName.pivSimAvg.v2dum, 'v2dsmMum') ;
                tubi.velocitySimpleAverage.v2dum = v2dsmMum ;
            end
            if any(strcmp(varargin, 'vn'))
                load(tubi.fileName.pivSimAvg.vn, 'vnsmM') ;
                tubi.velocitySimpleAverage.vn = vnsmM ;
            end
            if any(strcmp(varargin, 'vf'))
                load(tubi.fileName.pivSimAvg.vf, 'vfsmM') ;
                tubi.velocitySimpleAverage.vf = vfsmM ;
            end
            if any(strcmp(varargin, 'v2v'))
                load(tubi.fileName.pivSimAvg.vf, 'vvsmM') ;
                tubi.velocitySimpleAverage.vv = vvsmM ;
            end
        end
        function getVelocitySimpleAverage(tubi, varargin)
            % todo: check if all varargin are already loaded
            loadVelocitySimpleAverage(tubi, varargin{:})
        end
        % NOTE: the following have a simple option for averagingStyle
        % plotTimeAvgVelocities(tubi, options)
        % helmholtzHodge(tubi, options)
        
        %% compressible/incompressible flow on evolving surface
        measureMetricKinematics(tubi, options)
        plotMetricKinematics(tubi, options)
        measurePathlineMetricKinematics(tubi, options)
        plotPathlineMetricKinematics(tubi, options)
        
        %% infer stokes forces
        measureStokesForces(tubi, options)
        
        %% Strain RATE
        measureMetricStrainRate(tubi, options)
        measureStrainRate(tubi, options)
        plotStrainRate(tubi, options)
        plotStrainRate3DFiltered(tubi, options)
        measurePathlineStrainRate(tubi, options)
        measureDxDyStrainFiltered(tubi, options)
        % Also makes fund forms in regularlized (zeta, phi) t0 Lagrangian frame
        
        % strain RATE along pathlines
        plotPathlineStrainRate(tubi, options)
        
        % Strain by following flow lines
        function strain = getCurrentPathlineStrain(tubi, t0Pathlines, varargin)
            %
            % Parameters
            % ----------
            % tubi : current class instance
            % t0Pathlines : int or empty
            %   reference timepoint. If empty, set to tubi.t0set()
            % varargin : strings ('strain', 'beltrami')
            %   which strain measures to load from disk and attribute to
            %   self
            % 
            % Returns
            % -------
            % strain : tubi.currentStrain.pathline
            
            if isempty(tubi.currentTime)
                error('Must first set time with tubi.setTime()')
            end
            if nargin < 2 
                t0Pathlines = tubi.t0set() ;
            else
                if isempty(t0Pathlines)
                    t0Pathlines = tubi.t0set() ;
                end
                if nargin < 3
                    varargin = {'strain', 'beltrami'} ;
                end
            end
            
            tubi.currentStrain.pathline.t0Pathlines = t0Pathlines ;
            if any(contains(varargin, 'strain'))
                ffn = sprintf(tubi.fullFileBase.pathlines.strain, t0Pathlines, tubi.currentTime) ;
                tubi.currentStrain.pathline.strain = load(ffn) ;
            end

            % Grab all beltramis and assign to self since stored as single
            % file on disk
            if any(contains(varargin, 'beltrami'))
                tidx = tubi.xp.tIdx(tubi.currentTime) ;
                % stores beltramis 
                tubi.getubieltramiCoefficient() ;
                tubi.currentStrain.pathline.beltrami.mu_material = ...
                    tubi.pathlines.beltrami.mu_material(tidx, :);
                tubi.currentStrain.pathline.beltrami.mu_material_filtered = ...
                    tubi.pathlines.beltrami.mu_material_filtered(tidx, :);
                tubi.currentStrain.pathline.beltrami.mu_material_vertices = ...
                    tubi.pathlines.beltrami.mu_material_vertices(tidx, :);
                tubi.currentStrain.pathline.beltrami.filterOptions = ...
                    tubi.pathlines.beltrami.filterOptions ;
            end
            if nargout > 0
                strain = tubi.currentStrain.pathline ;
            end
        end
        measurePathlineStrain(tubi, options)
        function getPathlineStrain(tubi, options)
            for tidx = 1:length(tubi.xp.fileMeta.timePoints)
                tp = tubi.xp.fileMeta.timePoints(tidx) ;
                ffn = sprintf(tubi.fullFileBase.pathlines.strain, t0Pathlines, tp) ;
                strains(tidx, :) = load(ffn) ;
            end
            tubi.pathlines.strains = strains ;
        end
        % Note: plotPathlineStrain also plots pathline beltrami coeffs
        plotPathlineStrain(tubi, options)
        function plotPathlineBeltramiKymograph(tubi, t0Pathlines, options)
            % Example usage for 2021 gut paper:
            % options = struct('ylim', [0, 2] )
            % tubi.plotPathlineBeltramiKymograph([], options)
            % 
            if nargin < 2 
                t0Pathlines = tubi.t0set() ;
            else
                if isempty(t0Pathlines)
                    t0Pathlines = tubi.t0set() ;
                end
            end 
            
            convert2hrs = false ;
            timeunits = tubi.timeUnits ;
            timepoints = tubi.xp.fileMeta.timePoints - tubi.t0 ;
            if nargin < 3
                options = struct() ;
                if contains(lower(tubi.timeUnits), 'min')
                    convert2hrs = length(tubi.xp.fileMeta.timePoints) > 60 / tubi.timeInterval ;
                    timeunits = 'hr' ;
                    timepoints = timepoints/ 60 ;
                end
            else
                if isfield(options, 'convert2hrs')
                    convert2hrs = options.convert2hrs ;
                elseif contains(lower(tubi.timeUnits), 'min')
                    convert2hrs = length(tubi.xp.fileMeta.timePoints) > 60 / tubi.timeInterval ;
                    timeunits = 'hr' ;
                    timepoints = timepoints/ 60 ;
                end
            end
            
            muAPfn = sprintf(...
                tubi.fileName.pathlines.kymographs.mu, t0Pathlines) ;
            tmp = load(muAPfn) ;
            
            re = real(tmp.mu_apM) ;
            im = imag(tmp.mu_apM) ;
            magAP = sqrt(re.^2 + im .^2) ;
            phase = atan2(im, re) ;
            % Create pullback pathline grid as mesh
            mesh = struct() ;
            [uu, vv] = meshgrid(linspace(0,1,tubi.nU), timepoints) ;
            mesh.v = [uu(:), vv(:), 0*uu(:)] ;
            pOptions = struct('cbarlabel', '$\mu$', ...
                'mesh', mesh, 'xlabel', 'ap position, $\zeta/L$', ...
                'ylabel', ['time [' timeunits ']']) ;
            if isfield(options, 'ylim') 
                pOptions.ylim = options.ylim ; 
            end
            
            close all
            hf = figure('Position', [100 100 320 320], 'units', 'centimeters');
            [meshHandle, cbs] = plotNematicField(magAP, phase, pOptions)
            axis square 
            set(gcf, 'color','w')
            % set(gcf, 'renderer', 'painters')
            
            saveas(gcf, fullfile(sprintf(...
                tubi.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0.png', ...
                t0Pathlines)))
            
            % Save axis content as png image with accompanying pdf of
            % figure frame
            FF = getframe(gca) ;
            imwrite(FF.cdata, fullfile(sprintf(...
                tubi.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0_image.png', ...
                t0Pathlines)))
            
            cla ; grid off
            axes(cbs{1})
            % Save figure frame as PDF without image content (for fusion in
            % vector graphics software)
            saveas(gcf, fullfile(sprintf(...
                tubi.dir.pathlines.kymographs, t0Pathlines), ...
                sprintf(...
                'mu_apM_kumograph_pullbackPathlines_%06dt0_frame.pdf', ...
                t0Pathlines)))
            
            % Make colorwheel
            
        end
                               
        %% timepoint-specific coordinate transformations
        sf = interpolateOntoPullbackXY(tubi, XY, scalar_field, options)
                
        %% coordSysDemo
        function coordSystemDemo(tubi, options)
            % Image for publication/presentation on method & coordinate system
            % Create coordinate system charts visualization using smoothed meshes
            coordsys = 'spsm' ;
            if isfield(options, 'coordSys')
                coordsys = options.coordSys ;
            end
            
            if strcmpi(coordsys, 'spsm')
                mesh = tubi.loadCurrentSPCutMeshSmRS() ;
            elseif strcmpi(coordsys, 'sp')
                tmp = tubi.loadCurrentSPCutMesh() ;
                mesh = struct() ;
                mesh.u = tmp.sphi ;
                mesh.v = tubi.xyz2APDV(tmp.v) ;
                mesh.nU = tmp.nU ;
                mesh.nV = tmp.nV ;
            end
            fig = figure('units', 'centimeters', 'position', [0, 0, 13, 13]) ;
            uu = mesh.u(:, 1) ;
            vv = mesh.u(:, 2) ;
            xx = (mesh.v(:, 1)) ;
            yy = (mesh.v(:, 2)) ;
            zz = (mesh.v(:, 3)) ;
            nU = mesh.nU ;
            nV = mesh.nV ;
            colorsV = viridis(mesh.nV) ;
            colorsU = viridis(mesh.nU) ;

            % LONGITUDE 3D
            subplot(2, 3, [1,2])
            hold off
            for qq = 1:mesh.nV  
                plot3(xx(qq:nU:end), yy(qq:nU:end), zz(qq:nU:end), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis equal
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
            zlabel('dv position [$\mu$m]', 'interpreter', 'latex')

            % AZIMUTH 3D
            subplot(2, 3, [4,5])
            hold off
            for qq = 1:2:mesh.nU
                inds = (qq-1)*nU+1:qq*nU ;
                plot3(xx(inds), yy(inds), zz(inds), '-', ...
                    'color', colorsU(qq, :));
                hold on ;
            end
            axis equal
            xlabel('ap position [$\mu$m]', 'interpreter', 'latex')
            ylabel('lateral position [$\mu$m]', 'interpreter', 'latex')
            zlabel('dv position [$\mu$m]', 'interpreter', 'latex')

            % LONGITUDE
            subplot(2, 3, 3)
            hold off
            for qq = 1:mesh.nV  
                plot(uu(qq:nU:end), vv(qq:nU:end), '-', ...
                    'color', colorsV(qq, :));
                hold on ;
            end
            axis square
            xlabel('$s$ [$\mu$m]', 'interpreter', 'latex')
            ylabel('$\phi$ [1/$2\pi$]', 'interpreter', 'latex')

            % AZIMUTH
            subplot(2, 3, 6)
            hold off
            for qq = 1:mesh.nU
                inds = (qq-1)*nU+1:qq*nU ;
                plot(uu(inds), vv(inds), '-', ...
                    'color', colorsU(qq, :));
                hold on ;
            end
            axis square
            xlabel('$s$ [$\mu$m]', 'interpreter', 'latex')
            ylabel('$\phi$ [1/$2\pi$]', 'interpreter', 'latex')

            set(gcf, 'color', 'w')
            set(gcf, 'visible', 'on')
        end
    end
    
    methods (Static)
        function uv = XY2uv(im, XY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map pixel positions (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT MAP IS
            % [xesz, yesz] = [size(im, 1), size(im, 2)]
            % uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            %
            % NOTE THAT INVERSE MAP IS
            % x--> (xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (Ysz-1)) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            %
            % Parameters
            % ----------
            % im : NxM numeric array or 2 ints for size
            %   image in which pixel coordinates are defined or dimensions
            %   of the image (pullback image in pixels)
            % XY : Qx2 numeric array
            %   pixel coordinates to convert to pullback space
            % doubleCovered : bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % uv : Qx2 numeric array
            %   pullback coordinates of input pixel coordinates
            
            % Input defaults
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            % Input checking
            if size(XY, 2) ~= 2
                if size(XY, 1) == 2
                    XY = XY' ;
                else
                    error('XY must be passed as #pts x 2 numeric array')
                end
            end
            % size of extended image
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            % map extended image size to (0, 1), (-0.5, 1.5) if double
            % covered. 
            % subtract 1 since pixel positions range from (1, sizeIm)
            xesz = double(Xsz - 1) ;
            yesz = double(Ysz - 1) ;
            % map from pixel y to network y (sphi)
            uv = zeros(size(XY)) ;
            % convert x axis
            uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % convert y axis
            % subtract 1 since pixel positions range from (1, sizeIm)
            if doubleCovered
                uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            else
                uv(:, 2) = vmax * (XY(:, 2) - 1) / yesz ;
            end
        end
        
        function duv = dXY2duv(im, dXY, doubleCovered, umax, vmax)
            %XY2uv(im, XY, doubleCovered, umax, vmax)
            %   Map difference in pixel positions defined on 
            %   (1, sizeImX) and (1, sizeImY) to
            %   (0, umax) and (0, vmax) of pullback space if singleCover,
            %   or y coords are mapped to (-0.5, 1.5)*vmax if doubleCover
            % 
            % NOTE THAT FULL COORDINATE MAP IS
            % [xesz, yesz] = [size(im, 1), size(im, 2)]
            % uv(:, 1) = umax * (XY(:, 1) - 1) / xesz ;
            % uv(:, 2) = vmax * 2.0 * (XY(:, 2) - 1) / yesz - 0.5 ;
            %
            % NOTE THAT INVERSE MAP IS
            % x--> (xy(:, 1) * (Xsz-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (Ysz-1)) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            %
            % Parameters
            % ----------
            % im : NxM numeric array or 2 ints for size
            %   image in which pixel coordinates are defined or dimensions
            %   of the image (pullback image in pixels)
            % dXY : Qx2 numeric array
            %   difference in pixel coordinates as vector, to convert to 
            %   pullback space
            % doubleCovered : bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % uv : Qx2 numeric array
            %   pullback coordinates of input pixel coordinates
            
            % Input defaults
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            % Input checking
            if size(dXY, 2) ~= 2
                if size(dXY, 1) == 2
                    dXY = dXY' ;
                elsed
                    error('XY must be passed as #pts x 2 numeric array')
                end
            end
            % size of extended image
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            % map extended image size to (0, 1), (-0.5, 1.5) if double
            % covered. 
            % subtract 1 since pixel positions range from (1, sizeIm)
            xesz = double(Xsz - 1) ;
            yesz = double(Ysz - 1) ;
            % map from pixel y to network y (sphi)
            duv = zeros(size(dXY)) ;
            % convert x axis
            duv(:, 1) = umax * (dXY(:, 1)) / xesz ;
            % convert y axis
            % WOULD subtract 1 since pixel positions range from (1, sizeIm)
            if doubleCovered
                duv(:, 2) = vmax * 2.0 * (dXY(:, 2)) / yesz - 0.5 ;
            else
                duv(:, 2) = vmax * (dXY(:, 2)) / yesz ;
            end
        end
        
        function XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            % XY = uv2XY(im, uv, doubleCovered, umax, vmax) 
            %   Map from pullback uv u=(0,1), v=(0,1) to pixel XY
            % x--> (xy(:, 1) * (size(im, 2)-1)) / (1*umax) + 1 , ...
            % y--> (xy(:, 2) * (size(im, 1)-1)) / (2*vmax) + 0.75 + (size(im,1)-1)*0.25
            %
            % Parameters
            % ----------
            % im : NxM numeric array or length(2) int array
            %   2D image into whose pixel space to map or size(im)
            % uv : Q*2 numeric array
            %   mesh coordinates to convert to pullback pixel space (XY)
            % doubleCovered: bool
            %   the image is a double cover of the pullback (extended/tiled
            %   so that the "top" half repeats below the bottom and the
            %   "bottom" half repeats above the top. That is, 
            %   consider im to be a double cover in Y (periodic in Y and 
            %   covers pullback space twice (-0.5 * Ly, 1.5 * Ly)
            % umax : float
            %   extent of pullback mesh coordinates in u direction (X)
            % vmax : float 
            %   extent of pullback mesh coordinates in v direction (Y)
            %   before double covering/tiling
            %
            % Returns
            % -------
            % XY : N x 2 float array
            %   positions of uv coordinates in pullback pixel space
            %
            if nargin < 3
                doubleCovered = true ;
            end
            if nargin < 4
                umax = 1.0 ;
            end
            if nargin < 5
                vmax = 1.0 ;
            end
            
            if any(size(im) > 2) 
                Xsz = size(im, 2) ;
                Ysz = size(im, 1) ;
            else
                % Interpret im as imsize
                Xsz = im(1) ;
                Ysz = im(2) ;
            end
            XY = 0*uv ;
            XY(:, 1) = uv(:, 1) * (Xsz-1) / (1*umax) + 1 ;
            if doubleCovered
                % image is double cover of physical object (periodic
                % cylinder)
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (2*vmax) + 1 + (Ysz-1)*0.25 ;
            else
                % singleCover image of physical cylindrical object
                XY(:, 2) = uv(:, 2) * (Ysz-1) / (1*vmax) + 1  ;
            end        
        end
        
        function [xx, yy] = clipXY(xx, yy, Lx, Ly)
            % Clip x at (1, Lx) and clip Y as periodic (1=Ly, Ly=1), for
            % image that is periodic in Y. Consider Y in [1, Ly]. If the
            % pullback is a doubleCover, Ly = 2*mesh width in pixels
            %
            % Parameters
            % ----------
            % xx : 
            % yy : 
            % Lx : int
            %   number of pixels along x dimension
            % Ly : int
            %   number of pixels along y dimension
            % 
            % Returns
            % -------
            % [xx, yy] : x and y coordinates clipped to [1, Lx] and [1, Ly]
            
            % Note we use minimum values of 1 (in pixels)
            minX = 1 ;
            minY = 1 ;
            
            % Clip in X
            xx(xx > Lx) = Lx ;
            xx(xx < minX ) = 1 ;
            % modulo in Y
            yy(yy > Ly) = yy(yy > Ly) - Ly + minY;
            yy(yy < minY) = yy(yy < minY) + Ly ;
        end
        
        function XY = doubleToSingleCover(XY, Ly)
            % detect if XY is passed as a pair of grids
            if length(size(XY)) > 2 && size(XY, 2) > 2
                % XY is a pair of position grids each as 2D arrays. Clip Y
                tmp = XY(:, :, 2) ;
                tmp(tmp < Ly * .25) = tmp(tmp < Ly * .25) + Ly * 0.5 ;
                tmp(tmp > Ly * .75) = tmp(tmp > Ly * .75) - Ly * 0.5 ;
                XY(:, :, 2) = tmp ;
            elseif size(XY, 2) == 2
                % XY is input as Nx2 array
                XY(XY(:, 2) < Ly * .25, 2) = XY(XY(:, 2) < Ly * .25, 2) + Ly * 0.5 ;
                XY(XY(:, 2) > Ly * .75, 2) = XY(XY(:, 2) > Ly * .75, 2) - Ly * 0.5 ;
            else
                error('XY must be passed as Nx2 or QxRx2 array')
            end
        end
        
        % function uv2pix_old(im, aspect)
        %     % map from network xy to pixel xy
        %     % Note that we need to flip the image (Yscale - stuff) since saved ims had
        %     % natural ydirection.
        %     % Assume here a coord system xy ranging from (0, xscale) and (0, 1) 
        %     % maps to a coord system XY ranging from (0, Yscale * 0.5) and (0, Yscale)
        %     x2Xpix = @(x, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     % y2Ypix = @(y, h, Yscale) Yscale - (Yscale*0.5)*(y+0.5) + h ;
        %     y2Ypix = @(y, h, Yscale) (Yscale*0.5)*(y+0.5) + h ;
        % 
        %     dx2dX = @ (y, Yscale, xscale) (Yscale * 0.5) * aspect * x / xscale ;
        %     dy2dY = @ (y, Yscale) (Yscale*0.5)*y ;
        % 
        %     % Now map the coordinates
        % end
        
        [cutMesh, cutMeshC] = doubleResolution(cutMesh, preview)
        
        function [mag_ap, theta_ap] = dvAverageNematic(magnitude, theta)
            %[mag_ap, theta_ap] = DVAVERAGENEMATIC(magnitude, theta)
            % Given a nematic field defined on a rectilinear grid, with 
            % Q(i,j) being in the ith ap position and jth dv position,
            % average the nematic field over the dv positions (dimension 2)
            % 
            % Parameters
            % ----------
            % magnitude : nU x nV numeric array
            %   magnitude of nematic field in 2D rectilinear grid
            % theta : nU x nV float array, with values as angles in radians
            %   angle (mod pi) of nematic field in 2D rectilinear grid
            %
            % Returns
            % -------
            % mag_ap : nU x 1 float array
            %   average magnitude of nematic field along the ap dimension
            % theta_ap : nU x 1 float array
            %   average angle (mod pi) of nematic field along the ap 
            %   dimension
            ap_x = magnitude .* cos(2*theta) ;
            ap_y = magnitude .* sin(2*theta) ;
            ap_xy = [mean(ap_x, 2) , mean(ap_y, 2)] ;
            mag_ap = vecnorm(ap_xy, 2, 2) ;
            theta_averages = atan2(ap_xy(:, 2), ap_xy(:, 1)) ;
            theta_ap = 0.5 * mod(theta_averages, 2*pi) ;
        end
       
        function invRot = invertRotation(rot)
            % Obtain rotation matrix that undoes the APDV rotation 
            rotM = [rot(1, :), 0; rot(2,:), 0; rot(3,:), 0; 0,0,0,1] ;
            tform = affine3d(rotM) ;
            invtform = invert(tform) ;
            invRot = invtform.T(1:3,1:3) ;
        end
        
        function [dorsal, ventral, left, right] = quarterIndicesDV(nV)
            %[dorsal, ventral, left, right] = quarterIndicesDV(nV)
            % indices for each quarter of a DV section in grid coordinates
            if nargin < 1
                nV = tubi.nV ;
            end            
            q0 = round(nV * 0.125) ;
            q1 = round(nV * 0.375) ;
            q2 = round(nV * 0.625) ;
            q3 = round(nV * 0.875) ;
            left = q0:q1 ;
            ventral = q1:q2 ;
            right = q2:q3 ;
            dorsal = [q3:nV, 1:q1] ;
        end
    end
    
end