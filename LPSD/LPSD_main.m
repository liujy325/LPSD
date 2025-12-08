function [H,rmse,matchedPoints1,matchedPoints2,params] = LPSD_main(im1, im2, varargin)
% LPSD_v1.0
% Written by liujy

% Input validation
if nargin < 2
    error('LPSD_main:InsufficientInputs', 'At least two images are required');
end

% Parse input parameters
if nargin == 3 && isstruct(varargin{1})
    % Structure input
    params = varargin{1};
else
    % Name-value pairs or positional arguments
    if nargin >= 3 && isnumeric(varargin{1})
        % Legacy positional arguments
        params = struct();
        if nargin >= 3, params.nkpts = varargin{1}; end
        if nargin >= 4, params.patch_size = varargin{2}; end
        if nargin >= 5, params.nbin = varargin{3}; end
        if nargin >= 6, params.nblock = varargin{4}; end
        if nargin >= 7, params.r_num = varargin{5}; end
        if nargin >= 8, params.t_num = varargin{6}; end
    else
        % Name-value pairs
        params = struct();
        for i = 1:2:length(varargin)
            params.(varargin{i}) = varargin{i+1};
        end
    end
end

% Set default parameters
defaultParams = struct(...
    'nkpts', 5000, ...
    'patch_size', 128, ...
    'nbin', 6, ...
    'nblock', 6, ...
    'r_num', 36, ...
    't_num', 36);

% Merge with user parameters
fields = fieldnames(defaultParams);
for i = 1:length(fields)
    if ~isfield(params, fields{i})
        params.(fields{i}) = defaultParams.(fields{i});
    end
end

% Extract parameters
nkpts = params.nkpts;
patch_size = params.patch_size;
nbin = params.nbin;
nblock = params.nblock;
r_num = params.r_num;
t_num = params.t_num;

% Convert to grayscale if needed
if size(im1,3)==3
    im1 = rgb2gray(im1);
end
if size(im2,3)==3
    im2 = rgb2gray(im2);
end

%% Feature Detection
[kpts1,~] = detector(im1,nkpts);
[kpts2,~] = detector(im2,nkpts);

%% Feature Description
[des1,kpts1_] = descriptor_LPS(im1,kpts1',patch_size,nblock,nbin,r_num,t_num);
[des2,kpts2_] = descriptor_LPS(im2,kpts2',patch_size,nblock,nbin,r_num,t_num);

%% Feature Matching
[H,~,~,~] = LPSD_match(des1',kpts1_,des2',kpts2_,'similarity',10,1);

%% Rematching
[kpts1,kpts2] = kpts_elimination(kpts1,kpts2,H,10);
[des1,kpts1] = descriptor_LPS_ref(im1,kpts1',patch_size,nblock,nbin,r_num,t_num);
[des2,kpts2] = descriptor_LPS_sen(im2,kpts2',patch_size,nblock,nbin,r_num,t_num,H);
[H,rmse,matchedPoints1,matchedPoints2] = LPSD_match(des1',kpts1,des2',kpts2,'affine',3,0);