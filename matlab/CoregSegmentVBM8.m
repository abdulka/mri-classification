%% Coregister and Segment Structural MRI Images
%
%  This function requires SPM8 (r5236) and VBM8 (r435).
%
%  CoregSegment(T1imagefile,doCoreg,doSetOrigin)
%
%  Inputs:
%   T1imagefile : file name of structural T1 image to be processed
%   doCoreg     : logical variable that determines whether coregistration
%                 is done or not
%   doSetOrigin : logical variable that teremines wheter the origin is set
%                 to the center of the image or not
%
%  Outputs:
%   (realigned and segmented images output to file system)
%
% Autohr:
% Ahmed Abdulkadir

function CoregSegmentVBM8(T1,doCoreg,doSetOrigin)

% set default to no coregistration and no reset of the origin
if ~exist('doCoreg',    'var'),     doCoreg = false; end
if ~exist('doSetOrigin','var'), doSetOrigin = false; end
% make sure that the variables are logical
assert(islogical(doCoreg),    'variable doCoreg must be logical');
assert(islogical(doSetOrigin),'variable doSetOrigin must be logical');

% check if SPM is available
assert(~isempty(which('spm')),'ERROR: SPM not loaded. Abording.');

% check if input image file exists
assert(exist(T1,'file') > 0, ['ERROR: input T1 image file does '         ...
    'not exist %s'], T1);
% set reference (fixed) image
reference = fullfile(spm('dir'),'canonical','single_subj_T1.nii,1');
assert(exist(reference(1:end-2),'file') > 0, ['ERROR: desired reference'   ...
    ' image does'...
    ' not exist']);

% input image
[T1pth,T1fname,x] = fileparts(T1);

% center origin if requested
if doSetOrigin,
    setOriginToImageCenter(T1);
end

% run SPM batch
jobs = repmat({ defineSegmentJobVBM8(T1,doCoreg) },1,1);
spm_jobman('serial', jobs,'',{});


% _______________________________________________________________________ %

% _______________________________________________________________________ %
function matlabbatch = defineSegmentJobVBM8(T1,doCoreg)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.cfg_named_file.name = 'sMRI Images (T1)';
matlabbatch{1}.cfg_basicio.cfg_named_file.files = { {T1} }';
% Coregistration T1->MNI
if doCoreg,
    matlabbatch{end+1}.spm.spatial.coreg.estimate.ref = {                  ...
        fullfile(spm('dir'),'canonical','single_subj_T1.nii,1') };
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1) = cfg_dep;
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).tname = 'Source Image';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(1).name = 'class';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(1).value = 'cfg_files';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(2).name = 'strtype';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).tgt_spec{1}(2).value = 'e';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).sname = 'Named File Selector: sMRI Images (T1)(1) - Files';
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
    matlabbatch{end}.spm.spatial.coreg.estimate.source(1).src_output = substruct('.','files', '{}',{1});
    matlabbatch{end}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02  ...
        0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
end
% Segmentation of T1
matlabbatch{end+1}.spm.tools.vbm8.estwrite.data(1) = cfg_dep;
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).tname = 'Volumes';
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).tgt_spec{1}(1).name = 'class';
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).sname = ['Named File '    ...
    'Selector: sMRI Images (T1,FLAIR)(1) - Files'];
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).src_exbranch =            ...
    substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{end}.spm.tools.vbm8.estwrite.data(1).src_output =              ...
    substruct('.','files', '{}',{1});
spmversion = spm('ver');
switch spmversion
    case 'SPM8',
        matlabbatch{end}.spm.tools.vbm8.estwrite.opts.tpm = {
            fullfile(spm('dir'),'toolbox','Seg','TPM.nii') };
    case 'SPM12b',
        error('Incompatible SPM version: %s\nPlease use SPM8 instead.\n',  ...
            spmversion);
    otherwise,
        error('Unknown or incompatible SPM version: %s\n',spmversion);
end
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.ngaus = [2 2 2 3 4 2];
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.biasreg = 0.0001;
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.biasfwhm = 60;
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.affreg = 'mni';
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.warpreg = 4;
matlabbatch{end}.spm.tools.vbm8.estwrite.opts.samp = 3;
matlabbatch{end}.spm.tools.vbm8.estwrite.extopts.dartelwarp = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.extopts.ornlm = 0.7;
matlabbatch{end}.spm.tools.vbm8.estwrite.extopts.mrf = 0.15;
matlabbatch{end}.spm.tools.vbm8.estwrite.extopts.cleanup = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.extopts.print = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.GM.native = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.GM.warped = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.GM.modulated = 2;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.GM.dartel = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.WM.native = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.WM.warped = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.WM.modulated = 2;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.WM.dartel = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.CSF.native = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.CSF.warped = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.CSF.modulated = 2;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.CSF.dartel = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.bias.native = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.bias.warped = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.bias.affine = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.label.native = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.label.warped = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.label.dartel = 1;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.jacobian.warped = 0;
matlabbatch{end}.spm.tools.vbm8.estwrite.output.warps = [1 1];


% _______________________________________________________________________ %
function setOriginToImageCenter(fnames)
% check if argumant has been given
if ~exist('fnames','var');
    fnames = spm_select([1 Inf],'image','Select 3D Images ..','');
end
% transform char array of images into cell array of char
if ischar(fnames),
    fnames = cellstr(fnames);
end
% N: number of images
N = size(fnames,1);
% apply translation
for i=1:N, % for all images
    % V: memory-mapped volumes
    V = spm_vol(fnames{i});
    % center.vx: image center in voxel coordinates
    center.vx = [V.dim(1)/2 V.dim(2)/2 V.dim(3)/2 1];
    % center.vx: image center in world coordinates
    center.mm = (V.mat*center.vx')';
    % apply translation
    mat = spm_matrix(-center.mm(1:3))*V.mat;
    % write header
    spm_get_space(V.fname,mat);
end
% _______________________________________________________________________ %

