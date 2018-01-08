%% Extract structural features
function [voxelGM,voxelGMs,voxelGMf,voxelWM,voxelWMs,voxelWMf,...
    voxelCSF,voxelCSFs,voxelCSFf,...
    regionGM,regionWM,regionCSF,vol,LPBAGM,LPBAWM,  ...
    ROIGM, ...
    GMmask,WMmask,CSFmask,GMmaskLowRes,WMmaskLowRes,CSFmaskLowRes,GMAtlas,WMAtlas] = ...
    extractstructuralfeatures(VBM8folders,basenames)

% tGM: template GM
tGM = spm_vol(fullfile('S:/LPBA40/LPBA40.SPM5.nifti', ...
    '/tissue/lpba40.spm5.avg152T1.gm.nii'));
tWM = spm_vol(fullfile('S:/LPBA40/LPBA40.SPM5.nifti', ...
    '/tissue/lpba40.spm5.avg152T1.wm.nii'));
tCSF = spm_vol(fullfile('S:/LPBA40/LPBA40.SPM5.nifti', ...
    '/tissue/lpba40.spm5.avg152T1.csf.nii'));

% mask GM and WM images
GMmask  = tGM.private.dat(:)>0;
WMmask  = tWM.private.dat(:)>0;
CSFmask  = tCSF.private.dat(:)>0;
% low resolution GM and WM mask
GMmaskLowRes  = tGM.private.dat(1:2:end,1:2:end,1:2:end)>0;
WMmaskLowRes  = tWM.private.dat(1:2:end,1:2:end,1:2:end)>0;
CSFmaskLowRes  = tCSF.private.dat(1:2:end,1:2:end,1:2:end)>0;
GMmaskLowRes  = GMmaskLowRes(:);
WMmaskLowRes  = WMmaskLowRes(:);
CSFmaskLowRes  = CSFmaskLowRes(:);

% vxx,vxy,vxz: voxel positions
[vxx,vxy,vxz] = ndgrid(1:tGM.dim(1),1:tGM.dim(2),1:tGM.dim(3));
% pos_mm: mm positions (full resolution)
mmXYZ       = tGM.mat*[vxx(:),vxy(:),vxz(:),ones(numel(vxx),1)]';
% vxx,vxy,vxz: voxel positions
[vxx,vxy,vxz] = ndgrid(1:2:tGM.dim(1),1:2:tGM.dim(2),1:2:tGM.dim(3));
% mmXYZ: mm positions (half resolution)
mmXYZlowres = tGM.mat*[vxx(:),vxy(:),vxz(:),ones(numel(vxx),1)]';

clear vxx vxy vxz

% load atlas data
fprintf('Loading LONI probabilistic atlas .. ');
GMAtlas   = getGMAtlas;
LPBAGM    = reshape(permute(spm_read_vols(spm_vol(char(GMAtlas))),[4 1 2 3]),[56 prod(tGM.dim)]);
LPBAGM    = LPBAGM(:,GMmask);
LPBAGM(~isfinite(LPBAGM)) = 0;
sumLPBAGM = sum(LPBAGM,2);
WMAtlas   = getAtlas;
LPBAWM      = reshape(permute(spm_read_vols(spm_vol(char(WMAtlas))),[4 1 2 3]),[56 prod(tGM.dim)]);
LPBAWM      = LPBAWM(:,WMmask);
LPBAWM(~isfinite(LPBAWM))     = 0;
sumLPBAWM   = sum(LPBAWM,2);
CSFAtlas   = getAtlas;
LPBACSF      = reshape(permute(spm_read_vols(spm_vol(char(CSFAtlas))),[4 1 2 3]),[56 prod(tGM.dim)]);
LPBACSF      = LPBACSF(:,CSFmask);
LPBACSF(~isfinite(LPBACSF))     = 0;
sumLPBACSF   = sum(LPBACSF,2);
fprintf('Done!\n');

if ischar(basenames),
    basenames = { basenames };
end

% GM: gray matter probabilities segmentation
fprintf('Getting handles to %i GM images .. ',numel(VBM8folders));
GM = cellfun(@(SPM12folder,basename) spm_vol( ...
    fullfile(SPM12folder,['smwc1' basename '.nii'])),VBM8folders,basenames,...
    'ErrorHandler',@(err,~,~) [] ,'UniformOutput',false);
fprintf('Done!\n');

% WM: gray matter probabilities segmentation
fprintf('Getting handles to %i WM images .. ',numel(VBM8folders));
WM = cellfun(@(SPM12folder,basename) spm_vol( ...
    fullfile(SPM12folder,['smwc2' basename '.nii'])),VBM8folders,basenames,...
    'ErrorHandler',@(err,~,~) [],'UniformOutput',false);
fprintf('Done!\n');

% CSF: gray matter probabilities segmentation
fprintf('Getting handles to %i CSF images .. ',numel(VBM8folders));
CSF = cellfun(@(SPM12folder,basename) spm_vol( ...
    fullfile(SPM12folder,['smwc3' basename '.nii'])),VBM8folders,basenames,...
    'ErrorHandler',@(err,~,~) [],'UniformOutput',false);
fprintf('Done!\n');

% set sample positions
vxXYZ       = GM{1}.mat\mmXYZ;
vxXYZlowres = GM{1}.mat\mmXYZlowres;

% smooth images with 8mm FWHM
FWHM = 8;
GMs  = GM;
WMs  = WM;
CSFs  = CSF;
voxelGM  = NaN(numel(GM) ,sum(GMmask,1),'single');
voxelWM  = NaN(numel(WM) ,sum(WMmask,1),'single');
voxelCSF  = NaN(numel(WM) ,sum(CSFmask,1),'single');
voxelGMf = cell(numel(GM),1);
voxelWMf = cell(numel(GM),1);
voxelCSFf = cell(numel(GM),1);
voxelGMs = NaN(numel(GMs),sum(GMmaskLowRes,1),'single');
voxelWMs = NaN(numel(WMs),sum(WMmaskLowRes,1),'single');
voxelCSFs = NaN(numel(CSFs),sum(CSFmaskLowRes,1),'single');
vol = NaN(numel(basenames),3);
ROIGM = NaN(numel(basenames),7);

parfor i=1:numel(GMs)
    try
        % GM
        fprintf('%-04i %s\n',i,GM{i}.fname);
        Vs = NaN(GMs{i}.dim,'single');
        spm_smooth(GMs{i},Vs,[FWHM FWHM FWHM],'single');
        GMs{i}.dat      = Vs;
        GMs{i}.pinfo(3) = 0;        
        voxelGM(i,:)  = single(spm_sample_vol(GM {i},vxXYZ(1,GMmask),vxXYZ(2,GMmask),vxXYZ(3,GMmask),1));
        voxelGMs(i,:) = single(spm_sample_vol(GMs{i},vxXYZlowres(1,GMmaskLowRes),...
            vxXYZlowres(2,GMmaskLowRes),vxXYZlowres(3,GMmaskLowRes),1));
        GMs{i}.dat      = [];

        % WM
        fprintf('%-04i %s\n',i,WM{i}.fname);
        Vs = NaN(WMs{i}.dim,'single');
        spm_smooth(WMs{i},Vs,[FWHM FWHM FWHM],'single');
        WMs{i}.dat      = Vs;
        WMs{i}.pinfo(3) = 0;
        voxelWM (i,:) = single(spm_sample_vol(WM {i},vxXYZ(1,WMmask),vxXYZ(2,WMmask),vxXYZ(3,WMmask),1));
        voxelWMs(i,:) = single(spm_sample_vol(WMs{i},vxXYZlowres(1,WMmaskLowRes),...
            vxXYZlowres(2,WMmaskLowRes),vxXYZlowres(3,WMmaskLowRes),1));
        WMs{i}.dat      = [];

        % CSF
        fprintf('%-04i %s\n',i,CSF{i}.fname);
        Vs = NaN(CSFs{i}.dim,'single');
        spm_smooth(CSFs{i},Vs,[FWHM FWHM FWHM],'single');
        CSFs{i}.dat      = Vs;
        CSFs{i}.pinfo(3) = 0;
        voxelCSF (i,:) = single(spm_sample_vol(CSF {i},vxXYZ(1,CSFmask),vxXYZ(2,CSFmask),vxXYZ(3,CSFmask),1));
        voxelCSFs(i,:) = single(spm_sample_vol(CSFs{i},vxXYZlowres(1,CSFmaskLowRes),...
            vxXYZlowres(2,CSFmaskLowRes),vxXYZlowres(3,CSFmaskLowRes),1));
        CSFs{i}.dat      = [];
        
        % tissue volumes
        vol(i,:) = [sum(sum(sum(GM{i}.private.dat(:,:,:))))
                sum(sum(sum(WM{i}.private.dat(:,:,:))))
                sum(sum(sum(CSF{i}.private.dat(:,:,:))))]'/abs(det(GM{i}.mat))/10^5;
        ni = GM{i};
        % {
        ROIGM(i,:) = [ ...
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[  6 -48  66]'),@mean) ... Precuneus R
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[ 24 -14 -18]'),@mean) ... Parahyppocampal gyrus R
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[-24 -14 -16]'),@mean) ... Amygdala L
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[ -2 -68  54]'),@mean) ... Precuneus L
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[ 16  50 -20]'),@mean) ... Anterior cingulate cortex R
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[-10  -2   6]'),@mean) ... Thalamus L
        spm_summarise(ni,struct('def','sphere','spec',5.01,'xyz',[  4  14 -20]'),@mean)]; % Ventromedial frontal cortex R
        %}
    catch e
        warning('Error occured in iteration %i\n%s',i,e.message);
    end
end
voxelGM(~isfinite(voxelGM))   = 0;
voxelWM(~isfinite(voxelWM))   = 0;
voxelCSF(~isfinite(voxelCSF))   = 0;
voxelGMs(~isfinite(voxelGMs)) = 0;
voxelWMs(~isfinite(voxelWMs)) = 0;
voxelCSFs(~isfinite(voxelCSFs)) = 0;

% summarize GM and WM over regions
regionGM = cell2mat(arrayfun(@(i) ...
    sum(bsxfun(@times,voxelGM,LPBAGM(i,:)),2)./sumLPBAGM(i),(1:56),'UniformOutput',false));
regionWM = cell2mat(arrayfun(@(i) ...
    sum(bsxfun(@times,voxelWM,LPBAWM(i,:)),2)./sumLPBAWM(i),(1:56),'UniformOutput',false));
regionCSF = cell2mat(arrayfun(@(i) ...
    sum(bsxfun(@times,voxelCSF,LPBACSF(i,:)),2)./sumLPBACSF(i),(1:56),'UniformOutput',false));
end

function Atlas = getGMAtlas()
% /Users/abdulka/Research/Neuroimaging/LPBA40.SPM5.nifti/
%Atlas = cellfun(@(x) fullfile('/projects/aalib/var/LPBA40.SPM5.nifti',x), { ...
Atlas = cellfun(@(x) fullfile('S:/LPBA40/LPBA40.SPM5.nifti',x), { ...    
    'PDFGM/lpba40.spm5.avg152T1.brainstem.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.cerebellum.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.angular_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.caudate.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.cingulate_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.cuneus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.fusiform_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.gyrus_rectus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.hippocampus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.inferior_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.inferior_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.inferior_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.insular_cortex.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.lateral_orbitofrontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.lingual_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.middle_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.middle_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.middle_orbitofrontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.middle_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.parahippocampal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.postcentral_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.precentral_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.precuneus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.putamen.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.superior_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.superior_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.superior_parietal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.superior_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.L.supramarginal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.angular_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.caudate.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.cingulate_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.cuneus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.fusiform_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.gyrus_rectus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.hippocampus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.inferior_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.inferior_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.inferior_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.insular_cortex.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.lateral_orbitofrontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.lingual_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.middle_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.middle_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.middle_orbitofrontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.middle_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.parahippocampal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.postcentral_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.precentral_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.precuneus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.putamen.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.superior_frontal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.superior_occipital_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.superior_parietal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.superior_temporal_gyrus.gm.pdf.nii'
'PDFGM/lpba40.spm5.avg152T1.R.supramarginal_gyrus.gm.pdf.nii'},...
    'UniformOutput',false);

end

function Atlas = getAtlas()
%Atlas = cellfun(@(x) fullfile('/projects/aalib/var/LPBA40.SPM5.nifti',x), { ...
Atlas = cellfun(@(x) fullfile('S:/LPBA40/LPBA40.SPM5.nifti',x), { ...
'PDF/lpba40.spm5.avg152T1.brainstem.pdf.nii'
'PDF/lpba40.spm5.avg152T1.cerebellum.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.angular_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.caudate.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.cingulate_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.cuneus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.fusiform_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.gyrus_rectus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.hippocampus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.inferior_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.inferior_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.inferior_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.insular_cortex.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.lateral_orbitofrontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.lingual_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.middle_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.middle_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.middle_orbitofrontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.middle_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.parahippocampal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.postcentral_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.precentral_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.precuneus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.putamen.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.superior_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.superior_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.superior_parietal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.superior_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.L.supramarginal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.angular_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.caudate.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.cingulate_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.cuneus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.fusiform_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.gyrus_rectus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.hippocampus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.inferior_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.inferior_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.inferior_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.insular_cortex.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.lateral_orbitofrontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.lingual_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.middle_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.middle_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.middle_orbitofrontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.middle_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.parahippocampal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.postcentral_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.precentral_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.precuneus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.putamen.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.superior_frontal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.superior_occipital_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.superior_parietal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.superior_temporal_gyrus.pdf.nii'
'PDF/lpba40.spm5.avg152T1.R.supramarginal_gyrus.pdf.nii'},...
    'UniformOutput',false);
end
