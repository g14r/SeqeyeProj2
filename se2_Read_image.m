%% Preprocess the MRI
% read and align the MRI
mri             = ft_read_mri('preop_mri.nii.gz'); % or a single file of a DICOM series

cfg             = [];
cfg.method      = 'interactive';
cfg.coordsys    = 'tal'; % anterior & posterior commissure, positive midline, and right side
mri = ft_volumerealign(cfg, mri);

% write to file, for later processing
ft_write_mri('Subject_1_MR_preproc.nii', mri.anatomy, 'transform', mri.transform, 'dataformat', 'nifti_spm');

%% Preprocess the CT and co-register to the MR
ct = ft_read_mri('postop_ct.nii.gz'); % or a single file of a DICOM series

cfg             = [];
cfg.method      = 'interactive'; % nasion, left & right ear, and positive midline
ct = ft_volumerealign(cfg, ct);
ct = ft_convert_coordsys(ct, 'tal'); % convert ctf to spm/tal coordsys

% co-register the CT to the MRI
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'tal';
cfg.viewresult  = 'yes'; % view realignment result
ct = ft_volumerealign(cfg, ct, mri);

% write to file, as a backup
ft_write_mri('Subject_1_CT_coreg.nii', ct.anatomy, 'transform', ct.transform, 'dataformat', 'nifti_spm');

%% Localize electrodes in the CT
load('SubjectUCI29_hdr.mat'); % Contains the hdr
cfg             = [];
cfg.channel     = hdr.label; % these will be assigned to a location
elec = ft_electrodeplacement(cfg, ct);

% write to file, for later visualization
save('SubjectUCI29_elec.mat', 'elec');







%% Preprocess the CT and co-register to the MR
ct_mr = ft_read_mri('ct_to_mri.nii.gz'); % or a single file of a DICOM series

cfg             = [];
cfg.method      = 'interactive'; % nasion, left & right ear, and positive midline
ct_mr = ft_volumerealign(cfg, ct_mr);
ct = ft_convert_coordsys(ct_mr, 'tal'); % convert ctf to spm/tal coordsys

% co-register the CT to the MRI
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'tal';
cfg.viewresult  = 'yes'; % view realignment result
