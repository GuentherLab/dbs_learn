%%%% top level script for running audio syncing between behavioral audio files with CED
..... and add CED-synced timepoints to behavioral annotation tables
..... so that we can later get percept-responses synced to these timepoints


clc; clear

datetime

op.sub = 'sap022'; 
op.ses = 'stut'; 
op.grouping_columns = {'step','task'}; 
op.exclude_criteria = {'filename','sub-sap022_ses-stut_recording-thinkpad-lapelmic.wav'}; 
% op.include_criteria = {'step','Z001'}; 
% op.include_criteria = {'step',{'Z001','C001'}}; 

op.show_qc_figs  = 1; 


paths = setpaths_dbs_learn(op); 


sync_in = paths.landmarks_file; 


op.annot_dir = paths.annot; % contains annot files to give synced times to
op.sync_table_savename = paths.sync_ses; 
op.qc_figs_savedir = paths.sync_qc; 
op.annotation_file_suffixes = {'annot-speech'}; 
op.event_wise_audvid_dir = paths.event_audvid; 
    op.event_audvid_buffer_sec = [1 1]; 



[sync_out, op_out] = sync_audio_and_annot(sync_in, op)