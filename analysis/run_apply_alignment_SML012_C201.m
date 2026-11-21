%% 1) Local setup: keep code dependencies inside AnlysCode
% This runner should only need code from this repository plus an external
% FieldTrip installation already on the MATLAB path.
%
% Local code paths added here:
%   repo_root      - AnlysCode root
%   demo           - demo_apply_alignment.m
%   helper         - shared MATLAB helpers used by this project
%   packages       - local third-party MATLAB utilities bundled in repo


%% Demo: inspect aligned SML012 C201 MSL recording
% This script reads the preproc MAT file and creates diagnostics only.
% It does not modify the source preproc file or write outputs.

clearvars;
clc;

this_file = mfilename('fullpath');
demo_dir = fileparts(this_file);
repo_root = fileparts(demo_dir);
helper_dir = fullfile(repo_root, 'helper');
packages_dir = fullfile(repo_root, 'packages');

addpath(repo_root);
addpath(demo_dir);
if exist(helper_dir, 'dir') == 7
    addpath(helper_dir);
end
if exist(packages_dir, 'dir') == 7
    addpath(genpath(packages_dir));
end

fdir_dbs = "C:\Users\speech\Partners HealthCare Dropbox\Andrew Meier\R01-SML_data_shared\sourcedata\sml012\ses-multisyl\dbs-recordings\";
fdir_preproc = fdir_dbs;

fpath_preproc = fullfile( ...
    fdir_preproc, 'task-trainB\', ...
    'C201_MSL-Train-B_DBS-OFF_SURVEY_preproc_260812-230758.mat');

assert(exist(fpath_preproc, 'file') == 2, ...
    'Preproc file was not found:\n%s', fpath_preproc);
assert(exist('ft_preprocessing', 'file') == 2, ...
    ['FieldTrip is not on the MATLAB path. Add FieldTrip and run ft_defaults ', ...
     'before running this demo.']);

%% 2) What this file is and which fields the demo reads
% File name:
%   C201_MSL-Train-B_DBS-OFF_SURVEY_preproc_260812-230758.mat
%
% Subject / session association:
%   SML012                - subject ID (from this demo and parent folder)
%   ses-multisyl          - multisyllable session (from parent folder)
%   C201                  - step / recording ID for this block
%   MSL-Train-B           - multisyllable learning, Train-B block
%   DBS-OFF               - stimulation state during this recording
%   SURVEY                - Percept recording mode
%   260812-230758         - saved timestamp (YYMMDD-HHMMSS)
%
% Top-level variables expected by demo_apply_alignment:
%   CED                   - preprocessed CED-side signals and event tables
%   PCP                   - preprocessed Percept/PCP signals
%   m_Align               - alignment metadata and time-translation handles
%   bhv                   - behavior metadata, including bhv.evTrial when present
%
% CED fields confirmed from the MAT-file metadata:
%   SR, mic, voice, audStim, audioFR, fpath_voice, fpath_audStimHP,
%   fpath_audioFRHP, Sync, cfgE, evPCP, stimArt, EMG, KbEvent, evTask,
%   fileStartTime_dt, chanInfo
%
% Subfields used directly by this demo:
%   CED.evPCP.sync_t      - PCP sync pulse times in CED clock
%   CED.evTask.sync_t     - task sync pulse times in CED clock
%   CED.KbEvent           - keyboard/event table; demo uses Time_sec and Code
%   CED.stimArt           - CED-time reference and optional stim-artifact trace
%   CED.voice / mic / audStim
%                         - source streams for voice-envelope construction
%   CED.EMG               - optional EMG channels for plotting
%
% PCP / alignment fields used directly by this demo:
%   PCP.data              - FieldTrip-like continuous PCP/LFP struct
%   m_Align.fcn_translateTime
%   m_Align.TimeTrans_PCP_to_CED_byEvent
%   m_Align.TimeTrans_BHV_to_CED
%   m_Align.TimeTrans_TRL_to_CED
%   m_Align.t_tskSync_BHV
%   m_Align.t_pcpSync_CED
%
% BHV table expectations:
%   bhv.evTrial           - trial table with t_* event columns
%   timingRef             - original clock reference per trial row
%   t_go_aud_on,
%   t_stim_vis_on,
%   t_go_aud_off,
%   t_stim_aud_on         - landmark columns used for translation/plotting
%
% Notes:
%   The subject / session / task / DBS / mode items above are inferred from the
%   filename and parent folder. The struct names and confirmed CED fields were
%   recovered from the MAT-file metadata itself.

% Start with an overview that remains stable in the MATLAB desktop renderer.
% Full-session plotting with all yyaxis overlays can exhaust the graphics scene.
plot_dur_s = 120;

% demo_apply_alignment inputs:
%   fpath_preproc       - full path to the preprocessed MAT file to load
%   'msl'               - task label; enables multisyllable-learning trial
%                         event lines in the diagnostic figure
%   show_LFP            - plot one aligned LFP time-series panel per channel
%   show_LFP_PSD        - add a power-spectrum panel beside each LFP channel
%   show_voice          - plot the selected voice-envelope source
%   show_voice_raw      - overlay the unprocessed voice waveform (right axis)
%   show_EMG            - plot EMG channels when they are present in CED
%   show_stimArt        - overlay CED stimulation-artifact trace on LFP panels
%   show_KbEvents       - mark CED.KbEvent times and event-code labels
%   show_syncLines      - mark PCP/task/BHV synchronization events
%   voice_env_channel   - CED source used for the envelope: voice, mic, or audStim
%   psd_fmax            - highest frequency shown in the LFP PSD panels (Hz)
%   plot_dur_s          - initial CED-time display duration in seconds
%
% demo_apply_alignment outputs:
%   lfp_aligned         - FieldTrip LFP struct after filtering; .trial{1} is
%                         [channels x samples], and .time{1} is in CED seconds
%   t_PCP_inCED         - PCP/LFP sample-time vector translated into CED seconds
%   t_bhv_inCED         - first available BHV landmark per trial in CED seconds;
%                         uses t_go_aud_on, t_stim_vis_on, or t_stim_aud_on
%   ppCED_out           - CED struct returned by alignment, including voiceEnv
%                         when a usable voice/mic/audStim stream is available
%   tbl_BHV             - behavioral trial table with every available t_* event
%                         column translated from its source clock to CED seconds
[lfp_aligned, t_PCP_inCED, t_bhv_inCED, ppCED_out, tbl_BHV] = ...
    demo_apply_alignment( ...
        fpath_preproc, 'msl', ...
        'show_LFP',          true, ...
        'show_LFP_PSD',      true, ...
        'show_voice',        true, ...
        'show_voice_raw',    false, ...
        'show_EMG',          true, ...
        'show_stimArt',      false, ...
        'show_KbEvents',     true, ...
        'show_syncLines',    true, ...
        'voice_env_channel', 'voice', ...
        'psd_fmax',          150, ...
        'plot_dur_s',        plot_dur_s);

if 'print to command window'
    fprintf('\n--- SML012 C201 alignment summary ---\n');
    fprintf('LFP channels: %d\n', numel(lfp_aligned.label));
    fprintf('LFP samples:  %d at %.0f Hz\n', ...
        numel(lfp_aligned.time{1}), lfp_aligned.fsample);
    fprintf('PCP->CED time range: %.3f to %.3f s\n', ...
        t_PCP_inCED(1), t_PCP_inCED(end));
    fprintf('BHV landmark events: %d\n', numel(t_bhv_inCED));
    fprintf('BHV trial rows: %d\n', height(tbl_BHV));
    fprintf('Displayed overview: first %.0f seconds of CED time\n', plot_dur_s);

    disp('Available CED fields:');
    disp(fieldnames(ppCED_out));

    if ~isempty(tbl_BHV)
        disp('First BHV trial rows (all translated t_* columns use CED time):');
        disp(tbl_BHV(1:min(10, height(tbl_BHV)), :));
    end
end
%% If you wish to output the aligned LFP with CED to a MAT file
% Set this true to rerun alignment without creating diagnostic figures.
% Saved MAT-file variables:
%   lfp          - aligned LFP FieldTrip struct, with time in CED seconds
%   ced          - CED-side signals and events in their native CED clock
%   t_pcp_ced    - PCP/LFP sample times in CED seconds
%   t_bhv_ced    - BHV landmark times in CED seconds
%   bhv_trials   - BHV trial-event table with t_* columns in CED seconds
%   source_file  - full path of the input preproc MAT file
% Set fdir_output to an existing output directory, or leave it empty to save
% beside the source preproc file.
do_save_aligned = false;
fdir_output = "";

if do_save_aligned
    % Disable every plot panel. demo_apply_alignment will still load, align,
    % filter, and return data, but it will not create a MATLAB figure.
    [lfp, t_pcp_ced, t_bhv_ced, ced, bhv_trials] = demo_apply_alignment( ...
        fpath_preproc, 'msl', ...
        'show_LFP',          false, ...
        'show_LFP_PSD',      false, ...
        'show_voice',        false, ...
        'show_voice_raw',    false, ...
        'show_EMG',          false, ...
        'show_stimArt',      false, ...
        'show_KbEvents',     false, ...
        'show_syncLines',    false, ...
        'voice_env_channel', 'voice', ...
        'psd_fmax',          120, ...
        'plot_dur_s',        Inf);

    [source_dir, preproc_name] = fileparts(fpath_preproc);
    if strlength(fdir_output) == 0
        fdir_output = source_dir;
    end
    assert(exist(fdir_output, 'dir') == 7, ...
        'Output directory was not found:\n%s', fdir_output);

    fpath_aligned = fullfile(fdir_output, ...
        [preproc_name, '_desc-aligned-ephys.mat']);
    source_file = fpath_preproc;
    save(fpath_aligned, ...
        'lfp', 'ced', 't_pcp_ced', 't_bhv_ced', 'bhv_trials', ...
        'source_file', '-v7.3');
    fprintf('Saved aligned LFP and CED data:\n%s\n', fpath_aligned);
end


%% Lightweight fallback overview: no yyaxis or per-event graphics objects
% This figure is intentionally simple and should remain visible even if the
% richer diagnostic tiled layout has a MATLAB desktop-rendering issue.
t_lfp = lfp_aligned.time{1};
idx_lfp = t_lfp >= t_lfp(1) & t_lfp <= t_lfp(1) + plot_dur_s;
n_ch = numel(lfp_aligned.label);
fig_simple = figure( ...
    'Name', 'SML012 C201 aligned overview (fallback)', ...
    'Color', 'w', 'NumberTitle', 'off');
tl_simple = tiledlayout(fig_simple, n_ch + 1, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_simple, 'SML012 C201 MSL: aligned LFP and voice overview', ...
    'Interpreter', 'none');
ax_simple = gobjects(n_ch + 1, 1);
for i_ch = 1:n_ch
    ax_simple(i_ch) = nexttile;
    plot(ax_simple(i_ch), t_lfp(idx_lfp), ...
        lfp_aligned.trial{1}(i_ch, idx_lfp), 'k-', 'LineWidth', 0.4);
    grid(ax_simple(i_ch), 'on');
    ylabel(ax_simple(i_ch), lfp_aligned.label{i_ch}, 'Interpreter', 'none');
end
ax_simple(end) = nexttile;
if isfield(ppCED_out, 'voiceEnv') && ~isempty(ppCED_out.voiceEnv)
    t_voice = ppCED_out.voiceEnv.time{1};
    idx_voice = t_voice >= t_lfp(1) & t_voice <= t_lfp(1) + plot_dur_s;
    plot(ax_simple(end), t_voice(idx_voice), ...
        ppCED_out.voiceEnv.trial{1}(1, idx_voice), 'Color', [0.80 0.40 0.10]);
    ylabel(ax_simple(end), 'Voice env');
else
    text(ax_simple(end), 0.5, 0.5, 'No voice envelope available', ...
        'HorizontalAlignment', 'center');
end
grid(ax_simple(end), 'on');
xlabel(ax_simple(end), 'CED time (s)');
linkaxes(ax_simple, 'x');

% Useful workspace outputs:
%   lfp_aligned  - bandpass-filtered LFP with PCP times mapped to CED clock
%   t_PCP_inCED  - PCP sample time vector in CED time
%   t_bhv_inCED  - first available BHV landmark event vector in CED time
%   ppCED_out    - CED struct, including voiceEnv and available sync data
%   tbl_BHV      - all available BHV t_* event columns in CED time
%
% To inspect a later range, rerun the function with a shorter window and then
% use the figure's linked x-axes to navigate. Keep show_stimArt/show_voice_raw
% false initially; turn either on only after the base figure renders normally.
