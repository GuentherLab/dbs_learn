function [lfp_aligned, t_PCP_inCED, t_bhv_inCED, ppCED_out, tbl_BHV] = ...
        demo_apply_alignment(varargin)
% DEMO_APPLY_ALIGNMENT  Load, align, and inspect a BUSpeech preproc .mat.
%
%  Standalone diagnostic — no dependency on FileInfo excel, fileConfig, or
%  per-subject helper functions.  Provide the preproc file path directly.
%
%  Usage
%    demo_apply_alignment()
%    demo_apply_alignment(fpath_preproc)
%    demo_apply_alignment(fpath_preproc, task_name)
%    demo_apply_alignment(fpath_preproc, task_name, Name, Value)
%    [lfp, t_pcp, t_bhv, ppCED, tbl] = demo_apply_alignment(...)
%
%  Positional arguments
%    fpath_preproc  – full path to *_preproc*.mat   ('' → file-picker)
%    task_name      – task label string, controls default display panels:
%                       'ft'|'msl'|'ssl' → BHV trial event lines shown
%                       'shoc'|'titrate' → KbEvents prominent, no BHV lines
%                       '' (default)     → show everything available
%
%  Optional name-value pairs                              (default)
%    'show_LFP'           – LFP time-series (one tile per channel)   true
%    'show_LFP_PSD'       – PSD tile alongside each LFP row          true
%    'show_voice'         – Voice-envelope row                        true
%    'show_EMG'           – EMG row (silently skipped if absent)      true
%    'show_stimArt'       – stimArt overlay on LFP yyaxis right       true
%    'show_KbEvents'      – KbEvent vertical lines + code labels      true
%    'show_syncLines'     – evPCP / evTask / BHV trial event xlines   true
%    'show_voice_raw'     – Raw CED.voice trace on voice row (yyright)  true
%    'plot_dur_s'         – Seconds of data to display (Inf = all)    Inf
%    'voice_env_channel'  – Source channel: 'voice'|'mic'|'audStim' 'voice'
%    'psd_fmax'           – Upper freq limit for PSD plots (Hz)       150
%
%  Outputs
%    lfp_aligned   – FieldTrip struct, bandpass-filtered LFP in CED time
%    t_PCP_inCED   – [1×nSmp]  PCP sample times mapped to CED clock (s)
%    t_bhv_inCED   – [nTrials×1] key BHV landmark times in CED time ([] if none)
%    ppCED_out     – CED preprocessing struct (with .voiceEnv added)
%    tbl_BHV       – Trial event table with all t_* columns in CED time

    %% ── 0. Parse inputs ─────────────────────────────────────────────────
    p = inputParser();
    p.addOptional('fpath_preproc', '',       @(x) ischar(x) || isstring(x));
    p.addOptional('task_name',     '',       @(x) ischar(x) || isstring(x));
    p.addParameter('show_LFP',           true,    @islogical);
    p.addParameter('show_LFP_PSD',       true,    @islogical);
    p.addParameter('show_voice',         true,    @islogical);
    p.addParameter('show_EMG',           true,    @islogical);
    p.addParameter('show_stimArt',       true,    @islogical);
    p.addParameter('show_KbEvents',      true,    @islogical);
    p.addParameter('show_syncLines',     true,    @islogical);
    p.addParameter('show_voice_raw',     true,    @islogical);
    p.addParameter('plot_dur_s',         Inf,     @isnumeric);
    p.addParameter('voice_env_channel',  'voice', @(x) ischar(x) || isstring(x));
    p.addParameter('psd_fmax',           150,     @isnumeric);
    p.parse(varargin{:});
    r = p.Results;

    fpath_preproc = char(r.fpath_preproc);
    task_name     = lower(char(r.task_name));

    %% ── 1. File selection ───────────────────────────────────────────────
    if isempty(fpath_preproc) || exist(fpath_preproc, 'file') ~= 2
        [fn, fd] = uigetfile( ...
            {'*_preproc*.mat', 'Preproc files (*_preproc*.mat)'; ...
             '*.mat',          'All MAT files (*.mat)'}, ...
            'Select preproc .mat file');
        if isequal(fn, 0)
            error('demo_apply_alignment: no file selected — aborted.');
        end
        fpath_preproc = fullfile(fd, fn);
    end
    fprintf('Loading: %s\n', fpath_preproc);

    %% ── 2. Load, align, filter ──────────────────────────────────────────
    cfg_lad                      = [];
    cfg_lad.voice_env_channel    = r.voice_env_channel;
    cfg_lad.evpcp_reject_pad_sec = 1.0;

    [lfp_aligned, ~, tbl_BHV, ~, ~, ~, t_PCP_inCED, ~, ppCED_out] = ...
        load_and_align_data_local(fpath_preproc, cfg_lad);

    %% ── 3. Derive t_bhv_inCED ───────────────────────────────────────────
    % Return the first available landmark event column as a convenience vector.
    t_bhv_inCED = [];
    if ~isempty(tbl_BHV)
        for iC = 1:3
            landmark_cols = {'t_go_aud_on', 't_stim_vis_on', 't_stim_aud_on'};
            col = landmark_cols{iC};
            if ismember(col, tbl_BHV.Properties.VariableNames)
                t_bhv_inCED = tbl_BHV.(col);
                fprintf('t_bhv_inCED  ←  tbl_BHV.%s  (%d events)\n', ...
                    col, numel(t_bhv_inCED));
                break;
            end
        end
    end

    %% ── 4. Diagnostic plot ──────────────────────────────────────────────
    [~, fn_short] = fileparts(fpath_preproc);
    plot_sync_alignment(lfp_aligned, ppCED_out, tbl_BHV, task_name, fn_short, r);

    fprintf('\ndemo_apply_alignment: complete.\n');
    fprintf('  lfp_aligned   [%d ch × %d smp @ %.0f Hz, CED time]\n', ...
        numel(lfp_aligned.label), numel(lfp_aligned.time{1}), lfp_aligned.fsample);
    fprintf('  t_PCP_inCED   [1 × %d]\n',  numel(t_PCP_inCED));
    fprintf('  t_bhv_inCED   [%d × 1]\n',  numel(t_bhv_inCED));
end


% =========================================================================
function plot_sync_alignment(lfp, ppCED, tbl_BHV, task_name, fn_label, r)
% PLOT_SYNC_ALIGNMENT  Multi-panel diagnostic figure.
%
%  Layout (tiledlayout):
%    If show_LFP_PSD = true  → 4-column grid
%       • LFP channel rows : cols 1-3 = time series, col 4 = PSD
%       • Voice / EMG rows : span all 4 columns
%    If show_LFP_PSD = false → 1-column grid, all rows equal width
%
%  Colour legend
%    dark-grey dotted   evPCP sync pulses
%    magenta            evTask sync pulses
%    red                t_stim_vis_on   (BHV trial events)
%    green dashed       t_go_aud_on
%    lime dotted        t_go_aud_off
%    pink dashed        KbEvent codes

    % ── Availability flags ───────────────────────────────────────────────
    no_trial_tasks = {'shoc', 'titrate', 'tell'};
    is_shoc        = strcmp(task_name, 'shoc');
    is_titrate     = strcmp(task_name, 'titrate');
    show_bhv_lines = r.show_syncLines && ~isempty(tbl_BHV) && ...
                     ~any(strcmp(task_name, no_trial_tasks));

    has_voice     = r.show_voice     && isfield(ppCED,'voiceEnv') && ~isempty(ppCED.voiceEnv);
    has_voice_raw = r.show_voice_raw && isfield(ppCED,'voice')    && ~isempty(ppCED.voice);
    has_emg       = r.show_EMG       && isfield(ppCED,'EMG')      && ~isempty(ppCED.EMG);
    has_stim      = r.show_stimArt   && isfield(ppCED,'stimArt')  && ~isempty(ppCED.stimArt);
    has_kb        = r.show_KbEvents  && isfield(ppCED,'KbEvent')  && ...
                    istable(ppCED.KbEvent) && height(ppCED.KbEvent) > 0;
    has_psd       = r.show_LFP && r.show_LFP_PSD;
    has_stimAmp   = is_titrate && isfield(ppCED,'stimSettings') && ...
                    isfield(ppCED.stimSettings,'stimAmpL') && ...
                    ~isempty(ppCED.stimSettings.stimAmpL);

    nCh   = numel(lfp.label);
    nRows = r.show_LFP * nCh + has_voice + has_emg + has_stimAmp;
    if nRows == 0
        warning('plot_sync_alignment: no panels selected — nothing to show.');
        return;
    end

    % 1 column (no PSD) or 4 columns (LFP ts + PSD, voice/EMG span all)
    nCols = 1 + 3 * has_psd;

    % ── Time window ──────────────────────────────────────────────────────
    t_lfp = lfp.time{1};
    if isfinite(r.plot_dur_s)
        xl = [t_lfp(1), t_lfp(1) + r.plot_dur_s];
    else
        xl = [t_lfp(1), t_lfp(end)];
    end

    % ── Figure / layout ──────────────────────────────────────────────────
    fig = figure('Name', sprintf('Alignment | %s | %s', upper(task_name), fn_label), ...
        'Color','w', 'NumberTitle','off', ...
        'Units','normalized', 'Position',[0.02 0.02 0.96 0.94]);
    tl = tiledlayout(fig, nRows, nCols, 'TileSpacing','tight', 'Padding','compact');
    title(tl, sprintf('Alignment diagnostic  —  task: %s\n%s', ...
        upper(task_name), fn_label), 'Interpreter','none');

    ax_linked = gobjects(0);   % time-series axes to linkaxes

    % ── LFP rows ─────────────────────────────────────────────────────────
    if r.show_LFP
        idx_l  = t_lfp >= xl(1) & t_lfp <= xl(2);
        fs_lfp = lfp.fsample;

        for iCh = 1:nCh
            sig = lfp.trial{1}(iCh, :);

            % Time-series tile (spans 3 of 4 cols when PSD enabled)
            if has_psd
                ax = nexttile([1, 3]);
            else
                ax = nexttile;
            end
            hold(ax,'on');  grid(ax,'on');
            plot(ax, t_lfp(idx_l), sig(idx_l), ...
                'Color',[0.18 0.40 0.80], 'LineWidth',0.5);
            ylabel(ax, 'LFP (\muV)');
            xlabel(ax, 'CED time (s)');
            title(ax, lfp.label{iCh}, 'Interpreter','none');
            xlim(ax, xl);

            % stimArt on yyaxis right (z-scored for visual scale)
            if has_stim
                yyaxis(ax, 'right');
                t_sa  = ppCED.stimArt.time{1};
                s_sa  = ppCED.stimArt.trial{1}(1,:);
                idx_s = t_sa >= xl(1) & t_sa <= xl(2);
                plot(ax, t_sa(idx_s), zscore(s_sa(idx_s)), ...
                    'Color',[0.50 0.50 0.50 0.55], 'LineWidth',0.35);
                ylabel(ax, 'stimArt (z)');
                ylim(ax, [-8, 8]);
                yyaxis(ax, 'left');
            end

            local_add_event_lines(ax, xl, ppCED, tbl_BHV, ...
                show_bhv_lines, has_kb, r.show_syncLines, is_shoc);
            ax_linked(end+1) = ax; %#ok<AGROW>

            % PSD tile (4th column)
            if has_psd
                ax_p = nexttile;
                hold(ax_p,'on');  grid(ax_p,'on');
                sig_w = sig(idx_l);
                sig_w(isnan(sig_w)) = 0;
                nfft = min(2^nextpow2(round(fs_lfp * 2)), numel(sig_w));
                [pxx, f_ax] = pwelch(sig_w, nfft, round(nfft/2), nfft, fs_lfp);
                fmax_p = min(r.psd_fmax, fs_lfp / 2);
                fk     = f_ax <= fmax_p & f_ax > 0;
                semilogy(ax_p, f_ax(fk), pxx(fk), ...
                    'Color',[0.18 0.40 0.80], 'LineWidth',1.0);
                xlabel(ax_p, 'Freq (Hz)');
                ylabel(ax_p, 'PSD (\muV²/Hz)');
                title(ax_p, sprintf('%s  PSD', lfp.label{iCh}), 'Interpreter','none');
                xlim(ax_p, [0, fmax_p]);
            end
        end  % iCh
    end  % show_LFP

    % ── Voice envelope row ───────────────────────────────────────────────
    if has_voice
        ax = nexttile([1, nCols]);
        hold(ax,'on');  grid(ax,'on');
        t_v   = ppCED.voiceEnv.time{1};
        s_v   = ppCED.voiceEnv.trial{1}(1,:);
        idx_v = t_v >= xl(1) & t_v <= xl(2);
        % raw CED.voice on yyaxis right (z-scored, semi-transparent)
        if has_voice_raw
            yyaxis(ax, 'right');
            t_vr   = ppCED.voice.time{1};
            s_vr   = ppCED.voice.trial{1}(1,:);
            idx_vr = t_vr >= xl(1) & t_vr <= xl(2);
            plot(ax, t_vr(idx_vr), zscore(s_vr(idx_vr)), ...
                'Color',[0.25 0.65 0.90 0.45], 'LineWidth',0.3);
            ylabel(ax, 'Raw voice (z)');
            ylim(ax, [-8, 8]);
            yyaxis(ax, 'left');
        end
        plot(ax, t_v(idx_v), s_v(idx_v), ...
            'Color',[0.80 0.40 0.10], 'LineWidth',0.8);
        ylabel(ax, 'Voice Env (a.u.)');
        xlabel(ax, 'CED time (s)');
        title(ax, 'Voice  (envelope + raw)');
        xlim(ax, xl);
        local_add_event_lines(ax, xl, ppCED, tbl_BHV, ...
            show_bhv_lines, has_kb, r.show_syncLines, is_shoc);
        ax_linked(end+1) = ax; %#ok<AGROW>
    end

    % ── EMG row ──────────────────────────────────────────────────────────
    if has_emg
        ax = nexttile([1, nCols]);
        hold(ax,'on');  grid(ax,'on');
        emg   = ppCED.EMG;
        t_e   = emg.time{1};
        idx_e = t_e >= xl(1) & t_e <= xl(2);
        pal_e = lines(max(size(emg.trial{1}, 1), 1));
        for iE = 1:size(emg.trial{1}, 1)
            plot(ax, t_e(idx_e), emg.trial{1}(iE, idx_e), ...
                'Color', pal_e(iE,:), 'LineWidth',0.5);
        end
        if ~isempty(emg.label)
            legend(ax, emg.label, 'Interpreter','none', ...
                'Location','northwest', 'FontSize',7);
        end
        ylabel(ax, 'EMG (a.u.)');
        xlabel(ax, 'CED time (s)');
        title(ax, 'EMG');
        xlim(ax, xl);
        local_add_event_lines(ax, xl, ppCED, tbl_BHV, ...
            show_bhv_lines, has_kb, r.show_syncLines, is_shoc);
        ax_linked(end+1) = ax; %#ok<AGROW>
    end

    % ── Stim amplitude row (Titrate only) ──────────────────────────────────
    if has_stimAmp
        ax = nexttile([1, nCols]);
        hold(ax,'on');  grid(ax,'on');
        ss    = ppCED.stimSettings;
        t_ss  = ss.t(:)';
        idx_ss = t_ss >= xl(1) & t_ss <= xl(2);
        ampL  = double(ss.stimAmpL(:)');  ampL(~idx_ss) = [];
        ampR  = double(ss.stimAmpR(:)');  ampR(~idx_ss) = [];
        t_ss_w = t_ss(idx_ss);
        % plot as step-style so each sample value is visible
        stairs(ax, t_ss_w, ampL, 'Color',[0.10 0.45 0.85], 'LineWidth',1.2);
        stairs(ax, t_ss_w, ampR, 'Color',[0.85 0.25 0.10], 'LineWidth',1.2);
        legend(ax, {'stimAmpL','stimAmpR'}, 'Interpreter','none', ...
            'Location','northwest', 'FontSize',7);
        ylabel(ax, 'Stim Amp (mA)');
        xlabel(ax, 'CED time (s)');
        title(ax, 'Stimulation Amplitude  (L = blue, R = red)', 'Interpreter','none');
        xlim(ax, xl);
        local_add_event_lines(ax, xl, ppCED, tbl_BHV, ...
            show_bhv_lines, has_kb, r.show_syncLines, is_shoc);
        ax_linked(end+1) = ax; %#ok<AGROW>
    end

    % ── Link all time-series axes ─────────────────────────────────────────
    if numel(ax_linked) > 1
        linkaxes(ax_linked, 'x');
    end
end


% =========================================================================
function local_add_event_lines(ax, xl, ppCED, tbl_BHV, ...
        show_bhv_lines, has_kb, show_sync_lines, is_shoc)
% LOCAL_ADD_EVENT_LINES  Draw sync xlines and KbEvent labels on axes ax.
%
%  Draws onto whatever axes ax points to; called for each time-series panel.

    if nargin < 8, is_shoc = false; end

    % ── SHOC segment shading (drawn first so lines render on top) ────────
    if is_shoc && has_kb
        local_shade_shoc_segments(ax, xl, ppCED.KbEvent);
    end

    COL_PCP  = [0.25 0.25 0.25];   % dark grey : evPCP sync pulses
    COL_TSK  = [0.70 0.10 0.80];   % magenta   : evTask sync pulses
    COL_VIS  = [0.85 0.10 0.10];   % red       : t_stim_vis_on
    COL_GAUD = [0.10 0.65 0.10];   % green     : t_go_aud_on
    COL_GOFF = [0.50 0.80 0.15];   % lime      : t_go_aud_off
    COL_KB   = [0.85 0.30 0.65];   % pink      : KbEvent

    % ── Sync xlines ──────────────────────────────────────────────────────
    if show_sync_lines
        if isfield(ppCED,'evPCP') && isfield(ppCED.evPCP,'sync_t') && ...
                ~isempty(ppCED.evPCP.sync_t)
            t_s = ppCED.evPCP.sync_t;
            t_s = t_s(isfinite(t_s) & t_s >= xl(1) & t_s <= xl(2));
            if ~isempty(t_s)
                xline(ax, t_s, ':', 'Color',COL_PCP, 'LineWidth',0.7, 'Alpha',0.5);
            end
        end
        if isfield(ppCED,'evTask') && isfield(ppCED.evTask,'sync_t') && ...
                ~isempty(ppCED.evTask.sync_t)
            t_s = ppCED.evTask.sync_t;
            t_s = t_s(isfinite(t_s) & t_s >= xl(1) & t_s <= xl(2));
            if ~isempty(t_s)
                xline(ax, t_s, '-', 'Color',COL_TSK, 'LineWidth',0.8, 'Alpha',0.35);
            end
        end
    end

    % ── BHV trial event xlines ───────────────────────────────────────────
    if show_bhv_lines && ~isempty(tbl_BHV)
        ev_defs = { ...
            't_stim_vis_on', COL_VIS,  '-',  0.40; ...
            't_go_aud_on',   COL_GAUD, '--', 0.40; ...
            't_go_aud_off',  COL_GOFF, ':',  0.40 };
        for iE = 1:size(ev_defs,1)
            col_name = ev_defs{iE,1};
            if ~ismember(col_name, tbl_BHV.Properties.VariableNames), continue; end
            ev_t = tbl_BHV.(col_name);
            ev_t = ev_t(isfinite(ev_t) & ev_t >= xl(1) & ev_t <= xl(2));
            if isempty(ev_t), continue; end
            xline(ax, ev_t, ev_defs{iE,3}, ...
                'Color', ev_defs{iE,2}, 'LineWidth',0.6, 'Alpha',ev_defs{iE,4});
        end
    end

    % ── KbEvent vertical lines + rotated code labels ──────────────────────
    if has_kb
        KbEvt  = ppCED.KbEvent;
        yl_now = ylim(ax);
        y_text = yl_now(2) - 0.03 * diff(yl_now);   % just below top edge
        for jj = 1:height(KbEvt)
            t_ev = KbEvt.Time_sec(jj);
            if ~isfinite(t_ev) || t_ev < xl(1) || t_ev > xl(2), continue; end
            xline(ax, t_ev, '--', 'Color',COL_KB, 'LineWidth',1.2);
            code_v = KbEvt.Code(jj);
            if iscell(code_v),    code_v = code_v{1};    end
            if isstring(code_v),  code_v = char(code_v); end
            text(ax, t_ev, y_text, code_v, ...
                'Color', COL_KB, 'FontSize',7, ...
                'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                'Rotation',75, 'Clipping','on');
        end
    end
end


% =========================================================================
function [lfp, voiceEnv, tbl_BHV, fpass, cfg_filt, m_Align, t_PCP_inCED, KbEvt, ppCED] = ...
        load_and_align_data_local(fpath_pp, cfg)
% LOAD_AND_ALIGN_DATA_LOCAL  Load a preproc .mat, apply time-alignment,
%   filter LFP, build voiceEnv, translate BHV event times to CED.
%
%   Local self-contained copy of load_and_align_data() from
%   qwikCheck_alignment_and_evokedActivity.m — no inline plotting.
%   7th output is t_PCP_inCED (PCP sample times in CED clock) rather than
%   the raw CED time vector.
%
%  Outputs
%   lfp           – FieldTrip continuous struct (bandpass-filtered, CED time)
%   voiceEnv      – FieldTrip continuous struct (voice envelope, CED time)
%   tbl_BHV       – trial event table with all t_* cols in CED time
%   fpass         – [flo fhi] bandpass limits used (Hz)
%   cfg_filt      – ft_preprocessing config (reusable for baseline)
%   m_Align       – alignment struct from the .mat file
%   t_PCP_inCED   – [1×nSmp] PCP sample times mapped to CED clock (s)
%   KbEvt         – keyboard event table (Time_sec, Code) in CED time
%   ppCED         – full CED preprocessing struct (with .voiceEnv added)

    voice_env_channel    = ft_getopt(cfg, 'voice_env_channel',    'voice');
    evpcp_reject_pad_sec = ft_getopt(cfg, 'evpcp_reject_pad_sec', 1.0);

    %% Load
    S = load(fpath_pp);
    S = ft_checkconfig(S, 'renamed', {'preprocEMG', 'CED'});
    S = ft_checkconfig(S, 'renamed', {'preprocLFP', 'PCP'});
    ppCED   = S.CED;
    ppPCP   = S.PCP;
    m_Align = S.m_Align;

    pcp   = ppPCP.data;
    t_PCP = pcp.time{1};

    %% Sync / translation handles
    try
        t_tskSync_BHV = m_Align.t_tskSync_BHV;
    catch
        t_tskSync_BHV = [];
    end

    fcn_translateTime = m_Align.fcn_translateTime;
    tTrans_PCP_to_CED = m_Align.TimeTrans_PCP_to_CED_byEvent;
    tTrans_BHV_to_CED = m_Align.TimeTrans_BHV_to_CED;

    %% PCP → CED time
    if ~isempty(tTrans_PCP_to_CED)
        t_PCP_inCED = fcn_translateTime(t_PCP, tTrans_PCP_to_CED);
    else
        t_PCP_inCED = t_PCP;
        fprintf('load_and_align_data_local: no PCP→CED transform — using raw PCP time.\n');
    end
    pcp.oldtime = pcp.time;
    pcp.time{1} = t_PCP_inCED;

    if ~isempty(tTrans_BHV_to_CED) && ~isempty(t_tskSync_BHV)
        fcn_translateTime(t_tskSync_BHV, tTrans_BHV_to_CED); % side-effect: validates transform
    end

    %% Mask evPCP windows in CED channels
    evpcp_sync_t = local_get_evpcp_sync_times(ppCED, m_Align);
    if ~isempty(evpcp_sync_t) && evpcp_reject_pad_sec > 0
        for iF = 1:4
            fn_list = {'voice','mic','audStim','EMG'};
            fn = fn_list{iF};
            if isfield(ppCED, fn) && ~isempty(ppCED.(fn))
                ppCED.(fn) = mask_artifact_times(ppCED.(fn), evpcp_sync_t, evpcp_reject_pad_sec);
            end
        end
    end

    %% Translate BHV trial-table event times → CED
    if isfield(S,'bhv') && isfield(S.bhv,'evTrial') && ~isempty(S.bhv.evTrial)
        tbl_BHV = S.bhv.evTrial;
        if ismember('t_go_aud_on', tbl_BHV.Properties.VariableNames)
            tbl_BHV = rmmissing(tbl_BHV, 'DataVariables', {'t_go_aud_on'});
        end
        timeCols = tbl_BHV.Properties.VariableNames( ...
            startsWith(tbl_BHV.Properties.VariableNames, 't_'));

        for iCol = 1:numel(timeCols)
            col = timeCols{iCol};
            if ~ismember(col, tbl_BHV.Properties.VariableNames), continue; end

            if ismember('timingRef', tbl_BHV.Properties.VariableNames)
                refs = tbl_BHV.timingRef;
            else
                refs = repmat({'bhvFileStart'}, height(tbl_BHV), 1);
            end
            need_transform = ~all(ismember(refs, {'CEDFileStart'}));

            if all(ismember(refs, {'bhvFileStart'}))
                if isfield(m_Align,'TimeTrans_TRL_to_CED')
                    a = fcn_translateTime(tbl_BHV.(col)(:), m_Align.TimeTrans_TRL_to_CED)';
                    need_transform = false;
                else
                    error('load_and_align_data_local: TimeTrans_TRL_to_CED missing.');
                end
            elseif all(ismember(refs, {'fristPCP'}))
                a = tbl_BHV.(col)(:)';
            elseif all(ismember(refs, {'CEDFileStart'}))
                a = tbl_BHV.(col)';
            else
                error('load_and_align_data_local: inconsistent timingRef values.');
            end

            if need_transform
                a = fcn_translateTime(a, tTrans_BHV_to_CED);
            end
            tbl_BHV.(col) = a(:);
        end
        tbl_BHV.timingRef(:) = {'CEDFileStart'};
    else
        tbl_BHV = table();
    end

    %% Filter LFP
    fprintf('Filtering LFP …\n');
    fpass          = [0.5, min(140, pcp.fsample / 2 * 0.8)];
    cfg_filt           = [];
    cfg_filt.bsfilter  = 'no';
    cfg_filt.bpfilter  = 'yes';
    cfg_filt.bpfreq    = fpass;
    lfp = ft_preprocessing(cfg_filt, pcp);

    %% Voice envelope (downsampled to CED fs)
    t_CED       = ppCED.stimArt.time{1};
    fs_voice_DS = ppCED.stimArt.fsample;

    switch lower(voice_env_channel)
        case 'mic',     src_chan = ppCED.mic;
        case 'audstim', src_chan = ppCED.audStim;
        otherwise,      src_chan = ppCED.voice;   % default
    end

    factor_DS = src_chan.fsample / fs_voice_DS;
    % Replace NaN with 0 before envelope (NaN in RMS window → silent artifact)
    src_chan.trial{1}(isnan(src_chan.trial{1})) = 0;
    raw_voice  = src_chan.trial{1};
    envVoice   = downsample(envelope(raw_voice, 2000, 'rms'), factor_DS);
    t_voice    = downsample(src_chan.time{1}, factor_DS);

    % Trim/pad by ≤2 samples (1-sample rounding tolerance at boundary)
    nCED = numel(t_CED);
    nV   = numel(t_voice);
    if nV ~= nCED
        if abs(nV - nCED) <= 2
            nUse     = min(nV, nCED);
            t_voice  = t_voice(1:nUse);
            envVoice = envVoice(1:nUse);
            while numel(t_voice) < nCED
                dt       = t_voice(end) - t_voice(end-1);
                t_voice  = [t_voice,  t_voice(end)  + dt]; %#ok<AGROW>
                envVoice = [envVoice, envVoice(end)];       %#ok<AGROW>
            end
        else
            error('Voice envelope length (%d) ≠ CED length (%d)', nV, nCED);
        end
    end

    voiceEnv         = [];
    voiceEnv.fsample = fs_voice_DS;
    voiceEnv.time    = {t_voice};
    voiceEnv.trial   = {envVoice};
    voiceEnv.label   = {'voiceEnvelope'};

    ppCED.voiceEnv = voiceEnv;

    %% Extract stim amplitude from PCP.data.stimSettings (for Titrate)
    try
        ss = ppPCP.data.stimSettings;
        ppCED.stimSettings.stimAmpL = ss.stimAmpL;  % [1×nSmp] mA or steps
        ppCED.stimSettings.stimAmpR = ss.stimAmpR;
        ppCED.stimSettings.t        = t_PCP_inCED;  % aligned CED time
    catch
        % field absent — leave ppCED.stimSettings unset
    end

    %% KbEvt output
    KbEvt = ft_getopt(ppCED, 'KbEvent', table());
end


% =========================================================================
function lfp_out = mask_artifact_times(lfp_in, t_artifacts, pad_sec)
% MASK_ARTIFACT_TIMES  Set ±pad_sec samples around each artifact time to NaN.

    lfp_out = lfp_in;
    t       = lfp_in.time{1};
    nMasked = 0;

    for iTrl = 1:numel(lfp_out.trial)
        mask = false(1, numel(t));
        for iA = 1:numel(t_artifacts)
            mask = mask | (t >= (t_artifacts(iA) - pad_sec) & ...
                           t <= (t_artifacts(iA) + pad_sec));
        end
        lfp_out.trial{iTrl}(:, mask) = NaN;
        nMasked = nMasked + sum(mask);
    end

    fprintf('mask_artifact_times: %d artifact(s) × ±%.2f s  →  %d samples NaN\n', ...
        numel(t_artifacts), pad_sec, nMasked);
end


% =========================================================================
function evpcp_sync_t = local_get_evpcp_sync_times(CED, m_Align)
% LOCAL_GET_EVPCP_SYNC_TIMES  Prefer CED.evPCP.sync_t; fall back to m_Align.

    evpcp_sync_t = [];
    if nargin >= 1 && isstruct(CED) && isfield(CED,'evPCP') && ...
            isstruct(CED.evPCP) && isfield(CED.evPCP,'sync_t') && ...
            ~isempty(CED.evPCP.sync_t)
        evpcp_sync_t = double(CED.evPCP.sync_t(:));
    elseif nargin >= 2 && isstruct(m_Align) && ...
            isfield(m_Align,'t_pcpSync_CED') && ~isempty(m_Align.t_pcpSync_CED)
        evpcp_sync_t = double(m_Align.t_pcpSync_CED(:));
    end

    evpcp_sync_t = evpcp_sync_t(isfinite(evpcp_sync_t));
    evpcp_sync_t = unique(evpcp_sync_t(:)');
end


% =========================================================================
function local_shade_shoc_segments(ax, xl, KbEvt)
% LOCAL_SHADE_SHOC_SEGMENTS  Draw semi-transparent filled rectangles for
%   SHOC segments defined by consecutive on/off pairs of KbEvent codes.
%
%   SHOC codes (case-insensitive): b, h, r, d
%   Pairing rule: within each code, sort events by time and pair
%   consecutive rows as (onset=row 2k-1, offset=row 2k).
%
%   Each code gets a distinct near-transparent colour; labels appear at
%   the top of the patch.

    if isempty(KbEvt) || height(KbEvt) == 0, return; end

    % Code → colour map (RGBA-like: store RGB + FaceAlpha separately)
    shoc_codes = {'b',  'h',  'r',  'd'};
    shoc_cols  = [0.25 0.50 0.90;   % blue
                  0.95 0.55 0.10;   % orange
                  0.85 0.15 0.15;   % red
                  0.15 0.70 0.30];  % green
    face_alpha = 0.12;   % nearly transparent
    edge_alpha = 0.30;

    % Current y-limits (patches will span full y range)
    yl = ylim(ax);

    % Normalise code column to char cell array
    raw_codes = KbEvt.Code;
    if isstring(raw_codes),  raw_codes = cellstr(raw_codes); end
    if isnumeric(raw_codes), raw_codes = cellstr(num2str(raw_codes(:))); end
    % Strip whitespace
    raw_codes = strtrim(raw_codes);

    for iC = 1:numel(shoc_codes)
        code_str = shoc_codes{iC};
        col_rgb  = shoc_cols(iC,:);

        % Find all KbEvt rows whose code matches (case-insensitive, exact)
        mask = strcmpi(raw_codes, code_str);
        if ~any(mask), continue; end

        t_code = KbEvt.Time_sec(mask);
        t_code = sort(t_code(isfinite(t_code)));

        % Pair consecutive rows: (1,2), (3,4), …
        n_pairs = floor(numel(t_code) / 2);
        if n_pairs == 0, continue; end

        for iP = 1:n_pairs
            t_on  = t_code(2*iP - 1);
            t_off = t_code(2*iP);

            % Clip to visible window
            t_on  = max(t_on,  xl(1));
            t_off = min(t_off, xl(2));
            if t_on >= t_off, continue; end

            % Patch corners: [x1 x1 x2 x2], [y_lo y_hi y_hi y_lo]
            patch(ax, ...
                [t_on, t_on, t_off, t_off], ...
                [yl(1), yl(2), yl(2), yl(1)], ...
                col_rgb, ...
                'FaceAlpha', face_alpha, ...
                'EdgeColor', col_rgb, ...
                'EdgeAlpha', edge_alpha, ...
                'LineWidth', 0.8, ...
                'HitTest',   'off');

            % Label at top of patch
            text(ax, (t_on + t_off)/2, yl(2), upper(code_str), ...
                'Color',              col_rgb, ...
                'FontSize',           8, ...
                'FontWeight',         'bold', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment',  'top', ...
                'Clipping',           'on');
        end
    end
end
