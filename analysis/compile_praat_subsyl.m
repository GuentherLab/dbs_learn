% =========================================================================
%  Load Praat TextGrid phonetic annotation data into a MATLAB trial table.
%
%  Dependencies:
%    - mpraat toolbox (https://github.com/bbTomas/mPraat)
%    - MATLAB R2016b or later
% =========================================================================

function compile_praat_subsyl(op)

vardefault('op',struct);
op.ses = 'subsyl';

field_default('op','sub','sml003');
field_default('op','task','test1');
field_default('op','chan','headphone');

paths = setpaths_dbs_learn(op);
lndmks = readtable(paths.landmarks_file,'FileType','text','Delimiter','tab');
op.run = lndmks.run(find(contains(lndmks.filename,op.task),1));
paths = setpaths_dbs_learn(op);
paths.audiofiles_table_filename = [paths.trial_audio,filesep, paths.filestr, 'audiofiles-',op.chan, '.tsv'];
paths.trial_audio_task_chan = [paths.trial_audio, filesep, 'ses-',op.ses, '_task-',op.task,'_run-',num2str(op.run),'_',op.chan];

%% ---- CONFIGURATION ------------------------------------------------------
% Path to Praat executable -- opened automatically on consistency errors.
%   Windows: 'C:\Program Files\Praat\Praat.exe'
%   macOS:   '/Applications/Praat.app/Contents/MacOS/Praat'
%   Linux:   '/usr/bin/praat'
praatExe = [paths.code_dbs_learn, filesep, 'praat', filesep, 'praat.exe'];   % copy of praat in the repo, so we don't have to find the user's normal version
% -------------------------------------------------------------------------

%% Constants
EXPECTED_TIERS = { ...
    'cons.1.accuracy',    'cons.1.error.type', 'vow.accuracy',      ...
    'cons.2.accuracy',    'cons.2.error.type', 'transcript',        ...
    'disfluency',         'comments',          'unusable.trial',    ...
    'difficult.to.score', 'speech.epoch',      'vowel.epoch',       ...
    'nontarget.sounds.epoch'};
N_TIERS         = 13;
VALID_ERR_TYPES = [1 3 4 5 6 7 8];

% TextGrid-derived column names (used for column ordering below)
TG_VARS = {'cons1_accuracy','cons1_error_type','vow_accuracy','cons2_accuracy', ...
           'cons2_error_type','transcript','disfluency','comments', ...
           'unusable_trial','difficult_to_score','speech_onset','speech_offset', ...
           'vowel_onset','vowel_offset','nontarget_sounds_epoch'};

if ~isempty(paths.trial_audio_task_chan) && paths.trial_audio_task_chan(end) ~= filesep
    paths.trial_audio_task_chan = [paths.trial_audio_task_chan, filesep];
end

%% Load trial table --------------------------------------------------------
audiofiles = readtable(paths.audiofiles_table_filename,'FileType','text','Delimiter','tab');
fprintf('Loading trial table: %s\n', paths.audiofiles_table_filename);
nTrials = height(audiofiles);
fprintf('  %d trials found.\n\n', nTrials);

assert(ismember('filename', audiofiles.Properties.VariableNames), ...
    'Trial table must contain a column named "filename".');
assert(ismember('name', audiofiles.Properties.VariableNames), ...
    'Trial table must contain a column named "name".');

if iscell(audiofiles.filename), audiofiles.filename = string(audiofiles.filename); end
if iscell(audiofiles.name),     audiofiles.name     = string(audiofiles.name);     end

%% Build trials table from audiofiles -------------------------------------
trials = audiofiles;

% Remove unwanted columns if present
for rv = {'postproc_gain','adjust_start','adjust_end'}
    if ismember(rv{1}, trials.Properties.VariableNames)
        trials.(rv{1}) = [];
    end
end

% Add trialnum column extracted from filename (e.g. "trial-5_..." -> 5)
if ~ismember('trialnum', trials.Properties.VariableNames)
    trialNums = NaN(nTrials, 1);
    for i = 1:nTrials
        tok = regexp(char(trials.filename(i)), 'trial-(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            trialNums(i) = str2double(tok{1});
        end
    end
    trials.trialnum = trialNums;
end

%% Initialise textgrid annotation columns ----------------------------------
trials.cons1_accuracy         = NaN(nTrials, 1);
trials.cons1_error_type       = cell(nTrials, 1);
trials.vow_accuracy           = NaN(nTrials, 1);
trials.cons2_accuracy         = NaN(nTrials, 1);
trials.cons2_error_type       = cell(nTrials, 1);
trials.transcript             = repmat("", nTrials, 1);
trials.disfluency             = repmat("", nTrials, 1);
trials.comments               = repmat("", nTrials, 1);
trials.unusable_trial         = NaN(nTrials, 1);   % double: NaN=unscored, 0=usable, 1=unusable
trials.difficult_to_score     = false(nTrials, 1);
trials.speech_onset           = NaN(nTrials, 1);
trials.speech_offset          = NaN(nTrials, 1);
trials.vowel_onset            = NaN(nTrials, 1);
trials.vowel_offset           = NaN(nTrials, 1);
trials.nontarget_sounds_epoch = cell(nTrials, 1);

for i = 1:nTrials
    trials.cons1_error_type{i}       = [];
    trials.cons2_error_type{i}       = [];
    trials.nontarget_sounds_epoch{i} = zeros(2, 0);
end

%% Reorder columns: trialnum, stim_group, is_native, name, [tg], [rest] ---
FIXED_FIRST = {'trialnum','stim_group','is_native','name'};
allVars     = trials.Properties.VariableNames;

% Collect remaining original vars in their original order
remainingVars = {};
for v = allVars
    if ~ismember(v{1}, FIXED_FIRST) && ~ismember(v{1}, TG_VARS)
        remainingVars{end+1} = v{1}; %#ok<AGROW>
    end
end

desiredOrder = [FIXED_FIRST, TG_VARS, remainingVars];
finalOrder   = desiredOrder(ismember(desiredOrder, allVars));  % skip any that don't exist
trials       = trials(:, finalOrder);

%% Main loop ---------------------------------------------------------------
nMissing = 0;

for iTrial = 1:nTrials

    wavName  = char(audiofiles.filename(iTrial));
    stimName = char(audiofiles.name(iTrial));
    [~, baseName, ~] = fileparts(wavName);
    tgPath  = [paths.trial_audio_task_chan, baseName, '.TextGrid'];
    wavPath = [paths.trial_audio_task_chan, baseName, '.wav'];

    fprintf('Trial %d/%d: %s\n', iTrial, nTrials, baseName);

    % Skip trials that have no audio file (audfile_start = NaN means the
    % trial entry exists in the table but no clip was ever extracted)
    if isnan(audiofiles.audfile_start(iTrial))
        fprintf('  Skipping: no audio file (audfile_start = NaN).\n');
        continue;
    end

    % Skip if TextGrid has not been scored yet
    if ~exist(tgPath, 'file')
        warning('Trial %d: TextGrid not found -- %s. Skipping.', ...
                iTrial, tgPath);
        nMissing = nMissing + 1;
        continue;
    end

    % Stim-name sanity check
    if ~contains(baseName, stimName)
        error('Trial %d: stim "%s" expected in filename "%s" but not found.', ...
              iTrial, stimName, baseName);
    end

    % Load TextGrid.
    % Suppress "encoding 'UTF-16' not supported" warning: Praat saves
    % TextGrids as UTF-16 when they contain IPA/Unicode text; tgRead with
    % 'auto' handles this correctly but emits the warning on its first
    % fopen attempt regardless.
    prevWarn = warning('off', 'MATLAB:iofun:UnsupportedEncoding'); %%%%%%%%%% not working - doesn't actually suppress
    tg       = tgRead(tgPath, 'auto');
    warning(prevWarn);

    % Validate tier count and names
    nFound = tgGetNumberOfTiers(tg);
    if nFound ~= N_TIERS
        error('Trial %d (%s): expected %d tiers, found %d.', ...
              iTrial, baseName, N_TIERS, nFound);
    end
    for k = 1:N_TIERS
        gotName = tgGetTierName(tg, k);
        if ~strcmp(gotName, EXPECTED_TIERS{k})
            error('Trial %d (%s): tier %d expected "%s", got "%s".', ...
                  iTrial, baseName, k, EXPECTED_TIERS{k}, gotName);
        end
    end

    % --- Tier 1: cons.1.accuracy  (0, 1, or 2) --------------------------
    lbl = getSingleLabel(tg, 1, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1 2])
            error('Trial %d (%s), cons.1.accuracy: expected 0/1/2, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        trials.cons1_accuracy(iTrial) = val;
    end

    % --- Tier 2: cons.1.error.type  --------------------------------------
    lbl = getSingleLabel(tg, 2, iTrial, baseName);
    if ~isempty(lbl)
        trials.cons1_error_type{iTrial} = parseErrTypes(lbl, iTrial, baseName, ...
                                              'cons.1.error.type', VALID_ERR_TYPES);
    end

    % --- Tier 3: vow.accuracy  (0 or 1) ----------------------------------
    lbl = getSingleLabel(tg, 3, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), vow.accuracy: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        trials.vow_accuracy(iTrial) = val;
    end

    % --- Tier 4: cons.2.accuracy  (0, 1, or 2) --------------------------
    lbl = getSingleLabel(tg, 4, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1 2])
            error('Trial %d (%s), cons.2.accuracy: expected 0/1/2, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        trials.cons2_accuracy(iTrial) = val;
    end

    % --- Tier 5: cons.2.error.type  --------------------------------------
    lbl = getSingleLabel(tg, 5, iTrial, baseName);
    if ~isempty(lbl)
        trials.cons2_error_type{iTrial} = parseErrTypes(lbl, iTrial, baseName, ...
                                              'cons.2.error.type', VALID_ERR_TYPES);
    end

    % --- Tiers 6-8: free text strings ------------------------------------
    trials.transcript(iTrial) = string(getSingleLabel(tg, 6, iTrial, baseName));
    trials.disfluency(iTrial) = string(getSingleLabel(tg, 7, iTrial, baseName));
    trials.comments(iTrial)   = string(getSingleLabel(tg, 8, iTrial, baseName));

    % --- Tier 9: unusable.trial  (double: NaN=unscored, 0=usable, 1=unusable)
    lbl = getSingleLabel(tg, 9, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), unusable.trial: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        trials.unusable_trial(iTrial) = val;   % stored as double, not logical
    end

    % --- Tier 10: difficult.to.score  (logical) --------------------------
    lbl = getSingleLabel(tg, 10, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), difficult.to.score: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        trials.difficult_to_score(iTrial) = logical(val);
    end

    % --- Tier 11: speech.epoch -------------------------------------------
    [on, off] = getEpochBounds(tg, 11, 'speech.epoch', iTrial, baseName);
    trials.speech_onset(iTrial)  = on;
    trials.speech_offset(iTrial) = off;

    % --- Tier 12: vowel.epoch --------------------------------------------
    [on, off] = getEpochBounds(tg, 12, 'vowel.epoch', iTrial, baseName);
    trials.vowel_onset(iTrial)  = on;
    trials.vowel_offset(iTrial) = off;

    % --- Tier 13: nontarget.sounds.epoch ---------------------------------
    trials.nontarget_sounds_epoch{iTrial} = getNtEpochs(tg, 13, ...
        'nontarget.sounds.epoch', iTrial, baseName);

    % --- Consistency checks: error (not warning) on any issue -----------
    issues = collectIssues(trials, iTrial);
    if ~isempty(issues)
        fprintf('\n=== Consistency issues for trial %d (%s) ===\n', iTrial, baseName);
        for ii = 1:numel(issues)
            fprintf('  [!] %s\n', issues{ii});
        end
        fprintf('Opening trial in Praat for inspection...\n');
        openInPraat(wavPath, tgPath, praatExe);
        error('compile_praat_subsyl:consistencyFailed', ...
              'Consistency check failed for trial %d (%s). Fix the TextGrid and re-run.', ...
              iTrial, baseName);
    end

end  % for iTrial

%% Summary -----------------------------------------------------------------
fprintf('\n=== Processing complete ===\n');
fprintf('Trials processed: %d  |  Skipped (no TextGrid): %d\n', ...
        nTrials - nMissing, nMissing);

%% Save --------------------------------------------------------------------
if ~isempty(paths.beh_annot_table)
    Texport  = trials;
    errFmt   = @(v) strjoin(arrayfun(@num2str, v, 'UniformOutput', false), ' ');
    Texport.cons1_error_type = cellfun(errFmt, trials.cons1_error_type, ...
                                       'UniformOutput', false);
    Texport.cons2_error_type = cellfun(errFmt, trials.cons2_error_type, ...
                                       'UniformOutput', false);
    Texport.nontarget_sounds_epoch = cellfun(@serializeNtEpochs, ...
        trials.nontarget_sounds_epoch, 'UniformOutput', false);
    writetable(Texport, paths.beh_annot_table, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('Textgrid scoring compiled into table file: %s\n', paths.beh_annot_table);
end



end


%% =========================================================================
%%  LOCAL FUNCTIONS
%% =========================================================================

function lbl = getSingleLabel(tg, tierInd, iTrial, baseName)
% Return trimmed label from a tier expected to have a single label.
% If multiple non-empty labels exist, warn and concatenate them.
    nInt = tgGetNumberOfIntervals(tg, tierInd);
    if nInt == 1
        lbl = strtrim(char(tgGetLabel(tg, tierInd, 1)));
        return;
    end
    parts = {};
    for j = 1:nInt
        l = strtrim(char(tgGetLabel(tg, tierInd, j)));
        if ~isempty(l)
            parts{end+1} = l; %#ok<AGROW>
        end
    end
    if     isempty(parts),  lbl = '';
    elseif numel(parts)==1, lbl = parts{1};
    else
        warning('Trial %d (%s), tier %d: multiple non-empty labels; concatenating.', ...
                iTrial, baseName, tierInd);
        lbl = strjoin(parts, ' ');
    end
end


function errVec = parseErrTypes(lbl, iTrial, baseName, tierName, valid)
% Parse space/comma/semicolon-delimited error-type string -> validated row vector.
    parts  = strsplit(strtrim(lbl), {' ',',',';',char(9)}, ...
                      'CollapseDelimiters', true);
    errVec = zeros(1, numel(parts));
    for i = 1:numel(parts)
        v = str2double(parts{i});
        if isnan(v)
            error('Trial %d (%s), %s: cannot parse "%s" as a number.', ...
                  iTrial, baseName, tierName, parts{i});
        end
        if ~ismember(v, valid)
            error('Trial %d (%s), %s: error type %g not in valid set %s.', ...
                  iTrial, baseName, tierName, v, mat2str(valid));
        end
        errVec(i) = v;
    end
end


function [onset, offset] = getEpochBounds(tg, tierInd, tierName, iTrial, baseName)
% 1 interval  -> not marked -> [NaN, NaN]
% 3 intervals -> onset = T1(2), offset = T2(2)
% other count -> error
    nInt = tgGetNumberOfIntervals(tg, tierInd);
    switch nInt
        case 1
            onset  = NaN;
            offset = NaN;
        case 3
            onset  = tg.tier{tierInd}.T1(2);
            offset = tg.tier{tierInd}.T2(2);
        otherwise
            error('Trial %d (%s), %s: expected 1 or 3 intervals, found %d.', ...
                  iTrial, baseName, tierName, nInt);
    end
end


function epochs = getNtEpochs(tg, tierInd, tierName, iTrial, baseName)
% Returns 2 x nEpochs (row1=onsets, row2=offsets), or zeros(2,0) if empty.
% Internal boundaries are treated as sequential onset/offset pairs.
    nInt = tgGetNumberOfIntervals(tg, tierInd);
    if nInt == 1, epochs = zeros(2,0); return; end
    bounds = tg.tier{tierInd}.T2(1:nInt-1);
    bounds = bounds(:)';
    if mod(numel(bounds), 2) ~= 0
        error('Trial %d (%s), %s: odd number of internal boundaries (%d); must be even.', ...
              iTrial, baseName, tierName, numel(bounds));
    end
    nEp    = numel(bounds) / 2;
    epochs = zeros(2, nEp);
    for k = 1:nEp
        epochs(1,k) = bounds(2*k-1);
        epochs(2,k) = bounds(2*k);
    end
end


function issues = collectIssues(T, iTrial)
% Collect all logical consistency issues for one trial.
% Returns a cell array of description strings; empty = no issues.
    issues = {};

    unusable = ~isnan(T.unusable_trial(iTrial)) && T.unusable_trial(iTrial) == 1;
    c1acc    = T.cons1_accuracy(iTrial);
    c2acc    = T.cons2_accuracy(iTrial);
    vacc     = T.vow_accuracy(iTrial);
    c1err    = T.cons1_error_type{iTrial};
    c2err    = T.cons2_error_type{iTrial};
    trans    = strtrim(char(T.transcript(iTrial)));
    comments = strtrim(char(T.comments(iTrial)));
    sp_on    = T.speech_onset(iTrial);
    sp_off   = T.speech_offset(iTrial);
    vw_on    = T.vowel_onset(iTrial);
    vw_off   = T.vowel_offset(iTrial);
    nt       = T.nontarget_sounds_epoch{iTrial};

    % Unusable trial: only check that comments are filled in, then stop
    if unusable
        if isempty(comments)
            issues{end+1} = 'unusable_trial=1 but comments is empty';
        end
        return;
    end

    % Required tiers for usable trials
    if isnan(c1acc), issues{end+1} = 'cons.1.accuracy is empty (required for usable trials)'; end
    if isnan(vacc),  issues{end+1} = 'vow.accuracy is empty (required for usable trials)';    end
    if isnan(c2acc), issues{end+1} = 'cons.2.accuracy is empty (required for usable trials)'; end
    if isnan(sp_on), issues{end+1} = 'speech.epoch is empty (required for usable trials)';    end
    if isnan(vw_on), issues{end+1} = 'vowel.epoch is empty (required for usable trials)';     end

    % Cluster accuracy / error-type consistency
    issues = [issues, checkCluster(c1acc, c1err, 'cons.1')];
    issues = [issues, checkCluster(c2acc, c2err, 'cons.2')];

    % Transcript required when any accuracy is below maximum
    anyError = (~isnan(c1acc) && c1acc < 2) || ...
               (~isnan(vacc)  && vacc  < 1) || ...
               (~isnan(c2acc) && c2acc < 2);
    if anyError && isempty(trans)
        issues{end+1} = 'accuracy below maximum but transcript is empty';
    end

    % Vowel epoch must lie within speech epoch
    if ~isnan(sp_on) && ~isnan(vw_on)
        if vw_on < sp_on
            issues{end+1} = sprintf('vowel onset (%.4f s) precedes speech onset (%.4f s)', ...
                                    vw_on, sp_on);
        end
        if vw_off > sp_off
            issues{end+1} = sprintf('vowel offset (%.4f s) exceeds speech offset (%.4f s)', ...
                                    vw_off, sp_off);
        end
    end

    % Non-target epochs must not overlap speech epoch
    if ~isempty(nt) && ~isnan(sp_on)
        for k = 1:size(nt,2)
            if nt(1,k) < sp_off && nt(2,k) > sp_on
                issues{end+1} = sprintf( ...
                    'non-target epoch %d [%.4f, %.4f] overlaps speech epoch [%.4f, %.4f]', ...
                    k, nt(1,k), nt(2,k), sp_on, sp_off);
            end
        end
    end
end


function clusterIssues = checkCluster(acc, errVec, tag)
% Consistency checks for one consonant cluster. Returns cell of issue strings.
    clusterIssues = {};
    if isnan(acc) || acc == 2, return; end   % full score or unscored: nothing to check

    non78 = errVec(~ismember(errVec, [7 8]));

    if isempty(errVec)
        clusterIssues{end+1} = sprintf( ...
            '%s.accuracy=%g but %s.error.type is empty', tag, acc, tag);
        return;   % remaining checks need errVec to be non-empty
    end

    if acc == 1 && isempty(non78)
        clusterIssues{end+1} = sprintf( ...
            '%s.accuracy=1 but only epenthesis/prothesis errors (types 7/8) listed; at least one other error type required', ...
            tag);
    end

    if acc == 0 && ~ismember(3,non78) && ~ismember(5,non78) && numel(non78) <= 1
        clusterIssues{end+1} = sprintf( ...
            '%s.accuracy=0 but errors [%s] are insufficient (need type 3, type 5, or >1 non-{7,8} errors)', ...
            tag, num2str(errVec));
    end
end


function openInPraat(wavPath, tgPath, praatExe)
% Load wav and TextGrid into Praat's object list for inspection.
% After Praat opens: select both objects in the list, then click View & Edit.
    if isempty(praatExe) || ~exist(praatExe, 'file')
        warning(['Praat executable not found at "%s".\n' ...
                 'Inspect manually:\n  WAV: %s\n  TG:  %s'], ...
                praatExe, wavPath, tgPath);
        return;
    end
    if ispc
        cmd = sprintf('start "" "%s" --open "%s" "%s"', praatExe, wavPath, tgPath);
    else   % macOS / Linux
        cmd = sprintf('"%s" --open "%s" "%s" &', praatExe, wavPath, tgPath);
    end
    system(cmd);
    pause(2);   % give Praat time to launch before MATLAB throws the error
end


function s = serializeNtEpochs(mat)
% Serialise 2×n epoch matrix to "on1,off1;on2,off2;..." for TSV export.
    if isempty(mat), s = ''; return; end
    pairs = arrayfun(@(k) sprintf('%.6f,%.6f', mat(1,k), mat(2,k)), ...
                     1:size(mat,2), 'UniformOutput', false);
    s = strjoin(pairs, ';');
end