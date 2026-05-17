
% =========================================================================
%  Load Praat TextGrid phonetic annotation data into a MATLAB trial table.
%
%  Dependencies:
%    - mpraat toolbox (https://github.com/bbTomas/mPraat)
%      Add mpraat directory to MATLAB path before running.
%    - MATLAB R2016b or later
%
%  TextGrid tier order (all IntervalTier, 13 total):
%    1  cons.1.accuracy        2  cons.1.error.type
%    3  vow.accuracy           4  cons.2.accuracy
%    5  cons.2.error.type      6  transcript
%    7  disfluency             8  comments
%    9  unusable.trial        10  difficult.to.score
%   11  speech.epoch          12  vowel.epoch
%   13  nontarget.sounds.epoch
% =========================================================================

vardefault('op',struct);

op.ses = 'subsyl'; 

field_default('op','sub','sml003');
field_default('op','task','test1');
field_default('op','chan','headphone');

paths = setpaths_dbs_learn(op);
lndmks = readtable(paths.landmarks_file,'FileType','text','Delimiter','tab');
op.run = lndmks.run(find(contains(lndmks.filename,op.task),1)); 
paths = setpaths_dbs_learn(op); % run-specific paths
paths.audiofiles_table_filename = [paths.trial_audio,filesep, paths.filestr, 'audiofiles-',op.chan, '.tsv']; 

% folder containing .TextGrid files and trial audio clips for this run and this recording chan
paths.trial_audio_task_chan = [paths.trial_audio, filesep, 'ses-',op.ses, '_task-',op.task,'_run-',num2str(op.run),'_',op.chan]; 

% --------------------------------------------------------------------------

%% Constants
EXPECTED_TIERS = { ...
    'cons.1.accuracy',    'cons.1.error.type', 'vow.accuracy',      ...
    'cons.2.accuracy',    'cons.2.error.type', 'transcript',        ...
    'disfluency',         'comments',          'unusable.trial',    ...
    'difficult.to.score', 'speech.epoch',      'vowel.epoch',       ...
    'nontarget.sounds.epoch'};
N_TIERS         = 13;
VALID_ERR_TYPES = [1 3 4 5 6 7 8];

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

%% Initialise output columns -----------------------------------------------
audiofiles.cons1_accuracy         = NaN(nTrials, 1);
audiofiles.cons1_error_type       = cell(nTrials, 1);
audiofiles.vow_accuracy           = NaN(nTrials, 1);
audiofiles.cons2_accuracy         = NaN(nTrials, 1);
audiofiles.cons2_error_type       = cell(nTrials, 1);
audiofiles.transcript             = repmat("", nTrials, 1);
audiofiles.disfluency             = repmat("", nTrials, 1);
audiofiles.comments               = repmat("", nTrials, 1);
audiofiles.unusable_trial         = false(nTrials, 1);
audiofiles.difficult_to_score     = false(nTrials, 1);
audiofiles.speech_onset           = NaN(nTrials, 1);
audiofiles.speech_offset          = NaN(nTrials, 1);
audiofiles.vowel_onset            = NaN(nTrials, 1);
audiofiles.vowel_offset           = NaN(nTrials, 1);
audiofiles.nontarget_sounds_epoch = cell(nTrials, 1);

for i = 1:nTrials
    audiofiles.cons1_error_type{i}       = [];
    audiofiles.cons2_error_type{i}       = [];
    audiofiles.nontarget_sounds_epoch{i} = zeros(2, 0);
end

%% Main loop ---------------------------------------------------------------
nMissing   = 0;
nWarnTotal = 0;

for iTrial = 1:nTrials

    wavName  = char(audiofiles.filename(iTrial));
    stimName = char(audiofiles.name(iTrial));
    [~, baseName, ~] = fileparts(wavName);
    tgPath   = [paths.trial_audio_task_chan, baseName, '.TextGrid'];

    fprintf('Trial %d/%d: %s\n', iTrial, nTrials, baseName);

    % File existence
    if ~exist(tgPath, 'file')
        warning('Trial %d: TextGrid not found -- %s. Skipping.', ...
                iTrial, tgPath);
        nMissing = nMissing + 1;
        continue;
    end

    % Stim-name check: filename must contain expected stim label
    if ~contains(baseName, stimName)
        error('Trial %d: stim "%s" expected in filename "%s" but not found.', ...
              iTrial, stimName, baseName);
    end

    % Load TextGrid
    tg     = tgRead(tgPath,'auto');
    nFound = tgGetNumberOfTiers(tg);

    % Validate tier count
    if nFound ~= N_TIERS
        error('Trial %d (%s): expected %d tiers, found %d.', ...
              iTrial, baseName, N_TIERS, nFound);
    end

    % Validate tier names
    for k = 1:N_TIERS
        gotName = tgGetTierName(tg, k);
        if ~strcmp(gotName, EXPECTED_TIERS{k})
            error('Trial %d (%s): tier %d expected "%s", got "%s".', ...
                  iTrial, baseName, k, EXPECTED_TIERS{k}, gotName);
        end
    end

    % ------------------------------------------------------------------
    % Tier 1: cons.1.accuracy  (0, 1, or 2)
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 1, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1 2])
            error('Trial %d (%s), cons.1.accuracy: expected 0/1/2, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        audiofiles.cons1_accuracy(iTrial) = val;
    end

    % ------------------------------------------------------------------
    % Tier 2: cons.1.error.type  (subset of {1,3,4,5,6,7,8})
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 2, iTrial, baseName);
    if ~isempty(lbl)
        audiofiles.cons1_error_type{iTrial} = parseErrTypes(lbl, iTrial, baseName, ...
                                         'cons.1.error.type', VALID_ERR_TYPES);
    end

    % ------------------------------------------------------------------
    % Tier 3: vow.accuracy  (0 or 1)
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 3, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), vow.accuracy: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        audiofiles.vow_accuracy(iTrial) = val;
    end

    % ------------------------------------------------------------------
    % Tier 4: cons.2.accuracy  (0, 1, or 2)
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 4, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1 2])
            error('Trial %d (%s), cons.2.accuracy: expected 0/1/2, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        audiofiles.cons2_accuracy(iTrial) = val;
    end

    % ------------------------------------------------------------------
    % Tier 5: cons.2.error.type  (subset of {1,3,4,5,6,7,8})
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 5, iTrial, baseName);
    if ~isempty(lbl)
        audiofiles.cons2_error_type{iTrial} = parseErrTypes(lbl, iTrial, baseName, ...
                                         'cons.2.error.type', VALID_ERR_TYPES);
    end

    % ------------------------------------------------------------------
    % Tiers 6-8: free text strings
    % ------------------------------------------------------------------
    audiofiles.transcript(iTrial) = string(getSingleLabel(tg, 6, iTrial, baseName));
    audiofiles.disfluency(iTrial) = string(getSingleLabel(tg, 7, iTrial, baseName));
    audiofiles.comments(iTrial)   = string(getSingleLabel(tg, 8, iTrial, baseName));

    % ------------------------------------------------------------------
    % Tier 9: unusable.trial  (0 or 1, stored as logical)
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 9, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), unusable.trial: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        audiofiles.unusable_trial(iTrial) = logical(val);
    end

    % ------------------------------------------------------------------
    % Tier 10: difficult.to.score  (0 or 1, stored as logical)
    % ------------------------------------------------------------------
    lbl = getSingleLabel(tg, 10, iTrial, baseName);
    if ~isempty(lbl)
        val = str2double(lbl);
        if isnan(val) || ~ismember(val, [0 1])
            error('Trial %d (%s), difficult.to.score: expected 0/1, got "%s".', ...
                  iTrial, baseName, lbl);
        end
        audiofiles.difficult_to_score(iTrial) = logical(val);
    end

    % ------------------------------------------------------------------
    % Tier 11: speech.epoch
    %   Annotator places 2 boundaries => 3 intervals.
    %   onset  = start of interval 2
    %   offset = end   of interval 2
    %   1 interval (no boundaries placed) => not scored => NaN
    % ------------------------------------------------------------------
    [on, off] = getEpochBounds(tg, 11, 'speech.epoch', iTrial, baseName);
    audiofiles.speech_onset(iTrial)  = on;
    audiofiles.speech_offset(iTrial) = off;

    % ------------------------------------------------------------------
    % Tier 12: vowel.epoch  (same convention as speech.epoch)
    % ------------------------------------------------------------------
    [on, off] = getEpochBounds(tg, 12, 'vowel.epoch', iTrial, baseName);
    audiofiles.vowel_onset(iTrial)  = on;
    audiofiles.vowel_offset(iTrial) = off;

    % ------------------------------------------------------------------
    % Tier 13: nontarget.sounds.epoch
    %   Internal boundaries are sequential onset/offset pairs (even count).
    %   Stored as 2 x nEpochs: row 1 = onsets, row 2 = offsets.
    % ------------------------------------------------------------------
    audiofiles.nontarget_sounds_epoch{iTrial} = getNtEpochs(tg, 13, ...
        'nontarget.sounds.epoch', iTrial, baseName);

    % Consistency checks for this trial
    nWarnTotal = nWarnTotal + runChecks(audiofiles, iTrial, baseName);

end  % for iTrial

%% Summary -----------------------------------------------------------------
fprintf('\n=== Processing complete ===\n');
fprintf('Trials processed:     %d\n', nTrials - nMissing);
fprintf('Missing TextGrids:    %d\n', nMissing);
fprintf('Consistency warnings: %d\n', nWarnTotal);

%% Optional TSV save -------------------------------------------------------
if ~isempty(paths.beh_annot_table)
    % Cell columns must be serialised to strings for TSV output.
    % Workspace variable T retains the native cell/array types.
    Texport  = audiofiles;
    errFmt   = @(v) strjoin(arrayfun(@num2str, v, 'UniformOutput', false), ' ');
    Texport.cons1_error_type = cellfun(errFmt, audiofiles.cons1_error_type, ...
                                       'UniformOutput', false);
    Texport.cons2_error_type = cellfun(errFmt, audiofiles.cons2_error_type, ...
                                       'UniformOutput', false);
    Texport.nontarget_sounds_epoch = cellfun(@serializeNtEpochs, ...
        audiofiles.nontarget_sounds_epoch, 'UniformOutput', false);
    writetable(Texport, paths.beh_annot_table, 'Delimiter', '\t', 'FileType', 'text');
    fprintf('\nAugmented table saved to: %s\n', paths.beh_annot_table);
end


%% =========================================================================
%%  LOCAL FUNCTIONS
%% =========================================================================

function lbl = getSingleLabel(tg, tierInd, iTrial, baseName)
% Return trimmed label from a tier expected to have one interval.
% If the tier has multiple non-empty labels, warn and concatenate them.
    nInt = tgGetNumberOfIntervals(tg, tierInd);
    if nInt == 1
        lbl = strtrim(char(tgGetLabel(tg, tierInd, 1)));
        return;
    end
    parts = {};
    for j = 1:nInt
        l = strtrim(char(tgGetIntervalLabel(tg, tierInd, j)));
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
% Parse space/comma/semicolon-delimited error-type string into a
% validated numeric row vector.
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
% Extract epoch onset/offset from a Praat epoch tier.
%   1 interval  => not marked   => [NaN, NaN]
%   3 intervals => onset  = tg.tier{tierInd}.T1(2)
%                  offset = tg.tier{tierInd}.T2(2)
%   other count => error
% Direct struct access is used because mpraat stores boundary times
% as column vectors in tg.tier{i}.T1 and tg.tier{i}.T2.
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
% Extract non-target sound onset/offset pairs from tier 13.
% Returns 2 x nEpochs (row1 = onsets, row2 = offsets).
% Returns zeros(2,0) when no epochs were marked (1 interval only).
%
% Internal boundaries = T2 values of intervals 1 through nInt-1,
% treated as sequential onset/offset pairs (count must be even).
    nInt = tgGetNumberOfIntervals(tg, tierInd);
    if nInt == 1
        epochs = zeros(2, 0);
        return;
    end
    bounds = tg.tier{tierInd}.T2(1:nInt-1);
    bounds = bounds(:)';   % row vector
    nBound = numel(bounds);
    if mod(nBound, 2) ~= 0
        error('Trial %d (%s), %s: odd number of internal boundaries (%d); must be even (onset/offset pairs).', ...
              iTrial, baseName, tierName, nBound);
    end
    nEp    = nBound / 2;
    epochs = zeros(2, nEp);
    for k = 1:nEp
        epochs(1, k) = bounds(2*k - 1);   % onset
        epochs(2, k) = bounds(2*k);       % offset
    end
end


function nW = runChecks(T, iTrial, id)
% Run logical consistency checks for one trial.
% Issues MATLAB warnings (does not halt). Returns warning count.
    nW = 0;

    unusable = T.unusable_trial(iTrial);
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

    % -- UNUSABLE trial: only check that comments are filled in -----------
    if unusable
        if isempty(comments)
            warning('Trial %d (%s): unusable_trial=true but comments is empty.', ...
                    iTrial, id);
            nW = nW + 1;
        end
        return;  % skip all other checks
    end

    % -- Required tiers for usable trials ---------------------------------
    reqNames = {'cons.1.accuracy', 'vow.accuracy', 'cons.2.accuracy'};
    reqVals  = {c1acc,              vacc,           c2acc};
    for r = 1:3
        if isnan(reqVals{r})
            warning('Trial %d (%s): %s is empty (required for usable trials).', ...
                    iTrial, id, reqNames{r});
            nW = nW + 1;
        end
    end
    if isnan(sp_on)
        warning('Trial %d (%s): speech.epoch is empty (required for usable trials).', ...
                iTrial, id);
        nW = nW + 1;
    end
    if isnan(vw_on)
        warning('Trial %d (%s): vowel.epoch is empty (required for usable trials).', ...
                iTrial, id);
        nW = nW + 1;
    end

    % -- Cluster accuracy / error-type consistency ------------------------
    if ~isnan(c1acc) && c1acc < 2
        nW = nW + checkCluster(c1acc, c1err, 'cons.1', iTrial, id);
    end
    if ~isnan(c2acc) && c2acc < 2
        nW = nW + checkCluster(c2acc, c2err, 'cons.2', iTrial, id);
    end

    % -- Transcript required when any accuracy below maximum --------------
    anyError = (~isnan(c1acc) && c1acc < 2) || ...
               (~isnan(vacc)  && vacc  < 1) || ...
               (~isnan(c2acc) && c2acc < 2);
    if anyError && isempty(trans)
        warning('Trial %d (%s): errors scored but transcript is empty.', ...
                iTrial, id);
        nW = nW + 1;
    end

    % -- Vowel epoch must lie within speech epoch -------------------------
    if ~isnan(sp_on) && ~isnan(vw_on)
        if vw_on < sp_on
            warning('Trial %d (%s): vowel onset (%.4f s) precedes speech onset (%.4f s).', ...
                    iTrial, id, vw_on, sp_on);
            nW = nW + 1;
        end
        if vw_off > sp_off
            warning('Trial %d (%s): vowel offset (%.4f s) exceeds speech offset (%.4f s).', ...
                    iTrial, id, vw_off, sp_off);
            nW = nW + 1;
        end
    end

    % -- Non-target epochs must not overlap speech epoch ------------------
    if ~isempty(nt) && ~isnan(sp_on)
        for k = 1:size(nt, 2)
            nt_on  = nt(1, k);
            nt_off = nt(2, k);
            if nt_on < sp_off && nt_off > sp_on   % standard interval-overlap test
                warning(['Trial %d (%s): non-target epoch %d [%.4f, %.4f] ' ...
                         'overlaps speech epoch [%.4f, %.4f].'], ...
                        iTrial, id, k, nt_on, nt_off, sp_on, sp_off);
                nW = nW + 1;
            end
        end
    end
end


function nW = checkCluster(acc, errVec, tag, iTrial, id)
% Check accuracy/error-type logical consistency for one CC cluster.
% Returns number of warnings issued.
    nW    = 0;
    non78 = errVec(~ismember(errVec, [7 8]));

    % Error types must not be empty when accuracy < 2
    if isempty(errVec)
        warning('Trial %d (%s): %s.accuracy=%g but %s.error.type is empty.', ...
                iTrial, id, tag, acc, tag);
        nW = nW + 1;
        return;
    end

    if acc == 1
        % 1-point score requires at least one non-epenthesis/prothesis error
        if isempty(non78)
            warning(['Trial %d (%s): %s.accuracy=1 but only epenthesis/' ...
                     'prothesis errors (types 7/8) listed; need at least ' ...
                     'one other error type.'], iTrial, id, tag);
            nW = nW + 1;
        end

    elseif acc == 0
        % 0-point score requires deletion (3), substitution (5),
        % or more than one non-{7,8} error
        if ~ismember(3, non78) && ~ismember(5, non78) && numel(non78) <= 1
            warning(['Trial %d (%s): %s.accuracy=0 but error types [%s] ' ...
                     'do not satisfy zero-score criteria ' ...
                     '(need error 3, error 5, or >1 non-{7/8} errors).'], ...
                    iTrial, id, tag, num2str(errVec));
            nW = nW + 1;
        end
    end
end


function s = serializeNtEpochs(mat)
% Convert 2×n epoch matrix to "on1,off1;on2,off2;..." for TSV export.
    if isempty(mat)
        s = '';
        return;
    end
    pairs = arrayfun(@(k) sprintf('%.6f,%.6f', mat(1,k), mat(2,k)), ...
                     1:size(mat,2), 'UniformOutput', false);
    s = strjoin(pairs, ';');
end