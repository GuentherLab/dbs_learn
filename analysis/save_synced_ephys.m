 %% Top-level script for saving master-aligned (CED time) FieldTrip ephys files
clear;
clc;

op.sub = 'sap022';
op.ses = 'stut';
paths = setpaths_dbs_learn(op);

op.ephys_to_align_dir     = paths.src_dbs;
op.ephys_aligned_dir      = paths.fieldtrip;
op.ephys_to_align_suffix  = 'preproc';
op.time_block_labels      = 'step';
op.ephys_to_align_varname = 'PCP';
op.paths                  = paths; % saved into the new FieldTrip file

% Input sync-landmarks table (.xlsx, .tsv, or MATLAB table variable)
sync_in = paths.landmarks_file;

[sync_summary, op_out] = save_synced_ephys_(sync_in, op);





function [summary_tbl, op] = save_synced_ephys_(sync_in, op)
% SAVE_SYNCED_EPHYS  Export FieldTrip ephys files with time aligned to master (CED) clock.
%
%   [summary_tbl, op] = save_synced_ephys(sync_in, op)
%
%   Inputs
%     sync_in  - Path to sync-landmarks table (.tsv or .xlsx) or a MATLAB table
%                containing rows with a 'master' column (master == 1 for CED rows).
%     op       - Options struct with fields:
%                  .sub                    - Subject ID (e.g., 'sap022') [required]
%                  .ses                    - Session ID (e.g., 'stut')   [required]
%                  .ephys_to_align_dir     - Dir containing unaligned ephys .mat files
%                                            (e.g., paths.src_dbs)
%                  .ephys_aligned_dir      - Output dir for aligned FieldTrip .mat files
%                                            (e.g., paths.fieldtrip)
%                  .time_block_labels      - 'run', 'step', or cell array {'run','step'}
%                                            (default: 'step'; throws error otherwise)
%                  .ephys_to_align_suffix  - Suffix before timestamp in source .mat files
%                                            (default: 'preproc')
%                  .ephys_to_align_varname - Struct variable name in source .mat holding
%                                            FieldTrip .data struct (default: 'PCP')
%                  .paths                  - Optional paths struct to store inside saved op
%
%   Outputs
%     summary_tbl - Table summarizing status, source file, and saved output per master row.
%     op          - Options struct with defaults populated (.task is removed after loop).

    %% 1. Parse and validate options
    if nargin < 2 || isempty(op)
        op = struct();
    end

    assert(isfield(op, 'sub') && ~isempty(op.sub), 'save_synced_ephys:MissingSub', ...
        'op.sub must be specified.');
    assert(isfield(op, 'ses') && ~isempty(op.ses), 'save_synced_ephys:MissingSes', ...
        'op.ses must be specified.');
    op.sub = char(op.sub);
    op.ses = char(op.ses);

    % Defaults
    if ~isfield(op, 'ephys_to_align_suffix') || isempty(op.ephys_to_align_suffix)
        op.ephys_to_align_suffix = 'preproc';
    end
    op.ephys_to_align_suffix = regexprep(char(op.ephys_to_align_suffix), '^[_-]+', '');

    if ~isfield(op, 'ephys_to_align_varname') || isempty(op.ephys_to_align_varname)
        op.ephys_to_align_varname = 'PCP';
    end
    op.ephys_to_align_varname = char(op.ephys_to_align_varname);

    if ~isfield(op, 'time_block_labels') || isempty(op.time_block_labels)
        op.time_block_labels = 'step';
    end

    % Validate time_block_labels ('run' and/or 'step' only)
    if ischar(op.time_block_labels) || isstring(op.time_block_labels)
        block_labels = cellstr(op.time_block_labels);
    elseif iscell(op.time_block_labels)
        block_labels = cellfun(@char, op.time_block_labels, 'UniformOutput', false);
    else
        error('save_synced_ephys:InvalidTimeBlockLabels', ...
            'op.time_block_labels must be ''run'', ''step'', or a cell array of those strings.');
    end
    block_labels = block_labels(:)';
    if isempty(block_labels) || ~all(ismember(block_labels, {'run', 'step'}))
        error('save_synced_ephys:InvalidTimeBlockLabels', ...
            'op.time_block_labels only accepts ''run'' and/or ''step''.');
    end

    % Resolve directories
    if (~isfield(op, 'ephys_to_align_dir') || isempty(op.ephys_to_align_dir)) && ...
            isfield(op, 'paths') && isfield(op.paths, 'src_dbs')
        op.ephys_to_align_dir = op.paths.src_dbs;
    end
    assert(isfield(op, 'ephys_to_align_dir') && exist(op.ephys_to_align_dir, 'dir') == 7, ...
        'save_synced_ephys:MissingInputDir', ...
        'op.ephys_to_align_dir does not exist or was not provided.');
    op.ephys_to_align_dir = char(op.ephys_to_align_dir);

    if (~isfield(op, 'ephys_aligned_dir') || isempty(op.ephys_aligned_dir)) && ...
            isfield(op, 'paths') && isfield(op.paths, 'fieldtrip')
        op.ephys_aligned_dir = op.paths.fieldtrip;
    end
    assert(isfield(op, 'ephys_aligned_dir') && ~isempty(op.ephys_aligned_dir), ...
        'save_synced_ephys:MissingOutputDir', ...
        'op.ephys_aligned_dir must be specified.');
    op.ephys_aligned_dir = char(op.ephys_aligned_dir);

    if exist(op.ephys_aligned_dir, 'dir') ~= 7
        mkdir(op.ephys_aligned_dir);
    end

    %% 2. Load sync-landmarks table
    tbl = local_load_sync_table(sync_in);

    required_cols = [{'master', 'filename', 'task'}, block_labels];
    missing_cols = setdiff(required_cols, tbl.Properties.VariableNames);
    assert(isempty(missing_cols), 'save_synced_ephys:MissingColumns', ...
        'Sync table is missing required column(s): %s', strjoin(missing_cols, ', '));

    has_suffix_col = ismember('suffix', tbl.Properties.VariableNames);

    % Filter to master rows (master == 1)
    is_master = local_to_logical_master(tbl.master);
    tbl_master = tbl(is_master, :);

    fprintf('save_synced_ephys: found %d master row(s) in sync table.\n', height(tbl_master));

    % Preallocate summary tracking
    nRows = height(tbl_master);
    sum_task       = cell(nRows, 1);
    sum_block      = cell(nRows, 1);
    sum_status     = cell(nRows, 1);
    sum_source_mat = cell(nRows, 1);
    sum_output_mat = cell(nRows, 1);

    %% 3. Process each master row
    for iRow = 1:nRows
        % Note task for current iteration (removed via rmfield at end of iteration)
        task_val = strtrim(char(string(tbl_master.task(iRow))));
        op.task = task_val;
        sum_task{iRow} = task_val;

        % Build block label string (e.g., '_step-C001' or '_run-1_step-C001')
        block_parts = cell(1, numel(block_labels));
        for iB = 1:numel(block_labels)
            lbl = block_labels{iB};
            raw_val = tbl_master.(lbl)(iRow);
            if iscell(raw_val)
                raw_val = raw_val{1};
            end
            val_str = strtrim(char(string(raw_val)));
            if isempty(val_str) || strcmpi(val_str, 'nan') || strcmpi(val_str, '<missing>')
                error('save_synced_ephys:EmptyBlockValue', ...
                    'Master row %d has an empty/NaN value in column ''%s''.', iRow, lbl);
            end
            block_parts{iB} = sprintf('_%s-%s', lbl, val_str);
        end
        block_str = strjoin(block_parts, '');
        sum_block{iRow} = regexprep(block_str, '^_', '');

        % Determine master-audio basename (filename minus extension, suffix, and leading _/-)
        fname_raw = strtrim(char(string(tbl_master.filename(iRow))));
        if has_suffix_col
            sfx_raw = strtrim(char(string(tbl_master.suffix(iRow))));
        else
            sfx_raw = '';
        end
        master_basename = local_get_master_basename(fname_raw, sfx_raw);

        % Construct ephys-to-align basename
        ephys_basename = [master_basename, '_', op.ephys_to_align_suffix];

        % Search ephys_to_align_dir and immediate subdirectories (1 level below)
        matches = local_find_matching_mats(op.ephys_to_align_dir, ephys_basename);

        if numel(matches) > 1
            % Clean up op.task before throwing error
            if isfield(op, 'task')
                op = rmfield(op, 'task');
            end
            error('save_synced_ephys:MultipleMatches', ...
                'Multiple matching .mat files found for basename ''%s'':\n  %s', ...
                ephys_basename, strjoin(matches, '\n  '));
        elseif isempty(matches)
            fprintf('[%d/%d] SKIPPED (no matching %s*.mat found): %s\n', ...
                iRow, nRows, ephys_basename, fname_raw);
            sum_status{iRow}     = 'skipped_no_match';
            sum_source_mat{iRow} = '';
            sum_output_mat{iRow} = '';
        else
            src_fpath = matches{1};
            sum_source_mat{iRow} = src_fpath;

            % Check variables inside the .mat file before loading
            varname = op.ephys_to_align_varname;
            file_vars = {whos('-file', src_fpath).name};
            assert(ismember('m_Align', file_vars), 'save_synced_ephys:MissingMAlign', ...
                'Struct ''m_Align'' not found in %s', src_fpath);
            assert(ismember(varname, file_vars), 'save_synced_ephys:MissingEphysVar', ...
                'Variable ''%s'' not found in %s', varname, src_fpath);

            S = load(src_fpath, 'm_Align', varname);

            assert(isstruct(S.m_Align) && ...
                isfield(S.m_Align, 'fcn_translateTime') && ...
                isfield(S.m_Align, 'TimeTrans_PCP_to_CED_byEvent') && ...
                ~isempty(S.m_Align.TimeTrans_PCP_to_CED_byEvent), ...
                'save_synced_ephys:InvalidMAlign', ...
                'm_Align in %s is missing fcn_translateTime or TimeTrans_PCP_to_CED_byEvent.', ...
                src_fpath);

            ephys_struct = S.(varname);
            assert(isstruct(ephys_struct) && isfield(ephys_struct, 'data') && ...
                isstruct(ephys_struct.data), ...
                'save_synced_ephys:MissingDataField', ...
                'Variable ''%s'' in %s does not contain a .data struct.', varname, src_fpath);

            D = ephys_struct.data;
            ft_req = {'trial', 'time', 'fsample', 'label'};
            assert(all(isfield(D, ft_req)) && iscell(D.time) && ~isempty(D.time), ...
                'save_synced_ephys:InvalidFieldTripStruct', ...
                '%s.data in %s is missing required FieldTrip fields (trial, time, fsample, label).', ...
                varname, src_fpath);

            % Translate time vector(s) to master (CED) coordinates
            for iTrl = 1:numel(D.time)
                t_master_aligned = S.m_Align.fcn_translateTime( ...
                    D.time{iTrl}, S.m_Align.TimeTrans_PCP_to_CED_byEvent);
                D.time{iTrl} = t_master_aligned;
            end

            % Build BIDS output filename and save D and op
            out_fname = sprintf('sub-%s_ses-%s_task-%s%s_ft-ref.mat', ...
                op.sub, op.ses, op.task, block_str);
            out_fpath = fullfile(op.ephys_aligned_dir, out_fname);

            save(out_fpath, 'D', 'op', '-v7.3');

            fprintf('[%d/%d] SAVED: %s  (t = [%.3f, %.3f] s)\n', ...
                iRow, nRows, out_fname, D.time{1}(1), D.time{1}(end));
            sum_status{iRow}     = 'saved';
            sum_output_mat{iRow} = out_fpath;
        end

        % Remove iteration-specific op.task at the end of every iteration
        if isfield(op, 'task')
            op = rmfield(op, 'task');
        end
    end

    summary_tbl = table(sum_task, sum_block, sum_status, sum_source_mat, sum_output_mat, ...
        'VariableNames', {'task', 'time_block', 'status', 'source_file', 'saved_file'});
end


% =========================================================================
function tbl = local_load_sync_table(sync_in)
% Load .tsv, .xlsx, or accept a MATLAB table directly while preserving text IDs.
    if istable(sync_in)
        tbl = sync_in;
        return;
    end

    assert((ischar(sync_in) || isstring(sync_in)) && exist(sync_in, 'file') == 2, ...
        'save_synced_ephys:SyncFileNotFound', ...
        'Sync table file not found: %s', char(string(sync_in)));

    [~, ~, ext] = fileparts(char(sync_in));
    switch lower(ext)
        case {'.tsv', '.txt'}
            opts = detectImportOptions(sync_in, 'FileType', 'text', 'Delimiter', '\t');
        case {'.xlsx', '.xls'}
            opts = detectImportOptions(sync_in, 'FileType', 'spreadsheet');
        otherwise
            error('save_synced_ephys:UnsupportedTableFormat', ...
                'Unsupported sync table extension ''%s''. Use .tsv or .xlsx.', ext);
    end

    opts.VariableNamingRule = 'preserve';
    text_cols = intersect(opts.VariableNames, ...
        {'sub', 'ses', 'task', 'step', 'run', 'filename', 'suffix', 'filetype', 'filepath'});
    if ~isempty(text_cols)
        opts = setvartype(opts, text_cols, 'string');
    end
    tbl = readtable(sync_in, opts);
end


% =========================================================================
function is_m = local_to_logical_master(master_col)
% Robustly convert numeric, logical, or string/cell master column to logical.
    if islogical(master_col)
        is_m = master_col;
    elseif isnumeric(master_col)
        is_m = (master_col == 1);
    else
        s = strtrim(string(master_col));
        is_m = (s == "1") | strcmpi(s, "true");
    end
end


% =========================================================================
function base = local_get_master_basename(fname, sfx)
% Strip extension, suffix column string, and any trailing underscore or dash.
    [~, fname_stem, ~] = fileparts(fname);

    if ~isempty(sfx) && ~strcmpi(sfx, 'nan') && ~strcmpi(sfx, '<missing>')
        sfx_clean = regexprep(sfx, '^[_-]+', '');
        [~, sfx_stem, sfx_ext] = fileparts(sfx_clean);
        if isempty(sfx_stem)
            sfx_stem = sfx_clean;
        end

        % Try matching full suffix (if filename had double ext) or stem
        candidates = unique({sfx_clean, sfx_stem}, 'stable');
        for iC = 1:numel(candidates)
            cand = candidates{iC};
            if ~isempty(cand) && endsWith(fname_stem, cand)
                fname_stem = fname_stem(1:(end - length(cand)));
                break;
            end
        end
        % Also handle case where suffix had an extension like 'voice.wav'
        if ~isempty(sfx_ext) && endsWith(fname_stem, sfx_stem)
            fname_stem = fname_stem(1:(end - length(sfx_stem)));
        end
    end

    base = regexprep(fname_stem, '[_-]+$', '');
end


% =========================================================================
function matches = local_find_matching_mats(root_dir, basename_prefix)
% Search root_dir and immediate subdirectories (1 level deep) for <basename_prefix>*.mat.
    search_dirs = {root_dir};
    entries = dir(root_dir);
    for iE = 1:numel(entries)
        if entries(iE).isdir && ~ismember(entries(iE).name, {'.', '..'})
            search_dirs{end+1} = fullfile(root_dir, entries(iE).name); %#ok<AGROW>
        end
    end

    matches = {};
    for iD = 1:numel(search_dirs)
        d = dir(fullfile(search_dirs{iD}, [basename_prefix, '*.mat']));
        for iF = 1:numel(d)
            if ~d(iF).isdir && startsWith(d(iF).name, basename_prefix)
                matches{end+1} = fullfile(search_dirs{iD}, d(iF).name); %#ok<AGROW>
            end
        end
    end
    matches = unique(matches, 'stable');
end