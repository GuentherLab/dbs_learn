%%% load CED data and percept data 
% percept data will be synced to CED time coordinates
% all data in fieldtrip format
%
% 'physdata' contains CED and non-synced percept data
% 'pcp' contains synced CED data
%
% input cfg contains fields 'sub','ses','task'
%
%%% syncing/aligning code from Zeyang Yu


function [physdat, pcp] = load_ced_and_percept_data(cfg)

setpaths_dbsmulti()

physiopath = [paths.data, filesep, ['sub-',cfg.sub], filesep, ['ses-',num2str(cfg.ses)], filesep, 'physio'];
physiodirlist = struct2table(dir(physiopath)); 
last_preproc_file_idx = find(contains(physiodirlist.name,[task2dir(cfg.task),'_preproc']),1,'last'); % most updated preproc file
preprocfile = physiodirlist.name{last_preproc_file_idx,1}; 

% load preproc containing CED and un-aligned Percept data
physdat = load([physiopath, filesep, preprocfile]);

%%%%%%%%%%% translate Percept times into CED times
physdat = ft_checkconfig(physdat,'renamed',{'preprocEMG','CED'});
physdat = ft_checkconfig(physdat,'renamed',{'preprocLFP','PCP'});
ppCED = physdat.CED;
ppBHV = physdat.bhv;
ppPCP = physdat.PCP;
m_Align    = physdat.m_Align;
pcp = ppPCP.data;
t_PCP = pcp.time{1};
t_CED = ppCED.stimArt.time{1};
t_BHV = [];
% aligned times
fcn_translateTime = m_Align.fcn_translateTime; % function to apply the translation
tTrans_PCP_to_CED = m_Align.TimeTrans_PCP_to_CED_byEvent; % translating PCP time into CED time
tTrans_BHV_to_CED = m_Align.TimeTrans_BHV_to_CED; % translating BHV time into CED time
%
if ~isempty(tTrans_PCP_to_CED)
    t_PCP_inCED = fcn_translateTime(t_PCP,tTrans_PCP_to_CED);
else
    t_PCP_inCED = [];
end
pcp.oldtime = pcp.time;
pcp.time{1} = t_PCP_inCED;