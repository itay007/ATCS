function Data = affysnpemconvert(Data, EA_snpid, Mapping_snpid)
%AFFYSNPEMCONVERT A conversion function for the AFFYSNPCNVDEMO.
%
%   DATA = AFFYSNPEMCONVERT(DATA, EA_SNPID, MAPPING_SNPID) converts the SNP
%   probe structure DATA of Affymetrix Early Access array to commercial
%   version Mapping arrays. For example, Early Access 50 Xba array is the
%   same as Mapping 50K Xba array, only some of the probes on 50KXba array
%   are masked on the commercial array. DATA is a structure return from the
%   function CELINTENSITYREAD containing probe set IDs and probe indices
%   and PM intensities for allele A  and allele B.  EA_SNPID is a cell
%   array containing the probe set IDs for SNPs on an Early Access array,
%   and MAPPING_SNPID is a cell array containing probe IDs for SNPs on a
%   Mapping array. The return DATA contains data for probes that exist on
%   both arrays, and the probe set IDs in the ProbeSetIDs field are
%   converted to Mapping probe set IDs. 

%   Copyright 2008 The MathWorks, Inc.


%== Convert ProbesetIDs
% Usually the Early Access array ID read from XLS file are numerical.  
if isnumeric(EA_snpid)
    probesetids = str2double(Data.ProbeSetIDs);
elseif iscell(EA_snpid)
    probesetids = Data.ProbeSetIDs;
end

%== Find not rejected ID
[snp_tf, loc] = ismember(probesetids, EA_snpid);

%== Find the corresponding probes
probe_s = find(Data.ProbeIndices == 0); % The start of each probe set
probe_e = [probe_s(2:end)-1; numel(Data.ProbeIndices)];%The end of each probe set
probe_tf = false(numel(Data.ProbeIndices), 1);
for i = 1:length(snp_tf)
    probe_tf(probe_s(i):probe_e(i)) = snp_tf(i);
end

%== Update Data structure
Data.NumProbeSets = sum(snp_tf);
Data.NumProbes = sum(probe_tf);
Data.ProbeSetIDs = Mapping_snpid(loc(snp_tf));
Data.ProbeIndices = Data.ProbeIndices(probe_tf);
Data.GroupNumbers = Data.GroupNumbers(probe_tf);
Data.PMIntensities = Data.PMIntensities(probe_tf, :);

if isfield(Data, 'MMAIntensities')
    Data.MMAIntensities = Data.MMAIntensities(probe_tf, :);
    Data.MMBIntensities = Data.MMBIntensities(probe_tf, :);
end
