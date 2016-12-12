function [cmz po] = mspalign(p,varargin)
%MSPALIGN aligns mass spectrometry peaks from different scans.
%
%  [CMZ POUT] = MSPALIGN(P) aligns mass spectra from multiple peak lists
%  (centroided data), by first estimating CMZ, a vector of common
%  mass/charge (m/z) ratio values estimated by considering the peaks in all
%  spectra in P, a cell array of peak lists, where each element corresponds
%  to a spectrum or retention time. Each peak list in P is represented by a
%  matrix with two columns with peak information, the first column has the
%  mass-charge values and the second column has the ion intensities. It
%  then aligns the peaks in each spectrum to the values in CMZ, creating
%  POUT, a cell array of aligned peak lists with the same dimensions as P
%  but with the adjusted mass-charge values.
%
%  MSPALIGN works better with scans that are represented by all the
%  extracted peaks or the so called centroided peaks. 
%
%  MSPALIGN(...,'QUANTILEVALUE',Q) determines which peaks are selected by
%  the estimation method to create CMZ, the vector of common m/z values.
%  Choices are a scalar between 0 and 1. Default is 0.95.
%
%  MSPALIGN(...,'ESTIMATIONMETHOD',EM) sets the method used to find the
%  common mass-charge vector (CMZ). Options are:
%
%    'histogram'    - Peak locations are clustered using a kernel density
%     (default)       estimation approach. The peak ion intensity is used
%                     as a weighting factor. The center of all the found
%                     clusters conform the CMZ vector. 
%
%    'regression'   - Takes a sample of the distances between observed
%                     significant peaks and regresses the inter-peak
%                     distance to create a CMZ vector with a similar
%                     inter-element distances.
%
%  
%  MSPALIGN(...,'CORRECTIONMETHOD',CM) sets the method to align each scan
%  to the CMZ.
%
%    'nearest-neighbor' - For each common peak in the CMZ vector, its
%        (default)        counterpart in each peak list is the peak that is
%                         closest to the common peak's m/z value.
%
%    'shortest-path' - For each common peak in the CMZ vector, its
%                      counterpart in each peak list is selected using the
%                      shortest path algorithm. 
%
%  MSPALIGN(...,'SHOWESTIMATION',true) displays assessment plots respective
%  to the estimation method and the common mass-charge vector. Defaults is
%  FALSE unless MSPALIGN is called without output arguments.
% 
%   Example:
%
%      load lcmsdata
%
%      % Show the unaligned LCMS dataset:
%      [MZ,Y] = msppresample(ms_peaks,5000);
%      msheatmap(MZ,ret_time,log(Y))
%      msdotplot(ms_peaks,ret_time)
%
%      % Align the peaks using the default method:
%      [CMZ, aligned_peaks] = mspalign(ms_peaks);
%
%      % Show the aligned LCMS dataset:
%      [MZ2,Y2] = msppresample(aligned_peaks,5000);
%      msheatmap(MZ2,ret_time,log(Y2))
%      msdotplot(aligned_peaks,ret_time)
%            
%      % Link the axes and zoom in to observe the detail:
%      linkaxes(findobj(0,'Tag','MSHeatMap'))
%      axis([480 532 375 485])
%
%   See also DIFFPROTDEMO, LCMSDEMO, MSALIGN, MSDOTPLOT, MSHEATMAP,
%   MSPEAKS, MZXML2PEAKS.

%   Copyright 2006-2008 The MathWorks, Inc.


% References:
% [1] Jeffries, N. "Algorithms for alignment of mass spectrometry proteomic
%     data", Bioinformatics, 21(14):3066-3073, 2005.
%
% [2] Purvine, S. et.al. "Spectral Quality Assessment for High-Throughput
%     Tandem Mass Spectrometry Proteomics", OMICS A Journal of Integrative
%     Biology 8(3), 2004.

% check inputs
bioinfochecknargin(nargin,1,mfilename);

% Set defaults
quan_int = .95;
histResolution = .1;
correctionMethod = 'nearest-neighbor';
estimationMethod = 'histogram';
if nargout==0
    showEstimation = true;
else
    showEstimation = false;
end

if  nargin > 1
    if rem(nargin,2) ~= 1
        error(message('bioinfo:mspalign:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'quantilevalue','estimationmethod','correctionmethod','showestimation'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,length(pname)));
        if isempty(k)
            error(message('bioinfo:mspalign:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:mspalign:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % quantile
                    if ~isnumeric(pval) || ~isscalar(pval) || pval<0 || pval>1
                        error(message('bioinfo:mspalign:badLevels'))
                    end
                    quan_int = pval;
                case 2  %'estimationmethod'
                    estimationMethods = {'histogram','regression'};
                    estimationMethod = strmatch(lower(pval),estimationMethods); 
                    if isempty(estimationMethod) 
                        error(message('bioinfo:mspalign:NotValidEstimationMethod'))
                    end
                    estimationMethod = estimationMethods{estimationMethod};
                case 3  %'correctionmethod'
                    correctionMethods = {'nearest-neighbor','shortest-path'};
                    correctionMethod = strmatch(lower(pval),correctionMethods); 
                    if isempty(correctionMethod) 
                        error(message('bioinfo:mspalign:NotValidCorrectionMethod'))
                    end
                    correctionMethod = correctionMethods{correctionMethod};    
                case 4 %'showestimation'
                    showEstimation = bioinfoprivate.opttf(pval);
                    if isempty(showEstimation)
                        error(message('bioinfo:mspalign:showEstimatioNotLogical', upper( char( okargs( k ) ) )));
                    end 
            end
        end
    end
end

% put all P into a single array
P = cell2mat(p);

% First find a suitable threshold in the intensity values
res_int = 10000;
th = myQuantile(P,quan_int,res_int);

% limits of the MZ vector
mzmi = min(P(:,1));
mzma = max(P(:,1));
mz2idx =  @(x,r) round((x - mzmi) / (mzma-mzmi) * (r-1) + .5);
idx2mz = @(x,r) (x-.5) / (r-1) * (mzma-mzmi) + mzmi;

switch estimationMethod
    case 'histogram'        
        % Count all peaks above threshold and make a corse histogram just below
        % the resolution of the data
        res_mz = round((mzma-mzmi)/histResolution)+1;
        mzva = accumarray(mz2idx(P(:,1),res_mz),P(:,2)>th,[res_mz 1]);
        maxcount = max(mzva);
        
        if showEstimation
            mzva_bak = mzva;
        end

        % Smooth it and detect peaks, these will be the common masses
        K = 50; %increase resolution to improve peak centroid accuracy
        mzva = bioinfoprivate.masmooth(idx2mz(1:res_mz,res_mz),mzva,15,'loess');
        mzva = interp1(idx2mz(1:res_mz,res_mz),mzva,idx2mz(1:res_mz*K,res_mz*K),'pchip');
        cmz = idx2mz((find(mzva(3:end)<mzva(2:end-1) & mzva(2:end-1)>mzva(1:end-2))+1),res_mz*K);
       
        % correct cmz for long empty spaces
        dcmz = diff(cmz);
        h = max(0,round(dcmz./median(dcmz(:),1))-1);
        cmz = [cell2mat(arrayfun(@(q) cmz(q)+(0:h(q))*dcmz(q)/(h(q)+1),1:numel(h),'Uniform',false)) cmz(end)];
        
        if showEstimation
            figure
            hl3 = plot(cmz(ceil(1/3:1/3:numel(cmz))),...
                       repmat([min(mzva) maxcount NaN],1,numel(cmz)),'r:',...
                       'DisplayName','Common Mass/Charge','Tag','mzmarker');
            hold on
            hl2 = plot(idx2mz(1:res_mz*K,res_mz*K),mzva,'-','Color',[.3 .3 1],...
                       'DisplayName','Smoothed Counts','Tag','mzcounts');
            hl1 = plot(idx2mz(1:res_mz,res_mz),mzva_bak,'.','Color',[0 .5 0],...
                  'MarkerSize',10,'DisplayName','Peak Counts','Tag','mzcounts');
            title('Estimated CMZ Vector by Histogram Method');
            xlabel('Mass/Charge (M/Z)')
            ylabel('Counts')
            legend([hl3,hl2,hl1])
            axis([idx2mz(1,res_mz*K) idx2mz(res_mz*K,res_mz*K) min(mzva) maxcount])
            grid on
            hold off
            setAllowAxesRotate(rotate3d(gcf),gca,false)
            msDataCursor(gcf)
        end

    case 'regression'

        % Find all peaks above threshold and which are contiguous:
        h = find(P(:,2)>th);   
        h = h(diff(h)==1);
        Q = [P(h,1) P(h+1,1)-P(h,1)]; % [Peak loc, Peak distance]
        Q(Q(:,2)<=0,:)=[]; % remove observations from diff spectra
        
        if size(Q,1)>10000
            Q = Q(randsample(size(Q,1),10000),:);
        elseif size(Q,1)<1000
            error(message('bioinfo:mspalign:TooFewSamples'))
        end
        
        if showEstimation
            Q_bak = Q;
        end
        
        % regress a smooth curve
        Q = sortrows(Q);
        Q(:,2) = bioinfoprivate.masmooth(Q(:,1),Q(:,2),.5);
        
        % remove repeated observations in mz (so we can use interp1 later)
        Q(~diff(Q(:,1)),:)=[];
        
        % create the cmz vector
        cmz = zeros(1,round((mzma-mzmi)/min(Q(:,2))));
        cmz(1) = mzmi;
        i = 1;
        while cmz(i)<=mzma
            cmz(i+1) = cmz(i)+interp1(Q(:,1),Q(:,2),cmz(i),'pchip','extrap');
            i = i+1;
        end
        cmz(i:end)=[]; 
        
        if showEstimation
            figure
            hl3 = plot(cmz(ceil(1/3:1/3:numel(cmz))),...
                       repmat([min(Q_bak(:,2)) max(Q_bak(:,2)) NaN],1,numel(cmz)),'r:',...
                       'DisplayName','Common Mass/Charge','Tag','mzmarker');
            hold on
            hl1 = plot(Q_bak(:,1),Q_bak(:,2),'.','Color',[0 .5 0],...
                  'MarkerSize',10,'DisplayName','Observed Distance','Tag','mzdistance');
            hl2 = plot(Q(:,1),Q(:,2),'-','Color',[.3 .3 1],...
                       'DisplayName','Smoothed Distance','Tag','mzdistance');
            title('Estimated CMZ Vector by Regression Method');
            xlabel('Mass/Charge (M/Z)')
            ylabel('Inter Peak Distance')
            legend([hl3,hl2,hl1])
            axis([min(Q_bak(:,1)) max(Q_bak(:,1)) min(Q_bak(:,2)) max(Q_bak(:,2))])
            grid on
            hold off
            setAllowAxesRotate(rotate3d(gcf),gca,false)
            msDataCursor(gcf)
        end
end

if nargout<2
    if nargout == 0
        clear cmz
    end
    return
end

po = cell(size(p));
switch correctionMethod
    case 'nearest-neighbor'
        lutSize = 1000000; 
        lut = zeros(lutSize,1);
        lut([1,mz2idx(cmz(1:end-1)+(cmz(2:end)-cmz(1:end-1))/2,1000000)])=1;
        lut = cumsum(lut);
        for i = 1:numel(p)
            pks = p{i};
            pks(:,1) = cmz(lut(mz2idx(pks(:,1),1000000)))';
    
            % remove peaks assigned to the same spot
            pks(find(~diff(pks(:,1))),2) = max([pks(find(~diff(pks(:,1))),2), pks(find(~diff(pks(:,1)))+1,2)],[],2); %#ok
            pks(find(~diff(pks(:,1)))+1,:) = [];
    
            po{i} = pks;
        end
    case 'shortest-path'
        rpks = cmz'*[1 0];
        rpks(:,2) = 1;
        band = median(diff(cmz));
        for i = 1:numel(p)
            pks = sortrows(p{i});
            [j,k] = samplealign(double(rpks),double(bsxfun(@times,pks,1./[1 max(pks(:,2))])),'band',band);           
            po{i} = [rpks(j,1) pks(k,2)];
        end
end


function th = myQuantile(P,quan_int,res_int)
% computes a quick quantile using accumarray
imi = min(P(:,2));
ima = max(P(:,2));
in2idx =  @(x,r) round((x - imi) / (ima-imi) * (r-1) + .5);
idx2in = @(x,r) (x-.5) / (r-1) * (ima-imi) + imi;
inva = accumarray(in2idx(P(:,2),res_int),P(:,1),[res_int 1]);
th = idx2in(interp1q(cumsum(inva)/sum(inva),(1:res_int)',quan_int),res_int);
