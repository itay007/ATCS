%% Reconstructing the Origin and the Diffusion of the SARS Epidemic
% This example shows an analysis of the origin and diffusion of the SARS
% epidemic. It is based on the discussion of viral phylogeny
% presented in Chapter 7 of "Introduction to Computational Genomics. A Case
% Studies Approach" [1]. 

%   Copyright 2007-2012 The MathWorks, Inc.


%% Introduction
% SARS (Severe Acute Respiratory Syndrome) is a recently-emerged disease
% caused by a new type of coronavirus (SARS-CoV). It consists of a 29,571
% base-long, single-stranded RNA and displays a characteristic spiky
% envelope protein that resembles a crown. 
%
% The first cases of SARS appeared in late 2002 in the Chinese province of
% Guangdong and grew into a major outbreak in the next few months
% (January through February 2003). The majority of the infected individuals
% acquired the disease in the Guangzhou Hospital. A doctor who had
% worked in this hospital traveled to Hong Kong in February 2003 and
% stayed at the Metropole Hotel. The doctor and a number of other hotel
% guests all became infected with the virus and traveled to different
% destinations (Vietnam, Canada, Singapore, Taiwan) carrying the disease and the
% virus.
%
% By analyzing the phylogenetic relationships between the samples of
% SARS-CoV that were collected in late 2002 and in 2003, we can reconstruct
% the history of the SARS epidemic and understand how it spread throughout
% the world in such a short period of time.

%% Loading the Sequence Data of SARS Strains
% We consider the nucleotide sequences of 13 strains of human SARS
% coronaviruses for which the location and the date of collection are
% known. The sequences correspond to the spike S protein, which is
% responsible for binding to specific receptors and is considered a major
% antigenic determinant. Because the Himalayan palm civet is believed to be
% the source of the human SARS-CoV, we also consider the sample derived
% from the palm civet. For the sake of convenience, the sequence data is
% stored in a MATLAB(R) structure called |spike| consisting of |Header| and
% |Sequence| fields for each viral strain. The data can also be downloaded
% from the GenBank(R) database using the accession numbers stored in the
% structure |accNum|. 
%

% Load the genomic data for the human and palm civet SARS strains
load sarsdata.mat

% Display the accession numbers and collection dates of the sequence
% dataset. 
disp('Genbank Acc. No.    Collection Date       Location')
disp(sprintf('%10s %22s %25s\n',accNum{reshape(1:42,14,3)'}))

%% Computing the Sequence Pair-Wise Distances 
% Obtain the distance matrix needed to build the phylogenetic tree by
% computing a symmetric matrix that holds pair-wise genetic distances with
% Jukes-Cantor corrections. Ignore sequence sites representing gaps. 

JC_distances = seqpdist(spike,'method','jukes-cantor','alphabet','NT', ...
                           'indels','pairwise-delete','squareform',true); 
numSeq=size(JC_distances,1);

%%
% By plotting the distance matrix, we can appreciate the presence of a
% subset of sequences that are more closely related to each other (central
% cluster, represented by the darker tones). The last sequence, which is
% associated to the Himalayan palm civet, is the most distant to the
% majority of the members of the set. This is expected because it is a
% nonhuman coronavirus.

figure;
imagesc(JC_distances);
colormap(bone);
colorbar; 
title ('Pair-wise distances (spike protein nt sequences)');

%% Constructing a Neighbor-Joining Phylogenetic Tree
% Using the distances computed above, construct a phylogenetic tree
% using the neighbor-joining method. In this case, we assume equal
% variance and independence of evolutionary distance estimates.

tree1 = seqneighjoin(JC_distances,'equivar',spike);
plot(tree1,'orient','left');
title('Neighbor-joining tree using Jukes-Cantor model');

%%
% The tree depicts the story of the epidemic. The early infections all
% occurred in Guangzhou and Zhongshan, labelled as GZ and ZS respectively.
% The international cases (Hong Kong, Singapore, Hanoi, Taiwan, Toronto)
% are all related to each other and seem to branch from the case traced
% back to the Metropole Hotel in Hong Kong.

%% Estimating the Date of Origin of the Epidemic
% Because the date of collection of each SARS strain is known, we can
% observe the progression of the virus mutations over time. Consider the
% pair-wise distances according to the Kimura model, which distinguishes
% between transitional and tranversional mutation rates. Then, restrict
% your analysis to the distance of each human strain from the
% Himalayan palm civet's strain. Finally, plot the genetic distance versus the
% date of collection.

K_distances=seqpdist(spike,'method','Kimura','squareform',true, ...
                           'alphabet','NT','indels','pairwise-delete');

% sequence of the palm civet
civet=find(~cellfun(@isempty, strfind({spike.Header}, 'civet')));

d = regexp({spike.Header},'\d+/\d+/\d+','match','once');
for i=1:numSeq-1
    % genetic distances with respect to the palm civet's strain
    scores(i,1)=K_distances(civet,i);
    % convert the collection dates into numbers
    dates(i,1)=datenum(d{i});
end

refDate=datenum('01/01/03','mm/dd/yy'); % reference date

figure();
plot(dates-refDate,scores,'k*');
ylabel('Genetic distance (relative to the palm civet)');
xlabel('Time distance from 01/01/03 (days)');
hold on;

%%
% In relation to the sequence of the palm civet, we observe that the
% genetic distance increases approximately in a linear manner with time.
% Perform a polynomial fitting and a least-square interpolation to outline
% the progression of the viral mutations over time and estimate the
% approximate date for the origin of the epidemic. The start of the
% infection corresponds more or less to the root of the polynomial fit,
% i.e., any date that is at zero genetic distance from the palm civet's sequence.

[P,S] = polyfit(dates-refDate,scores,1);
x=[-max(dates-refDate):.1:max(dates-refDate)];
[y,delta] = polyconf(P,x,S); % estimate 95% prediction interval

plot(x,y,'b-');
plot(x,y+delta,'r-',x,y-delta,'r-'); % confidence interval
line([-max(dates-refDate) max(dates-refDate)],[0 0],'LineStyle', ':');
title('Estimate of origin of SARS epidemic');


originDist=roots(P);% estimated distance between origin and reference date
estimated_origin=datestr(floor(originDist+refDate))
%dates(civet)=originDist+refDate;
plot(originDist, 0,'*b');
annotation(gcf,'textarrow', [0.245 0.245], [0.30 0.35], ...
           'String', {'estimated origin'}, 'color', [0 0 1]);

%% Rerooting the Phylogenetic Tree
% Because the disease caused by the novel strain of human SARS-CoV appears
% to have originated in the palm civet, we can assume that the location of
% the root for the human strains' phylogenetic tree is next to the node
% associated with the Himalayan palm civet.

civetNode = getbyname(tree1,'civet');
tree2 = reroot(tree1,civetNode,0);
plot(tree2,'orient','left');
title('Rerooted Neighbor-joining tree using Jukes-Cantor model');

%%
% The rerooted tree better illustrates the progression of the SARS epidemic.
% Starting with the early infections in the Guangdong province in 2002 (GZ
% 12/16/02 and ZS 12/26/02), the virus spread in the Guangzhou Hospital in
% early 2003 (GZ Hospital 01/31/03) and reached Hong Kong via the doctor
% who worked in the mentioned hospital and who stayed at the Metropole
% Hotel (Metropole 02/21/03). The virus was then carried across the borders
% via infected guests of the Metropole Hotel.

%% Observing the Phylogenetic Tree as It Builds
% Assuming that the samples represent the SARS coronavirus at different
% points in time, we can observe the virus evolution as the phylogenetic
% tree (built on the basis of genetic distances) is created. We can simulate
% the various steps in the tree reconstruction. The |movie| function
% animates the tree building process.

d = regexp({spike.Header},'\d+/\d+/\d+','match','once');
d{end} = datestr(estimated_origin, 'mm/dd/yy');
allDates = datenum(d);
[dummy,order] = sort(allDates); % sort according to collection date

for i = 2:numSeq
    sp = order(1:i);
    tr1 = seqneighjoin(JC_distances(sp,sp),'equivar',spike(sp));
    tr2 = reroot(tr1,getbyname(tr1,'civet'),0);
    h = plot(tr2,'leaflabels',true,'terminallabels',false);
    set(findobj(h.leafNodeLabels,'string',spike(sp(i)).Header),'Color','r')
    axis([-.0002 .0045 0 15])
    fs(i-1) = get(h.axes,'Parent');
    M(i-1) = getframe(fs(i-1));
end
close(fs) % close figures
%%

% movie(figure,M,1,1) % <== uncomment this line to play the animation

%% Visualizing the Diffusion of the Virus via a Directed Graph 
% We can also visualize the diffusion of the virus using a directed graph,
% where each node represents an infected individual and weights of edges
% are associated to the genetic distance between sequences. First, create
% an adjacency matrix based on the date of collection of the samples, such
% that possible paths run through nodes that are compatible in terms of the
% collection dates. Then, use the previously computed Jukes-Cantor
% distances to assign weights to the edges between nodes. And finally,
% determine the shortest path from the node associated with the palm civet
% and every other node.

 % adjacency matrix based on collection dates
gValid = bsxfun(@lt,allDates,allDates');
% weight matrix for the graph 
g1= sparse((gValid .* JC_distances));

% find shortest paths from civet node to all nodes
[dist,paths,pred_tree] = graphshortestpath(g1,civet); 
% create adjacency matrix for the winning shortest path
g2 = sparse(pred_tree(1:13),1:13,1,14,14).*g1;
% plot the graph 
spikeGraph=view(biograph(g2,{spike.Header}));

% nodes relative to Guangdong province (GZ and ZS)
guangdong=find((~cellfun(@isempty, strfind({spike.Header}, 'GZ'))) |...
    (~cellfun(@isempty, strfind({spike.Header}, 'ZS'))));
% node relative to the Metropole Hotel
metropole=find(~cellfun(@isempty, strfind({spike.Header}, 'Metropole')));
% node relative to the Guangzhou Hospital 
hospital=find(~cellfun(@isempty, strfind({spike.Header}, 'Hospital')));

% highlight some of the important nodes
set(spikeGraph.Nodes(civet),'Color',[1 1 1]) % white (palm civet)
set(spikeGraph.Nodes(guangdong),'Color',[1 .7 .7]) % pink (Guangdong)
set(spikeGraph.Nodes(metropole),'Color',[0.8 0.8 1]) % lavander (Metropole Hotel)
set(spikeGraph.Nodes(hospital),'Color',[1 0.3 0.3]) % red (GZ Hospital)

%% 
% This graph highlights the crucial role played by some of the infection
% events: 
%
% * The Himalayan palm civet is the source of the infection
% * The Metropole Hotel is the root of the branching for the international
% epidemic
% * The Guangzhou Hospital represents the bridge connecting the province of
% Guangdong (GZ and ZS) and the Metropole Hotel in Hong Kong.
%

%% References
%
% [1] Nello Cristianini and Matthew W. Hahn,  "Introduction to
% Computational Genomics. A Case Studies Approach", Cambridge University
% Press, 2007.
%
%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20SARSDEMO%20in%20Bioinformatics%20Toolbox%204.2 Provide feedback for this example.>*

displayEndOfDemoMessage(mfilename)
