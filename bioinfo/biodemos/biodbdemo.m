%% Connecting to Local Databases 
% This example shows how to collect sequence data from GenBank(R) and
% stores it in a Microsoft(R) Access(TM) database, using the
% *Bioinformatics Toolbox(TM)*. Many of the computations done in MATLAB(R)
% on biological data can be stored in an organized queriable manner. The
% example uses SQL (Structured Query Language) at the MATLAB command line
% but you can do this exercise using the *Visual Query Builder* graphical
% interface from the *Database Toolbox(TM)*. NOTE: this example uses
% Access(TM) but the database toolbox can interface with many different
% databases (MySQL(R), Excel(R), Oracle(R), Postgres).

% Copyright 2005-2012 The MathWorks, Inc.

%% Check for Database Toolbox(TM)
% Check to see if the Database Toolbox is installed

dbTbxVer = ver('database');
if isempty(dbTbxVer)
    disp(sprintf('This example requires the Database Toolbox.'));
    return
end

%% Make a Backup Copy of the Original Database
% This example modifies the file *sequences.mdb*. If you do not want to
% change the original version of this file, you should make a copy of it
% before continuing with this example.

%% Configure Database for Use with MATLAB(R)
% In order for MATLAB to interface with your database, a driver must be
% installed on your machine. To do this, consult the Database Toolbox
% documentation, "Setting Up a Data Source for ODBC Drivers."
% Create a new user DSN that points to the *sequences.mdb* database and 
% name it *FASTAdatabase*. This allows the *Database Toolbox* to connect
% to the Microsoft Access database.

doc database;

%% Connect to the Database
% Connect to the database using the DSN name *FASTAdatabase*. No username
% or password is required for this example database. The variable
% *connection* is used for the rest of the example to interact with the
% database. This connection will be terminated at the end of the exercise.

connection = database('FASTAdatabase','','');

%% Get Information About This Database
% The metadata for the database contains information about tables and field
% names in the database. First create a metadata object (*dbmeta*) using
% the *dmd* function and get the information on the database. The *tables*
% function returns all the tables that are in the database. There are
% system tables listed, but they are not used in this example and should
% generally only be used by an experienced database manager. To get the
% field names from the fasta table, use the function *columns*.

dbmeta = dmd(connection); 
mydbmeta  = get(dbmeta); 
databasetables = tables(dbmeta, mydbmeta.Catalogs,mydbmeta.Schemas) 
fastacolumns = columns(dbmeta, mydbmeta.Catalogs, mydbmeta.Schemas,'fasta')

%% Collect Sequence Data from Genbank and Insert it into the Database
% The cell array acc holds a list of accession numbers for BRCA genes (the 
% genes related to breast cancer). To enter the fasta information into the
% database from GenBank use the *getgenpept* function to download the
% protein sequences in fasta format. The data is  inserted into the 
% database with the *insert* function from the *Database Toolbox*. 

acc = {'NP_036646', 'NP_033894', 'AAG43492', 'AAP12647', 'NP_001013434', ...
            'AAK71628','AAN77220', 'AAC39589', 'AAC39584', 'AAC39585'};
        
for i = 1:length(acc)
    % Download the sequences
   brca = getgenpept(acc{i},'FILEFORMAT','FASTA'); 
    % Insert them into the database
   insert(connection,'fasta',fastacolumns,{brca.Sequence,brca.Header,acc{i}})
end

%% Verify That the Sequences Were Entered into the Database
% To ensure that the data was entered into the database, you can use a SQL
% select statement on the fasta table to view its contents. The variable 
% *cursor* holds information about the transaction with the database 
% including the sequence information. You can use the *fetch* and *exec* 
% functions to get the information into MATLAB. You can also check the 
% table's contents by seeing how many rows are returned. To view a sequence
% in fasta format, concatenate the header information above the sequence, 
% using the *char* function and the *seqdisp* function. 

cursor = fetch(exec(connection, ...
    'select fasta_header, fasta_sequence from fasta'));
 % The number of data points that match out query
numrows = rows(cursor) 
 % Display a sequence as like a FASTA file entry
char(['>' cursor.Data{1,1}], ...
    seqdisp(cursor.Data{1,2},'SHOWNUMBERS',false,'COLUMN',60))

%% Update Data in the Database 
% After information is entered in the database, the *update* function is 
% used to update or modify the data. In this case MATLAB determines 
% the sequence length of a sequence retrieved from the database and stores 
% the length in the database in the *length* field. You can use the *length*
% function to obtain sequence length.

for i = 1:length(acc)
  cursor = fetch(exec(connection, ...
    ['select fasta_sequence from fasta where acc_number = ''' acc{i} '''']));
  seq_length = length(cursor.Data{1,1});
   % Update the length field in the database
  update(connection,'fasta',{'length'},{seq_length},['where acc_number = ''' acc{i} '''']); 
end

%% Add Alignment Information to the Database
% Instead of aligning sequences repeatedly, an alignment can be stored in
% the database for later. This reduces the number of unnecessary
% alignments. As an example a local alignment of the Rat and Mouse BRCA gene
% is stored in the database. Information on the alignment table fields
% from the database is obtained using the metadata object created
% earlier (*dbmeta*). Each alignment is separated into three strings of 
% characters. This process could easily be automated for multiple sequences.

aligncolumns = columns(dbmeta,mydbmeta.Catalogs,mydbmeta.Schemas,'alignment');
% Retrieve the sequences
cursor = fetch(exec(connection,...
 'select fasta_sequence from fasta where acc_number in (''NP_036646'',''NP_033894'')'));
% Align the two sequences
[scr aln] = swalign(deblank(cursor.Data{1}),deblank(cursor.Data{2}));
% Insert the alignment
insert(connection,'alignment',aligncolumns, ...
    {aln(2,:),'NP_036646','NP_033894',scr,aln(1,:),aln(3,:)});
  
%% Retrieve an Alignment
% To retrieve an alignment from the database, simply select the alignment
% information from the alignment table. The *showalignment* function is 
% used to view the alignment. 

cursor = fetch(exec(connection,...
  'select align_seq_1, alignment, align_seq_2 from alignment where sequence_1_acc = ''NP_036646'' and sequence_2_acc = ''NP_033894'''));
showalignment([cursor.data{1,1}; cursor.data{1,2}; cursor.data{1,3}])
    
%% Add BLAST Report Information to the Database
% Saving a BLAST report in a queriable database can be an advantage
% when there are numerous sequences to analyze. The *blastncbi* function 
% sends sequences from the database to NCBI and BLASTs them. The resulting
% information is parsed by the *getblast* function into sections that can be
% stored in the database.

cursor = fetch(exec(connection, ...
 ['select fasta_sequence, fasta_header from fasta where acc_number = ''' acc{10} '''']));
    % Put the sequence into a MATLAB struct
    seq.Sequence = cursor.Data{1,1};
    seq.Header =   cursor.Data{1,2};
    
%% BLAST the Sequence    
RID = blastncbi(seq, 'blastp');

%% Wait for NCBI to Finish BLASTing
blastreport = [];
while(isempty(blastreport))
    try
        blastreport = getblast(RID);
        found = true;
    catch theException
        [match t] = regexp(theException.message,'BLAST.* (\d)+\sseconds','match','tokens');
        if ~isempty(match)
            disp(sprintf('%s', match{1}));
            pause(str2num(cell2mat(t{1})));
        else
            disp('Waiting for BLAST to respond. Hit Ctrl-C to abort.');
            pause(20);
        end
    end
end

blastcols = columns(dbmeta, mydbmeta.Catalogs, mydbmeta.Schemas, 'blast');
hitscols = columns(dbmeta, mydbmeta.Catalogs, mydbmeta.Schemas, 'hits');

%% Insert the BLAST Results into the Database
insert(connection,'blast',blastcols, ...
    {acc{10},blastreport.RID,blastreport.Algorithm,blastreport.Database});

 for i = 1:size(blastreport.Hits,2)
  insert(connection,'hits',hitscols, ...
    {blastreport.RID,blastreport.Hits(i).Name,blastreport.Hits(i).Length});
 end
  
%% Put Information into MATLAB Using the Visual Query Builder 
% The Database Toolbox has a visual query builder to query the database and
% put information into MATLAB. This tool makes interfacing with databases
% easier for people who are not comfortable with SQL.  If the ODBC drivers 
% are working properly, the DSN that you created early should appear in the 
% Visual Query Builder.

querybuilder

%% Empty the Database and Close All Connections
% After the desired information is collected from the database, you can 
% close all connections to the database. Empty any information that you 
% entered into the database.

% Find all of the non-System tables 
index = strcmp(databasetables(:,2),'TABLE');
tablestoempty = databasetables(index,1);
for i = 1:length(tablestoempty)
    % use delete to empty them without changing the structure of the DB
cursor = exec(connection,['delete * from ' tablestoempty{i}]);
end

% close our connection with the database
close(connection);

%%
% *<mailto:bioinfo-feedback@mathworks.com?subject=Feedback%20for%20BIODBDEMO%20in%20Bioinformatics%20Toolbox%204.2 Suggest an enhancement for this example.>*

displayEndOfDemoMessage(mfilename)
