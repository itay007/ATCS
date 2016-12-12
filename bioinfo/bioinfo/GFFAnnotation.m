classdef GFFAnnotation < SequenceAnnotation
%GFFANNOTATION Class representing a collection of GFF annotations.
%
%  A GFFAnnotation is a class representing a collection of annotations to
%  one or more reference sequences conforming with the data available from
%  a GFF formatted file. Annotations may be filtered by feature, by
%  reference or by sequence position. Then data pertaining to one
%  annotation or a subset of annotations can be pulled out into an struct 
%  array.
%
%  GFFAnnotation properties:
%  FieldNames - available data fields.
%  NumEntries - number of annotations in the object.
%
%  GFFAnnotation methods:
%  GFFAnnotation     - create a GFFAnnotation object.
%  getData           - retrieve data from the object as an struct array.
%  getIndex          - get an index to specific annotations. 
%  getSubset         - create another object containing a subset of the annotations.
%  getRange          - get the extent of the annotations in the sequence.
%  getReferenceNames - get the references available in the object.
%  getFeatureNames   - get the types of annotations available in the object.
%
%  See also SEQUENCEANNOTATION, GTFANNOTATION.

%   Copyright 2010-2012 The MathWorks, Inc.

%  GFFAnnotation Hidden methods:
%  getBrowserData           - retrieve data for the NGSBrowser.
%  getStructArray           - retrieve data as an struct array.
%  getRowInCompactAlignment - calculate a compact alignment.

    properties (GetAccess='public', SetAccess='private')
        % FIELDNAMES available data fields.        
        FieldNames = {'Reference','Start','Stop','Feature','Source','Score','Strand','Frame','Attributes'};
    end
    properties (GetAccess='protected', SetAccess='private')
        EfficientAccessFields = {'Reference','Feature','Start','Stop'};
        StringSearchableFields = {'Reference','Feature'};
    end
    properties (Dependent)
        Reference
        Start
        Stop
        Feature
        Source
        Score
        Strand
        Frame
        Attributes
    end
    
    methods
        %==================================================================
        % Constructor
        %==================================================================
        function obj = GFFAnnotation(in)
%GFFANNOTATION creates a GFFAnnotation object.
%
%  GFFANNOTATION(FILENAME) creates a GFFAnnotation object from a GTF or GFF
%  formatted file FILENAME.
%
%  Examples:
%
%  % Create a GFFAnnotation object from a GFF formatted file.
%  a = GFFAnnotation('tair8_1.gff')
%
%  % Create a GFFAnnotation object from a GTF formatted file.
%  b = GFFAnnotation('hum37_2_1M.gtf')
%
%  See also GTFANNOTATION.
           
            % Check if input is a structure with the proper fields
            if isstruct(in)
                % Check that the input structure has all the required
                % fields for this object (it may have more, that's fine,
                % the StructAnnotationDataAdapter will take care of return
                % the appropriate fields and in the correct order)
                if all(isfield(in,obj.FieldNames))
                    % Construct the StructAnnotationDataAdapter
                    obj.Data =  bioinfoprivate.StructAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields);
                else
                    error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidStruct'))
                end
                
            % Check if input is a dataset with the proper variables
            elseif isa(in,'dataset')
                if all(ismember(obj.FieldNames,get(in,'VarNames')))
                    obj.Data =  bioinfoprivate.DatasetAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields);
                else
                    error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidDataset'))
                end
            
            % Check if input is a dataset with sufficient number of fields
            % in the entries
            elseif isa(in,'BioIndexedFile')
                if strcmp(in.FileFormat,'table')
                    fieldBinding = [1 4 5 3 2 6 7 8 9];
                    numericFields = {'Start','Stop'};
                    obj.Data =  bioinfoprivate.BIFAnnotationDataAdapter(in,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields,fieldBinding,numericFields);             
                else
                    error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidFileFormatinBIF'))
                end
                
            % Check if input is a file
            elseif (ischar(in) && (isrow(in) || isempty(in))) && exist(in,'file')
                try
                   fid = fopen(in,'rt');
                   if fid==-1
                      error(message('bioinfo:GFFAnnotation:GFFAnnotation:CannotOpenFile', in))
                   end
                   feat = textscan(fid,'%s%s%s%d%d%s%s%s%s%*[^\n]','delimiter','\t','CommentStyle','#');
                   fclose(fid);
                   
                   if numel(unique(cellfun(@numel,feat)))>1
                      error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidFileFormat'))
                   end
                catch ME 
                   error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidFileFormat'))
                end
                fieldBinding = [1 4 5 3 2 6 7 8 9];
                obj.Data =  bioinfoprivate.CellAnnotationDataAdapter(feat,obj.FieldNames,obj.EfficientAccessFields,obj.StringSearchableFields,fieldBinding);             

            elseif ischar(in) && (isrow(in) || isempty(in))
                error(message('bioinfo:GFFAnnotation:GFFAnnotation:FileNotFound', in))
            else
                error(message('bioinfo:GFFAnnotation:GFFAnnotation:InvalidInput'))
            end
        end
        
    end
end
