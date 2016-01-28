path_data = datapath;

xlsfile = 'Metadata.xlsx';
path_metadata = fullfile(path_data,xlsfile);
[~,sheets] = xlsfinfo(path_metadata);

X = xlsreadmetadata(path_metadata,{'pp1'},1,'Label',[1],'Condition',[1 11]);
