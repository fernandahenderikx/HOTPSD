%% make single cell particle zenodo

clear all; 
files = dir('\\jett\awlab\DATA\IFCB CRUISE DATA\HOT\HOT*\proc\CNN_TS4_2025-02-26_Group2_perParticle\');

tax = readtable('\\jett\AWlab\DATA\Fernanda\IFCB\taxonomy_CNN_TS4_2025-02-26.xlsx');
u = unique(tax.GroupsPaper);
c = [1:length(u)];
groups_paper_list = nan(length(c),1);


% assign number to each row in the excel spreadsheet
for i = 1:length(c)
    ind = strmatch(u{i},tax.GroupsPaper);
    groups_paper_list(ind,1)=c(i);
end

% assign a number to each of the files to make it easier to group them
IDX = nan(length(files),1);
for i = 1:length(tax.GroupsPaper)

    ind = strmatch(string(tax.CNN(i)),string({files.name}));

    IDX(ind,1) = groups_paper_list(i);
end

% order ot T_all is same as order of u
clear T_all L_all data fname
for i = 1:length(c)
    ind = find(IDX==i); % 2 is bad; so will have to get rid of it at some point below
    % load each of the things that  Iwan to merge together, put all into
    % one row
    T=[]; L=[];
    for j = 1:length(ind)
        load(fullfile(files(ind(j)).folder,files(ind(j)).name));
        data = [ifcbTime lat lon cruise ifcbData.ESD_summed ifcbData.summedBiovolume ifcbData.MajorAxisLength ifcbData.MinorAxisLength];
        fname = ifcbData.Picture;
        % choice: remove all but station aloha??? for paper, i think so.
        xxx = find(lat>=22.6 & lat<=22.9 & lon>=-158.2 & lon <=-157.9);
        data=data(xxx,:);
        fname=fname(xxx,:);
        if ~isempty(data)
            T = [T; data];
            L=[L; fname];
        end
    end
    T_all{i,:} = T;
    L_all{i,:}=L;
    disp(i)
end



%% make cvs of single cell particles:
% Number of groups
u = unique(tax.GroupsPaper);

u{7,1}='cyanobacteria-crocosphaera-like';
u{19,1}='dictyochophyta';
u{20,1}='dinoflagellata-large';
u{21,1}='dinoflagellata-medium';
u{22,1}='dinoflagellata-small';
u{23,1}='dinoflagellata-unknown';
u{25,1}='haptophyta-chrysochromulina-like';
u{26,1}='haptophyta-emiliania-like';
u{34,1}='haptophyta-unknown';

ng = numel(T_all);

Tall = table();   % initialize output table
for i = 1:ng
    
    A = T_all{i};   % numeric matrix (n x 6)
    
    % convert to table
    T = array2table(A, ...
        'VariableNames',{'time','lat','lon','cruise','ESD','biovolume','major_axis','minor_axis'});
    
    % add group column
    T.group = repmat(string(u{i}), height(T), 1);

    % add file links
    imgNames = string(L_all{i});
    
    % split into bin name and image number
    parts   = split(imgNames, "_");
    binName = parts(:,1) + "_" + parts(:,2);
    imgNum  = string(str2double(parts(:,3)));   % removes leading zeros
    
    T.filelink = "https://ifcbdb2.soest.hawaii.edu/image?dataset=HOT" + ...
                 string(T.cruise) + ...
                 "&bin=" + binName + ...
                 "&image=" + imgNum;
    T.filelink = cellstr(T.filelink);

    % append
    Tall = [Tall; T];
    
    
end

ind = find(T.ESD==0);
Tall(ind,:)=[];

ind = strmatch('bad',Tall.group);
Tall(ind,:)=[];

% format:
Tall.time = datetime(Tall.time,'ConvertFrom','datenum');
Tall.lat = round(Tall.lat,4); Tall.lon = round(Tall.lon,4);
Tall.ESD = round(Tall.ESD,2);
Tall.major_axis = round(Tall.major_axis,2);
Tall.minor_axis = round(Tall.minor_axis,2);
Tall.biovolume = round(Tall.biovolume,4);


% write csv
writetable(Tall,'single_particle_IFCB_HOT294_to_HOT361.csv')


outfile = 'single_particle_IFCB_HOT294_to_HOT361.csv';

fid = fopen(outfile,'w');

% write custom header lines
fprintf(fid,'Single particle dataset from IFCB processing\n');
fprintf(fid,'Columns: time (image collection time, UTC), lat (deg N), lon (deg W), HOT cruise number, equivalent spherical diameter (from biovolume, micron), biovolume (micrometer^3), major axis length (micron), minor axis length (micron), annotation group (from Convolution Neural Network), dashboard link\n');
fprintf(fid,'F1-scores quantifying accuracy of particle annotations for the different groups are provided separately.\n');
fprintf(fid,'IFCB Dashboard is available here: https://ifcbdb2.soest.hawaii.edu/datasets\n');
fprintf(fid,'Generated on: %s\n', datestr(now));
fprintf(fid,'\n');   % blank line

fclose(fid);

% append the table after the header
writetable(Tall, outfile, 'WriteMode','append');
