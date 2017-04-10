function [dataMat, idxUPCEnter, idxUPCExit, idxUPCEE] = ProcessYogurtData()

%% Load from SQL Server
conn=database.ODBCConnection('SQLSERVER_yx10','','');
curs = exec(conn, 'select COLUMN_NAME from INFORMATION_SCHEMA.COLUMNS where TABLE_NAME = ''yogurt_all_info'''); 
setdbprefs('DataReturnFormat','cellarray');
curs = fetch(curs);
headline = curs.Data'; % UPC,UPCDesc,EAN,flavor,LAUNCH_YEAR,aisle,brand,deptid,manufacturer,parent,
                        % product,UPC_2,WEEK_MOVED,num_nfp,num_range,num_claims,num_nutr_claims,num_flavor,
                        % num_other,num_style,num_type,Flag_dup,brand_type,style,style_codes,other,
                        % other_codes,type,type_codes,tot_units,tot_volume,flag_impute,category,
                        % freq_08,freq_09,freq_10,freq_11,freq_12,fsum,exit_year
clear curs;

curs = exec(conn, 'select * from yogurt_all_info'); 
setdbprefs('DataReturnFormat','cellarray');
curs = fetch(curs);
data = curs.Data;
clear curs;

%% Extract the numeric columns of interest
numCol = length(headline);
numRow = size(data, 1);

% Columns we need
columns = {'UPC', 'LAUNCH_YEAR', 'exit_year', 'freq_08','freq_09','freq_10','freq_11','freq_12'};
dataMat = zeros(numRow, length(columns));
% find their indices in the raw data
for k = 1 : length(columns)
    for j = 1 : numCol
        if strcmpi(headline{j}, columns{k}) == 1
            for i = 1 : numRow
                if k <= 2 % UPC and launch_year needs to be converted into number first
                    dataMat(i, k) = str2num(data{i, j});
                else
                    dataMat(i, k) = data{i, j};
                end
            end
            break;
        end
    end
end


% sum over 2008-2010
fsum_080910 = sum(dataMat(:, 4 : 6), 2);
dataMat = [dataMat, fsum_080910];

% remove sum = 0
dataMat = dataMat(fsum_080910 >= 100, :);
fsum_080910 = dataMat(:, end);

% frequency statistics 
xt = 0:10:max(fsum_080910);
yt = hist(fsum_080910, xt);
figure
bar(xt, yt);
xlabel('frequency');
ylabel('number of UPCs');

% entering market during 2008-2010
idxUPCEnter = (dataMat(:, 2) >= 2008) .* (dataMat(:, 2) <= 2010);
numUPCEnter = sum(idxUPCEnter);
fprintf(sprintf('%d UPCs entered market during 2008 - 2010\n', numUPCEnter));

% exitting market during 2008-2010
idxUPCExit = (dataMat(:, 3) >= 2008) .* (dataMat(:, 3) <= 2010);
numUPCExit = sum(idxUPCExit);
fprintf(sprintf('%d UPCs exited market during 2008 - 2010\n', numUPCExit));

% entering AND exitting market during 2008-2010
idxUPCEE = (dataMat(:, 2) >= 2008) .* (dataMat(:, 2) <= 2010) .* ...
    (dataMat(:, 3) >= 2008) .* (dataMat(:, 3) <= 2010) .* ...
    (dataMat(:, 2) <= dataMat(:, 3));
numUPCEE = sum(idxUPCEE);
fprintf(sprintf('%d UPCs entered AND exited market during 2008 - 2010\n', numUPCEE));

