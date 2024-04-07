% Base directory and pattern for the file names
baseDirectory = 'E:\UBC\MEng\Term2\MECH 503\Project\lammps_files\Jincong';
baseFilename = 'W_'; % Base name of the files
numFiles = 7; % Total number of files

% Initialize a cell array to store the elastic constants from all files
allElasticConstants = {};

% Loop over each file number
for fileNum = 2:numFiles+1
    % Construct the full path to the file
    currentFilePath = fullfile(baseDirectory, sprintf('%s%d.txt', baseFilename, fileNum));
    
    % Process the current file
    fileContent = fileread(currentFilePath);
    pattern = 'Elastic Constant (\w+) = ([\d\.-]+) GPa';
    matches = regexp(fileContent, pattern, 'tokens');
    
    % Temporary struct to hold constants from the current file
    elasticConstants = struct();
    
    % Extract constants
    for i = 1:length(matches)
        currentMatch = matches{i};
        constantName = currentMatch{1};
        constantValue = str2double(currentMatch{2});
        elasticConstants.(constantName) = constantValue;
    end
    
    % Store the extracted constants in the cell array
    allElasticConstants{end+1} = elasticConstants;
end

% Now `allElasticConstants` contains the data from all files

%% 
% 'up' values
up = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
hea = [1,2,3,4,5,6,7,8,9,10];

% Assuming the first struct has all field names we're interested in, to get a list of them
fieldNames = fieldnames(allElasticConstants{1});
numCValues = numel(fieldNames); % Number of C values

% Number of files should match the length of 'up'
if length(allElasticConstants) ~= length(hea)
    error('The number of files does not match the length of the "up" vector.');
end

% Preparing to plot - creating a new figure
fig = figure;
hold on; % Allows multiple plots in the same figure
grid on; % Enhances readability with a grid

% Colors or markers can be specified here to differentiate the plots
colors = lines(numCValues); % Generate a set of colors

% Loop over each C value
for i = 1:numCValues
    fieldName = fieldNames{i}; % Current C value's name
    cValues = zeros(1, length(allElasticConstants)); % Preallocate for efficiency
    
    % Extract this C value across all files
    for j = 1:length(allElasticConstants)
        cValues(j) = allElasticConstants{j}.(fieldName);
    end
    
    % Plotting this C value against 'up'
   
    loglog(hea, cValues, 'o-', 'Color', colors(i,:), 'DisplayName', fieldName);
end
set(gca, 'XScale');
set(fig, 'Position', [100, 100, 1200, 800]);
% Enhancing the plot
xlabel('Random Lattice Structure R1 to R10');
ylabel('C values (GPa)');
title('Components of the Elastic Constant Tensor vs. Up for HEA');
legend('show', 'Location', 'bestoutside'); % Show a legend for clarity
hold off; % No more plots in this figure

%%
% Initialize arrays to store the mean values and maximum percentage differences
meanCValues = zeros(1, numCValues);
maxPercentDiffCValues = zeros(1, numCValues);

% Loop over each C value
for i = 1:numCValues
    fieldName = fieldNames{i};
 % Current C value's name
    cValues = zeros(1, length(allElasticConstants)); % Preallocate for efficiency
    
    % Extract this C value across all files
    for j = 1:length(allElasticConstants)
        cValues(j) = allElasticConstants{j}.(fieldName);
    end
    
    % Calculate mean
    meanCValues(i) = mean(cValues);
    
    % Calculate the maximum percentage difference
    maxVal = max(cValues);
    minVal = min(cValues);
    maxPercentDiffCValues(i) = ((maxVal - minVal) / meanCValues(i)) * 100;
end

% Display the results
for i = 1:numCValues
    fprintf('Mean of %s: %f GPa\n', fieldNames{i}, meanCValues(i));
    fprintf('Max percentage difference of %s: %f%%\n', fieldNames{i}, maxPercentDiffCValues(i));
end
%%
% Extract field names (assuming all files have the same structure)
fieldNames = fieldnames(allElasticConstants{1});
numCValues = numel(fieldNames); % Number of C values

% Initialize matrix to store the C values for all files
cValueMatrix = zeros(numFiles, numCValues);

% Fill the matrix with C values
for fileIndex = 1:numFiles
    for cIndex = 1:numCValues
        cValueMatrix(fileIndex, cIndex) = allElasticConstants{fileIndex}.(fieldNames{cIndex});
    end
end

% Open a file to write the LaTeX code
fid = fopen('W.tex', 'w');

% LaTeX table header
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{|c%s|}\n', repmat('|c', 1, numCValues)); % Create enough 'c's for the columns
fprintf(fid, '\\hline\n');
fprintf(fid, 'File Number %s \\\\\n', strjoin(fieldNames, ' & ')); % Table headers
fprintf(fid, '\\hline\n');

% Loop over the files and create each row of the table
for fileIndex = 1:numFiles
    fprintf(fid, '%d ', fileIndex); % File number
    for cIndex = 1:numCValues
        fprintf(fid, '& %.2f ', cValueMatrix(fileIndex, cIndex)); % C value with 2 decimal places
    end
    fprintf(fid, '\\\\\n'); % End of row
    fprintf(fid, '\\hline\n');
end

% LaTeX table footer
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Elastic Constants for BCC W}\n');
fprintf(fid, '\\label{tab:elastic_constants}\n');
fprintf(fid, '\\end{table}\n');

% Close the file
fclose(fid);

% ...


