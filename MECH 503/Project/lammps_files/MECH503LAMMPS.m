% Base directory and pattern for the file names
baseDirectory = 'E:\UBC\MEng\Term2\MECH 503\Project\lammps_files';
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

% Assuming the first struct has all field names we're interested in, to get a list of them
fieldNames = fieldnames(allElasticConstants{1});
numCValues = numel(fieldNames); % Number of C values

% Number of files should match the length of 'up'
if length(allElasticConstants) ~= length(up)
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
   
    loglog(up, cValues, 'o-', 'Color', colors(i,:), 'DisplayName', fieldName);
end
set(gca, 'XScale', 'log', 'XDir', 'reverse');
set(fig, 'Position', [100, 100, 1200, 800]);
% Enhancing the plot
xlabel('Up values');
ylabel('C values (GPa)');
title('Components of the Elastic Constant Tensor vs. Up for bcc W');
legend('show', 'Location', 'bestoutside'); % Show a legend for clarity
hold off; % No more plots in this figure


