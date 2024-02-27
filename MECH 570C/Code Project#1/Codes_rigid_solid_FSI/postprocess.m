clc
clear all
close all
filename = 'rigidV3.oisd';
fileID = fopen(filename, 'r');

if fileID == -1
    error('File cannot be opened');
end

dragCoefficients = [];
liftCoefficients = [];

while ~feof(fileID)
    line = fgetl(fileID); 
    if contains(line, 'Drag and Lift coefficient')
        coeffsLine = fgetl(fileID);
        coeffs = textscan(coeffsLine, '%f %f');
        dragCoefficients = [dragCoefficients; coeffs{2}];
        liftCoefficients = [liftCoefficients; coeffs{1}];
    end
end


fclose(fileID);

timeSteps = 1:length(dragCoefficients); 

figure;
plot(timeSteps, dragCoefficients, '-o', 'MarkerSize', 4, 'MarkerFaceColor', 'blue');
legend('Drag Coefficient');
xlabel('Time Step');
ylabel('Drag Coefficient');
title('Drag Coefficient vs Time Step');

figure;
plot(timeSteps, liftCoefficients, '-x', 'MarkerSize', 4, 'MarkerEdgeColor', 'red');
legend('Lift Coefficient');
xlabel('Time Step');
ylabel('Lift Coefficient');
title('Lift Coefficient vs Time Step');


meanDrag = mean(dragCoefficients(1:500));
disp(['Mean Drag Coefficient: ', num2str(meanDrag)]);


maxAbsLift = max(abs(liftCoefficients(1:500)));
disp(['Max Absolute Lift Coefficient: ', num2str(maxAbsLift)]);

meanLift =mean(liftCoefficients(1:500));
disp(['Mean Lift Coefficient: ', num2str(meanLift)]);

rmsdrag = sqrt(sum((dragCoefficients(1:500,1)-meanLift).^2)/500);
disp(['RMS of Drag Coefficient: ', num2str(rmsdrag)]);

rmsLift = sqrt(sum((liftCoefficients(1:500,1)-meanLift).^2)/500);
disp(['RMS of Lift Coefficient: ', num2str(rmsLift)]);

