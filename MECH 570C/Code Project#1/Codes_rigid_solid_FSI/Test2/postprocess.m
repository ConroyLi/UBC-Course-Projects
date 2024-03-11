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

timeSteps = 1:length(dragCoefficients(500:end)); 

figure;
plot(timeSteps, dragCoefficients(500:end), '-o', 'MarkerSize', 4, 'MarkerFaceColor', 'blue');
legend('Drag Coefficient');
xlabel('Time Step');
ylabel('Drag Coefficient');
title('Drag Coefficient vs Time Step');

figure;
plot(timeSteps, liftCoefficients(500:end), '-x', 'MarkerSize', 4, 'MarkerEdgeColor', 'red');
legend('Lift Coefficient');
xlabel('Time Step');
ylabel('Lift Coefficient');
title('Lift Coefficient vs Time Step');


meanDrag = mean(dragCoefficients(500:end));
disp(['Mean Drag Coefficient: ', num2str(meanDrag)]);


maxAbsLift = max(abs(liftCoefficients));
disp(['Max Absolute Lift Coefficient: ', num2str(maxAbsLift)]);

meanLift =mean(liftCoefficients(500:end));
disp(['Mean Lift Coefficient: ', num2str(meanLift)]);

rmsdrag = sqrt(mean((dragCoefficients/2-meanDrag).^2));
disp(['RMS of Drag Coefficient: ', num2str(rmsdrag)]);

rmsLift = sqrt(mean((liftCoefficients(500:end)/2-meanLift).^2));
disp(['RMS of Lift Coefficient: ', num2str(rmsLift)]);

