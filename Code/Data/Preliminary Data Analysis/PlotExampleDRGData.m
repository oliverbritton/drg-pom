% PlotExampleDRGData
% 06/08/2015

% Plot example voltage clamp experiments from Anabios
% TODO
% Standardise y limits


filenames{1} = 'EP2_2014-01-21_03_Activation.csv';
filenames{2} = 'EP2_2014-01-21_03_FastInactivation.csv';
filenames{3} = 'EP2_2014-01-21_03_SlowInactivation.csv';

names = {'Act' 'Fast Inact' 'Slow Inact'};


for i = 1:3
    
    experiment{i} = importdata(filenames{i},',',2);
    
    
end

for i = 1:3
    fig
    for j = 2:size(experiment{i}.data,2)
        
        hold on % Comment out for subplots
         subplot(6,4,j-1) % Comment out for overlay
        plot(experiment{i}.data(:,1),experiment{i}.data(:,j))
        if i == 1
            maxCurrent(j-1) = min(experiment{i}.data(:,j));
        end
        if j == 2
            title([experiment{i}.colheaders{j} ' mV' ' ' names{i}]);
        else
            title([experiment{i}.colheaders{j} ' mV']);
        end
    end
    
    
    
end

%% Plot I-V

% Activation
fig
V = -90:5:10
plot(V,-maxCurrent)

