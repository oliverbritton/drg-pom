% PlotExampleDRGActionPotentials
% 29/09/2015

directory = 'E:\CLPC48\Neuron Project\Data\AP Recordings\Preliminary';

cd(directory)

files = dir('*.csv');

firstEvokedAP = [6 4 5 3 3 6 3 5];
fig
for i = 1:length(files)
    
    data = importdata(files(i).name,',',2);
%     fig;
    numRecordings = size((data.data),2)-1;
    
    for j = firstEvokedAP(i)%1:numRecordings
%        subplot(3,3,j)
       hold on
       plot(data.data(:,1),data.data(:,j+1)*1000)
       title(data.colheaders{j+1})
       ylim([-100 100])
    end
end