% 08/02/2016

fig('width',25,'height',25)

for i = 0:99
    
    subplot(10,10,i+1)
    data = importdata(['E:\CLPC48\Neuron Project\Simulations\Techlab\Test\FirstTest_' num2str(i) '.dat']);
    plot(data(:,1),data(:,2))
end
    
    