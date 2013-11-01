function [xdata,ydata]=gen_input(amps,pers,nperiods,npperiod,noise)
% This function generates the input data
ydata=0;
for i=1:length(pers)
    freq=1/pers(i); % Frequency
    xdata=1:min(pers)/npperiod:(max(pers)*nperiods)+1;
    ydata=ydata+amps(i)*sin(2*pi*xdata*freq); % Create a sine wave...
end
ydata=ydata+noise.*randn(1,length(xdata));
end
% ----------
