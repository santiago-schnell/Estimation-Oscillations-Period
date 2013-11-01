function [res]=autocorrelation(xdata,ydata,noscip,maxperm,func,plts,fign)
% ----------
% func is 'perms' or 'normal'
% maxperm is the number of permutations for the calculation of the null distribution 

n=0; % Plot counter
res=zeros(1,noscip*3); % Final Results

% ----------

n=n+1;
if ~isempty(find(plts==n,1))
    figure(fign);
    subplot(2,2,1);
    plot(xdata,ydata,'k');
    set(gca,'fontweight','b','fontsize',16);
    xlabel('Time (min)','fontweight','b','fontsize',16);
    ylabel('Ca^{2+} ','fontweight','b','fontsize',16);
    grid on;
end

nLags=floor(length(ydata)/2);
[ACF, ~, ~]=autocorr(ydata,nLags);
ACF=ACF(2:end);
xaxis=(1:nLags)*(xdata(2)-xdata(1));

n=n+1;
if ~isempty(find(plts==n,1))
    figure(fign);
    subplot(2,2,2);
    plot(xaxis,ACF,'k');
    set(gca,'fontweight','b','fontsize',16);
    xlabel('Lag (min)','fontweight','b','fontsize',16);
    ylabel('ACF','fontweight','b','fontsize',16);
    grid on;
end

% Find the extrema of the ACF function
[ymax,imax,xxx,yyy]=extrema(ACF);

if ~isempty(ymax)
    % Only consider ACF extrema values higher than 0.1 the maximum ACF value
    cand=find(ymax>(max(ACF)*0.1));
    ymax=ymax(cand);
    imax=imax(cand);
    
    % Remove the first point of ACF (always unity)
    cand=find(imax==1);
    imax(cand)=[];
    ymax(cand)=[];
    
    % Sort by ascending period
    [imax order]=sort(imax);
    ymax=ymax(order);
      
    % Remove multiples of periods
    l=1;
    while l<length(imax)
        posits=~mod(imax,imax(l))&imax~=imax(l);
        imax(posits)=[];
        ymax(posits)=[];
        l=l+1;
    end
    
    % Filter results by the number of oscillation ACF values we are interested in
    if (length(ymax)>noscip)
        ymax=ymax(1:noscip);
        imax=imax(1:noscip);
    end
end

if ~isempty(ymax)
    
    if ~isempty(find(plts==2,1))
        figure(fign);
        subplot(2,2,2);
        hold on;
        plot(xaxis(imax),ymax,'ko');
    end
    
    % Obtain the random permutations
    ydata_aux=zeros(maxperm,length(xdata));
    for m=1:maxperm
        if strcmp(func,'perms')
            perm=randperm(length(xdata));
            ydata_aux(m,:)=ydata(perm);
        elseif strcmp(func,'normal')
            ydata_aux(m,:)=mean(ydata)+std(ydata).*randn(1,length(xdata));
        end
    end
    
    % Obtain the ACF from random permutations
    ACF_tt=zeros(length(imax),maxperm);
    for m=1:maxperm
        [ACF_aux,~,~]=autocorr(ydata_aux(m,:),floor(length(ydata_aux(m,:))/2));
        ACF_aux=ACF_aux(2:end);
        for l=1:length(imax)
            ACF_tt(l,m)=ACF_aux(imax(l));
        end
    end
    
    % Calculate the significance
    sigpos=[];
    for l=1:length(imax)
        aux=length(find(ACF_tt(l,:)>ymax(l)));
        pvalue=aux/maxperm;
        index=3*(l-1)+1;
        if (pvalue < 0.05)   % 5% Significance level
            sigpos=[sigpos l];
            res(index)=xdata(imax(l)+1)-xdata(1); % stores the period
        else
            res(index)=(xdata(imax(l)+1)-xdata(1))*-1;
        end
        res(index+1)=ymax(l); % stores the power
        res(index+2)=pvalue; % stores the pvalue
    end
    
    if ~isempty(find(plts==2,1))
        sprintf('Significant Periods according to random data: ');
        for m=1:length(sigpos)
            figure(fign);
            subplot(2,2,2);
            hold on;
            plot(xaxis(imax(sigpos(m))),ymax(sigpos(m)),'k*');
        end              
    end
    
end
  
end