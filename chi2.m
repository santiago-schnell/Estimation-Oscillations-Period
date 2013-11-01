function [res]=chi2(xdata,ydata,max_per,noscip,maxperm,func,plts,fign)
% ----------
% func is 'perms' or 'normal'
% maxperm is the number of permutations for the calculation of the null distribution

n=0; % Plot counter
chi=0; % Change this variable to 1 if you would want to make chi2 calculations
min_per=min(xdata); % Minimum period to evaluate
res=zeros(1,noscip*3); % Final Results

% ----------

n=n+1;
if ~isempty(find(plts==n,1))
    figure(fign);
    plot(xdata,ydata,'k');
    xlabel('Time');
    ylabel('Value');
end

minpoints=length(find(xdata<=min_per));
maxpoints=length(find(xdata<=max_per));
Qp=calc_Qp(xdata,ydata,minpoints,maxpoints,var(ydata));
xaxis=(minpoints:maxpoints)*(xdata(2)-xdata(1));

n=n+2;
if ~isempty(find(plts==n,3))
    figure(fign);
    subplot(2,2,3);
    plot(xaxis,Qp,'k');
    set(gca,'fontweight','b','fontsize',16);
    xlabel('Period (min)','fontweight','b','fontsize',16);
    ylabel('Qp','fontweight','b','fontsize',16);
    grid on;
end

% n=n+1;
% if ~isempty(find(plts==n,1))
%     figure(n);
%     plot(xdata(minpoints:maxpoints),Qp/maxpoints);
%     xlabel('Period');
%     ylabel('Robustness Qp (Fraction of the maximum possible Qp)');
% end

% ----------

% Plots the maximum Qp with increasing number of repetitions
% n=n+1;
% if ~isempty(find(plts==n,1))
%     figure(n);
%     plot(maxim);
%     xlabel('Number of repetitions');
%     ylabel('Maximum Qp');
% end

% ----------

% Find the extrema of the Qp function
[ymax,imax,xxx,yyy]=extrema(Qp);

if ~isempty(ymax)
    % Only consider Qp extrema values higher than 0.1 the maximum Qp value
    cand=find(ymax>(max(Qp)*0.1));
    ymax=ymax(cand);
    imax=imax(cand);
    
    % Remove the last point of xdata if it shows up
    cand=find(imax==length(xdata));
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
    
    % Filter results by the number of oscillation Qp values we are interested in
    if (length(ymax)>noscip)
        ymax=ymax(1:noscip);
        imax=imax(1:noscip);
    end
end

if ~isempty(ymax)
    if ~isempty(find(plts==3,1))
        figure(fign);
        subplot(2,2,3);
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
    
    % Calculate the significance
    sigpos=[];
    variance=var(ydata);
    for l=1:length(imax)
        minpoints=imax(l);
        maxpoints=minpoints;
        Qp_tt=zeros(1,maxperm);
        for m=1:maxperm
            Qp_tt(m)=calc_Qp(xdata,ydata_aux(m,:),minpoints,maxpoints,variance);
        end
        aux=length(find(Qp_tt>ymax(l)));
        pvalue=aux/maxperm;
        index=3*(l-1)+1;
        if (pvalue < 0.05)   % 5% Significance level
            sigpos=[sigpos l];
            res(index)=xdata(imax(l)+1)-xdata(1);
        else
            res(index)=(xdata(imax(l)+1)-xdata(1))*-1;
        end
        res(index+1)=ymax(l); % stores the power
        res(index+2)=pvalue; % stores the pvalue
    end
    
    if ~isempty(find(plts==3,1))
        sprintf('Significant Periods according to random data: ');
        for m=1:length(sigpos)
            figure(fign);
            subplot(2,2,3);
            hold on;
            plot(xaxis(imax(sigpos(m))),ymax(sigpos(m)),'k*');
        end
    end
    
end
       
% ----------

% This function calculates the Qp statistic value using the chi-square 
% periodogram for time series data
    function Qp = calc_Qp(x,y,minpoints,maxpoints,variance)
        len_Qp=maxpoints-minpoints+1;
        Qp=zeros(1,len_Qp);
        
        npoints=minpoints:maxpoints;
        for k=1:length(npoints)
            nblocks=floor(length(x)/npoints(k));
            table=zeros(nblocks,npoints(k));
            for i=1:nblocks
                for j=1:npoints(k)
                    pos=j+npoints(k)*(i-1);
                    table(i,j)=y(pos);
                end
            end
            if (len_Qp==1)
                Qp=npoints(k)*nblocks*var(mean(table,1))/variance;
            else
                Qp(k)=npoints(k)*nblocks*var(mean(table,1))/variance;
            end
        end
    end

% ----------

end