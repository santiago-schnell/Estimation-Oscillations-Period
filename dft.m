function res=dft(xdata,ydata,fs,max_per,noscip,maxperm,func,plts,fign)
% ----------
% func is 'perms' or 'normal'
% maxperm is the number of permutations for the calculation of the null distribution

n=0;
res=zeros(1,noscip*3); % Final Results

% ----------

n=n+1;
if ~isempty(find(plts==n,1))
    figure(fign);
    plot(xdata,ydata,'k');
    xlabel('Time');
    ylabel('Value');
end

[f,mx]=calc_fft(ydata,fs,max_per);

n=n+3;
if ~isempty(find(plts==n,1))
    figure(fign);
    subplot(2,2,4);
    plot(f,mx,'k');
    set(gca,'fontweight','b','fontsize',16);
    xlabel('Frequency (cycles per min)','fontweight','b','fontsize',16);
    ylabel('Power','fontweight','b','fontsize',16);
    grid on;
end

% Find the extrema and sort by the highest power (descencent)
[ymax,imax,xxx,yyy]=extrema(mx);

if ~isempty(ymax)
    % Only consider extrema Power values higher than 0.1 the maximum Qp value
    cand=find(ymax>(max(mx)*0.1));
    ymax=ymax(cand);
    imax=imax(cand);
    
    % Filter results by the number of oscillation powers we are interested in
    if (length(ymax)>noscip)
        ymax=ymax(1:noscip);
        imax=imax(1:noscip);
    end
    
    % Remove zero frequencies if they exist
    freqs=find(f(imax)==0);
    imax(freqs)=[];
    ymax(freqs)=[];
 
%     % Sort results by the highest frequencies (lower periods)
%     [sr1 sr2]=sort(f(imax),'descend');
%     ymax=ymax(sr2);
%     imax=imax(sr2);
end

if ~isempty(ymax)
     
    if ~isempty(find(plts==4,1))
        figure(fign);
        subplot(2,2,4);
        hold on;
        plot(f(imax),ymax,'ko');
    end
     
    % ---------------
    
    count=0;
    for m=1:maxperm
        if strcmp(func,'perms')
            perm=randperm(length(xdata));
            ydata_aux=ydata(perm);
        elseif strcmp(func,'normal')
            ydata_aux=mean(ydata)+std(ydata).*randn(1,length(xdata));
        end
        [f_perm mx_perm]=calc_fft(ydata_aux,fs,max_per);
        [ymax_perm,imax_perm,xxx,yyy]=extrema(mx_perm);
        if ~isempty(ymax_perm)
            if ~(f_perm(imax_perm(1))==0)
                count=count+1;
                pow(count)=ymax_perm(1); % I am using the biggest power in the spectra
            end
        end
    end
    
    % --------------- 
        
    sigpos=[];
    for n=1:length(ymax)
        aux=length(find(pow>ymax(n)));
        pvalue=aux/count;
        index=3*(n-1)+1;
        if (pvalue < 0.05)   % 5% Significance level
            sigpos=[sigpos n];
            res(index)=1/f(imax(n));
        else
            res(index)=(1/f(imax(n)))*-1;
        end
        res(index+1)=ymax(n); % stores the power
        res(index+2)=pvalue; % stores the pvalue
    end
    
    if ~isempty(find(plts==4,1))
        sprintf('Significant Periods according to random data: ');
        if ~isempty(sigpos)
            figure(fign);
            subplot(2,2,4);
            hold on;
            plot(f(imax(sigpos)),ymax(sigpos),'k*');
        end
    end
end

% ----------

% This function calculates the frequency and the period of potential
% oscillations using the discrete fourier transform
    function [f,MX] = calc_fft(x,Fs,max_per)
        opt=1;
        
        Fn=Fs/2;                                % Nyquist frequency
        NFFT=2.^(ceil(log(length(x))/log(2)));  % Next highest power of 2
                                                % greater than length(x).
        
        if (opt==1)
            [MX,f]=periodogram(x,[],NFFT,Fs);   % The default window is used
            ind_lim=find(f>1/max_per);              % This imposes a lower limit
            f=f(ind_lim);                           % limit on the frequency
            MX=MX(ind_lim);     
        else
            FFTX=fft(x,NFFT);                       % Take FFT, padding with zeros.
                                                % length(FFTX)==NFFT
            NumUniquePts = ceil((NFFT+1)/2);
            FFTX=FFTX(1:NumUniquePts);              % FFT is symmetric, throw away
                                                % second half
            MX=abs(FFTX);                           % Take magnitude of X
            MX=MX*2;                                % Multiply by 2 to take into
                                                    % account the fact that we
                                                    % threw out second half of
                                                    % FFTX above
            MX(1)=MX(1)/2;                          % Account for endpoint
                                                    % uniqueness
            MX(length(MX))=MX(length(MX))/2;        % We know NFFT is even
            MX=MX/length(x);                        % Scale the FFT so that it is
                                                    % not a function of the length
                                                    % of x.
        
            f=(0:NumUniquePts-1)*2*Fn/NFFT;
        
            ind_lim=find(f>1/max_per);              % This imposes a lower limit
            f=f(ind_lim);                           % limit on the frequency
            MX=MX(ind_lim);     
        end
    end

end

