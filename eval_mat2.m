function eval_mat2()
% --------------

clc;
close all;
clear all;

per=1:20;
fs=0.2:0.2:3;
tt=5:5:100;
noise=0;

for i=1:1
    switch i
        case 1
            str='Auto';
        case 2
            str='Enright';
        case 3
            str='DFT';
    end
    
    mat=[];
    for j=1:2
        f=sprintf('simulations/mat_%s_per_%d_%d_noise_%d.mat',noise,str,10*j-9,10*j,noise);
        s=load(f);
        mat=[mat;s.mat];
    end
   
    mat(find(mat==0))=Inf;
        
    % Plot results obtained from simulations
    plot_np(i);
    %plot_nppp(i);
    %plot_Z(i,str);
    %plot_x_y_z_1(i,str);
    %plot_hist(i,1);
    %plot_x_frequency(i,1);
    %plot_x_time(i,1);
end

% Plot average error percentage and number of
% non significant points over number of periods 
    function plot_np(i)
        figure(i);
        ngroups=10;
        col=['b','r','g','k'];
        res1=zeros(1,ngroups);
        res2=zeros(1,ngroups);
        for j=1:4
            aux1=mat2(find((mat2(:,1)>=(j*5-4))&(mat2(:,1)<=(j*5))),:);
            ind=ceil(ngroups*tiedrank(aux1(:,2))/length(aux1(:,2)));
            middlep=zeros(1,ngroups);
            for k=1:ngroups
                aux2=aux1(find(ind==k),:);
                middlep(k)=(min(aux2(:,2))+max(aux2(:,2)))/2;
                aux3=min(aux2(:,4:8),[],2);
                aux4=abs(aux2(:,1)-aux3)./(aux2(:,1));
                res1(k)=mean(aux4(find(aux4~=Inf)));
                res2(k)=length(find(aux4==Inf))./length(aux4);
            end
            subplot(2,4,j);
            hold on;
            y1=plot(middlep,res1,col(j),'LineWidth',4,'MarkerSize',4);
            set(gca,'fontweight','b','fontsize',13);   
            xlabel('Number of periods','fontweight','b','fontsize',13);
            ylabel('Average Difference Error','fontweight','b','fontsize',13);
            axis([min(middlep) max(middlep) 0 1]);
            set(gca,'fontweight','b','fontsize',13);
            title(str,'fontweight','b','fontsize',13);
            grid on;
            subplot(2,4,j+4);
            hold on;
            y2=plot(middlep,res2,col(j),'LineWidth',4,'MarkerSize',4);
            set(gca,'fontweight','b','fontsize',13);   
            xlabel('Number of periods','fontweight','b','fontsize',13);
            ylabel('Fraction of NSP','fontweight','b','fontsize',13);
            axis([min(middlep) max(middlep) 0 0.2]);
            set(gca,'fontweight','b','fontsize',13);
            title(str,'fontweight','b','fontsize',13);
            grid on;
        end
    end

% Plot average error percentage and number of
% non significant points over number of points per period
    function plot_nppp(i)
        figure(i);
        ngroups=10;
        col=['b','r','g','k'];
        res1=zeros(1,ngroups);
        res2=zeros(1,ngroups);
        for j=1:4
            aux1=mat2(find((mat2(:,1)>=(j*5-4))&(mat2(:,1)<=(j*5))),:);
            ind=ceil(ngroups*tiedrank(aux1(:,3))/length(aux1(:,3)));
            middlep=zeros(1,ngroups);
            for k=1:ngroups
                aux2=aux1(find(ind==k),:);
                middlep(k)=(min(aux2(:,3))+max(aux2(:,3)))/2;
                aux3=min(aux2(:,4:8),[],2);
                aux4=abs(aux2(:,1)-aux3)./(aux2(:,1));
                res1(k)=mean(aux4(find(aux4~=Inf)));
                res2(k)=length(find(aux4==Inf))./length(aux4);
            end
            subplot(2,4,j);
            hold on;
            y1=plot(middlep,res1,col(j),'LineWidth',4,'MarkerSize',4);
            set(gca,'fontweight','b','fontsize',13);   
            xlabel('Number of points per period','fontweight','b','fontsize',13);
            ylabel('Average Difference Error','fontweight','b','fontsize',13);
            axis([min(middlep) max(middlep) 0 1]);
            set(gca,'fontweight','b','fontsize',13);
            title(str,'fontweight','b','fontsize',13);
            grid on;
            subplot(2,4,j+4);
            hold on;
            y2=plot(middlep,res2,col(j),'LineWidth',4,'MarkerSize',4);
            set(gca,'fontweight','b','fontsize',13);
            xlabel('Number of points per period','fontweight','b','fontsize',13);
            ylabel('Fraction of NSP','fontweight','b','fontsize',13);
            axis([min(middlep) max(middlep) 0 0.2]);
            set(gca,'fontweight','b','fontsize',13);
            title(str,'fontweight','b','fontsize',13);
            grid on;
        end
    end

% Plot Z_1
    function plot_Z(i,str)
        uscmat=unique(mat2(:,2))'; % unique second column of mat
        utcmat=unique(mat2(:,3))'; % unique third column of mat        
        Z=ones(length(per),length(fs)).*-1;
        count=0;
        for j=2:2:(length(tt)/2)
            count=count+1;
            for k=1:length(per)
                for l=1:length(fs)
                    pos=find((mat(:,1)==per(k)) & (mat(:,2)==fs(l)) & (mat(:,3)==tt(j)));
                    minim=min(mat(pos,4:8),[],2);
                    Z(k,l)=abs(minim-mat(pos,1))./mat(pos,1);
                end
            end
            figure(i);
            subplot(2,5,count);
            surf(fs,per,Z);
            axis([min(fs) max(fs) min(per) max(per) 0 Inf 0 0.2]);
            set(gca,'fontweight','b','FontSize',11)
            view(0,90);
            %shading interp;
            grid on;
            xlabel('SF','fontweight','b','fontsize',11);
            ylabel('Period','fontweight','b','fontsize',11);
            s=sprintf('TR-1:%d %s',tt(j),str);
            title(s,'fontweight','b','fontsize',11);
            colorbar;
            set(gca,'fontweight','b','FontSize',11)
        end
    end

% Plot x_y_z_1
    function plot_x_y_z_1(i,str)
        Z=ones(length(per),length(fs)).*-1;
        count=0;
        for j=2:2:(length(tt)/2)
            count=count+1;
            for k=1:length(per)
                for l=1:length(fs)
                    pos=find((mat(:,1)==per(k)) & (mat(:,2)==fs(l)) & (mat(:,3)==tt(j)));
                    minim=min(mat(pos,4:8),[],2);
                    Z(k,l)=abs(minim-mat(pos,1))./mat(pos,1);
                end
            end
            figure(i);
            subplot(2,5,count);
            surf(fs,per,Z);
            axis([min(fs) max(fs) min(per) max(per) 0 Inf 0 0.2]);
            set(gca,'fontweight','b','FontSize',11)
            view(0,90);
            %shading interp;
            grid on;
            xlabel('SF','fontweight','b','fontsize',11);
            ylabel('Period','fontweight','b','fontsize',11);
            s=sprintf('TR-1:%d %s',tt(j),str);
            title(s,'fontweight','b','fontsize',11);
            colorbar;
            set(gca,'fontweight','b','FontSize',11)
        end
    end

% Plot histogram of positions 
    function plot_hist(i,plt)
        maxh=3; % maximum number of harmonics
        res1=zeros(6,maxh);
        res2=zeros(6,maxh);
        for j=1:length(tt)
            aux1=mat(find(mat(:,3)==tt(j)),:);
            for k=1:length(fs)
                aux2=aux1(find(aux1(:,2)==fs(k)),:);
                for l=1:length(per)
                    aux3=aux2(find(aux2(:,1)==per(l)),:);
                    [aux4 ind1]=min(aux3(:,4:8),[],2);
                    if ~(aux4==Inf)
                        aux5=floor(tt(j)/per(l)); % number of times real period repeats
                        aux6=abs((1:aux5).*per(l)-(ones(1,aux5).*aux4)); % take differences to ff and harmonics 
                        [aux7 ind2]=min(aux6,[],2); % checks where minimum occurs
                        if ~(ind2>maxh)
                            res1(ind1+1,ind2)=(res1(ind1+1,ind2)*res2(ind1+1,ind2)+(aux7/(per(l)*ind2)))/(res2(ind1+1,ind2)+1);
                            res2(ind1+1,ind2)=res2(ind1+1,ind2)+1;
                        end
                    else
                        res2(1,1)=res2(1,1)+1;
                    end
                end
            end
        end
        figure(plt);
        subplot(3,2,2*i);
        bar(1:6,res1);
        ylim([0 1]);
        set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1';'2';'3';'4';'5']);
        xlabel('Position in the results vector','fontweight','b','fontsize',13);
        ylabel('Average Difference Error','fontweight','b','fontsize',13);
        legend('Fundamental frequency','First harmonic','Second harmonic');
        grid;
        title('Harmonics according to position error (20% noise)','fontweight','b','fontsize',13);
        figure(plt+1);
        subplot(3,2,2*i);
        res2=res2./sum(sum(res2)); % normalize result matrix
        bar(1:6,res2);
        ylim([0 1]);
        set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1';'2';'3';'4';'5']);
        xlabel('Position in the results vector','fontweight','b','fontsize',13);
        ylabel('Normalized Count','fontweight','b','fontsize',13);
        legend('Fundamental frequency','First harmonic','Second harmonic');
        grid;
        title('Harmonics according to position count (20% noise)','fontweight','b','fontsize',13);
    end

% Plot average error percentage over sampling frequency and number of
% non significant points
    function plot_x_frequency(i,plt)
        figure(plt);
        col=['r','b','g','c','k'];
        res1=zeros(1,length(fs));
        res2=zeros(1,length(fs));
        count=1;
        for j=4:4:length(tt)
            aux1=mat(find(mat(:,3)==tt(j)),:);
            for k=1:length(fs)
                aux2=aux1(find(aux1(:,2)==fs(k)),:);
                aux3=min(aux2(:,4:8),[],2);
                aux4=abs(aux2(:,1)-aux3)./(aux2(:,1));
                res1(k)=mean(aux4(find(aux4~=Inf)));
                res2(k)=length(find(aux4==Inf))./length(aux4);
            end
            subplot(3,2,2*i-1);
            hold on;
            y1=plot(fs,res1,col(count),'LineWidth',2,'MarkerSize',3);
            subplot(3,2,2*i);
            hold on;
            y2=plot(fs,res2,col(count),'LineWidth',2,'MarkerSize',3);
            count=count+1;
        end
        subplot(3,2,2*i-1);
        set(gca,'FontSize',14);
        xlabel('Sampling Frequency','fontweight','b','fontsize',13);
        ylabel('Average Error','fontweight','b','fontsize',13);
        axis([min(fs) max(fs) 0 3]);
        set(gca,'fontweight','b','fontsize',13);
        legend('TR-1:20','TR-1:40','TR-1:60','TR-1:80','TR-1:100');
        title(str);
        grid on;
        subplot(3,2,2*i);
        set(gca,'FontSize',14);
        xlabel('Sampling Frequency','fontweight','b','fontsize',13);
        ylabel('Fraction of NSP','fontweight','b','fontsize',13);
        axis([min(fs) max(fs) 0 1]);
        set(gca,'fontweight','b','fontsize',13);
        legend('TR-1:20','TR-1:40','TR-1:60','TR-1:80','TR-1:100');
        title(str,'fontweight','b','fontsize',13);
        grid on;
    end

% Plot average error percentage over time and number of non significant points
    function plot_x_time(i,plt)
        figure(plt);
        col=['b','g','r','c','m','y','k'];
        res1=zeros(1,length(tt));
        res2=zeros(1,length(tt));
        count=0;
        for j=2:2:15
            count=count+1;
            aux1=mat(find(mat(:,2)==fs(j)),:);
            for k=1:length(tt)
                aux2=aux1(find(aux1(:,3)==tt(k)),:);
                aux3=min(aux2(:,4:8),[],2);
                aux4=abs(aux2(:,1)-aux3)./(aux2(:,1));
                res1(k)=mean(aux4(find(aux4~=Inf)));
                res2(k)=length(find(aux4==Inf))./length(aux4);
            end
            subplot(3,2,2*i-1);
            hold on;
            y1=plot(tt,res1,col(count),'LineWidth',2,'MarkerSize',2);
            subplot(3,2,2*i);
            hold on;
            y2=plot(tt,res2,col(count),'LineWidth',2,'MarkerSize',2);
        end
        subplot(3,2,2*i-1);
        set(gca,'FontSize',14);
        xlabel('Time Range','fontweight','b','fontsize',13);
        ylabel('Average Error','fontweight','b','fontsize',13);
        axis([min(tt) max(tt) 0 2]);
        set(gca,'fontweight','b','fontsize',13);
        legend(num2str(fs(2:2:15)'));
        title(str);
        grid on;
        subplot(3,2,2*i);
        set(gca,'fontweight','b','fontsize',13);
        xlabel('Time Range','fontweight','b','fontsize',13);
        ylabel('Fraction of NSP','fontweight','b','fontsize',13);
        axis([min(tt) max(tt) 0 1]);
        set(gca,'fontweight','b','fontsize',13);
        legend(num2str(fs(2:2:15)'));
        title(str);
        grid on;
    end

end

