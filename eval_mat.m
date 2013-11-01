function eval_mat()
% --------------

clc;
close all;
clear all;

plts=[2];

per=1:0.5:20;
fs=0.2:0.2:3;
tt=1:50;

method='D';
mat_tt=zeros(length(per),length(fs),length(tt));

for i=1:(max(tt)/25)
    aux1=25*(i-1)+1;
    aux2=25*i;
    f=sprintf('simulations_03/mat_dft_tt_%d_%d_noise_0.mat',aux1,aux2);
    s=load(f);
    mat=s.mat;
    mat_tt(:,:,aux1:aux2)=mat;
end

[x y]=find(mat_tt==Inf);
sprintf('Infinits: %d',length(x))

% Some preprocessing
for i=1:length(per)
    mat_tt(i,:,:)=mat_tt(i,:,:)./per(i);
end

if ~isempty(find(plts==1))
    figure(1);
    mat_tt_aux=mat_tt;
    aux_tt=10:10:50;
    for i=1:length(aux_tt)
        subplot(1,5,i);
        surf(fs,per,mat_tt_aux(:,:,aux_tt(i)));
        axis([min(fs) max(fs) min(per) max(per) -0.00001 0.1 -0.00001 0.1]);
        %surf(0.2:0.2:2,1:0.5:10,mat_tt_aux(1:19,1:10,aux_tt(i)));
        %axis([0.2 2 1 10 0 aux3 0 aux3]);
        %axis([min(fs) max(fs) 1 100 0 aux3 0 aux3]);
        hold on;
        %shading interp;
        view(0,90);
        grid on;
        xlabel('Sampling frequency');
        ylabel('Period');
        s=sprintf('Time: %d',aux_tt(i));
        title(s);
        colorbar;
        myWait(0.01);
    end
end

% Plot average error percentage over sampling frequency
if ~isempty(find(plts==2))
    figure(2);
    mat_tt_aux=mat_tt;
    col=['r','b','g','c','k'];
    res1=zeros(1,length(fs));
    res2=zeros(1,length(fs));
    count=1;
    for j=10:10:50
        for k=1:length(fs)
            aux=mat_tt_aux(:,k,j);
            res1(k)=mean(aux(find(aux>=0)));
            res2(k)=length(find(aux<0));
        end
        subplot(2,2,1);
        hold on;
        y1=plot(fs,res1,col(count),'LineWidth',2,'MarkerSize',3);
        subplot(2,2,2);
        hold on;
        y2=plot(fs,res2,col(count),'LineWidth',2,'MarkerSize',3);
        count=count+1;
    end
    subplot(2,2,1);
    set(gca,'FontSize',14);
    xlabel('Sampling Frequency','fontsize',14,'fontweight','b','color','k');
    ylabel('Average Percentage Error','fontsize',14,'fontweight','b','color','k');
    axis([min(fs) max(fs) 0 0.15]);
    legend('10','20','30','40','50');
    if method=='E'
        title('Enright - APE');
    else
        title('DFT - APE');
    end
    grid on;
    subplot(2,2,2);
    set(gca,'FontSize',14);
    xlabel('Sampling Frequency','fontsize',14,'fontweight','b','color','k');
    ylabel('Number of non-significant points','fontsize',14,'fontweight','b','color','k');
    axis([min(fs) max(fs) 0 50]);
    legend('10','20','30','40','50');
    if method=='E'
        title('Enright - DSP');
    else
        title('DFT - DSP');
    end
    grid on;
end

% Plot average error percentage over time
if ~isempty(find(plts==2))
    figure(2);
    mat_tt_aux=mat_tt;
    col=['r','b','g','c','k'];
    res1=zeros(1,length(tt));
    res2=zeros(1,length(tt));
    count=1;
    for j=1:5
        for k=1:length(tt)
            aux=mat_tt_aux(:,j,k);
            res1(k)=mean(aux(find(aux>=0)));
            res2(k)=length(find(aux<0));
        end
        subplot(2,2,3);
        hold on;
        y1=plot(tt,res1,col(count),'LineWidth',2,'MarkerSize',3);
        subplot(2,2,4);
        hold on;
        y2=plot(tt,res2,col(count),'LineWidth',2,'MarkerSize',3);
        count=count+1;
    end
    subplot(2,2,3);
    set(gca,'FontSize',14);
    xlabel('Time','fontsize',14,'fontweight','b','color','k');
    ylabel('Average Percentage Error','fontsize',14,'fontweight','b','color','k');
    axis([min(tt) max(tt) 0 0.15]);
    legend('1','2','3','4','5');
    if method=='E'
        title('Enright - APE');
    else
        title('DFT - APE');
    end
    grid on;
    subplot(2,2,4);
    set(gca,'FontSize',14);
    xlabel('Time','fontsize',14,'fontweight','b','color','k');
    ylabel('Number of non-significant points','fontsize',14,'fontweight','b','color','k');
    axis([min(tt) max(tt) 0 50]);
    legend('1','2','3','4','5');
    if method=='E'
        title('Enright - DSP');
    else
        title('DFT - DSP');
    end
    grid on;
end

end