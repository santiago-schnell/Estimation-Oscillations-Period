function eval_mat3()
% --------------

clc;
close all;
clear all;

per1=1:10;
per2=1:10;
fs=0.2:0.2:3;
tt=5:5:100;

mat=[];
str='';

%correct();
build();

% Auxiliary code for making corrections on the data      
    function correct()
        for i=1:3
            switch i
                case 1
                    str='auto';
                case 2
                    str='enright';
                case 3
                    str='dft';
            end
            
            r1=16;
            r2=20;
            f1=sprintf('simulations_06/mat_%s_2pers_per1_%d_%d_noise_0.2.mat',str,r1,r2);
            s=load(f1);
            mat=s.mat;
            count=2;
            for j=r1:(r2-1)
                for k=count:j
                    pos=find((mat(:,1)==j) & (mat(:,2)==k));
                    mat(pos,:)=[];
                end
                count=count+1;
            end
            size(mat)
            str=strcat('simulations_06_corr/mat_',str,'_2pers_per1_',num2str(r1),'_',num2str(r2),'_noise_0.2.mat');
            save(str,'mat');      
        end
    end

% This function builds all the plots
    function build()
        for i=1:3
            switch i
                case 1
                    str='auto';
                case 2
                    str='enright';
                case 3
                    str='dft';
            end
            
            mat=[];
            for j=1:4
                f1=sprintf('simulations_05_corr/mat_%s_2pers_per1_%d_%d_noise_0.mat',str,5*j-4,5*j);
                s=load(f1);
                mat=[mat;s.mat];
            end
            % Add missing data
            f1=sprintf('simulations_05_corr/mat_%s_2pers_per1_%d_%d_noise_0.mat',str,5,6);
            s=load(f1);
            mat=[mat;s.mat];
            
            mat(find(mat==0))=Inf;
            
            % Plot results obtained from simulations
            %plot_x_y_z(i,mat);
            plot_hist(i,4);
            %plot_x_frequency(i,6);
            %plot_x_time(i,8);
        end
    end

% Plot x_y_z
    function plot_x_y_z(i)
        Z1=ones(45,length(fs)).*-1;
        Z2=ones(45,length(fs)).*-1;
        count=0;
        for j=2:4:length(tt)
            count=count+1;
            count2=0;
            for k=1:(length(per1)-1)
                for l=(k+1):length(per2)
                    count2=count2+1;
                    for m=1:length(fs)
                        pos=find((mat(:,1)==per1(k))&(mat(:,2)==per2(l))&(mat(:,3)==fs(m))&(mat(:,4)==tt(j)));
                        if not((per1(k)==10)|(per1(k)==15)) % lost values (obtain results for these)
                            res=sort(mat(pos,5:9));
                            Z1(count2,m)=abs(res(1)-mat(pos,1))./mat(pos,1);
                            Z2(count2,m)=abs(res(2)-mat(pos,2))./mat(pos,2);
                        end
                    end
                end
            end
            figure(i);
            for k=1:2
                switch k
                    case 1 
                        count_aux=count;Z=Z1;
                    case 2 
                        count_aux=count+5;Z=Z2;
                end
                subplot(2,5,count_aux);
                surf(fs,1:45,Z);
                axis([min(fs) max(fs) 1 45 0 Inf 0 0.2]);
                set(gca,'fontweight','b','fontsize',11);
                view(0,90);
                %shading interp;
                grid on;
                xlabel('Sampling frequency','fontweight','b','fontsize',11);
                ylabel('Period','fontweight','b','fontsize',11);
                s=sprintf('Time Range: %d',tt(j));
                title(s,'fontweight','b','fontsize',11);
                colorbar;
                set(gca,'fontweight','b','fontsize',11);
            end
        end
    end

% Plot histogram of positions 
    function plot_hist(i,plt)
        res1=[];
        res2=[];
        res3=[];
        res4=[];
        for j=1:length(tt)
            aux1=mat(find(mat(:,4)==tt(j)),:);
            for k=1:length(fs)
                aux2=aux1(find(aux1(:,3)==fs(k)),:);
                for l=1:(length(per1)-1)
                    for m=(l+1):length(per2)
                        if ~((per1(l)==10) | (per1(l)==15)) % lost values (obtain results for these)
                            aux3=aux2(find((aux2(:,1)==per1(l)) & (aux2(:,2)==per2(m))),:);
                            [res ind]=sort(aux3(:,5:9));
                            if ~(res(1)==Inf)
                                aux4=floor(tt(j)/per1(l)); % number of times the first real period repeats
                                aux6=abs((1:aux4).*per1(l)-(ones(1,aux4).*res(1))); % take differences to ff and harmonics 
                                [aux8 ind1]=min(aux6,[],2); % checks where the minimum occurs
                            else
                               ind(1)=0;
                               ind1=0;
                            end
                            if ~(res(2)==Inf)
                                aux5=floor(tt(j)/per2(m)); % number of times the second real period repeats
                                aux7=abs((1:aux5).*per2(m)-(ones(1,aux5).*res(2))); % take differences to ff and harmonics 
                                [aux8 ind2]=min(aux7,[],2); % checks where the minimum occurs
                            else
                               ind(2)=0;
                               ind2=0;
                            end
                            res1=[res1 ind(1)];                           
                            res2=[res2 ind(2)];
                            res3=[res3 ind1];
                            res4=[res4 ind2];
                        end
                    end
                end
            end
        end
        aux={res1 res3 res2 res4};
        for j=1:2
            figure(plt+j-1);
            subplot(3,2,2*i-1);
            N=hist(aux{2*j-1},max(aux{2*j-1})+1);
            bar(N./length(aux{2*j-1}));
            axis([1 6 0 1]);
            set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1';'2';'3';'4';'5']);
            xlabel('Position','fontweight','b','fontsize',13);
            ylabel('Normalized Count','fontweight','b','fontsize',13);
            title('Minimum position');
            subplot(3,2,2*i);
            N=hist(aux{2*j},max(aux{2*j})+1);
            bar(N./length(aux{2*j}));
            axis([1 6 0 1]);
            set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1';'2';'3';'4';'5']);
            xlabel('Harmonics','fontweight','b','fontsize',13);
            ylabel('Normalized Count','fontweight','b','fontsize',13);
            title('Minimum Harmonic position');
        end
    end

% Plot average error percentage over sampling frequency and number of
% non significant points
    function plot_x_frequency(i,plt)
        figure(plt);
        col=['r','b','g','c','k'];
        res1=zeros(1,length(fs));
        res2=zeros(1,length(fs));
        res3=zeros(1,length(fs));
        res4=zeros(1,length(fs));
        count=1;
        for j=4:4:length(tt)
            aux1=mat(find(mat(:,4)==tt(j)),:);
            for k=1:length(fs)
                aux2=aux1(find(aux1(:,3)==fs(k)),:);
                res=sort(aux2(:,5:9),2);
                aux3=abs(aux2(:,1)-res(:,1))./(aux2(:,1));
                aux4=abs(aux2(:,2)-res(:,2))./(aux2(:,2));
                res1(k)=mean(aux3(find(aux3~=Inf)));
                res2(k)=mean(aux4(find(aux4~=Inf)));
                res3(k)=length(find(aux3==Inf))./length(aux3);
                res4(k)=length(find(aux4==Inf))./length(aux4);
            end
            aux={res1 res3 res2 res4};
            for j=1:2
                figure(plt+j-1);
                subplot(3,2,2*i-1);
                hold on;
                y1=plot(fs,aux{2*j-1},col(count),'LineWidth',2,'MarkerSize',3);
                subplot(3,2,2*i);
                hold on;
                y2=plot(fs,aux{2*j},col(count),'LineWidth',2,'MarkerSize',3);
            end
            count=count+1;
        end
        for j=1:2
            figure(plt+j-1);
            subplot(3,2,2*i-1);
            xlabel('Sampling Frequency','fontweight','b','fontsize',13);
            ylabel('Average Error','fontweight','b','fontsize',13);
            axis([min(fs) max(fs) 0 5]);
            set(gca,'fontweight','b','fontsize',13);
            legend('20','40','60','80','100');
            title(str,'fontweight','b','fontsize',13);
            grid on;
            subplot(3,2,2*i);
            set(gca,'fontweight','b','fontsize',13);
            xlabel('Sampling Frequency','fontweight','b','fontsize',13);
            ylabel('Fraction of NSP','fontweight','b','fontsize',13);
            axis([min(fs) max(fs) 0 1]);
            set(gca,'fontweight','b','fontsize',13);
            legend('20','40','60','80','100');
            title(str,'fontweight','b','fontsize',13);
            grid on;
        end
    end

% Plot average error percentage over time and number of non significant points
    function plot_x_time(i,plt)
        figure(plt);
        col=['b','g','r','c','m','y','k'];
        res1=zeros(1,length(tt));
        res2=zeros(1,length(tt));
        res3=zeros(1,length(tt));
        res4=zeros(1,length(tt));
        count=1;
        for j=2:2:15
            aux1=mat(find(mat(:,3)==fs(j)),:);
            for k=1:length(tt)
                aux2=aux1(find(aux1(:,4)==tt(k)),:);
                minim=sort(aux2(:,5:9),2);
                aux3=abs(aux2(:,1)-minim(:,1))./(aux2(:,1));
                aux4=abs(aux2(:,2)-minim(:,2))./(aux2(:,2));
                res1(k)=mean(aux3(find(aux3~=Inf)));
                res2(k)=mean(aux4(find(aux4~=Inf)));
                res3(k)=length(find(aux3==Inf))./length(aux3);
                res4(k)=length(find(aux4==Inf))./length(aux4);                
            end
            aux={res1 res3 res2 res4};
            for j=1:2
                figure(plt+j-1); 
                subplot(3,2,2*i-1);
                hold on;
                y1=plot(tt,aux{2*j-1},col(count),'LineWidth',2,'MarkerSize',2);
                subplot(3,2,2*i);
                hold on;
                y2=plot(tt,aux{2*j},col(count),'LineWidth',2,'MarkerSize',2);
            end
            count=count+1;
        end
        for j=1:2
            figure(plt+j-1);
            subplot(3,2,2*i-1);
            xlabel('Time','fontweight','b','fontsize',13);
            ylabel('Average Error','fontweight','b','fontsize',13);
            axis([min(tt) max(tt) 0 5]);
            set(gca,'fontweight','b','FontSize',13);
            legend(num2str(fs(2:2:15)'));
            title(str,'fontweight','b','fontsize',13);
            grid on;
            subplot(3,2,2*i);
            xlabel('Time','fontweight','b','fontsize',13);
            ylabel('Fraction of NSP','fontweight','b','fontsize',13);
            axis([min(tt) max(tt) 0 1]);
            set(gca,'fontweight','b','fontsize',13);
            legend(num2str(fs(2:2:15)'));
            title(str,'fontweight','b','fontsize',13);
            grid on;
        end
    end

end

