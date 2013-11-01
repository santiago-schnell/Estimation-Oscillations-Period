function eval_mat4()
% --------------

clc;
clear all;
%close all;
warning('off');

per=1:20;
nperiods=2:20;
npperiod=1:20;
noise=60;
str={'Auto','Enright','DFT'};
marker={'d','s','o'};
linestyle={':','--','-'};
colors={'r','g','b'};

avg_res5=zeros(6,3);
std_res5=zeros(6,3);
res3=zeros(6,5);

for i=1:3
    mat=[];
    
    for j=1:2
        f=sprintf('simulations_ver3/mat_%s_per_%d_%d_noise_%d.mat',str{i},10*j-9,10*j,noise);
        s=load(f);
        mat=[mat;s.mat];
    end
    
    mat(find((mat(:,3)==1)),:)=[]; % removes mat positions with one point per period
    mat(find((mat(:,3)==2)),:)=[]; % removes mat positions with two points per period
    
    % Plot results obtained from simulations
    plot_results(2,i);
    %calc_period_variance(i);
    %plot_hist(i);
end

%plot_hist_aux;

% Plot average error as a function of the number of
% periods or the number of points per period
    function plot_results(xaxis,i)
        
        if (xaxis==1)
            init=1;
            final=19;
            array=nperiods;
            axis_label='Number of cycles (NC)';
            col=[3 2];
        else if (xaxis==2)
                init=3;
                final=20;
                array=npperiod;
                axis_label='Number of points per cycle (NPC)';
                col=[2 3];
            end
        end
        
        aux1=mat(find((mat(:,col(1))>=1)&(mat(:,col(1))<=20)),:);
        count3=1;
        for k=init:final
            aux2=aux1(find(aux1(:,col(2))==array(k)),:);
            aux3=abs(aux2(:,4:3:16)); % the period
            aux4=aux2(:,5:3:17); % the power
            aux5=aux2(:,6:3:18); % the pvalue
            
            % Obtains the average relative difference of periods
            aux3(find(aux3==0))=Inf;
            %aux6=abs((repmat(aux2(:,1),1,5))-aux3);
            %aux7=min(aux6,[],2);
            
            [x y]=size(aux3);
            for l=1:x
                xxx=min(aux5(l,aux3(l,:)~=Inf));
                xxx(find(xxx>=0.01))=[];
                if ~isempty(xxx)
                    yyy=find(aux5(l,:)==xxx);
                    aux7(l)=aux3(l,yyy(1));
                else
                    aux7(l)=Inf;
                end
            end
            aux7=abs(aux7(:,1)-aux2(:,1));
            
            aux8=aux7./(aux2(:,1));
            aux9=aux8(find(aux8~=Inf));
            
            % Stores results
            res_mean_avgper(count3)=mean(aux9);
            res_std_avgper(count3)=std(aux9)/sqrt(length(aux9));
            count3=count3+1;
        end
        
        % Produces the plots
        figure(1);
        subplot(3,2,6);
        hold on;
        h=errorbar(array(init:final),res_mean_avgper,res_std_avgper,'k');
        set(h(1),'Marker',marker{i},'LineStyle',linestyle{i},'DisplayName',str{i},'Color',colors{i});
        set(gca,'fontweight','b','fontsize',16);
        xlabel(axis_label,'fontweight','b','fontsize',16);
        ylabel('Average RD','fontweight','b','fontsize',16);
        axis([1 20 0 0.08]);
        set(gca,'fontweight','b','fontsize',16);
        aux=['Noise - ',num2str(noise),'%'];
        title(aux,'fontweight','b','fontsize',16);
        legend(str);
        grid on;
    end

% Calculates period variance
    function calc_period_variance(i)
        count1=1;
        for j=1:length(nperiods)
            for k=3:length(npperiod)
                aux1=mat(find((mat(:,2)==nperiods(j))&(mat(:,3)==npperiod(k))),:);
                count2=1;
                for l=1:length(per)
                    aux2=aux1(find(aux1(:,1)==per(l)),:);
                    aux3=abs(aux2(:,4:3:16)); % the period
                    aux4=aux2(:,5:3:17); % the power
                    aux5=aux2(:,6:3:18); % the pvalue
                    
                    % Obtains the average relative difference of periods
                    aux3(find(aux3==0))=Inf;
                    aux6=abs((repmat(aux2(:,1),1,5))-aux3);
                    [aux7 ind1]=min(aux6,[],2);
                    aux8=aux7./(aux2(:,1));
                    aux9=aux8(find(aux8~=Inf));
                    
                    % Only consider relevant points
                    if ~isempty(aux9)
                        res1(count2)=aux9;
                        count2=count2+1;
                    end
                end
                final_res1(count1)=var(res1);
                count1=count1+1;
            end
        end
        
        % Write final results
        mean(final_res1)
    end

% Plot histogram of positions
    function plot_hist(i)
        res1=zeros(6,6);
        res2=zeros(6,5);
        res4=zeros(6,2);
        res5=zeros(6,1);
        count=1;
        for j=1:length(nperiods)
            for k=3:length(npperiod)
                aux1=mat(find((mat(:,2)==nperiods(j))&(mat(:,3)==npperiod(k))),:);
                for l=1:length(per)
                    aux2=aux1(find(aux1(:,1)==per(l)),:);
                    aux3=abs(aux2(:,4:3:16)); % the period
                    aux4=aux2(:,5:3:17); % the power
                    aux5=aux2(:,6:3:18); % the pvalue
                    
                    % Obtains the average relative difference of periods
                    aux3(find(aux3==0))=Inf;
                    if ~isempty(find(aux3~=Inf))
                        aux6=abs((repmat(aux2(:,1),1,5))-aux3);
                        [aux7 ind1]=min(aux6,[],2); % closest period
                        [aux8 ind2]=max(aux4(aux3~=Inf),[],2); % maximum power
                        [aux9 ind3]=min(aux5(aux3~=Inf),[],2); % minimum pvalue
                        aux10=aux7./(aux2(:,1));
                        aux11=(abs(aux3(ind2)-(aux2(:,1))))./(aux2(:,1));
                        aux12=(abs(aux3(ind3)-(aux2(:,1))))./(aux2(:,1));
                        
                        % Take the closest period with pvalue < 0.01
                        aux13=sort(aux6(find(aux5<0.01)),'ascend');
                        if ~isempty(aux13)
                            ind4=find(aux6==aux13(1));
                            aux14=(abs(aux3(ind4)-(aux2(:,1))))./(aux2(:,1));
                        else
                            aux14=Inf;
                        end
                        
                        % Take the maximum power with pvalue < 0.01
                        aux15=sort(aux4(find(aux5<0.01)),'descend');
                        if ~isempty(aux15)
                            ind5=find(aux4==aux15(1));
                            aux16=(abs(aux3(ind5)-(aux2(:,1))))./(aux2(:,1));
                        else
                            aux16=Inf;
                        end
                        
                        % Take the minimum pvalues with pvalue < 0.01
                        aux17=sort(aux5(find(aux5<0.01)),'ascend');
                        if ~isempty(aux17)
                            ind6=find(aux5==aux17(1));
                            aux18=(abs(aux3(ind6)-(aux2(:,1))))./(aux2(:,1));
                        else
                            aux18=Inf;
                        end
                    else
                        aux10=Inf;
                        aux11=Inf;
                        aux12=Inf;
                    end
                    
                    % Sort final results
                    aux19=[ind1 ind4 ind2 ind5 ind3 ind6];
                    aux20=[aux10 aux14 aux11 aux16 aux12 aux18];
                    
                    for opt=1:6
                        ind=aux19(opt);
                        aux=aux20(opt);
                        if aux~=Inf
                            res1(opt,ind+1)=res1(opt,ind+1)+1;
                            res2(opt,ind,count)=aux;
                            res4(opt,2)=res4(opt,2)+1;
                            res5(opt,i,count)=aux;
                        else
                            res1(opt,1)=res1(opt,1)+1;
                            res4(opt,1)=res4(opt,1)+1;
                        end
                    end
                    count=count+1;
                    
                    if (aux7~=Inf)
                        aux21=find(sort(aux4,'descend')==aux4(ind1)); %gets power position of the closest period
                        aux22=find(sort(aux5,'ascend')==aux5(ind1)); %gets p-value position of the closest period
                        res3(2*i-1,aux21(1))=res3(2*i-1,aux21(1))+1;
                        res3(2*i,aux22(1))=res3(2*i,aux22(1))+1;
                    end
                end
            end
        end
        
        for opt=1:6
            avg_res2(opt,:)=mean(res2(opt,:,:),3);
            avg_res5(opt,i)=mean(res5(opt,i,:),3);
            res1(opt,:)=res1(opt,:)./sum(res1(opt,:)); % normalize results
            res4(opt,:)=res4(opt,:)./sum(res4(opt,:)); % normalize results
            for l=1:5
                std_res2(opt,l)=std(res2(opt,l,:));
            end
            std_res5(opt,i)=std(res5(opt,i,:))/sqrt(length(res5(opt,i,:)));
        end
        
        res3(2*i-1,:)=res3(2*i-1,:)./sum(res3(2*i-1,:)); % normalize results
        res3(2*i,:)=res3(2*i,:)./sum(res3(2*i,:)); % normalize results
        
        %         figure(1);
        %         subplot(3,2,2*i-1);
        %         title(str(i),'fontweight','b','fontsize',13);
        %         hold on;
        %         barweb(res1',zeros(6,6));
        %         xlim([0 7]);
        %         ylim([0 1]);
        %         set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1';'2';'3';'4';'5']);
        %         xlabel('Position in the results vector','fontweight','b','fontsize',13);
        %         ylabel('Normalized Count','fontweight','b','fontsize',13);
        %         legend('Closest','Closest Elm NS','Max Power','Max Power Elm NS','Min P-Value','Min P-Value Elm NS');
        %         grid on;
        %         subplot(3,2,2*i);
        %         title(str(i),'fontweight','b','fontsize',13);
        %         hold on;
        %         barweb(-log(avg_res2)',-log(std_res2)');
        %         hold on;
        %         xlim([0 6]);
        %         ylim([0 15]);
        %         set(gca,'fontweight','b','fontsize',13,'XTickLabel',['1';'2';'3';'4';'5']);
        %         xlabel('Position in results (increasing period)','fontweight','b','fontsize',13);
        %         ylabel('-log(Relative difference)','fontweight','b','fontsize',13);
        %         %legend('Closest','Closest Elm NS','Maximum Power','Maximum Power Elm NS','Minimum P-Value','Minimum P-Value Elm NS');
        %         grid on;
        %
        %         figure(3);
        %         subplot(3,2,2*i-1);
        %         title(str(i),'fontweight','b','fontsize',13);
        %         hold on;
        %         barweb(res4',zeros(2,6));
        %         xlim([0 3]);
        %         ylim([0 1]);
        %         set(gca,'fontweight','b','fontsize',13,'XTickLabel',['0';'1']);
        %         ylabel('Normalized Count','fontweight','b','fontsize',13);
        %         legend('Closest','Closest (p<0.01)','Max Power','Max Power (p<0.01)','Min P-Value','Min P-Value (p<0.01)');
        %         grid on;
        
    end

    function plot_hist_aux()
        %         figure(1)
        %         subplot(2,2,1);
        %         title('Noise - 0%','fontweight','b','fontsize',16);
        %         hold on;
        %         barweb(res3',zeros(5,6));
        %         xlim([0 6]);
        %         ylim([0 1]);
        %         set(gca,'fontweight','b','fontsize',13,'XTickLabel',['1';'2';'3';'4';'5']);
        %         xlabel('Highest to Lowest Power / Lowest to Highest P-Value','fontweight','b','fontsize',13);
        %         ylabel('Normalized Count','fontweight','b','fontsize',13);
        %         legend('Power of closest period-Auto','P-value of closest period-Auto','Power of closest period-Enright','P-value of closest period-Enright','Power of closest period-DFT','P-value of closest period-DFT');
        %         grid on;
        
        subplot(3,1,3);
        title('Noise - 60%','fontweight','b','fontsize',16);
        hold on;
        avg_res5(1:2,:)=[];
        avg_res5(2:3,:)=[];
        barweb(-log(avg_res5)',zeros(3,2));
        hold on;
        ylim([0 10]);
        set(gca,'fontweight','b','fontsize',13,'XTickLabel',['1';'2';'3']);
        ylabel('-log(RD)','fontweight','b','fontsize',13);
        %legend('Closest','Closest (p<0.01)','Max Power','Max Power (p<0.01)','Min P-Value','Min P-Value (p<0.01)');
        legend('Max Peak','Min p-value (p<0.01)');
        grid on;
    end

end

