function []=test(ntest)

close all;
clc;
n=0;
% ------------------------------ Test 1 ------------------------------

n=n+1;
if ~isempty(find(ntest==1,1))
    npperiod=15; % Number of points per period
    per=[10]; % Period time
    nperiods=10; % Number of repetitions
    noise=0; % Noise level
    [xdata,ydata]=gen_input(ones(1,length(per)),per,nperiods,npperiod,noise);
    %[ymax,imax,xxx,yyy]=extrema(ydata);
    %ydata(imax)=ydata(imax)+0.*randn(1,length(imax));
    %ydata=ydata./(ones(1,length(ydata))*var(ydata));
    [res1,res2,res3]=apply(xdata,ydata,5,[1 2 3],1)
    subplot(2,1,2);
    plot(xdata,ydata,'b','LineWidth',2);
    set(gca,'box','off');
    set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
    axis off;
end

% ------------------------------ Test 2 ------------------------------

n=n+1;
if ~isempty(find(ntest==2,1))
    display('Test 2');
    data=load('rythms.txt')
    [x y]=size(data);
    data=data(2:x,:);
    [x y]=size(data);
    xdata=data(:,3);
    ydata=data(:,4)-(ones(x,1)*mean(data(:,4)));
    [res1 res2 res3]=apply(xdata,ydata,5,[1 2 3],1)
    
    %     for i=0:6
    %         i
    %         str=strcat('test_data/00000',num2str(i),'.datasets2.dat');
    %         data=load(str);
    %         [x y]=size(data);
    %         xdata=data(:,1);
    %         ydata=data(:,2)-(ones(x,1)*mean(data(:,2)));
    %         [res1 res2 res3]=apply(xdata,ydata,5,[1 2 3],1)
    %     end
end

% ------ Results ------

% Auto (Perms) - 000000 (21,23,27), 000001 (2), 000002 (23,26,28), 000003 (12), 000005 (8)
% Enright (Perms) - 000006 (24)
% FFT (Perms) -

% ------------------------------ Test 3 ------------------------------

n=n+1;
if ~isempty(find(ntest==3,1))
    display('Test 3');
    data=load('test_data/insulin.csv');
    [x,y]=size(data);
    for i=2:5
        i
        xdata=data(:,1);
        ydata=data(:,i)-(ones(x,1)*mean(data(:,i)));
        [res1,res2,res3]=apply(xdata,ydata,5,[1 2 3],1);
    end
end

% ------ Results ------

% Auto (Perms) - 2 (21,23,27), 3 (2), 4 (23,26,28)
% Enright (Perms) -
% FFT (Perms) -

% ------------------------------ Test 4 ------------------------------

n=n+1;
if ~isempty(find(ntest==4,1))
    display('Test 4');
    data=load('sunspot.dat');
    [x,y]=size(data);
    xdata=data(:,1);
    ydata=data(:,2)-(ones(x,1)*mean(data(:,2)));
    [res1,res2,res3]=apply(xdata,ydata,5,[1 2 3],1)
end

% ------ Results -----

% Auto (Perms) - (11,32,43,168,179)
% Enright (Perms) - (11,20,30)
% FFT (Perms) - (10.0044,11.0918,102.0444,170.0741)

% ------------------------------ Test 5 ------------------------------

n=n+1;
if ~isempty(find(ntest==5,1))
    display('Test 5');
    data=load('test_data/Protusion/SUM149-Z-rowdata-extended.csv');
    [x y]=size(data);
    xdata=data(1,2:y);
    for i=2:x
        str=strcat('----- Analyzing:  ',num2str(i),' with identifier: ',num2str(data(i,1)),' in file: something.csv', ' -----');
        display(str);
        ydata=data(i,2:y)-(ones(1,y-1)*mean(data(i,2:y)));
        [res1,res2,res3]=apply(xdata,ydata,10,[1 2 3],i)
    end
end

% ------ Results -----

% Enright (Perms) - 24 (4.2), 28.2 (5.1), 11 (5.5 and 6), 12 (6.6)
% FFT (Perms) - 2 (3.1333), 5 (4.7), 9 (4.7), 11 (4.7)

% ------------------------------ Test 6 ------------------------------

n=n+1;
if ~isempty(find(ntest==6,1))
    display('Test 6');
    data=load('test_data/Diabetes/data4.csv');
    [x y]=size(data);
    xdata=data(:,1);
    for i=2:y
        str=strcat('----- Analyzing identifier: ',num2str(i),' in file: data4.csv', ' -----');
        display(str);
        ydata=data(:,i)-(ones(x,1)*mean(data(:,i)));
        %subplot(ceil(y/4),ceil(y/4),i);
        %plot(xdata,ydata);
        apply(xdata,ydata,10,[1 2 3],i)
    end
end

% ------------------------------ Test 7 ------------------------------

n=n+1;
if ~isempty(find(ntest==7,1))
    display('Test 7');
    data=load('test_data/Diabetes/data1.csv');
    %pos1=find(data==0.3552);
    %pos2=find(data==28.4557);
    %pos3=find(data==61.3679);
    [x y]=size(data);
    pos_min=1;
    pos_max=x;
    xdata=data(:,1);
    %fid=fopen('test_data/Diabetes/SR03-F01-Fura2_results.txt','a');
    %figure(1);
    for i=4:4
        str=strcat('----- Analyzing column: ',num2str(i),' in file: data1.csv', ' -----');
        display(str);
        ydata=data(pos_min:pos_max,i)-(ones(length(pos_min:pos_max),1)*mean(data(pos_min:pos_max,i)));
        %subplot(4,4,i);
        %plot(xdata,ydata);
        grid on;
        [res1,res2,res3]=apply(xdata,ydata,5,[1 2 3],i)
        %fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n',res3);
    end
    %fclose(fid);
end

% ------------------------------ Test 8 ------------------------------

n=n+1;
if ~isempty(find(ntest==8,1))
    display('Test 8');
    data=load('test_data/Circadian/hat7PCRnumber.txt');
    [x y]=size(data);
    xdata=data(:,1);
    for i=5:5
        str=strcat('----- Analyzing identifier: ',num2str(i),' in file: hat7PCRnumber.csv', ' -----');
        display(str);
        ydata=data(:,i)-(ones(x,1)*mean(data(:,i)));
        [res1 res2 res3]=apply(xdata,ydata,3,[1 2 3],i)
    end
end

% ------------------------------ Test 9 ------------------------------

n=n+1;
if ~isempty(find(ntest==9,1))
    display('Test 9');
    data=load('~/Desktop/points.csv');
    [x y]=size(data);
    xdata=data(:,1);
    for i=2:14
        str=strcat('----- Analyzing identifier: ',num2str(i),' in file: points.csv', ' -----');
        display(str);
        ind_neg_values=(find(data(:,i)==-1));
        if (isempty(ind_neg_values))
            ind_neg_values=x;
        end
        ydata=data(1:ind_neg_values(1)-1,i)-(ones(ind_neg_values(1)-1,1)*mean(data(1:ind_neg_values(1)-1,i)));
        [res1,res2,res3]=apply(xdata(1:ind_neg_values(1)-1),ydata,3,[1 2 3],i)
    end
end

% ------------------------------ Test 10 ------------------------------

n=n+1;
if ~isempty(find(ntest==10,1))
    display('Test 10');
    data=load('~/Desktop/oscillations.csv');
    [x,y]=size(data);
    xdata=data(:,1);
    for i=2:y
        str=strcat('----- Analyzing identifier: ',num2str(i),' in file: oscillations.csv', ' -----');
        display(str);
        ind_neg_values=(find(data(:,i)==-1));
        if (isempty(ind_neg_values))
            ind_neg_values=x;
        end
        ydata=data(1:ind_neg_values(1)-1,i)-(ones(ind_neg_values(1)-1,1)*mean(data(1:ind_neg_values(1)-1,i)));
        [res1 res2 res3]=apply(xdata(1:ind_neg_values(1)-1),ydata,5,[1 2 3],i)
    end
end

% --------------------------- Apply Function --------------------------

    function [res1,res2,res3]=apply(xdata,ydata,nop,methods,fign)
        for i=1:length(methods)
            method=methods(i);
            if method==1
                display('Applying the Autocorrelation method...');
                res1=autocorrelation(xdata,ydata,nop,10000,'perms',[1 2],fign);
            else
                if method==2
                    display('Applying the Enright Periodogram method...');
                    per_range=max(xdata)-min(xdata);
                    res2=chi2(xdata,ydata,min(xdata)+floor(per_range/2),nop,10000,'perms',[3],fign);
                else
                    display('Applying the FFT method...');
                    per_range=max(xdata)-min(xdata);
                    Fs_aux=length(xdata)/per_range;
                    res3=dft(xdata,ydata,Fs_aux,floor(per_range/2),nop,10000,'perms',[4],fign);
                end
            end
        end
    end

end