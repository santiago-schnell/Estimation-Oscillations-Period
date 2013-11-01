function evaluate2(min_per1,max_per1,min_per2,max_per2,min_nperiods,max_nperiods,min_npperiod,max_npperiod,noise,mode)

per1=str2double(min_per1):str2double(max_per1);
per2=str2double(min_per2):str2double(max_per2);
nperiods=str2double(min_nperiods):str2double(max_nperiods);
npperiod=str2double(min_npperiod):str2double(max_npperiod);
nop=5;

tam=sum(1:(length(per1)-1))*length(nperiods)*length(npperiod);
mat=ones(tam,4+nop*3)*-1;
count=1;
for count11=1:(length(per1)-1)
    for count12=(count11+1):length(per2)
        for count2=1:length(nperiods)
            for count3=1:length(npperiod)
                [xdata,ydata]=gen_input([1 1],[per1(count11) per2(count12)],nperiods(count2),npperiod(count3),str2double(noise));
                if str2double(mode)==1
                    res=autocorrelation(xdata,ydata,nop,10000,'perms',[],1);
                else if str2double(mode)==2
                    per_range=max(xdata)-min(xdata);    
                    res=chi2(xdata,ydata,min(xdata)+floor(per_range/2),nop,10000,'perms',[],1);
                    else
                        per_range=max(xdata)-min(xdata);
                        res=dft(xdata,ydata,npperiod(count3)/min([per1(count11) per2(count12)]),floor(per_range/2),nop,10000,'perms',[],1);
                    end
                end
                mat(count,:)=[per1(count11) per2(count12) nperiods(count2) npperiod(count3) res];
                count=count+1;
            end
        end
    end
end

if str2double(mode)==1
    mod='auto';
else
    if str2double(mode)==2
        mod='enright';
    else
        mod='dft';
    end
end

ftxt11=strcat('mat_',mod,'_2pers','_per1_',min_per1,'_',max_per1,'_noise_',noise,'.mat');
save(ftxt11,'mat');

end
