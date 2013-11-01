function evaluate(min_per,max_per,min_nperiods,max_nperiods,min_npperiod,max_npperiod,noise,mode)

per=str2double(min_per):str2double(max_per);
nperiods=str2double(min_nperiods):str2double(max_nperiods);
npperiod=str2double(min_npperiod):str2double(max_npperiod);
nop=5;

tam=length(per)*length(nperiods)*length(npperiod);
mat=ones(tam,3+nop*3)*-1;
count=1;
for count1=1:length(per)
    for count2=1:length(nperiods)
        for count3=1:length(npperiod)
            [xdata,ydata]=gen_input([1],[per(count1)],nperiods(count2),npperiod(count3),str2double(noise));
            if str2double(mode)==1
                res=autocorrelation(xdata,ydata,nop,10000,'perms',[],1);
            else if str2double(mode)==2
                    per_range=max(xdata)-min(xdata);
                    res=chi2(xdata,ydata,min(xdata)+floor(per_range/2),nop,10000,'perms',[],1);
                else
                    per_range=max(xdata)-min(xdata);
                    res=dft(xdata,ydata,npperiod(count3)/per(count1),floor(per_range/2),nop,10000,'perms',[],1);
                end
            end
            mat(count,:)=[per(count1) nperiods(count2) npperiod(count3) res];
            count=count+1;
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

ftxt11=strcat('mat_',mod,'_per_',min_per,'_',max_per,'_noise_',noise,'.mat');
save(ftxt11,'mat');

end
