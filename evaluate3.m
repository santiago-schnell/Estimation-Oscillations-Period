function [mat]=evaluate3(min_per,max_per,min_nperiods,max_nperiods,min_npperiod,max_npperiod)

per=str2double(min_per):str2double(max_per);
nperiods=str2double(min_nperiods):str2double(max_nperiods);
npperiod=str2double(min_npperiod):str2double(max_npperiod);

tam=length(per)*length(nperiods)*length(npperiod);
mat=ones(tam,3)*-1;
count=1;
for count1=1:length(per)
    for count2=1:length(nperiods)
        for count3=1:length(npperiod)
            mat(count,:)=[per(count1) nperiods(count2) npperiod(count3)];
            count=count+1;
        end
    end
end

end
