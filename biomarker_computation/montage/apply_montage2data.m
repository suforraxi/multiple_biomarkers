% apply a grid montage (passed as a function) 
% data    - fieldtrip data structure
% fun     - handler to one of the possible montage function
%           the signature of the generic function takes 
%           two input parameters (channel label [cell], matrix of data [ch X samples])
%        
% outdata - fieldtrip data structure with the data transformed according to the montage funtion fun 
function outdata = apply_montage2data(data,fun)

outdata = data;
ntrial = numel(data.trial);
for i = 1: ntrial
    
    [outdata.label,outdata.trial{i}] = fun(data.label,data.trial{i});
    if(isempty(outdata.trial{i}))
        outdata = [];
        return
    end
end
