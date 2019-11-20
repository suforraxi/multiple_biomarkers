% apply a montage passed as a function
% data - fieldtrip data structure (.vhdr as data format)
% fun  - handler to one of the possible montage function
%        the signature of the generic function takes 
%        two input parameters (channel label [cell], matrix of data [ch X samples])
%        and return 
%        two output paramenters (channel label [cell], matrix tranformed according to the montage of data [ch X samples] )           
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
