 % USAGE:
%   [C,n,p]=h2_m(S1,S2,n,N);
% FUNCTION:
%   gives the non-linear assymetric association h^2 index b/n sig S1 and S2 
%   H2(S1 <--| S2);var(S1|S2)/var(S1)
% INPUT:
%    S1, S2: 1D or 2D signals, [channel, time]
%    n : when S2 is 1D:
%           n>1, number of bins for S2, 
%           if 0<n<=1, the bin size is that fraction of the std(S2)
%           if n<0 integer clustered data from 1 to N is assumed
%           if n is a vector - indicates the bin-borders (not implemented) 
%    N : default []
%        if []: max(bin)/sum(bin)-1/n: the simplified bin entropy
%        if > 1 :number of iterations for p-value 
%        if <1 : the prctile value computed at [1000] itterations 
%        if N=[N1 N2] then N1 is th e p-value and N2 the # of itterations
% OUTPUT:
%  C - h^2 association index = (total.variance-explained.variance)/total.variance
%  n - the resulting number of bins
%  p - depending on N: if []: max(bin)/sum(bin)-1/n: the simplified bin entropy;
%                      N>1 the p is the p-palue of C
%                      N<1 the N-th percentile (for 1000 itteratons) 
%                      N=[N1<1 N2>1] is the N1 prctile for N2 itterations
%                      
% VERSION:
%  Stiliyan Kalitzin 13.08.18

function [C,n,p]=h2_m(S1,S2,n,N)
S1=reshape(S1,[1 numel(S1)]); 
S2=reshape(S2,[1 numel(S2)]); 
s=std(S2); 
C=[];
p=[];
if s==0
    n=[]; 
    return;
end
if nargin<4 
    N=[]; 
end
d1=size(S1); 
d2=size(S2);
n1=d1(2); 
n2=d2(2); 
m=min(n1,n2);
if  m==0 
    return; 
end
 S1=detrend(S1(:,1:m)','constant')';
 S2=detrend(S2(:,1:m)','constant')';
if nargin<3 
    n=[]; 
end
if isempty(n)
    n=0.1; 
end
if (n>0)&&(n<=1) 
  n=round((max(S2)-min(S2))/(n*s)); 
end 
[C,p]=h2_local(S1,S2,n); 
if nargout<3 
    return; 
end
if isempty(N)
    return;
end
if length(N)>1 
    N1=N(2); 
    p1=N(1); 
elseif N>1 
    N1=N; 
else
    N1=1000; 
    p1=N;
end
C0=zeros(1,N1); 
for k=1:N1
    J=randperm(m); 
    C0(k)=h2_local(S1(J),S2,n); 
end 
if N(1)>1 
    p=sum(C0>=C)/N1;
else
    p=prctile(C0,(1-p1)*100); 
end 

function [C,u]=h2_local(S1,S2,n) 
if n>0 
    S2=floor(n*(S2-min(S2))/(max(S2)-min(S2)))+1;
    S2(S2>n)=n;
else 
    n=max(S2); 
end
D=zeros(1,n); 
M=zeros(1,n); 
N=zeros(1,n); 
for k=1:length(S2) 
  D(S2(k))=D(S2(k))+S1(k).^2;
  M(S2(k))=M(S2(k))+S1(k);
  N(S2(k))=N(S2(k))+1; 
end
D(N>0)=D(N>0)-M(N>0).^2./N(N>0); 
D=sum(D); 
C=1-D/sum(sum(S1.^2)); 
u=max(N)/sum(N)-1/n; 




