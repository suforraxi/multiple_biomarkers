% USAGE
%    [Z,r,A] = ExtractPoles(S,n)
% FUNCTION:
%     gives the complex poles from auto-regressive analysis
% INPUT: 
%  S - signal [channel, sample]
%  n - requested order
%  OUTPUT:
%    Z - complex poles of the transfer function Z=exp(-u+2*Pii*w)
%    r - residual variation (noise)
%    A - complex amplitudes 
% VERSION:
%   Stiliyan 20.02.2013

function [Z,r,A] = ExtractPoles(S,n)
if nargin<2
    Z=[]; 
    r=[]; 
    A=[]; 
    return;
end;
for c=1:size(S,1)
 [B,r(c)]=arcov(S(c,:),n); 
 Z(:,c)=roots(B);
 R=Z(:,c)*ones(1,size(Z,1)); 
 R=R-transpose(R)+eye(size(Z,1)); 
 A(c,:)=1./prod(R)/B(1);
end;
Z=Z'; 
 return
end