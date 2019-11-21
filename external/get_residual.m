% Calculate poles (Z), amplitudes (A), and residuals (r)
% with model order (n) and windowlength (L) with 50% overlap

% Stiliyan Kalitzin
function [Z,r,A] = get_residual(D,n,L)

for c=1:size(D,1)
  E=embed(D(c,:),1,L,round(L/2));
  for k=1:size(E,2)
     [Z(:,c,k),r(c,k),A(:,c,k)] = ExtractPoles(E(:,k)',n);
  end;
end;