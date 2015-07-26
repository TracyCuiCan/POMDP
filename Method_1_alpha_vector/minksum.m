% This function finds the Minkowski sum of two sets
function S=minksum(A,B)
if nargin<2
    error('Requires two inputs');
end
szA=size(A); szB=size(B);
if (szA(2)~=szB(2)) || numel(szA)>2 || numel(szB)>2
    error('Input arrays must be MxN arrays, with equal N');
end

%Pre-allocate memory
S=zeros(szA(1)*szB(1),szA(2));

%Loop through each dimension
for k=1:szA(2)
   temp=bsxfun(@plus,A(:,k),B(:,k)');
   S(:,k)=temp(:);
end
%Take unique for each row
S=unique(S,'rows');
end


