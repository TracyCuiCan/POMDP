% this function finds cross points between all lines in set C, and only
% keeps the points between [0,1]. 
function d=node(C)
n=size(C,1);
m=1;
for i=1:n
    for j=(i+1):n
        d(m)=intersec([0, C(i,1)],[1,C(i,2)],[0,C(j,1)],[1,C(j,2)]);
        if d(m)<0
            d(m)=0;
        else if d(m)>1
                d(m)=1;
            end
        end
        m=m+1;
    end
end
d=unique(d);
end