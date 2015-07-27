% This function is to calculate the distance between two points

function dist=mydistance(p1,p2)
lg=length(p1);
temp1=0;
for i=1:lg
    temp1=temp1+(abs(p1(i)-p2(i)))^2;
end
dist=sqrt(temp1);