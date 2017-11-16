function z=pi_op(x,m)
    y=x;
    if(min(x.x)<m)
    z.prob=y.prob(y.x>=m);
    z.prob(1)=z.prob(1)+sum(y.prob(y.x<m));
    z.x=y.x(y.x>=m);
    
    else
        z=x;
    end;
%     i=find(y.x>=m); % index
%     z.x = [m y.x(i)];
%     z.prob = [sum(y.prob(1:i(1)-1)) y.prob(i)];
end





