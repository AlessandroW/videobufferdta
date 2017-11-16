function z=pi_op_minto(x,m,to)
    if(x.x(1)<m && ~(to<m))
        if x.x(end) < to
            dx = x.x(2)-x.x(1);
            n = ceil((to - x.x(end))/dx);
            x.x(end+1:end+n)=x.x(end)+(1:n)*dx;
            x.prob(end+1:end+n)=zeros(1,n);
        end
        z.prob=x.prob;
        z.x=x.x;
        z.prob(x.x<m)=0;
        toind = find(z.x<=to,1,'last');
        z.prob(toind)=z.prob(toind)+sum(x.prob(x.x<m));
        nonzero = z.prob>0;
        z.prob(1:(find(nonzero,1,'first')-1)) = [];
        z.x(1:(find(nonzero,1,'first')-1)) = [];
    else
        z=x;
    end;
%     i=find(y.x>=m); % index
%     z.x = [m y.x(i)];
%     z.prob = [sum(y.prob(1:i(1)-1)) y.prob(i)];
end