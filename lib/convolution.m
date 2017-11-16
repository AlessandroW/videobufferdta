function z=convolution(x,y)
        
    %z.x=x.x(1)+y.x(1):y.x(2)-y.x(1):x.x(end)+y.x(end);
    z.prob=conv(x.prob,y.prob);
    z.x=linspace(x.x(1)+y.x(1),x.x(end)+y.x(end),length(z.prob));
end