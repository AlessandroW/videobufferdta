function y=mirrorY(B);
    y.prob=B.prob(end:-1:1);
    y.x=-(B.x(end:-1:1));
    y.prob=y.prob./sum(y.prob);
end