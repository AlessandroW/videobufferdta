function [pst Uneg Upos err] = iterateBufferA(T, segd, ma, sa, p, q, init)

if T<Inf
    N = ceil(T/segd);

    pst = nan(N-1,1);
    Uneg = nan(N-1,1);
    Upos = nan(N-1,1);
end
    err = nan;
    
    overhead=@(x)10.6082*exp(-0.059*x);
    overhead=@(x)9.2308./x+0.7692;
    overhead=@(x)1;
    %overhead=@(x)4.1447*exp(-0.0355*x);
        
    mtarget = ma;
    starget = sa;
    
%     vd = starget*starget;
%     md = mtarget;
%     mud = log((md^2)/sqrt(vd+md^2));
%     sigmad = sqrt(log(vd/(md^2)+1));
%     %d.x = linspace(logninv(0.0001,mud,sigmad),logninv(0.9999,mud,sigmad),1e4);
%     %d.x = linspace(0,max(100,logninv(0.99999,mud,sigmad)),1e4);
%     d.x = linspace(0,200,1e4);
%     d.prob = lognpdf(d.x,mud,sigmad);
%     d.prob = d.prob/sum(d.prob);

%%%%%%

% [mv, sv] = rbmapget(mtarget, starget, 0.05);
% 
% if isnan(mv)
%     return;
% end
% 
%     d.x = linspace(0,200,1e4);
%     d.prob = ratio_of_normalpdf(d.x, mb*segd, sb*segd, mv, sv);
%     d.prob=d.prob./sum(d.prob);

%%%%%%

%     err = d.x*d.prob'-mb*segd/mv;
%     if sqrt(err^2) > 0.05;
%         return
%     end
    %err = 1e3;
    
%     m = mv;
%     s = sv;   
%     while sqrt(err^2)/(mb*segd/mv) > 0.05
%         %d.x = linspace(max(0,mb-2*sb)*segd/(mv+2*sv),(mb+2*sb)*segd/max(1e-3,mv-2*sv),1e3);
%         d.x = linspace(0,200,4e3);
%         d.prob = ratio_of_normalpdf(d.x, mb*segd, sb*segd, m, s);
%         %d = ratio_of_lognormalpdf(d.x, mb*segd, sb*segd, mv, sv);
%         d.prob=d.prob./sum(d.prob);
%         %d.prob(end)=1-sum(d.prob(1:end-1));
%         
%         err = d.x*d.prob'-mb*segd/mv;
%         m = m + err;
%     end

     m = mtarget;
     s = starget;   
     while isnan(err) || sqrt(err^2)/(ma) > 0.001
              vd = starget*starget;
              md = m;
              mud = log((md^2)/sqrt(vd+md^2));
              sigmad = sqrt(log(vd/(md^2)+1));
              d.x = linspace(0,240,2e4);
              d.prob = lognpdf(d.x,mud,sigmad);
              d.prob = d.prob/sum(d.prob);
              err = d.x*d.prob'-ma;
              m = m + err;
     end
    if sqrt(((d.x.^2*d.prob')-(d.x*d.prob').^2)-starget*starget) > 1e-3
        disp('Warning: error stdd of segment IA > 1e-3')
    end
%     m=mtarget;
%     while isnan(err) || sqrt(err^2)/(mb*segd/mv) > 0.001
%         P=m*(1e4/200)/(starget*1e4/200)^2;
%         R=m*1e4/200*P/(1-P);
%         d.x=0:1:1e4;
%         d.prob=nbinpdf(d.x,R,P);
%         d.x=d.x*200/1e4;
%         d.prob = d.prob/sum(d.prob);
%         err = d.x*d.prob'-mb*segd/mv;
%         m = m + err;
%     end

       
%     figure(111);clf;box on;hold all;
%     plot(d.x,d.prob);
%     line([d.x*d.prob' d.x*d.prob'], [0 max(d.prob)]);
%     line([mb*segd/mv mb*segd/mv], [0 max(d.prob)],'LineStyle','--');
    
    if d.prob(1) > 1e-3 || d.prob(end) > 1e-3
        disp('Warning: definition space of ratio pdf exceeded')
    end

  
    b_n.x = [0 (floor(init/segd)+1)*segd];
    b_n.prob = [0 1];
    
%     figure(1);clf;stairs(b_n.x,cumsum(b_n.prob),'color','red');
%     %figure(1);clf;plot(b_n.x,b_n.prob,'o','color','red');
%     xlim([-30 50]);ylim([0 1]);
%     xlabel('k');ylabel('F_{u_n}(k)');
    eps = 1e-9;
    problast = 0;
    it = 0;
    maxit=Inf;
    if T<Inf
        maxit=N-1;
        eps = -1;
    end        
   
%    while it<maxit
while it<maxit && (length(problast) ~= length(b_n.prob) || ((b_n.prob - problast)*(b_n.prob - problast)') > eps)%((sum(b_n.prob(b_n.x<=0)) - pstlast)*(sum(b_n.prob(b_n.x<=0)) - pstlast)') > eps*eps
%         if (length(problast) == length(b_n.prob))
%             (b_n.prob - problast)*(b_n.prob - problast)'
%         else
%             length(b_n.prob)
%         end
        it=it+1;
        problast = b_n.prob;
        pstlast = sum(b_n.prob(b_n.x<0));
        Uneglast = b_n.prob(b_n.x<0)*b_n.x(b_n.x<0)';
        %size(problast)

        tmp = b_n.prob(b_n.x>0)*b_n.x(b_n.x>0)';
        
                %length(b_n.prob)
        %figure(1);plot(b_n.x,b_n.prob,'o');
        b_n_l = b_n;
        b_n_l.prob(b_n.x>q) = 0;
        b_n_h = b_n;
        b_n_h.prob(b_n.x<=q) = 0;
        
        %b_n = convolution(b_n_l,b);
        b_n = b_n_l;
        
        probh = sum(b_n_h.prob);
        if(probh>0)
            toind = find(b_n.x<=p,1,'last');
            b_n.prob(toind)=b_n.prob(toind)+probh;
        end
        
        %TODO
        %[b_ntmp.x, iu, il] = union(upper.x, lower.x);

        %b_n = convolution(b_n, mirrorY(d));
        b_n = convolution(b_n, mirrorY(d));
        ind0 = find(b_n.prob>0,1,'last');
        if (~isempty(ind0) && ind0<length(b_n.x))
            b_n.x(ind0+1:end) = [];
            b_n.prob(ind0+1:end) = [];
        end
        
        if T<Inf
            pst(it) = sum(b_n.prob(b_n.x<0));
            Uneg(it) = b_n.prob(b_n.x<0)*b_n.x(b_n.x<0)';
            Upos(it) = (segd/(-Uneg(it)+segd))*(tmp + b_n.prob(b_n.x>0)*b_n.x(b_n.x>0)')/2;
        else
            pst = sum(b_n.prob(b_n.x<0));
            Uneg = b_n.prob(b_n.x<0)*b_n.x(b_n.x<0)';
            Upos = (segd/(-Uneg+segd))*(tmp + b_n.prob(b_n.x>0)*b_n.x(b_n.x>0)')/2;
        end
        %b_n_d = mirrorY(b_n);
        
        %Variant 1: start playback when buffer is positive
%        figure(1);clf;stairs(b_n.x,cumsum(b_n.prob),'color','red');hold all
        %figure(1);clf;plot(b_n.x,b_n.prob,'o','color','red');
%         xlim([-30 50]);ylim([0 1]);
%         xlabel('k');ylabel('F_{u_n}(k)');
        %b_n=pi_op(b_n,0);
        b_n=pi_op_minto(b_n,0,max(0,floor(init/segd)*segd));
        %figure(1);plot(b_n.x,b_n.prob,'o','color','red');
%        stairs(b_n.x,cumsum(b_n.prob),'color','blue');
        
        if T<Inf && it == (N-1)
            b.x = [0 T-(N-1)*segd];
            b.prob = [0 1];
        else
            b.x = [0 segd];
            b.prob = [0 1];
        end
        
        b_n = convolution(b_n,b);
        
%          hold all;
%          figure(1);stairs(b_n.x,cumsum(b_n.prob),'color','blue');
        %figure(1);plot(b_n.x,b_n.prob,'o','color','blue');
        

    
    %Upos(it) = b_n.prob(b_n.x>0)*b_n.x(b_n.x>0)';
    %Upos(it) = tmp;
end
    %figure(2);plot(b_n_d.x(b_n_d.x>0),cumsum(b_n_d.prob(b_n_d.x>0)/sum(b_n_d.prob(b_n_d.x>0))),'color','blue');
end