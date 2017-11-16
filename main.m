addpath('lib')
segd = 10;
rhov = [0.8 1 1.2 1.6 3.2];

mb = 500;
sb = 0.1*mb;

T = 240;

p = 30;
q = p+10;

cvarv = 10.^linspace(-1,0.7,100);

mpst = zeros(length(cvarv),length(rhov));
mUneg = zeros(length(cvarv),length(rhov));
mUpos = zeros(length(cvarv),length(rhov));
errv = zeros(length(cvarv),length(rhov));

for k=1:length(cvarv)
    cvar = cvarv(k)
    for j=1:length(rhov)
        mv = rhov(j)*mb;

        [pst, Uneg, Upos, err] = ...
            iterateBuffer(T, segd, mb, sb, mv, cvar*mv, 0, p, q, 0);

        mpst(k,j) = mean(pst);
        mUneg(k,j) = mean(Uneg);
        mUpos(k,j) = mean(Upos);
        errv(k,j) = err;
    end
end
%%
figure(52);clf;box on;hold all;
plot(cvarv,-mUneg,'LineWidth',2,'color','black')
xlabel('coefficient of variation of bandwidth c_\lambda')
ylabel('average stalling duration L [s]')
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 10],'xticklabel',{'0.1','1','10'})
%%
figure(53);clf;box on;hold all;
plot(cvarv,mpst,'LineWidth',2,'color','black')
xlabel('coefficient of variation of bandwidth c_\lambda')
ylabel('average stalling probability p_{st}')
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 10],'xticklabel',{'0.1','1','10'})
%%
figure(54);clf;box on;hold all;
plot(cvarv,1+4*qoe(0.15,0.2,-mUneg,mpst*N),'LineWidth',2,'color','black')
xlabel('coefficient of variation of bandwidth c_\lambda')
ylabel('QoE value')
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 10],'xticklabel',{'0.1','1','10'})
%%
figure(55);clf;box on;hold all;
plot(cvarv,mpst./segd*60,'LineWidth',2,'color','black')
xlabel('coefficient of variation of bandwidth c_\lambda')
ylabel('rate of stalling events \phi [min^{-1}]')
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 10],'xticklabel',{'0.1','1','10'})
%%
figure(55);clf;box on;hold all;
plot(cvarv,mUpos,'LineWidth',2,'color','black')
xlabel('coefficient of variation of bandwidth c_\lambda')
ylabel('average buffer level [s]')
set(gca,'xscale','log')
set(gca,'xtick',[0.1 1 10],'xticklabel',{'0.1','1','10'})