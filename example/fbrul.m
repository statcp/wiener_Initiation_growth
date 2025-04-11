function fx = fbrul(R,wrul,mu0,sig02,a,b,sigb2)
% function offtu
% solving monotone increasing function
% mu0 = mu0G(j),sig0 = sig0G(j),a=aG(j),b=bG(j),sigb = sigbG(j);
% R=R1, wrul = wrul(ii)

uprul = frul(10000,wrul,mu0,sig02,a,b,sigb2);
while uprul < R
    R = uprul;
end

    lb = 0; up = 20;% lower bound and upper bound
    yb = frul(up,wrul,mu0,sig02,a,b,sigb2)-R;
    
     %找右确界
    while yb <= 0
        up = up+10;
        yb = frul(up,wrul,mu0,sig02,a,b,sigb2)-R;
    end
    
    while abs(up-lb) > .000001
        c = (lb+up)/2;
        yc = frul(c,wrul,mu0,sig02,a,b,sigb2)-R;
        if yc > 0
            up = c;
        else
            lb = c;
        end
    end
    fx= (lb+up)/2;
    
end


