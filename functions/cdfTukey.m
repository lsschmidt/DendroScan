function prtrng=cdfTukey(q,v,r)
% The function returns the cumulative distribution function of the
% studentized range distribution used in Tukey's HSD test
% q - observed value of Tukey's test
% v - degrees of freedom (total number of elements minus number of groups
% r - number of samples
pcutj=0.00003;
pcutk=0.0001;
step=0.1125; % 0.45 in the original Fortran file
vmax=120;
cv1=0.193064705;
cv2=0.293525326;
cvmax=0.39894228;
cv=[0.318309886,-0.268132716e-2, 0.347222222e-2, 0.833333333e-1];
jmin=12; % 3 in the original Fortran file
jmax=60; % 15 in the original Fortran file
kmin=28; % 7 in the original Fortran file
kmax=60; % 15 in the original Fortran file
vw=zeros(121); % 30 in the original Fortran file
qw=zeros(121); % 30 in the original Fortran file
% Minimum and maximum number of steps are controlled by
% jmin, jmax, kmin and kmax.  Accuracy can be increased
% by use of a finer grid - Increase sizes of arrays vw
% and qw, and jmin, jmax, kmin, kmax and 1/step proportionally.
prtrng=0;
if not(q<=0 || v<1 || r<2)
    g=step*r^-0.2;
    gmid=0.5*log(r);
    r1=r-1;
    c=log(r*g*cvmax);
    if v<=vmax
        h=step*v^-0.2;
        v2=v*0.5;
        switch v
            case 1
                c=cv1;
            case 2
                c=cv2;
            otherwise
                c=sqrt(v2)*cv(1)/(1+((cv(2)/v2+cv(3))/v2+cv(4))/v2);
        end
        c=log(c*r*g*h);
    end
    gstep=g;
    qw(1)=-1;
    qw(jmax+1)=-1;
    pk1=1;
    pk2=1;
    k=1;
    while k<=kmax % loop 28
        gstep=gstep-g;
        notDone21=true;
        while notDone21 % loop 21
            gstep=-gstep;
            gk=gmid+gstep;
            pk=0;
            if not(pk2<=pcutk && k>kmin) % test 26
                w0=c-gk*gk*0.5;
                pz=1-normcdf(gk);
                x=1-normcdf(gk-q)-pz;
                if x>0, pk=exp(w0+r1*log(x));end
                if not(v>vmax) % test 26 again
                    jump=-jmax;
                    notDone22=true;
                    while notDone22 % loop 22
                        jump=jump+jmax;
                        j=1;
                        while j<=jmax % loop 24
                            jj=j+jump;
                            if not(qw(jj)>0) % test 23
                                hj=h*j;
                                if j>jmax, qw(jj+1)=-1;end
                                ehj=exp(hj);
                                qw(jj)=q*ehj;
                                vw(jj)= v*(hj+0.5-ehj*ehj*0.5);
                            end % test 23
                            pj=0;
                            x=1-normcdf(gk-qw(jj))-pz;
                            if x>0, pj=exp(w0+vw(jj)+r1*log(x)); end
                            pk=pk+pj;
                            if not(pj>pcutj) % test 24
                                if jj>jmin || k>kmin % if conditions are met, make sure to exit loop 24
                                    j=jmax+1;
                                end
                            end
                            j=j+1;
                        end % loop 24
                        h=-h;
                        notDone22=h<0;
                    end % loop 22
                end
            end
            prtrng=prtrng+pk;
            if k>kmin && pk<=pcutk && pk1<=pcutk % if the conditions are met, the program will step out of loop 28 and loop 21
                k=kmax+1; % make sure to exit loop 28
                notDone21=false; % make sure to exit loop 21
            else
                pk2=pk1;
                pk1=pk;
                notDone21=gstep>0;
            end
        end % loop 21
        k=k+1;
    end % loop 28
end % if not(q<=0 || ifault==1)
    
