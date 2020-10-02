function [h,p,stats]=games_howell(data,group,alpha)
%
% GAMES_HOWELL: post-hoc Games-Howell test for one-way ANOVA.
%
% [h,p,stats]=games_howell(data,group,alpha) performs the Games-Howell test
% for one-way ANOVA (analysis of variance). The Games-Howell test compares
% all individual group means to each other in a pairwise fashion. It
% accommodates groups with unequal sample sizes (with a recommended minimum
% number of 6 observations in any given group) and variances.
%
% INPUTS:
%
% The 'data' input is a vector with individual observations, and the
% 'group' input is a vector with labels that indicate to which group each
% observation in 'data' belongs. 'group' can be a numerical, string or cell
% array. Both inputs are required. The syntax is similar to that of ANOVA1
% in the Statistics Toolbox. The 'alpha' input optionally lets the user
% select the significance level for statistical testing; the default is
% 0.05.
%
% OUTPUTS:
%
% The 'h' output is a numerical array that is set to 1 if the Games-Howell
% test rejects the null hypothesis of no difference between group mean
% pairs at the alpha significance level. The 'p' output is a numerical
% array that returns the probability of observing a value as extreme as or
% more extreme than the test statistic under the null hypothesis. 'h' and
% 'p' are square matrices whose dimension reflects the number of different
% in 'group'. Their diagonal is set to NaN. The 'stats' structure contains:
%  - gnames, the labels for the different groups
%  - md, the differences between pairs of groups
%  - q, the test statistic for each pairwise comparison
%  - df, the degrees of freedom for each pairwise comparison
%
% EXAMPLE:
%
% data=[1+2.*randn(9,1);1+randn(14,1);-1+randn(11,1);randn(7,1)];
% group=[ones(9,1);2.*ones(14,1);3.*ones(11,1);4.*ones(7,1)];
% [h,p,stats]=games_howell(data,group);
%
% REFERENCES:
%
% Games PA, Howell JF. Pairwise Multiple Comparison Procedures with Unequal
% N's and/or Variances: A Monte Carlo Study. Journal of Educational and
% Behavioral Statistics 1976;1:113. doi: 10.3102/10769986001002113
%
% I used the equations in the following web pages to prepare this function:
% - David C. Howell's Advanced Statistical Methods:
%   https://www.uvm.edu/~dhowell/gradstat/psych341/labs/Lab1/Multcomp.html
% - Jon Starkweather's Introduction to Statistics for the Social Sciences:
%   http://www.unt.edu/rss/class/Jon/ISSS_SC/Module009/isss_m91_onewayanova
%   /node7.html
%
% I tested the output against GNU PSPP 0.8.3 for Windows
% (http://www.gnu.org/www.gnu.org/software/pspp/)
%
% This function calls the function CDFTUKEY, by Peter Nagy. 
% http://www.mathworks.com/matlabcentral/fileexchange/37450
% 
% See also ANOVA1, CDFTUKEY.
% Author of this function: pierre.megevand@gmail.com
% perform some checks on the inputs
if nargin<2
    error('The ''data'' and ''group'' inputs need to be defined. Type ''help games_howell'' for more information.');
end
if nargin<3
    alpha=0.05; % default alpha level
end
if ~isvector(data)||~isvector(group)
    error('The ''data'' and ''group'' inputs need to be vectors. Type ''help games_howell'' for more information.');
end
if numel(data)~=numel(group)
    error('The number of elements in ''data'' and ''group'' need to be the same. Type ''help games_howell'' for more information.');
end
[~,gnamesidx]=unique(group);
gnames=group(sort(gnamesidx));
ng=numel(gnames);
[~,gidx]=ismember(group,gnames);
h=NaN(ng);
p=NaN(ng);
md=NaN(ng);
q=NaN(ng);
df=NaN(ng);
for n1=1:ng
    for n2=1:ng
        if n1~=n2
            md(n1,n2)=mean(data(gidx==n1))-mean(data(gidx==n2));
            q(n1,n2)=compq(data(gidx==n1),data(gidx==n2));
            df(n1,n2)=compdf(data(gidx==n1),data(gidx==n2));
            
            p(n1,n2)=(1-cdfTukey(abs(q(n1,n2)),df(n1,n2),ng));

            if p(n1,n2)<=alpha
                h(n1,n2)=1;
            else
                h(n1,n2)=0;
            end
        end
    end
end
% compute test statistic
    function q=compq(g1,g2)
        se=sqrt((var(g1)./numel(g1)+var(g2)./numel(g2))./2);
        q=(mean(g1)-mean(g2))./se;
    end
% compute degrees of freedom
    function df=compdf(g1,g2)
        df=((var(g1)./numel(g1)+var(g2)./numel(g2)).^2) ./ ...
            (((var(g1)./numel(g1)).^2)./(numel(g1)-1) + ...
            ((var(g2)./numel(g2)).^2)./(numel(g2)-1));
    end
% prepare stats structure
stats=[];
stats.gnames=gnames;
stats.md=md;
stats.q=q;
stats.df=df;
end