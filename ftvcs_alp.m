function [U, out] = ftvcs_alp(A,b,p,q,opts)
%
% Goal: solve   min sum ||D_i u|| + Ohm||Fu||1   (with or without the constraint u>=0) 
%                  s.t. Au = b
%       to recover image/signal u from encoded b,
%       which is equivalent to solve       min sum ||w_i||+Ohm||w_f||
%                                          s.t. D_i u = w_i
%                                               F u = w_f
%                                               Au = b
% ftvcs_al solves Augmented Lagrangian function:
% 
% min_{u,w} sum ||w_i||+Ohm||w_f|| - sigma'(Du-w) - delta'(Au-b) - Ohm(Sigmaw'(Fu-w_f))
%                   + beta/2 ||Du-w||_2^2 + mu/2||Au-b||_2^2 +Ohm*upsilon/2||Fu-w_f)||_2^2,
%_________________________________
% by an alternating algorithm:
% i)  while norm(up-u)/norm(up) > tol_inn
%     1) Fix w^k, do Gradient Descent to 
%            - sigma'(Du-w^k) - delta'(Au-b) + beta/2||Du-w^k||^2 + mu/2||Au-f||^2;
%            u^k+1 is determined in the following way:
%         a) compute step length tau > 0 by BB formula
%         b) determine u^k+1 by
%                  u^k+1 = u^k - alpha*g^k,
%            where g^k = -D'sigma - A'delta + beta D'(Du^k - w^k) + mu A'(Au^k-f), 
%            and alpha is determined by Amijo-like nonmonotone line search;
%     2) Given u^k+1, compute w^k+1 by shrinkage
%                 w^k+1 = shrink(Du^k+1-sigma/beta, 1/beta);
%     end
% ii) update Lagrangian multipliers by
%             sigma^k+1 = sigma^k - beta(Du^k+1 - w^k+1)
%             delta^k+1 = delta^k - mu(Au^k+1 - b).
% iii)accept current u as the initial guess to run the loop again
%
% Inputs:
%       A        : either an matrix representing the measurement or a struct 
%                  with 2 function handles:
%                           A(x,1) defines @(x) A*x;
%                           A(x,2) defines @(x) A'*x;
%       b        :  either real or complex input vector representing the
%                   noisy observation of a grayscale image
%       p, q     :  size of original image
%       opts     :  structure to restore parameters
%
%
% variables in this code:
%
% lam1 = sum ||wi||
% lam2 = ||Du-w||^2 (at current w).
% lam3 = ||Au-f||^2
% lam4 = sigma'(Du-w)
% lam5 = delta'(Au-b)
% lam6 = Ohm*||w_f||
% lam7 = Ohm*sigmaw'(Fu-w_f)
% lam8 = Ohm*||Fu-w_f||^2
%
%   f  = lam1 + beta/2 lam2 + mu/2 lam3 - lam4 - lam5 +lam6 +upsilon/2
%   lam8 - lam7
% 
%
%   g  = A'(Au-f)
%   g2 = D'(Du-w) (coefficients beta and mu are not included)
%
%
% Numerical tests illustrates that this solver doestn't require large beta
% and mu. ( <100 usually)
%
%
% Written by: Chengbo Li
% Advisor: Prof. Yin Zhang and Wotao Yin
% Computational and Applied Mathematics department, Rice University
% May. 2, 2009
%------------------
% Revised by Chenyang 


global D Dt
[D,Dt] = defDDt;

% problem dimension
n = p*q;

% unify implementation of A
if ~isa(A,'function_handle')
    A = @(u,mode) f_handleA(A,u,mode); 
end

% get or check opts
opts = ftvcs_al_opts(opts); 

% mark important constants
mu = opts.mu;
beta = opts.beta;
upsilon = opts.upsilon;
tol_inn = opts.tol_inn;
tol_out = opts.tol;
gam = opts.gam;
Ohm = opts.Ohm;
% check if A*A'=I
tmp = rand(length(b),1);
if norm(A(A(tmp,2),1)-tmp,1)/norm(tmp,1) < 1e-3
    opts.scale_A = false;
end
clear tmp;

% check scaling A
if opts.scale_A
    [mu,A,b] = ScaleA(n,mu,A,b,opts.consist_mu);
end 

% check scaling b
if opts.scale_b
    [mu,b,scl] = Scaleb(mu,b,opts.consist_mu);
end

% calculate A'*b
Atb = A(b,2);

% initialize U, beta
muf = mu;
betaf = beta;     % final beta
upsilonf = upsilon;
[U,mu,beta,upsilon] = ftvcs_al_init(p,q,Atb,scl,opts);    % U: p*q
if mu > muf; mu = muf; end
if beta > betaf; beta = betaf; end
if upsilon > upsilonf; upsilon = upsilonf; end


muDbeta = mu/beta;                 %????????
upDbeta = upsilon/beta;

rcdU = U;

% initialize multiplers
sigmax = zeros(p,q);                       % sigmax, sigmay: p*q 
sigmay = zeros(p,q);
sigmaf = zeros(p,q);

delta = zeros(length(b),1);                % delta: m

% initialize D^T sigma + A^T delta +ohm*F^T sigmaf 
DtsAtd = zeros(p*q,1); 

% initialize out.n2re
if isfield(opts,'Ut')
    Ut = opts.Ut*scl;        %true U, just for computing the error
    nrmUt = norm(Ut,'fro');
else
    Ut = []; 
end
if ~isempty(Ut)
    out.n2re = norm(U - Ut,'fro')/nrmUt; 
end

% prepare for iterations
out.mus = mu; out.betas = beta; out.upsilons = upsilon;
out.res = []; out.itrs = []; out.f = []; out.obj = []; out.reer = [];
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.lam6 = [];out.lam7 = [];out.lam8 = [];
out.itr = Inf;
out.tau = []; out.alpha = []; out.C = []; gp = [];
out.cnt = [];

[Ux,Uy] = D(U);                   % Ux, Uy: p*q
Uf = fft2c(U);


if opts.TVnorm == 1
    Wx = max(abs(Ux) - 1/beta, 0).*sign(Ux);
    Wy = max(abs(Uy) - 1/beta, 0).*sign(Uy);
    Wf = max(abs(Uf) - 1/upsilon, 0).*sign(Uf);
    lam1 = sum(sum(abs(Wx) + abs(Wy)));
    lam6 = Ohm*sum(sum(abs(Wf)));
else
    V = sqrt(Ux.*conj(Ux) + Uy.*conj(Uy));        % V: p*q
    V(V==0) = 1;
    S = max(V - 1/beta, 0)./V;        % S: p*q
    Wx = S.*Ux;                       % Wx, Wy: p*q
    Wy = S.*Uy;
    lam1 = sum(sum(sqrt(Wx.*conj(Wx) + Wy.*conj(Wy))));  
end  

% [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,...
%     Atb,sigmax,sigmay,delta);
 [lam2,lam3,lam4,lam5,lam7,lam8,f,g2,g3,Au,g] = get_g(U,Ux,Uy,Uf,Wx,Wy,Wf,lam1,...
  lam6,beta,mu,upsilon,A,b,Atb,sigmax,sigmay,sigmaf,delta,Ohm);
%lam, f: constant      g2: pq        Au: m         g: pq

% compute gradient
d = g2 + muDbeta*g - DtsAtd + Ohm*upDbeta*g3;

count = 1;
Q = 1; C = f;                     % Q, C: costant
out.f = [out.f; f]; out.C = [out.C; C];
out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5]; out.lam6 = [out.lam6; lam6];
out.lam7 = [out.lam7; lam7]; out.lam8 = [out.lam8; lam8];

for ii = 1:opts.maxit
    if opts.disp
        fprintf('outer iter = %d, total iter = %d, normU = %4.2e; \n',count,ii,norm(U,'fro'));
    end
    
    % compute tau first
    if ~isempty(gp)
        dg = g - gp;                        % dg: pq
        dg2 = g2 - g2p;                     % dg2: pq
        dg3 = g3 - g3p;
        ss = uup'*uup;                      % ss: constant
        sy = uup'*(dg2 + muDbeta*dg + Ohm*upDbeta*dg3);       % sy: constant   ??????Ohm??????
        % sy = uup'*((dg2 + g2) + muDbeta*(dg + g));
        % compute BB step length
        tau = abs(ss/max(sy,eps));               % tau: constant
        
        fst_itr = false;
    else
        % do Steepest Descent at the 1st ieration
        %d = g2 + muDbeta*g - DtsAtd;         % d: pq
        [dx,dy] = D(reshape(d,p,q));                    %dx, dy: p*q
        df = fft2c(reshape(d,p,q));
        dDd = norm(dx,'fro')^2 + norm(dy,'fro')^2 + Ohm*upDbeta*norm(df,'fro')^2;      % dDd: cosntant  ????
        Ad = A(d,1);                        %Ad: m
        % compute Steepest Descent step length
        tau = abs((d'*d)/(dDd + muDbeta*Ad'*Ad));
        
        % mark the first iteration 
        fst_itr = true;
    end    
    
    % keep the previous values
    Up = U; gp = g; g2p = g2; g3p = gp; Aup = Au; Uxp = Ux; Uyp = Uy; Ufp = Uf;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONE-STEP GRADIENT DESCENT %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    taud = tau*d;
    U = U(:) - taud;
    % projected gradient method for nonnegtivity
    if opts.nonneg
        U = max(real(U),0);
    end
    U = reshape(U,p,q);                    % U: p*q (still)
    [Ux,Uy] = D(U);                        % Ux, Uy: p*q
    Uf = fft2c(U);
    
    [lam2,lam3,lam4,lam5,lam7,lam8,f,g2,g3,Au,g] = get_g(U,Ux,Uy,Uf,Wx,Wy,Wf,lam1,...
     lam6,beta,mu,upsilon,A,b,Atb,sigmax,sigmay,sigmaf,delta,Ohm);
   
    % Nonmonotone Line Search
    alpha = 1;
    du = U - Up;                          % du: p*q
    const = opts.c*beta*(d'*taud);

    % Unew = Up + alpha*(U - Up)
    cnt = 0; flag = true;
    while f > C - alpha*const
        if cnt == 5
            % shrink gam
            gam = opts.rate_gam*gam;

            % give up and take Steepest Descent step
            if opts.disp
                disp('    count of back tracking attains 5 ');
            end

            %d = g2p + muDbeta*gp - DtsAtd;
            [dx,dy] = D(reshape(d,p,q));
            df = fft2c(reshape(d,p,q));
            dDd = norm(dx,'fro')^2 + norm(dy,'fro')^2 + Ohm*upDbeta*norm(df,'fro')^2;      %5/12 REVISE
            Ad = A(d,1);
            tau = abs((d'*d)/(dDd + muDbeta*Ad'*Ad));
            U = Up(:) - tau*d;
            % projected gradient method for nonnegtivity
            if opts.nonneg
                U = max(real(U),0);
            end
            U = reshape(U,p,q);
            [Ux,Uy] = D(U);
            Uf = fft2c(U);
            Uxbar = Ux - sigmax/beta;
            Uybar = Uy - sigmay/beta; 
            Ufbar = Uf - sigmaf/upsilon;
            if opts.TVnorm == 1
                % ONE-DIMENSIONAL SHRINKAGE STEP
                Wx = max(abs(Uxbar) - 1/beta, 0).*sign(Uxbar);
                Wy = max(abs(Uybar) - 1/beta, 0).*sign(Uybar);
                Wf = max(abs(Ufbar) - 1/upsilon, 0).*sign(Ufbar);
                lam1 = sum(sum(abs(Wx) + abs(Wy)));
                lam6 = sum(sum(abs(Wf)));
            else
                % TWO-DIMENSIONAL SHRINKAGE STEP
                V = sqrt(Uxbar.*conj(Uxbar) + Uybar.*conj(Uybar)); % V: p*q
                V(V==0) = 1;
                S = max(V - 1/beta, 0)./V;                         % S: p*q
                Wx = S.*Uxbar;
                Wy = S.*Uybar;
                lam1 = sum(sum(sqrt(Wx.*conj(Wx) + Wy.*conj(Wy))));
            end
%             [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Wx,Wy,lam1,...
%                 beta,mu,A,b,Atb,sigmax,sigmay,delta);
           [lam2,lam3,lam4,lam5,lam7,lam8,f,g2,g3,Au,g] = get_g(U,Ux,Uy,Uf,Wx,Wy,Wf,lam1,...
               lam6,beta,mu,upsilon,A,b,Atb,sigmax,sigmay,sigmaf,delta,Ohm);
            alpha = 0; % remark the failure of back tracking
            break;
        end
        if flag
            dg = g - gp;
            dg2 = g2 - g2p;
            dg3 = g3 - g3p;
            dAu = Au - Aup;                 % dAu: m
            dUx = Ux - Uxp;
            dUy = Uy - Uyp;
            dUf = Uf - Ufp;
            flag = false;
        end
        alpha = alpha*opts.gamma;
        [U,lam2,lam3,lam4,lam5,lam7,lam8,f,Ux,Uy,Uf,Au,g,g2,g3] = update_g(p,q,lam1,lam6,...
           alpha,beta,mu,upsilon,Up,du,gp,dg,g2p,dg2,g3p,dg3,Aup,dAu,Wx,Wy,Wf,Uxp,dUx,Uyp,dUy,Ufp,dUf,b,...
           sigmax,sigmay,sigmaf,delta,Ohm);
        cnt = cnt + 1;
    end
    
    % if back tracking is succeceful, then recompute
    if alpha ~= 0
        Uxbar = Ux - sigmax/beta;
        Uybar = Uy - sigmay/beta;
        Ufbar = Uf - sigmaf/upsilon;
        if opts.TVnorm == 1
            % ONE-DIMENSIONAL SHRINKAGE STEP
            Wx = max(abs(Uxbar) - 1/beta, 0).*sign(Uxbar);
            Wy = max(abs(Uybar) - 1/beta, 0).*sign(Uybar);
            Wf = max(abs(Ufbar) - 1/upsilon, 0).*sign(Ufbar);
        else
            % TWO-DIMENSIONAL SHRINKAGE STEP
            V = sqrt(Uxbar.*conj(Uxbar) + Uybar.*conj(Uybar));
            V(V==0) = 1;
            S = max(V - 1/beta, 0)./V;
            Wx = S.*Uxbar;
            Wy = S.*Uybar;
        end
        
        % update parameters related to Wx, Wy;
        [lam1,lam2,lam4,lam6,lam7,lam8,f,g2,g3] = update_W(beta,upsilon,...
            Wx,Wy,Wf,Ux,Uy,Uf,sigmax,sigmay,sigmaf,lam1,lam2,lam4,lam6,lam7,lam8,f,opts.TVnorm,Ohm);
    end
    
    % update reference value
    Qp = Q; Q = gam*Qp + 1; C = (gam*Qp*C + f)/Q;
    uup = U - Up; uup = uup(:);           % uup: pq
    nrmuup = norm(uup,'fro');                   % nrmuup: constant
    
    out.res = [out.res; nrmuup];
    out.f = [out.f; f]; out.C = [out.C; C]; out.cnt = [out.cnt;cnt];
    out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
    out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];out.lam6 = [out.lam6; lam6];
    out.lam7 = [out.lam7; lam7]; out.lam8 = [out.lam8; lam8];
    out.tau = [out.tau; tau]; out.alpha = [out.alpha; alpha];

    if ~isempty(Ut), out.n2re = [out.n2re; norm(U - Ut,'fro')/norm(Ut,'fro')]; end

    % compute gradient
%     d = g2 + muDbeta*g - DtsAtd;
    d = g2 + muDbeta*g - DtsAtd + Ohm*upDbeta*g3;

    nrmup = norm(Up,'fro');
    RelChg = nrmuup/nrmup;

    if RelChg < tol_inn && ~fst_itr
        count = count + 1;
        RelChgOut = norm(U-rcdU,'fro')/nrmup;
        out.reer = [out.reer; RelChgOut];
        rcdU = U;
        out.obj = [out.obj; f + lam4 + lam5];
        if isempty(out.itrs)
            out.itrs = ii;
        else
            out.itrs = [out.itrs; ii - sum(out.itrs)];
        end

        % stop if already reached final multipliers
        if RelChgOut < tol_out || count > opts.maxcnt 
            if opts.isreal
                U = real(U);
            end
            if exist('scl','var')
                U = U/scl;
            end
            out.itr = ii;
            fprintf('Number of total iterations is %d. \n',out.itr);
            return
        end
        
        % update multipliers
%         [sigmax,sigmay,delta,lam4,lam5,f] = update_mlp(beta,mu, ...
%             Wx,Wy,Ux,Uy,Au,b,sigmax,sigmay,delta,lam4,lam5,f);
        [sigmax,sigmay,sigmaf,delta,lam4,lam5,lam7,f] = update_mlp(beta,mu,upsilon, ...
         Wx,Wy,Wf,Ux,Uy,Uf,Au,b,sigmax,sigmay,sigmaf,delta,lam4,lam5,lam7,f,Ohm);
 
        % update penality parameters for continuation scheme
        beta0 = beta;
        beta = beta*opts.rate_ctn;
        mu = mu*opts.rate_ctn;
        upsilon = upsilon*opts.rate_ctn;
        if beta > betaf; beta = betaf; end
        if mu > muf; mu = muf; end
        if upsilon > upsilonf; upsilon = upsilonf; end 
        muDbeta = mu/beta;
        upDbeta = upsilon/beta;
        
        out.mus = [out.mus; mu]; out.betas = [out.betas; beta]; out.upsilons = [out.upsilons; upsilon];

        % update function value, gradient, and relavent constant
        f = lam1 + lam6 + beta/2*lam2 + mu/2*lam3 + upsilon/2*lam8 - lam4 - lam5 - lam7;
        DtsAtd = -(beta0/beta)*d;     % DtsAtd should be divded by new beta instead of the old one for consistency!  
        d = g2 + muDbeta*g + Ohm*upDbeta*g3 - DtsAtd;
        
        %initialize the constants
        gp = [];
        gam = opts.gam; Q = 1; C = f;
    end

end

if opts.isreal
    U = real(U);
end
if exist('scl','var')
    fprintf('Attain the maximum of iterations %d. \n',opts.maxit);
    U = U/scl;
end




function [lam2,lam3,lam4,lam5,lam7,lam8,f,g2,g3,Au,g] = get_g(U,Ux,Uy,Uf,Wx,Wy,Wf,lam1,...
    lam6,beta,mu,upsilon,A,b,Atb,sigmax,sigmay,sigmaf,delta,Ohm)
global Dt

% A*u 
Au = A(U(:),1);

% g   
g = A(Au,2) - Atb;



% lam2
Vx = Ux - Wx;
Vy = Uy - Wy;
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));


% g2 = D'(Du-w)
g2 = Dt(Vx,Vy);

% lam3
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;

%lam4
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));

%lam5
lam5 = delta'*Aub;

%lam7
Vf = Uf - Wf;
lam7 = Ohm*sum(sum(conj(sigmaf).*Vf));

%lam8
lam8 = Ohm*sum(sum(Vf.*conj(Vf)));

% g3 = F'(Fu-Wf)
g3 = ifft2c(Vf);
g3 = g3(:);

% f
f = lam1 + lam6 + beta/2*lam2 + mu/2*lam3 + upsilon/2*lam8 - lam4 - lam5 - lam7;



function [U,lam2,lam3,lam4,lam5,lam7,lam8,f,Ux,Uy,Uf,Au,g,g2,g3] = update_g(p,q,lam1,lam6,...
    alpha,beta,mu,upsilon,Up,du,gp,dg,g2p,dg2,g3p,dg3,Aup,dAu,Wx,Wy,Wf,Uxp,dUx,Uyp,dUy,Ufp,dUf,b,...
    sigmax,sigmay,sigmaf,delta,Ohm)

g = gp + alpha*dg;
g2 = g2p + alpha*dg2;
g3 = g3p + alpha*dg3;
U = Up + alpha*reshape(du,p,q);
Au = Aup + alpha*dAu;
Ux = Uxp + alpha*dUx;
Uy = Uyp + alpha*dUy;
Uf = Ufp + alpha*dUf;

Vx = Ux - Wx;
Vy = Uy - Wy;
Vf = Uf - Wf;
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
lam5 = delta'*Aub;
lam7 = Ohm*sum(sum(conj(sigmaf).*Vf));
lam8 = Ohm*sum(sum(Vf.*conj(Vf)));

f = lam1 + lam6 + beta/2*lam2 + mu/2*lam3 + upsilon/2*lam8 - lam4 - lam5 - lam7;




function [lam1,lam2,lam4,lam6,lam7,lam8,f,g2,g3] = update_W(beta,upsilon,...
    Wx,Wy,Wf,Ux,Uy,Uf,sigmax,sigmay,sigmaf,lam1,lam2,lam4,lam6,lam7,lam8,f,option,Ohm)
global Dt

% update parameters because Wx, Wy were updated
tmpf = f -lam1 - beta/2*lam2 + lam4 -lam6 + lam7 - upsilon/2*lam8;
if option == 1
    lam1 = sum(sum(abs(Wx) + abs(Wy)));
else
    lam1 = sum(sum(sqrt(Wx.^2 + Wy.^2)));
end
lam6 = sum(sum(abs(Wf)));
Vx = Ux - Wx;
Vy = Uy - Wy;
Vf = Uf - Wf;
g2 = Dt(Vx,Vy);
g3 = ifft2c(Vf);
g3 = g3(:);
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
lam7 = Ohm*sum(sum(conj(sigmaf).*Vf));
lam8 = Ohm*sum(sum(Vf.*conj(Vf)));
f = tmpf +lam1 + beta/2*lam2 - lam4 + lam6 - lam7 + upsilon/2*lam8;



function [sigmax,sigmay,sigmaf,delta,lam4,lam5,lam7,f] = update_mlp(beta,mu,upsilon, ...
    Wx,Wy,Wf,Ux,Uy,Uf,Au,b,sigmax,sigmay,sigmaf,delta,lam4,lam5,lam7,f,Ohm)

Vx = Ux - Wx;
Vy = Uy - Wy;
Vf = Uf - Wf;
sigmax = sigmax - beta*Vx;
sigmay = sigmay - beta*Vy;
sigmaf = sigmaf - upsilon*Vf;

Aub = Au-b;
delta = delta - mu*Aub;

tmpf = f + lam4 + lam5 + lam7;
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
lam5 = delta'*Aub;
lam7 = Ohm*sum(sum(conj(sigmaf).*Vf));
f = tmpf - lam4 - lam5 - lam7;




function [U,mu,beta,upsilon] = ftvcs_al_init(p,q,Atb,scl,opts)

% initialize mu beta
if isfield(opts,'mu0')
    mu = opts.mu0;
else
    error('Initial mu is not provided.');
end
if isfield(opts,'beta0')
    beta = opts.beta0;
else
    error('Initial beta is not provided.');
end
if isfield(opts,'upsilon0')
    upsilon = opts.upsilon0;
else
    error('Initial upsilon is not provided.');
end

% initialize U
[mm,nn] = size(opts.init);
if max(mm,nn) == 1
    switch opts.init
        case 0, U = zeros(p,q);
        case 1, U = reshape(Atb,p,q);    %????????
    end
else
    U = opts.init*scl;  
    if mm ~= p || nn ~= q
        fprintf('Input initial guess has incompatible size! Switch to the default initial guess. \n');
        U = reshape(Atb,p,q);
    end
end

