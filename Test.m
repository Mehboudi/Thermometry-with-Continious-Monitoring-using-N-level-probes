Given a Fermionic bath, we know the FI for a trajectory. (T^2*F/(\gamma \tau)) where \gamma is the SD band (flat) and \tau is measurement time. 
(i) what is the opitmal gap if we take the 2-level structure? (numeric good too)
(ii) what is the overal optimal structure, if we do a proper optimisation (use matlab's very own numerical tools)---this is what we do in this code. see the 3lev_ansatz... for the above
function [opt_fi,x] = Fermionic_optimal_Ham_structure(DD,ub)
%DD=2;%for tests
%epsilon=sym('epsilon',[1,DD-1]);
DDm1=DD-1;

w_vec=@(epsilon)[0,epsilon];
w_mat=@(epsilon)ones(DD,1)*w_vec(epsilon);
%eg w_mat_{12}=epsilon_1-epsilon_2 (not the opposit):
w_mat=@(epsilon)...
    transpose(w_mat(epsilon))-w_mat(epsilon);

p_vec=@(epsilon)exp(-w_vec(epsilon));
p_vec=@(epsilon)...
    p_vec(epsilon)./(sum(p_vec(epsilon)));
%%
FI =@(epsilon) w_mat(epsilon).^2.*exp(2*w_mat(epsilon))...
    ./(exp(w_mat(epsilon))+1).^3;
FI=@(epsilon)...
    FI(epsilon)-diag(diag(FI(epsilon)));%Since the summation is over non-equal indices
FI=@(epsilon) nan2zero(FI(epsilon));%In the Bosonic case, there are devision by infty that we should remove
FI =@(epsilon) sum(FI(epsilon),1);
FI =@(epsilon) p_vec(epsilon)*FI(epsilon)';
%FI_sym(epsilon)=FI;
%FF=matlabFunction(-FI,'Vars',{epsilon})%There's a minus sign, we minimise this
FF=@(epsilon) -FI(epsilon);

% %eps0=rand*ones(1,DD-1)
% eps0=rand(1,DDm1);
% UB=ub*ones(1,DDm1);
% LB=0*UB;
% %FF(pi)%This is a test also
% %% Defien the problem and optimise (global)
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon','x0',eps0,...
%     'objective',FF,'lb',LB,'ub',[]);
% x = run(gs,problem);
% opt_fi=FF(x);
% x=sort(x);
% % FF(x)
% % x=sort(x,'descend')

options = optimoptions('patternsearch','Display','iter',...
    'PlotFcn',@psplotbestf,'MaxFunctionEvaluations',1e5);
eps0=rand(1,DDm1);
UB=5*ones(1,DDm1);
LB=0*UB;
x=patternsearch(FF,eps0,[],[],[],[],LB,UB,options);
opt_fi=FF(x);
x=sort(x);

end

Now, take the 2-level anstz. Optimise over gap. Will you get a better result?
% %epsilon=sym('epsilon',[1,DD-1]);
% w_vec=@(epsilon)[0,epsilon*ones(1,DDm1)]%all degenerate (and a g-state)
% w_mat=@(epsilon)ones(DD,1)*w_vec(epsilon);
% %eg w_mat_{12}=epsilon_1-epsilon_2 (not the opposit):
% w_mat=@(epsilon)...
%     transpose(w_mat(epsilon))-w_mat(epsilon);
% 
% p_vec=@(epsilon)exp(-w_vec(epsilon));
% p_vec=@(epsilon)...
%     p_vec(epsilon)./(sum(p_vec(epsilon)));
% %%
% FI =@(epsilon) w_mat(epsilon).^2.*exp(2*w_mat(epsilon))...
%     ./(exp(w_mat(epsilon))+1).^3;
% FI=@(epsilon)...
%     FI(epsilon)-diag(diag(FI(epsilon)));%Since the summation is over non-equal indices
% FI=@(epsilon) nan2zero(FI(epsilon));%In the Bosonic case, there are devision by infty that we should remove
% FI =@(epsilon) sum(FI(epsilon),1);
% FI =@(epsilon) p_vec(epsilon)*FI(epsilon)';
% %FI_sym(epsilon)=FI;
% %FF=matlabFunction(-FI,'Vars',{epsilon})%There's a minus sign, we minimise this
% FF_deg=@(epsilon) -FI(epsilon)

% eps0_deg=rand
% FF_deg(eps0_deg)
% UB=100;
% gs = GlobalSearch;
% UB=100;
% LB=0*UB;
% problem = createOptimProblem('fmincon','x0',eps0_deg,...
%     'objective',FF_deg,'lb',LB,'ub',[]);
% x = run(gs,problem)
% FF_deg(x)
% x=sort(x)
% FF_deg(x)

% %%%Test, do we match the NJP paper?
% Ferm_2_level_NJP=@(w)w.^2/8*cosh(w)./cosh(w/2).^4
% Ferm_2_level_NJP(pi)

function nonnan=nan2zero(M)
    M(isnan(M))=0;
    nonnan=M;
end
