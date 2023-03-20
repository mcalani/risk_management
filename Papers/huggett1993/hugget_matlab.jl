% Model Parameters: Page 84 %%%%%%%%%
e   = [1,.1];   % Endowments        %
sig = 1.5;      % Sigma             %
b   = .99322;   % Beta              %
g   = 200;      % Grid Size         %
x   = [.925,.075;.5,.5];            %
tol = 1e-5;                         %
one = ones(g,1);                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_grid  = linspace(-2,3.97,g);
qhi     = 1.2;
qlo     = b/2;
conv    = 100;
t       = 1;

while (conv > tol) & (t < 30)
    q   = (qhi + qlo) / 2;
    for i = 1:g
        A   = a_grid(1,i);
        for j = 1:g
            A_next  = a_grid(1,j);
            c1(i,j) = A + e(1,1) - q*A_next;
            c2(i,j) = A + e(1,2) - q*A_next;
            u1(i,j) = (c1(i,j)^(1-sig)) / (1-sig);
            u2(i,j) = (c2(i,j)^(1-sig)) / (1-sig);
            if c1(i,j) <= 0
                u1(i,j) = -inf;
            end
            if c2(i,j) <= 0
                u2(i,j) = -inf;
            end
        end
    end
    v1 = zeros(1,g); v2 = v1;
    pos     = zeros(g,2);
    conv1   = 100;
    conv2   = 100;
    while (conv1 > .001) | (conv2 ~= 0)
        vala = (u1 + b*one*( x(1,1)*v1 + x(1,2)*v2 ) );
        valb = (u2 + b*one*( x(2,1)*v1 + x(2,2)*v2 ) );
        [V1,pos1]  = max( vala' );
        [V2,pos2]  = max( valb' );
        V           = [V1,V2];
        v           = [v1,v2];
        Pos         = [pos1',pos2'];
        conv1       = abs(V - v);
        conv1       = max( (max(conv1))' );
        conv2       = max(any(Pos-pos));
        v1 = V1;    v2 = V2;
        pos         = Pos;
    end
    posa = zeros(g,g); posb = posa;
    for i = 1:g
        posa(i,pos1(1,i)) = 1; posb(i,pos2(1,i)) = 1;
    end
    t_matx  = [x(1,1)*posa,x(1,2)*posa;
               x(2,1)*posb,x(2,2)*posb];
    dist = (ones(2*g,1))/(2*g);
    conv1   = 100;
    while conv1 > tol
        dist_p     = t_matx'*dist;
        conv1   = max(abs(dist_p-dist));
        dist       = dist_p;
    end
    ExDD    = [a_grid,a_grid]*dist;
    ExcessDemand(t) = ExDD;
    Q(t)    = q;
    if (t>1) & (abs(Q(t) - Q(t-1))<tol)
        break
    end

    if ExDD < 0
        qhi = (qhi + qlo) / 2;
    else
        qlo = (qhi + qlo) / 2;
    end
    Qhi(t)  = qhi;
    Qlo(t)  = qlo;
    t       = 1 + t;
end


pol_1       = zeros(1,g);
pol_2       = pol_1;

for i = 1:g
    pol_1(1,i) = a_grid(1,pos(i,1));
    pol_2(1,i) = a_grid(1,pos(i,2));
end
