clc; clear all; close all;
format long;

%% Iterating over values of the angle of attack
for alpha=-10:1:30
    Vinf       = 30;
    alfa       = (alpha)*(pi/180);
    
    %% Loading points from data file
    data       = load('Data.txt');       
    X1          = data(:,1);
    Y1          = data(:,2);
    X = X1/1000.0;
    Y = Y1/1000.0;
    Y(1,1)     = 0;
    Y(end,1)   = 0;
    n          = size(X,1)-1;
    for index=1:n
        %angle of flow with tangent of panel
        phi(index) = -alfa+...
            atan2((Y(index+1)-Y(index)),(X(index+1)-X(index)));
        
        %% Calculating the angle of flow with normal of the panel,
        % the midpoints (control points) and the length of the panels
        beta(index) = phi(index)+pi/2;
        midpoint_x(index) = (X(index+1)+X(index))/2;
        midpoint_y(index) = (Y(index+1)+Y(index))/2;
        S(index) = sqrt((Y(index+1)-Y(index))^2+...
            (X(index+1)-X(index))^2);
    end

    %The Source Panel Method
    for p=1:n
        neighbors(:,p)=[1:p-1 p+1:n];
        xi=midpoint_x(p);
        yi=midpoint_y(p);
        for index=1:n-1
            m      = neighbors(index,p);
            Xj     = X(m);
            Yj     = Y(m);
            Xj1    = X(m+1);
            Yj1    = Y(m+1);
            A      = -(xi-Xj)*cos(phi(m))-(yi-Yj)*sin(phi(m));
            B      = (xi-Xj)^2+(yi-Yj)^2;
            C      = sin(phi(p)-phi(m));
            D      = (yi-Yj)*cos(phi(p))-(xi-Xj)*sin(phi(p));
            E      = sqrt(B-A^2);
            Sj     = S(m);
            I(p,m) = C/2*log((Sj^2+2*A*Sj+B)/B)+...
                       (D-A*C)/E*(atan2((Sj+A),E)-atan2(A,E));
            J(p,m) = (D-A*C)/2/E*log((Sj^2+2*A*Sj+B)/B)...
                       -C*(atan2((Sj+A),E)-atan2(A,E));
        end
        F(p,1)=Vinf*cos(beta(p));
    end

    M=I/2/pi+eye(n)/2;
    lambda=-inv(M)*F;
    V        = Vinf*sin(beta)+lambda'/2/pi*J';
    Cp       = 1-(V/Vinf).^2;
    angles   = min(beta):0.01:max(beta);
    Cp_exact = 1-4*sin(angles).^2;

    %% Calculating the circulation
    circulation = dot(V,S);
    fprintf('%d -> %.5f\n',alpha, circulation);  
end
