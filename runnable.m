% Pattern formation in Faraday Waves
% Kirsten Meeker
% kmeeker@cs.ucsb.edu
% ORIGINAL FILE http://www.cs.ucsb.edu/~kmeeker/Faraday.m
clear all;
% modifier : @cobwebkanamachi
% environment: gnu octave 3.8.0, mac osx
% if you did not read or know about Faraday Waves, misc docs bellow are.
% http://arxiv.org/pdf/cond-mat/9704060v2.pdf
% above paper is mostly important to run this.
% http://arxiv.org/pdf/1005.5257.pdf
% http://web.archive.org/web/20070806053342/http://haides.caltech.edu/~lifshitz/patterns.html
% http://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/content/toolbox_graph/write_off.m
% http://www.princeton.edu/~wbialek/rome/refs/cross+hohenberg_93.pdf
% input parameters
e = [0.4 0.1 0.001 0.05 0.05 3 1 1 1 1 1];
q = [2*cos(pi/12) 1.2 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12) 2*cos(pi/12)];
a = [0 1 1 1.832 1.8 1.5 3.5 3.3 2 3 4];
c = [5000 200 1 1 1 1 1 1 5000 5000 200]; % linear operator scaling mentioned in Lifshitz paper
p = 10;  % selects which set of input parameters to use
% if suddenly end, please change t,dt 0 to 10,31 to 51, 81 to 101 etc.

% increments
fN = [32 64];   % Number of grid points on a side (square grid)
L = 2*pi;    % System size
fdt = [.05 0.1 0.2 0.5 1.0];

% increments
fN = [32 64];   % Number of grid points on a side (square grid)
L = 2*pi;    % System size
fdt = [.05 0.1 0.2 0.5 1.0];

% loop for different grid resolution
for s=1:1
    N = fN(s);
    dx = L/N;  % Grid spacing for periodic boundary conditions
    
    % initial condition
    %rand('state',0);
    u = 2 * rand(N,N) - 1.0;    % uniform random [-1, 1]
    u = u - mean(mean(u));      % normalize
    
    % compute eigenvalues of linear operator
%http://www.cs.ucsb.edu/~kmeeker/Faraday.pdf (8) lambda!
    k = (0:N-1)*dx;
    for i=1:N
        for j=1:N
            A(i,j) = (-(k(i)^2) - (k(j)^2) + 1)^2 * (-(k(i)^2) - (k(j))^2 + q(p)^2)^2;
        end
    end
    
    % loop for different time steps
    for r = 1:1;
        dt = fdt(r);
        
        % invert linear oeprator
        P = 1./(1 + dt * c(p) * A);
        
        i = 0;
        j = 0;
for t=90:dt:100 % tweak these if fail (31 to 51, 81 to 101, etc)
        %       for t=0:dt:50
            i = i + 1;
            
            % find transform of rhs
            f = u + dt * (e(p)*u + a(p)*u*u - u^3);
            %fd = double ( f );
            fT = fft2(f);   % FFT
            
            uT = fT.*P;
            %uTd = double ( uT );
            u = ifft2(uT);  %inverse FFT
            
            n(i) = norm(u);
            
            % display at intervals
            if mod(i-1,20) == 0
%scatter3(u)
if(0)
                save temp8.txt un;
                save temp7.txt a(p);
                save temp6.txt e(p);
                save temp5.txt dt;
                save temp4.txt p;
                save temp3.txt u;
                save temp2.txt uTd;
                save temp1.txt fTd;
                save temp0.txt f;
end        
                meshc(real(u));
                j = j + 1;
                mov(j) = gcf;
                %Movie(M,3);                
            end
            
        end
%%org Movie here...
%%        fn(r,s) = norm(u);
%%        fm(r,s) = mean(mean(u)); 
%%        fv(r,s) = var(var(u));
%        movie2avi(mov, "out.avi");
    end % r loop
end % s loop
