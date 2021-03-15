N = 20;

u = solver(N);

fplot(u,[0,2]);
ylim([0,40]);


function u = solver(N)
    %początek i koniec przedziału;
    START = 0;
    END = 2;
    
    syms x;
    E(x) = piecewise(and(x>=0, x<=1), 3, and(x>1, x<=2), 5);
    
    % parametry równania
    a = E;
    b = 0;
    c = 0;
    
    % długość przedziału
    len = (END - START)/N;
    
    % macierz
    matrix = zeros(N,N);
    
    basis = basis_f(N, len);
    
    for i = 1:N
        for j = i:N
            u(x) = basis{j};
            v(x) = basis{i};
            %main diagonal
            if(and(i==1,j==1))
                matrix(i,j) = E(0)*u(0)*v(0) - my_quad(E*diff(u)*diff(v),0,i*len);
            elseif(and(i==j,i~=1))
                matrix(i,j) = E(0)*u(0)*v(0) - my_quad(E*diff(u)*diff(v),(i-2)*len, i*len);
            elseif(j-i==1)
                matrix(i,j) = E(0)*u(0)*v(0) - my_quad(E*diff(u)*diff(v),(i-1)*len, i*len);
            end
            matrix(j,i) = matrix(i,j);
         end
    end
    
    matrix
    
    L = zeros(N,1);
    L(1) = 30;
    X = linsolve(matrix,L);
    
    
    u = 0;
    for m = 1:N
        u = u + basis{m}*X(m);
    end
   
end
    
% funkcja zwracająca tablicę funkcji bazowych
function basis = basis_f(N,len)
    syms x
    for k = 1:N+1
        if(k==1)
            e(x) = piecewise(and(x>=0, x<=len), (len-x)/len, 0);
        elseif(k==N+1)
            e(x) = piecewise(and(x>=(N*len-len), x<=N*len), (x-(N*len-len))/len, 0);
        else
            e(x) = piecewise(and(x>=len*(k-2), x<=len*(k-1)), ((x-len*(k-2))/len), and(x>len*(k-1), x<=len*k), (len*k-x)/len, 0);
        end
        
        basis{k} = e(x);
    end
    
end
    
function res = my_quad(f,a,b)
    points = [1/sqrt(3) -1/sqrt(3)];
    weights = [1,1];
    
    func = f(((b-a)/2)*points + (a+b)/2);
    
    res = (b-a)/2 * dot(weights,func);
end


