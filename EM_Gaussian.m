clc;
close all;
clear;
mu1=[3,3]';                      %intializing mean vector of first distribution
mu2=[1,-3]';                     %intializing mean vector of second distribution
sigma1=[1,0;0,2];                %intializing covariance matrix of first distribution
sigma2=[2,0;0,1];                %intializing covariance matrix of second distribution

N=500;                           %total number of sample points to be generated
pi1=0.8;                         %probability of data point belonging to first distribution
pi2=0.2;                         %probability of data point belonging to second distribution
choose=zeros(1,N);               %intializing a vector
for i=1:N
    choose(i)= binornd(1,0.8);   %selecting first distribution randomly with pi1 probability
end
X=zeros(N,2);                     
i = 0;
iter=0;
while(i<N) 
    i=i+1;
    if choose(i)==1
        X(i,:)=mvnrnd(mu1,sigma1);
    else
        X(i,:)=mvnrnd(mu2,sigma2);
    end
end

plot(X(:,1), X(:,2), '.');
X=X';


mus = rand(2,2);                %randomly generating mean
sigmas = zeros(2,2,2);          
sigmas(:,:,1) = rand*eye(2);    %randomly generating covariance matrix 
sigmas(:,:,2) = rand*eye(2);
pi = [1/2;1/2];

iter = 0;
while(true)

    musold = mus;
    % updating qnk
    [K, N] = size(X);
    qnk = zeros(N, K);

    for i = 1:K
        for j = 1: N

            cons = 0;
            for l = 1: K
               cons =  cons + pi(l) * calc_gauss( X(:,j), mus(:,l), sigmas(:,:,l) );
      
            end

            top = pi(i) * calc_gauss( X(:,j), mus(:,i), sigmas(:,:,i));

            qnk(j,i) = top/cons;
            
            
        end
    end
    
    %updating pi
    N = length(qnk);
    pi = sum(qnk)/N;
    pi = pi';


    [ N, K ] = size(qnk);
  
    % updating means

    for k = 1:K

        res = qnk(:,k)' .* X;
        res = sum(res,2);
        den = sum(qnk(:,k));

        mus(:,k) = res./den;

        
    end


    [N, K] = size(qnk);

    for k = 1:K
        num = 0;
        for n = 1:N
            num = num + qnk(n,k) * ( X(:,n) - mus(:,k) ) * ( X(:,n) - mus(:,k) )' ;
        end
        den = sum( qnk(:,k) );

        sigmas(:,:,k) = num/den;
        
    end

    if( norm( musold - mus ) < 1e-4 )
        break;
    end

    iter = iter + 1;
    disp(['_____ Iteration:  ', num2str(iter), '_______']);
    % disp(qnk)



end


first_count = 0;
second_count = 0;

for i = 1:N
    if( qnk(i,1) > qnk(i,2) )
        first_count = first_count + 1;
        %indices(n) = 1;%corresponds to class 1.
        first(:, first_count) = X(:, i);
    else
        second_count = second_count + 1;
        %indices(n) = 2; %corresponds to class 2.
        second(:, second_count) = X(:, i);
    end
end

figure;
hold on;
scatter(first(1,:), first(2,:))
scatter(second(1,:), second(2,:))
hold off;




function gauss = calc_gauss(x, mus, sigmas)
    
    N = length(x);
    cons = ( 2 * pi * det(sigmas) )^(N/2);
    cons = 1/cons;

    expnum = ( x - mus )' * ( sigmas^(-1) ) * ( x - mus );
    expnum = (-1/2) * expnum;

    gauss = cons * exp(expnum);

end
