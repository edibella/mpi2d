%%%%% joint histogram of two matrices.
% A=[ 0 1 2 3 4 ; 0 2 1 3 4 ; 0 3 5 3 2; 0 1 3 4 6; 0 0 1 2 3 ];
% B=[ 0 2 1 3 4;  0 1 3 4 6;  0 3 5 3 2;  0 1 2 3 4;0 0 1 3 2];
%%%% Arvind
%%%% Takes two matrices or images of equal sizes and ouputs a joint
%%%% histogram.
function h=jointhist(A,B,N)
[r1 r2]=size(A);
Aratio = N/max(A(:));
Bratio = N/max(B(:));
h=zeros(N,N);
for i=1:r1
    for j=1:r2
        h(ceil(A(i,j)*Aratio),ceil(B(i,j)*Bratio))=h(ceil(A(i,j)*Aratio),ceil(B(i,j)*Bratio))+1;
    end
end
return
% disp(h);
% keyboard