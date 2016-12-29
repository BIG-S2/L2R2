% this program partitions the given matrix A at the given position k and
% returns the k-th column after removing the k-th element

function [aj] = SliceColumn(A, k)
  d = size(A);
  if (d(1) ~= d(2)); 
      error('Error in function mypartition: need a square matrix');
  elseif (k<1 || d(1)<k) ;
      error('Error in function mypartition: partition location out of range');
  else 
     %aj = A(:,k);
     %aj(k,:) = []; % deleting k-th row(/element)
     
     %change for coder
     idx = zeros(size(A,2)-1,1);
     if k==1
         idx(:,1) = 2:size(A,2);
     elseif k==size(A,2)
         idx(:,1) = 1:(size(A,2)-1);
     else
         idx(:,1) = [1:(k-1) (k+1):(size(A,2))];
             end
     
     aj = A(idx,k);
   
             
  end
end
