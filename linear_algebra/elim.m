function [Amod,ord]=elim(A,b)

% Amod=elim(A,b)
%
% This function perform elimination with partial pivoting and scaling as
% described in Section 1.3.2 in the Hoffman textbook (viz. it does Gaussian
% elimination).  Note that the ordering which preserves upper triangularity
% is stored in the ord output variable.  

%Allocation and setup
Amod=cat(2,A,b);          %make a copy of A and modify with RHS of system
n=size(A,1);              %number of unknowns
ord=[1:n]';               %ord is a mapping from input row being operated upon to the actual row that represents in the matrix ordering

%Elimination with scaled, partial pivoting for matrix Amod; note all row
%indices must be screen through ord mapping.
for ir1=1:n-1
    %check scaled pivot elements to see if reordering should be done
    pivmax=0;
    ipivmax=ir1;      %max pivot element should never be higher than my current position
    for ipiv=ir1:n    %look only below my current position in the matrix
      pivcurr=Amod(ord(ipiv),ir1)/max(Amod(ord(ir1),:));      %note that columns never get reordered...
      if (pivcurr>pivmax) 
          pivmax=pivcurr;
          ipivmax=ipiv;     %this stores the index into ord for row having largest pivot element
      end %if
    end %for
    
    %reorder if situation calls for it
    if (ipivmax ~= ir1)
        itmp=ord(ir1);
        ord(ir1)=ord(ipivmax);
        ord(ipivmax)=itmp;
    end %if
        
    %perform the elimination for this row, former references to ir1 are now
    %mapped through the ord array
    for ir2=ir1+1:n
        fact=Amod(ord(ir2),ir1);
        Amod(ord(ir2),:)=Amod(ord(ir2),:)-fact/Amod(ord(ir1),ir1).*Amod(ord(ir1),:);
    end %for
end %for

end %function
