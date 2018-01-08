 function S = getargs(defaultS, varglist);
 
 if mod(length(varglist),2) ~=0
     error('Odd number of variable parameters');
 end
 
 S = defaultS;
 i=1;
 while i <= length(varglist)
     if isfield(S, varglist{i})
         % for Matlab R12
         %S = setfield(S, varglist{i}, varglist{i+1});
         
         % for Matlab R13 and above
         S.(varglist{i}) = varglist{i+1};
     else
         warning_wrap('getargs:unknown_param', ...
                 ['Unknown parameter "' varglist{i} '"']);
     end
     i = i+2;
 end