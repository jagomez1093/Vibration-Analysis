function [ABCD]=OUTPUT4_rd3(filename);
%
          p0=fopen(filename,'r');
%
          error_m = 0;
%  
          linn=fgets(p0);
          ncr = str2num(linn(1:32));
          AA = repmat(0,ncr(2),ncr(1));
          while (linn~=-1);
              for i=1:ncr(1)+1;
                  linn=fgets(p0);
                  icr = str2num(linn);
                  if (icr(1) > ncr(1));
                      break;
                  else;
                     nfix = fix(icr(3)/3);
                     nmod = mod(icr(3),3);
                     kk=0;
                     for j=1:nfix;
                         linn=fgets(p0);
                         for k=1:3;
                             kk=kk+1;
                             tmp(kk) = str2num(linn(1+(k-1)*23:k*23));
                         end;
                     end;
                     if (nmod ~=0);
                        linn=fgets(p0);
                        for k=1:nmod;
                            kk=kk+1;
                            tmp(kk) = str2num(linn(1+(k-1)*23:k*23));
                        end;
                     end;
%                     
                     if (size(tmp) ~= icr(3));
                         error_m = 1;
                         break;
                     end;
%                         
                     if (ncr(4) == 1 | ncr(4) == 2);
                        AA(icr(2):icr(2)+kk-1,icr(1)) = tmp';
                     else;
                        AA(icr(2):icr(2)+kk/2-1,icr(1)) = complex(tmp(1:2:end-1)',tmp(2:2:end)');
                     end;
%
                     clear tmp
%
                  end; %if
%                  
              end; %for i
%              
              if (icr(1) > ncr(1));
                  break;
              else;
                  error_m = 1;
                  break;
              end;
%                                 
          end; %while
%
          if (error_m == 1);
             dLgname = 'ERROR';
             tdum    = sprintf('*** ERROR in %s file',filename);
             errordlg(tdum,dLgname);
%             break;
              return;
          end;
%          
          ABCD = AA;
%
      fclose(p0);
%
