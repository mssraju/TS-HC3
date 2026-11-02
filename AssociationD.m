function [Population,z,count] = AssociationD(R,W,z,wb,Global)
    if Global.FE/Global.maxFE>0.8
           [F,M] = CNDSort(R.objs,R.cons,Global.N);
           R     = R(F<=M);
           zr=min(R.objs,[],1);
  
    end
    z=min(R.objs,[],1);
    if Global.FE/Global.maxFE>0.8
        z=zr;
    end
    [PopObj,z] = Normalization5(R.objs,z);
     normP  = sqrt(sum(PopObj.^2,2));
     Cosine = 1 - pdist2(PopObj,W,'cosine');
     d1     = repmat(normP,1,size(W,1)).*Cosine;
     d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
     sr=[];
     count=0;
    %% Clustering
    [~,class] = min(d2,[],2);
    asso.sol  =[];
    for i = 1 : length(W)
        asso(i).sol = find(class==i);   
        if isempty(asso(i).sol)==1
            [~,id]=min(d2(:,i));
            Population(i)=R(id);
            if id>Global.N
                  count=count+1;
            end
            sr=[sr,id];
        else
            [F,M] = CNDSort(R([asso(i).sol]).objs,R([asso(i).sol]).cons,length(R([asso(i).sol])));
            if length(find(F==1))>1
                a=find(F==1);
                pbi=d1(a,i)+5*d2(a,i);
                [~,ind]=min(pbi);
                Population(i)=R(asso(i).sol(a(ind)));
                if (asso(i).sol(a(ind)))>Global.N
                    count=count+1;
                end 
                sr=[sr,[asso(i).sol(a(ind))]];
            else
               Population(i)=R([asso(i).sol(find(F==1))]);
               if asso(i).sol(find(F==1))>Global.N
                   count=count+1;
               end 
               sr=[sr,[asso(i).sol(find(F==1))]];
            end    
        end   
    end
    
end