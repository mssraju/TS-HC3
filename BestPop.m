function [Population,ind] = BestPop(Population1,Population2,Population3)
    MPop =[Population1,Population2,Population3];
    CV   = sum(max(0,MPop.cons),2);
    MPop = MPop(find(CV==0));
    PF   = [max(MPop.objs);min(MPop.objs)];
    if isempty(Population1(all(Population1.cons<=0,2)).objs)
       HV1=0;
    else
       [HV1,~]  = HVH(Population1(all(Population1.cons<=0,2)).objs,PF);
    end
    if isempty(Population2(all(Population2.cons<=0,2)).objs)
       HV2=0;
    else   
       [HV2,~]  = HVH(Population2(all(Population2.cons<=0,2)).objs,PF);
    end
    if isempty(Population3(all(Population3.cons<=0,2)).objs)
       HV3=0;
    else   
       [HV3,~]  = HVH(Population3(all(Population3.cons<=0,2)).objs,PF);
    end
   
   
    [~,ind]=max([HV1,HV2,HV3]);
    if ind==1
        Population=Population1;
    elseif ind==2
        Population=Population2;
    elseif ind==3
        Population=Population3;
    end      
end