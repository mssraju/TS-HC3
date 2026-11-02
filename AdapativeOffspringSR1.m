function Off = AdapativeOffspringSR1(Population,Global)
    
    Off     = []; 
    if length(Population)>4  
        CV = sum(max(0,Population.cons),2);  
        feas = find(CV==0);
        infeas = find(CV~=0);
        FPop=Population(feas).objs;
        InPop=Population(infeas).objs;
        PopObj=Population.objs;
        if length(FPop)<=4
            FPop=Population.objs;
        end
        if length(Population)>4    
            for i=1:length(Population)
                    Ang       = 1 - pdist2(PopObj(i,:),FPop,'cosine');
                    [~,ind] = sort(Ang);
                    ran  =  rand;
                    if rand<0.9
                        if i == ind(1) 
                            P = [i,ind(2),ind(3)];
                        elseif i == ind(2)
                            P = [i,ind(1),ind(3)];
                        else 
                            P = [i,ind(1),ind(2)];
                        end    
                    else
                        P = randperm(length(Population));
                    end
                    Offspring   = OperatorDE(Global,Population(P(1)),Population(P(2)),Population(P(3)),{1,0.5,1,20});
                    Off         = [Off,Offspring];     
            end  
        end
    end

        
end