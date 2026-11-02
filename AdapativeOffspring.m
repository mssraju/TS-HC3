function Off = AdapativeOffspring(TPop,Problem,B)
%      fit2 = Fit(TPop.objs,TPop.cons);
% 
%        MatingPool        = TournamentSelection(2,Problem.N,fit2);
%        Off         = OperatorGA(Problem,TPop(MatingPool));
    Off     = [];    
    for i = 1 : Problem.N
        % Choose the parents
            if rand<0.9
            	P = B(i,randperm(size(B,2)));
            else
                P = randperm(Problem.N);
            end
            ran  =  rand;
            if i == 1 || i == Problem.N 
                 Offspring    = OperatorDE(Problem,TPop(i),TPop(P(1)),TPop(P(2)),{1,0.8,1,20});
                 Off          = [Off,Offspring];
            else
                 if ran<1
                    Offspring = OperatorGAhalf(Problem,TPop(P(1:2)));
                 else   
                  Offspring   = OperatorDE(Problem,TPop(i),TPop(P(1)),TPop(P(2)),{1,0.5,1,20});
                 end 
                 Off          = [Off,Offspring];
            end   
            
    end  
end