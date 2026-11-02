function Off1 =  Adapative1(Pop,Problem,B,t)
        if t==1
                 Off1                     =   AdapativeOffspringSR1(Pop,Problem);
  
        else
            Off1=[];
            Pa=Pop;          
            for i = 1 : length(Pa)
                    % Choose the parents
                        if rand<0.9
                            if i>Problem.N
                                P = B(i-Problem.N,randperm(size(B,2)));
                            else
                                P = B(i,randperm(size(B,2)));
                                a=find(P>length(Pa));
                                if length(a)>=1
                                    P(a)=P(a)-length(Pa);
                                end
                            end
                        else
                            P = randperm(length(Pa));
                        end
                        if rand<0.5
                            Offspring = OperatorGAhalf(Problem,Pa(P(1:2)));
                        else   
                          Offspring   = OperatorDE(Problem,Pa(i),Pa(P(1)),Pa(P(2)),{1,0.5,1,20});
                        end
                        Off1            = [Off1,Offspring];   
            end 
        end

    
end