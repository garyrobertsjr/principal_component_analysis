function [fp, pop_value] = pop(frontier, fp)
% pop.m- Pops a value from a stack
% Author: Alexander Murray
% Input: current stack, and stack pointer
% Output: updated stack pointer

    % Acquire value from current stack position
    pop_value = frontier(1, fp); 
    
    % Frontier pointer decrements
    fp = fp - 1;    

end