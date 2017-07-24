function [mod_frontier, fp] = push(frontier, value, fp)
% push.m- Pushes a value onto a stack
% Author: Alexander Murray
% Input: current stack, value to place on stack, and stack pointer
% Output: modified stack with new value, updated stack pointer

    % Frontier pointer increments
    fp = fp + 1;
    % Place value onto stack 
    frontier(1, fp) = value; 
    mod_frontier = frontier; 
    

end

