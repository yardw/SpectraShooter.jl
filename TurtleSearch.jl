module TurtleSearch
    export Turtles, next!, Action
    @enum Action turn_left = 1 turn_right = -1 go_straight = 0 stop = 2
    @enum Direction ipos ineg jpos jneg

    function usualurge(lforwardsignal, rforwardsignal)
        if lforwardsignal > rforwardsignal
            return turn_left
        elseif lforwardsignal < rforwardsignal
            return turn_right
        else
            return go_straight
        end
    end
    mutable struct Turtles
        lpos::CartesianIndex{2} #left foot position
        rpos::CartesianIndex{2} #right foot position(always keep distance 1 from lpos)
        forwardantena::Function #return the signal of the forward antena, given the turtle and a map
        urge::Function #decide the next action, given the left and right signals of the forward antena
        function Turtles(forwardantena::Function, lpos::CartesianIndex{2}, rpos::CartesianIndex{2} = lpos+CartesianIndex(1, 0),  urge::Function = usualurge)
            new(lpos, rpos, forwardantena, urge)
        end
    end

```
    getcurrentdir

    given a turtle, return its current direction
```
    function getcurrentdir(t::Turtles)
        dpos = t.rpos - t.lpos
        @assert dpos[1] * dpos[2] == 0 "turtle not in the line $(@show dpos t)" #check if the turtle's feet are in the line
        @assert abs(dpos[1]) + abs(dpos[2]) == 1 "turtle's feet is not nearby" #check if the turtle's feet are nearby
        if dpos == CartesianIndex(1, 0)
            return jpos
        elseif dpos == CartesianIndex(-1, 0)
            return jneg
        elseif dpos == CartesianIndex(0, 1)
            return ineg
        elseif dpos == CartesianIndex(0, -1)
            return ipos
        end
    end

```
move

given a position and a direction, return the position of the next step
```
    function move(idx::CartesianIndex, dir::Direction)
        if dir == ipos
            return idx + CartesianIndex(1, 0)
        elseif dir == ineg
            return idx + CartesianIndex(-1, 0)
        elseif dir == jpos
            return idx + CartesianIndex(0, 1)
        elseif dir == jneg
            return idx + CartesianIndex(0, -1)
        end
    end
    function move!(t::Turtles, action::Action)
        dpos = [(t.rpos - t.lpos)[i] for i in 1:2]
        turnL = [0 -1; 1 0]
        # turnR = [0 1; -1 0]
        if action == turn_left
            t.rpos = t.lpos + CartesianIndex(turnL * dpos ...)
            return t.lpos, t.rpos
        elseif action == turn_right
            t.lpos = t.rpos + CartesianIndex(turnL * dpos ...) #equals to CartesianIndex(turnR * (-dpos) ...), note that dpos = rpos - lpos
            return t.lpos, t.rpos
        elseif action == go_straight
            current_dir = getcurrentdir(t)
            t.lpos = move(t.lpos, current_dir)
            t.rpos = move(t.rpos, current_dir)
            return t.lpos, t.rpos
        elseif action == stop
            return nothing
        end
    end
    function move(t::Turtles, action::Action)
        dpos = [(t.rpos - t.lpos)[i] for i in 1:2]
        turnL = [0 -1; 1 0]
        turnR = [0 1; -1 0]
        if action == turn_left
            t.rpos = t.lpos + CartesianIndex(turnL * dpos ...)
            return t.lpos, t.rpos
        elseif action == turn_right
            t.lpos = t.rpos + CartesianIndex(turnR * dpos ...)
            return t.lpos, t.rpos
        elseif action == go_straight
            lpos = move(t.lpos, getcurrentdir(t))
            rpos = move(t.rpos, getcurrentdir(t))
            return lpos, rpos
        elseif action == stop
            return nothing
        end
    end
```
    detectnextdir

given a turtle and a map, return the action of the next step
```
    function decidenextaction(t::Turtles, mapmat::Matrix)
        lforwardpos, rforwardpos = move(t, go_straight)
        #check if the turtle is going to be out of the map
        if (
            lforwardpos[1] < 1 || lforwardpos[1] > size(mapmat)[1] 
            || lforwardpos[2] < 1 || lforwardpos[2] > size(mapmat)[2]
            || rforwardpos[1] < 1 || rforwardpos[1] > size(mapmat)[1]
            || rforwardpos[2] < 1 || rforwardpos[2] > size(mapmat)[2]
        )
        return stop
    end
    # @show lforwardpos rforwardpos t
        #if the turtle is not going to be out of the map, then decide the next step
        lforwardsignal = t.forwardantena(mapmat, lforwardpos, t.lpos)
        rforwardsignal = t.forwardantena(mapmat, rforwardpos, t.rpos)
        return nextaction = t.urge(lforwardsignal, rforwardsignal)
    end
```
next!

update the status of the turtle with given map
```
    function next!(t::Turtles, mapmat::Matrix)
        nextaction = decidenextaction(t, mapmat)
        move!(t, nextaction)
    end       
end