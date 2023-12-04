module TurtleSearch
    export Turtles, next!, Action, findfirstseed!, bisearch
    """
    an enum type for the **action** of the turtle, which tells the turtle **where** to go next.
    - The action is defined by the relative position of the left foot and the right foot. See also [`usualantena`](@ref) and [`usualurge`](@ref).
    """
    @enum Action turn_left = 1 turn_right = -1 go_straight = 0 stop = 2

    """
        an enum type for the direction of the turtle

    - the direction is defined by the relative position of the left foot and the right foot.
    # e.g., 
        ↑ x x
        j L R  = jpos
        0 i →

        ↑ R L
        j x x  = jneg
        0 i →

        ↑ L x
        j R x  = ipos
        0 i →

        ↑ x R
        j x L  = ineg
        0 i →
    """
    @enum Direction ipos ineg jpos jneg

    """
        usualurge(lforwardsignal, rforwardsignal) decides the next action of the turtle. Return [`Action`](@ref) type.
    # Arguments
    - `lforwardsignal`: the signal of the left forward antena
    - `rforwardsignal`: the signal of the right forward antena 
    """
    function usualurge(lforwardsignal, rforwardsignal)
        if lforwardsignal > rforwardsignal
            return turn_left
        elseif lforwardsignal < rforwardsignal
            return turn_right
        else
            return go_straight
        end
    end
    """
        usualantena(m, newind, oldind) returns the signal of the forward antena, given the turtle and a map
    # Arguments
    - `m`: the map
    - `newind`: the position of the forward antena
    - `oldind`: the position of the turtle
    """
    function usualantena(m::AbstractMatrix, newind::CartesianIndex, oldind::CartesianIndex)
        new = m[newind]
        old = m[oldind]
        @assert new * old != 0 "new * old == 0"
        return new * old < 0
    end
    
    mutable struct Turtles
        lpos::CartesianIndex{2} #left foot position
        rpos::CartesianIndex{2} #right foot position(always keep distance 1 from lpos)
        forwardantena::Function #return the signal of the forward antena, given the turtle and a map
        urge::Function #decide the next action, given the left and right signals of the forward antena
    end
    """
    # Arguments
    - `lpos`: the position of the left foot
    - `rpos`: the position of the right foot
    - `forwardantena`: the function to return the signal of the forward antena, given the turtle and a map
    - `urge`: the function to decide the next action, given the left and right signals of the forward antena

    # Examples
    ```jldoctest
    using .TurtleSearch
    m = [sin(y-x- pi/51) for x in range(0, stop=pi, length=50), y in range(0, stop=pi, length=50)]
    t = Turtles(CartesianIndex(1,2), CartesianIndex(1,1))
    ```
    - by default, the turtle is facing to the positive `x` axis, i.e., the direction of the turtle is `ipos`. See also [`Direction`](@ref).
    """
    function Turtles(lpos::CartesianIndex{2} = CartesianIndex(1,2), rpos::CartesianIndex{2} = lpos+CartesianIndex(0, -1); forwardantena::Function = usualantena, urge::Function = usualurge)# by default, i is y, j is x, dir = ipos
        new(lpos, rpos, forwardantena, urge)
    end
    """
        getcurrentdir
    - given a turtle, return its current direction
    """
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
        # turnR = [0 1; -1 0]
        if action == turn_left
            t.rpos = t.lpos + CartesianIndex(turnL * dpos ...)
            return t.lpos, t.rpos
        elseif action == turn_right
            t.lpos = t.rpos + CartesianIndex(turnL * dpos ...) #equals to CartesianIndex(turnR * (-dpos) ...), note that dpos = rpos - lpos
            return t.lpos, t.rpos
        elseif action == go_straight
            lpos = move(t.lpos, getcurrentdir(t))
            rpos = move(t.rpos, getcurrentdir(t))
            return lpos, rpos
        elseif action == stop
            return nothing
        end
    end
    """
        detectnextdir

    given a turtle and a map, return the action of the next step
    """
    function decidenextaction(t::Turtles, mapmat::AbstractMatrix)
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
    """
        next!, iterate

    update the status of the turtle with given map
    """
    function next!(t::Turtles, mapmat::AbstractMatrix)
        nextaction = decidenextaction(t, mapmat)
        move!(t, nextaction)
    end  
    function bisearch(t::Turtles, m::AbstractMatrix)
        if getcurrentdir(t) in [ipos, ineg]
            var_ind_range = (t.lpos[2], t.rpos[2])
            var_ind_range = (min(var_ind_range...), max(var_ind_range...))
            bisearch_functor_i(x::Real) = m[t.lpos[1]::Int, x::Real]
            return t.lpos[1], fzero(bisearch_functor_i, var_ind_range)
        else
            var_ind_range = (t.lpos[1], t.rpos[1])
            var_ind_range = (min(var_ind_range...), max(var_ind_range...))
            bisearch_functor_j(x::Real) = m[x::Real, t.lpos[2]::Int]
            return fzero(bisearch_functor_j, var_ind_range), t.lpos[2]
        end
    end    
    function fzero(f::Function, var_ind_range::Tuple{Real, Real})
        var_ind_range = (min(var_ind_range...), max(var_ind_range...))
        if f(var_ind_range[1]) * f(var_ind_range[2]) > 0
            error("f(var_ind_range[1]) * f(var_ind_range[2]) > 0")
        end
        while var_ind_range[2] - var_ind_range[1] > 1e-3
            mid = (var_ind_range[1] + var_ind_range[2]) / 2
            if f(mid) * f(var_ind_range[1]) > 0
                var_ind_range = (mid, var_ind_range[2])
            else
                var_ind_range = (var_ind_range[1], mid)
            end
        end
        return var_ind_range[1]
    end 
    function findfirstseed!(m::AbstractMatrix, id1::CartesianIndex{2}, id2::CartesianIndex{2})
        if id1 > id2
            id1, id2 = id2, id1
        end
        did = 2(id2 - id1)
        #check whether id1 and id2 are out of bound
        while (1 <= id1[1] <= size(m, 1) ) && (1 <= id1[2] <= size(m, 2) ) && (1 <= id2[1] <= size(m, 1) ) && (1 <= id2[2] <= size(m, 2) )
            # println(id1, id2)
            if m[id1] * m[id2] > 0
                id1, id2 = id2, id1 + did
                # print(id1, id2)
                continue
            else
                if id1 > id2
                    id1, id2 = id2, id1
                end
                return id1, id2
            end
        end
        return nothing
    end
end