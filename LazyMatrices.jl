module LazyMatrices
    export LazyMatrix, getaxes
    """
        LazyMatrix{T} defines a matrix that is lazily initialized.

    - The initialization function `init` is called only when certain matrix elements are actually accessed. 
    - This helps when the initialization is expensive and only a small portion of the matrix needs to be used.

    > The motivation comes from the following problem: given a function `f(x,y)`, we want to find the roots of `f(x,y) = 0` on a 2D grid. 
    The function `f(x,y)` is expensive to evaluate(e.g. includes PDE), and all we need are just roots(1D spectrum on 2D grid) which only exist in a small portion of the grid.
    Some algorithms for finding roots can dynamically adjust steps to stay close to the roots(spectrum), so we don't need to evaluate `f(x,y)` on the whole grid.
    Therefore, the `LazyMatrix` is designed to reduce unnecessary function evaluation costs.
    """
    struct LazyMatrix{T} <: AbstractMatrix{T}
        data::Matrix{T}
        init::Function
        # a flag matrix of booleans to indicate whether the value at the corresponding position has been initialized
        initialized::Matrix{Bool}
        # the axis(metric) functions: given an index, return the corresponding value on the axis
        xaxis::Function
        yaxis::Function

        
    end

    """
        LazyMatrix(data, init, xrange, yrange; logscaled = true) creates a `LazyMatrix` with the given data, initialization function, and axis ranges.

    # Arguments
    - `data`: the matrix to store the initialized values later. 
    Its sizes decide **resolution of the axes**.
    The initial values are not used.
    - `init`: the initialization function for the matrix elements. It takes two parameter values(**not the matrix indices**), and returns the evaluated answer corresponding to the given parameters.
    - `xrange`: the range of the x axis
    - `yrange`: the range of the y axis
    - `logscaled`: whether the `x` and `y` axes are log scaled. Default is `true`.

    # Examples
    ```
    using .LazyMatrices
    m = LazyMatrix(zeros(10,10), (a,b)->a+b, (-3,3), (-3,3), logscaled = true)
    ```
    - by default, the x and y axis ranges are taken as log scaled value, e.g. `xrange = (-3, 3)` means the x axis ranges from 1e-3 to 1e3)
    - if `logscaled = false`, the x and y axis ranges are taken as linearly scaled value, e.g. `xrange = (-3, 3)` means the x axis ranges from -3 to 3)
    - given the indices `(i,j)`, `a=xaxis(i)` and `b=yaxis(j)`, then `init = (a,b)->a+b` means the matrix element at `(i,j)` will be initialized as `a+b`.
    """
    function LazyMatrix(data::Matrix{T}, init::Function, xrange, yrange; logscaled = true) where T
        xmin, xmax = min(xrange...), max(xrange...)
        ymin, ymax = min(yrange...), max(yrange...)
        # specify the axis(metric) functions for both linear and log scaled axes
        xaxis4linscaled(i) = (i-1)/(size(data, 1)-1)*(xmax-xmin)+xmin
        yaxis4linscaled(j) = (j-1)/(size(data, 2)-1)*(ymax-ymin)+ymin
        xaxis4logscaled(i) = exp10((i-1)/(size(data, 1)-1)*(xmax-xmin)+xmin)
        yaxis4logscaled(j) = exp10((j-1)/(size(data, 2)-1)*(ymax-ymin)+ymin)
        # at the beginning, all values are uninitialized
        initialized = falses(size(data))
        if logscaled
            return new{T}(data, init, initialized, xaxis4logscaled, yaxis4logscaled)
        else
            return new{T}(data, init, initialized, xaxis4linscaled, yaxis4linscaled)
        end
    end
    function Base.getindex(m::LazyMatrix, i::Int, j::Int)
        if !m.initialized[i,j]
            m.data[i,j] = m.init(m.xaxis(i), m.yaxis(j))
            m.initialized[i,j] = true
        end
        m.data[i,j]
    end
    function Base.getindex(m::LazyMatrix, i::Real, j::Real)
        return m.init(m.xaxis(i), m.yaxis(j))
    end
    function Base.setindex!(m::LazyMatrix, v, i, j)
        m.data[i,j] = v
        if !m.initialized[i,j]
            m.initialized = true
        end
    end
    function Base.size(m::LazyMatrix)
        size(m.data)
    end
    function Base.IndexStyle(::Type{<:LazyMatrix})
        IndexCartesian()
    end
    """
        getaxes(m::LazyMatrix) returns the x and y axes(as vectors) of the matrix `m`.
    - The resolution of the axes is determined by the size of the matrix `m`.
    - The start and end values of the axes are determined by the initialization ranges of the matrix `m`.
    """
    function getaxes(m::LazyMatrix)
        return (m.xaxis.(range(1,size(m.data,1))), m.yaxis.(range(1,size(m.data,2))))
    end
end