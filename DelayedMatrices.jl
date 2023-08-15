module DelayedMatrices
    export DelayedMatrix
    # , DelayedMatrixView, DelayedMatrixTranspose, DelayedMatrixDiagonal, DelayedMatrixRow, DelayedMatrixCol, DelayedMatrixRowView, DelayedMatrixColView, DelayedMatrixRowTranspose, DelayedMatrixColTranspose, DelayedMatrixRowDiagonal, DelayedMatrixColDiagonal, DelayedMatrixRowViewTranspose, DelayedMatrixColViewTranspose, DelayedMatrixRowViewDiagonal, DelayedMatrixColViewDiagonal, DelayedMatrixRowTransposeView, DelayedMatrixColTransposeView, DelayedMatrixRowDiagonalView, DelayedMatrixColDiagonalView, DelayedMatrixRowTransposeViewTranspose, DelayedMatrixColTransposeViewTranspose, DelayedMatrixRowDiagonalViewDiagonal, DelayedMatrixColDiagonalViewDiagonal, DelayedMatrixRowTransposeViewDiagonal, DelayedMatrixColTransposeViewDiagonal, DelayedMatrixRowDiagonalViewTranspose, DelayedMatrixColDiagonalViewTranspose, DelayedMatrixRowDiagonalViewDiagonal, DelayedMatrixColDiagonalViewDiagonal
    struct DelayedMatrix{T} <: AbstractMatrix{T}
        data::Matrix{T}
        init::Function
        initialized::Matrix{Bool}
        xaxis::Function
        yaxis::Function

        function DelayedMatrix(data::Matrix{T}, init::Function, xs, ys; logscaled = true) where T
            xmin, xmax = min(xs), max(xs)
            ymin, ymax = min(ys), max(ys)
            xexpinterp(i) = exp10((i-1)/(size(data, 1)-1)*(xmax-xmin)+xmin)
            yexpinterp(j) = exp10((j-1)/(size(data, 2)-1)*(ymax-ymin)+ymin)
            xinterp(i) = (i-1)/(size(data, 1)-1)*(xmax-xmin)+xmin
            yinterp(j) = (j-1)/(size(data, 2)-1)*(ymax-ymin)+ymin
            initialized = falses(size(data))
            if logscaled
                return new{T}(data, init, initialized, xexpinterp, yexpinterp)
            else
                return new{T}(data, init, initialized, xinterp, yinterp)
            end
        end
    end
    function Base.getindex(m::DelayedMatrix, i::Int, j::Int)
        if !m.initialized[i,j]
            m.data[i,j] = m.init(xaxis(i), yaxis(j))
            m.initialized[i,j] = true
        end
        m.data[i,j]
    end
    function Base.getindex(m::DelayedMatrix, i::Real, j::Real)
        return m.init(xaxis(i), yaxis(j))
    end
    function Base.setindex!(m::DelayedMatrix, v, i, j)
        m.data[i,j] = v
        if !m.initialized[i,j]
            m.initialized = true
        end
    end
    function Base.size(m::DelayedMatrix)
        size(m.data)
    end
    function Base.IndexStyle(::Type{<:DelayedMatrix})
        IndexCartesian()
    end
    function getaxes(m::DelayedMatrix)
        return (m.xaxis.(range(1,size(data,1))), m.yaxis.(range(1,size(data,2))))
    end
end