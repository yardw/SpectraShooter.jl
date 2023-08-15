module DelayedMatrices
    export DelayedMatrix
    # , DelayedMatrixView, DelayedMatrixTranspose, DelayedMatrixDiagonal, DelayedMatrixRow, DelayedMatrixCol, DelayedMatrixRowView, DelayedMatrixColView, DelayedMatrixRowTranspose, DelayedMatrixColTranspose, DelayedMatrixRowDiagonal, DelayedMatrixColDiagonal, DelayedMatrixRowViewTranspose, DelayedMatrixColViewTranspose, DelayedMatrixRowViewDiagonal, DelayedMatrixColViewDiagonal, DelayedMatrixRowTransposeView, DelayedMatrixColTransposeView, DelayedMatrixRowDiagonalView, DelayedMatrixColDiagonalView, DelayedMatrixRowTransposeViewTranspose, DelayedMatrixColTransposeViewTranspose, DelayedMatrixRowDiagonalViewDiagonal, DelayedMatrixColDiagonalViewDiagonal, DelayedMatrixRowTransposeViewDiagonal, DelayedMatrixColTransposeViewDiagonal, DelayedMatrixRowDiagonalViewTranspose, DelayedMatrixColDiagonalViewTranspose, DelayedMatrixRowDiagonalViewDiagonal, DelayedMatrixColDiagonalViewDiagonal
    struct DelayedMatrix{T} <: AbstractMatrix{T}
        data::Matrix{T}
        init::Function
        initialized::Matrix{Bool}

        function DelayedMatrix(data::Matrix{T}, init::Function) where T
            initialized = falses(size(data))
            new{T}(data, init, initialized)
        end
    end
    function Base.getindex(m::DelayedMatrix, i::Int, j::Int)
        if !m.initialized[i,j]
            m.data[i,j] = m.init(i, j)
            m.initialized[i,j] = true
        end
        m.data[i,j]
    end
    function Base.getindex(m::DelayedMatrix, i::Real, j::Real)
        return m.init(i, j)
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
end