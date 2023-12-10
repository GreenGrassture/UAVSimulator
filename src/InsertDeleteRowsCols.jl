#### Functions for inserting and deleting rows and columns in arrays ####

function insertRow(arr, idx)
    M,N = size(arr)
    newArr = zeros(M+1, N)
    newArr[1:idx-1,:] = arr[1:idx-1,:]
    newArr[idx+1:end,:] = arr[idx:end,:]
    return newArr
end

function insertRows(arr, idxVec)
    idxVec = sort(idxVec, rev=false)
    newArr = arr
    for idx in idxVec
        newArr = insertRow(newArr, idx)
    end
    return newArr
end

function insertCol(arr, idx)
    M,N = size(arr)
    newArr = [arr[:,1:idx-1] zeros(M, 1) arr[:,idx:end]]
    return newArr
end

function insertCols(arr, idxVec)
    idxVec = sort(idxVec, rev=false)
    newArr = arr
    for idx in idxVec
        newArr = insertCol(newArr, idx)
    end
    return newArr
end

function deleteRow(arr, idx)
    newArr = [arr[1:idx-1,:]
              arr[idx+1:end,:]]
    return newArr
end

function deleteRows(arr,idxVec)
    idxVec = sort(idxVec, rev=true)
    newArr = arr
    for idx in idxVec
        newArr = deleteRow(newArr, idx)
    end
    return newArr
end

function deleteCol(arr, idx)
    newArr = [arr[:,1:idx-1] arr[:,idx+1:end]]
    return newArr
end

function deleteCols(arr, idxVec)
    idxVec = sort(idxVec, rev=true)
    newArr = arr
    for idx in idxVec
        newArr = deleteRow(newArr, idx)
    end
    return newArr
end


function deleteRowAndCol(arr, idx)
    return deleteCol(deleteRow(arr, idx), idx)
end

function insertRowAndCol(arr, idx)
    return insertCol(insertRow(arr, idx), idx)
end

function deleteRowsAndCols(arr, idxVec)
    # Deletes the rows and columns of an array indicated by idxVec
    idxVec = sort(idxVec, rev=true)
    newArr = arr
    for idx in idxVec
        newArr = deleteRowAndCol(newArr, idx)
    end
    return newArr
end

function insertRowsAndCols(arr, idxVec)
    # Adds a row and column of zeros in the positions indicated by idxVec
    idxVec = sort(idxVec, rev=false)
    newArr = arr
    for idx in idxVec
        newArr = insertRowAndCol(newArr, idx)
    end
    return newArr
end

