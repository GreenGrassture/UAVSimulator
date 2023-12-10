# Models the process that the microcontroller uses to get sensor readings

function gyroscope(x, p, t)
    return @views x[11:13]
end

function accelerometer(x, p, t)
    return @views x[]
end