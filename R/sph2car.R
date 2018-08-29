sph2car = function(long, lat, radius=1, deg=TRUE){
    if (is.matrix(long) || is.data.frame(long)) {
        if (ncol(long) == 1) {
            long = long[, 1]
        }
        else if (ncol(long) == 2) {
            lat = long[, 2]
            long = long[, 1]
        }
        else if (ncol(long) == 3) {
            radius = long[, 3]
            lat = long[, 2]
            long = long[, 1]
        }
    }
    if (missing(long) | missing(lat)) {
        stop("Missing full spherical 3D input data.")
    }
    if (deg) {
        long = long * pi/180
        lat = lat * pi/180
    }
    return = cbind(x = radius * cos(long) * cos(lat), y = radius * 
        sin(long) * cos(lat), z = radius * sin(lat))
}

car2sph = function (x, y, z, deg = TRUE) 
{
    if (is.matrix(x) || is.data.frame(x)) {
        if (ncol(x) == 1) {
            x = x[, 1]
        }
        else if (ncol(x) == 2) {
            y = x[, 2]
            x = x[, 1]
        }
        else if (ncol(x) == 3) {
            z = x[, 3]
            y = x[, 2]
            x = x[, 1]
        }
    }
    if (missing(x) | missing(y) | missing(z)) {
        stop("Missing full cartesian 3D input data.")
    }
    radius = sqrt(x^2 + y^2 + z^2)
    long = atan2(y, x)
    lat = asin(z/radius)
    if (deg) {
        long = long * 180/pi
        lat = lat * 180/pi
    }
    lat[radius == 0] = 0
    return = cbind(long = long, lat = lat, radius = radius)
}