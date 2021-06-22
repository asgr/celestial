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

eq2gal = function(RA, Dec, pole_RA = 192.859508, pole_Dec = 27.128336, eta = 32.932){
    if(is.matrix(RA) || is.data.frame(RA)){
        if(ncol(RA) == 2){
            Dec = RA[, 2]
            RA = RA[, 1]
        }
    }
    
    RA = RA * (pi/180)
    Dec = Dec * (pi/180)
    pole_RA = pole_RA * (pi/180)
    pole_Dec = pole_Dec * (pi/180)
    eta = eta * (pi/180)
    
    gal_lat = asin(cos(Dec) * cos(pole_Dec) * cos((RA - pole_RA)) + sin(Dec) * sin(pole_Dec))
    gal_long = atan2(sin(Dec) - sin(gal_lat) * sin(pole_Dec),
                    cos(Dec) * cos(pole_Dec) * sin((RA - pole_RA))) + eta
    return(cbind(gal_long = gal_long * (180/pi) %% 360, gal_lat = gal_lat * (180/pi)))
}

gal2eq = function(gal_long, gal_lat, pole_RA = 192.859508, pole_Dec = 27.128336, eta = 32.932){
    if(is.matrix(gal_long) || is.data.frame(gal_long)){
        if(ncol(gal_long) == 2){
            gal_lat = gal_long[, 2]
            gal_long = gal_long[, 1]
        }
    }
    
    gal_lat = gal_lat * (pi/180)
    gal_long = gal_long * (pi/180)
    pole_RA = pole_RA * (pi/180)
    pole_Dec = pole_Dec * (pi/180)
    eta = eta * (pi/180)
    
    RA = atan2((cos(gal_lat) * cos(gal_long - eta)), (sin(gal_lat) * cos(pole_Dec) - cos(gal_lat) * sin(pole_Dec) * sin(gal_long - eta))) + pole_RA
    Dec = asin(cos(gal_lat) * cos(pole_Dec) * sin(gal_long - eta) + sin(gal_lat) * sin(pole_Dec))
    return(cbind(RA = RA * (180/pi) %% 360, Dec = Dec * (180/pi)))
}
