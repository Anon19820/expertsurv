
mean_survspline0 <- function(gamma0, gamma1, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline1 <- function(gamma0, gamma1, gamma2, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline2 <- function(gamma0, gamma1, gamma2, gamma3, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline3 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

mean_survspline4 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline5 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline6 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


mean_survspline7 <- function(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp"){
    mean_survspline(gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline)
}




rmst_survspline0 <- function(t, gamma0, gamma1, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline1 <- function(t, gamma0, gamma1, gamma2, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline2 <- function(t, gamma0, gamma1, gamma2, gamma3, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline3 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline4 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}

rmst_survspline5 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline6 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}


rmst_survspline7 <- function(t, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots=c(-10, 10), scale="hazard", timescale="log", spline="rp", start=0){
    rmst_survspline(t, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline, start=start)
}




dsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


dsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


dsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}

dsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}

dsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


dsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


dsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


dsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp", log=FALSE){
    dsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline, log=log)
}


 
psurvspline0 <- function(q, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
psurvspline1 <- function(q, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
psurvspline2 <- function(q, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}


psurvspline3 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

psurvspline4 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}


psurvspline5 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
psurvspline6 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

psurvspline7 <- function(q, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    psurvspline(q, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}





qsurvspline0 <- function(p, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}


qsurvspline1 <- function(p, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}


qsurvspline2 <- function(p, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
qsurvspline3 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
qsurvspline4 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
qsurvspline5 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
qsurvspline6 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}

 
qsurvspline7 <- function(p, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp", lower.tail=TRUE, log.p=FALSE){
    qsurvspline(p, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline, lower.tail=lower.tail, log.p=log.p)
}



 
rsurvspline0 <- function(n, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


rsurvspline1 <- function(n, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
rsurvspline2 <- function(n, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
rsurvspline3 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

rsurvspline4 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


rsurvspline5 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
rsurvspline6 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

rsurvspline7 <- function(n, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp"){
    rsurvspline(n, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline)
}



hsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
hsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
hsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


hsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

hsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


hsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


hsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

hsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp"){
    hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


 
Hsurvspline0 <- function(x, gamma0, gamma1, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
Hsurvspline1 <- function(x, gamma0, gamma1, gamma2, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


Hsurvspline2 <- function(x, gamma0, gamma1, gamma2, gamma3, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


Hsurvspline3 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


Hsurvspline4 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5), knots=knots, scale=scale, timescale=timescale, spline=spline)
}

 
Hsurvspline5 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5,  gamma6, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


Hsurvspline6 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7), knots=knots, scale=scale, timescale=timescale, spline=spline)
}


Hsurvspline7 <- function(x, gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8, knots, scale="hazard", timescale="log", spline="rp"){
    Hsurvspline(x, gamma=cbind(gamma0, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, gamma7, gamma8), knots=knots, scale=scale, timescale=timescale, spline=spline)
}
