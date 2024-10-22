


dweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    dweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}


pweibullPH <- function(q, shape, scale = 1,
                       lower.tail=TRUE, log.p=FALSE) {
    pweibull(q, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}


qweibullPH <- function(p, shape, scale = 1, lower.tail=TRUE, log.p=FALSE) {
    qweibull(p, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}


hweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}


HweibullPH <- function(x, shape, scale=1, log=FALSE) {
    Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}


rweibullPH <- function(n, shape, scale=1) {
    rweibull(n, shape=shape, scale=scale^{-1/shape})
}


rmst_weibullPH = function(t, shape, scale=1, start=0){
  rmst_generic(pweibullPH, t, start=start, shape=shape, scale=scale)
}


mean_weibullPH = function(shape, scale=1){
  mean_weibull(shape=shape, scale=scale^{-1/shape})
}
