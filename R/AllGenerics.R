setGeneric("update", function(obj, vargin) standardGeneric("update"))

setGeneric("updatewslow", function(obj) standardGeneric("updatewslow"))
setGeneric("updatehslow", function(obj) standardGeneric("updatehslow"))
setGeneric("updatethetaslow", function(obj) standardGeneric("updatethetaslow"))

setGeneric("updatew", function(obj) standardGeneric("updatew"))
setGeneric("updateh", function(obj) standardGeneric("updateh"))
setGeneric("updatetheta", function(obj) standardGeneric("updatetheta"))

setGeneric("goodk", function(obj, cutoff) standardGeneric("goodk"))
setGeneric("clearbadk", function(obj) standardGeneric("clearbadk"))

setGeneric("recomputeexpectations", 
           function(obj) standardGeneric("recomputeexpectations"))

setGeneric("figures", function(obj) standardGeneric("figures"))

setGeneric("score", function(obj) standardGeneric("score"))

setGeneric("xbar", function(obj, goodk) standardGeneric("xbar"))

setGeneric("xtwid", function(obj, goodk) standardGeneric("xtwid"))
