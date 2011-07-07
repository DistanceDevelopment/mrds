apex.gamma <- function(key.scale, key.shape)
{
	fr <- (1/gamma(key.shape)) * (((key.shape - 1)/exp(1))^(key.shape - 1))
	return(key.scale*fr*(key.shape-1))
}
