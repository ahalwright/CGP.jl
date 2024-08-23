# Function to compute the density function for the normal distribution with mean mu, standard deviation sigma.
# The default values for mu and sigma are 0.0 and 1.0 
# Examples:
# normal_density( 1.0)             # (mean 0.0, std 1.0) returns 0.24197072451914337
# normal_density( 2.0, 1.0 )       # (mean 1.0, std 1.0) returns 0.24197072451914337
# normal_density( 0.0, 0.0, 2.0 )  # (mean 0.0, std 2.0) returns 0.19947114020071635

normal_density(x::Real, mu::Real=0.0, sigma::Real=1.0 ) = exp( -(x-mu)^2/(2sigma^2) )/( sigma*sqrt(2pi) )
