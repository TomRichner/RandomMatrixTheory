If you have a system of differential equations 

\frac{dx_i}{dt}=-x_i+\sum_{j=1}^{n}J_{ij}\phi\left(x_j\right)+u_i(t)

where phi is ReLU

J is a Gaussian random matrix

and u_i(t) is a constant random vector for the first period of time such that the network reaches a fixed point x^* and then changes to a new constant random vector for the remainder and then reaches a new fixed point x^{**}

Could u_i(t) interact with phi to change the stability  of the Jacobian at these two fixed points?  

%% your reply below:
