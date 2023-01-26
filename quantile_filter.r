library(CVXR)

quantile_filter=function(y,Omega,lambda,tau){
n=nrow(y)
x=Variable(n)
p=Problem( Minimize( 0.5*p_norm(y-x,1)+(tau-0.5)*sum(y-x)+lambda*quad_form(x,Omega) ) )
result=solve(p)
xhat=result$getValue(x)
return(xhat)
}
