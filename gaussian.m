function N = gaussian(X, mu, sigma)
N = exp(-((X-mu)'*pinv(sigma)*(X-mu))/2)/((det(sigma))^(0.5)*(2*pi)^1.5);