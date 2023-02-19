#CGL HASH for primes p = 3 mod 4

import sys

p = 2^256 * 45 - 1



F.<i> = GF(p^2, modulus=x^2 + 1)
F_vec = F.vector_space(map=False)

def c0(x):
	return F_vec(x)[0]

def c1(x):
	return F_vec(x)[1]


INV_2 = c0(F(2)^(-1))
NON_RES = c0(i^2)
INV_NON_RES = (NON_RES)^(-1)

#tonelli shanks parameters
#this is odd since P = 3 mod 4
P_MINUS_ONE_OVER_TWO = (p-1)//2
P_PLUS_ONE_OVER_FOUR = (p+1)//4
NQR_TO_P_MINUS_ONE_OVER_TWO = INV_NON_RES^P_MINUS_ONE_OVER_TWO


def tonelli_shanks(x):
	if x == 0:
		return x
	return x^P_PLUS_ONE_OVER_FOUR
	
	

def det_sqrt(x):
	x0 = c0(x)
	x1 = c1(x)
	if x1 == 0:
		if x0.is_square():
			return F(tonelli_shanks(x0))
		else:
			x1 = tonelli_shanks(x0*INV_NON_RES)
			return F(0 + x1 * i)
	alpha = x.norm()	
	alpha = tonelli_shanks(alpha)
	delta = INV_2* (x0 + alpha)
	if not(delta.is_square()):
		delta = INV_2*(x0 - alpha)
	x0 = tonelli_shanks(delta)
	x1 = (INV_2*x1)/x0
	return F(x0 + x1*i)


def step(j_i, j_i_minus_one, sign):
	j_i_sqr = j_i^2
	j_i_minus_one_sqr = j_i_minus_one^2

	a_i = -j_i_sqr + F(1488)*j_i - F(162000)
	b_i = F(1488)*j_i_sqr + F(40773375)*j_i + F(8748000000)
	D_i = (a_i + j_i_minus_one)^2 - F(4)*(b_i + a_i*j_i_minus_one + j_i_minus_one_sqr)
	S_i = det_sqrt(D_i)
	return (INV_2*(-a_i -j_i_minus_one + sign * S_i), j_i)

def CGL(walk_string, start, previous):
	assert all([True for bit in walk_string if bit == "0" or bit == "1"])
	j_i = start
	j_i_minus_one = previous
	for bit in walk_string:
		sign = F(-1) if bit == "0" else F(1)
		(j_i, j_i_minus_one) = step(j_i, j_i_minus_one, sign)
	if c1(j_i) == 0:
		print(c0(j_i))
	else:
		print(str(c0(j_i)) + " + " + str(c1(j_i)) + "*i")
		

inp = sys.argv[1]
CGL(inp, F(1728), F(1728))
