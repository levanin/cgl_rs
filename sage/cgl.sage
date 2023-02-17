import sys
p = 2^256 * 45 - 1


F.<i> = GF(p^2, modulus=x^2 + 1)

inv_2 = F(2)^(-1)
def step(j_i, j_i_minus_one, sign):
	j_i_sqr = j_i^2
	j_i_minus_one_sqr = j_i_minus_one^2

	a_i = -j_i_sqr + F(1488)*j_i - F(162000)
	b_i = F(1488)*j_i_sqr + F(40773375)*j_i + F(8748000000)
	D_i = (a_i + j_i_minus_one)^2 - F(4)*(b_i + a_i*j_i_minus_one + j_i_minus_one_sqr)
	S_i = D_i.sqrt()
	return (inv_2*(-a_i -j_i_minus_one + sign * S_i), j_i)

def CGL(walk_string, start, previous):
	assert all([True for bit in walk_string if bit == "0" or bit == "1"])
	j_i = start
	j_i_minus_one = previous
	for bit in walk_string:
		sign = F(-1) if bit == "0" else F(1)
		(j_i, j_i_minus_one) = step(j_i, j_i_minus_one, sign)
	print("j invariant", j_i)
		

inp = sys.argv[1]
CGL(inp, F(1728), F(1728))
