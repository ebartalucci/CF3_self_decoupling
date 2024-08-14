from scipy.special import jv

def redor_bessel(max_order):
    def calculator(t, d_IS):
        if max_order > 5:
            pass  # for raising error

        z = t * 2**0.5 * d_IS

        redor_predict = 1 - jv(0, z)**2

        for order in range(1,max_order):
            redor_predict = redor_predict + 2 / (16*order**2-1) * jv(order, z)**2
    
        return redor_predict
    return calculator