import numpy as np
import matplotlib.pyplot as plt
plt.ion()

def virialFactors(ScaleFactor, OmegaM0=0.3, OmegaLambda0=0.7):
  Hratio2=OmegaM0 /ScaleFactor**3+(1 - OmegaM0 - OmegaLambda0)/ScaleFactor**2+ OmegaLambda0
  OmegaZ=OmegaM0/ScaleFactor**3/Hratio2
  x=OmegaZ-1
  virialF_tophat=18.0*3.1416*3.1416+82.0*x-39.0*x*x #<Rho_vir>/Rho_cri
  virialF_c200=200.
  virialF_b200=200.*OmegaZ #virialF w.r.t contemporary critical density 
  return virialF_tophat, virialF_c200, virialF_b200


a=np.logspace(-3, 0, 20)
ft,fc,fb=virialFactors(a)

plt.loglog(a, ft, 'r')
plt.loglog(a, np.ones_like(a)*fc, 'g')
plt.loglog(a, fb, 'b--')
plt.legend(('tophat','200c','200b'),loc=3)
plt.yscale('linear')