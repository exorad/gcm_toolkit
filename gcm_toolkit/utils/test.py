from clouds import cloud_opacities
import numpy as np
# ~ wave = np.array([0.30005903563735375])*1e-6
# ~ radi = np.array([5.757487933113008e-09])
# ~ abun = np.array([5.8115456995974085e-05])
# ~ rhoo = np.array([3584.7867774375186])
# ~ volu = {
# ~ 'TiO2[s]':np.asarray([5.04200118e-04]),   
# ~ 'Mg2SiO4[s]':np.asarray([1.50400305e-01]),
# ~ 'SiO[s]':np.asarray([1.03864361e-01]),   
# ~ 'SiO2[s]':np.asarray([7.78611511e-02]),  
# ~ 'Fe[s]':np.asarray([2.45392274e-02]),   
# ~ 'Al2O3[s]':np.asarray([2.16019393e-02]),  
# ~ 'CaTiO3[s]':np.asarray([9.04848529e-04]), 
# ~ 'FeO[s]':np.asarray([4.14777325e-02]),  
# ~ 'FeS[s]':np.asarray([6.29387015e-02]),   
# ~ 'Fe2O3[s]':np.asarray([5.23936554e-02]), 
# ~ 'MgO[s]':np.asarray([1.51615181e-01]),   
# ~ 'MgSiO3[s]':np.asarray([1.07987739e-01]), 
# ~ 'CaSiO3[s]':np.asarray([4.95579987e-02]),
# ~ 'Fe2SiO4[s]':np.asarray([1.17013927e-01]),
# ~ 'C[s]':np.asarray([6.13476907e-05]),   
# ~ 'KCl[s]':np.asarray([5.05006575e-03]),
# ~ }

# ~ print(cloud_opacities(wave, radi, abun, rhoo, volu))


class Dog:
    def bark(self, a):
        print('Woof!' + a)

        

border_collie = Dog()

import types

def bark(self, a, b):
    print(a+b)

border_collie = Dog()
border_collie.bark('test')
border_collie.bark  = types.MethodType(bark, border_collie)
border_collie.bark('test', '2')


