from matplotlib import pyplot as plt
import numpy as np
import numerical

ground = np.hypot(1.821, 0.4476)
l1 = ground # ground 1
l2 = 3 # crank 1
l3 = 4 # coupler 1
l4 =  3.06 # ternary 1

l5 = 3.07 # ternary 2
l6 = 3.99 # coupler 2
l7 = 3.0508 # end
l8 = ground

origin_y = 0.4476

omega_2 = -1

def getBoundedAngle(angle):
    if angle < 0:
        return 360+angle
    elif angle > 360:
        return angle-360
    else:
        return angle

def fourbarpos(a,b,c,d,th_1,th2,delta=-1):
    th_3_arr = []
    th_4_arr = []

    for i in (th2):
    
        # K constants based on link lengths
        K1 = d/a
        K2 = d/c
        K3 = (a ** 2 - b ** 2 + c ** 2 + d ** 2) / (2 * a * c)
        K4 = d/b
        K5 = (c ** 2 - d ** 2 - a ** 2 - b ** 2)/(2 * a * b)
        
        th_1 = np.deg2rad(th_1)
        th_2 = np.deg2rad(i)  # converting to rads
        A = -K1 - K2*np.cos(th_2) + K3 + np.cos(th_2)
        B = -2*np.sin(th_2)
        C = K1 - K2*np.cos(th_2) + K3 - np.cos(th_2)
        D = np.cos(th_2) - K1 +K4*np.cos(th_2) + K5
        E = -2*np.sin(th_2)
        F = K1 + (K4 - 1)*np.cos(th_2) + K5
        
        # We don't want to handle imaginary roots
        disc_4 = B**2-(4*A*C)
        disc_3 = E**2-(4*D*F)
        if disc_4 < 0 or disc_3 < 0:
            print (B,A,C,E,D,F)
            print (disc_3,disc_4)
            print('rip2')
            raise SystemExit('Error: This function does not handle imaginary roots')
        
        # Solve for thetas 
        th_4 = 2*np.arctan2((-B + delta*np.sqrt(B**2-4*A*C)),(2*A))
        th_3 = 2*np.arctan2((-E + delta*np.sqrt(E**2-4*D*F)),(2*D))

        th_3_arr.append(getBoundedAngle(np.rad2deg(th_3)))
        th_4_arr.append(getBoundedAngle(np.rad2deg(th_4)))

    return np.zeros(360), th2, np.array(th_3_arr), np.array(th_4_arr)

    
def fourbarvel(a,b,c,d,th_vec,omega_2):
    
    th_vec = [np.deg2rad(x) for x in th_vec]
    omega_3 = np.array([0.0 for _ in range(len(th_vec[1]))])
    omega_4 = np.array([0.0 for _ in range(len(th_vec[1]))])
    VA = np.zeros((len(th_vec[1]), 2))
    VBA = np.zeros((len(th_vec[1]), 2))
    VB = np.zeros((len(th_vec[1]), 2))
  
    for i in range(len(th_vec[1])):
      #print(omega_2[i]*(a/b)*(np.sin(th_vec[3][i]-th_vec[1][i]))/(np.sin(th_vec[2][i]-th_vec[3][i])))
      #print(omega_2[i]*(a/c)*(np.sin(th_vec[1][i]-th_vec[2][i]))/(np.sin(th_vec[3][i]-th_vec[2][i])))
      omega_3[i] = omega_2[i]*(a/b)*(np.sin(th_vec[3][i]-th_vec[1][i]))/(np.sin(th_vec[2][i]-th_vec[3][i]))
      omega_4[i] = omega_2[i]*(a/c)*(np.sin(th_vec[1][i]-th_vec[2][i]))/(np.sin(th_vec[3][i]-th_vec[2][i]))

      VA[i] = np.array([a*omega_2[i]*x for x in [-np.sin(th_vec[1][i]),np.cos(th_vec[1][i])]])
      VBA[i] = np.array([b*omega_3[i]*x for x in [-np.sin(th_vec[2][i]),np.cos(th_vec[2][i])]])
      VB[i] = np.array([c*omega_4[i]*x for x in [-np.sin(th_vec[3][i]),np.cos(th_vec[3][i])]])

      print(i, omega_2[i],omega_3[i],omega_4[i],VA[i],VBA[i],VB[i])

    omega_vec_out = [0,omega_2,omega_3,omega_4]
    return (omega_vec_out,VA,VBA,VB)

th1, th2, th3, th4 = (fourbarpos(l2,l3,l4,l1,0,np.linspace(0,361,360)))
th5, th6, th7, th8 = (fourbarpos(l5,l6,l7,l8,0,th4 - 39.01))

omega_2 = np.full(360, 1.0)

vel = fourbarvel(l2,l3,l4,l1,[th1, th2, th3, th4],omega_2)
print('second')
vel = fourbarvel(l5,l6,l7,l8,[th5, th6, th7, th8],vel[0][3])


'''
for i in th4:
    j, k, h = (numerical.numerical(i - 39.01,ground,3.07181611,3.99108686,3.05088475))


#print('yo')
#print(th7,th8)


plt.plot(th6, th7)
plt.plot(th6, th8)
plt.legend('a','b')
plt.show()
'''
