import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm


# update latex font library
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


# define evolution of theta_g 
def func(kappa_ast,r_cur, G_ctx):
    
    theta_g_dot = G_ctx*(1.0 - r_cur* np.tanh(kappa_ast))
    
    return theta_g_dot


# curvature-sensitive parameter 
r_cur = np.linspace(-1.0,1.0,5) 


# normalized curvature 
kappa_ast = np.linspace(-5.0,5.0,100) 

# baseline growth rate 
G_ctx = 0.5 



# length of array
length_kappa_ast = len(kappa_ast)
length_r_cur = len(r_cur)


# initialize theta_g 
theta_g = np.zeros((length_kappa_ast, length_r_cur))


# loop over kappa_ast and r_cur
for i in range(length_kappa_ast): 
    
    
    for j in range(length_r_cur): 
        
            theta_g[i][j]  = func(kappa_ast[i],r_cur[j],G_ctx)


            
# plot the results        

# plot the results 
plt.plot(kappa_ast,theta_g[:,0],'b-')
plt.plot(kappa_ast,theta_g[:,1],'b-')
plt.plot(kappa_ast,theta_g[:,2],'k-')
plt.plot(kappa_ast,theta_g[:,3],'r-')
plt.plot(kappa_ast,theta_g[:,4],'r-')

plt.text(-4.5, .025, r'$r^{\mathrm{cur}}=-1$')
plt.text(-4.5, .2, r'$r^{\mathrm{cur}}=-0.5$')
plt.text(-4.5, .42, r'$r^{\mathrm{cur}}=0$')
plt.text(-4.5, .7, r'$r^{\mathrm{cur}}=0.5$')
plt.text(-4.5, 0.95, r'$r^{\mathrm{cur}}=1$')
plt.text(3, 0.42, r'$G^{\mathrm{ctx}}$')

plt.text(-5, -0.1, 'gyri')
plt.text(5, -0.1, 'sulci')



plt.xlim((-5.0, 5.0)) 
plt.ylim((0, 1.005)) 

plt.xlabel('$\kappa^{*}$',fontsize=23)
plt.ylabel('$\dot{\\vartheta}^{\mathrm{g}}$',fontsize=23)

ax = plt.gca()
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
plt.savefig("Fig2.png")


                        
            

              
            
            
            
            
            