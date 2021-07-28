import numpy as np
import scipy.io
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



# Read in human brain data
brain_data = np.loadtxt("human_t_ratios.txt", dtype='f')
# array stores coverage 
coverage = brain_data[:,0]
# get rid of the first col and now matrix stores only brain data 
brain_data = brain_data[:,1:]
# transpose
brain_data = brain_data.T
    
# data
#t = np.linspace(0,100,100)
#y_ = 5 * np.sin(t/10) + 4*np.random.randn(100*10).reshape(10, 100)

t = coverage
y_ = brain_data




# here is the tsplot function for plotting the data with percentile bands

def tsplot(x, y, n=20, percentile_min=1, percentile_max=99, color='r', plot_mean=True, plot_median=False, line_color='k', **kwargs):
    # calculate the lower and upper percentile groups, skipping 50 percentile
    perc1 = np.percentile(y, np.linspace(percentile_min, 50, num=n, endpoint=False), axis=0)
    perc2 = np.percentile(y, np.linspace(50, percentile_max, num=n+1)[1:], axis=0)

    if 'alpha' in kwargs:
        alpha = kwargs.pop('alpha')
    else:
        alpha = 1/n
    # fill lower and upper percentile groups
    for p1, p2 in zip(perc1, perc2):
        plt.fill_between(x, p1, p2, alpha=alpha, color=color, edgecolor=None)


    if plot_mean:
        plt.plot(x, np.mean(y, axis=0), color=line_color)



    if plot_median:
        plt.plot(x, np.median(y, axis=0), color=line_color)
    
    return plt.gca()    


# data from simulation (different stiffness ratio @ thetag = 3)
# offset parameters
fac1 = 0.1
fac2 = 0.05
fac3 = 0.1
fac4 = 0.05
fac5 = 0.02

# transparency 
light = 0.3
med = 0.6 
dark = 1.0 


# beta=1, thetag = 3.0
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac1,np.array([1.372,1.247,1.181,1.144]),'bo',alpha=light)# A=-2
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac2,np.array([1.4122,1.2621,1.1813,1.1414]),'bv',alpha=light)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac3,np.array([1.623,1.42,1.291,1.221]),'ko',alpha=light)# A=0
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac4,np.array([1.845,1.577,1.404,1.294]),'rv',alpha=light)# A=1
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac5,np.array([2.3040,2.0357,1.7463,1.5445]),'ro',alpha=light)# A=2

# beta=3, thetag = 3.0
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac1,np.array([1.0761,1.0429,1.0217,1.0168]),'bo',alpha=med) # A=-2
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac2,np.array([1.0694,1.0631,1.0483,1.0477]),'bv',alpha=med)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac3,np.array([1.1944,1.1618,1.1268,1.1086]),'ko',alpha=med)# A=0
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac4,np.array([1.1705,1.1861,1.1521,1.1309]),'rv',alpha=med)# A=1
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac5,np.array([1.7366,1.5167,1.3035,1.2498]),'ro',alpha=med)# A=2

# beta=5, thetag = 3.0
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac1,np.array([0.8927,0.9104,0.9205,0.94]),'bo',alpha=dark) # A=-2
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac2,np.array([1.0712,1.0130,1.0087,1.0129]),'bv',alpha=dark)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac3,np.array([1.035,1.0243,1.03,1.0218]),'ko',alpha=dark)# A=0
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac4,np.array([1.3227,1.1935,1.1651,1.121]),'rv',alpha=dark)# A=1
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac5,np.array([1.6653,1.402,1.3413,1.2469]),'ro',alpha=dark)# A=2

plt.xlim((0.0, 1.05)) 
plt.ylim((0.8, 2.5)) 


plt.xlabel('Coverage',fontsize=18)
plt.ylabel('$t_{\mathrm{g}}/t_{\mathrm{s}}$',fontsize=18)
plt.title('$\mu_{\mathrm{c}}/\mu_{\mathrm{s}} = [1,3,5]\quad$ at $\quad\\vartheta^{\mathrm{g}}_{\mathrm{ctx}}=3$')

# here I call tsplot to plot the percentile band 
tsplot(t, y_, n=5, percentile_min=0.0, percentile_max=100.0, plot_median=True, plot_mean=False, color='g', line_color='g')





# data from simulation (different differerent growth levels @ stiffness ratio = 1)
# offset parameters
fac1 = -0.1
fac2 = -0.05
fac3 = 0.05
fac4 = 0.03
fac5 = 0.03

# transparency 
light = 0.2
med_light = 0.4
med_dark = 0.6
dark = 1.0 

# thetag = 2.3 
plt.plot(np.array([0.07,0.38,0.69,1.0])+2.0/4.0*fac1,np.array([1.5568,1.3860,1.2732,1.2123]),'bo',alpha=light)# A=-2 
plt.plot(np.array([0.07,0.38,0.69,1.0])+2.0/4.0*fac2,np.array([1.4967,1.3219,1.2234,1.1728]),'bv',alpha=light)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])+2.0/4.0*fac3,np.array([1.5206,1.3511,1.2392,1.1829]),'ko',alpha=light)# A=0 
plt.plot(np.array([0.07,0.38,0.69,1.0])+2.0/4.0*fac4,np.array([1.6971,1.4769,1.316,1.2387]),'rv',alpha=light)# A=1 
plt.plot(np.array([0.07,0.38,0.69,1.0])+2.0/4.0*fac5,np.array([1.7948,1.5664,1.3980,1.2935]),'ro',alpha=light)# A=2 

# thetag = 2.5
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac1,np.array([1.4646,1.2931,1.2122,1.1644]),'bo',alpha=med_light)# A=-2 
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac2,np.array([1.4528,1.3044,1.2115,1.1644]),'bv',alpha=med_light)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac3,np.array([1.5557,1.3739,1.2511,1.1926]),'ko',alpha=med_light)# A=0 
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac4,np.array([1.7368,1.5016,1.3404,1.2533]),'rv',alpha=med_light)# A=1 
plt.plot(np.array([0.07,0.38,0.69,1.0])+1.0/4.0*fac5,np.array([1.9186,1.7108,1.5325,1.3894]),'ro',alpha=med_light)# A=2 


# thetag = 2.7
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac1,np.array([1.4216,1.2712,1.1986,1.1557]),'bo',alpha=med_dark)# A=-2 
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac2,np.array([1.3932,1.2532,1.1748,1.1371]),'bv',alpha=med_dark)# A=-1 
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac3,np.array([1.6073,1.4118,1.2755,1.2113]),'ko',alpha=med_dark)# A=0 
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac4,np.array([1.7747,1.5257,1.3678,1.2702]),'rv',alpha=med_dark)# A=1 
plt.plot(np.array([0.07,0.38,0.69,1.0])-1.0/4.0*fac5,np.array([2.0656,1.8293,1.6110,1.4488]),'ro',alpha=med_dark)# A=2 


# thetag = 3.0
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac1,np.array([1.372,1.247,1.181,1.144]),'bo',alpha=dark)# A=-2 
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac2,np.array([1.4122,1.2621,1.1813,1.1414]),'bv',alpha=dark)# A=-1
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac3,np.array([1.623,1.42,1.291,1.221]),'ko',alpha=dark)# A=0
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac4,np.array([1.845,1.577,1.404,1.294]),'rv',alpha=dark)# A=1
plt.plot(np.array([0.07,0.38,0.69,1.0])-2.0/4.0*fac5,np.array([2.3040,2.0357,1.7463,1.5445]),'ro',alpha=dark)# A=2




plt.xlim((0.0, 1.07)) 
plt.ylim((0.8, 2.5)) 


plt.xlabel('Coverage',fontsize=18)
plt.ylabel('$t_{\mathrm{g}}/t_{\mathrm{s}}$',fontsize=18)
plt.title('$\\vartheta^{\mathrm{g}}_{\mathrm{ctx}}=[2.3, 2.5, 2.7, 3.0]\quad$ at $\quad\mu_{\mathrm{c}}/\mu_{\mathrm{s}} = 1$')




# here I call tsplot to plot the percentile band 
tsplot(t, y_, n=5, percentile_min=0.0, percentile_max=100.0, plot_median=True, plot_mean=False, color='g', line_color='g')

plt.savefig("Fig9.png")