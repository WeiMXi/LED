import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from scipy.interpolate import griddata
from scipy.stats import chi2
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator

# 读取数据文件
# T2KSMdata = np.loadtxt('T2KSMprob.dat')
T2K4Ddata = np.loadtxt('T2Kprob4D.dat')
T2Kprobcr42data = np.loadtxt('T2Kprobcir2.dat')
T2Kprobcr45data = np.loadtxt('T2Kprobcir5.dat')
T2Kprobmuir01data = np.loadtxt('T2Kprobmuir0.01.dat')
T2Kprobmuir02data = np.loadtxt('T2Kprobmuir0.02.dat')
# x0=T2KSMdata[:,0]
# y0=T2KSMdata[:,1]
x1=T2K4Ddata[:,0]
y1=T2K4Ddata[:,1]
x2=T2Kprobcr42data[:,0]
y2=T2Kprobcr42data[:,1]
x3=T2Kprobcr45data[:,0]
y3=T2Kprobcr45data[:,1]
x4=T2Kprobmuir01data[:,0]
y4=T2Kprobmuir01data[:,1]
x5=T2Kprobmuir02data[:,0]
y5=T2Kprobmuir02data[:,1]
plt.figure(figsize=(8, 6)) 
plt.semilogx(x1,y1)
plt.semilogx(x2,y2)
plt.semilogx(x3,y3)
plt.semilogx(x4,y4)
plt.semilogx(x5,y5)
# plt.semilogx(x0,y0,linewidth=3.5,label=r"Standard Oscillation")
plt.semilogx(x1,y1,linewidth=2,label=r"4D finite")
plt.semilogx(x2,y2,linewidth=2,label=r"$|C_iRs|+2$")
plt.semilogx(x3,y3,linewidth=2,label=r"$|C_iRs|+5$")
plt.semilogx(x4,y4,linewidth=2,label=r"$\mu_iRs+0.01$")
plt.semilogx(x5,y5,linewidth=2,label=r"$\mu_iRs+0.02$")

plt.xlim(0.2,2.5)
# plt.ylim(0.6,1.05)
plt.xlabel("Energy/GeV",fontsize=13)
plt.ylabel("Oscillation Probability",fontsize=13)
plt.legend(prop={'size': 7})
ax = plt.gca()

# ax.xaxis.set_major_locator(MultipleLocator(0.1))
# ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.set_xticks([0.5,1,1.5,2],['0.5', '1', '0.5', '2'])
plt.tick_params(direction='in', which='both', top=True, bottom=True, left=True, right=True)
# 保存或显示结果
plt.savefig("T2Kprob.pdf", dpi=800)  # 保存为图片
plt.show()  # 显示图像
