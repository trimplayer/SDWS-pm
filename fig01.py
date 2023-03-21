import numpy as np
import matplotlib.pyplot as plt

a=1
A=1
w=3
h=10

x1=np.arange(0, w, 0.01)
n1=x1.size
y1=A*np.cos(2*np.pi/a*x1)
plt.plot(x1, y1, '-k', linewidth=2)

x2=np.arange(-0.2, w+0.5, 0.01)
n2=x2.size
y2=np.arange(n2)
y2.fill(0)
plt.plot(x2, y2, '-k', linewidth=0.7)

y3=np.arange(n1)
y3.fill(-h)
plt.plot(x1, y3, '-w', linewidth=1)
plt.fill_between(x1, y1, y3, facecolor='0.6')

y4=np.arange(-4, 3, 0.01)
n4=y4.size
x4=np.arange(n4, dtype=np.float32)
x4.fill(w/2)
plt.plot(x4, y4, '-k', linewidth=0.7)

y5=np.arange(-2.5, 0, 0.01)
n5=y5.size
x5=np.arange(n5, dtype=np.float32)
x5.fill(w/2+a)
plt.plot(x5, y5, '-k', linewidth=0.7)

x6=np.arange(w/2, w/2+a, 0.01)
n6=x6.size
y6=np.arange(n6)
y6.fill(-2)
plt.plot(x6, y6, '-k', linewidth=0.7)

plt.annotate('', xy=(w+0.3,-h+0.5),xytext=(w,-h+0.5), color='black', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),size=18)
plt.annotate('', xy=(-0.3,-h+0.5),xytext=(0,-h+0.5), color='black', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),size=18)
plt.annotate('', xy=((w+a)/2+0.3,a+2.1),xytext=((w+a)/2-0.05,a-0.05), color='black', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),size=18)

plt.annotate('', xy=(w-0.4,A-0.6),xytext=(w-0.05,A-0.05), color='black', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),size=18)
plt.annotate('', xy=(w-0.15,A+1.5),xytext=(w-0.05,A-0.1), color='black', arrowprops=dict(arrowstyle="->",connectionstyle="arc3"),size=18)

plt.text(0.1, -h+1, '$\Omega$', fontsize=18)
plt.text(0.1, 1.5, '$\Gamma$', fontsize=18)
plt.text(w/2+a/2, -3, '$a$',fontsize=18)
plt.text(w+0.4, -1.4, '$x_{1}$',fontsize=18)
plt.text(w/2-0.3, 2.5, '$x_{2}$',fontsize=18)
plt.text(w+0.15, -h+1.5, '$T$',fontsize=18)
plt.text(-0.3, -h+1.5, '$T$',fontsize=18)
plt.text(w/2-0.5, -1.8, '$-\epsilon a$', fontsize=18)
plt.text((w+a)/2+0.3, a+2, '$\mathbf{q}^{s}$', fontsize=18)
plt.text(w-0.5, A-0.25, '$\mathbf{t}$', fontsize=18, rotation=20)
plt.text(w-0.2, A+1.5, '$\mathbf{n}$', fontsize=18, rotation=20)


plt.axis('off')
#plt.savefig('KGAfig01.pdf', format='pdf', dpi=300, bbox_inches="tight")
#plt.savefig('fig01.eps', format='eps', dpi=300, bbox_inches="tight")
plt.savefig('fig01.jpg', format='jpg', dpi=300, bbox_inches="tight")
plt.show()